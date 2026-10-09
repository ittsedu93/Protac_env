#!/usr/bin/env python
"""
Troca o recrutador da E3 com as verificações encadeadas, de uma vez.

Roda no env `mdtools` (RDKit). Instantâneo.

    python definir_recrutador.py --e3 VHL \\
        --ligante .../6GFZ_ref_ligand.sdf \\
        --receptor .../6GFZ_receptor.pdb \\
        --out .../recruiter_VHL_6GFZ_cristal.sdf

Por que um script e não cinco comandos
--------------------------------------
Trocar o recrutador toca cinco coisas, e cada uma tem um modo de falhar
silencioso. Em cinco comandos manuais, esquecer um não dá erro — dá resultado
errado com cara de certo, que é o que custou sete dias a este projeto.

    1. quimiotipo        o ligante recruta a E3 que declaramos?
    2. vetor de saída    há ponto de conjugação FORA do farmacóforo?
    3. ÁTOMO 0           o WP2 conjuga o linker no átomo de índice 0. Um SDF de
                         cristal tem a ordem do PDB, não a nossa convenção.
                         Se o átomo 0 por acaso tiver um H, o WP2 NÃO reclama:
                         ele pendura o linker ali — possivelmente no meio da
                         hidroxiprolina, que é o farmacóforo da VHL — e os 120
                         PROTACs saem errados sem uma única mensagem
    4. coordenadas       renumerar não pode mover átomo nenhum: a pose é
                         cristalográfica e é o ativo mais valioso aqui
    5. marcadores        as fases 5 a 9 já estão marcadas como feitas. Sem
                         apagar os marcadores, o driver PULA tudo e roda a
                         fase 10 sobre os PROTACs antigos

Os três primeiros são portões: falhando, o script para. Os dois últimos ele
confere e relata.

Nada é reimplementado: o quimiotipo vem do checar_recrutador.py, o vetor de
saída do wp1_select_recruiter.py e a renumeração do generate_pcsk9_warheads.py.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parent))
from checar_recrutador import ASSINATURAS, motivos  # noqa: E402
from generate_pcsk9_warheads import reorder_attachment_first  # noqa: E402
from wp1_select_recruiter import exit_vector_do_recrutador  # noqa: E402

RDLogger.DisableLog("rdApp.*")

DESVIO_MAX_A = 1e-3     # renumerar não move átomo; isto é folga de float

# A QUE GRUPO o ponto de conjugação pertence. O filtro de farmacóforo do WP1 é
# genérico — ele procura doadores/aceptores de H e cargas — e por isso NÃO
# distingue uma metila de tiazol (ponto de conjugação usado em PROTACs de VHL
# publicados) de uma metila de terc-butilo (elemento de ligação enterrado no
# bolso). As duas são CH3 sem farmacóforo nenhum, e escolher a segunda destrói
# o reconhecimento sem disparar guarda alguma.
# Classifico só o que sei classificar. O resto sai como NAO RECONHECIDO, e isso
# NÃO é um passe: a lista de pontos de conjugação válidos não é enumerável por
# mim, então o que eu não reconheço fica para o olho humano no SMILES marcado.
# Tratar desconhecido como "ok" seria um passe silencioso sobre justamente o
# que eu não previ — o mesmo mecanismo que deixou um ligante de CRBN passar por
# recrutador de VHL.
GRUPOS_DO_PONTO = [
    ("metila de terc-butilo", "[CH3]C([CH3])([CH3])", "PROIBIDO"),
    ("metila substituinte de anel aromático", "[CH3]c", "ok"),
    ("carbono alifático periférico", "[CH3,CH2][CH2,CH1,NX3,OX2]", "ok"),
    ("metila de acetamida (N-cap do scaffold)", "[CH3]C(=O)N", "ATENCAO"),
    ("hidroxila", "[OX2H]", "ATENCAO"),
    ("posição C-H de anel aromático", "[cH]", "ATENCAO"),
]


def ambiente_do_atomo0(mol):
    """(linhas, pior_risco) descrevendo a que grupo o átomo 0 pertence."""
    a0 = mol.GetAtomWithIdx(0)
    linhas = [f"átomo 0: {a0.GetSymbol()}, {a0.GetTotalNumHs()} H, "
              f"grau {a0.GetDegree()}, aromático={a0.GetIsAromatic()}"]
    for n in a0.GetNeighbors():
        linhas.append(f"  vizinho: {n.GetSymbol()} (grau {n.GetDegree()}, "
                      f"anel={n.IsInRing()}, aromático={n.GetIsAromatic()})")
    marc = Chem.Mol(mol)
    marc.GetAtomWithIdx(0).SetAtomMapNum(1)
    linhas.append(f"  SMILES com o ponto marcado [:1]:")
    linhas.append(f"    {Chem.MolToSmiles(marc)}")
    pior, achados = None, []
    for nome, sma, risco in GRUPOS_DO_PONTO:
        q = Chem.MolFromSmarts(sma)
        if q is None:
            continue
        if any(0 in m for m in mol.GetSubstructMatches(q)):
            achados.append((nome, risco))
            if risco == "PROIBIDO":
                pior = "PROIBIDO"
            elif risco == "ATENCAO" and pior != "PROIBIDO":
                pior = "ATENCAO"
            elif pior is None:
                pior = "ok"
    if achados:
        linhas.append("  grupos que contêm o átomo 0: "
                      + ", ".join(f"{n} [{r}]" for n, r in achados))
    else:
        pior = "NAO RECONHECIDO"
        linhas.append("  o átomo 0 não casa com nenhum grupo que eu saiba"
                      " classificar.")
        linhas.append("  Isto NÃO é aprovação: leia o SMILES marcado acima e"
                      " confirme que o")
        linhas.append("  ponto é periférico e não participa do reconhecimento"
                      " pela E3.")
    return linhas, pior


def coords_receptor(pdb: Path) -> np.ndarray:
    xyz = [[float(l[30:38]), float(l[38:46]), float(l[46:54])]
           for l in pdb.read_text(errors="ignore").splitlines()
           if l.startswith("ATOM")]
    if not xyz:
        raise SystemExit(f"nenhum registro ATOM em {pdb}")
    return np.array(xyz)


def listar_pontos(mol, rec, raio=6.0):
    """Todo átomo pesado com H, classificado e com o enterramento MEDIDO.

    Por que isto existe
    -------------------
    A heurística do WP1 escolhe o átomo exposto mais distante do núcleo
    enterrado, excluindo farmacóforo. Para o ligante do cristal da VHL ela
    escolheu uma metila do TERC-BUTILO — que é elemento de ligação do scaffold
    VH032, e que nenhum SMARTS de farmacóforo marca, porque uma metila não tem
    farmacóforo.

    Num scaffold cuja química é conhecida, quem decide o vetor de saída é a
    química, não o enterramento. Então aqui não se escolhe: MEDE-SE e lista-se,
    e a escolha é explícita (--ponto-idx). O enterramento entra como dado da
    pose cristalográfica, não como critério único.
    """
    pos = mol.GetConformer().GetPositions()
    d = np.linalg.norm(pos[:, None, :] - rec[None, :, :], axis=2)
    nb = (d <= raio).sum(axis=1)
    centro = pos[np.argsort(nb)[-max(2, len(nb) // 3):]].mean(axis=0)

    linhas = []
    for a in mol.GetAtoms():
        if a.GetTotalNumHs() < 1:
            continue
        i = a.GetIdx()
        risco, grupo = "ok", "—"
        for nome, sma, r in GRUPOS_DO_PONTO:
            q = Chem.MolFromSmarts(sma)
            if q is not None and any(i in mm for mm in mol.GetSubstructMatches(q)):
                if r == "PROIBIDO" or (r == "ATENCAO" and risco == "ok"):
                    risco, grupo = r, nome
                elif risco == "ok":
                    grupo = nome
        viz = ",".join(n.GetSymbol() + ("(ar)" if n.GetIsAromatic() else "")
                       for n in a.GetNeighbors())
        linhas.append({
            "idx": i, "el": a.GetSymbol(), "nH": a.GetTotalNumHs(),
            "arom": a.GetIsAromatic(), "vizinhos": viz,
            "vizinhos_proteicos": int(nb[i]),
            "dist_nucleo": float(np.linalg.norm(pos[i] - centro)),
            "grupo": grupo, "risco": risco,
        })
    # mais exposto primeiro; proibidos no fim
    linhas.sort(key=lambda r: (r["risco"] == "PROIBIDO",
                              r["vizinhos_proteicos"], -r["dist_nucleo"]))
    return linhas


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--e3", required=True)
    ap.add_argument("--ligante", type=Path, required=True)
    ap.add_argument("--receptor", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--e3-chain", default=None,
                    help="cadeia da E3 no receptor, para a linha da config")
    ap.add_argument("--pipeline-out", type=Path, default=None,
                    help="$PIPELINE_OUT: lista os marcadores .done_* obsoletos")
    ap.add_argument("--burial", type=int, default=20)
    ap.add_argument("--so-ambiente", type=Path, default=None,
                    help="só descreve a que grupo pertence o átomo 0 de um SDF "
                         "já escrito, e sai")
    ap.add_argument("--aceitar-ponto-proibido", action="store_true")
    ap.add_argument("--listar-pontos", action="store_true",
                    help="lista os candidatos a ponto de conjugação com o "
                         "enterramento medido na pose, e sai")
    ap.add_argument("--ponto-smarts", default=None,
                    help="escolhe o ponto de conjugação por QUÍMICA, não por "
                         "posição: o SMARTS tem de casar exatamente um átomo "
                         "com H. Reprodutível entre SDFs de ordem diferente, "
                         "que é o que a tese precisa poder citar")
    ap.add_argument("--ponto-idx", type=int, default=None,
                    help="escolhe o ponto de conjugação EXPLICITAMENTE, pelo "
                         "índice no SDF de entrada. A direção de saída passa a "
                         "ser o vetor da ligação vizinho->átomo, que é a "
                         "definição química de vetor de saída")
    args = ap.parse_args()

    if args.so_ambiente:
        m = next(iter(Chem.SDMolSupplier(str(args.so_ambiente.expanduser()),
                                         removeHs=True)), None)
        if m is None:
            raise SystemExit(f"não li {args.so_ambiente}")
        print(f"{args.so_ambiente.name}")
        linhas, pior = ambiente_do_atomo0(m)
        for l in linhas:
            print("  " + l)
        print(f"\n  risco: {pior}")
        return 0 if pior != "PROIBIDO" else 2

    e3 = args.e3.upper()
    lig_p = args.ligante.expanduser()
    mol = next(iter(Chem.SDMolSupplier(str(lig_p), removeHs=True)), None)
    if mol is None or mol.GetNumConformers() == 0:
        raise SystemExit(f"não li {lig_p} com coordenadas")
    print(f"ligante : {lig_p.name} — {mol.GetNumAtoms()} átomos pesados")
    print(f"receptor: {args.receptor.name}")

    # --- PORTÃO 1: quimiotipo ---------------------------------------------
    ach = motivos(mol)
    sufs = {a for a, l in ach.items() if any(s for _, s in l)}
    print(f"\n[1/5] quimiotipo")
    for alvo, lista in sorted(ach.items()):
        print(f"      {alvo}: " + (", ".join(n for n, _ in lista) or "nenhuma"))
    if sufs - {e3}:
        raise SystemExit(
            f"\n*** recruta {', '.join(sorted(sufs - {e3}))}, e declaramos {e3}."
            f"\n*** Foi exatamente este erro que anulou sete dias de PRosettaC.")
    if e3 not in sufs:
        raise SystemExit(
            f"\n*** sem assinatura suficiente de {e3}. Pode ser um ligante"
            f" legítimo fora\n*** do quimiotipo canônico, mas isso precisa de"
            f" justificativa explícita\n*** antes de custar dias de máquina.")
    print(f"      OK: assinatura de {e3} e de mais nenhuma E3")

    # --- PORTÃO 2: vetor de saída -----------------------------------------
    rec = coords_receptor(args.receptor.expanduser())

    if args.listar_pontos:
        print("\ncandidatos a ponto de conjugação, mais EXPOSTO primeiro")
        print("(enterramento = átomos de proteína a 6 Å, medido NESTA pose)\n")
        print(f"  {'idx':>4} {'el':<4} {'nH':>2} {'viz.prot':>9} "
              f"{'d.núcleo':>9}  risco      grupo / vizinhos")
        for r in listar_pontos(mol, rec):
            mark = "  <-- PROIBIDO" if r["risco"] == "PROIBIDO" else ""
            print(f"  {r['idx']:>4} {r['el'] + ('(ar)' if r['arom'] else ''):<4} "
                  f"{r['nH']:>2} {r['vizinhos_proteicos']:>9} "
                  f"{r['dist_nucleo']:>9.2f}  {r['risco']:<10} "
                  f"{r['grupo']} [{r['vizinhos']}]{mark}")
        print("\n  Escolha um com --ponto-idx N. O mais exposto e quimicamente")
        print("  periférico é o candidato natural; a decisão é química, e por")
        print("  isso ela é sua e não da heurística.")
        return 0

    if args.ponto_smarts:
        q = Chem.MolFromSmarts(args.ponto_smarts)
        if q is None:
            raise SystemExit(f"SMARTS inválido: {args.ponto_smarts}")
        # O primeiro átomo do SMARTS é o ponto de conjugação, e ele tem de ter
        # H. Casar vários átomos com H seria ambíguo, e ambiguidade aqui vira
        # um PROTAC ligado num lugar que ninguém escolheu.
        alvos = sorted({m[0] for m in mol.GetSubstructMatches(q)
                        if mol.GetAtomWithIdx(m[0]).GetTotalNumHs() >= 1})
        if not alvos:
            raise SystemExit(
                f"o SMARTS '{args.ponto_smarts}' não casou nenhum átomo com H."
                f"\n  Rode --listar-pontos para ver os candidatos.")
        if len(alvos) > 1:
            raise SystemExit(
                f"o SMARTS '{args.ponto_smarts}' casou {len(alvos)} átomos com"
                f" H: {alvos}.\n  Precisa ser exatamente um — torne o padrão"
                f" mais específico, ou use --ponto-idx.")
        args.ponto_idx = alvos[0]
        print(f"\n      SMARTS '{args.ponto_smarts}' -> átomo {alvos[0]}")

    if args.ponto_idx is not None:
        i = args.ponto_idx
        if not (0 <= i < mol.GetNumAtoms()):
            raise SystemExit(f"índice {i} fora de 0..{mol.GetNumAtoms()-1}")
        a = mol.GetAtomWithIdx(i)
        if a.GetTotalNumHs() < 1:
            raise SystemExit(
                f"o átomo {i} ({a.GetSymbol()}) não tem H para ceder ao linker")
        pos = mol.GetConformer().GetPositions()
        vz = [n.GetIdx() for n in a.GetNeighbors()]
        if not vz:
            raise SystemExit(f"o átomo {i} não tem vizinho: direção indefinida")
        # A definição química de vetor de saída: a direção da ligação que o
        # linker vai substituir, do átomo pesado vizinho para o ponto de
        # conjugação. É mais precisa que (centroide exposto - centroide
        # enterrado), que é uma média de toda a molécula.
        v = pos[i] - pos[vz[0]]
        v = v / float(np.linalg.norm(v))
        nb = (np.linalg.norm(pos[:, None, :] - rec[None, :, :], axis=2)
              <= 6.0).sum(axis=1)
        ev = {"atom_idx": i, "atom_symbol": a.GetSymbol(),
              "atom_n_hs": a.GetTotalNumHs(),
              "exit_point": [round(float(x), 3) for x in pos[i]],
              "exit_direction": [round(float(x), 4) for x in v],
              "dist_ao_nucleo_A": 0.0,
              "vizinhos_proteicos": int(nb[i]),
              "n_atomos_enterrados": 0, "n_atomos_expostos": 0,
              "burial_min_max": [int(nb.min()), int(nb.max())],
              "excluidos_por_farmacoforo": [],
              "origem": f"escolha explícita --ponto-idx {i}; direção = vetor "
                        f"da ligação {a.GetSymbol()}{i} <- vizinho {vz[0]}"}
        print(f"\n[2/5] vetor de saída (ESCOLHA EXPLÍCITA)")
        print(f"      {ev['origem']}")
    else:
        ev = exit_vector_do_recrutador(mol, rec, burial=args.burial,
                                       mol_ref=None)
        print(f"\n[2/5] vetor de saída (heurística do WP1)")
    print(f"      átomo {ev['atom_idx']} ({ev['atom_symbol']},"
          f" {ev['atom_n_hs']} H), {ev['vizinhos_proteicos']} vizinhos proteicos")
    exc = ev.get("excluidos_por_farmacoforo") or []
    if exc:
        print(f"      {len(exc)} átomo(s) excluído(s) por farmacóforo:")
        for e in exc[:5]:
            print(f"          átomo {e['atom_idx']}: {e['motivo']}")
    if ev["atom_symbol"] == "H":
        raise SystemExit(
            "\n*** escolheu um HIDROGÊNIO como ponto de conjugação, o que nunca"
            "\n*** é válido: o linker entra no lugar de um H, num átomo pesado."
            "\n*** Nenhum átomo pesado exposto tem H livre fora do farmacóforo.")
    if ev["atom_n_hs"] < 1:
        raise SystemExit(
            f"\n*** o átomo {ev['atom_idx']} não tem H para ceder ao linker.")
    print(f"      OK: átomo pesado, com H, fora do farmacóforo")

    # --- PORTÃO 3 + 4: átomo 0 e coordenadas ------------------------------
    idx = int(ev["atom_idx"])
    xyz_antes = mol.GetConformer().GetPositions()
    novo = reorder_attachment_first(mol, idx)
    print(f"\n[3/5] renumeração: átomo {idx} -> índice 0")
    a0 = novo.GetAtomWithIdx(0)
    if a0.GetTotalNumHs() < 1:
        raise SystemExit(
            f"\n*** depois de renumerar, o átomo 0 ({a0.GetSymbol()}) não tem H."
            f"\n*** É exatamente a condição que o WP2 exige, e ela falhou.")
    print(f"      átomo 0 = {a0.GetSymbol()}, {a0.GetTotalNumHs()} H — "
          f"é o que o WP2 exige")
    linhas, pior = ambiente_do_atomo0(novo)
    for l in linhas:
        print("      " + l)
    if pior == "PROIBIDO" and not args.aceitar_ponto_proibido:
        raise SystemExit(
            "\n*** o ponto de conjugação está num TERC-BUTILO. Esse grupo é"
            " elemento\n*** de ligação do scaffold VH032, enterrado no bolso da"
            " VHL — pendurar o\n*** linker nele destrói o reconhecimento, e o"
            " filtro de farmacóforo não\n*** pega isso porque uma metila não"
            " tem farmacóforo nenhum.\n*** \n"
            "*** Escolha outro ponto, ou --aceitar-ponto-proibido se souber"
            " por quê.")
    if pior in ("ATENCAO", "NAO RECONHECIDO"):
        print(f"      [ATENÇÃO] risco {pior}: confira o SMILES marcado acima"
              f" antes de gastar dias de máquina")

    xyz_depois = novo.GetConformer().GetPositions()
    # a ordem mudou, então compara CONJUNTO de posições, não posição a posição
    perm = [idx] + [i for i in range(mol.GetNumAtoms()) if i != idx]
    desvio = float(np.abs(xyz_depois - xyz_antes[perm]).max())
    print(f"\n[4/5] coordenadas: desvio máximo {desvio:.2e} Å")
    if desvio > DESVIO_MAX_A:
        raise SystemExit(
            f"\n*** renumerar MOVEU átomos ({desvio:.3f} Å). A pose é"
            f" cristalográfica —\n*** é o ativo mais valioso aqui, e não pode"
            f" mudar.")
    p0 = xyz_depois[0]
    if float(np.abs(p0 - np.array(ev["exit_point"])).max()) > 1e-2:
        raise SystemExit(
            f"\n*** o átomo 0 não está no exit point. {p0} vs"
            f" {ev['exit_point']}")
    print(f"      OK: nada se moveu, e o átomo 0 está no exit point")

    out = args.out.expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    with Chem.SDWriter(str(out)) as w:
        w.write(novo)
    (out.with_suffix(".exit_vector.json")).write_text(
        json.dumps(ev, indent=2, ensure_ascii=False))

    # --- 5: o que a config precisa, e o que está obsoleto -----------------
    print(f"\n[5/5] escrito: {out}")
    print(f"       {out.with_suffix('.exit_vector.json').name}")
    print("\n" + "=" * 66)
    print("AS LINHAS DA CONFIG (substitua as existentes em pipeline_vhl.conf)")
    print("=" * 66)
    print(f'E3_NAME="{e3}"')
    print(f'E3_RECRUITER_SDF="{out}"')
    print(f'E3_RECEPTOR_PDB="{args.receptor.expanduser()}"')
    if args.e3_chain:
        print(f'E3_CHAIN="{args.e3_chain}"')
    print(f'E3_EXIT_POINT="{",".join(str(x) for x in ev["exit_point"])}"')
    print(f'E3_EXIT_DIRECTION="{",".join(str(x) for x in ev["exit_direction"])}"')
    print("=" * 66)
    print("Preencher a config NÃO é opcional: com elas vazias o driver lê o")
    print("recrutador do wp1_recruiter.json — que tem o ligante errado. O `:=`")
    print("do bash substitui em variável vazia, não só em variável ausente.")

    if args.pipeline_out:
        po = args.pipeline_out.expanduser()
        obsoletos = [f".done_{n}" for n in ("5", "6", "6a", "6b", "7", "8", "9", "10")
                     if (po / f".done_{n}").exists()]
        print(f"\nMARCADORES OBSOLETOS em {po}:")
        if not obsoletos:
            print("  nenhum — nada a apagar")
        else:
            for m in obsoletos:
                print(f"  {m}")
            print("\n  Todos dependem do recrutador. Sem apagá-los o driver")
            print("  PULA essas fases e roda a MD sobre os PROTACs antigos.")
            print("  A 6a entra na lista porque a restrição do PatchDock é")
            print("  entre os dois ÁTOMOS DE CONJUGAÇÃO, e o do lado da E3")
            print("  mudou de lugar — a curva do vão muda com ele.")
            print("\n  rm -f " + " ".join(f'"{po}/{m}"' for m in obsoletos))
        rj = po / "wp1_recruiter.json"
        if rj.exists():
            print(f"\n  E saia da frente do json obsoleto:")
            print(f"    mv \"{rj}\" \"{rj}.recrutador_errado\"")
    return 0


if __name__ == "__main__":
    sys.exit(main())
