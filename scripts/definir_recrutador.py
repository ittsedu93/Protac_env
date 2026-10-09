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


def coords_receptor(pdb: Path) -> np.ndarray:
    xyz = [[float(l[30:38]), float(l[38:46]), float(l[46:54])]
           for l in pdb.read_text(errors="ignore").splitlines()
           if l.startswith("ATOM")]
    if not xyz:
        raise SystemExit(f"nenhum registro ATOM em {pdb}")
    return np.array(xyz)


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
    args = ap.parse_args()

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
    ev = exit_vector_do_recrutador(mol, rec, burial=args.burial, mol_ref=None)
    print(f"\n[2/5] vetor de saída")
    print(f"      átomo {ev['atom_idx']} ({ev['atom_symbol']},"
          f" {ev['atom_n_hs']} H), {ev['dist_ao_nucleo_A']} Å do núcleo")
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
