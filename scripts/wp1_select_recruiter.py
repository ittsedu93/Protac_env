#!/usr/bin/env python
"""
Fase 0b — escolhe o recrutador da E3 ligase e extrai o exit vector da pose.

Roda no env `mdtools`. Minutos.

Esta é a peça que liga o WP1 ao WP2. O WP1 produziu rankings de triagem
(`screening/<E3>/<E3>_round2_results.csv`) e poses, mas nada disso vira entrada
do WP2 sem duas coisas que precisam ser derivadas da pose vencedora:

  1. o **ponto de conjugação** — o átomo do recrutador por onde o linker sai,
     que precisa estar no índice 0 do SDF (mesma convenção dos warheads, que é
     o que `ReplaceSubstructs(replacementConnectionPoint=0)` espera);
  2. o **exit vector** — ponto e direção, no referencial do receptor, para
     posicionar os confôrmeros do linker.

Ambos saem do mesmo cálculo de enterramento usado no ligante `063` da PCSK9:
átomos com muitos vizinhos proteicos formam o núcleo ancorado, átomos sem
vizinho nenhum apontam para o solvente.

    python wp1_select_recruiter.py \\
        --screening ~/PRosettaC_runs/vhl_crbn_pcsk9_protac/screening \\
        --prep ~/PRosettaC_runs/vhl_crbn_pcsk9_protac/prep \\
        --out ~/PRosettaC_runs/vhl_crbn_pcsk9_protac/pipeline/wp1_recruiter.json

ATENÇÃO: a escolha do átomo de conjugação é uma HEURÍSTICA geométrica (o átomo
exposto mais afastado do núcleo, com hidrogênio para ceder). Ela não sabe
química de acoplamento. O script grava um .cxc para você conferir no ChimeraX,
e o número escolhido deve ser validado antes de comprometer o WP2 inteiro.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

OBABEL_EXE = "/home/soberano/miniconda3/envs/obabel_env/bin/obabel"

# Átomos que NÃO podem receber o linker: são o farmacóforo que ancora o
# recrutador na E3 ligase. Conjugar neles destrói o recrutamento, e a
# heurística geométrica não tem como saber disso — ela vê um átomo exposto com
# hidrogênio e o escolhe.
#
# O caso que motivou esta lista: num análogo de talidomida o ÚNICO N com
# hidrogênio é o NH da glutarimida, que faz as ligações de hidrogênio no bolsão
# tri-triptofano da CRBN. A heurística o elegeu ponto de conjugação.
FARMACOFOROS = [
    ("NH da glutarimida (CRBN — bolsão tri-Trp)", "[NX3;H1;R](C(=O))C(=O)"),
    ("N imídico (CRBN)",                          "[NX3;R](C(=O))C(=O)"),
    ("hidroxila da hidroxiprolina (VHL)",         "[OX2;H1][CX4;R][CX4;R]"),
    ("amida da terc-leucina (VHL)",               "[NX3;H1]C(=O)[CX4][CX4](C)(C)C"),
]


def atomos_farmacoforicos(mol):
    """Índices que não devem receber o linker, com o motivo."""
    proibidos = {}
    for nome, smarts in FARMACOFOROS:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue
        for match in mol.GetSubstructMatches(patt):
            proibidos.setdefault(match[0], nome)
    return proibidos


def _esqueleto(mol):
    """Cópia com todas as ligações simples e sem aromaticidade."""
    m = Chem.RWMol(mol)
    for b in m.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    for a in m.GetAtoms():
        a.SetIsAromatic(False)
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(0)
        a.SetFormalCharge(0)
    out = m.GetMol()
    try:
        Chem.SanitizeMol(out, Chem.SanitizeFlags.SANITIZE_SYMMRINGS |
                              Chem.SanitizeFlags.SANITIZE_ADJUSTHS)
    except Exception:
        pass
    return out


def farmacoforos_da_pose(pose, ref):
    """Farmacóforos da pose, calculados na molécula de referência.

    A pose vem de PDBQT, formato que não guarda ordem de ligação: na volta
    para SDF as carbonilas viram ligações simples e os SMARTS de farmacóforo
    deixam de casar. Pior que ficar cega, a guarda passa a produzir FALSOS
    POSITIVOS — num análogo de talidomida com as ordens perdidas ela "detecta"
    a hidroxila da hidroxiprolina da VHL.

    Então os padrões são aplicados ao SDF preparado, que tem as ordens
    corretas, e os índices são transferidos para a pose pelo casamento das
    duas estruturas.
    """
    if ref is None or ref.GetNumAtoms() != pose.GetNumAtoms():
        return atomos_farmacoforicos(pose), "pose (sem referência)", None

    proibidos_ref = atomos_farmacoforicos(ref)
    if not proibidos_ref:
        return {}, "referência preparada", 0

    # ref -> pose: match[i] é o índice na pose do átomo i da referência.
    # O casamento é pelo ESQUELETO (conectividade e elemento, ordens de ligação
    # zeradas), porque é exatamente a ordem de ligação que a pose perdeu.
    match = (pose.GetSubstructMatch(ref)
             or _esqueleto(pose).GetSubstructMatch(_esqueleto(ref)))
    if not match:
        return (atomos_farmacoforicos(pose),
                "pose (as estruturas não casaram — guarda pouco confiável)",
                len(proibidos_ref))
    return ({match[i]: motivo for i, motivo in proibidos_ref.items()
             if i < len(match)}, "referência preparada", len(proibidos_ref))

# nomes prováveis das colunas nos CSVs do WP1, em ordem de preferência
COLS_ID = ("ligand_id", "ligand", "name", "id", "mol_id", "molecule")
COLS_SCORE = ("best_score", "score", "affinity", "vina_score", "energy")


def achar_coluna(df: pd.DataFrame, candidatas) -> str | None:
    baixo = {c.lower(): c for c in df.columns}
    for c in candidatas:
        if c in baixo:
            return baixo[c]
    return None


def ler_ranking(csv_path: Path):
    df = pd.read_csv(csv_path)
    col_id = achar_coluna(df, COLS_ID)
    col_sc = achar_coluna(df, COLS_SCORE)
    if col_id is None or col_sc is None:
        raise SystemExit(
            f"não reconheci as colunas de {csv_path.name}.\n"
            f"  colunas presentes: {list(df.columns)}\n"
            f"  esperava um id entre {COLS_ID} e um score entre {COLS_SCORE}.\n"
            f"  Renomeie no CSV ou me diga os nomes corretos.")
    return df.sort_values(col_sc), col_id, col_sc


def achar_pose(raiz: Path, lig_id: str) -> Path | None:
    """As poses do WP1 podem estar em vários layouts; procura os prováveis."""
    padroes = [
        f"**/round2/{lig_id}/*.pdbqt", f"**/round2/{lig_id}*.pdbqt",
        f"**/{lig_id}/run*.pdbqt", f"**/{lig_id}*.pdbqt",
        f"**/{lig_id}/*.sdf", f"**/{lig_id}*.sdf",
    ]
    for p in padroes:
        achados = sorted(Path(raiz).glob(p))
        if achados:
            return melhor_por_score(achados)
    return None


def melhor_por_score(arquivos: list[Path]) -> Path:
    melhor, melhor_sc = arquivos[0], np.inf
    for f in arquivos:
        if f.suffix != ".pdbqt":
            continue
        try:
            for linha in f.read_text().splitlines():
                if linha.startswith("REMARK VINA RESULT"):
                    sc = float(linha.split()[3])
                    if sc < melhor_sc:
                        melhor_sc, melhor = sc, f
                    break
        except (OSError, ValueError, IndexError):
            continue
    return melhor


def corrigir_ordens(pose, ref):
    """Devolve a molécula com a química do SDF preparado e as coordenadas
    do docking.

    O PDBQT não guarda ordem de ligação nem hidrogênio. A pose lida de volta
    tem os átomos com valência travada, e daí saem dois estragos encadeados:
    um nitrogênio de amida vira amina e estoura a valência ao receber o linker
    (as 75 falhas de montagem do WP3), e os carbonos ficam sem hidrogênio
    nenhum — radicais, que o antechamber rejeita com "the number of electrons
    is odd".

    Consertar a pose não funciona: AssignBondOrdersFromTemplate sozinho deixa
    33 elétrons radicalares, e liberar o hidrogênio implícito depois disso
    destrói a aromaticidade. Então o caminho é o inverso: parte-se da
    molécula PREPARADA, que tem a química certa, e transplantam-se as
    coordenadas da pose para ela.
    """
    if ref is None:
        return pose, "sem referência (química NÃO conferida)"
    try:
        # O obabel preserva os hidrogênios polares do PDBQT, então a pose
        # pode ter mais átomos que a molécula preparada. Remove H dos dois
        # lados antes de qualquer comparação.
        alvo = Chem.RemoveAllHs(Chem.Mol(ref), sanitize=False)
        pose = Chem.RemoveAllHs(Chem.Mol(pose), sanitize=False)
        if alvo.GetNumAtoms() != pose.GetNumAtoms():
            return pose, (f"contagem de átomos pesados difere (preparada "
                          f"{alvo.GetNumAtoms()}, pose {pose.GetNumAtoms()})")

        # 1) casamento exato, quando a pose preservou a química
        match = alvo.GetSubstructMatch(pose)

        # 2) ligações genéricas: é a ordem de ligação que a pose perdeu, então
        #    a comparação precisa ignorá-la
        if not match:
            try:
                par = Chem.AdjustQueryParameters.NoAdjustments()
                par.makeBondsGeneric = True
                match = alvo.GetSubstructMatch(
                    Chem.AdjustQueryProperties(pose, par))
            except Exception:
                match = ()

        # 3) identidade: pose e preparada são a MESMA molécula, e a conversão
        #    sdf -> pdbqt -> sdf preserva a ordem dos átomos pesados. Só é
        #    aceita depois de conferir elemento por elemento.
        if not match:
            if ([a.GetSymbol() for a in alvo.GetAtoms()]
                    == [a.GetSymbol() for a in pose.GetAtoms()]):
                match = tuple(range(alvo.GetNumAtoms()))

        if not match or len(match) != pose.GetNumAtoms():
            return pose, "estruturas não casaram (química NÃO conferida)"

        pconf = pose.GetConformer()
        novo_conf = Chem.Conformer(alvo.GetNumAtoms())
        for i_pose, i_alvo in enumerate(match):
            novo_conf.SetAtomPosition(i_alvo, pconf.GetAtomPosition(i_pose))
        alvo.RemoveAllConformers()
        alvo.AddConformer(novo_conf, assignId=True)
        Chem.SanitizeMol(alvo)

        rad = sum(a.GetNumRadicalElectrons() for a in alvo.GetAtoms())
        if rad:
            return pose, f"transplante deixou {rad} radicais — descartado"
        return alvo, "coordenadas transplantadas para a molécula preparada"
    except Exception as exc:
        return pose, f"falha ({type(exc).__name__}: {exc})"


def carregar_mol(caminho: Path, tmp: Path):
    caminho = Path(caminho)
    if caminho.suffix == ".pdbqt":
        tmp.mkdir(parents=True, exist_ok=True)
        sdf = tmp / (caminho.stem + ".sdf")
        if not sdf.exists():
            subprocess.run([OBABEL_EXE, str(caminho), "-O", str(sdf)],
                           check=True, capture_output=True)
        caminho = sdf
    if caminho.suffix in (".sdf", ".mol"):
        return Chem.MolFromMolFile(str(caminho), removeHs=True)
    if caminho.suffix == ".pdb":
        return Chem.MolFromPDBFile(str(caminho), removeHs=True)
    return None


def identificar_no_catalogo(lig_id: str, mol, anchors_sdf: Path | None,
                            prep: Path | None = None):
    """Descobre QUAL composto do catálogo é o recrutador escolhido.

    Sem isto o resultado é "ligand_042", que não serve para encomendar reagente
    nem para escrever na tese.

    A correspondência por POSIÇÃO no SDF não é confiável: o WP1 pode ter
    filtrado a biblioteca antes de numerar, e foi o que aconteceu — o índice 87
    do catálogo tem 22 átomos pesados contra 37 da pose de ligand_087. Então a
    identificação é feita pela ESTRUTURA: lê-se o SDF preparado em
    prep/ligands/<lig_id>.sdf e casa-se o InChIKey contra o catálogo.
    """
    if not anchors_sdf or not Path(anchors_sdf).exists():
        return {"catalogo": None}

    from rdkit.Chem import inchi

    # estrutura de referência: o SDF preparado é mais fiel que a pose,
    # que passou por pdbqt e perdeu ordens de ligação
    ref = None
    if prep:
        for cand in (Path(prep) / "ligands" / f"{lig_id}.sdf",
                     Path(prep) / f"{lig_id}.sdf"):
            if cand.exists():
                ref = Chem.MolFromMolFile(str(cand), removeHs=True)
                if ref is not None:
                    break
    if ref is None:
        ref = mol

    try:
        chave_ref = inchi.MolToInchiKey(ref).split("-")[0]   # só o esqueleto
    except Exception:
        chave_ref = None

    supp = Chem.SDMolSupplier(str(anchors_sdf), removeHs=True)
    achado, idx_achado = None, None
    for i, cand in enumerate(supp):
        if cand is None:
            continue
        if cand.GetNumHeavyAtoms() != ref.GetNumHeavyAtoms():
            continue
        try:
            if chave_ref and inchi.MolToInchiKey(cand).split("-")[0] == chave_ref:
                achado, idx_achado = cand, i
                break
        except Exception:
            continue

    if achado is None:
        return {"catalogo": None,
                "catalogo_nota": (
                    f"não achei {lig_id} no catálogo por estrutura "
                    f"({ref.GetNumHeavyAtoms()} átomos pesados). A biblioteca "
                    f"triada no WP1 pode não ser exatamente este SDF.")}

    props = {k: achado.GetProp(k) for k in achado.GetPropNames()}
    nome = (achado.GetProp("_Name") if achado.HasProp("_Name") else "") or ""
    for chave in ("IDNUMBER", "ID", "Catalog ID", "CSID", "Chemspace ID"):
        if not nome and chave in props:
            nome = props[chave]

    return {
        "catalogo": {
            "indice_no_sdf": idx_achado,
            "nome": nome or f"(sem nome, índice {idx_achado})",
            "smiles": Chem.MolToSmiles(achado),
            "n_heavy": achado.GetNumHeavyAtoms(),
            "casado_por": "InChIKey do esqueleto",
            "propriedades": {k: v for k, v in list(props.items())[:10]},
        },
        "catalogo_nota": "",
    }


def receptor_coords(receptor_pdb: Path) -> np.ndarray:
    return np.array([(float(l[30:38]), float(l[38:46]), float(l[46:54]))
                     for l in Path(receptor_pdb).read_text().splitlines()
                     if l.startswith("ATOM") and (l[76:78].strip() or "C") != "H"])


# ---------------------------------------------------------------------------
def exit_vector_do_recrutador(mol, rec_xyz: np.ndarray, burial: int = 20,
                              raio: float = 6.0, mol_ref=None):
    """Ponto de conjugação e direção de saída, pelo enterramento dos átomos.

    A divisão núcleo/exposto usa percentis em vez de cortes absolutos: um
    recrutador raso pode não ter nenhum átomo com 20 vizinhos, e um corte fixo
    faria os dois conjuntos coincidirem — o centroide de um igual ao do outro,
    e a direção indefinida.
    """
    pos = mol.GetConformer().GetPositions()
    d = np.linalg.norm(pos[:, None, :] - rec_xyz[None, :, :], axis=2)
    nb = (d <= raio).sum(axis=1)

    if nb.max() == 0:
        raise SystemExit(
            f"a pose não toca o receptor (nenhum átomo com vizinho a {raio} Å). "
            f"Pose e receptor estão no mesmo referencial?")

    if nb.max() == nb.min():
        raise SystemExit("enterramento uniforme: não dá para definir a saída")

    # terço mais enterrado contra terço mais exposto, sempre disjuntos
    alto = np.percentile(nb, 66)
    baixo = np.percentile(nb, 33)
    enterrados = nb >= max(alto, nb.min() + 1)
    expostos = nb <= baixo
    if enterrados.sum() < 2 or expostos.sum() < 2 or (enterrados & expostos).all():
        ordem = np.argsort(nb)
        k = max(2, len(nb) // 3)
        expostos = np.zeros(len(nb), bool); expostos[ordem[:k]] = True
        enterrados = np.zeros(len(nb), bool); enterrados[ordem[-k:]] = True

    centro_nucleo = pos[enterrados].mean(axis=0)
    centro_exposto = pos[expostos].mean(axis=0)
    v = centro_exposto - centro_nucleo
    n = np.linalg.norm(v)
    if n < 1e-6:
        raise SystemExit("núcleo e porção exposta coincidem: pose degenerada?")
    direcao = v / n

    # ponto de conjugação: entre os átomos expostos COM hidrogênio para ceder,
    # o mais afastado do núcleo ancorado — excluindo o farmacóforo
    proibidos, origem_guarda, n_ref = farmacoforos_da_pose(mol, mol_ref)
    candidatos = [int(i) for i in np.where(expostos)[0]
                  if mol.GetAtomWithIdx(int(i)).GetTotalNumHs() > 0]
    excluidos = [(i, proibidos[i]) for i in candidatos if i in proibidos]
    candidatos = [i for i in candidatos if i not in proibidos]

    if not candidatos:
        livres = [int(i) for i in np.where(expostos)[0] if int(i) not in proibidos]
        if not livres:
            raise SystemExit(
                "todos os átomos expostos fazem parte do farmacóforo: este "
                "recrutador não tem vetor de saída utilizável sem destruir o "
                "reconhecimento pela E3 ligase")
        candidatos = livres
    idx = max(candidatos, key=lambda i: float(np.linalg.norm(pos[i] - centro_nucleo)))

    atomo = mol.GetAtomWithIdx(int(idx))
    return {
        "atom_idx": int(idx),
        "atom_symbol": atomo.GetSymbol(),
        "atom_n_hs": atomo.GetTotalNumHs(),
        "exit_point": [round(float(x), 3) for x in pos[idx]],
        "exit_direction": [round(float(x), 4) for x in direcao],
        "dist_ao_nucleo_A": round(float(np.linalg.norm(pos[idx] - centro_nucleo)), 2),
        "vizinhos_proteicos": int(nb[idx]),
        "n_atomos_enterrados": int(enterrados.sum()),
        "n_atomos_expostos": int(expostos.sum()),
        "burial_min_max": [int(nb.min()), int(nb.max())],
        "excluidos_por_farmacoforo": [{"atom_idx": i, "motivo": m}
                                      for i, m in excluidos],
        "farmacoforos_no_recrutador": len(proibidos),
        "guarda_calculada_em": origem_guarda,
    }


def reordenar_conjugacao_primeiro(mol, idx: int):
    ordem = [idx] + [i for i in range(mol.GetNumAtoms()) if i != idx]
    return Chem.RenumberAtoms(mol, ordem)


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--screening", type=Path, required=True)
    ap.add_argument("--prep", type=Path, required=True)
    ap.add_argument("--e3", nargs="*", default=["VHL", "CRBN"])
    ap.add_argument("--burial", type=int, default=20)
    ap.add_argument("--anchors-sdf", type=Path,
                    help="biblioteca de anchors usada no WP1, para identificar "
                         "o composto de catálogo correspondente")
    ap.add_argument("--top-n", type=int, default=3,
                    help="quantos do topo inspecionar por E3")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    screening = args.screening.expanduser()
    prep = args.prep.expanduser()
    out = args.out.expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.parent / "_tmp_recruiter"

    resultados = {}
    for e3 in args.e3:
        print(f"\n{'=' * 60}\n{e3}\n{'=' * 60}")
        dir_e3 = screening / e3
        csvs = sorted(dir_e3.glob("*round2*.csv")) or sorted(dir_e3.glob("*.csv"))
        if not csvs:
            print(f"  [pulado] nenhum CSV de triagem em {dir_e3}")
            continue
        df, col_id, col_sc = ler_ranking(csvs[0])
        print(f"  ranking: {csvs[0].name} ({len(df)} ligantes, "
              f"id='{col_id}', score='{col_sc}')")

        # receptor preparado no WP1
        subdirs = [d for d in prep.iterdir() if d.is_dir() and d.name.startswith(e3)]
        if not subdirs:
            print(f"  [pulado] nenhuma pasta {e3}* em {prep}")
            continue
        rec_pdbs = sorted(subdirs[0].glob("*_receptor.pdb"))
        if not rec_pdbs:
            print(f"  [pulado] receptor não encontrado em {subdirs[0]}")
            continue
        rec_pdb = rec_pdbs[0]
        rec_xyz = receptor_coords(rec_pdb)
        print(f"  receptor: {rec_pdb.name} ({len(rec_xyz)} átomos pesados)")

        for _, linha in df.head(args.top_n).iterrows():
            lig_id = str(linha[col_id])
            score = float(linha[col_sc])
            pose = achar_pose(dir_e3, lig_id) or achar_pose(prep, lig_id)
            if pose is None:
                print(f"  [{lig_id}] score {score:.2f} — pose não encontrada")
                continue
            mol = carregar_mol(pose, tmp)
            if mol is None:
                print(f"  [{lig_id}] pose ilegível: {pose}")
                continue
            ref_prep = None
            for cand in (prep / "ligands" / f"{lig_id}.sdf", prep / f"{lig_id}.sdf"):
                if cand.exists():
                    ref_prep = Chem.MolFromMolFile(str(cand), removeHs=True)
                    if ref_prep is not None:
                        break
            mol, nota_ordens = corrigir_ordens(mol, ref_prep)
            print(f"      ordens de ligação: {nota_ordens}")
            try:
                ev = exit_vector_do_recrutador(mol, rec_xyz, args.burial,
                                               mol_ref=ref_prep)
            except SystemExit as exc:
                print(f"  [{lig_id}] {exc}")
                continue

            print(f"  [{lig_id}] score {score:.2f} | conjugação pelo átomo "
                  f"{ev['atom_idx']} ({ev['atom_symbol']}, {ev['atom_n_hs']} H) "
                  f"a {ev['dist_ao_nucleo_A']} Å do núcleo")
            print(f"      farmacóforos: {ev['farmacoforos_no_recrutador']} "
                  f"(via {ev['guarda_calculada_em']})")
            for ex in ev.get("excluidos_por_farmacoforo", []):
                print(f"      [protegido] átomo {ex['atom_idx']}: {ex['motivo']}")
            if ev["farmacoforos_no_recrutador"] == 0:
                print(f"      [ATENÇÃO] nenhum farmacóforo reconhecido neste "
                      f"recrutador: a guarda não está protegendo nada. "
                      f"Confira manualmente antes de aceitar.")

            if e3 in resultados:
                continue                       # guarda só o melhor por E3

            rad = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
            if rad:
                print(f"      [DESCARTADO] {rad} elétrons radicalares: a "
                      f"química desta pose não pôde ser reconstruída, e ela "
                      f"quebraria o antechamber lá na MD")
                continue

            sdf = out.parent / f"recruiter_{e3}_{lig_id}.sdf"
            with Chem.SDWriter(str(sdf)) as w:
                w.write(reordenar_conjugacao_primeiro(mol, ev["atom_idx"]))

            cxc = out.parent / f"conferir_{e3}.cxc"
            cxc.write_text(
                f"open {rec_pdb}\nopen {pose}\n"
                f"show #2 atoms\nstyle #2 stick\nlabel #2 atoms\n"
                f"# o átomo de conjugação escolhido é o índice "
                f"{ev['atom_idx']} do SDF ({ev['atom_symbol']})\n"
                f"# confira se ele aponta para o solvente e se é um ponto de\n"
                f"# acoplamento quimicamente razoável\n")

            cat = identificar_no_catalogo(lig_id, mol, args.anchors_sdf, prep)
            if cat.get("catalogo"):
                c = cat["catalogo"]
                print(f"      catálogo: {c['nome']} (índice {c['indice_no_sdf']}, "
                      f"casado por {c['casado_por']})")
            if cat.get("catalogo_nota"):
                print(f"      {cat['catalogo_nota']}")

            resultados[e3] = {
                "e3": e3, "ligand_id": lig_id, "score": score,
                **cat,
                "pose": str(pose), "recruiter_sdf": str(sdf),
                "receptor_pdb": str(rec_pdb),
                "receptor_pdbqt": str(rec_pdb.with_suffix(".pdbqt")),
                "ref_ligand_sdf": next(
                    (str(p) for p in subdirs[0].glob("*ref_ligand*.sdf")), ""),
                "redock_poses": str(subdirs[0] / "redock*.pdbqt"),
                "conferir_cxc": str(cxc),
                **ev,
            }

    if not resultados:
        raise SystemExit(
            "\nNenhum recrutador pôde ser selecionado. Veja as mensagens acima: "
            "normalmente é o CSV com nomes de coluna diferentes dos esperados, "
            "ou as poses do round 2 num layout que o script não procurou.")

    # escolhe a E3 com melhor score
    melhor = min(resultados.values(), key=lambda d: d["score"])
    saida = {"escolhido": melhor, "por_e3": resultados}
    out.write_text(json.dumps(saida, indent=2))

    print(f"\n{'=' * 60}")
    print(f"ESCOLHIDO: {melhor['e3']} / {melhor['ligand_id']} "
          f"(score {melhor['score']:.2f})")
    print(f"{'=' * 60}")
    if melhor.get("catalogo"):
        print(f"  composto   : {melhor['catalogo']['nome']}")
        print(f"  SMILES     : {melhor['catalogo']['smiles']}")
    print(f"  recrutador : {melhor['recruiter_sdf']}")
    print(f"  exit point : {melhor['exit_point']}")
    print(f"  direção    : {melhor['exit_direction']}")
    print(f"\n  CONFIRA a escolha do átomo de conjugação antes do WP2:")
    print(f"    chimerax {melhor['conferir_cxc']}")
    print(f"\n  A heurística é geométrica: o átomo exposto mais afastado do")
    print(f"  núcleo ancorado, com hidrogênio para ceder. Ela não sabe se esse")
    print(f"  ponto é sinteticamente acessível.")
    print(f"\n{out}")


if __name__ == "__main__":
    main()
