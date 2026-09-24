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


def identificar_no_catalogo(lig_id: str, mol, anchors_sdf: Path | None):
    """Descobre QUAL composto do catálogo é o recrutador escolhido.

    Sem isto o resultado é "ligand_042", que não serve para encomendar nada
    nem para escrever na tese. O WP1 nomeou os ligantes por posição na
    biblioteca (ligand_NNN), então NNN é o índice no SDF; a conferência por
    contagem de átomos pesados garante que o índice não escorregou.
    """
    if not anchors_sdf or not Path(anchors_sdf).exists():
        return {"catalogo": None}
    try:
        idx = int(str(lig_id).rsplit("_", 1)[-1])
    except ValueError:
        return {"catalogo": None, "catalogo_nota": "id não termina em número"}

    supp = Chem.SDMolSupplier(str(anchors_sdf), removeHs=True)
    if idx >= len(supp):
        return {"catalogo": None,
                "catalogo_nota": f"índice {idx} fora da biblioteca ({len(supp)})"}
    cand = supp[idx]
    if cand is None:
        return {"catalogo": None, "catalogo_nota": f"molécula {idx} ilegível"}

    props = {k: cand.GetProp(k) for k in cand.GetPropNames()}
    nome = (cand.GetProp("_Name") if cand.HasProp("_Name") else "") or ""
    for chave in ("IDNUMBER", "ID", "Catalog ID", "CSID", "Chemspace ID", "SMILES"):
        if chave in props and not nome:
            nome = props[chave]

    bate = cand.GetNumHeavyAtoms() == mol.GetNumHeavyAtoms()
    return {
        "catalogo": {
            "indice_no_sdf": idx,
            "nome": nome or f"(sem nome, índice {idx})",
            "smiles": Chem.MolToSmiles(cand),
            "n_heavy_catalogo": cand.GetNumHeavyAtoms(),
            "n_heavy_pose": mol.GetNumHeavyAtoms(),
            "confere": bool(bate),
            "propriedades": {k: v for k, v in list(props.items())[:10]},
        },
        "catalogo_nota": ("" if bate else
                          "ATENÇÃO: contagem de átomos pesados difere entre a "
                          "pose e a molécula do catálogo nesse índice — a "
                          "correspondência por posição pode estar errada"),
    }


def receptor_coords(receptor_pdb: Path) -> np.ndarray:
    return np.array([(float(l[30:38]), float(l[38:46]), float(l[46:54]))
                     for l in Path(receptor_pdb).read_text().splitlines()
                     if l.startswith("ATOM") and (l[76:78].strip() or "C") != "H"])


# ---------------------------------------------------------------------------
def exit_vector_do_recrutador(mol, rec_xyz: np.ndarray, burial: int = 20,
                              raio: float = 6.0):
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
    # o mais afastado do núcleo ancorado
    candidatos = [int(i) for i in np.where(expostos)[0]
                  if mol.GetAtomWithIdx(int(i)).GetTotalNumHs() > 0]
    if not candidatos:
        candidatos = [int(i) for i in np.where(expostos)[0]]
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
            try:
                ev = exit_vector_do_recrutador(mol, rec_xyz, args.burial)
            except SystemExit as exc:
                print(f"  [{lig_id}] {exc}")
                continue

            print(f"  [{lig_id}] score {score:.2f} | conjugação pelo átomo "
                  f"{ev['atom_idx']} ({ev['atom_symbol']}, {ev['atom_n_hs']} H) "
                  f"a {ev['dist_ao_nucleo_A']} Å do núcleo")

            if e3 in resultados:
                continue                       # guarda só o melhor por E3

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

            cat = identificar_no_catalogo(lig_id, mol, args.anchors_sdf)
            if cat.get("catalogo"):
                c = cat["catalogo"]
                print(f"      catálogo: {c['nome']}  "
                      f"({'confere' if c['confere'] else 'NÃO CONFERE'}: "
                      f"{c['n_heavy_pose']} vs {c['n_heavy_catalogo']} átomos)")
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
