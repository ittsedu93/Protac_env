#!/usr/bin/env python
"""
Etapa 3.0d — análise da triagem dos warheads na PCSK9.

Roda no env `mdtools` (RDKit + pandas). Segundos.

    python scripts/analyze_warhead_docking.py \\
        --docking ~/PCSK9_docking/docking \\
        --site ~/PCSK9_docking/pcsk9_site.json \\
        --warheads-csv ~/PCSK9_warheads/pcsk9_warheads.csv

Três perguntas que o score sozinho não responde
------------------------------------------------

1. **O ranking é artefato de lipofilia?** A função de scoring do Vina soma
   contribuições por átomo, então favorece sistematicamente ligantes maiores e
   mais lipofílicos. Numa série que varia de 229 a 420 Da, ordenar por score
   bruto tende a ranquear por tamanho. A eficiência de ligante
   (LE = -score / átomos pesados) corrige isso, e comparar os dois rankings
   mostra o tamanho do viés.

2. **A pose é estável?** Score bom com desvio alto entre runs independentes é
   candidato a falso positivo: o motor não converge para a mesma solução.

3. **O linker tem por onde sair?** O N de conjugação precisa apontar para o
   solvente. Se ficar enterrado no bolsão, o candidato morre no WP2 — e é
   melhor descobrir antes de montar PROTAC com ele.

Saídas
------
    warhead_ranking.csv   tabela completa com score, LE, consistência e
                          acessibilidade do vetor de saída
    (stdout)              enriquecimento por grupo R e recomendação
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


# ---------------------------------------------------------------------------
def receptor_heavy_coords(receptor_pdb: Path) -> np.ndarray:
    return np.array([(float(l[30:38]), float(l[38:46]), float(l[46:54]))
                     for l in Path(receptor_pdb).read_text().splitlines()
                     if l.startswith("ATOM") and (l[76:78].strip() or "C") != "H"])


def find_conjugation_nitrogen(mol):
    """N da etilamina — o ponto por onde o linker entra.

    As poses vêm de PDBQT convertido, onde a percepção de aromaticidade e de
    hidrogênios pode divergir do SDF original, então a busca é por camadas, da
    mais específica para a mais tolerante.
    """
    for smarts in ("c[CH2][CH2][NX3]", "[#6][CH2][CH2][NX3]",
                   "[#6][#6][#6][NX3;!R]"):
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue
        m = mol.GetSubstructMatches(patt)
        if m:
            return m[0][-1]
    # último recurso: qualquer N fora de anel
    fora = [a.GetIdx() for a in mol.GetAtoms()
            if a.GetSymbol() == "N" and not a.IsInRing()]
    return fora[0] if fora else None


def exit_vector_check(head_sdf: Path, site: dict, rec_xyz: np.ndarray,
                      min_clearance: float = 3.0):
    mol = Chem.MolFromMolFile(str(head_sdf), removeHs=True)
    if mol is None:
        return {"n_encontrado": False}
    idx = find_conjugation_nitrogen(mol)
    if idx is None:
        return {"n_encontrado": False}

    n_pos = np.array(mol.GetConformer().GetAtomPosition(idx))
    clearance = float(np.linalg.norm(rec_xyz - n_pos, axis=1).min())

    v = n_pos - np.array(site["center"])
    nv = np.linalg.norm(v)
    cos = (float(np.dot(v / nv, np.array(site["exit_vector"])))
           if nv > 1e-6 and site.get("exit_vector") else np.nan)

    return {"n_encontrado": True,
            "n_clearance_A": round(clearance, 2),
            "cos_vetor_saida": round(cos, 3) if np.isfinite(cos) else np.nan,
            "n_exposto": clearance >= min_clearance,
            "alinhado_com_063": bool(np.isfinite(cos) and cos > 0.3)}


# ---------------------------------------------------------------------------
def enriquecimento(df: pd.DataFrame, coluna_r: str, top_n: int = 20):
    """Compara a composição do topo por score e por LE com a do conjunto."""
    base = df[coluna_r].value_counts(normalize=True)
    top_score = df.nsmallest(top_n, "best_score")[coluna_r].value_counts(normalize=True)
    top_le = df.nlargest(top_n, "ligand_efficiency")[coluna_r].value_counts(normalize=True)

    tab = pd.DataFrame({
        "no_conjunto_%": (base * 100).round(1),
        f"top{top_n}_score_%": (top_score * 100).round(1),
        f"top{top_n}_LE_%": (top_le * 100).round(1),
    }).fillna(0.0)
    tab["viés_score"] = (tab[f"top{top_n}_score_%"] - tab["no_conjunto_%"]).round(1)
    return tab.sort_values("viés_score", ascending=False)


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--docking", type=Path, required=True)
    ap.add_argument("--site", type=Path, required=True)
    ap.add_argument("--warheads-csv", type=Path, required=True)
    ap.add_argument("--sd-max", type=float, default=0.5,
                    help="desvio máximo entre runs para a pose contar como "
                         "estável (kcal/mol, default 0.5)")
    ap.add_argument("--top-n", type=int, default=20)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    docking = args.docking.expanduser()
    site = json.loads(args.site.expanduser().read_text())

    r2 = pd.read_csv(docking / "round2_scores.csv")
    heads = pd.DataFrame(json.loads((docking / "heads_manifest.json").read_text()))
    wh = pd.read_csv(args.warheads_csv.expanduser())

    df = (r2.merge(heads, left_on="ligand_id", right_on="warhead_id",
                   how="inner", suffixes=("", "_head"))
            .merge(wh, on="warhead_id", how="left"))
    print(f"{len(df)} warheads com pose no round 2\n")

    # --- 1. eficiência de ligante ---
    df["ligand_efficiency"] = (-df["best_score"] / df["heavy_atoms"]).round(3)

    # --- 2. consistência ---
    df["pose_estavel"] = df["std_score"] <= args.sd_max

    # --- 3. vetor de saída ---
    rec_xyz = receptor_heavy_coords(Path(site["receptor_pdb"]))
    checks = [exit_vector_check(Path(p), site, rec_xyz) for p in df["head_sdf"]]
    df = pd.concat([df.reset_index(drop=True),
                    pd.DataFrame(checks)], axis=1)

    nao_lidos = int((~df["n_encontrado"]).sum())
    if nao_lidos:
        print(f"[aviso] {nao_lidos} poses sem N de conjugação identificável\n",
              file=sys.stderr)

    # ---------------- relatório ----------------
    print("=" * 66)
    print("1. O ranking por score é artefato de tamanho/lipofilia?")
    print("=" * 66)
    corr = df["best_score"].corr(df["heavy_atoms"])
    print(f"\ncorrelação score x átomos pesados: {corr:+.2f}")
    if corr < -0.5:
        print("  -> FORTE: score melhora com o tamanho. Ordenar por score bruto")
        print("     ranqueia por massa molecular, não por complementaridade.")
    elif corr < -0.3:
        print("  -> MODERADA: o viés existe e deve ser reportado.")
    else:
        print("  -> fraca: o ranking por score não é dominado por tamanho.")

    for coluna in ("r1", "r2", "r3"):
        if coluna in df.columns:
            print(f"\n--- grupo {coluna.upper()} ---")
            print(enriquecimento(df, coluna, args.top_n).to_string())

    print("\n" + "=" * 66)
    print("2. Consistência das poses")
    print("=" * 66)
    print(f"\n{df['pose_estavel'].sum()}/{len(df)} com desvio <= "
          f"{args.sd_max} kcal/mol entre os runs")
    instaveis = df[~df["pose_estavel"]].nsmallest(5, "best_score")
    if len(instaveis):
        print("\nmelhores scores COM pose instável (candidatos a falso positivo):")
        print(instaveis[["warhead_id", "best_score", "std_score"]].to_string(index=False))

    print("\n" + "=" * 66)
    print("3. Vetor de saída do linker")
    print("=" * 66)
    ok_exp = df["n_exposto"].fillna(False)
    ok_ali = df["alinhado_com_063"].fillna(False)
    print(f"\n{ok_exp.sum()}/{len(df)} com o N de conjugação exposto ao solvente")
    print(f"{ok_ali.sum()}/{len(df)} alinhados com o braço do 063 (cos > 0,3)")
    print(f"{(ok_exp & ok_ali).sum()}/{len(df)} satisfazem os dois")
    if (~ok_exp).sum():
        print(f"\n{(~ok_exp).sum()} candidatos têm o N enterrado: não há por onde")
        print("o linker sair, e eles morreriam no WP2.")

    # ---------------- priorização ----------------
    viaveis = df[ok_exp & ok_ali & df["pose_estavel"]].copy()
    print("\n" + "=" * 66)
    print(f"4. Priorização — {len(viaveis)}/{len(df)} passam nos três critérios")
    print("=" * 66)

    if len(viaveis):
        viaveis = viaveis.sort_values("ligand_efficiency", ascending=False)
        cols = [c for c in ("warhead_id", "r1", "r2", "r3", "mw",
                            "best_score", "ligand_efficiency", "std_score",
                            "n_clearance_A", "cos_vetor_saida", "clogp",
                            "net_charge_ph74", "filter_flags")
                if c in viaveis.columns]
        print("\nTop 15 por eficiência de ligante:")
        print(viaveis[cols].head(15).to_string(index=False))

        print("\nPara comparação, top 15 por score bruto:")
        print(viaveis.nsmallest(15, "best_score")[cols].to_string(index=False))

        sobrep = len(set(viaveis.head(15)["warhead_id"])
                     & set(viaveis.nsmallest(15, "best_score")["warhead_id"]))
        print(f"\nsobreposição entre os dois top-15: {sobrep}/15")
        if sobrep < 8:
            print("  -> os dois critérios discordam bastante. A escolha entre")
            print("     eles muda o que vai para o WP2, então registre qual")
            print("     você usou e por quê.")

    out = args.out or (docking / "warhead_ranking.csv")
    df.sort_values("ligand_efficiency", ascending=False).to_csv(out, index=False)
    print(f"\ntabela completa: {out}")


if __name__ == "__main__":
    main()
