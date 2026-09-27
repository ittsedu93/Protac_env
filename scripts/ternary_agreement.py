#!/usr/bin/env python
"""
Os modelos aprovados apontam para a MESMA geometria ternária?

Roda no env `mdtools` (MDAnalysis). Segundos.

    python ternary_agreement.py --candidato SC0006__WH023

A proporção de aprovados não basta
-----------------------------------
3 de 25 modelos passando o corte de interface pode significar duas coisas
opostas:

  - três geometrias DIFERENTES  -> ruído; o modelo não sabe onde as proteínas
                                   se encontram e acertou o critério por acaso
  - três vezes a MESMA          -> o modelo encontra aquela interface
                                   raramente, mas de forma reprodutível — o
                                   que é evidência, não sorte

A diferença é medível: superpõe-se cada modelo pela E3 (que é a âncora
experimental, do cristal) e mede-se o quanto o ALVO se desloca. Modelos que
concordam na interface ficam próximos; modelos que discordam, longe.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


def ler_ca(pdb: Path):
    """CA por cadeia, na ordem do arquivo."""
    por_cadeia = {}
    for l in pdb.read_text().splitlines():
        if not l.startswith("ATOM") or l[12:16].strip() != "CA":
            continue
        if l[16] not in (" ", "A"):
            continue
        cad = l[21]
        por_cadeia.setdefault(cad, []).append(
            (float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return {c: np.array(v) for c, v in por_cadeia.items()}


def kabsch(movel, ref):
    """Rotação+translação que leva `movel` a `ref`; devolve a matriz aplicada."""
    cm, cr = movel.mean(0), ref.mean(0)
    H = (movel - cm).T @ (ref - cr)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, cm, cr


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidato", required=True)
    ap.add_argument("--work", type=Path,
                    default=Path.home() / "PRosettaC_runs" / "vhl_crbn_pcsk9_protac")
    ap.add_argument("--so-aprovados", action="store_true",
                    help="compara só os que passaram o corte de interface")
    args = ap.parse_args()

    base = args.work.expanduser() / "boltz_ternario" / args.candidato
    rank = base / "ranking_poses.csv"
    if not rank.exists():
        raise SystemExit(f"rode antes o select_ternary_pose.py "
                         f"(não achei {rank.name})")
    linhas = list(csv.DictReader(open(rank)))
    if args.so_aprovados:
        linhas = [l for l in linhas if l["aprovado"] == "True"]
    if len(linhas) < 2:
        raise SystemExit("preciso de pelo menos 2 modelos para comparar")

    ids = json.loads((base / "entradas.json").read_text())["boltz_ids"]
    letras = list(ids.keys())
    papel = {c: ("E3" if ids[c].startswith("E3")
                 else "PROTAC" if ids[c].startswith("PROTAC") else "ALVO")
             for c in letras}
    cad_e3 = [c for c in letras if papel[c] == "E3"]
    cad_alvo = [c for c in letras if papel[c] == "ALVO"]

    # acha o pdb de cada modelo
    def pdb_de(nome):
        n = nome.split("/")[-1]
        for p in base.rglob(f"{n}.pdb"):
            if "/lote_" in str(p) or "/" not in nome:
                if nome.startswith("lote_") and f"/{nome.split('/')[0]}/" not in str(p):
                    continue
                if not nome.startswith("lote_") and "/lote_" in str(p):
                    continue
                return p
        return None

    estruturas = {}
    for l in linhas:
        p = pdb_de(l["modelo"])
        if p is None:
            print(f"  [pulado] não achei o .pdb de {l['modelo']}")
            continue
        estruturas[l["modelo"]] = ler_ca(p)

    nomes = list(estruturas)
    if len(nomes) < 2:
        raise SystemExit("não consegui ler .pdb suficientes")

    print(f"candidato: {args.candidato}")
    print(f"  E3: cadeias {cad_e3} | alvo: cadeias {cad_alvo}")
    print(f"  {len(nomes)} modelos comparados"
          + ("  (só os aprovados)" if args.so_aprovados else ""))
    print()

    def concat(d, cads):
        return np.vstack([d[c] for c in cads if c in d])

    # matriz: superpõe pela E3, mede o deslocamento do alvo
    M = np.zeros((len(nomes), len(nomes)))
    for i, a in enumerate(nomes):
        for j, b in enumerate(nomes):
            if i >= j:
                continue
            ea, eb = concat(estruturas[a], cad_e3), concat(estruturas[b], cad_e3)
            ta, tb = concat(estruturas[a], cad_alvo), concat(estruturas[b], cad_alvo)
            if ea.shape != eb.shape or ta.shape != tb.shape:
                M[i, j] = M[j, i] = np.nan
                continue
            R, cm, cr = kabsch(ea, eb)
            ta_al = (ta - cm) @ R.T + cr
            M[i, j] = M[j, i] = float(np.sqrt(np.mean(
                np.sum((ta_al - tb) ** 2, axis=1))))

    larg = max(len(n) for n in nomes) + 2
    print(" " * larg + "  ".join(f"{i:>6d}" for i in range(len(nomes))))
    for i, n in enumerate(nomes):
        print(f"{i:>2d} {n:<{larg-3}}" +
              "  ".join("     ." if i == j else
                        ("   n/a" if np.isnan(M[i, j]) else f"{M[i, j]:6.1f}")
                        for j in range(len(nomes))))

    fora = M[~np.eye(len(nomes), dtype=bool)]
    fora = fora[~np.isnan(fora)]
    print(f"\n  RMSD do alvo após superpor pela E3 (Å):")
    print(f"    menor {fora.min():.1f} | mediana {np.median(fora):.1f} "
          f"| maior {fora.max():.1f}")
    proximos = int((fora <= 5.0).sum() / 2)
    total = len(nomes) * (len(nomes) - 1) // 2
    print(f"    {proximos} de {total} pares a menos de 5 Å")
    print()
    if fora.min() > 10:
        print("  -> NENHUM par concorda. Cada modelo põe o alvo num lugar")
        print("     diferente: o Boltz não encontrou uma interface, encontrou")
        print("     várias incompatíveis. O critério de ipTM foi satisfeito")
        print("     por acaso.")
    elif proximos >= total / 2:
        print("  -> Os modelos CONCORDAM na geometria. Uma interface")
        print("     reprodutível entre amostras independentes é evidência real,")
        print("     mesmo que a proporção de aprovados seja baixa.")
    else:
        print("  -> Concordância PARCIAL: há um subgrupo que converge e outros")
        print("     que divergem. O subgrupo que concorda é o que vale — os")
        print("     pares a menos de 5 Å acima dizem quais são.")

    with open(base / "concordancia_poses.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([""] + nomes)
        for i, n in enumerate(nomes):
            w.writerow([n] + [f"{M[i, j]:.2f}" for j in range(len(nomes))])
    print(f"\n  {base / 'concordancia_poses.csv'}")


if __name__ == "__main__":
    main()
