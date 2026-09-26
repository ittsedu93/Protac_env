#!/usr/bin/env python
"""
A série temporal por trás das médias: deriva ou platô?

Só usa numpy/pandas. Segundos, lê md_series.json do disco.

O problema que isto resolve
--------------------------
Uma média sobre 200 ns inclui a relaxação inicial, quando o sistema ainda está
saindo das restrições de posição do equilíbrio. Um RMSD que sobe e estabiliza
em 3 Å e um que sobe sem parar até 5 Å podem ter a MESMA média — e significam
coisas opostas. A média não distingue; a série distingue.

    python md_series_report.py --md-dir ~/pipeline/md
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def blocos(v, n=8):
    """Média em n blocos iguais, para ver a forma da curva num relance."""
    v = np.asarray(v, dtype=float)
    corte = np.array_split(v, n)
    return [float(np.mean(c)) for c in corte if len(c)]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--md-dir", type=Path, required=True)
    ap.add_argument("--ns", type=float, default=200.0)
    ap.add_argument("--metade", action="store_true",
                    help="também resume só a segunda metade da trajetória")
    args = ap.parse_args()

    md = args.md_dir.expanduser()
    series = json.loads((md / "md_series.json").read_text())

    for chave, rotulo in (("rmsd_prot", "RMSD da proteína (Å)"),
                          ("rmsd_lig", "RMSD do PROTAC (Å)"),
                          ("c_rec", "contatos recrutador"),
                          ("c_wh", "contatos warhead")):
        if not any(chave in s for s in series.values()):
            continue
        print(f"\n{rotulo} — média por bloco de {args.ns / 8:.0f} ns")
        larguras = " ".join(f"{(i + 1) * args.ns / 8:6.0f}" for i in range(8))
        print(f"{'réplica':>8}  {larguras}   {'tendência':>22}")
        for rep in sorted(series):
            v = series[rep].get(chave)
            if not v:
                continue
            b = blocos(v)
            linha = " ".join(f"{x:6.2f}" for x in b)
            # a tendência é a diferença entre o último quarto e o segundo
            inicio, fim = np.mean(b[1:3]), np.mean(b[-2:])
            delta = fim - inicio
            if abs(delta) < 0.1 * max(abs(inicio), 1e-9):
                tend = "PLATÔ"
            else:
                tend = f"{'sobe' if delta > 0 else 'desce'} {abs(delta):.2f}"
            print(f"{rep:>8}  {linha}   {tend:>22}")

    if args.metade:
        print(f"\n{'=' * 66}\nSó a SEGUNDA METADE ({args.ns / 2:.0f}–{args.ns:.0f} ns)"
              f"\n{'=' * 66}")
        for chave, rotulo in (("rmsd_prot", "RMSD da proteína"),
                              ("rmsd_lig", "RMSD do PROTAC")):
            vals = []
            for rep in sorted(series):
                v = series[rep].get(chave)
                if not v:
                    continue
                meio = len(v) // 2
                m = float(np.mean(v[meio:]))
                vals.append(m)
                print(f"  {rep}: {rotulo} = {m:.2f} Å")
            if vals:
                print(f"  -> média das réplicas: {np.mean(vals):.2f} Å\n")

    print("\nComo ler: PLATÔ significa que a média sobre os 200 ns inteiros")
    print("está inflada pelo transiente inicial e não descreve o estado")
    print("amostrado. 'sobe' persistente significa que o sistema não")
    print("convergiu, e aí a média está certa em reprovar.")


if __name__ == "__main__":
    main()
