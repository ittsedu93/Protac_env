#!/usr/bin/env python
"""
Estruturas finais da MD: o PDB representativo e os quadros do filme.

Roda no env `mdtools` (MDAnalysis). Minutos.

    python md_export_structures.py --md-dir ~/pipeline/md

Escreve em <md-dir>/structures/:

    final_complex_<rep>.pdb     o quadro REPRESENTATIVO da segunda metade
    movie_<rep>.pdb             multi-modelo, para o filme
    representative_frames.csv   qual quadro foi escolhido e por quê

Por que o quadro representativo e não o último
----------------------------------------------
O último quadro é um instante arbitrário — pode ter pegado uma flutuação. O
quadro representativo é o que está mais perto da estrutura MÉDIA da segunda
metade da trajetória, que é o estado efetivamente amostrado depois do
equilíbrio. É ele que se mostra como "a estrutura do complexo".

O alinhamento é pelo SÍTIO do recrutador, não pela proteína inteira: assim o
bolso fica parado na tela e o que se vê mexer é o movimento de domínio,
relativo a ele. Alinhar pela proteína inteira repartiria o giro entre tudo.

Nomenclatura: este é o complexo BINÁRIO E3-PROTAC (CRBN + recrutador-linker-
warhead). NÃO é ternário — a PCSK9 não está no sistema.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def sitio(u, corte=12.0):
    """Seleção dos CA perto do ligante no primeiro quadro."""
    from MDAnalysis.analysis.distances import distance_array
    u.trajectory[0]
    lig = u.select_atoms("not protein")
    ca = u.select_atoms("protein and name CA")
    if len(lig) == 0 or len(ca) == 0:
        return "protein and name CA"
    d = distance_array(ca.positions, lig.positions).min(axis=1)
    resids = [int(a.resid) for a, dd in zip(ca, d) if dd <= corte]
    if len(resids) < 10:
        return "protein and name CA"
    return "protein and name CA and resid " + " ".join(str(r) for r in resids)


def exportar(md: Path, rep: str, out: Path, n_quadros: int, ns: float):
    import MDAnalysis as mda
    from MDAnalysis.analysis import align, rms

    gro, xtc = md / rep / "solutos.gro", md / rep / "solutos.xtc"
    if not (gro.exists() and xtc.exists()):
        print(f"  {rep}: sem solutos.* — rode md_fix_pbc.sh antes")
        return None

    u = mda.Universe(str(gro), str(xtc))

    # Sem caixa, o transfer_to_memory do AlignTraj morre com um TypeError
    # sobre dtype que não diz nada sobre a causa. Uma caixa grande o bastante
    # para não envolver nada é inócua aqui: a trajetória já passou pelo
    # md_fix_pbc.sh e a análise não usa imagens periódicas.
    if u.dimensions is None or not np.all(np.isfinite(np.asarray(
            u.dimensions, dtype=float))):
        print(f"  {rep}: sem caixa no .gro — usando uma caixa nominal "
              f"(a trajetória já está tratada, então isto é inócuo)")
        ext = u.atoms.positions.ptp(axis=0).max() * 3 + 100.0
        for ts in u.trajectory:
            ts.dimensions = [ext, ext, ext, 90.0, 90.0, 90.0]

    sel = sitio(u)

    # alinha tudo pelo sítio, em memória
    media = align.AverageStructure(u, u, select=sel, ref_frame=0).run()
    align.AlignTraj(u, media.results.universe, select=sel,
                    in_memory=True).run()

    n = len(u.trajectory)
    meia = n // 2

    # estrutura média da segunda metade, e o quadro mais perto dela
    todos = u.select_atoms("all")
    soma = np.zeros_like(todos.positions)
    for ts in u.trajectory[meia:]:
        soma += todos.positions
    ref = soma / (n - meia)

    melhor, melhor_d = None, None
    for i, ts in enumerate(u.trajectory[meia:], start=meia):
        d = float(np.sqrt(np.mean(np.sum((todos.positions - ref) ** 2, axis=1))))
        if melhor_d is None or d < melhor_d:
            melhor, melhor_d = i, d

    u.trajectory[melhor]
    alvo = out / f"final_complex_{rep}.pdb"
    todos.write(str(alvo))
    print(f"  {rep}: quadro {melhor}/{n} (a {melhor_d:.2f} Å da média) "
          f"-> {alvo.name}")

    # multi-modelo para o filme, subamostrado
    passo = max(1, n // n_quadros)
    filme = out / f"movie_{rep}.pdb"
    with mda.Writer(str(filme), todos.n_atoms, multiframe=True) as w:
        for ts in u.trajectory[::passo]:
            w.write(todos)
    print(f"  {rep}: {len(range(0, n, passo))} quadros -> {filme.name}")

    return {"replicate": rep, "n_frames": n, "representative_frame": melhor,
            "time_ns": round(melhor / max(n - 1, 1) * ns, 1),
            "rmsd_to_mean_A": round(melhor_d, 3),
            "movie_frames": len(range(0, n, passo))}


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--md-dir", type=Path, required=True)
    ap.add_argument("--quadros", type=int, default=200,
                    help="quadros no PDB multi-modelo do filme")
    ap.add_argument("--ns", type=float, default=200.0)
    ap.add_argument("--replica", default=None,
                    help="só uma réplica (default: todas)")
    args = ap.parse_args()

    md = args.md_dir.expanduser()
    out = md / "structures"
    out.mkdir(exist_ok=True)

    reps = [args.replica] if args.replica else \
        sorted(d.name for d in md.glob("rep*") if d.is_dir())

    print("Estruturas do complexo binário E3-PROTAC:")
    linhas = [r for r in (exportar(md, rep, out, args.quadros, args.ns)
                          for rep in reps) if r]
    if linhas:
        import pandas as pd
        pd.DataFrame(linhas).to_csv(out / "representative_frames.csv",
                                    index=False)
    print(f"\n{out}")


if __name__ == "__main__":
    main()
