#!/usr/bin/env python
"""
Publication-quality figures for the level-(ii) MD, in English.

Runs in the `mdtools` env (MDAnalysis + matplotlib). Minutes.

    python md_figures.py --md-dir ~/pipeline/md

Writes to <md-dir>/figures/:

    fig1_rmsd.png/.pdf        RMSD vs time: site, core, whole protein
    fig2_rmsf.png/.pdf        per-residue RMSF, site region shaded
    fig3_contacts.png/.pdf    recruiter and warhead contacts vs time
    fig4_protac.png/.pdf      PROTAC RMSD vs time
    *.csv                     source data for every panel

Every figure ships its numbers as CSV next to it: a figure whose data cannot be
re-plotted by someone else is a picture, not a result.

Note on naming: this is the **binary E3-PROTAC complex** (CRBN + recruiter-
linker-warhead), NOT a ternary complex — PCSK9 is not in the system. Calling it
ternary in a manuscript would be wrong.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

# --- design tokens --------------------------------------------------------
# Categorical slots in fixed order, never cycled (validated for CVD against the
# chart surface). Replicate identity always maps to the same hue.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#8d8c86"
GRID = "#e3e2dd"
SURFACE = "#ffffff"
CRIT = "#b0403a"          # criterion line — status colour, never a series


def estilo():
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 300,
        "savefig.bbox": "tight", "savefig.facecolor": SURFACE,
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "font.family": "sans-serif",
        "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 10,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
        "axes.edgecolor": MUTED, "axes.linewidth": 0.8,
        "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": INK_2, "ytick.color": INK_2,
        "xtick.direction": "out", "ytick.direction": "out",
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
        "axes.axisbelow": True, "legend.frameon": False,
        "lines.linewidth": 1.4, "lines.solid_capstyle": "round",
    })
    return plt


def limpar(ax, grid_y=True):
    """Recessive axes: no top/right spines, grid on one axis only."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y" if grid_y else "x")
    ax.grid(axis="x" if grid_y else "y", visible=False)


def rotular_fim(ax, x, y, texto, cor):
    """Direct label at the end of a line — identity never by colour alone."""
    ax.annotate(texto, xy=(x[-1], y[-1]), xytext=(4, 0),
                textcoords="offset points", color=cor, fontsize=8,
                va="center", ha="left")


def salvar(fig, out: Path, nome: str):
    for ext in ("png", "pdf"):
        fig.savefig(out / f"{nome}.{ext}")
    print(f"      {nome}.png / .pdf")


# ==========================================================================
def fig_rmsd(plt, series, out: Path, ns: float):
    """Three panels, shared y: the message is that the site is flat and low."""
    import pandas as pd

    chaves = [("rmsd_sitio", "Recruiter site\n(CA within 12 Å)"),
              ("rmsd_nucleo", "Structured core\n(termini excluded)"),
              ("rmsd_prot", "Whole protein\n(all CA)")]
    chaves = [(k, r) for k, r in chaves
              if any(k in s for s in series.values())]
    if not chaves:
        return

    fig, axes = plt.subplots(1, len(chaves), figsize=(3.1 * len(chaves), 2.9),
                             sharey=True)
    axes = np.atleast_1d(axes)
    linhas = []
    for ax, (chave, rotulo) in zip(axes, chaves):
        for i, rep in enumerate(sorted(series)):
            v = np.asarray(series[rep].get(chave, []), dtype=float)
            if v.size == 0:
                continue
            t = np.linspace(0, ns, v.size)
            ax.plot(t, v, color=SERIES[i % len(SERIES)],
                    label=f"replicate {i + 1}")
            linhas += [{"metric": chave, "replicate": i + 1,
                        "time_ns": round(float(a), 3), "rmsd_A": round(float(b), 4)}
                       for a, b in zip(t, v)]
        ax.axhline(3.5, color=CRIT, lw=1.0, ls=(0, (4, 3)), zorder=1)
        ax.set_title(rotulo, color=INK, pad=6)
        ax.set_xlabel("Time (ns)")
        ax.set_xlim(0, ns)
        limpar(ax)
    axes[0].set_ylabel("RMSD (Å)")
    axes[0].set_ylim(bottom=0)
    axes[-1].annotate("acceptance criterion 3.5 Å", xy=(ns * 0.98, 3.5),
                      xytext=(0, 4), textcoords="offset points",
                      color=CRIT, fontsize=7.5, ha="right")
    alças, rots = axes[0].get_legend_handles_labels()
    fig.legend(alças, rots, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.12))
    fig.suptitle("Backbone RMSD of the E3 ligase during 200 ns "
                 "(three independent replicates)",
                 fontsize=10.5, color=INK, y=1.16, x=0.02, ha="left")
    salvar(fig, out, "fig1_rmsd")
    pd.DataFrame(linhas).to_csv(out / "fig1_rmsd.csv", index=False)
    plt.close(fig)


def fig_contatos(plt, series, out: Path, ns: float):
    import pandas as pd
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.9), sharey=True)
    linhas = []
    for ax, (chave, rotulo) in zip(
            axes, [("c_rec", "Recruiter – E3"),
                   ("c_wh", "Warhead – E3")]):
        for i, rep in enumerate(sorted(series)):
            v = np.asarray(series[rep].get(chave, []), dtype=float)
            if v.size == 0:
                continue
            t = np.linspace(0, ns, v.size)
            # a light rolling mean on top of the raw trace: the raw data is
            # noisy per frame and the trend is the claim being made
            # `mode="same"` divide as pontas por uma janela cheia de zeros e
            # cria quedas que não existem nos dados. Estender pelas bordas
            # antes de convoluir mantém as extremidades honestas.
            jan = max(1, v.size // 50)
            pad = jan // 2
            suave = np.convolve(np.pad(v, pad, mode="edge"),
                                np.ones(jan) / jan, mode="same")[pad:pad + v.size]
            ax.plot(t, v, color=SERIES[i % len(SERIES)], lw=0.5, alpha=0.30)
            ax.plot(t, suave, color=SERIES[i % len(SERIES)],
                    label=f"replicate {i + 1}")
            linhas += [{"contact_set": chave, "replicate": i + 1,
                        "time_ns": round(float(a), 3), "contacts": int(b)}
                       for a, b in zip(t, v)]
        ax.set_title(rotulo, color=INK, pad=6)
        ax.set_xlabel("Time (ns)")
        ax.set_xlim(0, ns)
        limpar(ax)
    axes[0].set_ylabel("Heavy-atom contacts (< 4.5 Å)")
    axes[0].set_ylim(bottom=0)
    alças, rots = axes[0].get_legend_handles_labels()
    fig.legend(alças, rots, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.10))
    fig.suptitle("Contacts with the E3 ligase are maintained throughout",
                 fontsize=10.5, color=INK, y=1.06, x=0.02, ha="left")
    salvar(fig, out, "fig3_contacts")
    pd.DataFrame(linhas).to_csv(out / "fig3_contacts.csv", index=False)
    plt.close(fig)


def fig_protac(plt, series, out: Path, ns: float):
    import pandas as pd
    fig, ax = plt.subplots(figsize=(4.2, 2.9))
    linhas = []
    for i, rep in enumerate(sorted(series)):
        v = np.asarray(series[rep].get("rmsd_lig", []), dtype=float)
        if v.size == 0:
            continue
        t = np.linspace(0, ns, v.size)
        ax.plot(t, v, color=SERIES[i % len(SERIES)])
        rotular_fim(ax, t, v, f"rep {i + 1}", SERIES[i % len(SERIES)])
        linhas += [{"replicate": i + 1, "time_ns": round(float(a), 3),
                    "rmsd_A": round(float(b), 4)} for a, b in zip(t, v)]
    ax.axhline(5.0, color=CRIT, lw=1.0, ls=(0, (4, 3)), zorder=1)
    ax.annotate("acceptance criterion 5.0 Å", xy=(ns * 0.98, 5.0),
                xytext=(0, 4), textcoords="offset points",
                color=CRIT, fontsize=7.5, ha="right")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("RMSD (Å)")
    ax.set_xlim(0, ns * 1.10)
    ax.set_ylim(bottom=0)
    limpar(ax)
    ax.set_title("PROTAC heavy-atom RMSD", color=INK, pad=6)
    salvar(fig, out, "fig4_protac")
    pd.DataFrame(linhas).to_csv(out / "fig4_protac.csv", index=False)
    plt.close(fig)


def fig_rmsf(plt, dados, out: Path):
    """Per-residue RMSF. Localises the motion — the hinge claim lives here."""
    import pandas as pd
    if not dados:
        return
    fig, ax = plt.subplots(figsize=(7.2, 2.9))
    linhas = []
    for i, (rep, d) in enumerate(sorted(dados.items())):
        ax.plot(d["resids"], d["rmsf"], color=SERIES[i % len(SERIES)],
                label=f"replicate {i + 1}", lw=1.1)
        linhas += [{"replicate": i + 1, "resid": int(r),
                    "rmsf_A": round(float(v), 4)}
                   for r, v in zip(d["resids"], d["rmsf"])]
    prim = dados[sorted(dados)[0]]
    # site residues shaded: the reader must see WHERE the motion is not
    for ini, fim in prim.get("faixas_sitio", []):
        ax.axvspan(ini, fim, color="#2a78d6", alpha=0.10, lw=0, zorder=0)
    # Quatro tracejados marcam os cortes; quarenta viram hachura e escondem
    # os dados. Acima de um punhado, a informação vira ruído e sai.
    cortes = prim.get("cortes", [])
    if 0 < len(cortes) <= 10:
        for corte in cortes:
            ax.axvline(corte, color=MUTED, lw=0.8, ls=(0, (2, 3)), zorder=0)
    elif len(cortes) > 10:
        print(f"      {len(cortes)} cortes detectados — não marcados na figura "
              f"(viraria hachura); estão no CSV")
    ax.set_xlabel("Residue number")
    ax.set_ylabel("RMSF (Å)")
    ax.set_ylim(bottom=0)
    limpar(ax)
    alças, rots = ax.get_legend_handles_labels()
    fig.legend(alças, rots, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.14))
    sub = ("; chain breaks dashed" if 0 < len(cortes) <= 10 else "")
    ax.set_title(f"Per-residue fluctuation: the recruiter site (shaded) "
                 f"is rigid{sub}", color=INK, pad=6, loc="left")
    salvar(fig, out, "fig2_rmsf")
    pd.DataFrame(linhas).to_csv(out / "fig2_rmsf.csv", index=False)
    plt.close(fig)


# ==========================================================================
def calcular_rmsf(md: Path, reps):
    """RMSF por resíduo, alinhando a trajetória pelo sítio.

    Alinhar pelo SÍTIO e não pela proteína inteira é o ponto: alinhado pelo
    sítio, o que aparece como flutuação é o movimento RELATIVO ao bolso do
    recrutador — que é a pergunta. Alinhado pela proteína inteira, o giro de
    domínio se reparte entre todos os resíduos e não se localiza nada.
    """
    import MDAnalysis as mda
    from MDAnalysis.analysis import align, rms

    dados = {}
    for rep in reps:
        gro, xtc = md / rep / "solutos.gro", md / rep / "solutos.xtc"
        if not (gro.exists() and xtc.exists()):
            print(f"      {rep}: sem solutos.* — rode md_fix_pbc.sh antes")
            continue
        u = mda.Universe(str(gro), str(xtc))
        lig = u.select_atoms("not protein")
        u.trajectory[0]
        from MDAnalysis.analysis.distances import distance_array
        ca = u.select_atoms("protein and name CA")
        d = distance_array(ca.positions, lig.positions).min(axis=1)
        resids_sitio = [int(a.resid) for a, dd in zip(ca, d) if dd <= 12.0]
        sel = ("protein and name CA and resid "
               + " ".join(str(r) for r in resids_sitio)) if resids_sitio \
            else "protein and name CA"

        media = align.AverageStructure(u, u, select=sel, ref_frame=0).run()
        align.AlignTraj(u, media.results.universe, select=sel,
                        in_memory=True).run()
        ca = u.select_atoms("protein and name CA")
        R = rms.RMSF(ca).run()

        resids = np.array([int(a.resid) for a in ca])
        cortes = []
        u.trajectory[0]
        res = u.select_atoms("protein").residues
        for i in range(len(res) - 1):
            c = res[i].atoms.select_atoms("name C")
            n = res[i + 1].atoms.select_atoms("name N")
            if len(c) and len(n) and \
                    np.linalg.norm(c.positions[0] - n.positions[0]) > 2.5:
                cortes.append(int(res[i].resid))

        faixas, atual = [], None
        s = set(resids_sitio)
        for r in resids:
            if r in s and atual is None:
                atual = r
            elif r not in s and atual is not None:
                faixas.append((atual, r)); atual = None
        if atual is not None:
            faixas.append((atual, resids[-1]))

        dados[rep] = {"resids": resids, "rmsf": R.results.rmsf,
                      "faixas_sitio": faixas, "cortes": cortes}
        print(f"      {rep}: RMSF de {len(resids)} resíduos "
              f"(máx {R.results.rmsf.max():.2f} Å)")
    return dados


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--md-dir", type=Path, required=True)
    ap.add_argument("--ns", type=float, default=200.0)
    ap.add_argument("--sem-rmsf", action="store_true",
                    help="pula o RMSF (que precisa reler as trajetórias)")
    args = ap.parse_args()

    md = args.md_dir.expanduser()
    out = md / "figures"
    out.mkdir(exist_ok=True)
    series = json.loads((md / "md_series.json").read_text())
    plt = estilo()

    print("Figures (English, 300 dpi, PNG + PDF + source CSV):")
    fig_rmsd(plt, series, out, args.ns)
    fig_contatos(plt, series, out, args.ns)
    fig_protac(plt, series, out, args.ns)

    if not args.sem_rmsf:
        print("  RMSF (relê as trajetórias, alguns minutos)")
        dados = calcular_rmsf(md, sorted(series))
        fig_rmsf(plt, dados, out)

    print(f"\n{out}")


if __name__ == "__main__":
    main()
