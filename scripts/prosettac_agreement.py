#!/usr/bin/env python
"""
"20 clusters" é vinte lugares diferentes, ou uma região amostrada grosso?

Roda no env `mdtools`. Segundos.

    python prosettac_agreement.py --candidato SC0013__WH022

O que isto resolve
------------------
O PRosettaC fechou assim o `SC0013__WH022`:

    658 local docking solutions were generated.
    21 final models with energy below the threshold (0) were generated.
    20 clusters were generated.
    Out of them 0 have at least 5 members.

Pelo critério da própria ferramenta isso é negativo: nenhum cluster com 5
membros. Mas "20 clusters de 21 modelos" admite duas leituras opostas, e elas
levam a decisões diferentes:

  - 20 interfaces INCOMPATÍVEIS      -> o método não encontrou geometria
  - UMA região, amostrada mais grosso
    que o limiar de 2 Å do clustering -> há sinal, e o que falta é amostragem

A distinção é medível do mesmo jeito que foi para o Boltz-2: superpõe-se cada
modelo pela E3 — que é a âncora experimental, vinda do cristal — e mede-se o
quanto o ALVO se desloca entre modelos. Dez modelos a 3 Å um do outro são uma
interface; dez modelos a 40 Å são dez hipóteses.

O pareamento de átomos é por (cadeia, número de resíduo), não por ordem no
arquivo: modelos do Rosetta podem perder resíduo, e comparar por posição
alinharia átomos de resíduos diferentes sem reclamar.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from ternary_agreement import kabsch


def ler_ca_por_residuo(pdb: Path):
    """{(cadeia, num_residuo): (x, y, z)} dos CA."""
    d = {}
    for l in pdb.read_text().splitlines():
        if not l.startswith("ATOM") or l[12:16].strip() != "CA":
            continue
        if l[16] not in (" ", "A"):          # só a altloc principal
            continue
        try:
            chave = (l[21], int(l[22:26]))
        except ValueError:
            continue
        d[chave] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    return d


def cadeias_do_config(cfg: Path):
    """(cadeia da E3, cadeia do alvo), na ordem em que o config as declara."""
    for linha in cfg.read_text().splitlines():
        if linha.startswith("Chains: "):
            partes = linha.split(": ", 1)[1].split()
            if len(partes) >= 2:
                return partes[0], partes[1]
    raise SystemExit(f"não achei a linha Chains: em {cfg}")


def comuns(a: dict, b: dict, cadeia: str):
    """Coordenadas pareadas por resíduo, só da cadeia pedida."""
    ch = sorted(k for k in a if k[0] == cadeia and k in b)
    if not ch:
        return None, None
    return (np.array([a[k] for k in ch]), np.array([b[k] for k in ch]))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidato", required=True)
    ap.add_argument("--work", type=Path,
                    default=Path.home() / "PRosettaC_runs"
                    / "vhl_crbn_pcsk9_protac" / "pipeline" / "wp3"
                    / "prosettac")
    ap.add_argument("--corte", type=float, default=5.0,
                    help="Å abaixo do qual dois modelos concordam")
    args = ap.parse_args()

    base = args.work.expanduser() / args.candidato
    res = base / "Results"
    if not res.is_dir():
        raise SystemExit(f"não achei {res} — o PRosettaC terminou?")

    cad_e3, cad_alvo = cadeias_do_config(base / "prosetta_config.txt")
    print(f"candidato: {args.candidato}")
    print(f"  E3: cadeia {cad_e3} | alvo: cadeia {cad_alvo}")

    # um representante por cluster, e quantos membros cada um tem
    clusters = sorted(
        [d for d in res.iterdir() if d.is_dir()],
        key=lambda d: int("".join(c for c in d.name if c.isdigit()) or 0))
    modelos, membros = {}, {}
    for d in clusters:
        pdbs = sorted(d.rglob("*.pdb"))
        if not pdbs:
            continue
        modelos[d.name] = ler_ca_por_residuo(pdbs[0])
        membros[d.name] = len(pdbs)
    if len(modelos) < 2:
        raise SystemExit(f"achei {len(modelos)} modelo(s) — preciso de 2")

    nomes = list(modelos)
    print(f"  {len(nomes)} clusters, {sum(membros.values())} modelos no total")
    print(f"  membros por cluster: "
          f"{', '.join(str(membros[n]) for n in nomes)}")
    print()

    M = np.full((len(nomes), len(nomes)), np.nan)
    for i, a in enumerate(nomes):
        for j, b in enumerate(nomes):
            if i >= j:
                continue
            ea, eb = comuns(modelos[a], modelos[b], cad_e3)
            ta, tb = comuns(modelos[a], modelos[b], cad_alvo)
            if ea is None or ta is None:
                continue
            R, cm, cr = kabsch(ea, eb)
            ta_al = (ta - cm) @ R.T + cr
            M[i, j] = M[j, i] = float(np.sqrt(np.mean(
                np.sum((ta_al - tb) ** 2, axis=1))))

    larg = max(len(n) for n in nomes) + 2
    print(" " * (larg + 3) + " ".join(f"{i:>6d}" for i in range(len(nomes))))
    for i, n in enumerate(nomes):
        print(f"{i:>2d} {n:<{larg}}" + " ".join(
            "     ." if i == j else
            ("   n/a" if np.isnan(M[i, j]) else f"{M[i, j]:6.1f}")
            for j in range(len(nomes))))

    fora = M[~np.eye(len(nomes), dtype=bool)]
    fora = fora[~np.isnan(fora)]
    if not len(fora):
        raise SystemExit("\n  nenhum par comparável — as cadeias não casaram")

    total = len(nomes) * (len(nomes) - 1) // 2
    proximos = int((fora <= args.corte).sum() / 2)
    print(f"\n  RMSD do alvo após superpor pela E3 (Å):")
    print(f"    menor {fora.min():.1f} | mediana {np.median(fora):.1f} "
          f"| maior {fora.max():.1f}")
    print(f"    {proximos} de {total} pares a menos de {args.corte:.0f} Å")

    # o maior grupo mutuamente próximo: é ele que seria "o cluster" se o
    # limiar do PRosettaC fosse o desta medida
    perto = (M <= args.corte) & ~np.isnan(M)
    maior, quem = 0, []
    for i in range(len(nomes)):
        grupo = [i] + [j for j in range(len(nomes)) if j != i and perto[i, j]]
        # mutuamente próximos, não só próximos do centro
        grupo = [k for k in grupo
                 if all(k == m or perto[k, m] for m in grupo)]
        if len(grupo) > maior:
            maior, quem = len(grupo), grupo
    print(f"    maior grupo mutuamente a menos de {args.corte:.0f} Å: "
          f"{maior} modelo(s)")
    if maior > 1:
        print(f"      {', '.join(nomes[k] for k in quem)}")

    print()
    if fora.min() > 10:
        print("  -> NENHUM par concorda. Os 20 clusters são 20 interfaces")
        print("     incompatíveis, não uma região amostrada grosso. O")
        print("     PRosettaC não encontrou a geometria ternária deste par —")
        print("     e isso confirma, por um método independente, o que o")
        print("     Boltz-2 já havia indicado.")
    elif maior >= 5 or maior == len(nomes):
        # `maior == len(nomes)` cobre o caso de poucos clusters: com 3 modelos
        # todos concordando, exigir 5 daria "parcial" para uma concordância
        # completa.
        print(f"  -> Há um grupo de {maior} modelos concordando a menos de")
        print(f"     {args.corte:.0f} Å. O clustering do PRosettaC usa limiar mais")
        print("     apertado e por isso os separou — mas a região é uma só, e")
        print("     ela é candidata a pose para a MD nível (iii).")
    else:
        print(f"  -> Concordância PARCIAL: o maior grupo tem {maior} modelos.")
        print("     Não sustenta uma pose única, e sustenta menos ainda um")
        print("     veredito de ausência. O que falta é amostragem — e a causa")
        print("     provável está no linker no limite do alcance: modelo que só")
        print("     fecha a ponte estendido paga energia, e foi a energia que")
        print("     cortou 658 soluções para 21.")

    saida = base / "concordancia_clusters.csv"
    with open(saida, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([""] + nomes)
        for i, n in enumerate(nomes):
            w.writerow([n] + ["" if np.isnan(M[i, j]) else f"{M[i, j]:.2f}"
                              for j in range(len(nomes))])
    print(f"\n  {saida}")


if __name__ == "__main__":
    main()
