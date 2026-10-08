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
    # Com 105 clusters de um membro cada, comparar todos é comparar ruído — e
    # imprimir 105x105 é ilegível. Entram os maiores, que são os que importam.
    ap.add_argument("--limiar-clustering", type=float, default=4.0,
                    help="o limiar que o clustering.py do PRosettaC usou (o "
                         "3o argumento dele; main.py passa 4)")
    ap.add_argument("--max-clusters", type=int, default=30,
                    help="quantos clusters comparar, os de mais membros antes")
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

    total_clusters, total_modelos = len(modelos), sum(membros.values())
    # os maiores primeiro: um cluster de 8 membros diz mais que oito de um
    nomes = sorted(modelos, key=lambda n: (-membros[n], n))
    cortados = 0
    if len(nomes) > args.max_clusters:
        cortados = len(nomes) - args.max_clusters
        nomes = nomes[:args.max_clusters]
    print(f"  {total_clusters} clusters, {total_modelos} modelos no total")
    print(f"  membros dos maiores: "
          f"{', '.join(str(membros[n]) for n in nomes[:12])}"
          + (" ..." if len(nomes) > 12 else ""))
    com_massa = [n for n in modelos if membros[n] >= 5]
    if com_massa:
        print(f"  cluster(s) com 5+ membros: "
              + ", ".join(f"{n} ({membros[n]})" for n in com_massa))
    if cortados:
        print(f"  comparando os {len(nomes)} maiores "
              f"({cortados} de 1 membro ficaram de fora — use "
              f"--max-clusters para incluir)")
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

    if len(nomes) <= 24:
        larg = max(len(n) for n in nomes) + 2
        print(" " * (larg + 3) + " ".join(f"{i:>6d}"
                                          for i in range(len(nomes))))
        for i, n in enumerate(nomes):
            print(f"{i:>2d} {n:<{larg}}" + " ".join(
                "     ." if i == j else
                ("   n/a" if np.isnan(M[i, j]) else f"{M[i, j]:6.1f}")
                for j in range(len(nomes))))
    else:
        print(f"  ({len(nomes)}x{len(nomes)} não cabe na tela — está no CSV"
              f" do fim)")

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

    # O maior grupo MUTUAMENTE próximo — o que seria "o cluster" se o limiar
    # do PRosettaC fosse o desta medida. Isto é clique máxima, e a primeira
    # versão daqui fazia busca gulosa a partir de cada vértice: com 20 modelos
    # a busca exaustiva custa milissegundos e não erra, então não há motivo
    # para aproximar.
    def maior_grupo(corte: float):
        A = (M <= corte) & ~np.eye(len(nomes), dtype=bool) & ~np.isnan(M)
        melhor: list[int] = []

        # Orçamento de chamadas: a busca exata é exponencial, e com 30 nós num
        # grafo denso ela não termina. Estourando, devolve o melhor encontrado
        # até ali e AVISA — um número silenciosamente aproximado num veredito
        # seria pior que um número menor declarado.
        orcamento = [400_000]
        estourou = [False]

        def expandir(atual, cands):
            nonlocal melhor
            if orcamento[0] <= 0:
                estourou[0] = True
                return
            orcamento[0] -= 1
            if len(atual) > len(melhor):
                melhor = list(atual)
            for i, v in enumerate(cands):
                if len(atual) + len(cands) - i <= len(melhor):
                    return
                expandir(atual + [v], [u for u in cands[i + 1:] if A[v, u]])

        expandir([], list(range(len(nomes))))
        return melhor, estourou[0]

    # Um corte só responde à pergunta errada: a 5 Å tudo parece disperso, a
    # 20 Å tudo parece junto. A varredura mostra a FORMA da concordância.
    print(f"\n  maior grupo mutuamente próximo, por corte:")
    grupos = {}
    for corte in (args.corte, args.corte * 2, args.corte * 3, args.corte * 4):
        g, aprox = maior_grupo(corte)
        grupos[corte] = g
        print(f"    {corte:>5.0f} Å   {len(g):>2} modelo(s)"
              + ("  (mínimo: busca truncada)" if aprox else "")
              + (f"   {', '.join(nomes[k] for k in g)}" if 1 < len(g) <= 8
                 else ""))
    maior = len(grupos[args.corte])
    quem = grupos[args.corte]
    no_dobro = len(grupos[args.corte * 2])

    print()
    # ------------------------------------------------------------------
    # O CORTE DESTA MEDIDA NÃO PODE SER COMPARÁVEL AO DO CLUSTERING.
    #
    # Esta é a correção de um erro de projeto que já produziu um veredito
    # falso. O clustering.py do PRosettaC roda DBSCAN com eps = 4 Å sobre a
    # cadeia móvel — o MESMO deslocamento que esta análise mede. Dois clusters
    # distintos do DBSCAN, por construção, não são densamente conectados a
    # 4 Å: os representantes deles estão a MAIS de 4 Å um do outro.
    #
    # Exigir 5 representantes mutuamente a menos de 5 Å é, então, pedir que o
    # DBSCAN tenha produzido clusters que ele mesmo teria fundido. É quase
    # logicamente impossível, e o "NÃO CONVERGIU" que sai daí é artefato do
    # critério, não achado sobre a geometria.
    #
    # Isto vale a priori, lendo os dois códigos, e não porque o portão
    # reprovou. Um portão consertado DEPOIS de falhar só é honesto quando o
    # defeito se demonstra sem olhar o resultado — e este se demonstra.
    # ------------------------------------------------------------------
    piso_valido = args.limiar_clustering * 1.5
    if args.corte <= piso_valido:
        print(f"  [CORTE INCOMPARÁVEL] este corte ({args.corte:.1f} Å) está no"
              f" mesmo patamar do")
        print(f"  limiar do clustering do PRosettaC ({args.limiar_clustering:.1f} Å,"
              f" DBSCAN sobre a mesma")
        print(f"  medida). Representantes de clusters distintos estão a mais de"
              f" {args.limiar_clustering:.1f} Å")
        print(f"  POR CONSTRUÇÃO, então um grupo grande aqui seria contradição"
              f" com o")
        print(f"  próprio agrupamento. Nenhum veredito de convergência sai"
              f" desta medida")
        print(f"  com este corte — ela fica como DESCRIÇÃO da dispersão, acima.")
        print(f"  Para um veredito desta medida, use --corte >"
              f" {piso_valido:.1f}.")
    elif maior <= 1 and no_dobro <= 3:
        print(f"  -> NÃO CONVERGIU. Nenhum par a menos de {args.corte:.0f} Å, e mesmo")
        print(f"     dobrando o corte o maior grupo tem {no_dobro} de {len(nomes)} modelos.")
    elif maior >= 5 or maior == len(nomes):
        print(f"  -> Há um grupo de {maior} modelos concordando a menos de")
        print(f"     {args.corte:.0f} Å: a região é uma só, e é candidata a pose.")
    else:
        print(f"  -> Concordância PARCIAL: o maior grupo tem {maior} de {len(nomes)}.")

    # ------------------------------------------------------------------
    # A medida que de fato responde "há uma pose dominante?" é a CONCENTRAÇÃO
    # DA POPULAÇÃO, e não a coincidência de representantes. O critério usado
    # aqui é o do PRÓPRIO PRosettaC, que imprime no result_summary.txt
    # "Out of them N have at least 5 members" — ou seja, 5 membros é a régua
    # que os autores da ferramenta escolheram. Adotá-la não é mover o alvo:
    # ela existia antes deste resultado e não foi calibrada por ele.
    # ------------------------------------------------------------------
    mais_membros = max(membros.values())
    com_massa_n = sum(v for v in membros.values() if v >= 5)
    singletons = sum(1 for v in membros.values() if v == 1)
    print()
    print("  CONCENTRAÇÃO DA POPULAÇÃO (critério do próprio PRosettaC: 5+ membros)")
    print(f"    maior cluster ............ {mais_membros} de {total_modelos} "
          f"modelos ({100 * mais_membros / total_modelos:.0f}%)")
    print(f"    em clusters de 5+ membros  {com_massa_n} de {total_modelos} "
          f"({100 * com_massa_n / total_modelos:.0f}%)")
    print(f"    clusters de 1 membro ..... {singletons} de {total_clusters}")
    if mais_membros >= 5:
        print(f"    -> existe pose dominante a testar: o maior cluster. Ela é")
        print(f"       CANDIDATA, não confirmada — {100 - 100 * mais_membros / total_modelos:.0f}%"
              f" dos modelos estão em")
        print(f"       outras {total_clusters - 1} poses, e os representantes delas estão a")
        print(f"       {np.median(fora):.0f} Å de mediana. Quem decide é o método ORTOGONAL")
        print(f"       (Boltz-2/AlphaFold 3), como a metodologia do WP3 pede:")
        print(f"       se ele cair sobre esta geometria, a pose está confirmada.")
    else:
        print(f"    -> NENHUMA pose dominante: o maior cluster tem"
              f" {mais_membros} membros,")
        print(f"       abaixo da régua de 5 da própria ferramenta. Não há o que")
        print(f"       levar ao método ortogonal nem à MD.")

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
