#!/usr/bin/env python
"""
Fase 2.5 — o portão que faltava: o linker ALCANÇA o vão entre os dois sítios?

Roda no env `mdtools` (RDKit). Minutos para o catálogo inteiro.

    python linker_span.py \\
        --protacs ~/pipeline/wp3/protac_candidates.csv \\
        --e3-head ~/pipeline/wp1/recruiter_CRBN_ligand_116.sdf \\
        --warheads-sdf ~/PCSK9_docking/poses \\
        --requisito 18 \\
        --out ~/pipeline/wp3/linker_span.csv

O que isto resolve, e por que ele vem ANTES da MD
-------------------------------------------------
O PRosettaC restringe o docking global a uma distância máxima entre os dois
átomos de conjugação — o do recrutador, dentro do bolsão da CRBN, e o da
warhead, dentro do sítio da PCSK9. Essa distância é o ALCANCE do linker, e ela
é medida amostrando confôrmeros do PROTAC inteiro.

No `SC0006__WH023` o alcance foi 11 Å, e a varredura do PatchDock
(`patchdock_span_scan.sh`) mostrou que o par CRBN/PCSK9 só produz soluções a
partir de 14 Å, com regime utilizável em 18–20 Å:

     10 Å ->   0 transformadas      17 Å ->  19
     12 Å ->   0                    18 Å ->  40
     14 Å ->   1                    20 Å ->  77
     15 Å ->   7                    25 Å -> 295
     16 Å ->  14                   100 Å -> 8312   (controle: preparação OK)

O passo fino importa: 14 Å dá UMA transformada, mas 15 Å já dá 7 — e 7 pontos
de partida, com Full: True (50 refinamentos cada), rendem ~350 modelos. É pouco
e não é nada. Daí a faixa `marginal`: tentável, sabendo que o cluster pode não
se formar.

Ou seja: **o candidato estava geometricamente impossibilitado antes de a MD
começar**, e a MD de 46 h não tinha como perceber — ela simula E3 + PROTAC, sem
a PCSK9 na caixa. O teste que decide custa segundos e estava no fim da fila.

O ALCANCE EXIGIDO não depende do candidato: depende do par de proteínas e das
duas poses. Mede-se uma vez (é o `patchdock_span_scan.sh`) e vale para o
catálogo inteiro. Este script aplica o número.

O que ele mede, por candidato
-----------------------------
    alcance_A         o ALCANCE: maior distância âncora–âncora que a
                      topologia permite (limite de distance-geometry, 3 ms)
    n_ligacoes        caminho de ligações entre as duas âncoras
    ligacoes_faltando quanto falta para o requisito, pela taxa Å/ligação
                      medida no próprio conjunto
    alcance_mediana_A só com --amostrar: a mediana sobre confôrmeros. O
                      alcance diz se é POSSÍVEL; a mediana, se é confortável —
                      um PROTAC que só vence o vão estendido paga entropia.

NÃO use máximo amostrado como alcance: ele satura muito abaixo do real, e essa
foi a versão errada deste script. Os números estão em `limite_de_alcance`.

A âncora vem do `prosettac_anchors.ancora_do_head` — a mesma função que o
PRosettaC consome. Duas implementações de "qual átomo liga ao linker" acabariam
divergindo, e a distância medida aqui deixaria de ser a distância restringida
lá.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdDistGeom, rdmolops

from prosettac_anchors import ancora_do_head

RDLogger.DisableLog("rdApp.*")

ANGSTROM_POR_LIGACAO = 1.27   # projeção de uma ligação C–C numa cadeia anti


def limite_de_alcance(protac: Chem.Mol, i_a: int, i_b: int) -> float:
    """O ALCANCE: maior distância âncora–âncora que a topologia permite, em Å.

    Sai da matriz de limites de distance-geometry do RDKit, que resolve as
    restrições de comprimento de ligação, ângulo e anel — não é amostragem, é o
    teto, e sai em 3 ms.

    Por que não o máximo amostrado
    ------------------------------
    A primeira versão deste script usava o máximo sobre confôrmeros ETKDG. Ele
    SATURA muito abaixo do real, porque confôrmero gerado por aleatoriedade
    quase nunca cai na conformação estendida. Medido no PROTAC do
    `SC0013__WH022`:

        limite de distance-geometry .......... 18,03 Å   (3 ms)
        histograma do próprio PRosettaC ...... 18 Å
        amostrado com 200 confôrmeros ........ 14,44 Å   (27 s)
        amostrado com 800 confôrmeros ........ 14,44 Å   (108 s)

    Quadruplicar a amostra não move o máximo em 0,01 Å. E o erro não era
    acadêmico: com o máximo amostrado, este script classificou como `curto` o
    candidato que o PRosettaC aceitou e levou ao refinamento do Rosetta.

    O PRosettaC amostra até 200 confôrmeros POR BIN de distância — busca
    dirigida, não amostragem cega — e por isso ele encosta no teto. Como o que
    este portão precisa prever é o que o PRosettaC vai fazer, o teto é a medida
    certa.

    A ressalva fica no `alcance_mediana_A` (com --amostrar): um PROTAC que só
    alcança o vão na conformação estendida paga entropia para fazê-lo. O teto
    diz se é possível; a mediana, se é confortável.
    """
    b = rdDistGeom.GetMoleculeBoundsMatrix(protac)
    return float(max(b[i_a][i_b], b[i_b][i_a]))


def amostra_de_alcance(protac: Chem.Mol, i_a: int, i_b: int, n_confs: int,
                       semente: int, otimizar: bool):
    """Distribuição da distância entre as âncoras sobre confôrmeros.

    NÃO serve para o máximo (ver acima) — serve para a mediana, que diz quanto
    o alcance é confortável. Os H entram para a geometria ficar correta, e como
    o `AddHs` os acrescenta no FIM, os índices dos pesados não se movem.
    """
    mol = Chem.AddHs(protac)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = semente
    # Medido num PROTAC desta série (105 átomos com H, 120 confôrmeros):
    # 31 s numa thread, 8,6 s em todas. É o custo dominante do script, e num
    # catálogo de cem candidatos a diferença é de 50 min para 15.
    ps.numThreads = 0
    # Sem poda por RMSD: ela não acelerou nada aqui (31,4 s contra 31,2 s) e só
    # reduziria a amostra de onde sai o máximo, que é justamente o que se mede.
    ids = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=ps)
    if not ids:
        return None
    if otimizar:
        AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0, maxIters=200)
    d = []
    for cid in ids:
        c = mol.GetConformer(cid)
        pa, pb = c.GetAtomPosition(i_a), c.GetAtomPosition(i_b)
        d.append(pa.Distance(pb))
    return np.array(d)


def caminho_de_ligacoes(mol: Chem.Mol, i_a: int, i_b: int) -> int:
    caminho = rdmolops.GetShortestPath(mol, i_a, i_b)
    return max(len(caminho) - 1, 0)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--protacs", type=Path, required=True,
                    help="protac_candidates.csv do WP3")
    ap.add_argument("--e3-head", type=Path, required=True,
                    help=".sdf do recrutador na pose da E3")
    ap.add_argument("--warheads-sdf", type=Path, required=True,
                    help="diretório com <warhead_id>*.sdf nas poses do alvo")
    ap.add_argument("--requisito", type=float, default=18.0,
                    help="alcance exigido, do patchdock_span_scan.sh "
                         "(padrão 18: o regime com 40+ transformadas)")
    ap.add_argument("--minimo-absoluto", type=float, default=14.0,
                    help="onde aparece a PRIMEIRA solução; entre este valor e "
                         "o requisito o candidato é marginal, não aprovado")
    # A amostragem NÃO entra no portão (ela satura muito abaixo do alcance
    # real) e custa ~27 s por candidato, contra 3 ms do limite. Fica opcional,
    # para quando interessar saber se o alcance é confortável ou só possível.
    ap.add_argument("--amostrar", action="store_true",
                    help="medir também a mediana sobre confôrmeros (lento)")
    ap.add_argument("--confs", type=int, default=200)
    ap.add_argument("--semente", type=int, default=0xC0FFEE)
    ap.add_argument("--otimizar", action="store_true",
                    help="minimizar com MMFF (mais lento, pouco muda o máximo)")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.protacs.expanduser())
    print(f"{len(df)} PROTACs montados")

    head_e3 = Chem.MolFromMolFile(str(args.e3_head.expanduser()))
    if head_e3 is None:
        raise SystemExit(f"não consegui ler o head da E3: {args.e3_head}")

    dir_wh = args.warheads_sdf.expanduser()
    cache_wh: dict[str, Chem.Mol] = {}

    def head_warhead(wid: str):
        if wid in cache_wh:
            return cache_wh[wid]
        # o nome do arquivo varia com a etapa que o gerou (`WH023.sdf`,
        # `WH023_in_pcsk9_bo.sdf`); procurar é mais honesto que adivinhar
        achados = sorted(dir_wh.glob(f"{wid}*.sdf"))
        mol = None
        for a in achados:
            mol = Chem.MolFromMolFile(str(a))
            if mol is not None:
                break
        cache_wh[wid] = mol
        return mol

    linhas = []
    for n, r in enumerate(df.itertuples(), 1):
        protac = Chem.MolFromSmiles(r.protac_smiles)
        wh = head_warhead(r.warhead_id)
        base = {"candidate_id": r.candidate_id, "warhead_id": r.warhead_id}
        if protac is None or wh is None:
            linhas.append({**base, "situacao": "sem molécula legível"})
            continue

        _, i_e3, _, diag_e3 = ancora_do_head(head_e3, protac)
        _, i_wh, _, diag_wh = ancora_do_head(wh, protac)
        if i_e3 is None or i_wh is None:
            linhas.append({**base, "situacao": "âncora não determinada",
                           "diagnostico": (diag_e3 or "") + " | " +
                                          (diag_wh or "")})
            continue
        if i_e3 == i_wh:
            linhas.append({**base, "situacao": "as duas âncoras casaram no "
                                              "mesmo átomo"})
            continue

        plano = Chem.RemoveHs(protac)
        n_lig = caminho_de_ligacoes(plano, i_e3, i_wh)
        alcance = limite_de_alcance(protac, i_e3, i_wh)

        if alcance >= args.requisito:
            situacao = "alcança"
        elif alcance >= args.minimo_absoluto:
            situacao = "marginal"
        else:
            situacao = "curto"
        linha = {
            **base,
            "alcance_A": round(alcance, 2),
            "n_ligacoes": n_lig,
            "situacao": situacao,
            "protac_mw": getattr(r, "protac_mw", None),
            "protac_rotb": getattr(r, "protac_rotb", None),
        }
        if args.amostrar:
            d = amostra_de_alcance(protac, i_e3, i_wh, args.confs,
                                   args.semente, args.otimizar)
            if d is not None:
                linha["alcance_mediana_A"] = round(float(np.median(d)), 2)
                linha["alcance_amostrado_max_A"] = round(float(d.max()), 2)
                linha["n_confs"] = len(d)
        linhas.append(linha)
        if n % 25 == 0:
            print(f"  {n}/{len(df)}...", flush=True)

    res = pd.DataFrame(linhas)

    # Quantas ligações faltam, pela taxa MEDIDA e não pela teórica. A cadeia
    # toda anti rende 1,27 Å por ligação, mas confôrmero amostrado se enrola:
    # nesta série a taxa observada é ~0,98 Å/ligação. Usar 1,27 aqui subestima
    # o que falta — e o número existe justamente para escolher linker no
    # catálogo, onde errar por baixo devolve o mesmo zero do PatchDock.
    if "alcance_A" in res.columns:
        val = res.dropna(subset=["alcance_A", "n_ligacoes"])
        taxa = float(np.median(val["alcance_A"] / val["n_ligacoes"])) \
            if len(val) else ANGSTROM_POR_LIGACAO
        res["ligacoes_faltando"] = np.where(
            res["alcance_A"].notna(),
            np.ceil((args.requisito - res["alcance_A"]).clip(lower=0)
                    / taxa), np.nan)
        print(f"\n  taxa medida: {taxa:.2f} Å por ligação de cadeia "
              f"(teórica, cadeia anti: {ANGSTROM_POR_LIGACAO})")
        print(f"  ('ligacoes_faltando' usa essa taxa; é estimativa para "
              f"garimpar o\n  catálogo — os candidatos novos se medem, em 3 ms "
              f"cada.)")

    if "alcance_A" in res.columns:
        res = res.sort_values("alcance_A", ascending=False, na_position="last")
    out = args.out.expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out, index=False)

    print(f"\n  requisito: alcance ≥ {args.requisito:.0f} Å "
          f"(marginal a partir de {args.minimo_absoluto:.0f} Å)")
    contagem = res["situacao"].value_counts().to_dict()
    for k in ("alcança", "marginal", "curto"):
        if k in contagem:
            print(f"    {k:>9}: {contagem[k]}")
    for k, v in contagem.items():
        if k not in ("alcança", "marginal", "curto"):
            print(f"    {k}: {v}")

    mostrar = [c for c in ("candidate_id", "alcance_A", "n_ligacoes",
                           "alcance_mediana_A", "ligacoes_faltando",
                           "protac_mw", "protac_rotb", "situacao")
               if c in res.columns]
    print(f"\n  os 15 de maior alcance:")
    print(res[mostrar].head(15).to_string(index=False))
    print(f"\n  {out}")

    passam = res[res["situacao"] == "alcança"] if "situacao" in res else res
    print()
    if len(passam):
        print(f"  {len(passam)} candidato(s) já montado(s) ALCANÇAM o vão.")
        print(f"  Não é preciso voltar ao catálogo: re-ranqueie com o portão e")
        print(f"  o primeiro da fila é o próximo a ir para o PRosettaC.")
        print(f"    python rank_protacs.py ... --span {out} "
              f"--span-min {args.requisito:.0f}")
    else:
        marg = res[res["situacao"] == "marginal"] if "situacao" in res else []
        if len(marg):
            print(f"  Nenhum alcança {args.requisito:.0f} Å, mas {len(marg)} "
                  f"chega(m) à faixa marginal.")
            print(f"  Marginal quer dizer POUCAS transformadas no PatchDock — "
                  f"dá para tentar,")
            print(f"  sabendo que o cluster pode não se formar.")
        else:
            falta = int(res["ligacoes_faltando"].min()) \
                if "ligacoes_faltando" in res else None
            print(f"  NENHUM candidato montado alcança o vão.")
            if falta:
                print(f"  O melhor deles precisa de ~{falta} ligação(ões) a "
                      f"mais na cadeia.")
            print(f"  A conclusão é do WP2, não do WP3: o pool de linkers do")
            print(f"  Chemspace foi escolhido curto para este par de sítios.")
            print(f"  Refiltre o catálogo exigindo cadeia mais longa e remonte.")


if __name__ == "__main__":
    main()
