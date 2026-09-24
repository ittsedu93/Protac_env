#!/usr/bin/env python
"""
Fase 7 — ranqueia os PROTACs montados e escolhe o candidato para a MD.

Roda no env `mdtools`. Segundos.

Não existe um número que diga "este PROTAC vai funcionar". O que existe são
indicadores que, combinados, separam candidatos plausíveis de implausíveis — e
o ranking abaixo é explícito sobre o peso de cada um, para poder ser defendido
ou contestado na tese.

Critérios, e por que cada um
----------------------------
1. **Eficiência de ligante do warhead** — o score bruto do Vina correlacionou
   -0,92 com o número de átomos pesados nesta série: ordenar por ele é ordenar
   por massa molecular. LE corrige, e concorda com os filtros de propriedade.
2. **Geometria do linker** (`frac_sem_clash`) — fração dos confôrmeros que
   acomoda o exit vector sem colidir com a E3 ligase. Um linker que só encaixa
   num confôrmero está forçado, por melhor que pontue.
3. **Estabilidade da pose do warhead** (`std_score` entre runs independentes).
4. **Penalidade de tamanho do PROTAC** — massa e ligações rotacionáveis. Um
   PROTAC de 900 Da com 25 torções é uma molécula difícil por construção, e o
   custo disso aparece na MD e na permeabilidade, não no docking.

Os quatro entram normalizados (0–1) com pesos declarados em PESOS.

    python rank_protacs.py \\
        --protacs ~/pipeline/wp3/protac_candidates.csv \\
        --warhead-ranking ~/pipeline/warhead_ranking.csv \\
        --linker-ranking ~/pipeline/wp2/linker_ranking.csv \\
        --subcomplexes ~/pipeline/wp2/subcomplexes_manifest.json \\
        --out ~/pipeline/protac_ranking.csv
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

PESOS = {
    "le_warhead": 0.35,       # afinidade corrigida por tamanho
    "geometria_linker": 0.30,  # acomodação do exit vector
    "estabilidade_pose": 0.20,  # consistência entre runs
    "tamanho_protac": 0.15,   # penalidade de massa/flexibilidade
}

# faixas para a penalidade de tamanho: acima destes valores o PROTAC entra em
# território onde a permeabilidade celular é problema conhecido
MW_CONFORTAVEL = 800.0
MW_LIMITE = 1000.0
ROTB_CONFORTAVEL = 18
ROTB_LIMITE = 28


def normalizar(s: pd.Series, maior_melhor: bool = True) -> pd.Series:
    s = pd.to_numeric(s, errors="coerce")
    lo, hi = s.min(), s.max()
    if not np.isfinite(lo) or not np.isfinite(hi) or hi - lo < 1e-12:
        return pd.Series(0.5, index=s.index)
    z = (s - lo) / (hi - lo)
    return z if maior_melhor else 1.0 - z


def penalidade_tamanho(mw: pd.Series, rotb: pd.Series) -> pd.Series:
    p_mw = ((mw - MW_CONFORTAVEL) / (MW_LIMITE - MW_CONFORTAVEL)).clip(0, 1)
    p_rb = ((rotb - ROTB_CONFORTAVEL) / (ROTB_LIMITE - ROTB_CONFORTAVEL)).clip(0, 1)
    return 1.0 - (0.5 * p_mw + 0.5 * p_rb)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--protacs", type=Path, required=True)
    ap.add_argument("--warhead-ranking", type=Path, required=True)
    ap.add_argument("--linker-ranking", type=Path)
    ap.add_argument("--subcomplexes", type=Path)
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.protacs.expanduser())
    wh = pd.read_csv(args.warhead_ranking.expanduser())
    print(f"{len(df)} PROTACs montados")

    cols_wh = [c for c in ("warhead_id", "ligand_efficiency", "best_score",
                           "std_score", "r1", "r2", "r3", "mw", "clogp",
                           "net_charge_ph74", "cos_vetor_saida")
               if c in wh.columns]
    df = df.merge(wh[cols_wh], on="warhead_id", how="left",
                  suffixes=("", "_wh"))

    # geometria do linker, via manifesto dos sub-complexos
    if args.subcomplexes and args.subcomplexes.expanduser().exists():
        sub = pd.DataFrame(json.loads(args.subcomplexes.expanduser().read_text()))
        if "frac_sem_clash" in sub.columns:
            df = df.merge(sub[["subcomplex_id", "frac_sem_clash", "linker_smiles"]],
                          on="subcomplex_id", how="left")
    if "frac_sem_clash" not in df.columns:
        print("  [aviso] sem frac_sem_clash: geometria do linker entra neutra")
        df["frac_sem_clash"] = 0.5

    # --- componentes normalizados ---
    df["_le"] = normalizar(df.get("ligand_efficiency", pd.Series(0.5, index=df.index)))
    df["_geo"] = normalizar(df["frac_sem_clash"])
    df["_est"] = normalizar(df.get("std_score", pd.Series(0.5, index=df.index)),
                            maior_melhor=False)
    df["_tam"] = penalidade_tamanho(
        pd.to_numeric(df["protac_mw"], errors="coerce"),
        pd.to_numeric(df["protac_rotb"], errors="coerce"))

    df["score_composto"] = (
        PESOS["le_warhead"] * df["_le"]
        + PESOS["geometria_linker"] * df["_geo"]
        + PESOS["estabilidade_pose"] * df["_est"]
        + PESOS["tamanho_protac"] * df["_tam"]
    ).round(4)

    df = df.sort_values("score_composto", ascending=False).reset_index(drop=True)
    df.insert(0, "rank", df.index + 1)

    out = args.out.expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)

    mostrar = [c for c in ("rank", "candidate_id", "r1", "r2", "r3",
                           "protac_mw", "protac_rotb", "ligand_efficiency",
                           "frac_sem_clash", "std_score", "score_composto")
               if c in df.columns]
    print(f"\nPesos: {PESOS}")
    print(f"\nTop {args.top}:")
    print(df[mostrar].head(args.top).to_string(index=False))

    melhor = df.iloc[0]
    print(f"\n{'=' * 66}")
    print(f"CANDIDATO PARA A MD: {melhor['candidate_id']}")
    print(f"{'=' * 66}")
    print(f"  MW {melhor['protac_mw']:.0f} Da | {melhor['protac_rotb']} torções")
    print(f"  warhead {melhor['warhead_id']} (R1={melhor.get('r1')}, "
          f"R2={melhor.get('r2')}, R3={melhor.get('r3')})")
    print(f"  linker  {melhor.get('linker_smiles', '?')}")
    print(f"  score composto {melhor['score_composto']:.3f}")
    print(f"\n  SMILES: {melhor['protac_smiles']}")

    # o próximo da fila, para quando a MD reprovar o primeiro
    (out.parent / "md_candidato.json").write_text(json.dumps({
        "escolhido": {k: (v.item() if hasattr(v, "item") else v)
                      for k, v in melhor.to_dict().items()},
        "proximos": [{k: (v.item() if hasattr(v, "item") else v)
                      for k, v in r.to_dict().items()}
                     for _, r in df.iloc[1:6].iterrows()],
    }, indent=2, default=str))
    print(f"\n  se a MD reprovar, os próximos da fila estão em "
          f"{out.parent / 'md_candidato.json'}")
    print(f"\n{out}")


if __name__ == "__main__":
    main()
