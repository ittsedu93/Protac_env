#!/usr/bin/env python
"""
Fase 4b — quanto o par E3/alvo EXIGE de alcance, e quanto exige com folga.

Roda no env `mdtools`. Segundos (a varredura do PatchDock é que custa minutos).

    python span_requirement.py --scan <span_scan_*.csv> --out span_requirement.json

O que isto decide
-----------------
A varredura do `patchdock_span_scan.sh` devolve, para cada distância máxima
permitida entre as duas âncoras, quantas colocações das duas proteínas
sobrevivem. Este script lê essa curva e extrai três números:

    vao_minimo_A      onde aparece a PRIMEIRA transformada. É o piso absoluto,
                      e não serve de critério: com 14 Å o par CRBN/PCSK9 dava
                      UMA solução, e uma solução não sustenta cluster.

    vao_util_A        onde há transformadas suficientes para o Rosetta refinar
                      e o clustering ter material. Padrão: 20 transformadas.

    alcance_exigido_A o vão útil dividido pela razão mediana/teto — o critério
                      que a comparação entre dois candidatos reais produziu:

                        WH022: ponte no LIMITE do alcance (teto 18 Å, vão 18 Å)
                               -> 21 de 658 modelos passam a energia  (3,2%)
                        WH023: ponte com FOLGA  (teto 18 Å, vão 14 Å)
                               -> 134 de 1184 modelos passam         (11,3%)

                      Conformação estendida paga entropia e tensão, e isso
                      aparece como energia desfavorável: o filtro corta 97% dos
                      modelos quando a ponte se fecha no extremo. Exigir que o
                      TETO do linker exceda o vão pela razão mediana/teto põe a
                      ponte perto do meio da distribuição conformacional em vez
                      do extremo dela.

A razão mediana/teto medida na série (11,6 / 18,03) é 0,64. Ela é UM ponto
experimental, não uma constante da natureza — e é por isso que ela é um
parâmetro aqui, com o valor e a origem no JSON de saída, em vez de um número
enterrado no código.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

RAZAO_MEDIANA_TETO = 0.64   # medido em SC0013__WH022: mediana 11,6 / teto 18,03


def ler_varredura(csv_path: Path):
    """[(distancia_A, transformadas)] ordenado por distância."""
    pontos = []
    with open(csv_path) as fh:
        for r in csv.DictReader(fh):
            try:
                d = float(r["dist_thr_A"])
                t = int(float(r["transforms"]))
            except (KeyError, ValueError):
                continue
            pontos.append((d, t))
    return sorted(pontos)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scan", type=Path, required=True,
                    help="CSV do patchdock_span_scan.sh")
    ap.add_argument("--transformadas-uteis", type=int, default=20,
                    help="quantas transformadas fazem um vão ser 'útil'")
    ap.add_argument("--razao", type=float, default=RAZAO_MEDIANA_TETO,
                    help="mediana/teto do alcance; converte vão em teto exigido")
    ap.add_argument("--controle", type=float, default=100.0,
                    help="a distância que serve de controle na varredura")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    pontos = ler_varredura(args.scan.expanduser())
    if not pontos:
        raise SystemExit(f"não consegui ler pontos de {args.scan}")

    print(f"varredura: {len(pontos)} pontos de "
          f"{pontos[0][0]:.0f} a {pontos[-1][0]:.0f} Å")
    for d, t in pontos:
        marca = "   <- controle" if d >= args.controle else ""
        print(f"  {d:>6.0f} Å   {t:>6d} transformadas{marca}")

    # O controle vem primeiro: sem restrição efetiva o docking tem de devolver
    # milhares. Zero ali invalida a curva inteira, e ler a curva como química
    # nesse caso seria concluir sobre o linker a partir de um erro de preparação.
    controle = [t for d, t in pontos if d >= args.controle]
    if controle and max(controle) < 100:
        raise SystemExit(
            f"\n  [CONTROLE FALHOU] a {args.controle:.0f} Å o docking devolveu "
            f"{max(controle)} transformadas.\n"
            f"  Sem restrição efetiva era esperado milhares. O problema está na\n"
            f"  preparação (Init0/Init1) ou nos parâmetros de superfície, não no\n"
            f"  alcance do linker. A curva não significa nada assim.")

    uteis = [(d, t) for d, t in pontos if d < args.controle]
    minimo = next((d for d, t in uteis if t > 0), None)
    util = next((d for d, t in uteis if t >= args.transformadas_uteis), None)

    # O primeiro ponto já com transformadas significa que o piso REAL está
    # abaixo da faixa varrida: o número reportado é um limite superior, não a
    # medida. Conservador na direção certa (exige mais alcance do que talvez
    # baste), mas dizer "vão mínimo = 10 Å" quando não se olhou abaixo de 10
    # seria afirmar o que não foi medido.
    if uteis and uteis[0][1] > 0:
        piso_nao_medido = True
        print(f"\n  [ATENÇÃO] o ponto mais baixo da varredura ({uteis[0][0]:.0f} Å) "
              f"já tem {uteis[0][1]} transformadas.")
        print(f"  O piso real está ABAIXO disso e não foi medido. Os números a"
              f" seguir são\n  limites superiores: varra distâncias menores para"
              f" achar o piso — um vão\n  menor admite linker mais curto, e"
              f" linker mais curto é massa molecular\n  que não se gasta.")
    else:
        piso_nao_medido = False

    if minimo is None:
        raise SystemExit(
            "\n  Nenhuma distância abaixo do controle produziu transformada.\n"
            "  O par não admite colocação com restrição de distância: a decisão\n"
            "  é de E3 ou de sítio, e não de linker.")
    if util is None:
        util = uteis[-1][0]
        print(f"\n  [ATENÇÃO] nenhum ponto chegou a {args.transformadas_uteis} "
              f"transformadas; usando o maior medido ({util:.0f} Å), e o vão\n"
              f"  útil pode estar acima da faixa varrida.")

    exigido = util / args.razao

    d = {
        "vao_minimo_A": round(minimo, 1),
        "piso_abaixo_da_varredura": piso_nao_medido,
        "vao_util_A": round(util, 1),
        "transformadas_uteis": args.transformadas_uteis,
        "razao_mediana_teto": args.razao,
        "alcance_exigido_A": round(exigido, 1),
        "ligacoes_minimas": int(round(exigido / 1.20)),
        "origem": {
            "varredura": str(args.scan),
            "razao": "mediana/teto medida em SC0013__WH022 (11,6/18,03); um "
                     "ponto experimental, não constante — ajuste com --razao "
                     "quando houver mais candidatos medidos",
            "por_que_nao_o_minimo": "com 14 Å o par CRBN/PCSK9 dava UMA "
                                    "transformada; ponte no limite do alcance "
                                    "reprovou 97% dos modelos por energia",
        },
        "curva": [{"dist_thr_A": d_, "transforms": t_} for d_, t_ in pontos],
    }
    args.out.expanduser().parent.mkdir(parents=True, exist_ok=True)
    args.out.expanduser().write_text(json.dumps(d, indent=2,
                                                ensure_ascii=False) + "\n")

    print()
    print(f"  vão mínimo (1a transformada) .......... {d['vao_minimo_A']} Å"
          + ("   <- limite superior; o piso não foi medido" if piso_nao_medido
             else ""))
    print(f"  vão útil ({args.transformadas_uteis}+ transformadas) ......... "
          f"{d['vao_util_A']} Å")
    print(f"  ALCANCE EXIGIDO do linker (com folga) . {d['alcance_exigido_A']} Å")
    print(f"  ou seja, cadeia de ~{d['ligacoes_minimas']} ligações entre os "
          f"pontos de conjugação")
    print()
    print(f"  {args.out}")


if __name__ == "__main__":
    main()
