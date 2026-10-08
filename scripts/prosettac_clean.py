#!/usr/bin/env python
"""
Devolve o diretório de trabalho do PRosettaC ao estado de entrada.

Só biblioteca padrão. Instantâneo.

    python prosettac_clean.py --dir <run_dir>          # lista o que apagaria
    python prosettac_clean.py --dir <run_dir> --aplicar

O problema que isto resolve
---------------------------
Uma tentativa que morre no meio deixa produtos pela metade, e o PRosettaC
decide o que fazer pela EXISTÊNCIA do arquivo, não pelo conteúdo:

    -rw-rw-r--  1 ...      0 Sep 27 09:30 random_sampling.sdf
    -rw-rw-r--  1 ...   6582 Sep 27 09:30 PT0.params
    -rw-rw-r--  1 ... 476284 Sep 27 09:30 Init0.pdb

Zero byte é o pior caso: o arquivo existe, a etapa que o produziria é pulada, e
o que vem depois opera sobre nada. A primeira versão da limpeza deste pipeline
apagava só `*_[A-Z].pdb`, `*.fasta`, `log.txt` e `*_H.sdf` — os produtos que eu
tinha visto naquele dia. Nomear produtos um a um é uma lista que sempre está
incompleta.

A regra aqui é a inversa, e não depende de eu conhecer a ferramenta: o que fica
é o que é ENTRADA — os arquivos que o próprio config declara, o config, e os
backups. Todo o resto é produto e sai. Um diretório de trabalho não tem
conteúdo permanente além das entradas.

    MANTÉM  o config, seu .bak, e tudo que está em Structures/Heads/Protac
            (mais `*.bak` e `*.velho`, que são cópias nossas do original)
    APAGA   qualquer outra coisa, arquivo ou diretório

Por isso o `--seco` é o padrão: quem apaga por regra e não por lista tem de
poder ver a lista antes.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

CAMPOS_DE_ENTRADA = ("Structures", "Heads", "Protac", "Linkers")
SUFIXOS_PRESERVADOS = (".bak", ".velho")
# Arquivos que NÓS escrevemos. Saem junto, mas não são etapa pulada do
# PRosettaC — dizer que um log de zero byte fez pular uma etapa seria inventar
# um mecanismo, e a saída deste script é lida como diagnóstico.
NOSSOS = ("run_prosettac.log", ".heads.out", ".ancoras.out")

# ACIMA DISTO A LIMPEZA PARA E PERGUNTA.
#
# A regra inversa ("fica a entrada, sai o resto") é certa, e é justamente por
# ser certa que ela é perigosa: ela apaga sem precisar entender a ferramenta.
# A guarda que existia era "há algo rodando?", com três sinais — processo vivo,
# fila com jobs, arquivo escrito nos últimos 15 min. Em 08/10 os três estavam
# negativos: o orquestrador tinha morrido 3 dias antes, a fila havia drenado
# sozinha (o SLURM não depende dele), e nada era escrito desde então. Mesmo
# assim o diretório tinha 271 mil arquivos: os 2682 jobs de conformação do
# linker TODOS concluídos, sete dias de máquina, faltando só o agrupamento
# final. "Nada rodando" não é "nada feito".
#
# Então a guarda nova não pergunta se algo roda — pergunta QUANTO TRABALHO
# está aqui. Uma tentativa que morreu no início deixa um punhado de arquivos;
# uma que morreu no fim deixa centenas de milhares. A primeira é resto, a
# segunda é patrimônio, e a diferença é contável.
LIMITE_DE_RESTO = 500


def entradas_declaradas(cfg: Path) -> set[str]:
    """Os nomes que o config aponta como entrada, e mais nada."""
    nomes = {cfg.name, cfg.name + ".bak", cfg.with_suffix(".txt.bak").name}
    for linha in cfg.read_text().splitlines():
        if ": " not in linha:
            continue
        chave, valores = linha.split(": ", 1)
        if chave.strip() not in CAMPOS_DE_ENTRADA:
            continue
        for v in valores.split():
            nomes.add(Path(v).name)
    return nomes


def a_apagar(d: Path, manter: set[str]) -> list[Path]:
    fora = []
    for p in sorted(d.iterdir()):
        if p.name in manter:
            continue
        if p.suffix in SUFIXOS_PRESERVADOS:
            continue
        fora.append(p)
    return fora


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--config", default="prosetta_config.txt")
    ap.add_argument("--aplicar", action="store_true")
    ap.add_argument("--aceitar-perda", action="store_true",
                    help=f"apagar mesmo com mais de {LIMITE_DE_RESTO} arquivos "
                         f"de produto (sete dias de máquina caberiam aqui)")
    ap.add_argument("--limite", type=int, default=LIMITE_DE_RESTO,
                    help="quantos arquivos ainda contam como 'resto'")
    args = ap.parse_args()

    d = args.dir.expanduser().resolve()
    cfg = d / args.config
    if not cfg.exists():
        raise SystemExit(f"não achei {cfg} — sem config eu não sei o que é"
                         f" entrada, e sem isso não apago nada")

    manter = entradas_declaradas(cfg)
    fora = a_apagar(d, manter)

    if not fora:
        print("  diretório já limpo — só as entradas")
        return

    print(f"  {len(fora)} resto(s) de tentativa anterior:")
    for p in fora:
        if p.is_dir():
            n = sum(1 for _ in p.rglob("*"))
            print(f"    {p.name}/  ({n} arquivos)")
        else:
            tam = p.stat().st_size
            if p.name in NOSSOS:
                aviso = "   <- log da tentativa anterior, nosso"
            elif tam == 0:
                aviso = "   <- ZERO BYTES, e é por isso que a etapa era pulada"
            else:
                aviso = ""
            print(f"    {p.name}  ({tam} bytes){aviso}")

    # Quanto trabalho está aqui, contado antes de qualquer decisão.
    n_arq, n_bytes, mais_novo = 0, 0, 0.0
    for p in fora:
        alvos = p.rglob("*") if p.is_dir() else [p]
        for q in alvos:
            if q.is_file():
                n_arq += 1
                st = q.stat()
                n_bytes += st.st_size
                mais_novo = max(mais_novo, st.st_mtime)
    print(f"\n  total: {n_arq} arquivos, {n_bytes / 1e9:.2f} GB")
    if mais_novo:
        import time
        print(f"  escrita mais recente: "
              f"{time.strftime('%d/%m %H:%M', time.localtime(mais_novo))}"
              f"  (há {(time.time() - mais_novo) / 3600:.0f} h)")

    if not args.aplicar:
        print(f"\n  para apagar:  python {Path(__file__).name} "
              f"--dir {d} --aplicar")
        return

    if n_arq > args.limite and not args.aceitar_perda:
        raise SystemExit(
            f"\n*** NÃO VOU APAGAR: {n_arq} arquivos é produto, não resto.\n"
            f"*** O limite é {args.limite}; acima dele isto é trabalho feito.\n"
            f"***\n"
            f"*** Uma execução que morreu no INÍCIO deixa um punhado de\n"
            f"*** arquivos. Esta deixou {n_arq} ({n_bytes / 1e9:.2f} GB), e um\n"
            f"*** diretório assim costuma estar a uma etapa do fim — com as\n"
            f"*** etapas caras já pagas e só o agrupamento faltando.\n"
            f"***\n"
            f"*** Para RETOMAR sem apagar (o PRosettaC pula etapa cujo\n"
            f"*** produto existe, e aqui isso é o que se quer):\n"
            f"***     PROSETTAC_SEM_LIMPAR=1 bash scripts/run_prosettac.sh <cand>\n"
            f"***\n"
            f"*** Para apagar de propósito e começar do zero, sabendo o custo:\n"
            f"***     python {Path(__file__).name} --dir {d} --aplicar"
            f" --aceitar-perda")

    for p in fora:
        shutil.rmtree(p) if p.is_dir() else p.unlink()
    print(f"  apagados — o PRosettaC começa do zero, sem pular etapa por"
          f" causa de arquivo pela metade")


if __name__ == "__main__":
    main()
