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
            aviso = "   <- ZERO BYTES, e é por isso que a etapa era pulada" \
                if tam == 0 else ""
            print(f"    {p.name}  ({tam} bytes){aviso}")

    if not args.aplicar:
        print(f"\n  para apagar:  python {Path(__file__).name} "
              f"--dir {d} --aplicar")
        return

    for p in fora:
        shutil.rmtree(p) if p.is_dir() else p.unlink()
    print(f"  apagados — o PRosettaC começa do zero, sem pular etapa por"
          f" causa de arquivo pela metade")


if __name__ == "__main__":
    main()
