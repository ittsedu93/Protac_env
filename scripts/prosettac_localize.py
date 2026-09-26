#!/usr/bin/env python
"""
Traz as entradas do PRosettaC para dentro do diretório de trabalho.

Só biblioteca padrão. Instantâneo.

    python prosettac_localize.py --dir <run_dir>

O problema que isto resolve
---------------------------
O PRosettaC espera as entradas NO diretório de trabalho, com nomes relativos —
é assim que o exemplo da própria ferramenta as declara. Com caminho absoluto
apontando para outra pasta, ele se parte ao meio: o `clean_pdb` do Rosetta
escreve `<struct>_<chain>.fasta` no diretório ATUAL, mas o `rosetta.py` o
procura ao lado do PDB de ENTRADA. Os dois lugares deixam de coincidir:

    4TZ4_receptor C   361 --- --- MOD --- OK      <- a limpeza funcionou
    FileNotFoundError: .../prep/CRBN_4TZ4/4TZ4_receptor_C.fasta

O arquivo existe — no diretório de trabalho. Copiar as entradas para lá e usar
nomes relativos faz os dois lugares serem o mesmo lugar.

Idempotente: um config que já esteja com nomes relativos não é tocado.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

CAMPOS = ("Structures:", "Heads:", "Protac:")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--config", default="prosetta_config.txt")
    args = ap.parse_args()

    d = args.dir.expanduser().resolve()
    cfg = d / args.config
    if not cfg.exists():
        raise SystemExit(f"não achei {cfg}")

    linhas, copiados, ja_ok = [], [], 0
    for linha in cfg.read_text().splitlines():
        campo = linha.split(" ", 1)[0] if " " in linha else linha
        if campo not in CAMPOS:
            linhas.append(linha)
            continue
        chave, valores = linha.split(": ", 1)
        novos = []
        for v in valores.split():
            p = Path(v)
            if not p.is_absolute():
                ja_ok += 1
                novos.append(v)
                continue
            if not p.exists():
                raise SystemExit(f"entrada não existe: {p}")
            alvo = d / p.name
            if alvo.resolve() != p.resolve():
                shutil.copy2(p, alvo)
                copiados.append(f"{p.name}  <- {p.parent}")
            novos.append(p.name)
        linhas.append(f"{chave}: {' '.join(novos)}")

    if not copiados:
        print(f"  entradas já locais ({ja_ok} arquivos) — config intacto")
        return

    shutil.copy2(cfg, cfg.with_suffix(".txt.bak"))
    cfg.write_text("\n".join(linhas) + "\n")
    print(f"  {len(copiados)} entrada(s) copiada(s) para o diretório de trabalho:")
    for c in copiados:
        print(f"    {c}")
    print(f"  config reescrito com nomes relativos "
          f"(original em {cfg.with_suffix('.txt.bak').name})")


if __name__ == "__main__":
    main()
