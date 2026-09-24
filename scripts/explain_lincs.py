#!/usr/bin/env python
"""
Diz QUAL resíduo está explodindo, a partir dos índices do LINCS.

Só usa a biblioteca padrão. Instantâneo.

O problema que isto resolve
--------------------------
O GROMACS relata instabilidade por número de átomo:

    Step 29  LINCS WARNING
    rms 0.009537, max 0.520117 (between atoms 539 and 540)

`539` não diz nada. Saber que é o H de uma cadeia lateral reconstruída, ou um
átomo do ligante, ou uma água encravada, é a diferença entre corrigir a causa e
mexer em parâmetros de integração no escuro. O mapeamento está no .gro, e é
trabalho de um segundo — mas só se alguém o fizer.

    python explain_lincs.py --gro complexo.gro --log mdrun_nvt.out
"""

from __future__ import annotations

import argparse
import re
from collections import Counter
from pathlib import Path


def ler_gro(gro: Path):
    """Átomo (1-based) -> (resíduo, nome do átomo). Formato fixo do .gro."""
    linhas = gro.read_text().splitlines()
    n = int(linhas[1])
    atomos = {}
    for i, l in enumerate(linhas[2:2 + n], start=1):
        atomos[i] = (f"{l[5:10].strip()}{l[0:5].strip()}", l[10:15].strip())
    return atomos, n


def indices_do_log(log: Path):
    """Todos os átomos citados nos avisos do LINCS, na ordem de aparição."""
    texto = log.read_text(errors="ignore")
    idx = []
    # "(between atoms 539 and 540)"
    for a, b in re.findall(r"between atoms\s+(\d+)\s+and\s+(\d+)", texto):
        idx += [int(a), int(b)]
    # as linhas da tabela "atom 1 atom 2 angle ..."
    for bloco in re.findall(r"bonds that rotated more than 30 degrees:\n"
                            r"\s*atom 1 atom 2.*\n((?:\s*\d+\s+\d+\s+[\d.]+.*\n)+)",
                            texto):
        for a, b in re.findall(r"^\s*(\d+)\s+(\d+)\s+[\d.]+", bloco, re.M):
            idx += [int(a), int(b)]
    # "on atom 570" (Fmax da minimização)
    for a in re.findall(r"on atom\s*=?\s*(\d+)", texto):
        idx.append(int(a))
    return idx


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gro", type=Path, required=True)
    ap.add_argument("--log", type=Path, required=True)
    ap.add_argument("--n-proteina", type=int, default=None,
                    help="último átomo da proteína, para separar do ligante")
    args = ap.parse_args()

    gro, log = args.gro.expanduser(), args.log.expanduser()
    if not (gro.exists() and log.exists()):
        return

    atomos, n = ler_gro(gro)
    idx = indices_do_log(log)
    if not idx:
        return

    print("      Quem está explodindo, pelos índices do LINCS:")
    contagem = Counter(idx)
    residuos = Counter()
    for i, vezes in contagem.most_common():
        if i not in atomos:
            print(f"        átomo {i}: fora do .gro ({n} átomos) — "
                  f"o .gro e o .tpr não são do mesmo sistema?")
            continue
        res, nome = atomos[i]
        residuos[res] += vezes
        print(f"        átomo {i:6d} = {nome:<5s} de {res:<10s} "
              f"({vezes}x nos avisos)")

    print("      Resíduos envolvidos, do mais citado ao menos:")
    for res, vezes in residuos.most_common(5):
        print(f"        {res}  ({vezes}x)")

    # o que fazer depende de QUEM é
    alvos = set(residuos)
    if any(r.startswith(("SOL", "HOH", "WAT")) for r in alvos):
        print("      >>> Há ÁGUA entre os culpados: o solvate encaixou uma")
        print("      >>> molécula perto demais. Raios de VdW são adivinhados")
        print("      >>> por nome de átomo, e nomes não padrão atrapalham.")
    if any(r.startswith(("PTC", "UNL", "LIG", "MOL")) for r in alvos):
        print("      >>> O LIGANTE está entre os culpados: a geometria de")
        print("      >>> partida do PROTAC é o suspeito, não o protocolo.")
        print("      >>> Reembeber com mais tentativas, ou relaxá-lo isolado")
        print("      >>> antes de montar o complexo, é o caminho.")
    else:
        print("      >>> Só PROTEÍNA entre os culpados. Se for um resíduo")
        print("      >>> reconstruído (veja o log do fix_receptor_for_md.py),")
        print("      >>> o rotâmero escolhido caiu num mínimo local que a")
        print("      >>> minimização não consegue deixar.")


if __name__ == "__main__":
    main()
