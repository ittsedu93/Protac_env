#!/usr/bin/env python
"""
Monta o topol.top do complexo proteína + ligante, na ordem que o GROMACS exige.

Roda no env `mdtools`. Segundos.

O problema que isto resolve
---------------------------
O `.itp` que o acpype gera traz `[ atomtypes ]` e `[ moleculetype ]` no mesmo
arquivo. O GROMACS exige que TODO `[ atomtypes ]` venha antes do primeiro
`[ moleculetype ]` — e o topol.top do pdb2gmx já define a proteína logo depois
do include do campo de força. Um `#include` único do arquivo do acpype, em
qualquer posição, viola a ordem:

    Fatal error: Syntax error - File PTC_GMX.itp, line 3
    Last line read: '[ atomtypes ]'
    Invalid order for directive atomtypes

A solução é separar as duas seções em arquivos distintos e incluir cada uma no
lugar certo: os tipos de átomo logo após o campo de força, a molécula depois da
proteína e antes da água.

    python build_topology.py --topol topol.top --ligand-itp PTC.acpype/PTC_GMX.itp \\
        --resname PTC --ligand-gro PTC.acpype/PTC_GMX.gro
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

# seções que pertencem ao bloco global, antes de qualquer moleculetype
GLOBAIS = {"atomtypes", "bondtypes", "angletypes", "dihedraltypes",
           "constrainttypes", "nonbond_params", "pairtypes", "implicit_genborn_params"}
# seções que o campo de força já define e que não podem ser redeclaradas
DESCARTAR = {"defaults"}


def separar_itp(itp: Path):
    """Divide o .itp do acpype em (bloco global, bloco da molécula)."""
    texto = Path(itp).read_text()
    linhas = texto.splitlines()

    blocos, atual, nome = [], [], None
    for l in linhas:
        m = re.match(r"\s*\[\s*([a-zA-Z_0-9]+)\s*\]", l)
        if m:
            if nome is not None or atual:
                blocos.append((nome, atual))
            nome, atual = m.group(1).lower(), [l]
        else:
            atual.append(l)
    if nome is not None or atual:
        blocos.append((nome, atual))

    globais, molecula = [], []
    for nome, corpo in blocos:
        if nome is None:                      # cabeçalho/comentários
            globais.append("\n".join(corpo))
            continue
        if nome in DESCARTAR:
            continue
        (globais if nome in GLOBAIS else molecula).append("\n".join(corpo))
    return "\n".join(globais).strip(), "\n".join(molecula).strip()


def posre_ligante(gro: Path, resname: str, out: Path, fc: int = 1000):
    """Restrições de posição nos átomos pesados do ligante, para o equilíbrio.

    Sem isto o ligante fica solto durante o NVT/NPT com -DPOSRES enquanto a
    proteína está travada, e pode escorregar do sítio antes da produção.
    """
    linhas = Path(gro).read_text().splitlines()
    n = int(linhas[1])
    corpo = ["; restrições do ligante, geradas por build_topology.py",
             "[ position_restraints ]", "; atom  type      fx      fy      fz"]
    contados = 0
    for i, l in enumerate(linhas[2:2 + n], start=1):
        nome_atomo = l[10:15].strip()
        if nome_atomo.upper().startswith("H"):
            continue
        corpo.append(f"{i:6d}     1 {fc:7d} {fc:7d} {fc:7d}")
        contados += 1
    out.write_text("\n".join(corpo) + "\n")
    return contados


def montar(topol: Path, itp: Path, resname: str, gro: Path | None):
    texto = topol.read_text()
    if f"{resname}_atomtypes.itp" in texto:
        print("  topol.top já montado — nada a fazer")
        return

    globais, molecula = separar_itp(itp)
    if not molecula:
        raise SystemExit(f"não achei [ moleculetype ] em {itp}")

    dir_ = topol.parent
    f_glob = dir_ / f"{resname}_atomtypes.itp"
    f_mol = dir_ / f"{resname}_moleculetype.itp"
    f_glob.write_text(globais + "\n")
    f_mol.write_text(molecula + "\n")
    print(f"  {f_glob.name}: {len(globais.splitlines())} linhas")
    print(f"  {f_mol.name}: {len(molecula.splitlines())} linhas")

    # restrições do ligante, referenciadas dentro do bloco da molécula
    if gro and Path(gro).exists():
        f_posre = dir_ / f"posre_{resname}.itp"
        n = posre_ligante(Path(gro), resname, f_posre)
        f_mol.write_text(
            molecula
            + f'\n\n#ifdef POSRES\n#include "{f_posre.name}"\n#endif\n')
        print(f"  posre_{resname}.itp: {n} átomos pesados restringidos")

    linhas = texto.splitlines()

    # 1) tipos de átomo: logo após o include do campo de força
    idx_ff = next((i for i, l in enumerate(linhas)
                   if "forcefield.itp" in l and l.strip().startswith("#include")),
                  None)
    if idx_ff is None:
        raise SystemExit("não achei o #include do campo de força em topol.top")
    linhas.insert(idx_ff + 1, f'\n; tipos de atomo do ligante, antes de '
                              f'qualquer molecula\n'
                              f'#include "{f_glob.name}"')

    # 2) a molécula: antes da água, já depois da proteína
    idx_agua = next((i for i, l in enumerate(linhas)
                     if l.strip().startswith("#include")
                     and re.search(r"(tip3p|spce?|water)\.itp", l)), None)
    if idx_agua is None:
        idx_agua = next((i for i, l in enumerate(linhas)
                         if l.strip().startswith("[ system ]")), len(linhas))
    linhas.insert(idx_agua, f'\n; ligante\n#include "{f_mol.name}"\n')

    # 3) contagem em [ molecules ], antes de qualquer SOL
    idx_mols = next((i for i, l in enumerate(linhas)
                     if l.strip().startswith("[ molecules ]")), None)
    if idx_mols is None:
        raise SystemExit("não achei [ molecules ] em topol.top")
    fim = len(linhas)
    idx_sol = next((i for i in range(idx_mols + 1, fim)
                    if linhas[i].split() and linhas[i].split()[0] == "SOL"), fim)
    linhas.insert(idx_sol, f"{resname}                 1")

    topol.write_text("\n".join(linhas) + "\n")
    print(f"  topol.top montado")


def verificar(topol: Path):
    """Confere a ordem: nenhum atomtypes depois do primeiro moleculetype."""
    texto = topol.read_text()
    dir_ = topol.parent
    # expande os includes locais para ver a ordem real
    def expandir(t, prof=0):
        if prof > 3:
            return t
        saida = []
        for l in t.splitlines():
            m = re.match(r'\s*#include\s+"([^"]+)"', l)
            if m and (dir_ / m.group(1)).exists():
                saida.append(expandir((dir_ / m.group(1)).read_text(), prof + 1))
            else:
                saida.append(l)
        return "\n".join(saida)

    # tira os comentários antes de procurar diretivas: o GROMACS ignora tudo
    # depois de ';', e uma diretiva citada num comentário não é uma diretiva
    plano = "\n".join(l.split(";", 1)[0] for l in expandir(texto).splitlines())
    primeiro_mol = None
    for m in re.finditer(r"\[\s*(atomtypes|moleculetype)\s*\]", plano):
        if m.group(1) == "moleculetype" and primeiro_mol is None:
            primeiro_mol = m.start()
        elif m.group(1) == "atomtypes" and primeiro_mol is not None:
            return False, ("[ atomtypes ] aparece depois do primeiro "
                           "[ moleculetype ] — o GROMACS vai recusar")
    return True, "ordem das diretivas OK"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--topol", type=Path, required=True)
    ap.add_argument("--ligand-itp", type=Path, required=True)
    ap.add_argument("--ligand-gro", type=Path)
    ap.add_argument("--resname", default="PTC")
    args = ap.parse_args()

    montar(args.topol.expanduser(), args.ligand_itp.expanduser(),
           args.resname, args.ligand_gro)
    ok, msg = verificar(args.topol.expanduser())
    print(f"  verificação: {msg}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
