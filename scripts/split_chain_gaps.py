#!/usr/bin/env python
"""
Insere TER nas lacunas da cadeia, antes do pdb2gmx.

Só usa a biblioteca padrão. Segundos.

O problema que isto resolve
--------------------------
Cristais têm alças não resolvidas. O pdb2gmx, vendo uma cadeia contínua no
arquivo, cria uma ligação peptídica entre os resíduos que *ladeiam* a lacuna —
resíduos que estão a nanômetros de distância:

    Long Bond (3358-3360 = 1.39835 nm)

Uma ligação peptídica mede 0,133 nm. Uma mola harmônica esticada dez vezes
gera força suficiente para desmontar o sistema nos primeiros passos da
dinâmica, e o estouro raramente se apresenta como estouro:

    Step 5  LINCS WARNING ... 539 540  0.9983 863453.1250  0.1090
    CUDA error #700 (cudaErrorIllegalAddress)

O erro de CUDA é consequência de coordenadas que viraram lixo, não a causa. Na
CPU o mesmo sistema estouraria com outra mensagem.

A correção é dizer ao pdb2gmx que ali a cadeia termina: ele divide em cadeias
nos registros TER, e cada segmento ganha seus próprios términos. A carga total
não muda, porque cada corte acrescenta um NH3+ (+1) e um COO- (-1).

    python split_chain_gaps.py --receptor receptor_fixed.pdb \\
        --out receptor_capped.pdb [--ligante protac.sdf]
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

# Uma ligação peptídica C-N mede 1,33 Å. Acima de 2,5 Å não há ligação: há
# lacuna. O corte é generoso de propósito — geometria ruim de cristal pode
# esticar uma ligação real até ~1,8 Å sem que ela deixe de ser uma ligação.
CORTE_A = 2.5


def ler_residuos(pdb: Path):
    """Resíduos na ordem do arquivo: (chave, rótulo, {nome_atomo: (x,y,z)})."""
    residuos, ordem = {}, []
    for linha in pdb.read_text().splitlines():
        if not linha.startswith(("ATOM", "HETATM")):
            continue
        chave = (linha[21], linha[22:27])           # cadeia + resSeq + iCode
        if chave not in residuos:
            residuos[chave] = {"rot": f"{linha[17:20].strip()}{linha[22:27].strip()}",
                               "cadeia": linha[21], "atomos": {}}
            ordem.append(chave)
        residuos[chave]["atomos"][linha[12:16].strip()] = (
            float(linha[30:38]), float(linha[38:46]), float(linha[46:54]))
    return ordem, residuos


def dist(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def achar_lacunas(ordem, residuos, corte: float = CORTE_A):
    """Pares consecutivos cujo C-N é longo demais para ser ligação."""
    lacunas = []
    for i in range(len(ordem) - 1):
        k1, k2 = ordem[i], ordem[i + 1]
        r1, r2 = residuos[k1], residuos[k2]
        if r1["cadeia"] != r2["cadeia"]:
            continue                                 # já são cadeias distintas
        c, n = r1["atomos"].get("C"), r2["atomos"].get("N")
        if c is None or n is None:
            # sem backbone não há como medir; cortar é o lado seguro do erro
            lacunas.append((i, k1, k2, None))
            continue
        d = dist(c, n)
        if d > corte:
            lacunas.append((i, k1, k2, d))
    return lacunas


def centro_do_ligante(sdf: Path):
    """Centro geométrico dos átomos pesados do .sdf, sem depender do RDKit."""
    linhas = sdf.read_text().splitlines()
    if len(linhas) < 4:
        return None
    try:
        n = int(linhas[3][:3])
    except ValueError:
        return None
    pts = []
    for l in linhas[4:4 + n]:
        try:
            x, y, z = float(l[:10]), float(l[10:20]), float(l[20:30])
        except ValueError:
            continue
        if l[31:34].strip() != "H":
            pts.append((x, y, z))
    if not pts:
        return None
    return tuple(sum(p[i] for p in pts) / len(pts) for i in range(3))


def escrever(pdb: Path, out: Path, chaves_corte: set):
    """Copia o PDB inserindo TER depois de cada resíduo marcado."""
    saida, atual, serial_ter = [], None, 0
    for linha in pdb.read_text().splitlines():
        if linha.startswith(("ATOM", "HETATM")):
            chave = (linha[21], linha[22:27])
            if atual is not None and chave != atual and atual in chaves_corte:
                serial_ter += 1
                saida.append("TER")
            atual = chave
            saida.append(linha)
        elif linha.startswith("TER"):
            saida.append("TER")
            atual = None
        elif linha.startswith("END"):
            continue
        else:
            saida.append(linha)
    saida.append("TER")
    saida.append("END")
    out.write_text("\n".join(saida) + "\n")
    return serial_ter


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--receptor", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--ligante", type=Path,
                    help="para avisar se alguma lacuna cai perto do sítio")
    ap.add_argument("--corte", type=float, default=CORTE_A)
    args = ap.parse_args()

    pdb = args.receptor.expanduser()
    ordem, residuos = ler_residuos(pdb)
    print(f"  {len(ordem)} resíduos lidos de {pdb.name}")

    lacunas = achar_lacunas(ordem, residuos, args.corte)

    if not lacunas:
        print("  nenhuma lacuna na cadeia — nada a cortar")
        args.out.expanduser().write_text(pdb.read_text())
        return

    centro = centro_do_ligante(args.ligante.expanduser()) \
        if args.ligante and args.ligante.expanduser().exists() else None

    print(f"  {len(lacunas)} lacuna(s) na cadeia:")
    perto = []
    for _, k1, k2, d in lacunas:
        med = f"C-N = {d:.2f} Å" if d is not None else "sem backbone para medir"
        extra = ""
        if centro is not None:
            ds = [dist(p, centro) for p in residuos[k1]["atomos"].values()]
            ds += [dist(p, centro) for p in residuos[k2]["atomos"].values()]
            dmin = min(ds)
            extra = f", a {dmin:.1f} Å do ligante"
            if dmin < 12.0:
                perto.append((residuos[k1]['rot'], residuos[k2]['rot'], dmin))
        print(f"    {residuos[k1]['rot']} | {residuos[k2]['rot']}  ({med}{extra})")

    n = escrever(pdb, args.out.expanduser(), {k1 for _, k1, _, _ in lacunas})
    print(f"  {n} TER inserido(s) -> {args.out}")
    print("  a carga total não muda: cada corte soma um NH3+ (+1) e um COO- (-1)")

    if perto:
        print("\n  [ATENÇÃO] lacuna(s) a menos de 12 Å do ligante:")
        for a, b, dmin in perto:
            print(f"    {a} | {b}  ({dmin:.1f} Å)")
        print("  Os términos criados ali são carregados e ficam perto do sítio.")
        print("  Registre na tese; se a interação medida depender dessa região,")
        print("  o certo é modelar a alça em vez de cortá-la.")


if __name__ == "__main__":
    main()
