#!/usr/bin/env python
"""
Monta o grupos.ndx do GROMACS descobrindo os números dos grupos.

Só usa a biblioteca padrão — roda em qualquer env. Segundos.

O problema que isto resolve
--------------------------
Os `tc-grps` dos .mdp precisam nomear grupos que existam no index. Escrever as
fusões à mão exige adivinhar DOIS valores frágeis: o número do grupo, que muda
com a composição do sistema, e o nome do resíduo do ligante. E o nome do
ligante é uma armadilha: o acpype batiza a *moleculetype* como se pede em `-b`
(PTC), mas as coordenadas escritas pelo RDKit trazem o resíduo como `UNL`,
então o make_ndx cria um grupo `UNL` e um `Protein_PTC` nunca aparece:

    Group    12 (          Other) has   105 elements
    Group    13 (            UNL) has   105 elements

Aqui nada é adivinhado. Roda-se `make_ndx` só para listar, leem-se os números
pelos nomes, as fusões são montadas com os números reais, e o índice resultante
é verificado antes de ser aceito: os dois grupos têm de existir e a soma dos
dois tem de cobrir o sistema inteiro (átomo sem grupo de acoplamento = o
grompp recusa a simulação).

    python build_index.py --gro neutro.gro --out grupos.ndx
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

GMX = "/usr/local/gromacs/bin/gmx"

# grupos que o GROMACS cria por conta própria e que NÃO são o ligante
NAO_LIGANTE = {
    "system", "protein", "protein-h", "c-alpha", "backbone", "mainchain",
    "mainchain+cb", "mainchain+h", "sidechain", "sidechain-h", "prot-masses",
    "non-protein", "water", "sol", "non-water", "ion", "ions",
    "water_and_ions", "other", "na", "cl", "k", "mg", "ca", "zn", "dna", "rna",
}

RE_GRUPO = re.compile(
    r"^\s*(?:Group\s+)?(\d+)\s*\(\s*(\S+)\s*\)\s+has\s+(\d+)\s+elements", re.M)


def rodar_make_ndx(gro: Path, gmx: str, entrada: str, out: str):
    """make_ndx com um roteiro de comandos na entrada padrão."""
    return subprocess.run([gmx, "make_ndx", "-f", str(gro), "-o", out],
                          input=entrada, capture_output=True, text=True)


def listar_grupos(gro: Path, gmx: str) -> dict[int, tuple[str, int]]:
    """Números, nomes e tamanhos dos grupos, lidos da saída do próprio gmx."""
    r = rodar_make_ndx(gro, gmx, "q\n", "/dev/null")
    saida = (r.stdout or "") + (r.stderr or "")
    grupos = {int(n): (nome, int(tam)) for n, nome, tam in RE_GRUPO.findall(saida)}
    if not grupos:
        raise SystemExit(f"o make_ndx não listou grupos.\n{saida[-1500:]}")
    return grupos


def achar(grupos, *nomes):
    """Primeiro grupo cujo nome esteja entre `nomes` (sem diferenciar caixa)."""
    alvo = [n.lower() for n in nomes]
    for pref in alvo:                       # respeita a ordem pedida
        for num in sorted(grupos):
            if grupos[num][0].lower() == pref:
                return num, grupos[num][0]
    return None, None


def achar_ligante(grupos, n_protein: int):
    """O ligante é o grupo pequeno que não é proteína, água nem íon.

    Pega o menor grupo cujo nome não está na lista conhecida — que é como um
    resíduo batizado UNL, PTC, LIG ou MOL aparece, qualquer que seja o nome.
    Não achando nenhum, cai no `Other`, que num complexo proteína+ligante é
    exatamente o ligante.
    """
    candidatos = [(tam, num, nome) for num, (nome, tam) in grupos.items()
                  if nome.lower() not in NAO_LIGANTE and 0 < tam < n_protein]
    if candidatos:
        _, num, nome = min(candidatos)
        return num, nome
    num, nome = achar(grupos, "Other")
    if num is not None and 0 < grupos[num][1] < n_protein:
        return num, nome
    return None, None


def tamanhos_do_ndx(texto: str) -> dict[str, int]:
    """Nome -> quantidade de átomos, lendo o .ndx escrito."""
    tam = {}
    nome = None
    for linha in texto.splitlines():
        m = re.match(r"\s*\[\s*(.+?)\s*\]\s*$", linha)
        if m:
            nome = m.group(1)
            tam[nome] = 0
        elif nome is not None:
            tam[nome] += len(linha.split())
    return tam


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gro", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--gmx", default=GMX)
    ap.add_argument("--nome-solutos", default="Protein_LIG",
                    help="nome do grupo proteína+ligante (default: Protein_LIG)")
    ap.add_argument("--nome-solvente", default="Water_and_ions")
    args = ap.parse_args()

    gro = args.gro.expanduser()
    out = args.out.expanduser()
    if not gro.exists():
        raise SystemExit(f"não achei {gro}")

    grupos = listar_grupos(gro, args.gmx)
    print(f"  {len(grupos)} grupos no sistema")

    n_prot, nome_prot = achar(grupos, "Protein")
    if n_prot is None:
        raise SystemExit("não achei o grupo Protein — o .gro tem proteína?")
    tam_prot = grupos[n_prot][1]

    n_lig, nome_lig = achar_ligante(grupos, tam_prot)
    n_agua, _ = achar(grupos, "Water", "SOL")
    n_ion, nome_ion = achar(grupos, "Ion", "Ions", "CL", "NA")

    print(f"  proteína : {n_prot} ({nome_prot}, {tam_prot} átomos)")
    if n_lig is None:
        raise SystemExit(
            "não identifiquei o grupo do ligante entre: "
            + ", ".join(f"{n}={grupos[n][0]}" for n in sorted(grupos)))
    print(f"  ligante  : {n_lig} ({nome_lig}, {grupos[n_lig][1]} átomos)")
    if n_agua is None:
        raise SystemExit("não achei o grupo da água")
    print(f"  água     : {n_agua} ({grupos[n_agua][1]} átomos)")
    print(f"  íons     : "
          + (f"{n_ion} ({nome_ion}, {grupos[n_ion][1]} átomos)"
             if n_ion is not None else "nenhum"))

    # Monta o roteiro. Um grupo já existente com o nome pedido é reaproveitado:
    # num sistema neutralizado o GROMACS já cria Water_and_ions, e criar um
    # segundo com o mesmo nome deixaria o grompp escolhendo entre homônimos.
    proximo = max(grupos) + 1
    cmds, criados = [], []

    n_sol, _ = achar(grupos, args.nome_solutos)
    if n_sol is None:
        cmds += [f"{n_prot} | {n_lig}", f"name {proximo} {args.nome_solutos}"]
        criados.append(f"{args.nome_solutos} = {n_prot} | {n_lig}")
        proximo += 1
    else:
        print(f"  {args.nome_solutos} já existe (grupo {n_sol})")

    n_slv, _ = achar(grupos, args.nome_solvente)
    if n_slv is None:
        fusao = f"{n_agua} | {n_ion}" if n_ion is not None else f"{n_agua}"
        cmds += [fusao, f"name {proximo} {args.nome_solvente}"]
        criados.append(f"{args.nome_solvente} = {fusao}")
        proximo += 1
    else:
        print(f"  {args.nome_solvente} já existe (grupo {n_slv})")

    cmds.append("q")
    print("  fusões   : " + (" ; ".join(criados) if criados
                             else "nenhuma (só copiando os grupos padrão)"))

    r = rodar_make_ndx(gro, args.gmx, "\n".join(cmds) + "\n", str(out))
    saida = (r.stdout or "") + (r.stderr or "")
    texto = out.read_text() if out.exists() else ""

    faltam = [n for n in (args.nome_solutos, args.nome_solvente)
              if f"[ {n} ]" not in texto]
    if faltam:
        out.unlink(missing_ok=True)   # não deixa índice pela metade no disco
        raise SystemExit(f"grupos ausentes no índice: {faltam}\n{saida[-1200:]}")

    tam = tamanhos_do_ndx(texto)
    n_sistema = grupos[0][1] if 0 in grupos else None
    soma = tam.get(args.nome_solutos, 0) + tam.get(args.nome_solvente, 0)
    print(f"  {args.nome_solutos}: {tam.get(args.nome_solutos)} átomos | "
          f"{args.nome_solvente}: {tam.get(args.nome_solvente)} átomos")

    if n_sistema and soma != n_sistema:
        out.unlink(missing_ok=True)
        raise SystemExit(
            f"os dois grupos somam {soma} átomos mas o sistema tem "
            f"{n_sistema}: {n_sistema - soma} ficariam sem acoplamento "
            f"térmico e o grompp recusaria a simulação")

    print(f"  índice verificado: {out}")


if __name__ == "__main__":
    main()
