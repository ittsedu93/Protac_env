#!/usr/bin/env python
"""
Conserta os .sdf dos heads para o PRosettaC conseguir casá-los.

Roda no env `mdtools` (RDKit). Segundos.

    python prosettac_fix_heads.py --dir <run_dir>            # diagnostica
    python prosettac_fix_heads.py --dir <run_dir> --aplicar   # reescreve

O problema que isto resolve
---------------------------
O `translate_anchors` do PRosettaC faz:

    OldSdf = Chem.SDMolSupplier(old_sdf, sanitize=False)[0]
    try:    Chem.SanitizeMol(OldSdf)
    except: pass                      # <- a falha é engolida
    NewMatch = NewSdf.GetSubstructMatch(OldSdf)
    if len(NewMatch) == 0: return -1

Falhando a sanitização, o mol fica com átomos sem contagem de hidrogênio
implícito — o RDKit os escreve entre colchetes:

    Oc1[c]c([C][C]Nc2nnn[nH]2)[c]c(Oc2[c][c][c][c][c]2)[c]1

Um átomo de consulta nesse estado nunca casa com um carbono aromático que tem
H. O match sai vazio, a função devolve -1, e o -1 só aparece cinquenta linhas
adiante como índice negativo num GetAtomPosition.

Este script reproduz exatamente o que o translate_anchors faz — mesma leitura,
mesmo casamento — e, falhando, reescreve o .sdf a partir de uma versão que
sanitiza, preservando as coordenadas. Depois reconfere.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def como_o_prosettac_le(caminho: Path):
    """A mesma leitura do translate_anchors, e o mesmo except silencioso."""
    mol = Chem.SDMolSupplier(str(caminho), sanitize=False)[0]
    if mol is None:
        return None, "o RDKit não leu o arquivo nem sem sanitizar"
    erro = None
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        erro = str(e).split("\n")[0][:110]
    return mol, erro


def casamento_do_prosettac(caminho: Path):
    """O teste que importa: o GetSubstructMatch do PRosettaC casa ou não?

    Procurar `[c]` no SMILES era heurística, e errada — um tetrazol sem o
    [nH] sai como `nnnn`, sem colchete nenhum, e quebra igual. O critério é
    fazer o próprio casamento que a ferramenta faz.
    """
    velho, erro = como_o_prosettac_le(caminho)
    novo = Chem.SDMolSupplier(str(caminho), sanitize=True)[0]
    if velho is None:
        return 0, erro, None
    if novo is None:
        # nem sanitizado o RDKit lê: o match contra o _H também falharia
        return 0, erro or "não lê sanitizado", Chem.MolToSmiles(velho)
    match = novo.GetSubstructMatch(velho)
    return len(match), erro, Chem.MolToSmiles(velho)


def consertar(caminho: Path):
    """Reescreve o .sdf de modo que ele sanitize, preservando coordenadas.

    O caso real é o tetrazol da warhead: as ligações vêm marcadas aromáticas,
    nenhum N carrega hidrogênio explícito, e o RDKit não consegue kekulizar —
    ele não tem como adivinhar qual dos quatro N é o NH. A reconstrução testa
    cada N candidato e fica com o primeiro que sanitiza. É determinístico e
    não inventa: ou existe um tautômero válido, ou o arquivo está errado de
    outro jeito.
    """
    bruto = Chem.SDMolSupplier(str(caminho), sanitize=False)[0]
    if bruto is None:
        return None, "ilegível mesmo sem sanitizar"

    # 1) o caminho simples: já sanitiza?
    mol = Chem.Mol(bruto)
    try:
        Chem.SanitizeMol(mol)
        return mol, None
    except Exception as e:
        primeiro_erro = str(e).split("\n")[0][:110]

    # 2) candidatos a NH: nitrogênios de anel, com duas ligações e sem H.
    #
    # NÃO dá para exigir GetIsAromatic() aqui: num mol lido sem sanitizar o
    # RDKit ainda não percebeu aromaticidade nos ÁTOMOS (só nas ligações, a
    # partir do tipo 4 do mol block). Exigir o flag descartava todos os
    # candidatos e a correção nunca acontecia. O que é observável sem
    # sanitizar: símbolo, grau, H explícitos e pertencer a um anel.
    bruto.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(bruto)
    candidatos = [a.GetIdx() for a in bruto.GetAtoms()
                  if a.GetSymbol() == "N" and a.GetDegree() == 2
                  and a.GetNumExplicitHs() == 0 and a.IsInRing()]
    for idx in candidatos:
        mol = Chem.Mol(bruto)
        at = mol.GetAtomWithIdx(idx)
        at.SetNumExplicitHs(1)
        at.SetNoImplicit(True)
        try:
            Chem.SanitizeMol(mol)
            return mol, f"NH atribuído ao N {idx} (de {len(candidatos)} candidatos)"
        except Exception:
            continue

    return None, primeiro_erro


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
    heads = []
    for l in cfg.read_text().splitlines():
        if l.startswith("Heads: "):
            heads = [d / h for h in l.split(": ", 1)[1].split()]
    if not heads:
        raise SystemExit(f"não achei a linha Heads: em {cfg}")

    problemas = 0
    for h in heads:
        print(f"\n{h.name}")
        if not h.exists():
            print("  [ERRO] não existe")
            problemas += 1
            continue
        n, erro, smi = casamento_do_prosettac(h)
        if smi:
            print(f"  como o PRosettaC lê: {smi[:72]}")
        if erro:
            print(f"  sanitização FALHOU: {erro}")
        if n > 0:
            print(f"  OK — o GetSubstructMatch casa {n} átomos")
            continue

        print("  [PROBLEMA] o GetSubstructMatch sai VAZIO.")
        print("             O translate_anchors devolve -1, e o -1 vira o")
        print("             OverflowError cinquenta linhas adiante.")
        problemas += 1
        if not args.aplicar:
            continue
        novo, nota = consertar(h)
        if novo is None:
            print(f"  [ERRO] não consegui consertar: {nota}")
            continue
        shutil.copy2(h, h.with_suffix(".sdf.bak"))
        with Chem.SDWriter(str(h)) as w:
            w.write(novo)
        n2, erro2, smi2 = casamento_do_prosettac(h)
        if n2 > 0:
            print(f"  CORRIGIDO: {smi2[:72]}")
            if nota:
                print(f"             ({nota})")
            print(f"  match de {n2} átomos | original em "
                  f"{h.with_suffix('.sdf.bak').name}, coordenadas preservadas")
            problemas -= 1
        else:
            print(f"  [ERRO] ainda vazio depois de reescrever: {smi2}")

    print()
    if problemas and not args.aplicar:
        print(f"  {problemas} head(s) com problema. Para corrigir:")
        print(f"    python {Path(__file__).name} --dir {d} --aplicar")
        raise SystemExit(1)
    if problemas:
        raise SystemExit(f"  {problemas} head(s) seguem com problema")
    print("  todos os heads em ordem para o PRosettaC")


if __name__ == "__main__":
    main()
