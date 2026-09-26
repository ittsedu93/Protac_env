#!/usr/bin/env python
"""
Descobre os `Anchor atoms` corretos, a partir dos arquivos que o PRosettaC lê.

Roda no env `mdtools` (RDKit). Segundos.

    python prosettac_anchors.py --dir <run_dir>          # só informa
    python prosettac_anchors.py --dir <run_dir> --aplicar # corrige o config

O problema que isto resolve
---------------------------
O PRosettaC morre assim quando a âncora está errada:

    anchor_b = HeadB.GetConformer().GetAtomPosition(Anchors[1])
    OverflowError: can't convert negative value to unsigned int

O main.py lê a âncora 1-based e converte para 0-based (`int(i) - 1`), depois
manda o `translate_anchors` remapear o índice para a versão reescrita do head.
Falhando o remapeamento, ele devolve -1 — e o -1 só aparece lá adiante, como
um índice negativo num `GetAtomPosition`. A mensagem não diz "âncora errada".

A âncora é o átomo do HEAD que se liga ao LINKER. Isso é calculável: casando
cada head contra o PROTAC completo, o átomo do head cuja contraparte no PROTAC
tem vizinho FORA do head é o ponto de conjugação. Nada de tentar 1, depois 2.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import rdFMCS

RDLogger.DisableLog("rdApp.*")


def ler_config(cfg: Path):
    d = {}
    for l in cfg.read_text().splitlines():
        if ": " in l:
            k, v = l.split(": ", 1)
            d[k.strip()] = v.strip()
    return d


def ancora_do_head(head: Chem.Mol, protac: Chem.Mol):
    """Índice 0-based do átomo do head que se liga ao linker.

    Devolve (indice, quantos_candidatos, diagnostico).
    """
    head = Chem.RemoveHs(head)
    protac = Chem.RemoveHs(protac)

    mcs = rdFMCS.FindMCS([protac, head], timeout=30,
                         ringMatchesRingOnly=True, completeRingsOnly=False,
                         atomCompare=rdFMCS.AtomCompare.CompareElements,
                         bondCompare=rdFMCS.BondCompare.CompareOrderExact)
    if mcs.numAtoms == 0:
        return None, 0, "nenhuma subestrutura comum entre head e PROTAC"
    patt = Chem.MolFromSmarts(mcs.smartsString)
    m_pro = protac.GetSubstructMatch(patt)
    m_head = head.GetSubstructMatch(patt)
    if not m_pro or not m_head:
        return None, 0, "o MCS não casou nos dois"

    # head -> protac, pelos índices correspondentes do MCS
    para_protac = dict(zip(m_head, m_pro))
    dentro = set(m_pro)

    candidatos = []
    for i_head, i_pro in para_protac.items():
        a = protac.GetAtomWithIdx(i_pro)
        fora = [v.GetIdx() for v in a.GetNeighbors() if v.GetIdx() not in dentro]
        if fora:
            candidatos.append((i_head, head.GetAtomWithIdx(i_head).GetSymbol(),
                               len(fora)))

    if not candidatos:
        return None, 0, (f"o head casou {mcs.numAtoms} átomos, mas nenhum tem "
                         f"vizinho fora dele — o head parece ser o PROTAC inteiro")
    # havendo mais de um, o que tem mais ligações para fora é o ponto de saída
    candidatos.sort(key=lambda c: -c[2])
    return candidatos[0][0], len(candidatos), \
        f"MCS de {mcs.numAtoms} átomos; candidatos: " + \
        ", ".join(f"{i}({s})" for i, s, _ in candidatos)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--config", default="prosetta_config.txt")
    ap.add_argument("--aplicar", action="store_true",
                    help="reescreve a linha 'Anchor atoms' do config")
    args = ap.parse_args()

    d = args.dir.expanduser().resolve()
    cfg = d / args.config
    p = ler_config(cfg)

    heads = [d / h for h in p["Heads"].split()]
    smi_path = d / (p.get("Protac") or p.get("Linkers"))
    if not smi_path.exists():
        raise SystemExit(f"não achei o SMILES do PROTAC: {smi_path}")
    smi = smi_path.read_text().split()[0]
    protac = Chem.MolFromSmiles(smi)
    if protac is None:
        raise SystemExit(f"SMILES do PROTAC inválido: {smi[:60]}")

    print(f"PROTAC: {Chem.MolToSmiles(protac)[:70]}...")
    print(f"        {protac.GetNumAtoms()} átomos pesados\n")

    atuais = p.get("Anchor atoms", "").split()
    novos = []
    ok = True
    for i, h in enumerate(heads):
        mol = Chem.SDMolSupplier(str(h))[0] if h.exists() else None
        if mol is None:
            print(f"  [ERRO] não consegui ler {h.name}")
            ok = False
            novos.append(atuais[i] if i < len(atuais) else "1")
            continue
        idx, n_cand, diag = ancora_do_head(mol, protac)
        atual = atuais[i] if i < len(atuais) else "?"
        print(f"  {h.name}  ({Chem.RemoveHs(mol).GetNumAtoms()} átomos pesados)")
        print(f"    {diag}")
        if idx is None:
            print(f"    [ERRO] não achei o ponto de conjugação")
            ok = False
            novos.append(atual)
            continue
        um_based = idx + 1
        simbolo = Chem.RemoveHs(mol).GetAtomWithIdx(idx).GetSymbol()
        marca = "  <- igual ao config" if atual == str(um_based) else \
                f"  <- config diz {atual}, DIFERENTE"
        print(f"    âncora = átomo {um_based} (1-based), elemento {simbolo}"
              f"{marca}")
        if n_cand > 1:
            print(f"    [ATENÇÃO] {n_cand} átomos do head tocam o linker; "
                  f"escolhi o de mais ligações para fora")
        novos.append(str(um_based))
        print()

    linha = "Anchor atoms: " + " ".join(novos)
    print(f"  linha correta:  {linha}")

    if not ok:
        raise SystemExit("\n  não aplico com erro acima — resolva antes.")

    if args.aplicar:
        texto = re.sub(r"^Anchor atoms: .*$", linha, cfg.read_text(),
                       flags=re.M)
        cfg.write_text(texto)
        print(f"  config atualizado: {cfg}")
    else:
        print(f"\n  para aplicar:  python {Path(__file__).name} "
              f"--dir {d} --aplicar")


if __name__ == "__main__":
    main()
