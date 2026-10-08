#!/usr/bin/env python
"""
O recrutador é quimicamente compatível com a E3 que o projeto declara?

Roda no env `mdtools` (RDKit). Instantâneo.

    python checar_recrutador.py --e3 VHL --sdf recruiter_VHL_ligand_087.sdf
    python checar_recrutador.py --e3 VHL --smiles "O=C1CC..."
    python checar_recrutador.py --e3 VHL --protac-smi protac.smi

Por que esta verificação não existia, e o que custou
---------------------------------------------------
O WP1 escolheu o recrutador por SCORE DE DOCKING no sítio da E3, mais
enterramento. Nenhuma das duas medidas sabe que E3 a molécula de fato recruta:
uma biblioteca de blocos de PROTAC (Chemspace) contém ligantes de CRBN
(análogos de talidomida) E de VHL, e um análogo de talidomida pode pontuar bem
encaixado no bolso da VHL sem ter razão nenhuma para se ligar nela.

Foi o que aconteceu. O `protac.smi` do candidato SC0013__WH023 — o PROTAC que
consumiu sete dias de PRosettaC no "track da VHL" — carrega glutarimida +
ftalimida fluorada, ou seja, uma talidomida: o recrutador canônico da CRBN.
Nenhuma das três assinaturas de ligante de VHL está presente. A corrida inteira
dockeou a proteína VHL contra um ligante de CRBN.

E isso explica o resultado sem precisar de mais nada: um "recrutador" que não
se liga à E3 não produz interface cooperativa, e 103 clusters com mediana de
41 Å é o que a geometria devolve quando a única coisa que prende as duas
proteínas é a restrição de distância do linker.

Calibração das evidências
-------------------------
A identificação POSITIVA é a forte: glutarimida ligada a ftalimida é
talidomida, e talidomida liga CRBN. Isso basta.

A ausência das assinaturas de VHL é evidência forte mas não absoluta: a
esmagadora maioria dos ligantes de VHL deriva do VH032 e tem
(2S,4R)-4-hidroxiprolina, mas existem exceções na literatura. Por isso o script
reporta os dois lados em vez de decidir por ausência.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

# SMARTS por E3. A chave é a E3 que a molécula RECRUTA.
ASSINATURAS = {
    "CRBN": [
        ("glutarimida",   "O=C1CCC(C(=O)N1)N",  True),
        ("ftalimida",     "O=C1c2ccccc2C(=O)N1", False),
        ("isoindolinona", "O=C1Cc2ccccc2N1",     False),
    ],
    "VHL": [
        ("4-hidroxiprolina",     "OC1CC(NC1)C(=O)N",            True),
        ("terc-leucina amida",   "CC(C)(C)C(NC(=O)C)C(=O)N",    False),
        ("tiazol 4-metil",       "Cc1scnc1",                    False),
        ("benzilamida terminal", "c1ccccc1CNC(=O)",             False),
    ],
}
# o terceiro campo diz se a assinatura é SUFICIENTE sozinha para identificar


def motivos(mol):
    achados = {}
    for e3, lista in ASSINATURAS.items():
        achados[e3] = [(nome, suf) for nome, sma, suf in lista
                       if (q := Chem.MolFromSmarts(sma)) is not None
                       and mol.HasSubstructMatch(q)]
    return achados


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--e3", required=True, help="a E3 que o projeto declara")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--sdf", type=Path)
    g.add_argument("--smiles")
    g.add_argument("--protac-smi", type=Path,
                   help="o protac.smi do job: confere a molécula inteira")
    args = ap.parse_args()

    if args.sdf:
        mol = next(iter(Chem.SDMolSupplier(str(args.sdf.expanduser()),
                                           removeHs=True)), None)
        fonte = args.sdf.name
    elif args.protac_smi:
        txt = args.protac_smi.expanduser().read_text().split()
        mol = Chem.MolFromSmiles(txt[0]) if txt else None
        fonte = args.protac_smi.name
    else:
        mol = Chem.MolFromSmiles(args.smiles)
        fonte = "--smiles"
    if mol is None:
        raise SystemExit(f"não consegui ler a molécula de {fonte}")

    e3 = args.e3.upper()
    print(f"molécula: {fonte} — {mol.GetNumAtoms()} átomos pesados")
    print(f"E3 declarada pelo projeto: {e3}")
    ach = motivos(mol)
    print()
    for alvo, lista in ach.items():
        if lista:
            print(f"  assinaturas de {alvo}: "
                  + ", ".join(n + (" [suficiente]" if s else "")
                              for n, s in lista))
        else:
            print(f"  assinaturas de {alvo}: NENHUMA")

    identificadas = {a for a, l in ach.items() if any(s for _, s in l)}
    print()
    if e3 in identificadas and len(identificadas) == 1:
        print(f"  OK: a molécula tem assinatura suficiente de {e3} e de mais"
              f" nenhuma E3.")
        return 0
    outras = identificadas - {e3}
    if outras:
        print(f"  *** INCOMPATÍVEL: assinatura suficiente de"
              f" {', '.join(sorted(outras))},")
        print(f"  *** e o projeto declara {e3}.")
        print(f"  ***")
        print(f"  *** Score de docking não sabe qual E3 a molécula recruta. Uma")
        print(f"  *** biblioteca de blocos de PROTAC tem ligantes das duas, e um")
        print(f"  *** análogo de talidomida encaixa no bolso da VHL sem ter razão")
        print(f"  *** para se ligar nela. Seguir daqui modela um complexo que a")
        print(f"  *** química não sustenta.")
        return 2
    print(f"  *** SEM ASSINATURA de {e3}, e sem assinatura de outra E3.")
    print(f"  *** Isto não prova que a molécula não se liga à {e3} — há")
    print(f"  *** ligantes fora dos quimiotipos canônicos — mas exige")
    print(f"  *** justificativa explícita antes de gastar dias de máquina.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
