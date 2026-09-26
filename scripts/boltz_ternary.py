#!/usr/bin/env python
"""
Monta a entrada do Boltz-2 para o complexo ternário PCSK9 + PROTAC + CRBN.

Só biblioteca padrão. Segundos.

    python boltz_ternary.py --candidato SC0006__WH023

O Boltz-2 recebe SEQUÊNCIA, não estrutura. Isso tem uma consequência boa: ele
modela a CRBN com as 4 alças que o cristal 4TZ4 não resolve — as mesmas que
tivemos de cortar para a MD de nível (ii). O modelo ternário sai sem os 8
términos artificiais, e o movimento de domínio que apareceu lá não se repete
por essa causa.

As sequências saem dos PDB que o pipeline já usou, não de um FASTA à parte:
assim o que se prediz é o mesmo objeto que se docou.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

TRES_PARA_UM = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
    # variantes de protonação que o pdb2gmx/ff14SB escrevem
    "HID": "H", "HIE": "H", "HIP": "H", "HISD": "H", "HISE": "H", "HISP": "H",
    "CYX": "C", "CYM": "C", "ASH": "D", "GLH": "E", "LYN": "K",
    "MSE": "M",
}


def sequencia(pdb: Path, cadeia: str | None = None):
    """Sequência de uma letra, na ordem do arquivo, sem depender de bibliotecas.

    Lê resíduos por (cadeia, resSeq, iCode) para não repetir o mesmo resíduo a
    cada átomo, e ignora altLoc secundário.
    """
    vistos, seq, faltantes = set(), [], set()
    for l in pdb.read_text().splitlines():
        if not l.startswith(("ATOM", "HETATM")):
            continue
        if l[16] not in (" ", "A"):          # altLoc secundário
            continue
        cad = l[21]
        if cadeia and cad != cadeia:
            continue
        chave = (cad, l[22:27])
        if chave in vistos:
            continue
        vistos.add(chave)
        nome = l[17:20].strip().upper()
        if nome in TRES_PARA_UM:
            seq.append(TRES_PARA_UM[nome])
        elif nome in ("HOH", "WAT", "NA", "CL", "SOL", "ZN", "MG"):
            continue
        else:
            faltantes.add(nome)
    return "".join(seq), sorted(faltantes)


def achar(*candidatos):
    for c in candidatos:
        if c and Path(c).expanduser().exists():
            return Path(c).expanduser()
    return None


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidato", required=True)
    ap.add_argument("--work", type=Path,
                    default=Path.home() / "PRosettaC_runs" / "vhl_crbn_pcsk9_protac")
    ap.add_argument("--e3-pdb", type=Path, default=None)
    ap.add_argument("--e3-chain", default=None)
    ap.add_argument("--pcsk9-pdb", type=Path, default=None)
    ap.add_argument("--pcsk9-chain", default="B")
    ap.add_argument("--smiles", default=None)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    work = args.work.expanduser()
    out = (args.out or work / "boltz_ternario" / args.candidato).expanduser()
    out.mkdir(parents=True, exist_ok=True)

    # --- receptores -------------------------------------------------------
    e3 = args.e3_pdb or achar(
        work / "prep" / "CRBN_4TZ4" / "4TZ4_receptor.pdb",
        work / "pipeline" / "md" / "receptor_fixed.pdb")
    pcsk9 = args.pcsk9_pdb or achar(
        Path.home() / "PCSK9_docking" / "receptor.pdb",
        Path.home() / "PCSK9_docking" / "6U26_receptor.pdb",
        Path.home() / "structures" / "6U26.pdb")
    if e3 is None or pcsk9 is None:
        raise SystemExit(
            "não achei os PDB dos receptores.\n"
            f"  E3   : {e3}\n  PCSK9: {pcsk9}\n"
            "Passe --e3-pdb e --pcsk9-pdb explicitamente.")

    # --- SMILES do PROTAC -------------------------------------------------
    smi = args.smiles
    if smi is None:
        for p in (work / "pipeline" / "md" / "md_sistema.json",
                  work / "pipeline" / "md_candidato.json"):
            if p.exists():
                d = json.loads(p.read_text())
                if d.get("protac_smiles") and \
                        d.get("candidate_id") == args.candidato:
                    smi = d["protac_smiles"]
                    break
    if smi is None:
        csv = work / "pipeline" / "protac_candidates.csv"
        if csv.exists():
            import csv as _csv
            for row in _csv.DictReader(open(csv)):
                if row.get("candidate_id") == args.candidato:
                    smi = row.get("protac_smiles")
                    break
    if not smi:
        raise SystemExit(f"não achei o SMILES de {args.candidato}. "
                         f"Passe --smiles.")

    # --- sequências -------------------------------------------------------
    seq_e3, falta_e3 = sequencia(e3, args.e3_chain)
    seq_pc, falta_pc = sequencia(pcsk9, args.pcsk9_chain)

    print(f"candidato: {args.candidato}")
    print(f"  E3    : {e3.name}"
          + (f" cadeia {args.e3_chain}" if args.e3_chain else " (todas)")
          + f" -> {len(seq_e3)} resíduos")
    print(f"  PCSK9 : {pcsk9.name} cadeia {args.pcsk9_chain} "
          f"-> {len(seq_pc)} resíduos")
    print(f"  PROTAC: {smi[:70]}{'...' if len(smi) > 70 else ''}")
    for rot, falta in (("E3", falta_e3), ("PCSK9", falta_pc)):
        if falta:
            print(f"  [ATENÇÃO] resíduos não reconhecidos em {rot}: {falta}")
            print(f"            eles NÃO entraram na sequência — confira se "
                  f"são heteroátomos (ok) ou aminoácidos (não ok)")
    if len(seq_e3) < 50 or len(seq_pc) < 50:
        raise SystemExit(
            "sequência curta demais para ser uma proteína — cadeia errada?\n"
            "Liste as cadeias com: grep '^ATOM' <pdb> | cut -c22 | sort -u")

    # --- YAML -------------------------------------------------------------
    yaml = out / "ternary.yaml"
    yaml.write_text(
        "version: 1\n"
        "sequences:\n"
        "  - protein:\n"
        "      id: A\n"
        f"      sequence: \"{seq_e3}\"\n"
        "  - protein:\n"
        "      id: B\n"
        f"      sequence: \"{seq_pc}\"\n"
        "  - ligand:\n"
        "      id: L\n"
        f"      smiles: \"{smi}\"\n")

    (out / "entradas.json").write_text(json.dumps({
        "candidate_id": args.candidato,
        "e3_pdb": str(e3), "e3_chain": args.e3_chain,
        "e3_residues": len(seq_e3),
        "pcsk9_pdb": str(pcsk9), "pcsk9_chain": args.pcsk9_chain,
        "pcsk9_residues": len(seq_pc),
        "protac_smiles": smi,
    }, indent=2))

    print(f"\n  {yaml}")
    print(f"  cadeia A = E3 (CRBN) | cadeia B = PCSK9 | L = PROTAC")
    print(f"\nPróximo: bash scripts/run_boltz_ternary.sh {args.candidato}")


if __name__ == "__main__":
    main()
