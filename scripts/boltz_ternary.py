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


def sequencias_por_cadeia(pdb: Path, cadeias=None):
    """Uma sequência POR CADEIA, nunca concatenadas.

    Concatenar cadeias distintas numa única entrada de proteína inventa uma
    ligação peptídica entre elas. A PCSK9 é o caso típico: o prodomínio
    permanece ligado NÃO covalentemente depois da clivagem autocatalítica, e
    emendá-lo à cadeia catalítica produziria uma molécula que não existe.
    """
    por_cadeia, faltantes, vistos = {}, set(), set()
    for l in pdb.read_text().splitlines():
        if not l.startswith(("ATOM", "HETATM")):
            continue
        if l[16] not in (" ", "A"):
            continue
        cad = l[21]
        if cadeias and cad not in cadeias:
            continue
        chave = (cad, l[22:27])
        if chave in vistos:
            continue
        vistos.add(chave)
        nome = l[17:20].strip().upper()
        if nome in TRES_PARA_UM:
            por_cadeia.setdefault(cad, []).append(TRES_PARA_UM[nome])
        elif nome in ("HOH", "WAT", "NA", "CL", "SOL", "ZN", "MG", "CA"):
            continue
        else:
            faltantes.add(nome)
    # cadeias curtas demais são peptídeos de cristalização, não a proteína
    return ({c: "".join(s) for c, s in por_cadeia.items() if len(s) >= 20},
            sorted(faltantes))


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
    ap.add_argument("--pcsk9-chain", default=None,
                    help="default: TODAS as cadeias do receptor preparado — "
                         "a PCSK9 madura é de duas cadeias e o docking manteve "
                         "as duas (PCSK9_KEEP_CHAINS)")
    ap.add_argument("--smiles", default=None)
    ap.add_argument("--job-dir", type=Path, default=None,
                    help="diretório do job do PRosettaC deste candidato. O "
                         "prosetta_config.txt dele diz as estruturas, as "
                         "cadeias e o SMILES — é a fonte autoritativa, e é "
                         "independente de track.")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    work = args.work.expanduser()
    out = (args.out or work / "boltz_ternario" / args.candidato).expanduser()
    out.mkdir(parents=True, exist_ok=True)

    # --- receptores, do CONFIG DO JOB -------------------------------------
    # A primeira versão disto tinha os caminhos da CRBN escritos no código
    # (`prep/CRBN_4TZ4/4TZ4_receptor.pdb`, `pipeline/md/receptor_fixed.pdb`) e
    # o SMILES vinha de `pipeline/`, que é a saída DAQUELE track. Num projeto
    # que abandonou a CRBN, isso não é só inútil: é um script que só sabe
    # trabalhar no track errado.
    #
    # A fonte autoritativa é o prosetta_config.txt do job: ele nomeia as
    # estruturas, as cadeias e o arquivo do PROTAC, foi o que o PRosettaC de
    # fato usou, e não sabe nem precisa saber de qual track veio.
    cfg_cadeias = {}
    if args.job_dir:
        jd = args.job_dir.expanduser()
        cfg = jd / "prosetta_config.txt"
        if not cfg.exists():
            raise SystemExit(f"não achei {cfg}")
        campos = {}
        for linha in cfg.read_text().splitlines():
            if ": " in linha:
                k, v = linha.split(": ", 1)
                campos[k.strip()] = v.split()
        estruturas = campos.get("Structures", [])
        cadeias = campos.get("Chains", [])
        heads = campos.get("Heads", [])
        if not (len(estruturas) == len(cadeias) == len(heads) == 2):
            raise SystemExit(
                f"o config não tem 2 estruturas/cadeias/heads:\n"
                f"  Structures: {estruturas}\n  Chains: {cadeias}\n"
                f"  Heads: {heads}")

        # QUAL LADO É A E3, decidido pelo CONTEÚDO e não pela posição. A
        # convenção do PRosettaC é que a estrutura i casa com o head i, e o
        # head da E3 é o recrutador. Procurar "recruiter" no nome é mais
        # verificável que confiar na ordem — e se não der para decidir, o
        # script para em vez de inverter E3 e alvo silenciosamente, o que
        # trocaria o significado de toda a comparação adiante.
        idx_e3 = [i for i, h in enumerate(heads) if "recruit" in h.lower()]
        if len(idx_e3) != 1:
            raise SystemExit(
                f"não consegui decidir qual head é o recrutador da E3: {heads}."
                f"\nPasse --e3-pdb/--e3-chain e --pcsk9-pdb/--pcsk9-chain.")
        i_e3 = idx_e3[0]
        i_pc = 1 - i_e3
        e3 = args.e3_pdb or (jd / estruturas[i_e3])
        pcsk9 = args.pcsk9_pdb or (jd / estruturas[i_pc])
        cfg_cadeias = {"e3": cadeias[i_e3], "pcsk9": cadeias[i_pc]}
        print(f"do config do job: E3 = {estruturas[i_e3]} cadeia "
              f"{cadeias[i_e3]} (head {heads[i_e3]})")
        print(f"                  alvo = {estruturas[i_pc]} cadeia "
              f"{cadeias[i_pc]} (head {heads[i_pc]})")
    else:
        # Caminho antigo, do track da CRBN. Mantido só para reproduzir um
        # número daquele track; para a VHL, use --job-dir.
        e3 = args.e3_pdb or achar(
            work / "prep" / "CRBN_4TZ4" / "4TZ4_receptor.pdb",
            work / "pipeline" / "md" / "receptor_fixed.pdb")
        pcsk9 = args.pcsk9_pdb or achar(
            Path.home() / "PCSK9_docking" / "receptor" / "6U26_receptor.pdb",
            Path.home() / "PCSK9_docking" / "receptor.pdb",
            Path.home() / "structures" / "6U26.pdb")
    if e3 is None or pcsk9 is None or not Path(e3).exists() \
            or not Path(pcsk9).exists():
        raise SystemExit(
            "não achei os PDB dos receptores.\n"
            f"  E3   : {e3}\n  PCSK9: {pcsk9}\n"
            "Passe --job-dir <dir do job do PRosettaC> — o config dele nomeia\n"
            "as estruturas — ou --e3-pdb e --pcsk9-pdb explicitamente.")

    # --- SMILES do PROTAC -------------------------------------------------
    smi = args.smiles
    if smi is None and args.job_dir:
        # O protac.smi do job é o SMILES que o PRosettaC de fato usou. Os
        # JSON/CSV abaixo vivem em `pipeline/`, a saída do track da CRBN.
        ps = args.job_dir.expanduser() / "protac.smi"
        if ps.exists() and ps.read_text().split():
            smi = ps.read_text().split()[0]
            print(f"SMILES: do protac.smi do job ({len(smi)} caracteres)")
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

    # --- sequências, uma entrada por cadeia -------------------------------
    cad_e3 = [args.e3_chain] if args.e3_chain else (
        [cfg_cadeias["e3"]] if cfg_cadeias else None)
    # A PCSK9 madura é de duas cadeias e o receptor preparado manteve as duas,
    # então aqui o default continua sendo TODAS — restringir à cadeia do
    # config daria ao Boltz-2 metade do alvo.
    cad_pc = [args.pcsk9_chain] if args.pcsk9_chain else None
    seqs_e3, falta_e3 = sequencias_por_cadeia(e3, cad_e3)
    seqs_pc, falta_pc = sequencias_por_cadeia(pcsk9, cad_pc)

    print(f"candidato: {args.candidato}")
    print(f"  E3    : {e3.name}")
    for c, s in sorted(seqs_e3.items()):
        print(f"            cadeia {c or '(sem id)'}: {len(s)} resíduos")
    print(f"  PCSK9 : {pcsk9.name}")
    for c, s in sorted(seqs_pc.items()):
        print(f"            cadeia {c or '(sem id)'}: {len(s)} resíduos")
    print(f"  PROTAC: {smi[:70]}{'...' if len(smi) > 70 else ''}")
    for rot, falta in (("E3", falta_e3), ("PCSK9", falta_pc)):
        if falta:
            print(f"  [nota] resíduos fora da sequência em {rot}: {falta}")
            print(f"         (heteroátomos e ligantes do cristal são esperados "
                  f"aqui; aminoácido nesta lista não é)")
    if not seqs_e3 or not seqs_pc:
        raise SystemExit(
            "não extraí sequência de um dos receptores — cadeia errada?\n"
            "Liste as cadeias com: grep '^ATOM' <pdb> | cut -c22 | sort -u")

    total = sum(len(s) for s in seqs_e3.values()) + \
        sum(len(s) for s in seqs_pc.values())
    print(f"\n  total: {total} resíduos em "
          f"{len(seqs_e3) + len(seqs_pc)} cadeias")
    if total > 1400:
        print(f"  [ATENÇÃO] sistema grande para o Boltz-2. Se estourar a "
              f"memória da GPU, rode com 1 amostra.")

    # --- YAML: uma entrada de proteína POR CADEIA -------------------------
    ids = "ABCDEFGHIJKMNOPQRSTUVWXYZ"          # L fica para o ligante
    linhas_yaml, mapa = ["version: 1", "sequences:"], {}
    i = 0
    for rotulo, seqs in (("E3", seqs_e3), ("PCSK9", seqs_pc)):
        for c, s in sorted(seqs.items()):
            linhas_yaml += ["  - protein:", f"      id: {ids[i]}",
                            f'      sequence: "{s}"']
            mapa[ids[i]] = f"{rotulo} cadeia {c or '?'} ({len(s)} res)"
            i += 1
    linhas_yaml += ["  - ligand:", "      id: L", f'      smiles: "{smi}"']
    mapa["L"] = "PROTAC"

    yaml = out / "ternary.yaml"
    yaml.write_text("\n".join(linhas_yaml) + "\n")

    (out / "entradas.json").write_text(json.dumps({
        "candidate_id": args.candidato,
        "e3_pdb": str(e3),
        "e3_chains": {c: len(s) for c, s in seqs_e3.items()},
        "pcsk9_pdb": str(pcsk9),
        "pcsk9_chains": {c: len(s) for c, s in seqs_pc.items()},
        "boltz_ids": mapa,
        "total_residues": total,
        "protac_smiles": smi,
    }, indent=2))

    print(f"\n  {yaml}")
    for k, v in mapa.items():
        print(f"    id {k} = {v}")
    print(f"\nPróximo: bash scripts/run_boltz_ternary.sh {args.candidato}")


if __name__ == "__main__":
    main()
