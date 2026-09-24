#!/usr/bin/env python
"""
Fase 8a — prepara o sistema de MD do PROTAC escolhido.

Roda no env `mdtools`. Minutos (o antechamber/AM1-BCC é o passo longo).

Monta o nível (ii) da metodologia: **E3 ligase + recrutador–linker–warhead,
sem a PCSK9**. É o nível que detecta interação espúria entre o warhead e a
própria E3 ligase — uma orientação improdutiva do linker que condena o
candidato antes de gastar dias no ternário completo.

O nível (iii), ternário completo, exige a estrutura do PRosettaC ou do
AlphaFold 3 e por isso não entra aqui.

Como as coordenadas são obtidas
--------------------------------
O sub-complexo do WP2 é montado por substituição de subestrutura e não tem
geometria 3D utilizável. Aqui o PROTAC é embebido com o **recrutador fixo nas
coordenadas do docking** (ConstrainedEmbed), de modo que a parte ancorada na
E3 ligase fica onde o docking a colocou e o resto se acomoda a partir dali.

    python md_prepare.py \\
        --candidato ~/pipeline/md_candidato.json \\
        --recruiter ~/pipeline/recruiter_CRBN_ligand_116.sdf \\
        --receptor ~/.../4TZ4_receptor.pdb \\
        --outdir ~/pipeline/md
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors

RDLogger.DisableLog("rdApp.*")

LIG_RESNAME = "PTC"          # 3 letras, como o GROMACS espera


def embutir_com_recrutador_fixo(protac_smiles: str, recruiter,
                                n_tentativas: int = 20, seed: int = 0xC0FFEE):
    """PROTAC em 3D com o recrutador nas coordenadas do docking."""
    protac = Chem.MolFromSmiles(protac_smiles)
    if protac is None:
        raise SystemExit("SMILES do PROTAC inválido")
    protac = Chem.AddHs(protac)

    core = Chem.RemoveHs(Chem.Mol(recruiter))
    match = protac.GetSubstructMatch(core)
    if not match:
        # o recrutador pode ter perdido um H ao se conjugar: tenta sem os
        # átomos terminais do core
        core = Chem.DeleteSubstructs(core, Chem.MolFromSmarts("[#1]"))
        match = protac.GetSubstructMatch(core)
    if not match:
        raise SystemExit(
            "o recrutador não casa como subestrutura do PROTAC — a montagem do "
            "WP3 alterou o recrutador? Confira o SDF passado em --recruiter.")

    mapa = {idx_protac: core.GetConformer().GetAtomPosition(i)
            for i, idx_protac in enumerate(match)}

    melhor, melhor_e = None, np.inf
    params = AllChem.ETKDGv3()
    for tentativa in range(n_tentativas):
        params.randomSeed = seed + tentativa
        m = Chem.Mol(protac)
        if AllChem.EmbedMolecule(m, coordMap=mapa, randomSeed=seed + tentativa,
                                 useRandomCoords=True) != 0:
            continue
        ff = AllChem.MMFFGetMoleculeForceField(
            m, AllChem.MMFFGetMoleculeProperties(m))
        if ff is None:
            continue
        for idx in match:                      # trava o recrutador
            ff.AddFixedPoint(idx)
        ff.Minimize(maxIts=2000)
        e = ff.CalcEnergy()
        if e < melhor_e:
            melhor_e, melhor = e, m

    if melhor is None:
        raise SystemExit("não consegui embeber o PROTAC com o recrutador fixo")

    # conferência: o recrutador ficou onde o docking o colocou?
    pos = melhor.GetConformer().GetPositions()
    desvio = float(np.sqrt(np.mean([
        np.sum((pos[idx] - np.array(core.GetConformer().GetAtomPosition(i))) ** 2)
        for i, idx in enumerate(match)])))
    return melhor, melhor_e, desvio, len(match)


def nomear_residuo(mol, resname: str):
    """Batiza o resíduo do ligante no PDB.

    Sem isto o RDKit escreve `UNL`, o acpype propaga `UNL` para o .gro, e o
    nome viaja até o fim do pipeline: o make_ndx cria um grupo `UNL` em vez do
    esperado, e a análise procura `resname PTC` na trajetória e não acha nada.
    Batizar aqui, na única vez que as coordenadas são escritas, resolve nos
    dois lugares.
    """
    for i, atom in enumerate(mol.GetAtoms()):
        info = Chem.AtomPDBResidueInfo()
        info.SetResidueName(f"{resname:>3s}")
        info.SetResidueNumber(1)
        info.SetChainId("X")
        info.SetName(f" {atom.GetSymbol()}{i}"[:4].ljust(4))
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)
    return mol


def escrever_complexo(receptor_pdb: Path, ligante, out_pdb: Path,
                      resname: str = LIG_RESNAME):
    """Receptor + ligante num único PDB, o ligante como HETATM."""
    linhas = [l for l in Path(receptor_pdb).read_text().splitlines()
              if l.startswith(("ATOM", "TER"))]
    serial = 90000
    conf = ligante.GetConformer()
    for i, atom in enumerate(ligante.GetAtoms()):
        p = conf.GetAtomPosition(i)
        el = atom.GetSymbol()
        nome = f"{el}{i}"[:4]
        linhas.append(
            f"HETATM{serial:5d} {nome:<4s}{resname:>3s} X 900    "
            f"{p.x:8.3f}{p.y:8.3f}{p.z:8.3f}  1.00  0.00          {el:>2s}")
        serial += 1
    linhas.append("END")
    out_pdb.write_text("\n".join(linhas) + "\n")
    return out_pdb


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidato", type=Path, required=True)
    ap.add_argument("--recruiter", type=Path, required=True)
    ap.add_argument("--receptor", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--rank", type=int, default=1,
                    help="1 = melhor; use 2, 3... para os próximos da fila "
                         "quando a MD reprovar o primeiro")
    ap.add_argument("--n-tentativas", type=int, default=20)
    args = ap.parse_args()

    out = args.outdir.expanduser()
    out.mkdir(parents=True, exist_ok=True)

    dados = json.loads(args.candidato.expanduser().read_text())
    cand = (dados["escolhido"] if args.rank == 1
            else dados["proximos"][args.rank - 2])
    print(f"[1/4] Candidato rank {args.rank}: {cand['candidate_id']}")
    print(f"      MW {cand['protac_mw']} Da, {cand['protac_rotb']} torções")

    recruiter = Chem.MolFromMolFile(str(args.recruiter.expanduser()))
    if recruiter is None:
        raise SystemExit(f"não consegui ler {args.recruiter}")

    print(f"[2/4] Embebendo o PROTAC com o recrutador fixo "
          f"({args.n_tentativas} tentativas)")
    mol, energia, desvio, n_fixos = embutir_com_recrutador_fixo(
        cand["protac_smiles"], recruiter, args.n_tentativas)
    print(f"      {n_fixos} átomos do recrutador travados nas coordenadas "
          f"do docking (desvio {desvio:.3f} Å)")
    print(f"      energia MMFF {energia:.1f} kcal/mol")
    if desvio > 0.5:
        print(f"      [ATENÇÃO] desvio acima de 0,5 Å: o recrutador saiu da "
              f"pose do docking, e a MD partiria de uma geometria que a "
              f"triagem não validou")

    lig_sdf = out / "protac.sdf"
    with Chem.SDWriter(str(lig_sdf)) as w:
        w.write(mol)
    lig_pdb = out / "protac.pdb"
    nomear_residuo(mol, LIG_RESNAME)
    Chem.MolToPDBFile(mol, str(lig_pdb))

    print("[3/4] Escrevendo o complexo E3 + PROTAC")
    complexo = escrever_complexo(args.receptor.expanduser(), mol,
                                 out / "complexo.pdb")
    print(f"      {complexo}")

    carga = Chem.GetFormalCharge(Chem.RemoveHs(mol))
    info = {
        "candidate_id": cand["candidate_id"],
        "rank": args.rank,
        "protac_smiles": cand["protac_smiles"],
        "protac_mw": round(Descriptors.MolWt(mol), 1),
        "carga_formal": carga,
        "n_atomos_recrutador_fixos": n_fixos,
        "desvio_do_docking_A": round(desvio, 3),
        "energia_mmff": round(energia, 1),
        "lig_sdf": str(lig_sdf), "lig_pdb": str(lig_pdb),
        "complexo_pdb": str(complexo),
        "receptor_pdb": str(args.receptor.expanduser()),
        "resname": LIG_RESNAME,
    }
    (out / "md_sistema.json").write_text(json.dumps(info, indent=2))

    print(f"[4/4] Carga formal do PROTAC: {carga:+d}")
    print(f"      É ela que vai em `acpype -n`. Carga errada aqui contamina")
    print(f"      a trajetória inteira sem dar erro.")
    print(f"\n{out / 'md_sistema.json'}")


if __name__ == "__main__":
    main()
