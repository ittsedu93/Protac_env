#!/usr/bin/env python
"""
O vetor de saída de UM recrutador, sem passar pela triagem do WP1.

Roda no env `mdtools` (RDKit). Instantâneo.

    python exit_vector_de.py --ligante 6GFZ_ref_ligand.sdf \\
                             --receptor 6GFZ_receptor.pdb

Para que serve
--------------
A triagem do WP1 escolhe um recrutador entre 159 e extrai o vetor de saída do
vencedor. Quando o recrutador já está DECIDIDO — porque é o ligante
co-cristalizado do próprio cristal que o projeto usa, e portanto tem pose
experimental em vez de pose dockeada — não há o que triar: falta só o vetor.

E falta responder a única coisa que pode desqualificar um recrutador
cristalográfico: ele tem ponto de conjugação utilizável, ou todo átomo exposto
dele faz parte do farmacóforo que a E3 reconhece? Pendurar o linker no
farmacóforo destrói o reconhecimento, e descobrir isso depois de montar os
PROTACs é descobrir tarde.

A conta é a MESMA do WP1 — `exit_vector_do_recrutador` importado de lá, não
reimplementado. Duas implementações do mesmo cálculo divergem, e neste projeto
já divergiram três vezes.

A referência de farmacóforo fica em None de propósito: ela existe no WP1 porque
as poses vêm de PDBQT, formato que perde ordem de ligação, e os SMARTS passam a
dar falso positivo. Um SDF de cristal tem as ordens corretas, então aplicar os
padrões direto na molécula é o caminho certo aqui.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parent))
from wp1_select_recruiter import exit_vector_do_recrutador  # noqa: E402

RDLogger.DisableLog("rdApp.*")


def coords_receptor(pdb: Path) -> np.ndarray:
    xyz = [[float(l[30:38]), float(l[38:46]), float(l[46:54])]
           for l in pdb.read_text(errors="ignore").splitlines()
           if l.startswith("ATOM")]
    if not xyz:
        raise SystemExit(f"nenhum registro ATOM em {pdb}")
    return np.array(xyz)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ligante", type=Path, required=True)
    ap.add_argument("--receptor", type=Path, required=True)
    ap.add_argument("--burial", type=int, default=20)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    # removeHs=True, e isto foi achado no teste: com os H explícitos, o
    # conjunto de átomos expostos inclui hidrogênios, e quando nenhum átomo
    # pesado exposto tem H para ceder o `exit_vector_do_recrutador` cai no
    # fallback e escolhe um H como ponto de conjugação. Ponto de conjugação é
    # átomo PESADO por definição — é nele que o linker se liga, no lugar de um
    # H. Com os H implícitos, `GetTotalNumHs()` continua respondendo certo e o
    # conjunto candidato é o correto.
    mol = next(iter(Chem.SDMolSupplier(str(args.ligante.expanduser()),
                                       removeHs=True)), None)
    if mol is None or mol.GetNumConformers() == 0:
        raise SystemExit(f"não li {args.ligante} com coordenadas")
    rec = coords_receptor(args.receptor.expanduser())

    print(f"ligante : {args.ligante.name} — {mol.GetNumAtoms()} átomos")
    print(f"receptor: {args.receptor.name} — {len(rec)} átomos")

    ev = exit_vector_do_recrutador(mol, rec, burial=args.burial, mol_ref=None)

    if ev["atom_symbol"] == "H":
        raise SystemExit(
            f"\n*** o cálculo escolheu um HIDROGÊNIO (átomo {ev['atom_idx']})"
            f" como ponto\n*** de conjugação, o que nunca é válido: o linker"
            f" se liga a um átomo\n*** pesado, no lugar de um H.\n"
            f"*** Isto indica que nenhum átomo pesado exposto tem H para ceder"
            f" —\n*** ou seja, este recrutador não tem ponto de conjugação"
            f" livre fora do\n*** farmacóforo. Veja os excluídos no --out e"
            f" escolha outro recrutador.")

    print()
    print(f"  ponto de conjugação: átomo {ev['atom_idx']} "
          f"({ev['atom_symbol']}, {ev['atom_n_hs']} H)")
    print(f"  exit point ....... {ev['exit_point']}")
    print(f"  exit direction ... {ev['exit_direction']}")
    print(f"  distância ao núcleo ancorado: {ev['dist_ao_nucleo_A']} Å")
    print(f"  vizinhos proteicos do átomo escolhido: {ev['vizinhos_proteicos']}")
    print(f"  enterrados/expostos: {ev['n_atomos_enterrados']}/"
          f"{ev['n_atomos_expostos']}  (burial {ev['burial_min_max']})")
    exc = ev.get("excluidos_por_farmacoforo") or []
    if exc:
        print(f"  excluídos por farmacóforo: {len(exc)} átomo(s)")
        for e in exc[:6]:
            print(f"      átomo {e['atom_idx']}: {e['motivo']}")
    print()
    print("  As DUAS linhas para a config, se este for o recrutador:")
    print(f"    E3_EXIT_POINT=\"{','.join(str(x) for x in ev['exit_point'])}\"")
    print(f"    E3_EXIT_DIRECTION=\"{','.join(str(x) for x in ev['exit_direction'])}\"")

    if args.out:
        args.out.expanduser().write_text(json.dumps(ev, indent=2,
                                                    ensure_ascii=False))
        print(f"\n  {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
