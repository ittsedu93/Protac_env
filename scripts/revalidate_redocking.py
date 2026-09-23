#!/usr/bin/env python
"""
Recalcula o portão de validação por redocking com o RMSD correto.

Use para reconferir o WP1: as poses de redocking VHL/CRBN que passaram pelo
corte de 2,0 Å foram medidas por `redocking_rmsd()` (célula 1.6), que usa
`rdMolAlign.GetBestRMS`. Essa função SUPERPÕE as moléculas antes de medir, então
reporta RMSD ~0 mesmo para uma pose ancorada no bolsão errado. Este script
remede as mesmas poses sem alinhar.

    conda activate mdtools
    python revalidate_redocking.py \\
        --receptor-pdb  .../VHL_6GFZ_receptor.pdb \\
        --ref-ligand    .../6GFZ_ref_ligand.sdf \\
        --poses         ".../redock*/*.pdbqt" \\
        --out           .../wp1_revalidation_VHL.csv

Não roda docking: só remede o que já está no disco. Se nada sobreviver, o
protocolo do WP1 precisa ser refeito antes de qualquer conclusão da triagem.
"""

from __future__ import annotations

import argparse
import csv
import glob
import subprocess
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rmsd_inplace import rmsd_inplace, buried_atom_indices  # noqa: E402

RDLogger.DisableLog("rdApp.warning")
OBABEL_EXE = "/home/soberano/miniconda3/envs/obabel_env/bin/obabel"
DEFAULT_CUTOFF = 2.0


def receptor_heavy_coords(receptor_pdb: Path) -> np.ndarray:
    xyz = []
    for line in Path(receptor_pdb).read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            resname = line[17:20].strip()
            if resname in ("HOH", "WAT"):
                continue
            el = (line[76:78].strip() or line[12:16].strip()[0])
            if el != "H":
                xyz.append((float(line[30:38]), float(line[38:46]),
                            float(line[46:54])))
    return np.array(xyz)


def load_pose(path: Path, tmp_dir: Path):
    """Lê pose .sdf/.mol/.pdb direto; converte .pdbqt via obabel."""
    p = Path(path)
    if p.suffix.lower() == ".pdbqt":
        tmp_dir.mkdir(parents=True, exist_ok=True)
        sdf = tmp_dir / (p.stem + ".sdf")
        if not sdf.exists():
            subprocess.run([OBABEL_EXE, str(p), "-O", str(sdf)],
                           check=True, capture_output=True)
        p = sdf
    if p.suffix.lower() in (".sdf", ".mol"):
        return Chem.MolFromMolFile(str(p), removeHs=True)
    if p.suffix.lower() == ".pdb":
        return Chem.MolFromPDBFile(str(p), removeHs=True)
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--receptor-pdb", type=Path, required=True)
    ap.add_argument("--ref-ligand", type=Path, required=True,
                    help="ligante co-cristalizado (SDF/MOL/PDB)")
    ap.add_argument("--poses", required=True,
                    help="glob das poses de redocking (entre aspas)")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--cutoff", type=float, default=DEFAULT_CUTOFF)
    ap.add_argument("--core-threshold", type=int, default=20,
                    help="vizinhos a 6 Å para um átomo contar como enterrado; "
                         "0 desliga o RMSD de núcleo")
    ap.add_argument("--label", default="", help="rótulo (ex.: VHL_6GFZ)")
    args = ap.parse_args()

    ref = load_pose(args.ref_ligand, args.out.parent / "_tmp")
    if ref is None:
        raise SystemExit(f"não consegui ler {args.ref_ligand}")

    core = None
    if args.core_threshold > 0:
        rec = receptor_heavy_coords(args.receptor_pdb)
        core = buried_atom_indices(ref, rec, threshold=args.core_threshold)
        print(f"núcleo enterrado: {len(core)}/{ref.GetNumAtoms()} átomos pesados")
        if len(core) < 3:
            print("  poucos átomos enterrados — desligando o RMSD de núcleo")
            core = None

    files = sorted(glob.glob(args.poses))
    if not files:
        raise SystemExit(f"nenhuma pose casou com {args.poses}")
    print(f"{len(files)} poses para remedir\n")

    rows, failed = [], 0
    for f in files:
        probe = load_pose(Path(f), args.out.parent / "_tmp")
        if probe is None:
            failed += 1
            continue
        try:
            full = rmsd_inplace(probe, ref)
            cored = rmsd_inplace(probe, ref, atom_indices=core) if core else None
        except ValueError as exc:
            print(f"  [pulado] {Path(f).name}: {exc}")
            failed += 1
            continue
        decisive = cored if cored is not None else full
        rows.append({"label": args.label, "pose": f,
                     "rmsd_core_A": round(cored, 3) if cored is not None else "",
                     "rmsd_full_A": round(full, 3),
                     "passa_corte": decisive <= args.cutoff})

    if not rows:
        raise SystemExit("nenhuma pose mensurável")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    passed = [r for r in rows if r["passa_corte"]]
    key = "rmsd_core_A" if core else "rmsd_full_A"
    best = min(float(r[key]) for r in rows)

    print(f"\n{'=' * 58}")
    print(f"{args.label or 'protocolo'}: {len(passed)}/{len(rows)} poses "
          f"<= {args.cutoff} Å  (melhor: {best:.2f} Å)")
    if failed:
        print(f"  {failed} poses ilegíveis")
    if passed:
        print("  -> protocolo VALIDADO pelo critério correto")
    else:
        print("  -> protocolo NAO VALIDADO.")
        print("     As conclusões da triagem do WP1 que dependiam deste portão")
        print("     precisam ser refeitas: nenhuma pose reproduz o cristal")
        print("     quando o RMSD é medido sem superposição.")
    print(f"{'=' * 58}")
    print(f"\nCSV: {args.out}")


if __name__ == "__main__":
    main()
