#!/usr/bin/env python
"""
Etapa 3.0b — docking dos warheads PCSK9 e produção dos `Heads` do PRosettaC.

Roda no env `pf_vs` (precisa de `vina`; `pip install vina==1.2.7`).
Lê o `pcsk9_site.json` escrito por `prep_pcsk9_receptor.py`.

    conda activate pf_vs
    # 1) valida o protocolo por redocking do ligante co-cristalizado
    python dock_warheads_pcsk9.py --site .../pcsk9_site.json \\
        --warheads-pdbqt ~/PCSK9_warheads/pdbqt --outdir .../docking --validate-only
    # 2) só depois, a triagem em dois estágios
    python dock_warheads_pcsk9.py --site .../pcsk9_site.json \\
        --warheads-pdbqt ~/PCSK9_warheads/pdbqt --outdir .../docking \\
        --warheads-csv ~/PCSK9_warheads/pcsk9_warheads.csv --series A

Processo longo: rode com `nohup ... &` ou via SLURM, nunca dentro do kernel do
Jupyter. O script é retomável — candidatos já concluídos são pulados.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rmsd_inplace import rmsd_inplace, buried_atom_indices  # noqa: E402

OBABEL_EXE = "/home/soberano/miniconda3/envs/obabel_env/bin/obabel"

REDOCK_RMSD_CUTOFF = 2.0
ROUND1_N_RUNS = 30
ROUND2_N_RUNS = 100
EXHAUSTIVENESS = 32


# ---------------------------------------------------------------------------
def run_vina_once(receptor_pdbqt, ligand_pdbqt, center, box_size, seed,
                  exhaustiveness, n_poses, out_pose_pdbqt):
    from vina import Vina
    v = Vina(sf_name="vina", seed=seed, cpu=0, verbosity=0)
    v.set_receptor(str(receptor_pdbqt))
    v.set_ligand_from_file(str(ligand_pdbqt))
    v.compute_vina_maps(center=list(center), box_size=list(box_size))
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    v.write_poses(str(out_pose_pdbqt), n_poses=n_poses, overwrite=True)
    return v.energies(n_poses=n_poses)[:, 0].tolist()


def independent_runs(receptor_pdbqt, ligand_pdbqt, center, box_size, n_runs,
                     exhaustiveness, out_dir, seed_start=1):
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for i in range(n_runs):
        seed = seed_start + i
        pose = out_dir / f"run{seed:03d}.pdbqt"
        if pose.exists():                       # retomada
            continue
        scores = run_vina_once(receptor_pdbqt, ligand_pdbqt, center, box_size,
                               seed, exhaustiveness, 1, pose)
        rows.append({"run": i + 1, "seed": seed, "best_score": scores[0],
                     "pose": str(pose)})
    return rows


def summarize(rows):
    s = [r["best_score"] for r in rows]
    if not s:
        return {}
    return {"best_score": min(s), "mean_score": float(np.mean(s)),
            "std_score": float(np.std(s)), "n_runs": len(s)}


def pdbqt_to_sdf(pdbqt: Path, sdf: Path):
    subprocess.run([OBABEL_EXE, str(pdbqt), "-O", str(sdf)], check=True,
                   capture_output=True)
    return sdf


def receptor_heavy_coords(receptor_pdb: Path) -> np.ndarray:
    xyz = []
    for line in Path(receptor_pdb).read_text().splitlines():
        if line.startswith("ATOM"):
            el = (line[76:78].strip() or line[12:16].strip()[0])
            if el != "H":
                xyz.append((float(line[30:38]), float(line[38:46]),
                            float(line[46:54])))
    return np.array(xyz)


# ---------------------------------------------------------------------------
def validate_protocol(site: dict, out_dir: Path, exhaustiveness: int,
                      n_runs: int = 5) -> bool:
    """Redocking do ligante co-cristalizado.

    Dois números: RMSD do núcleo enterrado (o que decide) e RMSD do ligante
    inteiro (informativo). O `063` tem um braço de ~15 Å no solvente que o
    docking não tem como reproduzir e que não faz parte do que se está
    validando — reprovar por causa dele seria reprovar o protocolo errado.
    """
    from rdkit import Chem

    ref_sdf = site.get("ref_ligand_sdf")
    ref_pdbqt = site.get("ref_ligand_pdbqt")
    if not ref_sdf or not ref_pdbqt:
        print("Sem ligante de referência no site.json — o protocolo NÃO pode")
        print("ser validado por redocking. Rode prep_pcsk9_receptor.py com")
        print("--site-mode ligand, ou assuma explicitamente o risco na tese.")
        return False

    out_dir.mkdir(parents=True, exist_ok=True)
    ref = Chem.MolFromMolFile(str(ref_sdf), removeHs=True)
    if ref is None:
        raise SystemExit(f"não consegui ler {ref_sdf}")

    rec_xyz = receptor_heavy_coords(Path(site["receptor_pdb"]))
    core = buried_atom_indices(ref, rec_xyz,
                               threshold=site.get("burial_threshold", 20))
    print(f"  núcleo enterrado: {len(core)}/{ref.GetNumAtoms()} átomos pesados")

    results = []
    for seed in range(1, n_runs + 1):
        pose_pdbqt = out_dir / f"redock_seed{seed}.pdbqt"
        scores = run_vina_once(site["receptor_pdbqt"], ref_pdbqt,
                               site["center"], site["box_size"], seed,
                               exhaustiveness, 1, pose_pdbqt)
        pose_sdf = pdbqt_to_sdf(pose_pdbqt, pose_pdbqt.with_suffix(".sdf"))
        probe = Chem.MolFromMolFile(str(pose_sdf), removeHs=True)
        if probe is None:
            print(f"  seed {seed}: não consegui ler a pose")
            continue
        try:
            rmsd_core = rmsd_inplace(probe, ref, atom_indices=core)
            rmsd_full = rmsd_inplace(probe, ref)
        except ValueError as exc:
            print(f"  seed {seed}: RMSD falhou — {exc}")
            continue
        results.append((seed, scores[0], rmsd_core, rmsd_full))
        print(f"  seed {seed}: score {scores[0]:7.2f} kcal/mol | "
              f"RMSD núcleo {rmsd_core:5.2f} Å | ligante inteiro {rmsd_full:5.2f} Å")

    if not results:
        print("\nNenhum redocking mensurável — protocolo NÃO validado.")
        return False

    best = min(r[2] for r in results)
    ok = best <= REDOCK_RMSD_CUTOFF
    print(f"\n  melhor RMSD de núcleo: {best:.2f} Å "
          f"(corte {REDOCK_RMSD_CUTOFF} Å) -> "
          f"{'VALIDADO' if ok else 'NAO VALIDADO'}")
    if not ok:
        print("  Não siga para a triagem. Reveja o box, a protonação do")
        print("  receptor e o exhaustiveness antes de dockar a série.")

    (out_dir / "validation.json").write_text(json.dumps(
        {"cutoff": REDOCK_RMSD_CUTOFF, "best_core_rmsd": best, "validated": ok,
         "n_core_atoms": len(core),
         "runs": [{"seed": s, "score": sc, "rmsd_core": rc, "rmsd_full": rf}
                  for s, sc, rc, rf in results]}, indent=2))
    return ok


# ---------------------------------------------------------------------------
def staged_screening(site, ligands: dict, out_dir: Path, exhaustiveness: int,
                     round1_cutoff: float):
    out_dir.mkdir(parents=True, exist_ok=True)
    r1 = []
    print(f"\nRound 1 — {ROUND1_N_RUNS} runs x {len(ligands)} warheads")
    for i, (lid, lpath) in enumerate(ligands.items(), 1):
        rows = independent_runs(site["receptor_pdbqt"], lpath, site["center"],
                                site["box_size"], ROUND1_N_RUNS, exhaustiveness,
                                out_dir / "round1" / lid)
        s = summarize(rows)
        if s:
            s["ligand_id"] = lid
            r1.append(s)
            print(f"  [{i:3d}/{len(ligands)}] {lid}: best {s['best_score']:.2f} "
                  f"(sd {s['std_score']:.2f})")

    survivors = [r["ligand_id"] for r in r1 if r["best_score"] <= round1_cutoff]
    print(f"\nRound 1: {len(survivors)}/{len(ligands)} abaixo de {round1_cutoff}")

    r2 = []
    print(f"\nRound 2 — {ROUND2_N_RUNS} runs x {len(survivors)} sobreviventes")
    for i, lid in enumerate(survivors, 1):
        rows = independent_runs(site["receptor_pdbqt"], ligands[lid],
                                site["center"], site["box_size"],
                                ROUND2_N_RUNS, exhaustiveness,
                                out_dir / "round2" / lid)
        s = summarize(rows)
        if s:
            s["ligand_id"] = lid
            r2.append(s)
            print(f"  [{i:3d}/{len(survivors)}] {lid}: best {s['best_score']:.2f} "
                  f"(sd {s['std_score']:.2f})")
    return r1, r2


def export_heads(r2, out_dir: Path, heads_dir: Path):
    """Melhor pose de cada sobrevivente -> SDF posicionado = HeadB do PRosettaC."""
    heads_dir.mkdir(parents=True, exist_ok=True)
    made = []
    for row in r2:
        lid = row["ligand_id"]
        runs = sorted((out_dir / "round2" / lid).glob("run*.pdbqt"))
        if not runs:
            continue
        best_pose, best_score = None, np.inf
        for p in runs:
            for line in p.read_text().splitlines():
                if line.startswith("REMARK VINA RESULT"):
                    sc = float(line.split()[3])
                    if sc < best_score:
                        best_score, best_pose = sc, p
                    break
        if best_pose is None:
            continue
        sdf = heads_dir / f"{lid}_in_pcsk9.sdf"
        pdbqt_to_sdf(best_pose, sdf)
        made.append({"warhead_id": lid, "head_sdf": str(sdf),
                     "best_score": best_score, "source_pose": str(best_pose)})
        print(f"  {lid}: {best_score:.2f} kcal/mol -> {sdf.name}")
    return made


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--site", type=Path, required=True)
    ap.add_argument("--warheads-pdbqt", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--warheads-csv", type=Path,
                    help="manifesto, para filtrar por série")
    ap.add_argument("--series", choices=["A", "B", "all"], default="A",
                    help="A = R3 H (N primário, recomendada); B = R3 mantido")
    ap.add_argument("--validate-only", action="store_true")
    ap.add_argument("--skip-validation", action="store_true",
                    help="pula o redocking; assume o risco explicitamente")
    ap.add_argument("--exhaustiveness", type=int, default=EXHAUSTIVENESS)
    ap.add_argument("--round1-cutoff", type=float, default=-7.0)
    args = ap.parse_args()

    site = json.loads(args.site.read_text())
    out = args.outdir.expanduser()
    out.mkdir(parents=True, exist_ok=True)

    print(f"Receptor : {Path(site['receptor_pdbqt']).name}")
    print(f"Centro   : {site['center']}")
    print(f"Box (Å)  : {site['box_size']}")
    print(f"Origem   : {site['provenance']}\n")

    if not args.skip_validation:
        print("[validação] redocking do ligante co-cristalizado")
        ok = validate_protocol(site, out / "validation", args.exhaustiveness)
        if args.validate_only:
            return
        if not ok:
            raise SystemExit(
                "\nProtocolo não validado — parando. Use --skip-validation "
                "apenas se for registrar essa decisão por escrito.")
    elif args.validate_only:
        raise SystemExit("--validate-only com --skip-validation não faz sentido")

    # seleção dos warheads
    ids = None
    if args.warheads_csv and args.series != "all":
        rows = list(csv.DictReader(open(args.warheads_csv)))
        keep_r3 = ["H"] if args.series == "A" else ["Ms", "Tet", "CH2Tet"]
        ids = {r["warhead_id"] for r in rows
               if r["filter_flags"] == "ok" and r["r3"] in keep_r3}
        print(f"\nSérie {args.series}: {len(ids)} warheads do manifesto")

    ligands = {p.stem: p for p in sorted(args.warheads_pdbqt.glob("*.pdbqt"))
               if ids is None or p.stem in ids}
    if not ligands:
        raise SystemExit(f"nenhum PDBQT em {args.warheads_pdbqt} "
                         f"(rodou o prepare_pdbqt.sh?)")
    print(f"{len(ligands)} warheads para dockar")

    r1, r2 = staged_screening(site, ligands, out, args.exhaustiveness,
                              args.round1_cutoff)

    for name, data in (("round1_scores.csv", r1), ("round2_scores.csv", r2)):
        if data:
            with open(out / name, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(data[0].keys()))
                w.writeheader()
                w.writerows(data)

    print("\nExportando Heads posicionados para o PRosettaC:")
    heads = export_heads(r2, out, out / "heads_pcsk9")
    (out / "heads_manifest.json").write_text(json.dumps(heads, indent=2))

    print(f"\n{len(heads)} Heads em {out / 'heads_pcsk9'}")
    print("Inspecione as poses antes de seguir para o WP3:")
    print(f"  chimerax {site['receptor_pdb']} {site.get('ref_ligand_sdf','')} "
          f"{out / 'heads_pcsk9'}/*.sdf")


if __name__ == "__main__":
    main()
