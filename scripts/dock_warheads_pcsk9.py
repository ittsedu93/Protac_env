#!/usr/bin/env python
"""
Etapa 3.0b — docking dos warheads PCSK9 e produção dos `Heads` do PRosettaC.

Roda no env `pf_vs`. Dois motores, mesma função de scoring:

  --engine unidock  (PADRÃO) usa o Uni-Dock já instalado no pf_vs, em GPU e em
                    lote. Nada a instalar — e nesta workstation os envs
                    pertencem a outro usuário, então instalar não é opção.
  --engine vina     AutoDock Vina via API Python, se você tiver um env próprio
                    onde deu para instalá-lo.

Confira o motor antes de gastar horas:
    python scripts/docking_engines.py

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
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rmsd_inplace import (rmsd_inplace, buried_atom_indices,  # noqa: E402
                          truncate_to_core)
from docking_engines import (dock_batch, check_engine, parse_pdbqt_score,  # noqa: E402
                             VINA, UNIDOCK)

OBABEL_EXE = "/home/soberano/miniconda3/envs/obabel_env/bin/obabel"

REDOCK_RMSD_CUTOFF = 2.0
ROUND1_N_RUNS = 30
ROUND2_N_RUNS = 100
EXHAUSTIVENESS = 32


# ---------------------------------------------------------------------------
def dock_one(engine, site, ligand_pdbqt, seed, exhaustiveness, out_dir,
             lig_id="lig", **kw):
    """Uma molécula, um seed. Usado na validação por redocking."""
    res = dock_batch(engine, receptor_pdbqt=site["receptor_pdbqt"],
                     ligands={lig_id: ligand_pdbqt}, center=site["center"],
                     box_size=site["box_size"], seed=seed,
                     exhaustiveness=exhaustiveness, out_dir=Path(out_dir), **kw)
    return res[lig_id]


def runs_by_seed(engine, site, ligands: dict, n_runs, exhaustiveness, out_dir,
                 seed_start=1, **kw):
    """N runs independentes de TODOS os ligantes.

    O laço externo é o seed, não o ligante: assim cada seed vira uma única
    chamada em lote na GPU (Uni-Dock) em vez de um processo por molécula.
    Retomável — um seed cuja pasta já tem as poses é pulado pelo próprio motor.

    Imprime progresso por seed com flush: sem isso um lote de horas não produz
    saída nenhuma até terminar, e um `tail -f` do log parece travado.

    Devolve {ligand_id: [{"seed", "best_score", "pose"}, ...]}.
    """
    out_dir = Path(out_dir)
    por_ligante = {lid: [] for lid in ligands}
    t0 = time.time()
    for i in range(n_runs):
        seed = seed_start + i
        t1 = time.time()
        res = dock_batch(engine, receptor_pdbqt=site["receptor_pdbqt"],
                         ligands=ligands, center=site["center"],
                         box_size=site["box_size"], seed=seed,
                         exhaustiveness=exhaustiveness,
                         out_dir=out_dir / f"seed{seed:03d}", **kw)
        n_ok = 0
        for lid, (score, pose) in res.items():
            if score is not None:
                por_ligante[lid].append({"seed": seed, "best_score": score,
                                         "pose": str(pose)})
                n_ok += 1

        dt = time.time() - t1
        decorrido = time.time() - t0
        falta = decorrido / (i + 1) * (n_runs - i - 1)
        print(f"  seed {i + 1:3d}/{n_runs}  {n_ok}/{len(ligands)} poses  "
              f"{dt:6.1f}s  (decorrido {decorrido / 60:.0f} min, "
              f"faltam ~{falta / 60:.0f} min)", flush=True)
    return por_ligante


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
def fraction_in_box(coords, center, box_size) -> float:
    """Fração dos átomos que cabe no box de busca."""
    c = np.asarray(center, dtype=float)
    b = np.asarray(box_size, dtype=float) / 2.0
    dentro = ((coords >= c - b) & (coords <= c + b)).all(axis=1)
    return float(dentro.mean())


def validate_protocol(site: dict, out_dir: Path, exhaustiveness: int,
                      engine: str = UNIDOCK, n_runs: int = 5, **kw) -> bool:
    """Redocking do ligante co-cristalizado, recortado ao núcleo ancorado.

    Por que recortar: o box é dimensionado para o NÚCLEO do ligante, que é o
    tamanho certo para os warheads (230-420 Da). Redockar o ligante INTEIRO
    nesse box é um teste impossível — no 6U26, 19 dos 76 átomos da pose
    cristalográfica caem fora do box, então a pose correta não é sequer
    alcançável e o docking é forçado a dobrar o braço de solvente para dentro.
    O resultado é um RMSD alto que não diz nada sobre o protocolo, depois de
    uma busca longa em 31 torções.

    Recortado ao núcleo (anéis inteiros, valências capeadas), o fragmento cabe
    no box, tem poucas torções e mede exatamente o que interessa: se o
    protocolo reproduz a ancoragem observada no cristal.
    """
    from rdkit import Chem

    ref_sdf = site.get("ref_ligand_sdf")
    if not ref_sdf:
        print("Sem ligante de referência no site.json — o protocolo NÃO pode")
        print("ser validado por redocking. Rode prep_pcsk9_receptor.py com")
        print("--site-mode ligand, ou assuma explicitamente o risco na tese.")
        return False

    out_dir.mkdir(parents=True, exist_ok=True)
    ref_full = Chem.MolFromMolFile(str(ref_sdf), removeHs=True)
    if ref_full is None:
        raise SystemExit(f"não consegui ler {ref_sdf}")

    rec_xyz = receptor_heavy_coords(Path(site["receptor_pdb"]))
    frac = fraction_in_box(ref_full.GetConformer().GetPositions(),
                           site["center"], site["box_size"])
    print(f"  ligante de referência: {ref_full.GetNumHeavyAtoms()} átomos "
          f"pesados, {frac * 100:.0f}% dentro do box")

    core = buried_atom_indices(ref_full, rec_xyz,
                               threshold=site.get("burial_threshold", 20))
    if len(core) < 3:
        raise SystemExit("núcleo enterrado pequeno demais para validar")

    ref = truncate_to_core(ref_full, core)
    ref_heavy = Chem.RemoveHs(ref)
    frac_core = fraction_in_box(ref_heavy.GetConformer().GetPositions(),
                                site["center"], site["box_size"])
    print(f"  recortado ao núcleo  : {ref_heavy.GetNumHeavyAtoms()} átomos "
          f"pesados, {frac_core * 100:.0f}% dentro do box")

    if frac_core < 1.0:
        raise SystemExit(
            f"mesmo recortado, {(1 - frac_core) * 100:.0f}% do núcleo fica fora "
            f"do box — aumente --box-pad no prep_pcsk9_receptor.py")

    core_sdf = out_dir / "ref_core.sdf"
    with Chem.SDWriter(str(core_sdf)) as w:
        w.write(ref)
    core_pdbqt = out_dir / "ref_core.pdbqt"
    subprocess.run([OBABEL_EXE, str(core_sdf), "-O", str(core_pdbqt)],
                   check=True, capture_output=True)

    results = []
    for seed in range(1, n_runs + 1):
        score, pose_pdbqt = dock_one(engine, site, core_pdbqt, seed,
                                     exhaustiveness,
                                     out_dir / f"redock_seed{seed}",
                                     lig_id="ref_core", **kw)
        if score is None:
            print(f"  seed {seed}: o motor não produziu pose")
            continue
        pose_sdf = pdbqt_to_sdf(pose_pdbqt, pose_pdbqt.with_suffix(".sdf"))
        probe = Chem.MolFromMolFile(str(pose_sdf), removeHs=True)
        if probe is None:
            print(f"  seed {seed}: não consegui ler a pose")
            continue
        try:
            rmsd = rmsd_inplace(probe, ref_heavy)
        except ValueError as exc:
            print(f"  seed {seed}: RMSD falhou — {exc}")
            continue
        results.append((seed, score, rmsd))
        print(f"  seed {seed}: score {score:7.2f} kcal/mol | RMSD {rmsd:5.2f} Å")

    if not results:
        print("\nNenhum redocking mensurável — protocolo NÃO validado.")
        return False

    best = min(r[2] for r in results)
    ok = best <= REDOCK_RMSD_CUTOFF
    print(f"\n  melhor RMSD do núcleo: {best:.2f} Å "
          f"(corte {REDOCK_RMSD_CUTOFF} Å) -> "
          f"{'VALIDADO' if ok else 'NAO VALIDADO'}")
    if not ok:
        print("  Não siga para a triagem. Reveja o box, a protonação do")
        print("  receptor e o exhaustiveness antes de dockar a série.")

    (out_dir / "validation.json").write_text(json.dumps(
        {"cutoff": REDOCK_RMSD_CUTOFF, "best_rmsd": best, "validated": ok,
         "ref_full_heavy_atoms": ref_full.GetNumHeavyAtoms(),
         "ref_core_heavy_atoms": ref_heavy.GetNumHeavyAtoms(),
         "ref_full_fraction_in_box": round(frac, 3),
         "runs": [{"seed": s, "score": sc, "rmsd": r}
                  for s, sc, r in results]}, indent=2))
    return ok


# ---------------------------------------------------------------------------
def staged_screening(site, ligands: dict, out_dir: Path, exhaustiveness: int,
                     round1_cutoff: float, engine: str = UNIDOCK, **kw):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nRound 1 — {ROUND1_N_RUNS} runs x {len(ligands)} warheads", flush=True)
    bruto1 = runs_by_seed(engine, site, ligands, ROUND1_N_RUNS, exhaustiveness,
                          out_dir / "round1", **kw)
    r1 = []
    for lid, rows in bruto1.items():
        s = summarize(rows)
        if s:
            s["ligand_id"] = lid
            r1.append(s)
    r1.sort(key=lambda d: d["best_score"])
    for s in r1[:10]:
        print(f"  {s['ligand_id']}: best {s['best_score']:.2f} "
              f"(sd {s['std_score']:.2f}, n={s['n_runs']})")

    survivors = {s["ligand_id"]: ligands[s["ligand_id"]] for s in r1
                 if s["best_score"] <= round1_cutoff}
    print(f"\nRound 1: {len(survivors)}/{len(ligands)} abaixo de {round1_cutoff}", flush=True)
    if not survivors:
        print("  Nenhum sobrevivente. Corte frouxo demais ou sítio errado —")
        print("  reveja antes de afrouxar --round1-cutoff por conveniência.")
        return r1, []

    print(f"\nRound 2 — {ROUND2_N_RUNS} runs x {len(survivors)} sobreviventes", flush=True)
    bruto2 = runs_by_seed(engine, site, survivors, ROUND2_N_RUNS,
                          exhaustiveness, out_dir / "round2", **kw)
    r2 = []
    for lid, rows in bruto2.items():
        s = summarize(rows)
        if s:
            s["ligand_id"] = lid
            r2.append(s)
    r2.sort(key=lambda d: d["best_score"])
    for s in r2[:10]:
        print(f"  {s['ligand_id']}: best {s['best_score']:.2f} "
              f"(sd {s['std_score']:.2f}, n={s['n_runs']})")
    return r1, r2


def export_heads(r2, out_dir: Path, heads_dir: Path):
    """Melhor pose de cada sobrevivente -> SDF posicionado = HeadB do PRosettaC."""
    heads_dir.mkdir(parents=True, exist_ok=True)
    made = []
    for row in r2:
        lid = row["ligand_id"]
        poses = sorted((Path(out_dir) / "round2").glob(f"seed*/{lid}*_out.pdbqt"))
        melhor, melhor_score = None, np.inf
        for p in poses:
            sc = parse_pdbqt_score(p)
            if sc is not None and sc < melhor_score:
                melhor_score, melhor = sc, p
        if melhor is None:
            print(f"  [sem pose] {lid}")
            continue
        sdf = heads_dir / f"{lid}_in_pcsk9.sdf"
        pdbqt_to_sdf(melhor, sdf)
        made.append({"warhead_id": lid, "head_sdf": str(sdf),
                     "best_score": melhor_score, "source_pose": str(melhor)})
        print(f"  {lid}: {melhor_score:.2f} kcal/mol -> {sdf.name}")
    return made


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--site", type=Path, required=True)
    ap.add_argument("--warheads-pdbqt", type=Path,
                    help="pasta dos PDBQT; dispensável com --validate-only, "
                         "que só usa o ligante co-cristalizado")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--warheads-csv", type=Path,
                    help="manifesto, para filtrar por série")
    ap.add_argument("--series", choices=["A", "B", "all"], default="A",
                    help="A = R3 H (N primário, recomendada); B = R3 mantido")
    ap.add_argument("--validate-only", action="store_true")
    ap.add_argument("--skip-validation", action="store_true",
                    help="pula o redocking; assume o risco explicitamente")
    ap.add_argument("--engine", choices=[UNIDOCK, VINA], default=UNIDOCK,
                    help="unidock (padrão, GPU, já instalado) ou vina")
    ap.add_argument("--unidock-exe", help="caminho do executável, se não no PATH")
    ap.add_argument("--batch-size", type=int, default=40,
                    help="ligantes por lote na GPU (default 40); reduza se a "
                         "memória da GPU estourar")
    ap.add_argument("--search-mode", choices=["fast", "balance", "detail"],
                    help="preset do Uni-Dock; substitui --exhaustiveness")
    ap.add_argument("--exhaustiveness", type=int, default=EXHAUSTIVENESS)
    ap.add_argument("--round1-cutoff", type=float, default=-7.0)
    args = ap.parse_args()

    print("[motor]", check_engine(args.engine, args.unidock_exe))
    kw = {}
    if args.engine == UNIDOCK:
        kw = dict(unidock_exe=args.unidock_exe, batch_size=args.batch_size,
                  search_mode=args.search_mode)

    site = json.loads(args.site.read_text())
    out = args.outdir.expanduser()
    out.mkdir(parents=True, exist_ok=True)

    print(f"Receptor : {Path(site['receptor_pdbqt']).name}")
    print(f"Centro   : {site['center']}")
    print(f"Box (Å)  : {site['box_size']}")
    print(f"Origem   : {site['provenance']}\n")

    if not args.skip_validation:
        print("[validação] redocking do ligante co-cristalizado")
        ok = validate_protocol(site, out / "validation", args.exhaustiveness,
                               engine=args.engine, **kw)
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

    if args.warheads_pdbqt is None:
        raise SystemExit("--warheads-pdbqt é obrigatório para a triagem")
    ligands = {p.stem: p for p in sorted(args.warheads_pdbqt.glob("*.pdbqt"))
               if ids is None or p.stem in ids}
    if not ligands:
        raise SystemExit(f"nenhum PDBQT em {args.warheads_pdbqt} "
                         f"(rodou o prepare_pdbqt.sh?)")
    print(f"{len(ligands)} warheads para dockar")

    r1, r2 = staged_screening(site, ligands, out, args.exhaustiveness,
                              args.round1_cutoff, engine=args.engine, **kw)

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
