#!/usr/bin/env python
"""
Fase 5 (WP2) — filtra linkers, avalia a geometria no exit vector do recrutador
e monta os sub-complexos recrutador–linker.

Roda no env `mdtools`. Minutos a ~1 h, conforme o tamanho da biblioteca.

    python wp2_build_subcomplexes.py \\
        --linkers Chemspace_PROTACs_linkers_SDF.sdf \\
        --recruiter recrutador.sdf \\
        --receptor-pdb VHL_receptor.pdb \\
        --exit-point 12.3,45.6,7.8 --exit-direction 0.1,0.9,0.4 \\
        --outdir ~/pipeline/wp2

Três correções em relação ao WP2 original estão em `wp2_linker_tools.py`:
bifuncionalidade exigida, confôrmero POSICIONADO no exit vector antes de medir
clash, e pontos de conjugação rotulados para a montagem do WP3 não falhar em
silêncio. Aqui eles são aplicados em lote.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parent))
from wp2_linker_tools import (filter_linkers, evaluate_linker,  # noqa: E402
                              label_attachment_points, cap_for_md,
                              find_dummy, RECRUITER_ISOTOPE, WARHEAD_ISOTOPE)

RDLogger.DisableLog("rdApp.*")


def carregar_sdf(caminhos: list[Path]):
    mols = []
    for c in caminhos:
        c = Path(c).expanduser()
        if not c.exists():
            print(f"  [aviso] não encontrei {c}", file=sys.stderr)
            continue
        n0 = len(mols)
        for m in Chem.SDMolSupplier(str(c), removeHs=True):
            if m is not None:
                mols.append(m)
        print(f"  {len(mols) - n0:5d} moléculas de {c.name}")
    return mols


def receptor_coords(receptor_pdb: Path) -> np.ndarray:
    return np.array([(float(l[30:38]), float(l[38:46]), float(l[46:54]))
                     for l in Path(receptor_pdb).read_text().splitlines()
                     if l.startswith("ATOM") and (l[76:78].strip() or "C") != "H"])


def vetor(txt: str) -> list[float]:
    partes = [p for p in txt.replace(";", ",").split(",") if p.strip()]
    if len(partes) != 3:
        raise SystemExit(f"esperava 3 números em '{txt}'")
    return [float(p) for p in partes]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--linkers", nargs="+", required=True)
    ap.add_argument("--linkers-extra", nargs="*", default=[])
    ap.add_argument("--recruiter", type=Path, required=True,
                    help="SDF do recrutador; átomo 0 = ponto de conjugação")
    ap.add_argument("--receptor-pdb", type=Path, required=True)
    ap.add_argument("--exit-point", required=True)
    ap.add_argument("--exit-direction", required=True)
    ap.add_argument("--atom-range", nargs=2, type=int, default=[5, 40])
    ap.add_argument("--rotb-max", type=int, default=15)
    ap.add_argument("--n-confs", type=int, default=50)
    ap.add_argument("--n-spins", type=int, default=12)
    ap.add_argument("--top-n", type=int, default=10)
    ap.add_argument("--min-frac-sem-clash", type=float, default=0.1,
                    help="fração mínima de confôrmeros sem clash (default 0.1); "
                         "um linker que encaixa em 1 de 50 está forçado")
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()

    out = args.outdir.expanduser()
    (out / "subcomplexes").mkdir(parents=True, exist_ok=True)

    print("[1/4] Carregando bibliotecas de linkers")
    linkers = carregar_sdf([*args.linkers, *args.linkers_extra])
    if not linkers:
        raise SystemExit("nenhum linker carregado")
    print(f"  total: {len(linkers)}")

    print("\n[2/4] Filtrando (bifuncionalidade exigida)")
    ok = filter_linkers(linkers, atom_count_range=tuple(args.atom_range),
                        rotb_max=args.rotb_max)
    if not ok:
        raise SystemExit(
            "nenhum linker passou. Se a maioria caiu por monofuncionalidade, "
            "o subconjunto exportado do catálogo provavelmente veio errado.")

    print(f"\n[3/4] Geometria no exit vector ({len(ok)} linkers x "
          f"{args.n_confs} confôrmeros x {args.n_spins} rotações)")
    rec = receptor_coords(args.receptor_pdb.expanduser())
    ep, ed = vetor(args.exit_point), vetor(args.exit_direction)

    linhas = []
    for i, m in enumerate(ok, 1):
        attach = int(m.GetProp("attach_1_idx"))
        _, res = evaluate_linker(m, attach, ep, ed, rec,
                                 n_confs=args.n_confs, n_spins=args.n_spins)
        if not res:
            continue
        melhor = res[0]
        sem_clash = [r for r in res if r["n_clashes"] == 0]
        linhas.append({
            "idx": i - 1,
            "smiles": Chem.MolToSmiles(m),
            "n_heavy": m.GetNumHeavyAtoms(),
            "attach_1_type": m.GetProp("attach_1_type"),
            "attach_2_type": m.GetProp("attach_2_type"),
            "extension_A": melhor["extension_A"],
            "cos_to_exit": melhor["cos_angle_to_exit"],
            "min_dist_A": melhor["min_dist_to_receptor_A"],
            "frac_sem_clash": round(len(sem_clash) / len(res), 3),
        })
        if i % 25 == 0 or i == len(ok):
            print(f"  {i}/{len(ok)} avaliados", flush=True)

    df = pd.DataFrame(linhas)
    if df.empty:
        raise SystemExit("nenhum linker gerou confôrmero avaliável")

    # frac_sem_clash é o critério principal: um linker que só encaixa num
    # confôrmero está geometricamente forçado, por melhor que ele pontue
    df = df.sort_values(["frac_sem_clash", "cos_to_exit"], ascending=False)
    df.to_csv(out / "linker_ranking.csv", index=False)

    viaveis = df[df["frac_sem_clash"] >= args.min_frac_sem_clash]
    print(f"\n  {len(viaveis)}/{len(df)} com frac_sem_clash >= "
          f"{args.min_frac_sem_clash}")
    if viaveis.empty:
        raise SystemExit(
            "nenhum linker acomoda o exit vector. Confira se --exit-point e "
            "--exit-direction estão no referencial do receptor.")
    escolhidos = viaveis.head(args.top_n)
    print(f"  seguindo com os {len(escolhidos)} melhores")
    print(escolhidos[["smiles", "n_heavy", "extension_A", "cos_to_exit",
                      "frac_sem_clash"]].to_string(index=False))

    print("\n[4/4] Montando sub-complexos recrutador–linker")
    # removeHs: o ponto de conjugação precisa de valência aberta. Com H
    # explícito no átomo 0 a montagem estoura com "Explicit valence ... is
    # greater than permitted", que não diz o que está errado de verdade.
    recruiter = Chem.MolFromMolFile(str(args.recruiter.expanduser()),
                                    removeHs=True)
    if recruiter is None:
        raise SystemExit(f"não consegui ler {args.recruiter}")
    a0 = recruiter.GetAtomWithIdx(0)
    if a0.GetTotalNumHs() < 1:
        raise SystemExit(
            f"o átomo 0 do recrutador ({a0.GetSymbol()}) não tem hidrogênio "
            f"para ceder ao linker. O SDF do recrutador precisa ter o PONTO DE "
            f"CONJUGAÇÃO no índice 0, com valência aberta — mesma convenção "
            f"dos warheads (veja reorder_attachment_first em "
            f"generate_pcsk9_warheads.py).")
    print(f"  recrutador: {recruiter.GetNumHeavyAtoms()} átomos pesados, "
          f"conjugação pelo átomo 0 ({a0.GetSymbol()})")

    from rdkit.Chem import AllChem
    feitos = []
    for r in escolhidos.itertuples():
        m = ok[r.idx]
        sc_id = f"SC{r.idx:04d}"
        try:
            rotulado = label_attachment_points(m, int(m.GetProp("attach_1_idx")),
                                               int(m.GetProp("attach_2_idx")))
            sub = AllChem.ReplaceSubstructs(
                rotulado, Chem.MolFromSmarts(f"[{RECRUITER_ISOTOPE}#0]"),
                recruiter, replacementConnectionPoint=0)[0]
            Chem.SanitizeMol(sub)
            if find_dummy(sub, WARHEAD_ISOTOPE) is None:
                raise ValueError("o [2*] sumiu ao conjugar o recrutador")
        except Exception as exc:
            print(f"  [falha] {sc_id}: {exc}")
            continue

        d = out / "subcomplexes" / sc_id
        d.mkdir(parents=True, exist_ok=True)
        with Chem.SDWriter(str(d / f"{sc_id}.sdf")) as w:
            w.write(sub)                      # -> WP3, [2*] intacto
        with Chem.SDWriter(str(d / f"{sc_id}_capped_md.sdf")) as w:
            w.write(cap_for_md(sub))          # -> MD nível (i)
        feitos.append({"subcomplex_id": sc_id, "linker_smiles": r.smiles,
                       "smiles": Chem.MolToSmiles(sub),
                       "frac_sem_clash": r.frac_sem_clash,
                       "path": str(d / f"{sc_id}.sdf")})

    (out / "subcomplexes_manifest.json").write_text(json.dumps(feitos, indent=2))
    print(f"\n{len(feitos)} sub-complexos em {out / 'subcomplexes'}")
    if not feitos:
        raise SystemExit("nenhum sub-complexo montado")


if __name__ == "__main__":
    main()
