#!/usr/bin/env python
"""
Fase 6 (WP3) — monta os PROTACs completos e emite os jobs de PRosettaC/AF3.

Roda no env `mdtools`. Minutos.

    python wp3_assemble_protacs.py \\
        --subcomplexes ~/pipeline/wp2/subcomplexes \\
        --warhead-ranking ~/pipeline/warhead_ranking.csv \\
        --warheads-sdf ~/PCSK9_warheads/sdf \\
        --heads-pcsk9 ~/PCSK9_docking/docking/heads_pcsk9 \\
        --e3-structure VHL_receptor.pdb --e3-head recrutador.sdf \\
        --pcsk9-structure 6U26_receptor.pdb \\
        --outdir ~/pipeline/wp3

A montagem usa `assemble_protac`, que procura o `[2*]` rotulado e LEVANTA
EXCEÇÃO se ele não existir. A versão original do notebook recebia de volta a
molécula inalterada sem erro, e o lote inteiro viraria recrutador-linker sem
warhead nenhum, em silêncio.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors

sys.path.insert(0, str(Path(__file__).resolve().parent))
from wp2_linker_tools import assemble_protac, WARHEAD_ISOTOPE  # noqa: E402

RDLogger.DisableLog("rdApp.*")


def selecionar_warheads(ranking_csv: Path, top_n: int, cos_min: float):
    df = pd.read_csv(ranking_csv)
    faltando = [c for c in ("n_exposto", "alinhado_com_063", "pose_estavel",
                            "ligand_efficiency") if c not in df.columns]
    if faltando:
        raise SystemExit(f"colunas ausentes em {ranking_csv}: {faltando}")

    viaveis = df[df["n_exposto"].fillna(False)
                 & (df["cos_vetor_saida"].fillna(-1) >= cos_min)
                 & df["pose_estavel"].fillna(False)].copy()
    if viaveis.empty:
        raise SystemExit(
            f"nenhum warhead passa nos critérios (cos >= {cos_min}). "
            f"Afrouxe --cos-min e registre a escolha.")

    # eficiência de ligante, não score bruto: o score do Vina correlaciona
    # fortemente com o tamanho, então ordenar por ele ordena por massa
    viaveis = viaveis.sort_values("ligand_efficiency", ascending=False)
    return viaveis.head(top_n)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--subcomplexes", type=Path, required=True)
    ap.add_argument("--warhead-ranking", type=Path, required=True)
    ap.add_argument("--warheads-sdf", type=Path, required=True)
    ap.add_argument("--heads-pcsk9", type=Path, required=True)
    ap.add_argument("--top-warheads", type=int, default=15)
    ap.add_argument("--cos-min", type=float, default=0.3)
    ap.add_argument("--e3-structure", type=Path)
    ap.add_argument("--e3-chain", default="A")
    ap.add_argument("--e3-head", type=Path)
    ap.add_argument("--e3-anchor-serial", type=int, default=1)
    ap.add_argument("--pcsk9-structure", type=Path)
    ap.add_argument("--pcsk9-chain", default="B")
    ap.add_argument("--anchor-serial", type=int, default=1)
    ap.add_argument("--prosettac-dir", default="/mnt/hd2tb/Documentos/PRosettaC")
    ap.add_argument("--prosettac-full", default="False")
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()

    out = args.outdir.expanduser()
    (out / "protacs").mkdir(parents=True, exist_ok=True)

    print("[1/3] Selecionando warheads")
    wh = selecionar_warheads(args.warhead_ranking.expanduser(),
                             args.top_warheads, args.cos_min)
    print(f"  {len(wh)} warheads (por eficiência de ligante)")
    print(wh[[c for c in ("warhead_id", "r1", "r2", "r3", "mw",
                          "ligand_efficiency") if c in wh.columns]]
          .to_string(index=False))

    subs = sorted(p for p in Path(args.subcomplexes).expanduser().rglob("*.sdf")
                  if "_capped_md" not in p.name)
    if not subs:
        raise SystemExit(f"nenhum sub-complexo em {args.subcomplexes}")
    print(f"\n  {len(subs)} sub-complexos")

    print(f"\n[2/3] Montando {len(subs)} x {len(wh)} PROTACs")
    linhas, falhas = [], 0
    for sp in subs:
        sc_mol = Chem.MolFromMolFile(str(sp))
        if sc_mol is None:
            print(f"  [ilegível] {sp.name}")
            continue
        sc_id = sp.stem
        for r in wh.itertuples():
            wid = r.warhead_id
            wsdf = Path(args.warheads_sdf).expanduser() / f"{wid}.sdf"
            try:
                # MolFromMolFile levanta OSError em arquivo ausente ou ilegível,
                # não devolve None: sem o try, um único SDF faltando mata o lote
                wmol = Chem.MolFromMolFile(str(wsdf))
            except OSError as exc:
                print(f"  [SDF ilegível] {wid}: {exc}")
                continue
            if wmol is None:
                print(f"  [sem SDF] {wid}")
                continue
            try:
                full = assemble_protac(sc_mol, wmol)
            except ValueError as exc:
                falhas += 1
                print(f"  [falha] {sc_id} x {wid}: {exc}")
                continue
            cand = f"{sc_id}__{wid}"
            smi = Chem.MolToSmiles(full)
            d = out / "protacs" / cand
            d.mkdir(parents=True, exist_ok=True)
            (d / "protac.smi").write_text(f"{smi}\t{cand}\n")
            linhas.append({
                "candidate_id": cand, "subcomplex_id": sc_id, "warhead_id": wid,
                "protac_smiles": smi,
                "protac_mw": round(Descriptors.MolWt(full), 1),
                "protac_rotb": rdMolDescriptors.CalcNumRotatableBonds(full),
                "protac_clogp": round(Descriptors.MolLogP(full), 2),
                "protac_smi_path": str(d / "protac.smi"),
            })

    if not linhas:
        raise SystemExit("nenhum PROTAC montado")
    df = pd.DataFrame(linhas)
    df.to_csv(out / "protac_candidates.csv", index=False)
    print(f"\n  {len(df)} PROTACs montados ({falhas} falhas)")
    print(f"  MW: {df['protac_mw'].min():.0f}–{df['protac_mw'].max():.0f} Da | "
          f"rotB: {df['protac_rotb'].min()}–{df['protac_rotb'].max()}")

    print("\n[3/3] Emitindo jobs do PRosettaC")
    if not (args.e3_structure and args.e3_head and args.pcsk9_structure):
        print("  [pulado] faltam --e3-structure / --e3-head / --pcsk9-structure")
        return

    root = out / "prosettac"
    root.mkdir(parents=True, exist_ok=True)
    cmds = []
    for r in df.itertuples():
        head_b = Path(args.heads_pcsk9).expanduser() / f"{r.warhead_id}_in_pcsk9.sdf"
        if not head_b.exists():
            print(f"  [sem Head] {r.warhead_id}")
            continue
        # o Head veio de PDBQT e perdeu as ordens de ligação; o PRosettaC
        # precisa da química certa, então corrige pelo SDF do gerador
        corrigido = head_b.parent / f"{r.warhead_id}_in_pcsk9_bo.sdf"
        if not corrigido.exists():
            try:
                from rdkit.Chem import AllChem
                tmpl = Chem.MolFromMolFile(
                    str(Path(args.warheads_sdf).expanduser() / f"{r.warhead_id}.sdf"))
                posed = Chem.MolFromMolFile(str(head_b))
                if tmpl is not None and posed is not None:
                    fix = AllChem.AssignBondOrdersFromTemplate(tmpl, posed)
                    with Chem.SDWriter(str(corrigido)) as w:
                        w.write(fix)
            except Exception as exc:
                print(f"  [ordens não corrigidas] {r.warhead_id}: {exc}")
        if corrigido.exists():
            head_b = corrigido
        d = root / r.candidate_id
        d.mkdir(parents=True, exist_ok=True)
        (d / "prosetta_config.txt").write_text(
            f"Structures: {args.e3_structure} {args.pcsk9_structure}\n"
            f"Chains: {args.e3_chain} {args.pcsk9_chain}\n"
            f"Heads: {args.e3_head} {head_b}\n"
            f"Anchor atoms: {args.e3_anchor_serial} {args.anchor_serial}\n"
            f"Protac: {r.protac_smi_path}\n"
            f"Full: {args.prosettac_full}\n"
            f"ClusterName: SLURM\n")     # sem isto o default é PBS
        cmds.append(f"cd {d} && nohup {args.prosettac_dir}/run_prosettac.sh "
                    f"{d} prosetta_config.txt > run_prosettac.log 2>&1 &")

    (root / "launch_all.sh").write_text(
        "#!/usr/bin/env bash\nset -u\n"
        "# ATENÇÃO: rode UM job antes do lote e confira o 'Anchor atoms'.\n"
        "# A convenção 0-based vs 1-based varia por build do PRosettaC.\n"
        + "\n".join(cmds) + "\n")
    (root / "launch_all.sh").chmod(0o755)

    (root / "launch_one.sh").write_text(
        "#!/usr/bin/env bash\nset -u\n" + (cmds[0] if cmds else "") + "\n")
    (root / "launch_one.sh").chmod(0o755)

    print(f"  {len(cmds)} jobs em {root}")
    print(f"\n  PRIMEIRO:  bash {root}/launch_one.sh")
    print(f"  confira o config e o resultado, DEPOIS: bash {root}/launch_all.sh")


if __name__ == "__main__":
    main()
