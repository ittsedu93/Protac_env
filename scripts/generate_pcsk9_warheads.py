#!/usr/bin/env python
"""
Enumeração combinatória dos warheads PCSK9 (scaffold feniletilamina) para o WP3
do pipeline PROTAC (protac_wp1_wp3_pipeline.ipynb, Seção 3).

Scaffold (Figura 3):

        R1
        |
       / \\
      |   |          R1 = -OAr, ciclopropil, ciclobutil  (+ F, Cl, OH)
       \\ /           R2 = -CF3, naftaleno                (+ F, Cl, OH)
      /   \\          R3 = metilsulfonil, tetrazol        (+ H = amina livre)
    R2     CH2-CH2-NH-R3

Anel 1,3,5-trissubstituído: R1 e R2 nas posições meta, cadeia etilamina na
terceira posição meta.

O que o script produz (tudo em --outdir):

  pcsk9_warheads.csv              manifesto: id, SMILES, R1/R2/R3, descritores,
                                  átomo do vetor de saída, carga em pH 7,4
  pcsk9_warheads.smi              SMILES + id (entrada rápida para outros tools)
  pcsk9_warheads_3d.sdf           todos os confôrmeros de menor energia, 3D
  sdf/<id>.sdf                    um SDF 3D por warhead  -> `Heads` do PRosettaC
  sdf_dummy/<id>_dummy.sdf        mesma molécula com um átomo dummy [*] no ponto
                                  de conjugação, átomo 0 = ponto de ancoragem
                                  -> entrada de assemble_full_protac() no WP3
  prepare_pdbqt.sh                comando Meeko (env pf_vs) para gerar PDBQT

Uso:
    python generate_pcsk9_warheads.py --outdir ~/PCSK9_warheads
    python generate_pcsk9_warheads.py --outdir ... --no-3d        # só 2D/CSV
    python generate_pcsk9_warheads.py --outdir ... --n-confs 100
    python generate_pcsk9_warheads.py --outdir ... --only-figure  # sem F/Cl/OH extras

Requer RDKit (env `mdtools`).
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.warning")

# ---------------------------------------------------------------------------
# 1. Definição dos grupos R
# ---------------------------------------------------------------------------
# Cada entrada: rotulo_curto -> (SMILES do substituinte, origem)
#   origem "figura" = explicitamente listado na Figura 3
#   origem "quadro" = variação F/Cl/OH marcada nos quadros magenta do esquema

R1_GROUPS = {
    # --- explicitamente propostos como R1 na figura ---
    "OPh":      ("Oc1ccccc1",        "figura"),   # -OAr (aril = fenil)
    "O-4F-Ph":  ("Oc1ccc(F)cc1",     "figura"),   # -OAr (aril = 4-fluorofenil)
    "O-4Cl-Ph": ("Oc1ccc(Cl)cc1",    "figura"),   # -OAr (aril = 4-clorofenil)
    "cPr":      ("C1CC1",            "figura"),   # ciclopropil
    "cBu":      ("C1CCC1",           "figura"),   # ciclobutil
    # --- variações do quadro magenta ---
    "F":        ("F",                "quadro"),
    "Cl":       ("Cl",               "quadro"),
    "OH":       ("O",                "quadro"),
}

R2_GROUPS = {
    # --- explicitamente propostos como R2 na figura ---
    "CF3":      ("C(F)(F)F",          "figura"),
    "2-Naph":   ("c1ccc2ccccc2c1",    "figura"),  # naftalen-2-il
    "1-Naph":   ("c1cccc2ccccc12",    "figura"),  # naftalen-1-il
    # --- variações do quadro magenta ---
    "F":        ("F",                 "quadro"),
    "Cl":       ("Cl",                "quadro"),
    "OH":       ("O",                 "quadro"),
}

# R3 vai no nitrogênio da etilamina.
R3_GROUPS = {
    "H":        ("",                  "quadro"),  # amina primária (produto direto)
    "Ms":       ("S(=O)(=O)C",        "figura"),  # metilsulfonil (sulfonamida)
    "Tet":      ("c1nnn[nH]1",        "figura"),  # 1H-tetrazol-5-il
    "CH2Tet":   ("Cc1nnn[nH]1",       "figura"),  # (1H-tetrazol-5-il)metil
}

# Núcleo 1,3,5-trissubstituído com pontos de ligação mapeados.
# [*:1] = R1, [*:2] = R2, [*:3] = R3 (no nitrogênio da etilamina).
# A montagem usa Chem.molzip (e NÃO concatenação de strings), porque os SMILES
# dos substituintes têm anéis próprios e o dígito de fechamento colidiria com o
# do anel central.
SCAFFOLD_CORE = "[*:1]c1cc([*:2])cc(CCN[*:3])c1"

# Filtros (informativos: marcam a linha, NÃO descartam — warhead é fragmento)
FILTERS = dict(mw_max=400.0, hbd_max=3, hba_max=6, rotb_max=8, clogp_max=4.5)


# ---------------------------------------------------------------------------
# 2. Enumeração
# ---------------------------------------------------------------------------
def build_mol(r1_smi: str, r2_smi: str, r3_smi: str):
    """Monta o warhead ligando os três substituintes ao núcleo via molzip."""
    core = Chem.MolFromSmiles(SCAFFOLD_CORE)
    combo = Chem.RWMol(core)
    for map_num, sub in ((1, r1_smi), (2, r2_smi), (3, r3_smi)):
        frag = Chem.MolFromSmiles(f"[*:{map_num}]{sub}" if sub else f"[*:{map_num}][H]")
        if frag is None:
            return None
        combo.InsertMol(frag)
    zipped = Chem.molzip(combo.GetMol())
    if zipped is None:
        return None
    zipped = Chem.RemoveHs(zipped)
    Chem.SanitizeMol(zipped)
    return zipped


def enumerate_warheads(only_figure: bool = False):
    """Gera todas as combinações R1 x R2 x R3, deduplicadas por SMILES canônico."""

    def keep(entry):
        return entry[1] == "figura" if only_figure else True

    r1s = {k: v for k, v in R1_GROUPS.items() if keep(v)}
    r2s = {k: v for k, v in R2_GROUPS.items() if keep(v)}
    r3s = {k: v for k, v in R3_GROUPS.items() if keep(v) or k == "H"}

    seen: dict[str, str] = {}   # canonical smiles -> id do primeiro que apareceu
    records = []
    idx = 0

    for r1_name, (r1_smi, r1_src) in r1s.items():
        for r2_name, (r2_smi, r2_src) in r2s.items():
            for r3_name, (r3_smi, r3_src) in r3s.items():
                mol = build_mol(r1_smi, r2_smi, r3_smi)
                if mol is None:
                    print(f"  [aviso] montagem falhou: "
                          f"R1={r1_name} R2={r2_name} R3={r3_name}", file=sys.stderr)
                    continue
                can = Chem.MolToSmiles(mol)
                if can in seen:
                    # R1/R2 são ambos meta: pares simétricos colapsam aqui
                    continue
                idx += 1
                wid = f"WH{idx:03d}"
                seen[can] = wid
                records.append(
                    dict(
                        warhead_id=wid,
                        smiles=can,
                        r1=r1_name,
                        r2=r2_name,
                        r3=r3_name,
                        source=("figura" if {r1_src, r2_src, r3_src} == {"figura"}
                                else "figura+quadro"),
                    )
                )
    return records


# ---------------------------------------------------------------------------
# 3. Descritores, carga em pH 7,4 e vetor de saída
# ---------------------------------------------------------------------------
AMINE_PATTERNS = [
    Chem.MolFromSmarts("[NX3;H2;!$(N[!#6]);!$(N=*)][CX4]"),          # -CH2-NH2
    Chem.MolFromSmarts("[NX3;H1;!$(N[S,C]=O);!$(Nc)][CX4][CX4]"),     # -CH2-NH-CH2-
]
SULFONAMIDE_N = Chem.MolFromSmarts("[NX3][SX4](=O)(=O)")
TETRAZOLE = Chem.MolFromSmarts("c1nnn[nH]1")


def describe(mol) -> dict:
    return dict(
        mw=round(Descriptors.MolWt(mol), 2),
        clogp=round(Descriptors.MolLogP(mol), 2),
        tpsa=round(rdMolDescriptors.CalcTPSA(mol), 2),
        hbd=rdMolDescriptors.CalcNumHBD(mol),
        hba=rdMolDescriptors.CalcNumHBA(mol),
        rotb=rdMolDescriptors.CalcNumRotatableBonds(mol),
        heavy_atoms=mol.GetNumHeavyAtoms(),
        formal_charge=Chem.GetFormalCharge(mol),
        fsp3=round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
    )


def flag_filters(d: dict) -> str:
    """Não descarta nada; só aponta o que sai da faixa típica de warhead."""
    bad = []
    if d["mw"] > FILTERS["mw_max"]:
        bad.append("MW")
    if d["hbd"] > FILTERS["hbd_max"]:
        bad.append("HBD")
    if d["hba"] > FILTERS["hba_max"]:
        bad.append("HBA")
    if d["rotb"] > FILTERS["rotb_max"]:
        bad.append("RotB")
    if d["clogp"] > FILTERS["clogp_max"]:
        bad.append("cLogP")
    return ",".join(bad) if bad else "ok"


def charge_at_ph74(mol) -> int:
    """Estimativa grosseira: amina alifática protona (+1), tetrazol desprotona (-1),
    sulfonamida N-H fica neutra. Serve para `acpype -n` e para o preparo Meeko."""
    q = 0
    if any(mol.HasSubstructMatch(p) for p in AMINE_PATTERNS if p is not None):
        q += 1
    if mol.HasSubstructMatch(TETRAZOLE):
        q -= 1
    return q


def find_exit_vector(mol) -> tuple[int, str]:
    """Átomo pelo qual o linker do PROTAC será conjugado.

    Regra: o nitrogênio da etilamina é o ponto de conjugação (é exatamente a
    posição que a Figura 3 chama de R3 — o linker entra no lugar de R3, ou no
    N-H remanescente). Retorna (índice do átomo, justificativa).
    """
    # nitrogênio da cadeia -CH2-CH2-N
    patt = Chem.MolFromSmarts("c[CH2][CH2][NX3]")
    matches = mol.GetSubstructMatches(patt)
    if not matches:
        raise ValueError("nitrogênio da etilamina não encontrado")
    n_idx = matches[0][-1]
    n = mol.GetAtomWithIdx(n_idx)
    heavy = n.GetDegree()
    if heavy == 1:
        why = ("N primario (R3=H): linker entra direto no lugar de um H; "
               "serie preferencial para conjugacao")
    else:
        why = ("N ja substituido por R3: o linker ocupa o H restante e o N vira "
               "terciario — confirme que R3 nao e essencial para a ligacao a PCSK9")
    return n_idx, why


def reorder_attachment_first(mol, attach_idx: int):
    """Renumera para que o átomo de conjugação seja o índice 0.

    Isso é requisito de assemble_full_protac() no notebook: RDKit
    ReplaceSubstructs(..., replacementConnectionPoint=0) liga o linker ao
    átomo 0 da molécula de substituição.
    """
    order = [attach_idx] + [i for i in range(mol.GetNumAtoms()) if i != attach_idx]
    return Chem.RenumberAtoms(mol, order)


def make_dummy_variant(mol):
    """Versão com átomo dummy [*] no N, para casar com CAP_SMARTS='[#0]' do WP2."""
    rw = Chem.RWMol(mol)
    attach_idx, _ = find_exit_vector(rw)
    dummy = Chem.Atom(0)
    d_idx = rw.AddAtom(dummy)
    rw.AddBond(attach_idx, d_idx, Chem.BondType.SINGLE)
    n = rw.GetAtomWithIdx(attach_idx)
    if n.GetNumExplicitHs() > 0:
        n.SetNumExplicitHs(n.GetNumExplicitHs() - 1)
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return out


# ---------------------------------------------------------------------------
# 4. 3D
# ---------------------------------------------------------------------------
def embed_best_conformer(mol, n_confs: int, seed: int = 0xC0FFEE):
    """Gera n confôrmeros (ETKDGv3), otimiza em MMFF94s e devolve o de menor energia."""
    molh = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.useSmallRingTorsions = True
    cids = AllChem.EmbedMultipleConfs(molh, numConfs=n_confs, params=params)
    if not cids:
        return None, None
    # MMFF94s não tipa átomos dummy ([*]); nesses casos cai para UFF.
    if AllChem.MMFFHasAllMoleculeParams(molh):
        res = AllChem.MMFFOptimizeMoleculeConfs(molh, maxIters=2000, mmffVariant="MMFF94s")
    else:
        res = AllChem.UFFOptimizeMoleculeConfs(molh, maxIters=2000)
    energies = [(e, cid) for cid, (conv, e) in zip(cids, res)]
    energies.sort()
    best_e, best_cid = energies[0]
    # mantém apenas o melhor confôrmero
    keep = Chem.Mol(molh)
    keep.RemoveAllConformers()
    keep.AddConformer(molh.GetConformer(best_cid), assignId=True)
    return keep, best_e


# ---------------------------------------------------------------------------
# 5. main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", type=Path, required=True,
                    help="pasta de saída (será criada)")
    ap.add_argument("--n-confs", type=int, default=50,
                    help="confôrmeros por molécula na busca 3D (default 50)")
    ap.add_argument("--no-3d", action="store_true",
                    help="pula a geração 3D (só CSV/SMI, rápido)")
    ap.add_argument("--only-figure", action="store_true",
                    help="apenas os R explicitamente listados na figura "
                         "(sem as variações F/Cl/OH dos quadros)")
    ap.add_argument("--seed", type=int, default=0xC0FFEE)
    args = ap.parse_args()

    out = args.outdir.expanduser()
    (out / "sdf").mkdir(parents=True, exist_ok=True)
    (out / "sdf_dummy").mkdir(parents=True, exist_ok=True)

    print(f"Enumerando warheads PCSK9 (scaffold feniletilamina)...")
    records = enumerate_warheads(only_figure=args.only_figure)
    print(f"  {len(records)} moléculas únicas após deduplicação por SMILES canônico")

    rows = []
    writer_all = None if args.no_3d else Chem.SDWriter(str(out / "pcsk9_warheads_3d.sdf"))
    n_3d_ok = 0

    for rec in records:
        mol = Chem.MolFromSmiles(rec["smiles"])
        desc = describe(mol)
        attach_idx, attach_why = find_exit_vector(mol)
        mol_ord = reorder_attachment_first(mol, attach_idx)

        row = dict(rec)
        row.update(desc)
        row["filter_flags"] = flag_filters(desc)
        row["net_charge_ph74"] = charge_at_ph74(mol)
        row["exit_vector_atom_idx0"] = 0            # após renumeração
        row["exit_vector_atom_serial_1based"] = 1   # PRosettaC "Anchor atoms"
        row["exit_vector_note"] = attach_why
        row["smiles_reordered"] = Chem.MolToSmiles(mol_ord, canonical=False)
        row["smiles_dummy"] = Chem.MolToSmiles(make_dummy_variant(mol))

        if not args.no_3d:
            mol3d, energy = embed_best_conformer(mol_ord, args.n_confs, args.seed)
            if mol3d is None:
                print(f"  [aviso] embedding 3D falhou para {rec['warhead_id']}", file=sys.stderr)
                row["mmff_energy_kcal"] = ""
            else:
                n_3d_ok += 1
                row["mmff_energy_kcal"] = round(energy, 2)
                for k, v in row.items():
                    mol3d.SetProp(str(k), str(v))
                mol3d.SetProp("_Name", rec["warhead_id"])
                writer_all.write(mol3d)
                with Chem.SDWriter(str(out / "sdf" / f"{rec['warhead_id']}.sdf")) as w:
                    w.write(mol3d)

                # O átomo dummy [*] não tem tipo nem em MMFF nem em UFF, então
                # o RDKit despeja um UFFTYPER por confôrmero. É esperado e não
                # afeta a geometria dos átomos reais — silenciado para não
                # afogar avisos que importam.
                try:
                    from rdkit import rdBase
                    _silencio = rdBase.BlockLogs()
                except (ImportError, AttributeError):
                    _silencio = None
                dummy3d, _ = embed_best_conformer(
                    make_dummy_variant(mol_ord), max(10, args.n_confs // 5), args.seed)
                del _silencio
                if dummy3d is not None:
                    dummy3d.SetProp("_Name", f"{rec['warhead_id']}_dummy")
                    dummy3d.SetProp("warhead_id", rec["warhead_id"])
                    with Chem.SDWriter(str(out / "sdf_dummy" /
                                           f"{rec['warhead_id']}_dummy.sdf")) as w:
                        w.write(dummy3d)

        rows.append(row)

    if writer_all is not None:
        writer_all.close()

    # --- CSV ---
    csv_path = out / "pcsk9_warheads.csv"
    cols = list(rows[0].keys())
    with open(csv_path, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=cols)
        wr.writeheader()
        wr.writerows(rows)

    # --- SMI ---
    smi_path = out / "pcsk9_warheads.smi"
    with open(smi_path, "w") as fh:
        for r in rows:
            fh.write(f"{r['smiles']}\t{r['warhead_id']}\n")

    # --- script auxiliar de PDBQT (Meeko, env pf_vs) ---
    sh = out / "prepare_pdbqt.sh"
    sh.write_text(
        "#!/usr/bin/env bash\n"
        "# Converte os warheads 3D em PDBQT para docking (Vina/Uni-Dock).\n"
        "# Rodar no env pf_vs:  conda activate pf_vs && bash prepare_pdbqt.sh\n"
        "set -euo pipefail\n"
        f'OUT="{out}/pdbqt"\n'
        'mkdir -p "$OUT"\n'
        f'for f in "{out}"/sdf/*.sdf; do\n'
        '  base=$(basename "$f" .sdf)\n'
        '  mk_prepare_ligand.py -i "$f" -o "$OUT/$base.pdbqt"\n'
        'done\n'
        'echo "PDBQT em $OUT"\n'
    )
    sh.chmod(0o755)

    # --- resumo ---
    print(f"\nSaídas em {out}:")
    print(f"  pcsk9_warheads.csv        {len(rows)} linhas")
    print(f"  pcsk9_warheads.smi")
    if not args.no_3d:
        print(f"  pcsk9_warheads_3d.sdf     {n_3d_ok} confôrmeros 3D")
        print(f"  sdf/                      {n_3d_ok} arquivos (Heads do PRosettaC)")
        print(f"  sdf_dummy/                warheads com [*] no ponto de conjugação")
    print(f"  prepare_pdbqt.sh          (rodar no env pf_vs)")

    flagged = [r for r in rows if r["filter_flags"] != "ok"]
    print(f"\n{len(rows) - len(flagged)}/{len(rows)} dentro da faixa típica de warhead "
          f"{FILTERS}")
    if flagged:
        print(f"  {len(flagged)} marcados (coluna filter_flags) — não foram descartados")
    charged = [r for r in rows if r["net_charge_ph74"] != 0]
    print(f"{len(charged)}/{len(rows)} com carga líquida ≠ 0 em pH 7,4 "
          f"(use em `acpype -n` e no genion)")


if __name__ == "__main__":
    main()
