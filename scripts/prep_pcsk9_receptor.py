#!/usr/bin/env python
"""
Etapa 3.0a — preparo do receptor PCSK9 e definição do sítio, a partir do
co-cristal 6U26 (PCSK9 + composto 16, 1,53 Å, ligante HET `063`).

Roda no env `mdtools` (numpy; ChimeraX e obabel como executáveis). Não precisa
do Vina.

Por que este script existe
--------------------------
Os warheads da Figura 3 são moléculas novas, sem pose conhecida na PCSK9, e o
PRosettaC exige `Heads` já posicionados. O co-cristal 6U26 resolve as três
coisas de que o docking precisa:

  1. ONDE dockar — o sítio é o bolsão ocupado pelo ligante `063`, medido da
     estrutura, não estimado por detecção de cavidade.
  2. COMO VALIDAR — há ligante co-cristalizado, então o redocking com corte
     RMSD <= 2,0 Å da metodologia do WP1 passa a ser aplicável aqui também.
  3. PARA ONDE SAI O LINKER — o braço PEG/guanidina do `063` está totalmente
     exposto ao solvente e define experimentalmente o vetor de saída.

Modos de sítio
--------------
`--site-mode ligand` (default)  box a partir do ligante co-cristalizado.
    Usa só os átomos ENTERRADOS do ligante (>= `--burial-threshold` vizinhos
    proteicos a 6 Å), porque o `063` tem um braço de ~14 Å no solvente que, se
    entrasse na conta, inflaria o box para ~28 Å e transformaria o docking em
    busca cega.
`--site-mode interface`  box na interface com uma cadeia parceira (ex.: EGF(A)
    do LDLR num complexo como 3BPS).
`--site-mode residues`   box em volta de uma lista de resíduos.

Saídas (em --outdir)
--------------------
    chains.txt              composição de cadeias
    receptor/*_receptor.pdb/.pdbqt    receptor limpo e protonado
    receptor/ref_ligand.pdb/.pdbqt/.sdf   ligante de referência (redocking)
    pcsk9_site.json         centro, box, vetor de saída e proveniência

Uso:
    python prep_pcsk9_receptor.py --pdb-file 6U26.pdb --outdir ~/PCSK9_docking --inspect-only
    python prep_pcsk9_receptor.py --pdb-file 6U26.pdb --outdir ~/PCSK9_docking \\
        --target-chain B --keep-chains A B --ref-ligand-resname 063
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import urllib.request
from collections import Counter, OrderedDict
from pathlib import Path

import numpy as np

CHIMERAX_EXE = "/usr/bin/chimerax"
OBABEL_EXE = "/home/soberano/miniconda3/envs/obabel_env/bin/obabel"

DEFAULT_INTERFACE_RESIDUES = [153, 154, 155, 194, 238, 239, 367, 369,
                              372, 374, 375, 377, 378, 379]
SOLVENT = {"HOH", "WAT", "DOD"}
IONS = {"NA", "CL", "MG", "ZN", "CA", "K", "SO4", "PO4", "EDO", "GOL", "PEG"}


# ---------------------------------------------------------------------------
def fetch_pdb(pdb_id: str, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{pdb_id.upper()}.pdb"
    if not out_path.exists():
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        print(f"  baixando {url}")
        urllib.request.urlretrieve(url, out_path)
    return out_path


def parse_atoms(pdb_path: Path):
    atoms = []
    for line in Path(pdb_path).read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            atoms.append(dict(
                record=line[:6].strip(),
                name=line[12:16].strip(),
                resname=line[17:20].strip(),
                chain=line[21],
                resid=int(line[22:26]),
                xyz=(float(line[30:38]), float(line[38:46]), float(line[46:54])),
                raw=line,
            ))
        except ValueError:
            continue
    return atoms


def describe_chains(pdb_path: Path, out_txt: Path) -> str:
    atoms = parse_atoms(pdb_path)
    prot, het = OrderedDict(), Counter()
    het_loc = OrderedDict()
    for a in atoms:
        if a["record"] == "ATOM":
            prot.setdefault(a["chain"], set()).add(a["resid"])
        elif a["resname"] not in SOLVENT:
            het[a["resname"]] += 1
            het_loc.setdefault(a["resname"], set()).add(f"{a['chain']}{a['resid']}")

    lines = [f"Estrutura: {pdb_path.name}", "", "Cadeias de proteína:"]
    for c, res in prot.items():
        lines.append(f"  cadeia {c}: {len(res):4d} residuos, "
                     f"numeracao {min(res)}-{max(res)}")
    lines += ["", "Heteroátomos (sem água):"]
    if het:
        for rn, n in het.most_common():
            tag = " [ION/CRIOPROTETOR]" if rn in IONS else ""
            lines.append(f"  {rn}: {n} atomos em {sorted(het_loc[rn])}{tag}")
    else:
        lines.append("  nenhum")
    lines += ["",
              "Na PCSK9 madura o pró-domínio fica ~61-152 e o catalítico+CHRD",
              "~153-692. Em 6U26: cadeia A = pró-domínio, cadeia B = catalítico",
              "+CHRD, ligante 063 = composto 16, ancorado na cadeia B."]
    text = "\n".join(lines)
    out_txt.write_text(text + "\n")
    return text


# ---------------------------------------------------------------------------
# Sítio
# ---------------------------------------------------------------------------
def _burial(lig_xyz: np.ndarray, prot_xyz: np.ndarray, radius: float = 6.0):
    d = np.linalg.norm(lig_xyz[:, None, :] - prot_xyz[None, :, :], axis=2)
    return (d <= radius).sum(axis=1), d


def site_from_ligand(pdb_path: Path, resname: str, pad: float,
                     burial_threshold: int, max_edge: float):
    """Box sobre a porção ENTERRADA do ligante co-cristalizado, mais o vetor de
    saída apontado pela porção exposta ao solvente."""
    atoms = parse_atoms(pdb_path)
    lig = [a for a in atoms if a["resname"] == resname]
    prot = [a for a in atoms if a["record"] == "ATOM"]
    if not lig:
        raise SystemExit(f"ligante '{resname}' não encontrado em {pdb_path.name}")

    L = np.array([a["xyz"] for a in lig])
    P = np.array([a["xyz"] for a in prot])
    nb, d = _burial(L, P)

    buried = nb >= burial_threshold
    if buried.sum() < 3:
        raise SystemExit(
            f"só {buried.sum()} átomos com >= {burial_threshold} vizinhos; "
            f"baixe --burial-threshold (máximo observado: {nb.max()})")

    core = L[buried]
    center = core.mean(axis=0)
    size = (core.max(axis=0) - core.min(axis=0)) + 2 * pad
    size = np.minimum(size, max_edge)          # trava a aresta máxima

    # vetor de saída: do núcleo enterrado para os átomos sem nenhum vizinho
    exposed = nb == 0
    exit_vec, exit_len, exposed_names = None, None, []
    if exposed.sum() >= 2:
        ec = L[exposed].mean(axis=0)
        v = ec - center
        exit_len = float(np.linalg.norm(v))
        exit_vec = (v / exit_len).tolist()
        exposed_names = [lig[i]["name"] for i in np.where(exposed)[0]]

    contacts = sorted({(a["chain"], a["resid"], a["resname"])
                       for i, a in enumerate(prot) if d[:, i].min() <= 4.5},
                      key=lambda t: (t[0], t[1]))

    return dict(
        center=[round(float(x), 3) for x in center],
        box_size=[round(float(x), 2) for x in size],
        ref_ligand_resname=resname,
        ref_ligand_n_atoms=len(lig),
        n_buried_atoms=int(buried.sum()),
        burial_threshold=burial_threshold,
        exit_vector=[round(x, 4) for x in exit_vec] if exit_vec else None,
        exit_vector_origin=[round(float(x), 3) for x in center],
        exit_arm_length_A=round(exit_len, 1) if exit_len else None,
        exit_arm_atoms=exposed_names,
        contact_residues=[f"{c}:{rn}{ri}" for c, ri, rn in contacts],
        provenance=(f"atomos do ligante {resname} com >= {burial_threshold} "
                    f"vizinhos proteicos a 6 A em {pdb_path.name}, "
                    f"padding {pad} A, aresta limitada a {max_edge} A"),
    )


def site_from_interface(pdb_path, target_chain, partner_chain, cutoff, pad):
    atoms = parse_atoms(pdb_path)
    tgt = [a for a in atoms if a["record"] == "ATOM" and a["chain"] == target_chain]
    prt = [a for a in atoms if a["record"] == "ATOM" and a["chain"] == partner_chain]
    if not tgt or not prt:
        raise SystemExit("cadeia alvo ou parceira não encontrada")
    tc = np.array([a["xyz"] for a in tgt])
    pc = np.array([a["xyz"] for a in prt])
    mask = np.linalg.norm(tc[:, None] - pc[None], axis=2).min(axis=1) <= cutoff
    if not mask.any():
        raise SystemExit(f"nenhum contato a <= {cutoff} Å entre as cadeias")
    sel = tc[mask]
    return dict(
        center=[round(float(x), 3) for x in sel.mean(axis=0)],
        box_size=[round(float(x), 2) for x in (sel.max(0) - sel.min(0)) + 2 * pad],
        interface_residues=sorted({tgt[i]["resid"] for i in np.where(mask)[0]}),
        provenance=(f"interface cadeia {target_chain} <= {cutoff} A da cadeia "
                    f"{partner_chain} em {pdb_path.name}"),
    )


def site_from_residues(pdb_path, target_chain, residues, pad):
    atoms = parse_atoms(pdb_path)
    sel = [a for a in atoms if a["record"] == "ATOM"
           and a["chain"] == target_chain and a["resid"] in set(residues)]
    if not sel:
        raise SystemExit(f"nenhum dos resíduos encontrado na cadeia {target_chain}")
    found = sorted({a["resid"] for a in sel})
    missing = sorted(set(residues) - set(found))
    if missing:
        print(f"  [aviso] resíduos ausentes: {missing}", file=sys.stderr)
    c = np.array([a["xyz"] for a in sel])
    return dict(
        center=[round(float(x), 3) for x in c.mean(axis=0)],
        box_size=[round(float(x), 2) for x in (c.max(0) - c.min(0)) + 2 * pad],
        interface_residues=found, residues_missing=missing,
        provenance=f"residuos {found} da cadeia {target_chain} em {pdb_path.name}",
    )


def sanity_check_box(site: dict):
    e = site["box_size"]
    if max(e) > 30.0:
        print(f"  [ATENÇÃO] aresta {max(e):.0f} Å > 30 Å: o Vina vai varrer quase")
        print("            a superfície toda. Suba --burial-threshold.")
    if min(e) < 14.0:
        print(f"  [ATENÇÃO] aresta {min(e):.0f} Å < 14 Å: pode não caber os")
        print("            warheads maiores (série naftaleno, ~400 Da).")


# ---------------------------------------------------------------------------
# Receptor e ligante de referência
# ---------------------------------------------------------------------------
def extract_ref_ligand(pdb_path: Path, resname: str, out_dir: Path):
    """Extrai o ligante co-cristalizado. Vai servir de referência do redocking,
    então a mesma perceção de ligações do obabel é usada nos dois lados (pose e
    referência) — é o que garante que o GetBestRMS compare grafos iguais."""
    out_dir.mkdir(parents=True, exist_ok=True)
    lig_pdb = out_dir / "ref_ligand.pdb"
    lines = [a["raw"] for a in parse_atoms(pdb_path) if a["resname"] == resname]
    if not lines:
        raise SystemExit(f"ligante {resname} não encontrado para extração")
    lig_pdb.write_text("\n".join(lines) + "\nEND\n")

    lig_sdf = out_dir / "ref_ligand.sdf"
    lig_pdbqt = out_dir / "ref_ligand.pdbqt"
    subprocess.run([OBABEL_EXE, str(lig_pdb), "-O", str(lig_sdf), "-h"], check=True)
    subprocess.run([OBABEL_EXE, str(lig_sdf), "-O", str(lig_pdbqt)], check=True)
    print(f"    ligante de referência: {lig_sdf.name}, {lig_pdbqt.name}")
    return lig_pdb, lig_sdf, lig_pdbqt


def prep_receptor(pdb_path: Path, keep_chains: list[str], out_dir: Path,
                  drop_resnames: list[str]):
    """Remove cadeias extras/solvente/ligante, protona e converte para PDBQT.

    Sintaxe do ChimeraX: cadeias vão num único `/` separadas por vírgula
    (`/A,B`), e a negação é `~/A,B`. A forma `~(/A,/B)` é rejeitada com
    "invalid atoms specifier". O resíduo é selecionado por `::name="063"` —
    a própria forma que o ChimeraX sugere no log — porque `:063` começando
    com dígito é ambíguo com número de resíduo.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    receptor_pdb = out_dir / f"{pdb_path.stem}_receptor.pdb"
    receptor_pdbqt = receptor_pdb.with_suffix(".pdbqt")
    script = out_dir / f"{pdb_path.stem}_prep.cxc"

    todas = sorted({a["chain"] for a in parse_atoms(pdb_path)
                    if a["record"] == "ATOM"})
    remover = [c for c in todas if c not in keep_chains]

    cmds = [f"open {pdb_path}"]
    if remover:                      # spec vazio também é erro no ChimeraX
        cmds.append(f"delete ~/{','.join(keep_chains)}")
    cmds += ["delete solvent", "delete ions"]
    for rn in drop_resnames:
        cmds.append(f'delete ::name="{rn}"')
    cmds += ["addh", f"save {receptor_pdb}", "exit"]
    script.write_text("\n".join(cmds) + "\n")

    print(f"    ChimeraX: mantendo {keep_chains}"
          + (f", removendo cadeias {remover}" if remover else "")
          + f", removendo {drop_resnames}")

    # --exit garante que ele saia mesmo se um comando falhar, em vez de ficar
    # parado no prompt cmd> esperando entrada que não vem.
    r = subprocess.run([CHIMERAX_EXE, "--nogui", "--exit", str(script)],
                       capture_output=True, text=True)

    if not receptor_pdb.exists():
        saida = (r.stdout or "") + (r.stderr or "")
        raise SystemExit(
            f"ChimeraX não gerou {receptor_pdb}.\n"
            f"Script: {script}\n"
            f"Últimas linhas do ChimeraX:\n"
            + "\n".join(saida.strip().splitlines()[-25:]))

    subprocess.run([OBABEL_EXE, str(receptor_pdb), "-O", str(receptor_pdbqt),
                    "-xr", "-p", "7.4"], check=True)
    print(f"    receptor: {receptor_pdb.name}, {receptor_pdbqt.name}")
    return receptor_pdb, receptor_pdbqt


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--pdb-file", type=Path, help="PDB local (ex.: 6U26_1.pdb)")
    src.add_argument("--pdb-id", help="baixa do RCSB")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--inspect-only", action="store_true")
    ap.add_argument("--target-chain", default="B")
    ap.add_argument("--keep-chains", nargs="*", default=None,
                    help="cadeias no receptor (default: A e B, pró-domínio junto)")
    ap.add_argument("--site-mode", choices=["ligand", "interface", "residues"],
                    default="ligand")
    ap.add_argument("--ref-ligand-resname", default="063")
    ap.add_argument("--burial-threshold", type=int, default=20,
                    help="vizinhos proteicos a 6 Å para um átomo contar como "
                         "enterrado (default 20)")
    ap.add_argument("--partner-chain")
    ap.add_argument("--interface-cutoff", type=float, default=5.0)
    ap.add_argument("--residues", type=int, nargs="*",
                    default=DEFAULT_INTERFACE_RESIDUES)
    ap.add_argument("--box-pad", type=float, default=4.0)
    ap.add_argument("--max-box-edge", type=float, default=24.0)
    ap.add_argument("--skip-receptor-prep", action="store_true")
    args = ap.parse_args()

    out = args.outdir.expanduser()
    out.mkdir(parents=True, exist_ok=True)

    print("[1/4] Estrutura")
    if args.pdb_file:
        src_path = args.pdb_file.expanduser()
        if not src_path.exists():
            raise SystemExit(f"não encontrei {src_path}")
        pdb_path = out / "pdb" / src_path.name
        pdb_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(src_path, pdb_path)
        print(f"    {pdb_path}")
    else:
        pdb_path = fetch_pdb(args.pdb_id, out / "pdb")

    print("[2/4] Composição de cadeias")
    report = describe_chains(pdb_path, out / "chains.txt")
    print("\n".join("    " + l for l in report.splitlines()))

    if args.inspect_only:
        print("\n--inspect-only: parando aqui.")
        return

    print(f"\n[3/4] Sítio (modo {args.site_mode})")
    if args.site_mode == "ligand":
        site = site_from_ligand(pdb_path, args.ref_ligand_resname, args.box_pad,
                                args.burial_threshold, args.max_box_edge)
    elif args.site_mode == "interface":
        if not args.partner_chain:
            raise SystemExit("--site-mode interface exige --partner-chain")
        site = site_from_interface(pdb_path, args.target_chain,
                                   args.partner_chain, args.interface_cutoff,
                                   args.box_pad)
    else:
        site = site_from_residues(pdb_path, args.target_chain,
                                  args.residues, args.box_pad)

    site["source_structure"] = pdb_path.name
    site["target_chain"] = args.target_chain
    print(f"    centro : {site['center']}")
    print(f"    box (Å): {site['box_size']}")
    if site.get("exit_vector"):
        print(f"    vetor de saída : {site['exit_vector']} "
              f"(braço exposto de {site['exit_arm_length_A']} Å)")
    if site.get("contact_residues"):
        print(f"    {len(site['contact_residues'])} resíduos de contato: "
              f"{', '.join(site['contact_residues'][:8])}...")
    sanity_check_box(site)

    print("\n[4/4] Receptor e ligante de referência")
    if args.skip_receptor_prep:
        print("    --skip-receptor-prep: pulando")
    else:
        keep = args.keep_chains if args.keep_chains is not None else ["A", "B"]
        rec_dir = out / "receptor"
        if args.site_mode == "ligand":
            _, lig_sdf, lig_pdbqt = extract_ref_ligand(
                pdb_path, args.ref_ligand_resname, rec_dir)
            site["ref_ligand_sdf"] = str(lig_sdf)
            site["ref_ligand_pdbqt"] = str(lig_pdbqt)
        rec_pdb, rec_pdbqt = prep_receptor(
            pdb_path, keep, rec_dir, [args.ref_ligand_resname])
        site["receptor_pdb"] = str(rec_pdb)
        site["receptor_pdbqt"] = str(rec_pdbqt)

    site_json = out / "pcsk9_site.json"
    site_json.write_text(json.dumps(site, indent=2))
    print(f"\nSítio salvo em {site_json}")
    print("\nConfira o box no ChimeraX antes de dockar a série inteira:")
    print(f"  chimerax {site.get('receptor_pdb', '<receptor>')} "
          f"{site.get('ref_ligand_sdf', '')}")
    print("\nDepois, no env pf_vs:")
    print(f"  python scripts/dock_warheads_pcsk9.py --site {site_json} \\")
    print(f"      --warheads-pdbqt ~/PCSK9_warheads/pdbqt \\")
    print(f"      --outdir {out}/docking --validate-only")


if __name__ == "__main__":
    main()
