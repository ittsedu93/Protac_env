#!/usr/bin/env python
"""
Fase 8a-bis — completa resíduos incompletos do receptor antes do pdb2gmx.

Roda no env `mdtools`. Segundos a minutos.

Estruturas cristalográficas têm cadeias laterais parcialmente resolvidas:
o átomo existe na proteína mas não aparece no mapa de densidade, e o PDB o
omite. Para o docking isso raramente importa; para a MD, importa sempre. O
GROMACS recusa a estrutura com mensagens como

    Fatal error: Incomplete ring in HIS68

que vêm da análise de protonação das histidinas — um anel imidazólico com
átomo faltando não pode ser protonado.

Duas rotas, nesta ordem:
  1. **ChimeraX `swapaa`** reconstrói a cadeia lateral no mesmo tipo de
     resíduo, usando rotâmero de biblioteca. É o conserto de verdade.
  2. **Truncar para alanina** quando a reconstrução falha. É uma mutação, e
     por isso só se aplica a resíduos longe do sítio — o script recusa
     truncar resíduo próximo do ligante e diz por quê.

    python fix_receptor_for_md.py --receptor 4TZ4_receptor.pdb \\
        --ligante protac.sdf --out 4TZ4_receptor_fixed.pdb
"""

from __future__ import annotations

import argparse
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np

CHIMERAX_EXE = "/usr/bin/chimerax"

# átomos pesados esperados na cadeia lateral de cada resíduo
SIDECHAIN = {
    "ARG": "CB CG CD NE CZ NH1 NH2", "ASN": "CB CG OD1 ND2",
    "ASP": "CB CG OD1 OD2", "CYS": "CB SG", "GLN": "CB CG CD OE1 NE2",
    "GLU": "CB CG CD OE1 OE2", "HIS": "CB CG ND1 CD2 CE1 NE2",
    "ILE": "CB CG1 CG2 CD1", "LEU": "CB CG CD1 CD2", "LYS": "CB CG CD CE NZ",
    "MET": "CB CG SD CE", "PHE": "CB CG CD1 CD2 CE1 CE2 CZ",
    "PRO": "CB CG CD", "SER": "CB OG", "THR": "CB OG1 CG2",
    "TRP": "CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2",
    "TYR": "CB CG CD1 CD2 CE1 CE2 CZ OH", "VAL": "CB CG1 CG2",
    "ALA": "CB", "GLY": "",
}
BACKBONE = {"N", "CA", "C", "O"}


def ler_residuos(pdb: Path):
    res = defaultdict(lambda: {"atomos": set(), "xyz": []})
    for l in Path(pdb).read_text().splitlines():
        if not l.startswith("ATOM"):
            continue
        nome = l[12:16].strip()
        if nome.startswith("H") or l[76:78].strip() == "H":
            continue
        chave = (l[21], int(l[22:26]), l[17:20].strip())
        res[chave]["atomos"].add(nome)
        res[chave]["xyz"].append((float(l[30:38]), float(l[38:46]),
                                  float(l[46:54])))
    return res


def incompletos(pdb: Path):
    faltas = []
    for (cad, num, tipo), dados in ler_residuos(pdb).items():
        esperado = SIDECHAIN.get(tipo)
        if esperado is None:
            continue
        falta_sc = set(esperado.split()) - dados["atomos"]
        falta_bb = BACKBONE - dados["atomos"]
        if falta_sc or falta_bb:
            faltas.append({"cadeia": cad, "num": num, "tipo": tipo,
                           "falta_sidechain": sorted(falta_sc),
                           "falta_backbone": sorted(falta_bb),
                           "xyz": np.array(dados["xyz"])})
    return sorted(faltas, key=lambda d: (d["cadeia"], d["num"]))


def dist_ao_ligante(res, lig_xyz):
    if lig_xyz is None or len(lig_xyz) == 0:
        return float("inf")
    return float(np.linalg.norm(
        res["xyz"][:, None, :] - lig_xyz[None, :, :], axis=2).min())


def ler_ligante_xyz(caminho: Path | None):
    if caminho is None or not Path(caminho).exists():
        return None
    try:
        from rdkit import Chem
        m = (Chem.MolFromMolFile(str(caminho))
             if str(caminho).endswith((".sdf", ".mol"))
             else Chem.MolFromPDBFile(str(caminho)))
        return m.GetConformer().GetPositions() if m else None
    except Exception:
        return None


def reconstruir_chimerax(receptor: Path, alvos, out_pdb: Path):
    """swapaa para o mesmo tipo de resíduo: recria a cadeia lateral."""
    script = out_pdb.parent / "fix_receptor.cxc"
    cmds = [f"open {receptor}"]
    for r in alvos:
        cmds.append(f"swapaa /{r['cadeia']}:{r['num']} {r['tipo']} "
                    f"criteria highest")
    cmds += [f"save {out_pdb}", "exit"]
    script.write_text("\n".join(cmds) + "\n")
    r = subprocess.run([CHIMERAX_EXE, "--nogui", "--exit", str(script)],
                       capture_output=True, text=True)
    return out_pdb.exists(), (r.stdout or "") + (r.stderr or "")


def truncar_para_ala(receptor: Path, alvos, out_pdb: Path):
    """Remove a cadeia lateral além do CB e renomeia para ALA.

    É uma MUTAÇÃO. Só se aplica a resíduos longe do sítio, e cada uma é
    registrada para constar na descrição do sistema.
    """
    manter = BACKBONE | {"CB"}
    chaves = {(r["cadeia"], r["num"]) for r in alvos}
    saida = []
    for l in Path(receptor).read_text().splitlines():
        if l.startswith("ATOM"):
            chave = (l[21], int(l[22:26]))
            if chave in chaves:
                if l[12:16].strip() not in manter:
                    continue
                l = l[:17] + "ALA" + l[20:]
        saida.append(l)
    out_pdb.write_text("\n".join(saida) + "\n")
    return out_pdb


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--receptor", type=Path, required=True)
    ap.add_argument("--ligante", type=Path,
                    help="para medir a distância ao sítio antes de truncar")
    ap.add_argument("--dist-min-truncar", type=float, default=8.0,
                    help="não trunca resíduo a menos desta distância do "
                         "ligante (Å, default 8)")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    rec = args.receptor.expanduser()
    out = args.out.expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)

    faltas = incompletos(rec)
    if not faltas:
        print("nenhum resíduo incompleto — copiando o receptor como está")
        out.write_text(rec.read_text())
        return

    print(f"{len(faltas)} resíduos incompletos em {rec.name}:")
    lig = ler_ligante_xyz(args.ligante)
    for r in faltas[:20]:
        d = dist_ao_ligante(r, lig)
        r["dist"] = d
        print(f"  {r['cadeia']}:{r['tipo']}{r['num']:<4d} falta "
              f"{r['falta_sidechain'] + r['falta_backbone']}"
              + (f"  ({d:.1f} Å do ligante)" if np.isfinite(d) else ""))
    if len(faltas) > 20:
        print(f"  ... e mais {len(faltas) - 20}")

    print(f"\n[1] ChimeraX swapaa — reconstruindo as cadeias laterais")
    ok, log = reconstruir_chimerax(rec, faltas, out)
    if ok:
        restantes = incompletos(out)
        if not restantes:
            print(f"    reconstruído: {out}")
            return
        print(f"    ainda restam {len(restantes)} incompletos")
        faltas = restantes
    else:
        print("    ChimeraX não gerou a estrutura:")
        print("\n".join("      " + l for l in log.strip().splitlines()[-8:]))

    print(f"\n[2] Truncando para alanina o que sobrou")
    for r in faltas:
        if "dist" not in r:
            r["dist"] = dist_ao_ligante(r, lig)
    perto = [r for r in faltas if r["dist"] < args.dist_min_truncar]
    if perto:
        print(f"    RECUSADO: {len(perto)} resíduos incompletos a menos de "
              f"{args.dist_min_truncar} Å do ligante:")
        for r in perto:
            print(f"      {r['cadeia']}:{r['tipo']}{r['num']} "
                  f"({r['dist']:.1f} Å)")
        raise SystemExit(
            "Truncar resíduo do sítio mudaria a interação que se quer medir. "
            "Reconstrua manualmente (ChimeraX: swapaa) ou escolha outra "
            "estrutura do receptor.")

    fonte = out if out.exists() else rec
    truncar_para_ala(fonte, faltas, out)
    print(f"    {len(faltas)} resíduos truncados para ALA (longe do sítio)")
    for r in faltas:
        print(f"      {r['cadeia']}:{r['tipo']}{r['num']} -> ALA "
              f"({r['dist']:.1f} Å do ligante)")
    print(f"\n    ATENÇÃO: são mutações. Registre-as na descrição do sistema.")
    print(f"\n{out}")


if __name__ == "__main__":
    main()
