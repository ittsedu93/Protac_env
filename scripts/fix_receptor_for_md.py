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


def _tentar_chimerax(receptor: Path, alvos, out_pdb: Path, comando, rotulo: str):
    """Roda um script ChimeraX e devolve (gerou_arquivo, log)."""
    script = out_pdb.parent / f"fix_receptor_{rotulo}.cxc"
    cmds = [f"open {receptor}"] + [comando(r) for r in alvos]
    cmds += [f"save {out_pdb}", "exit"]
    script.write_text("\n".join(cmds) + "\n")
    if out_pdb.exists():
        out_pdb.unlink()
    r = subprocess.run([CHIMERAX_EXE, "--nogui", "--exit", str(script)],
                       capture_output=True, text=True)
    return out_pdb.exists(), (r.stdout or "") + (r.stderr or "")


def reconstruir(receptor: Path, alvos, out_pdb: Path):
    """Reconstrói as cadeias laterais, tentando as rotas disponíveis em ordem.

    A sintaxe do `swapaa` do ChimeraX aceita `criteria` como uma sequência de
    LETRAS (d, c, h, p) ou um número de rotâmero. Passar "highest" faz o
    parser ler 'h', depois 'i', e falhar com "Unknown criteria: 'i'" — foi o
    que derrubou a primeira versão. Sem o argumento, ele usa o critério
    padrão, que é o que se quer.

    Cada rota é verificada pelo mesmo detector que apontou o problema: só
    conta como sucesso se os resíduos deixarem de estar incompletos.
    """
    rotas = [
        ("swapaa padrão",
         lambda r: f"swapaa /{r['cadeia']}:{r['num']} {r['tipo']}"),
        ("swapaa criteria dchp",
         lambda r: f"swapaa /{r['cadeia']}:{r['num']} {r['tipo']} criteria dchp"),
        ("swapaa rotâmero 1",
         lambda r: f"swapaa /{r['cadeia']}:{r['num']} {r['tipo']} criteria 1"),
    ]
    for rotulo, comando in rotas:
        if not Path(CHIMERAX_EXE).exists():
            break
        ok, log = _tentar_chimerax(receptor, alvos, out_pdb, comando, rotulo)
        if ok:
            restantes = incompletos(out_pdb)
            print(f"    [{rotulo}] {len(alvos)} -> {len(restantes)} incompletos")
            if not restantes:
                return out_pdb, []
            if len(restantes) < len(alvos):
                return out_pdb, restantes     # progresso parcial já ajuda
        else:
            erro = [l for l in log.splitlines() if "ERROR" in l or "Unknown" in l]
            print(f"    [{rotulo}] não gerou o arquivo"
                  + (f": {erro[-1][:80]}" if erro else ""))

    # PDBFixer, se estiver disponível neste env
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
        print("    [pdbfixer] reconstruindo átomos faltantes")
        fx = PDBFixer(filename=str(receptor))
        fx.findMissingResidues()
        fx.findMissingAtoms()
        fx.addMissingAtoms()
        with open(out_pdb, "w") as fh:
            PDBFile.writeFile(fx.topology, fx.positions, fh, keepIds=True)
        restantes = incompletos(out_pdb)
        print(f"    [pdbfixer] {len(alvos)} -> {len(restantes)} incompletos")
        return out_pdb, restantes
    except ImportError:
        print("    [pdbfixer] não instalado neste env")
    except Exception as exc:
        print(f"    [pdbfixer] falhou: {type(exc).__name__}: {exc}")

    return None, alvos


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
    ap.add_argument("--truncar-perto", action="store_true",
                    help="aceita truncar resíduo incompleto do sítio; é uma "
                         "mutação, e precisa constar na descrição do sistema")
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

    print(f"\n[1] Reconstruindo as cadeias laterais")
    gerado, faltas = reconstruir(rec, faltas, out)
    if gerado and not faltas:
        print(f"    todos reconstruídos: {out}")
        return

    for r in faltas:
        if "dist" not in r:
            r["dist"] = dist_ao_ligante(r, lig)
    perto = [r for r in faltas if r["dist"] < args.dist_min_truncar]

    print(f"\n[2] Restam {len(faltas)} incompletos "
          f"({len(perto)} a menos de {args.dist_min_truncar} Å do ligante)")

    if perto and not args.truncar_perto:
        print("\n    Resíduos incompletos no sítio:")
        for r in perto:
            print(f"      {r['cadeia']}:{r['tipo']}{r['num']} "
                  f"({r['dist']:.1f} Å) falta "
                  f"{r['falta_sidechain'] + r['falta_backbone']}")
        print("\n    Truncar um resíduo do sítio muda a interação que a MD")
        print("    existe para medir, então isso não é feito por conta própria.")
        print("\n    Três saídas, em ordem de preferência:")
        print("      1. reconstruir à mão no ChimeraX e salvar por cima:")
        for r in perto:
            print(f"           swapaa /{r['cadeia']}:{r['num']} {r['tipo']}")
        print(f"           save {out}")
        print("      2. usar outra estrutura do receptor, com o sítio completo")
        print("      3. aceitar a truncagem conscientemente, com")
        print("           --truncar-perto  (registre na descrição do sistema)")
        raise SystemExit("\nreceptor não pôde ser completado no sítio")

    fonte = out if (gerado and out.exists()) else rec
    truncar_para_ala(fonte, faltas, out)
    print(f"\n    {len(faltas)} resíduos truncados para ALA:")
    for r in faltas:
        marca = "  <-- NO SÍTIO" if r["dist"] < args.dist_min_truncar else ""
        print(f"      {r['cadeia']}:{r['tipo']}{r['num']} -> ALA "
              f"({r['dist']:.1f} Å){marca}")
    print("\n    São MUTAÇÕES. Registre-as na descrição do sistema.")

    sobrando = incompletos(out)
    if sobrando:
        raise SystemExit(f"ainda restam {len(sobrando)} incompletos após a "
                         f"truncagem: {[(r['tipo'], r['num']) for r in sobrando[:5]]}")
    print(f"\n{out}")


if __name__ == "__main__":
    main()
