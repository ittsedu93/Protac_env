#!/usr/bin/env python
"""
Fase 8c — analisa a MD e decide se o candidato passa.

Roda no env `mdtools` (MDAnalysis). Minutos.

    python md_analyze.py --md-dir ~/pipeline/md

Mede, por réplica e no conjunto:

  RMSD da proteína           o sistema é estável?
  RMSD do PROTAC             o ligante fica onde foi colocado?
  contatos recrutador–E3     a ancoragem na E3 ligase se mantém?
  contatos warhead–E3        interação ESPÚRIA: é o que o nível (ii) existe
                             para detectar. Se o warhead gruda na própria E3
                             ligase na ausência da PCSK9, a orientação do
                             linker é improdutiva e o candidato cai.
  RMSF por resíduo           que partes se mexem
  torções do linker          o linker explora conformações ou trava?

O veredito é uma regra explícita, não um número mágico — os cortes estão em
CRITERIOS e devem ser justificados na tese.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

CRITERIOS = {
    "rmsd_proteina_max_A": 3.5,      # acima disso o sistema não equilibrou
    "rmsd_protac_max_A": 5.0,        # o PROTAC é flexível; o corte é generoso
    "contatos_recrutador_min": 20,   # ancoragem mantida na média da trajetória
    "contatos_warhead_max": 15,      # acima disso há interação espúria com a E3
    "fracao_frames_ancorado_min": 0.7,
}


GMX = "/usr/local/gromacs/bin/gmx"


def topologia_legivel(md_dir: Path, rep: str, gmx: str = GMX):
    """Arquivo de topologia que o MDAnalysis consiga ler.

    O .tpr do GROMACS 2026 está na versão tpx 138, e o TPRParser do
    MDAnalysis para nela:

        ValueError: Your tpx version is 138, which this parser does not
        support, yet

    Isso não é problema dos dados: o .tpr é pedido apenas por causa dos nomes
    de átomo, nomes de resíduo e numeração, e tudo isso está igualmente num
    .gro do MESMO sistema. Ligações não são usadas em lugar nenhum desta
    análise. Tenta-se o .tpr primeiro para não mudar o comportamento onde ele
    funciona; não dando, usa-se um .gro, e em último caso pede-se ao próprio
    gmx que converta o .tpr — ele sempre entende o formato que ele mesmo
    escreveu.
    """
    import MDAnalysis as mda

    rep_dir = md_dir / rep
    tpr = rep_dir / "prod.tpr"

    if tpr.exists():
        try:
            mda.Universe(str(tpr))
            return tpr, "prod.tpr"
        except Exception as e:
            motivo = str(e).splitlines()[-1][:90]
            print(f"      o MDAnalysis não lê {rep}/prod.tpr ({motivo})")

    # um .gro do mesmo sistema tem a mesma ordem de átomos que a trajetória
    for alt, rotulo in ((rep_dir / "prod.gro", f"{rep}/prod.gro"),
                        (md_dir / "npt.gro", "npt.gro")):
        if alt.exists():
            print(f"      usando {rotulo} como topologia")
            return alt, rotulo

    # último recurso: o próprio gmx converte o .tpr que só ele entende
    gro = rep_dir / "topol_mda.gro"
    if tpr.exists() and not gro.exists():
        subprocess.run([gmx, "editconf", "-f", str(tpr), "-o", str(gro)],
                       capture_output=True, text=True)
    if gro.exists():
        print(f"      usando {gro.name} (convertido do .tpr pelo gmx)")
        return gro, gro.name

    raise SystemExit(
        f"não achei topologia legível para {rep}: nem um .tpr que o "
        f"MDAnalysis leia, nem um .gro do sistema solvatado.")


def carregar(md_dir: Path, rep: str):
    import MDAnalysis as mda
    xtc = md_dir / rep / "prod.xtc"
    if not xtc.exists():
        return None
    top, _ = topologia_legivel(md_dir, rep)
    u = mda.Universe(str(top), str(xtc))

    # Um .gro não traz elementos, e as seleções dependem de nomes. Conferir
    # aqui evita um RMSD calculado sobre zero átomo mais adiante.
    n_ca = len(u.select_atoms("protein and name CA"))
    if n_ca == 0:
        nomes = sorted({r.resname for r in u.residues})[:25]
        raise SystemExit(
            f"a seleção 'protein and name CA' não achou nada em {rep}.\n"
            f"Resíduos vistos: {nomes}\n"
            f"Se os nomes não são os padrão do campo de força, a topologia "
            f"({top.name}) não corresponde a esta trajetória.")
    return u


def descobrir_resname(u, pedido: str) -> str:
    """Nome do resíduo do ligante, conferido contra a trajetória.

    O nome escrito nas coordenadas nem sempre é o que o md_sistema.json diz:
    o acpype batiza a moleculetype como se pede em `-b`, mas propaga para o
    .gro o resíduo que veio do PDB de entrada (`UNL`, o default do RDKit).
    Em vez de parar por causa disso, acha o resíduo que sobra depois de tirar
    proteína, água e íons — que é o PROTAC, qualquer que seja seu nome.
    """
    if len(u.select_atoms(f"resname {pedido}")) > 0:
        return pedido
    resto = u.select_atoms(
        "not protein and not resname SOL WAT HOH TIP3 NA CL K MG CA ZN")
    nomes = sorted({r.resname for r in resto.residues})
    if len(nomes) == 1:
        print(f"      resíduo do ligante é '{nomes[0]}', não '{pedido}' "
              f"(nome herdado do PDB de entrada) — seguindo com ele")
        return nomes[0]
    if not nomes:
        raise SystemExit(
            f"não achei o ligante na trajetória: nem '{pedido}' nem nenhum "
            f"resíduo fora de proteína, água e íons")
    raise SystemExit(
        f"'{pedido}' não está na trajetória e há mais de um candidato a "
        f"ligante: {nomes}. Rode com --resname <nome>.")


def analisar_replica(u, resname: str):
    from MDAnalysis.analysis import rms
    from MDAnalysis.analysis.distances import distance_array

    prot = u.select_atoms("protein and name CA")
    lig = u.select_atoms(f"resname {resname}")
    if len(lig) == 0:
        raise SystemExit(f"resíduo {resname} não encontrado na trajetória")

    R = rms.RMSD(u, select="protein and name CA", ref_frame=0)
    R.run()
    rmsd_prot = R.results.rmsd[:, 2]

    Rl = rms.RMSD(u, select=f"resname {resname}",
                  groupselections=[f"resname {resname}"], ref_frame=0)
    Rl.run()
    rmsd_lig = Rl.results.rmsd[:, 3] if Rl.results.rmsd.shape[1] > 3 else \
        Rl.results.rmsd[:, 2]

    # o PROTAC é dividido pela geometria: a metade mais próxima da E3 no
    # primeiro frame é o lado do recrutador, a outra é o lado do warhead
    u.trajectory[0]
    d0 = distance_array(lig.positions, u.select_atoms("protein").positions)
    perto = d0.min(axis=1)
    mediana = np.median(perto)
    idx_recrutador = np.where(perto <= mediana)[0]
    idx_warhead = np.where(perto > mediana)[0]

    c_rec, c_wh = [], []
    for _ in u.trajectory:
        p = u.select_atoms("protein").positions
        d = distance_array(lig.positions, p)
        c_rec.append(int((d[idx_recrutador] <= 4.5).sum()))
        c_wh.append(int((d[idx_warhead] <= 4.5).sum()))

    return {
        "n_frames": len(rmsd_prot),
        "rmsd_proteina_media": float(np.mean(rmsd_prot)),
        "rmsd_proteina_final": float(rmsd_prot[-1]),
        "rmsd_protac_media": float(np.mean(rmsd_lig)),
        "rmsd_protac_final": float(rmsd_lig[-1]),
        "contatos_recrutador_media": float(np.mean(c_rec)),
        "contatos_warhead_media": float(np.mean(c_wh)),
        "fracao_frames_ancorado": float(np.mean(
            np.array(c_rec) >= CRITERIOS["contatos_recrutador_min"])),
        "_series": {"rmsd_prot": rmsd_prot.tolist(),
                    "rmsd_lig": np.asarray(rmsd_lig).tolist(),
                    "c_rec": c_rec, "c_wh": c_wh},
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--md-dir", type=Path, required=True)
    ap.add_argument("--resname", default=None)
    args = ap.parse_args()

    md = args.md_dir.expanduser()
    info = json.loads((md / "md_sistema.json").read_text())
    resname = args.resname or info.get("resname", "PTC")
    print(f"candidato: {info['candidate_id']} (rank {info.get('rank', '?')})\n")

    reps = sorted(d.name for d in md.glob("rep*") if d.is_dir())
    if not reps:
        raise SystemExit(f"nenhuma réplica em {md}")

    linhas, series = [], {}
    for rep in reps:
        u = carregar(md, rep)
        if u is None:
            print(f"  [{rep}] trajetória incompleta, pulada")
            continue
        resname = descobrir_resname(u, resname)
        r = analisar_replica(u, resname)
        series[rep] = r.pop("_series")
        r["replica"] = rep
        linhas.append(r)
        print(f"  [{rep}] {r['n_frames']} frames | "
              f"RMSD prot {r['rmsd_proteina_media']:.2f} Å | "
              f"RMSD PROTAC {r['rmsd_protac_media']:.2f} Å | "
              f"contatos recrut {r['contatos_recrutador_media']:.0f} / "
              f"warhead {r['contatos_warhead_media']:.0f}")

    if not linhas:
        raise SystemExit("nenhuma réplica analisável")

    df = pd.DataFrame(linhas)
    df.to_csv(md / "md_resultados.csv", index=False)
    (md / "md_series.json").write_text(json.dumps(series))

    m = df.mean(numeric_only=True)
    print(f"\n{'=' * 66}\nVEREDITO (média de {len(df)} réplicas)\n{'=' * 66}")

    checks = [
        ("RMSD da proteína", m["rmsd_proteina_media"],
         CRITERIOS["rmsd_proteina_max_A"], "<=", "Å"),
        ("RMSD do PROTAC", m["rmsd_protac_media"],
         CRITERIOS["rmsd_protac_max_A"], "<=", "Å"),
        ("contatos recrutador–E3", m["contatos_recrutador_media"],
         CRITERIOS["contatos_recrutador_min"], ">=", ""),
        ("contatos warhead–E3 (espúrio)", m["contatos_warhead_media"],
         CRITERIOS["contatos_warhead_max"], "<=", ""),
        ("frações ancoradas", m["fracao_frames_ancorado"],
         CRITERIOS["fracao_frames_ancorado_min"], ">=", ""),
    ]
    passou_tudo = True
    for nome, valor, corte, op, un in checks:
        ok = valor <= corte if op == "<=" else valor >= corte
        passou_tudo &= ok
        print(f"  {'OK ' if ok else 'NAO'}  {nome:32s} {valor:7.2f} {un} "
              f"({op} {corte})")

    print()
    if passou_tudo:
        print("  -> CANDIDATO APROVADO no nível (ii).")
        print("     Próximo passo: ternário completo (PRosettaC/AF3) e MD nível (iii).")
    else:
        print("  -> CANDIDATO REPROVADO no nível (ii).")
        if m["contatos_warhead_media"] > CRITERIOS["contatos_warhead_max"]:
            print("     O warhead está interagindo com a PRÓPRIA E3 ligase na")
            print("     ausência da PCSK9. É orientação improdutiva do linker:")
            print("     no ternário ele competiria com a ligação ao alvo.")
        print("\n     Rode o próximo da fila:")
        prox = json.loads((md.parent / "md_candidato.json").read_text())
        for i, c in enumerate(prox.get("proximos", [])[:3], start=2):
            print(f"       rank {i}: {c['candidate_id']}")
        print(f"\n       python scripts/md_prepare.py --rank 2 ...  "
              f"e depois md_run.sh")

    (md / "md_veredito.json").write_text(json.dumps({
        "candidate_id": info["candidate_id"], "aprovado": bool(passou_tudo),
        "criterios": CRITERIOS,
        "medias": {k: float(v) for k, v in m.items()},
    }, indent=2))
    print(f"\n{md / 'md_veredito.json'}")


if __name__ == "__main__":
    main()
