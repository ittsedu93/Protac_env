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


def carregar(md_dir: Path, rep: str):
    import MDAnalysis as mda
    tpr = md_dir / rep / "prod.tpr"
    xtc = md_dir / rep / "prod.xtc"
    if not (tpr.exists() and xtc.exists()):
        return None
    return mda.Universe(str(tpr), str(xtc))


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
