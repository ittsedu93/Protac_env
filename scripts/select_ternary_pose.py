#!/usr/bin/env python
"""
Lê os modelos ternários do Boltz-2 e escolhe o que vai para a MD de nível (iii).

Só biblioteca padrão + pandas. Segundos.

    python select_ternary_pose.py --candidato SC0006__WH023

O que ele mede, e por quê
-------------------------
O `confidence_score` do Boltz é um composto da própria ferramenta. Ele é útil,
mas não é a pergunta: um modelo pode ter score alto com as duas proteínas
confiantes e a INTERFACE entre elas frouxa — e é a interface que decide se o
ternário existe.

Então saem três números por modelo, e o ranking usa os três:

  iptm_E3_alvo    a interface que o PROTAC precisa criar. Abaixo de ~0,5 o
                  modelo não sabe onde as duas proteínas se encontram.
  ligand_iptm     o PROTAC foi colocado com confiança, ou jogado no meio?
  complex_plddt   a estrutura como um todo

Os índices de cadeia do Boltz (0, 1, 2...) seguem a ordem do YAML; o
entradas.json guarda o mapa, então nada aqui depende de eu lembrar a ordem.
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path


def papel(descricao: str) -> str:
    """'E3 cadeia C (361 res)' -> 'E3'."""
    d = descricao.upper()
    if d.startswith("PROTAC"):
        return "PROTAC"
    return descricao.split()[0]


def mapa_de_cadeias(entradas: Path):
    """Índice do Boltz -> papel, pela ordem em que o YAML foi escrito."""
    d = json.loads(entradas.read_text())
    ids = d.get("boltz_ids", {})
    return {i: papel(desc) for i, desc in enumerate(ids.values())}, d


def par(m, a, b):
    """ipTM do par, simetrizado: a matriz do Boltz não é simétrica."""
    try:
        return (m[str(a)][str(b)] + m[str(b)][str(a)]) / 2
    except (KeyError, TypeError):
        return None


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidato", required=True)
    ap.add_argument("--work", type=Path,
                    default=Path.home() / "PRosettaC_runs" / "vhl_crbn_pcsk9_protac")
    ap.add_argument("--iptm-min", type=float, default=0.5,
                    help="corte da interface E3-alvo (default 0.5)")
    args = ap.parse_args()

    base = args.work.expanduser() / "boltz_ternario" / args.candidato
    entradas = base / "entradas.json"
    if not entradas.exists():
        raise SystemExit(f"não achei {entradas}")
    cadeias, info = mapa_de_cadeias(entradas)

    preds = sorted(base.glob("boltz_results_*/predictions/*/confidence_*.json"))
    if not preds:
        raise SystemExit(f"não achei confidence_*.json sob {base}")

    idx_e3 = [i for i, p in cadeias.items() if p == "E3"]
    idx_lig = [i for i, p in cadeias.items() if p == "PROTAC"]
    idx_alvo = [i for i, p in cadeias.items() if p not in ("E3", "PROTAC")]

    print(f"candidato: {args.candidato}")
    print("  cadeias: " + ", ".join(f"{i}={p}" for i, p in cadeias.items()))
    print()

    linhas = []
    for p in preds:
        d = json.loads(p.read_text())
        m = d.get("pair_chains_iptm", {})
        # a interface que precisa existir: o PIOR par E3 x alvo, não a média.
        # Uma interface boa com uma cadeia e ruim com a outra não é interface.
        pares = [par(m, a, b) for a in idx_e3 for b in idx_alvo]
        pares = [x for x in pares if x is not None]
        pares_lig = [par(m, a, b) for a in idx_lig for b in (idx_e3 + idx_alvo)]
        pares_lig = [x for x in pares_lig if x is not None]
        linhas.append({
            "modelo": p.stem.replace("confidence_", ""),
            "confidence": round(d.get("confidence_score", float("nan")), 3),
            "iptm_E3_alvo": round(min(pares), 3) if pares else None,
            "ligand_iptm": round(d.get("ligand_iptm", float("nan")), 3),
            "lig_pior_par": round(min(pares_lig), 3) if pares_lig else None,
            "complex_plddt": round(d.get("complex_plddt", float("nan")), 3),
            "pdb": str(p.with_name(p.name.replace("confidence_", "")
                                   .replace(".json", ".pdb"))),
        })

    # Tabela e CSV com a biblioteca padrão: este script precisa rodar em
    # qualquer env, e uma dependência a menos é uma falha a menos.
    for l in linhas:
        l["aprovado"] = (l["iptm_E3_alvo"] is not None
                         and l["iptm_E3_alvo"] >= args.iptm_min)
    linhas.sort(key=lambda l: (l["aprovado"], l["iptm_E3_alvo"] or -1,
                               l["confidence"] or -1), reverse=True)

    cols = ["modelo", "confidence", "iptm_E3_alvo", "ligand_iptm",
            "lig_pior_par", "complex_plddt", "aprovado"]
    larg = {c: max(len(c), *(len(str(l[c])) for l in linhas)) for c in cols}
    print("  " + "  ".join(c.rjust(larg[c]) for c in cols))
    for l in linhas:
        print("  " + "  ".join(str(l[c]).rjust(larg[c]) for c in cols))

    print(f"\n  iptm_E3_alvo = o PIOR par entre cadeias da E3 e do alvo.")
    print(f"  Corte: {args.iptm_min}. Abaixo disso o modelo não sabe onde as")
    print(f"  duas proteínas se encontram, e o resto do número não importa.")

    import csv as _csv
    with open(base / "ranking_poses.csv", "w", newline="") as fh:
        w = _csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for l in linhas:
            w.writerow({c: l[c] for c in cols})

    melhores = [l for l in linhas if l["aprovado"]]
    if not melhores:
        print("\n  -> NENHUM modelo passa o corte da interface.")
        print("     O Boltz não convergiu para uma geometria ternária. Espere")
        print("     o PRosettaC, que parte das duas poses já validadas.")
        return

    top = melhores[0]
    alvo = base / "melhor_ternario.pdb"
    shutil.copy2(top["pdb"], alvo)
    print(f"\n  -> ESCOLHIDO: {top['modelo']}")
    print(f"     interface E3-alvo {top['iptm_E3_alvo']} | "
          f"ligante {top['ligand_iptm']} | plDDT {top['complex_plddt']}")
    print(f"     {alvo}")

    (base / "melhor_ternario.json").write_text(json.dumps({
        "candidate_id": args.candidato,
        "fonte": "Boltz-2",
        "modelo": top["modelo"],
        "pdb": str(alvo),
        "confidence_score": float(top["confidence"]),
        "iptm_E3_alvo_pior_par": float(top["iptm_E3_alvo"]),
        "ligand_iptm": float(top["ligand_iptm"]),
        "complex_plddt": float(top["complex_plddt"]),
        "cadeias": {str(k): v for k, v in cadeias.items()},
        "criterio_iptm_min": args.iptm_min,
    }, indent=2))

    print(f"\n  ATENÇÃO: estas são as métricas de confiança do PRÓPRIO modelo,")
    print(f"  não validação experimental. Modelos de estrutura são conhecidos")
    print(f"  por superestimar a colocação de ligantes grandes e flexíveis em")
    print(f"  interfaces proteína-proteína — e este tem 19 torções. É por isso")
    print(f"  que o WP3 pede DOIS métodos: a concordância com o PRosettaC vale")
    print(f"  mais que qualquer um dos dois sozinho.")


if __name__ == "__main__":
    main()
