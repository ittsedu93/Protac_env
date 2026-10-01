#!/usr/bin/env python
"""
Quanto falta, medido nos carimbos de hora dos arquivos — não estimado.

Só biblioteca padrão. Instantâneo.

    python prosettac_progress.py --candidato SC0013__WH023

Por que medir em vez de estimar
-------------------------------
A duração do PRosettaC escala com as transformadas que o PatchDock devolveu, e
elas variam por duas ordens de grandeza entre pares de proteínas: a CRBN deu 34
e a VHL 2468. Uma estimativa vinda da CRBN errou por um fator de 30 no track da
VHL, e errar assim tem custo real — foi o que fez o driver quase declarar falha
numa execução que estava indo bem.

O que este script faz é diferente: conta as soluções de docking local já
escritas, olha o intervalo entre a primeira e a última, e divide. A taxa é
medida NESTA execução, nesta máquina, com esta carga.

O alvo tem uma suposição, e ela está declarada: o PRosettaC refina no máximo as
1000 melhores transformadas, e a razão medida na CRBN foi ~35 soluções de
docking por transformada (1184 soluções de 34 transformadas). O alvo é
`min(transformadas, 1000) x 35`, e `--alvo` o sobrescreve.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import time
from pathlib import Path

SOLUCOES_POR_TRANSFORMADA = 34.8   # 1184 / 34, medido no track da CRBN
TETO_DO_PROSETTAC = 1000           # `Global` com Full: True, no main.py dele


def transformadas(log: Path) -> int | None:
    if not log.exists():
        return None
    achados = re.findall(r"Transforms remain after interface clustering[^:]*:\s*(\d+)",
                         log.read_text(errors="ignore"))
    return int(achados[-1]) if achados else None


ETAPAS = [          # na ordem do main.py do PRosettaC
    ("Cleaning structures", "preparando receptores (clean + relax)"),
    ("Sampling the distance", "amostrando a distância entre as âncoras"),
    ("Running PatchDock", "docking global com a restrição"),
    ("Run Rosetta local docking", "refinamento local do Rosetta"),
    ("Generating up to 100 constrained", "conformações do linker"),
    ("Clustering the top results", "agrupamento final"),
    ("run has finished", "TERMINADO"),
]


def etapa_atual(log: Path):
    """A última etapa que o PRosettaC anunciou, e se o docking já acabou.

    Sem isto o script mede só os arquivos de docking local e, quando essa etapa
    acaba, reporta "92,5%, falta 1h54" para sempre — porque o trabalho mudou de
    lugar e ele continua olhando o lugar antigo. Percentual de uma etapa que
    terminou é pior que nenhum percentual.
    """
    if not log.exists():
        return None, False
    txt = log.read_text(errors="ignore")
    atual, passou_docking = None, False
    for chave, nome in ETAPAS:
        if chave in txt:
            atual = nome
            if chave in ("Generating up to 100 constrained",
                         "Clustering the top results", "run has finished"):
                passou_docking = True
    return atual, passou_docking


def progresso_pela_fila(d: Path):
    """Progresso da etapa atual pela FILA, não pelos arquivos.

    Cada etapa do PRosettaC escreve arquivos diferentes, e um contador amarrado
    a um padrão de nome mede uma etapa e cega nas outras. O que vale para todas
    é o número de jobs: o log anuncia quantos foram emitidos, e o `squeue` diz
    quantos ainda não acabaram.

    A taxa precisa de duas medidas. Em vez de pedir ao operador que rode duas
    vezes e faça a conta, guardamos a leitura anterior ao lado do log.
    """
    emitidos = 0
    txt = (d / "log.txt").read_text(errors="ignore") if (d / "log.txt").exists() else ""
    listas = re.findall(r"jobs: \[(.*?)\]", txt)
    if listas:
        emitidos = len([x for x in listas[-1].split(",") if x.strip()])
    try:
        saida = subprocess.run(["squeue", "-h", "-u", os.environ.get("USER", "")],
                               capture_output=True, text=True, timeout=20)
        restantes = len([l for l in saida.stdout.splitlines() if l.strip()])
        tem_fila = saida.returncode == 0
    except Exception:
        restantes, tem_fila = 0, False

    if not tem_fila:
        print("\n  (sem `squeue` nesta máquina — acompanhe por "
              "`pgrep -u $USER -c -f PRosettaC`)")
        return
    if not emitidos:
        print(f"\n  jobs ainda na fila: {restantes}")
        return

    feitos = max(emitidos - restantes, 0)
    print()
    print(f"  jobs desta etapa: {emitidos} emitidos, {restantes} na fila, "
          f"{feitos} concluídos")
    print(f"  progresso da etapa: {100 * feitos / emitidos:.1f}%")

    # taxa a partir da leitura anterior
    estado = d / ".progresso.json"
    agora = time.time()
    anterior = None
    if estado.exists():
        try:
            anterior = json.loads(estado.read_text())
        except Exception:
            anterior = None
    try:
        estado.write_text(json.dumps({"t": agora, "restantes": restantes,
                                      "emitidos": emitidos}))
    except Exception:
        pass

    if anterior and anterior.get("emitidos") == emitidos:
        dt = agora - anterior["t"]
        df = anterior["restantes"] - restantes
        if dt > 300 and df > 0:
            taxa = df / dt
            print(f"  desde a leitura anterior ({humano(dt)}): {df} jobs "
                  f"-> {taxa * 3600:.0f} jobs/h")
            print(f"  falta: ~{humano(restantes / taxa)}  (fim por volta de "
                  f"{time.strftime('%d/%m %H:%M', time.localtime(agora + restantes / taxa))})")
        elif dt > 300:
            print(f"  desde a leitura anterior ({humano(dt)}): nenhum job "
                  f"concluiu — se isso se repetir, há travamento")
    else:
        print("  rode de novo em ~30 min: com duas leituras sai a taxa e o fim")

    print()
    print("  Enquanto a fila baixa, é só esperar. O fim da etapa é o "
          "diretório\n  Results/ aparecer.")


def humano(seg: float) -> str:
    if seg < 0:
        return "?"
    h, m = divmod(int(seg // 60), 60)
    d, h = divmod(h, 24)
    if d:
        return f"{d}d {h}h {m}min"
    return f"{h}h {m}min" if h else f"{m}min"


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidato", required=True)
    ap.add_argument("--work", type=Path,
                    default=Path.home() / "PRosettaC_runs"
                    / "vhl_crbn_pcsk9_protac" / "pipeline_vhl" / "wp3"
                    / "prosettac")
    ap.add_argument("--alvo", type=int,
                    help="quantas soluções esperar (sobrescreve a suposição)")
    args = ap.parse_args()

    d = args.work.expanduser() / args.candidato
    if not d.is_dir():
        raise SystemExit(f"não achei {d}")

    print(f"candidato: {args.candidato}")
    if (d / "Results").is_dir():
        print("  Results/ já existe — TERMINOU")
        resumo = d / "result_summary.txt"
        if resumo.exists():
            print()
            print("\n".join("  " + l for l in resumo.read_text().splitlines()))
        return

    pd_dir = d / "Patchdock_Results"
    arquivos = sorted(pd_dir.glob("*_docking_????.pdb"),
                      key=lambda p: p.stat().st_mtime) if pd_dir.is_dir() else []
    n = len(arquivos)
    n_pd = len(list(pd_dir.iterdir())) if pd_dir.is_dir() else 0
    tr = transformadas(d / "run_prosettac.log")

    etapa, passou = etapa_atual(d / "log.txt")
    if etapa:
        print(f"  etapa atual: {etapa}")
    print(f"  Patchdock_Results: {n_pd} arquivos")
    print(f"  transformadas do PatchDock: {tr if tr is not None else '?'}"
          + (f"  (o PRosettaC refina as {TETO_DO_PROSETTAC} melhores)"
             if tr and tr > TETO_DO_PROSETTAC else ""))
    print(f"  soluções de docking local escritas: {n}")

    if n < 2:
        print("\n  poucas soluções para medir taxa — rode de novo em 30 min.")
        return

    if passou:
        ultimo = time.strftime("%d/%m %H:%M",
                               time.localtime(arquivos[-1].stat().st_mtime))
        parado_h = (time.time() - arquivos[-1].stat().st_mtime) / 3600
        print()
        print(f"  O refinamento local terminou (último arquivo: {ultimo}, "
              f"há {parado_h:.0f} h).")
        print(f"  A etapa atual é '{etapa}', que não escreve nesses arquivos.")
        progresso_pela_fila(d)
        return

    t0, t1 = arquivos[0].stat().st_mtime, arquivos[-1].stat().st_mtime
    decorrido = t1 - t0
    if decorrido <= 0:
        print("\n  intervalo nulo entre a primeira e a última — sem taxa ainda.")
        return
    taxa = n / decorrido                      # soluções por segundo
    print(f"  primeira às {time.strftime('%d/%m %H:%M', time.localtime(t0))}, "
          f"última às {time.strftime('%d/%m %H:%M', time.localtime(t1))}")
    print(f"  decorrido: {humano(decorrido)} | taxa: {taxa * 60:.1f} soluções/min")

    alvo = args.alvo
    suposto = False
    if alvo is None and tr:
        alvo = int(min(tr, TETO_DO_PROSETTAC) * SOLUCOES_POR_TRANSFORMADA)
        suposto = True
    if not alvo:
        print("\n  sem alvo: passe --alvo para estimar o fim.")
        return

    falta = max(alvo - n, 0)
    print()
    print(f"  alvo: ~{alvo} soluções"
          + ("  (suposição: min(transformadas, 1000) x 34,8 — a razão medida"
             " na CRBN)" if suposto else ""))
    print(f"  progresso: {100 * n / alvo:.1f}%")
    print(f"  falta: ~{falta} soluções -> ~{humano(falta / taxa)}")
    print(f"  fim estimado: "
          f"{time.strftime('%d/%m às %H:%M', time.localtime(time.time() + falta / taxa))}")
    print()
    print("  A taxa é a MEDIDA desta execução; o alvo é que é suposto. Se a"
          " razão\n  soluções/transformada desta E3 for diferente da CRBN, o fim"
          " muda — o\n  progresso em % é o número a acompanhar, não a data.")


if __name__ == "__main__":
    main()
