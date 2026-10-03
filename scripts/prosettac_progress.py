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


def serie_do_log_do_driver(pipeline_out: Path):
    """[(minutos, fila)] lidos dos heartbeats do driver, em ordem.

    O `+N min` é monotônico desde o início da espera, então não há ambiguidade
    de data nem de virada de meia-noite — o problema clássico de reconstruir
    tempo a partir de `[HH:MM]`.
    """
    logs = sorted(pipeline_out.glob("pipeline_*.log"),
                  key=lambda f: f.stat().st_mtime, reverse=True)
    for log in logs[:3]:
        pontos = []
        for m in re.finditer(r"\+(\d+) min \| fila: (\d+)",
                             log.read_text(errors="ignore")):
            pontos.append((int(m.group(1)), int(m.group(2))))
        if len(pontos) >= 2:
            return pontos, log
    return [], None


def taxa_da_serie(pontos):
    """(taxa_longa, taxa_recente) em jobs/h, ou (None, None).

    Duas janelas de propósito. A longa é a que projeta o fim; a recente mostra
    se a execução desacelerou — e é a diferença entre as duas que informa.
    Uma média de 5 dias esconde uma parada de 12 h; a janela curta sozinha
    confunde flutuação com tendência.
    """
    if len(pontos) < 2:
        return None, None
    t0, f0 = pontos[0]
    t1, f1 = pontos[-1]
    longa = (f0 - f1) / ((t1 - t0) / 60) if t1 > t0 and f0 > f1 else None
    # janela recente: os pontos dos últimos 180 min
    recentes = [p for p in pontos if p[0] >= t1 - 180]
    curta = None
    if len(recentes) >= 2:
        ta, fa = recentes[0]
        tb, fb = recentes[-1]
        if tb > ta and fa >= fb:
            curta = (fa - fb) / ((tb - ta) / 60)
    return longa, curta


def progresso_pela_fila(d: Path):
    """Progresso da etapa atual pela FILA, não pelos arquivos.

    Cada etapa do PRosettaC escreve arquivos diferentes, e um contador amarrado
    a um padrão de nome mede uma etapa e cega nas outras. O que vale para todas
    é o número de jobs: o log anuncia quantos foram emitidos, e o `squeue` diz
    quantos ainda não acabaram.

    A taxa precisa de duas medidas no tempo, e a melhor fonte delas NÃO é este
    script: é o log do driver. O esperar_results.sh já grava, a cada 30 min,

        [09:15] +2790 min | fila: 1785 | arquivos: 115487 | processos: 1

    ou seja, uma série temporal da fila com um contador de minutos monotônico —
    que não depende de data, de fuso, nem de quando o operador rodou isto pela
    última vez. Dezenas de pontos, de graça.

    O arquivo de estado ao lado do log continua como reserva, para quando o
    driver não estiver registrando. Mas ele tinha dois defeitos, e o primeiro
    explica uma pergunta sem resposta:

      silêncio   o `if dt > 300 and df > 0 / elif dt > 300` não tinha `else`.
                 Com a leitura anterior feita menos de 5 min antes, o script
                 não dizia NADA sobre a taxa — nem o número, nem por que não
                 havia número. Quem pergunta "está ok?" fica sem resposta e sem
                 saber que perguntou errado
      amnésia    o estado era sobrescrito em TODA execução. Rodar duas vezes
                 seguidas apagava a leitura antiga e destruía a base de
                 comparação justamente de quem estava conferindo com atenção
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

    # --- a taxa: primeiro do log do driver, que tem dezenas de pontos -----
    taxa_h = None
    medido = False          # o log do driver respondeu (mesmo que "parada")
    pontos, log_drv = serie_do_log_do_driver(d.parent.parent.parent)
    if pontos:
        longa, curta = taxa_da_serie(pontos)
        span_h = (pontos[-1][0] - pontos[0][0]) / 60
        print(f"  medido no log do driver ({log_drv.name}, {len(pontos)} "
              f"leituras em {humano(span_h * 3600)}):")
        if longa:
            print(f"      média: {longa:.0f} jobs/h")
            taxa_h = longa
            medido = True
        if curta is not None:
            nota = ""
            if longa and longa > 0 and curta > 0:
                if curta < 0.6 * longa:
                    nota = "  <- desacelerou"
                elif curta > 1.4 * longa:
                    nota = "  <- acelerou"
            print(f"      últimas 3 h: {curta:.0f} jobs/h{nota}")
            # A projeção usa a MENOR das duas. Projetar pela média quando a
            # execução desacelerou promete uma data que não vai acontecer, e
            # uma data otimista é pior que nenhuma: é com ela que se decide
            # desligar o computador e voltar na quinta.
            if taxa_h and curta > 0:
                taxa_h = min(taxa_h, curta)
            elif curta == 0:
                # Com a fila PARADA não existe projeção. A versão anterior caía
                # na média e anunciava "falta ~4d 18h" logo abaixo de "0 jobs/h
                # nas últimas 3 h" — duas linhas que se contradizem, e a data é
                # a que o operador leva embora.
                taxa_h = None
                medido = True      # medimos: a medida é "parada"
                print("      nenhum job concluído nas últimas 3 h: NÃO há data"
                      " a projetar.")
                print("      Se repetir na próxima leitura, é travamento —"
                      " confira a fila e o disco:")
                print("          squeue -u $USER | head")
                print("          df -h $(dirname $PWD)")

    # --- reserva: o arquivo de estado, para quando não há log do driver ---
    estado = d / ".progresso.json"
    agora = time.time()
    hist = []
    if estado.exists():
        try:
            dados = json.loads(estado.read_text())
            if dados.get("emitidos") == emitidos:
                hist = [tuple(x) for x in dados.get("leituras", [])]
        except Exception:
            hist = []
    # Guarda HISTÓRICO, não só a última: sobrescrever a cada execução apagava a
    # base de comparação de quem rodava duas vezes seguidas para conferir.
    hist.append((agora, restantes))
    hist = hist[-60:]
    try:
        estado.write_text(json.dumps({"emitidos": emitidos, "leituras": hist}))
    except Exception:
        pass

    # A reserva só entra se o log do driver NÃO respondeu. "Fila parada" é uma
    # resposta, não ausência de fonte: sem esta distinção o script dizia "não há
    # data a projetar" e logo abaixo "confira se o driver está vivo" — duas
    # linhas que se contradizem, sobre um driver que acabou de ser lido.
    if taxa_h is None and not medido:
        # a leitura mais antiga que dê uma janela de pelo menos 10 min
        base = next((h for h in hist if agora - h[0] >= 600), None)
        if base:
            dt = agora - base[0]
            df = base[1] - restantes
            if df > 0:
                taxa_h = df / dt * 3600
                print(f"  medido entre duas leituras deste script "
                      f"({humano(dt)}): {df} jobs -> {taxa_h:.0f} jobs/h")
            else:
                print(f"  em {humano(dt)} de leituras deste script nenhum job "
                      f"concluiu — se repetir, é travamento")
        else:
            # O else que faltava. Sem ele, uma leitura anterior recente fazia o
            # script não dizer NADA sobre a taxa — nem o número, nem por que
            # não havia número.
            quando = humano(agora - hist[0][0]) if len(hist) > 1 else "agora"
            print(f"  ainda sem taxa: a leitura mais antiga é de {quando}, e a"
                  f" janela\n  mínima é 10 min. O driver também registra a fila"
                  f" a cada 30 min —\n  se este número não aparecer, confira se"
                  f" o driver está vivo:\n      pgrep -af esperar_results.sh")

    if taxa_h and taxa_h > 0:
        falta_s = restantes / taxa_h * 3600
        print(f"  falta: ~{humano(falta_s)}  (fim por volta de "
              f"{time.strftime('%d/%m %H:%M', time.localtime(agora + falta_s))})")

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

    # A guarda de "poucas soluções" ficava AQUI, antes do desvio para a fila —
    # e bloqueava a única medida que funciona nesta etapa. As soluções de
    # docking local são a matéria-prima de UMA etapa; a fila vale para todas.
    # Com o Patchdock_Results limpo ou renomeado, o script dizia "poucas
    # soluções para medir taxa" de uma execução com 2682 jobs medíveis na fila.
    # Guarda de uma etapa não pode barrar a medida de outra.
    if passou:
        print()
        if arquivos:
            ultimo = time.strftime("%d/%m %H:%M",
                                   time.localtime(arquivos[-1].stat().st_mtime))
            parado_h = (time.time() - arquivos[-1].stat().st_mtime) / 3600
            print(f"  O refinamento local terminou (último arquivo: {ultimo}, "
                  f"há {parado_h:.0f} h).")
        else:
            print("  O refinamento local já passou, e não há arquivos "
                  "*_docking_*.pdb aqui.")
        print(f"  A etapa atual é '{etapa}', que não escreve nesses arquivos.")
        progresso_pela_fila(d)
        return

    if n < 2:
        print("\n  poucas soluções de docking para medir taxa — rode de novo"
              " em 30 min.")
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
