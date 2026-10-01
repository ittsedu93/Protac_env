#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Espera o Results/ do PRosettaC, desistindo por FALTA DE PROGRESSO.
#
#   bash scripts/esperar_results.sh <dir_do_job>
#
#   códigos: 0 = Results apareceu
#            1 = o PRosettaC terminou SEM Results (conclusão sobre o candidato)
#            2 = TRAVAMENTO (conclusão sobre a máquina, não sobre o candidato)
#
# Por que isto é um script e não um laço dentro do driver
# -------------------------------------------------------
# Porque ele já errou duas vezes, e das duas o erro só apareceria depois de
# horas ou dias:
#
#   1) o contador de progresso estava amarrado ao padrão de nome de UMA etapa,
#      e congelou em "92,5%, falta 1h54" de uma etapa já concluída
#   2) o limite era de RELÓGIO (168 h), e a etapa de conformações do linker foi
#      medida em 12 jobs/h com 2532 na fila — 8,6 dias. O driver declararia
#      falha em 08/10 de uma execução que termina em 10/10, jogando fora 7 dias
#      de máquina com uma mensagem sobre "ausência de geometria ternária" que
#      seria simplesmente falsa
#
# Um laço que só pode ser exercitado em 12 h não é um laço testado. Aqui ele
# roda sozinho, com TERNARIO_ESPERA_S=1 no teste e 300 em produção, e os três
# desfechos se verificam em segundos.
#
# O que conta como progresso
# --------------------------
# Três medidas, e qualquer uma delas mexendo basta:
#
#   a fila do SLURM        o sinal desta etapa; os contadores de arquivo ficam
#                          parados nas conformações do linker porque ela não
#                          escreve em Patchdock_Results
#   arquivos no diretório  o sinal das etapas de docking
#   tamanho do log.txt     cobre as transições entre etapas, quando a fila
#                          esvazia e nenhum arquivo nasce por alguns minutos
#
# Enquanto um dos três se move, a execução está trabalhando, e esperar é o
# certo por mais dias que leve. Relógio não distingue devagar de parado — e
# essa é exatamente a pergunta.
# ---------------------------------------------------------------------------
set -uo pipefail

DIRC="${1:?uso: esperar_results.sh <dir_do_job>}"
[[ -d "$DIRC" ]] || { echo "não é um diretório: $DIRC"; exit 3; }

ESPERA_S="${TERNARIO_ESPERA_S:-300}"          # 5 min em produção
SEM_PROG_H="${TERNARIO_SEM_PROGRESSO_H:-12}"  # sem mudança por tanto tempo
LIMITE_H="${TERNARIO_LIMITE_H:-336}"          # teto absoluto, só como rede
BATIDAS="${TERNARIO_BATIDAS:-6}"              # 1 relatório a cada 6 voltas
# Janela da quarta medida (abaixo): o dobro do intervalo entre relatórios, para
# que uma escrita feita entre dois relatórios não passe sem ser vista.
JANELA_MIN=$(( 2 * BATIDAS * ESPERA_S / 60 )); (( JANELA_MIN < 1 )) && JANELA_MIN=1
U="${USER:-$(id -un)}"

# O PADRÃO IMPORTA, e descobri isso no teste: `pgrep -f "PRosettaC"` casa o
# PRÓPRIO processo deste script, porque o diretório do job que ele recebe como
# argumento é .../PRosettaC_runs/.../prosettac/SC0013__WH023 — a string está
# na linha de comando. A consequência não é teórica: a execução nunca pareceria
# ter terminado, o script esperaria as 12 h do travamento, e o diagnóstico
# diria TRAVAMENTO onde a verdade é "terminou sem Results" — conclusões
# diferentes, uma sobre a máquina e outra sobre o candidato.
#
# `PRosettaC/main\.py` é o orquestrador, vive enquanto a execução viver, e a
# barra depois do nome não casa com `PRosettaC_runs`. O $$ e o $PPID ficam
# excluídos de qualquer forma, porque isto é barato e definitivo.
PADRAO_ORQ='PRosettaC/main\.py'
PADRAO_TRAB='PRosettaC/(main|constraint_generation|distance_constraints)\.py'

# Nem o padrão nem excluir $$/$PPID bastam: no teste quem casou foi o shell
# AVÔ, cuja linha de comando continha a string. Excluir parentesco é uma lista
# que nunca fecha. O que fecha é olhar o CONTEÚDO: um orquestrador de verdade é
# um interpretador python cujo primeiro argumento é o main.py de dentro da
# instalação do PRosettaC. Um shell que menciona o caminho não é isso.
eh_orquestrador() {  # eh_orquestrador <pid>
  local pid="$1" c0 c1
  [[ -r "/proc/$pid/cmdline" ]] || return 1
  mapfile -d '' -t argv < "/proc/$pid/cmdline" 2>/dev/null || return 1
  c0="${argv[0]:-}"; c1="${argv[1]:-}"
  [[ "$c0" == *python* ]] || return 1
  [[ "$c1" == */PRosettaC/main.py || "$c1" == */PRosettaC//main.py ]] || return 1
  return 0
}

orquestrador_vivo() {
  local pid
  for pid in $(pgrep -u "$U" -f "$PADRAO_ORQ" 2>/dev/null); do
    [[ "$pid" == "$$" ]] && continue
    eh_orquestrador "$pid" && return 0
  done
  return 1
}

# Os trabalhadores são python rodando um .py de dentro da instalação. Mesma
# verificação de conteúdo: argv[0] é python, argv[1] é um .py sob PRosettaC/.
n_trabalhadores() {
  local n=0 pid c0 c1
  for pid in $(pgrep -u "$U" -f "$PADRAO_TRAB" 2>/dev/null); do
    [[ "$pid" == "$$" ]] && continue
    [[ -r "/proc/$pid/cmdline" ]] || continue
    mapfile -d '' -t argv < "/proc/$pid/cmdline" 2>/dev/null || continue
    c0="${argv[0]:-}"; c1="${argv[1]:-}"
    [[ "$c0" == *python* ]] || continue
    [[ "$c1" == *PRosettaC/*.py || "$c1" == *PRosettaC//*.py ]] || continue
    n=$((n + 1))
  done
  echo "$n"
}

fila_tem_algo() { squeue -h -u "$U" 2>/dev/null | grep -q .; }

echo "  esperando o Results em $DIRC"
echo "  desiste por FALTA DE PROGRESSO ($SEM_PROG_H h sem mudança), não por"
echo "  tempo total; teto absoluto de $LIMITE_H h apenas como rede"

voltas=0; ref_prog=""; ref_voltas=0
total=$(( LIMITE_H * 3600 / ESPERA_S )); [[ $total -lt 1 ]] && total=1

for _ in $(seq 1 "$total"); do
  sleep "$ESPERA_S"
  voltas=$((voltas + 1))
  [[ -d "$DIRC/Results" ]] && { echo "  [$(date +%H:%M)] Results apareceu"; exit 0; }

  if (( voltas % BATIDAS == 0 )); then
    # `ls dir/*.pdb` com dezenas de milhares de arquivos estoura o limite de
    # argumentos e devolve zero — indistinguível de "não produziu".
    n_pd=$(find "$DIRC/Patchdock_Results" -maxdepth 1 -type f 2>/dev/null | wc -l)
    n_fila=$(squeue -h -u "$U" 2>/dev/null | wc -l)
    tam_log=$(stat -c %s "$DIRC/log.txt" 2>/dev/null || echo 0)
    n_proc=$(n_trabalhadores)
    etapa=$(tail -1 "$DIRC/log.txt" 2>/dev/null | cut -c1-58)
    echo "  [$(date +%H:%M)] +$((voltas * ESPERA_S / 60)) min | fila: $n_fila | arquivos: $n_pd | processos: $n_proc"
    [[ -n "$etapa" ]] && echo "            etapa: $etapa"

    # QUARTA MEDIDA, e ela é um OVERRIDE e não parte da tupla. Uma etapa que
    # reescreve os MESMOS arquivos deixa a contagem parada e o log do mesmo
    # tamanho — trabalho real que as três primeiras medidas não veem. Se algo
    # foi escrito aqui na janela, isso é progresso, e ponto.
    # (Não pode entrar na tupla: "escreveu recentemente" fica constante em
    #  "sim" enquanto a etapa trabalha, a tupla pararia de mudar, e o
    #  travamento seria declarado justamente durante o trabalho.)
    houve_escrita=0
    [[ -n "$(find "$DIRC" -maxdepth 2 -type f -newermt "-${JANELA_MIN} minutes" \
                ! -name 'prosetta_config.txt*' ! -name '.*' 2>/dev/null | head -1)" ]] \
      && houve_escrita=1

    prog="${n_fila}:${n_pd}:${tam_log}"
    if [[ "$prog" != "$ref_prog" || $houve_escrita -eq 1 ]]; then
      ref_prog="$prog"; ref_voltas=$voltas
    else
      paradas_min=$(( (voltas - ref_voltas) * ESPERA_S / 60 ))
      if (( paradas_min >= SEM_PROG_H * 60 )); then
        echo "  [$(date +%H:%M)] nada mudou em $((paradas_min / 60)) h: fila, arquivos e log.txt parados"
        exit 2
      fi
      (( paradas_min >= 60 )) && \
        echo "            (sem mudança há $((paradas_min / 60)) h de $SEM_PROG_H h)"
    fi
  fi

  # Morreu no meio: só é fim quando o PRosettaC diz que terminou. Processo
  # ausente E fila vazia é fim — de execução bem ou malsucedida, e quem
  # distingue as duas é a existência do Results, conferida no topo do laço.
  if ! orquestrador_vivo && ! fila_tem_algo; then
    echo "  [$(date +%H:%M)] nenhum processo do PRosettaC vivo e fila vazia"
    [[ -d "$DIRC/Results" ]] && exit 0
    exit 1
  fi
done

echo "  [$(date +%H:%M)] teto absoluto de $LIMITE_H h alcançado"
[[ -d "$DIRC/Results" ]] && exit 0
exit 2
