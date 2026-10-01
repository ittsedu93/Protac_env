#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Quais processos do PRosettaC estão vivos. Para SOURCE, não para executar:
#
#     source "$(dirname "$0")/proc_prosettac.sh"
#     prosettac_orquestrador_vivo && echo "rodando"
#     echo "trabalhadores: $(prosettac_n_trabalhadores)"
#
# Existe num arquivo só porque dois scripts precisam da mesma resposta — o
# esperar_results.sh, para decidir se a execução morreu, e o status.sh, para
# dizer quantos processos por job. Duas implementações da mesma pergunta
# divergem; neste projeto já divergiram três vezes, e na última as duas estavam
# erradas de maneiras DIFERENTES, ao mesmo tempo:
#
#     status.sh         "processos do PRosettaC: 19"   contava a si mesmo
#     esperar_results   "processos: 1"                 via 1 de 19
#
# Duas armadilhas, as duas achadas em execução real
# -------------------------------------------------
# 1) BARRA DUPLA. O PRosettaC monta os caminhos dos trabalhadores assim:
#
#        python /mnt/hd2tb/Documentos/PRosettaC//constraint_generation.py
#                                            ^^
#
#    Um padrão com `PRosettaC/` casa a primeira barra e então exige o nome do
#    script onde está a SEGUNDA barra — e falha. Eram 18 processos trabalhando
#    reportados como zero, com "processos: 1" no log: a cara exata de uma
#    execução quase morta. `PRosettaC/+` aceita uma ou mais barras.
#
# 2) O PADRÃO CASA QUEM O MENCIONA. `pgrep -f` olha a linha de comando inteira,
#    e o diretório do job é .../PRosettaC_runs/.../prosettac/SC0013__WH023 — a
#    string está no argumento de quem pergunta. Excluir $$ não basta (casou o
#    shell avô), e excluir parentesco é uma lista que nunca fecha. O que fecha é
#    olhar o CONTEÚDO: um processo do PRosettaC é um interpretador python cujo
#    primeiro argumento é um .py de dentro da instalação. Um shell que menciona
#    o caminho não é isso.
# ---------------------------------------------------------------------------

# `/+` = uma ou mais barras (ERE), por causa da armadilha 1.
PROSETTAC_PADRAO_ORQ='PRosettaC/+main\.py'
PROSETTAC_PADRAO_TRAB='PRosettaC/+[A-Za-z0-9_]+\.py'

# _prosettac_eh <pid> [so_main] -> 0 se for processo do PRosettaC
_prosettac_eh() {
  local pid="$1" so_main="${2:-0}" c0 c1
  local -a argv
  [[ "$pid" == "$$" ]] && return 1
  [[ -r "/proc/$pid/cmdline" ]] || return 1
  mapfile -d '' -t argv < "/proc/$pid/cmdline" 2>/dev/null || return 1
  c0="${argv[0]:-}"; c1="${argv[1]:-}"
  [[ "$c0" == *python* ]] || return 1
  # glob, não regex: `*` casa `/` também, então a barra dupla não incomoda aqui
  [[ "$c1" == *PRosettaC/*.py ]] || return 1
  if [[ "$so_main" == "1" ]]; then
    [[ "$c1" == */main.py ]] || return 1
  fi
  return 0
}

prosettac_orquestrador_vivo() {
  local pid u="${USER:-$(id -un)}"
  for pid in $(pgrep -u "$u" -f "$PROSETTAC_PADRAO_ORQ" 2>/dev/null); do
    _prosettac_eh "$pid" 1 && return 0
  done
  return 1
}

prosettac_n_trabalhadores() {
  local n=0 pid u="${USER:-$(id -un)}"
  for pid in $(pgrep -u "$u" -f "$PROSETTAC_PADRAO_TRAB" 2>/dev/null); do
    _prosettac_eh "$pid" && n=$((n + 1))
  done
  echo "$n"
}

prosettac_fila_tem_algo() {
  squeue -h -u "${USER:-$(id -un)}" 2>/dev/null | grep -q .
}
