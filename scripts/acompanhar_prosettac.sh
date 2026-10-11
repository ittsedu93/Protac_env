#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Um contador de progresso que ESCREVE NUM LOG, para acompanhar dias de
# PRosettaC com um `tail -f`.
#
#   bash scripts/acompanhar_prosettac.sh SC0018__WH023 [minutos]
#
# Para rodar destacado, junto da corrida:
#
#   nohup bash scripts/acompanhar_prosettac.sh SC0018__WH023 15 \
#       > /dev/null 2>&1 &
#
# O relatório completo (prosettac_progress.py sem --uma-linha) responde "está
# ok?" UMA vez. Uma série de linhas responde "está andando?", que é outra
# pergunta — e é a que importa quando a resposta leva dias. Cada linha traz a
# etapa, as transformadas, o contador, a taxa MEDIDA nesta execução, o que
# falta, o consumo de disco e o livre da máquina.
#
# Para parar: `touch <log>.pare` ou mate o processo. Ele também para sozinho
# quando o Results/ aparece — ou seja, quando o PRosettaC termina.
# ---------------------------------------------------------------------------
set -uo pipefail

AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONF="${PIPELINE_CONF:-$AQUI/../config/pipeline.conf}"
[[ -f "$CONF" ]] || [[ "$CONF" = /* ]] || CONF="$AQUI/../$CONF"
[[ -f "$CONF" ]] || { echo "config não encontrada: $CONF"; exit 1; }
# shellcheck disable=SC1090
source "$CONF"
export PIPELINE_CONF="$CONF"
# shellcheck disable=SC1090
source "$AQUI/achar_job.sh"

CAND="${1:?uso: acompanhar_prosettac.sh <candidate_id> [minutos]}"
MIN="${2:-15}"
# Piso de 1 min: com 0 (ou um valor não numérico, que o bash lê como 0) este
# laço gira sem dormir e passa dias queimando um núcleo da máquina do
# laboratório — justamente o que ele existe para ajudar a não fazer.
[[ "$MIN" =~ ^[0-9]+$ ]] || { echo "minutos deve ser um número: '$MIN'"; exit 1; }
(( MIN < 1 )) && { echo "  [aviso] minutos=$MIN é curto demais; usando 1"; MIN=1; }
INTERVALO=$(( MIN * 60 ))

DIRC="$(achar_job_dir "$PIPELINE_OUT" "$CAND")" || {
  echo "não achei o diretório de job de $CAND em $PIPELINE_OUT"; exit 1; }
WORKDIR="$(dirname "$DIRC")"
LOGP="$PIPELINE_OUT/progresso_${CAND}.log"

{
  echo "=============================================================="
  echo "acompanhamento de $CAND — iniciado $(date -Is)"
  echo "  job em: $DIRC"
  echo "  uma linha a cada $MIN min; para parar: touch $LOGP.pare"
  echo "=============================================================="
} >> "$LOGP"

echo "escrevendo em: $LOGP"
echo "acompanhe com: tail -f $LOGP"

while true; do
  if [[ -e "$LOGP.pare" ]]; then
    echo "$(date '+%d/%m %H:%M') | PARADO a pedido ($LOGP.pare)" >> "$LOGP"
    rm -f "$LOGP.pare"
    break
  fi

  # O `du` é o instrumento que atribui o consumo a NÓS; o `df` diz o que sobra
  # para o resto do laboratório. Os dois, porque respondem perguntas
  # diferentes, e a segunda é a que decide se estamos incomodando alguém.
  GB_NOSSO=$(du -sBM "$DIRC" 2>/dev/null | cut -f1 | tr -dc '0-9')
  GB_LIVRE=$(df -BM --output=avail "$PIPELINE_OUT" 2>/dev/null | tail -1 | tr -dc '0-9')
  # python3 do sistema, sem conda: o prosettac_progress.py é só biblioteca
  # padrão (está escrito no topo dele). Um `conda run` a cada 15 min, por
  # dias, custaria segundos de CPU de graça e quebraria aqui se o env mudasse
  # de nome — e este laço precisa sobreviver a dias sem vigilância.
  LINHA="$(python3 "$AQUI/prosettac_progress.py" \
             --candidato "$CAND" --work "$WORKDIR" --uma-linha 2>/dev/null \
           | grep -v '^$' | tail -1)"
  [[ -n "$LINHA" ]] || LINHA="$(date '+%d/%m %H:%M') | (sem leitura do progresso)"
  printf '%s | nosso:%sGB livre:%sGB\n' "$LINHA" \
     "$(awk -v m="${GB_NOSSO:-0}" 'BEGIN{printf "%.1f", m/1024}')" \
     "$(awk -v m="${GB_LIVRE:-0}" 'BEGIN{printf "%.0f", m/1024}')" >> "$LOGP"

  if [[ -d "$DIRC/Results" ]]; then
    echo "$(date '+%d/%m %H:%M') | Results/ apareceu — o PRosettaC TERMINOU." >> "$LOGP"
    echo "  veredicto:  PIPELINE_CONF=$CONF bash $AQUI/status.sh" >> "$LOGP"
    break
  fi
  sleep "$INTERVALO"
done
