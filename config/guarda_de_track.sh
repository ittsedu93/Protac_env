# ---------------------------------------------------------------------------
# Sourceável DEPOIS da config: recusa rodar no track abandonado.
#
#     source "$AQUI/../config/guarda_de_track.sh"
#
# Por que isto existe
# -------------------
# A CRBN foi abandonada como decisão de projeto: o PRosettaC não produziu
# geometria ternária com ela, e o Boltz-2 concordou. Mas o default de TODO
# script deste repositório era `config/pipeline.conf`, que é a config DA CRBN —
# e quatro deles ignoravam a variável PIPELINE_CONF, então nem passá-la
# ajudava.
#
# Isso não dá erro: dá números, do track errado, com a mesma cara de certo. Já
# aconteceu uma vez (o status.sh respondendo sobre a CRBN quando se pedia a
# VHL), e o custo de acontecer de novo cresce à medida que o track da VHL
# acumula resultados que irão para a tese.
#
# A guarda transforma "lembrar de passar PIPELINE_CONF" em "não é possível
# esquecer". Para rodar no track antigo de propósito — reproduzir um número da
# CRBN para a comparação da tese, por exemplo — PROTAC_TRACK_OK=1 libera, e
# fica registrado no comando que foi deliberado.
# ---------------------------------------------------------------------------
if [[ "${TRACK:-CRBN}" != "${PROTAC_TRACK_ESPERADO:-VHL}" ]]; then
  if [[ "${PROTAC_TRACK_OK:-0}" != "1" ]]; then
    echo "" >&2
    echo "*** TRACK ERRADO: esta config é do track '${TRACK:-CRBN}' e o" >&2
    echo "*** projeto trabalha no '${PROTAC_TRACK_ESPERADO:-VHL}'." >&2
    echo "***" >&2
    echo "*** config carregada: ${CONF:-?}" >&2
    echo "***" >&2
    echo "*** A CRBN foi abandonada: o PRosettaC não produziu geometria" >&2
    echo "*** ternária com ela e o Boltz-2 concordou. Rodar nela agora" >&2
    echo "*** escreveria em \$WORK/pipeline e misturaria os dois resultados." >&2
    echo "***" >&2
    echo "*** Use:" >&2
    echo "***   PIPELINE_CONF=\$HOME/Protac_env/config/pipeline_vhl.conf \\" >&2
    echo "***     <o comando que você ia rodar>" >&2
    echo "***" >&2
    echo "*** Para rodar na CRBN de propósito (reproduzir um número para a" >&2
    echo "*** comparação da tese): PROTAC_TRACK_OK=1 antes do comando." >&2
    echo "" >&2
    exit 3
  fi
  echo "  [TRACK ${TRACK:-CRBN}] liberado por PROTAC_TRACK_OK=1" >&2
fi
