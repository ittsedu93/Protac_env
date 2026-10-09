#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Põe de lado tudo o que foi construído com o recrutador ANTIGO.
#
#   bash scripts/arquivar_recrutador_antigo.sh             # mostra o plano
#   bash scripts/arquivar_recrutador_antigo.sh --aplicar   # move
#
# Por que isto é necessário, e não higiene
# ----------------------------------------
# O recrutador da VHL foi trocado: o antigo era um ligante de CRBN (glutarimida
# + fluorftalimida — talidomida), o novo é o VH101/VH298 cristalográfico do
# 6GFZ. Tudo o que o pipeline produziu entre as fases 5 e 9 foi construído
# sobre o antigo.
#
# E as fases NÃO refazem sozinhas. Medido na fonte:
#
#   prosettac_prepare.sh:60    if [[ -s Init0.pdb && -s Init1.pdb && ... ]]
#                              echo "  preparação já existe — nada a fazer"
#
# A fase 6a chama esse script antes de medir o vão. Se o diretório do candidato
# ainda tiver os Init0/Init1 do recrutador antigo, ela os REUSA e devolve a
# mesma medição de antes — e os nomes dos candidatos (SC0013__WH023) são
# derivados de numeração, então um candidato novo pode cair num diretório
# velho. O número que sairia (vão útil 10,0 Å -> alcance exigido 15,6 Å) é
# exatamente o que o preflight ainda mostra: ele é do recrutador errado.
#
# Mover, não apagar: o resultado do PRosettaC antigo (Results/, score.sc) é o
# registro de uma semana de máquina e entra na tese como tentativa descartada,
# com o motivo. E `mv` dentro do mesmo sistema de arquivos é renomeação —
# ZERO bytes escritos, o que importa num disco compartilhado.
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
source "$(dirname "$CONF")/guarda_de_track.sh"

APLICAR=0
[[ "${1:-}" == "--aplicar" ]] && APLICAR=1
[[ -n "${1:-}" && "$APLICAR" == "0" ]] && { echo "opção desconhecida: $1"; exit 1; }

echo "=============================================================="
echo "Arquivar o que foi construído com o recrutador antigo"
echo "=============================================================="
echo "  config:   $CONF   (track ${TRACK:-?})"
echo "  saída:    $PIPELINE_OUT"

[[ -d "$PIPELINE_OUT" ]] || { echo "  $PIPELINE_OUT não existe — nada a fazer"; exit 0; }
if [[ -z "${E3_RECRUITER_SDF:-}" ]]; then
  echo "*** E3_RECRUITER_SDF está vazio na config. Sem saber qual é o"
  echo "*** recrutador ATUAL eu não sei o que é velho — e não mexo no disco"
  echo "*** às cegas. Rode o definir_recrutador.py antes."
  exit 2
fi
ATUAL="$(basename "$E3_RECRUITER_SDF")"
echo "  recrutador atual: $ATUAL"

# --- nada se move com o PRosettaC vivo ------------------------------------
# shellcheck disable=SC1090
source "$AQUI/proc_prosettac.sh"
N_VIVOS="$(prosettac_n_trabalhadores 2>/dev/null || echo 0)"
if [[ "${N_VIVOS:-0}" -gt 0 ]] || prosettac_orquestrador_vivo 2>/dev/null; then
  echo ""
  echo "*** Há PRosettaC rodando ($N_VIVOS trabalhadores). Mover o diretório"
  echo "*** debaixo de um processo que escreve nele é como se perde uma corrida."
  echo "*** Pare primeiro, ou espere."
  exit 3
fi

# --- quem referencia qual recrutador --------------------------------------
echo ""
echo "  diretórios de job e o recrutador de cada um:"
VELHOS=0; NOVOS=0; SEM=0
while IFS= read -r cfg; do
  [[ -n "$cfg" ]] || continue
  h="$(awk '/^Heads:/{print $2; exit}' "$cfg" 2>/dev/null)"
  b="$(basename "${h:-?}")"
  if [[ "$b" == "$ATUAL" ]]; then NOVOS=$((NOVOS+1))
  elif [[ -z "$h" ]]; then SEM=$((SEM+1))
  else
    VELHOS=$((VELHOS+1))
    [[ "$VELHOS" -le 3 ]] && echo "      [velho] $(basename "$(dirname "$cfg")")  ->  $b"
  fi
  # o que já foi arquivado não conta como "velho em uso": senão uma segunda
  # execução acharia o recrutador antigo dentro do próprio arquivo morto.
done < <(find "$PIPELINE_OUT" -maxdepth 5 -type f -name prosetta_config.txt \
              -not -path "*/_recrutador_antigo_*" 2>/dev/null | sort)
[[ "$VELHOS" -gt 3 ]] && echo "      ... e outros $((VELHOS-3)) com o mesmo recrutador velho"
echo "      total: $VELHOS velhos, $NOVOS já no atual, $SEM sem linha 'Heads:'"

# --- o que sai da frente --------------------------------------------------
ALVOS=(
  "$PIPELINE_OUT/wp2"
  "$PIPELINE_OUT/wp3"
  "$PIPELINE_OUT/span_requirement.json"
  "$PIPELINE_OUT/protac_ranking.csv"
  "$PIPELINE_OUT/ortogonal_confirmado.json"
)
# Os marcadores das fases construídas sobre o recrutador antigo. Sem apagá-los
# o driver diz "[fase 5 pulada]" e reusa tudo. A 1-4 ficam: warheads, docking e
# a escolha do recrutador não dependem do SDF que foi trocado.
for f in 5 6 6a 6b 7 8 9 10; do ALVOS+=("$PIPELINE_OUT/.done_$f"); done

DESTINO="$PIPELINE_OUT/_recrutador_antigo_$(date +%Y%m%d_%H%M%S)"
echo ""
echo "  destino: $DESTINO"
echo ""
PRESENTES=()
for a in "${ALVOS[@]}"; do
  [[ -e "$a" ]] || continue
  PRESENTES+=("$a")
  if [[ -d "$a" ]]; then
    printf "      %-34s %8s  (%s arquivos)\n" "$(basename "$a")/" \
      "$(du -sh "$a" 2>/dev/null | cut -f1)" \
      "$(find "$a" -type f 2>/dev/null | wc -l)"
  else
    printf "      %-34s %8s\n" "$(basename "$a")" \
      "$(du -sh "$a" 2>/dev/null | cut -f1)"
  fi
done
if [[ ${#PRESENTES[@]} -eq 0 ]]; then
  echo "      (nada presente — já está limpo)"
  exit 0
fi

if [[ "$VELHOS" -eq 0 && "$NOVOS" -gt 0 ]]; then
  echo ""
  echo "  Todos os jobs já referenciam o recrutador ATUAL. Não há o que"
  echo "  arquivar — e eu não mexo numa saída que já é a boa."
  exit 0
fi

# `mv` só é gratuito dentro do mesmo sistema de arquivos; entre dois, é cópia.
DEV_ORI="$(stat -c %d "$PIPELINE_OUT" 2>/dev/null)"
mkdir -p "$DESTINO" 2>/dev/null || { echo "*** não consegui criar $DESTINO"; exit 1; }
DEV_DST="$(stat -c %d "$DESTINO" 2>/dev/null)"
if [[ "$DEV_ORI" != "$DEV_DST" ]]; then
  rmdir "$DESTINO" 2>/dev/null
  echo "*** $DESTINO está em outro sistema de arquivos: o mv copiaria os"
  echo "*** bytes em vez de renomear. Escolha um destino no mesmo disco."
  exit 1
fi
echo ""
echo "  mesmo sistema de arquivos ($DEV_ORI): o mv renomeia, 0 bytes escritos"

if [[ "$APLICAR" == "0" ]]; then
  rmdir "$DESTINO" 2>/dev/null
  echo ""
  echo "  ESTE FOI O PLANO. Para executar:"
  echo "      bash ${BASH_SOURCE[0]##*/} --aplicar"
  exit 0
fi

for a in "${PRESENTES[@]}"; do
  mv "$a" "$DESTINO/" && echo "      movido: $(basename "$a")" \
    || { echo "*** falhou mover $a"; exit 1; }
done
echo ""
echo "  pronto. O driver volta a rodar as fases 5 em diante do zero, com"
echo "  $ATUAL. O que foi arquivado continua em:"
echo "      $DESTINO"
