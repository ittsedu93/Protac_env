#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Só a PREPARAÇÃO do PRosettaC: Init0.pdb, Init1.pdb e Patchdock_params.txt.
#
#   bash scripts/prosettac_prepare.sh <candidate_id> [minutos_max]
#
# Para que isto existe
# --------------------
# O vão que o par E3/alvo exige — a distância mínima entre os dois átomos de
# conjugação para existir alguma colocação sem clash — é o número que deveria
# governar a escolha de linker. Ele se mede com o PatchDock sobre o Init0 e o
# Init1, que são a E3 com o recrutador na pose e o alvo com a warhead na pose.
#
# Só que quem produz esses dois arquivos é o próprio PRosettaC, nos primeiros
# minutos de uma execução que depois leva horas. Este script aproveita esses
# minutos: lança, espera os arquivos aparecerem, e para o processo.
#
# Assim o vão exigido é medido ANTES de escolher linker, que é a ordem certa —
# e não depois de 46 h de MD num candidato impossível, que foi a ordem que a
# CRBN nos impôs.
#
# O PROTAC do config aqui é irrelevante: ele só afeta a restrição de distância
# que o `patchdock_span_scan.sh` sobrescreve ponto a ponto.
# ---------------------------------------------------------------------------
set -uo pipefail

AQUI="$(cd "$(dirname "$0")" && pwd)"
_OUT_ENV="${PIPELINE_OUT:-}"
_PROS_ENV="${PROSETTAC_DIR:-}"
CONF="$AQUI/../config/pipeline.conf"
[[ -f "$CONF" ]] && source "$CONF"
OUT="${_OUT_ENV:-${PIPELINE_OUT:-$HOME/PRosettaC_runs/vhl_crbn_pcsk9_protac/pipeline}}"
PROSETTAC="${_PROS_ENV:-${PROSETTAC_DIR:-/mnt/hd2tb/Documentos/PRosettaC}}"
CFG_NOME="prosetta_config.txt"

CAND="${1:?uso: prosettac_prepare.sh <candidate_id> [minutos_max]}"
MAX_MIN="${2:-20}"

echo "=============================================================="
echo "Preparação do PRosettaC (só Init0/Init1) — $CAND"
echo "=============================================================="

DIR=""
while IFS= read -r achado; do DIR="$(dirname "$achado")"; break
done < <(find "$OUT" -maxdepth 5 -type f -name "$CFG_NOME" -path "*/$CAND/*" \
              2>/dev/null | sort)
[[ -n "$DIR" ]] || { echo "não achei $CFG_NOME para $CAND em $OUT"; exit 1; }
cd "$DIR" || exit 1
echo "  diretório: $DIR"

if [[ -s Init0.pdb && -s Init1.pdb && -s Patchdock_params.txt ]]; then
  echo "  preparação já existe — nada a fazer"
  exit 0
fi

# A validação e o conserto são os mesmos do lançamento de verdade: rodar a
# preparação com head torto ou âncora errada produziria um Init inútil, e o
# erro só apareceria na medição seguinte.
python "$AQUI/prosettac_localize.py" --dir . --config "$CFG_NOME" || exit 1
python "$AQUI/prosettac_clean.py" --dir . --config "$CFG_NOME" --aplicar || exit 1
python "$AQUI/prosettac_fix_heads.py" --dir . --config "$CFG_NOME" --aplicar \
  | grep -E "OK —|CORRIGIDO|casa [0-9]+ átomos" | sed 's/^/  /'
python "$AQUI/prosettac_anchors.py" --dir . --config "$CFG_NOME" --aplicar \
  | grep -E "âncora =|config atualizado" | sed 's/^/  /'

export SETVARS_CALL="${SETVARS_CALL:-}"
echo
echo "[$(date -Is)] preparando (limite de $MAX_MIN min)..."
setsid "$PROSETTAC/run_prosettac.sh" "$DIR" "$CFG_NOME" \
    > prepare.log 2>&1 &
disown -a

pronto=0
for _ in $(seq 1 $((MAX_MIN * 6))); do
  sleep 10
  if [[ -s Init0.pdb && -s Init1.pdb && -s Patchdock_params.txt ]]; then
    pronto=1; break
  fi
  # morreu antes de produzir: não faz sentido esperar o resto do limite
  if ! pgrep -f "$PROSETTAC" > /dev/null \
     && ! pgrep -u "${USER:-$(id -un)}" -f "PRosettaC.*main\.py" > /dev/null; then
    echo "  o processo terminou sem produzir os Init — log:"
    tail -25 prepare.log 2>/dev/null | sed 's/^/      /'
    exit 1
  fi
done

# O resto da execução são horas de Rosetta que não interessam aqui. Parar é
# deliberado: o que este script entrega são os dois Init e o params.
pkill -f "$PROSETTAC/run_prosettac.sh" 2>/dev/null
pkill -u "${USER:-$(id -un)}" -f "PRosettaC.*main\.py" 2>/dev/null
sleep 2

if [[ $pronto -eq 0 ]]; then
  echo "  [TEMPO ESGOTADO] $MAX_MIN min sem Init0/Init1. Últimas linhas:"
  tail -25 prepare.log 2>/dev/null | sed 's/^/      /'
  exit 1
fi

echo "  pronto:"
ls -la Init0.pdb Init1.pdb Patchdock_params.txt | sed 's/^/    /'
echo
echo "  agora a medição do vão exigido:"
echo "    bash scripts/patchdock_span_scan.sh $CAND"
