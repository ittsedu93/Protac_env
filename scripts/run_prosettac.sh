#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Dispara o PRosettaC para o candidato aprovado no nível (ii).
#
#   bash scripts/run_prosettac.sh <candidate_id> [--lote]
#
# Sem --lote roda UM job e para. Isso é deliberado e está no RUNBOOK desde o
# começo: a convenção de `Anchor atoms` (0-based vs 1-based) varia por build do
# PRosettaC, e um lote inteiro com o átomo errado é um lote perdido. O primeiro
# job existe para conferir isso.
#
# Retomável: um candidato cujo diretório já tem resultado é pulado.
# ---------------------------------------------------------------------------
set -uo pipefail

CONF="$(dirname "$0")/../config/pipeline.conf"
[[ -f "$CONF" ]] && source "$CONF"
OUT="${PIPELINE_OUT:-$HOME/PRosettaC_runs/vhl_crbn_pcsk9_protac/pipeline}"
PROSETTAC="${PROSETTAC_DIR:-/mnt/hd2tb/Documentos/PRosettaC}"

CAND="${1:?uso: run_prosettac.sh <candidate_id> [--lote]}"
LOTE="${2:-}"
# A fase 6 escreve os jobs em $PIPELINE_OUT/wp3/<candidato> (o --outdir do
# wp3_assemble_protacs.py), não na raiz do pipeline.
DIR="$OUT/wp3/$CAND"
[[ -d "$DIR" ]] || DIR="$OUT/$CAND"

[[ -d "$PROSETTAC" ]] || { echo "não achei o PRosettaC em $PROSETTAC"; exit 1; }
[[ -x "$PROSETTAC/run_prosettac.sh" ]] || {
  echo "não achei $PROSETTAC/run_prosettac.sh (executável)"; exit 1; }
[[ -s "$DIR/prosetta_config.txt" ]] || {
  echo "não achei $DIR/prosetta_config.txt"
  echo "ele é escrito pela fase 6 (wp3_assemble_protacs.py)."
  echo "Candidatos com job emitido:"
  ls -1 "$OUT/wp3" 2>/dev/null | sed 's/^/    /' | head -20
  exit 1; }

echo "=============================================================="
echo "PRosettaC — $CAND"
echo "=============================================================="
echo
echo "config que será usado:"
sed 's/^/    /' "$DIR/prosetta_config.txt"
echo

# --- conferências antes de gastar horas ------------------------------------
falta=0
while read -r chave valores; do
  case "$chave" in
    Structures:|Heads:|Protac:)
      for f in $valores; do
        if [[ ! -s "$f" ]]; then echo "  [FALTA] $chave $f"; falta=1; fi
      done ;;
  esac
done < "$DIR/prosetta_config.txt"
[[ $falta -eq 0 ]] || { echo -e "\n*** arquivos de entrada ausentes"; exit 1; }

ANCHORS=$(awk -F': ' '/^Anchor atoms:/{print $2}' "$DIR/prosetta_config.txt")
echo "  Anchor atoms = $ANCHORS"
echo "  >>> CONFIRA no resultado deste job que os âncoras são os átomos que"
echo "  >>> ligam cada head ao linker. 0-based vs 1-based varia por build,"
echo "  >>> e o lote inteiro depende disto."
echo

if [[ -d "$DIR/Results" || -d "$DIR/results" ]]; then
  echo "  já há resultado em $DIR — nada a fazer"
  exit 0
fi

cd "$DIR"
echo "[$(date -Is)] lançando..."
setsid "$PROSETTAC/run_prosettac.sh" "$DIR" prosetta_config.txt \
    > "$DIR/run_prosettac.log" 2>&1 &
disown -a
sleep 5
echo "  PID $(pgrep -f "run_prosettac.sh $DIR" | head -1) — pode desligar o notebook"
echo "  log: tail -f $DIR/run_prosettac.log"

if [[ "$LOTE" == "--lote" ]]; then
  echo
  echo "*** --lote foi ignorado de propósito nesta execução."
  echo "*** Confira o Anchor atoms deste job PRIMEIRO; depois rode"
  echo "***   bash $OUT/wp3/launch_all.sh"
fi
