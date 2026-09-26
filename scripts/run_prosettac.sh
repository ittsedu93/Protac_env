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
# Onde a fase 6 escreveu os jobs é coisa que se PROCURA, não que se adivinha:
# já errei este caminho duas vezes seguidas supondo o layout em vez de olhar.
# O arquivo que importa é o prosetta_config.txt do candidato; achá-lo pelo
# nome funciona em qualquer layout, hoje e depois de qualquer refatoração.
DIR=""
while IFS= read -r achado; do
  DIR="$(dirname "$achado")"; break
done < <(find "$OUT" -maxdepth 5 -type f -name prosetta_config.txt \
             -path "*/$CAND/*" 2>/dev/null | sort)
[[ -n "$DIR" ]] || DIR="$OUT/wp3/prosettac/$CAND"

[[ -d "$PROSETTAC" ]] || { echo "não achei o PRosettaC em $PROSETTAC"; exit 1; }
[[ -x "$PROSETTAC/run_prosettac.sh" ]] || {
  echo "não achei $PROSETTAC/run_prosettac.sh (executável)"; exit 1; }
[[ -s "$DIR/prosetta_config.txt" ]] || {
  echo "não achei $DIR/prosetta_config.txt"
  echo "ele é escrito pela fase 6 (wp3_assemble_protacs.py)."
  echo
  echo "Candidatos com job emitido (procurados por prosetta_config.txt):"
  find "$OUT" -maxdepth 5 -type f -name prosetta_config.txt 2>/dev/null \
    | sed 's|/prosetta_config.txt$||' | xargs -r -n1 basename \
    | sort -u | sed 's/^/    /' | head -20
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
sleep 8

# O setsid forka, então $! não é o processo final e um pgrep pelo comando
# completo erra. Quem responde "está vivo?" é o log crescer.
if pgrep -f "prosetta_config.txt" > /dev/null || [[ -s "$DIR/run_prosettac.log" ]]; then
  echo "  rodando — pode desligar o notebook"
else
  echo "  [ATENÇÃO] nada rodando e o log está vazio."
  echo "  O PRosettaC submete a um gerenciador de filas; o config pede SLURM."
  echo "  Se esta máquina não tem SLURM, ele não tem onde submeter:"
  command -v sbatch > /dev/null \
    && echo "    sbatch existe: $(command -v sbatch)" \
    || echo "    sbatch NÃO existe nesta máquina — é provável que seja isso"
fi
echo "  log: tail -f $DIR/run_prosettac.log"

if [[ "$LOTE" == "--lote" ]]; then
  echo
  echo "*** --lote foi ignorado de propósito nesta execução."
  echo "*** Confira o Anchor atoms deste job PRIMEIRO; depois rode"
  echo "***   bash $(dirname "$DIR")/launch_all.sh"
fi
