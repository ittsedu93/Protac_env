#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Predição do complexo ternário com Boltz-2 (GPU), destacada do terminal.
#
#   bash scripts/run_boltz_ternary.sh <candidate_id> [n_amostras]
#
# Por que Boltz-2 e não AlphaFold 3: o AF3 público (alphafoldserver.com) aceita
# ligantes de uma lista fechada, não SMILES arbitrário, e o AF3 local exige
# pedir os pesos. O Boltz-2 está instalado aqui, aceita SMILES, e dá a mesma
# ortogonalidade que o WP3 pede do AF3 — previsão por aprendizado profundo
# contra a amostragem física do PRosettaC.
#
# As FLAGS são conferidas contra o `--help` da versão instalada antes de
# lançar. Rodar 40 minutos para descobrir que uma flag não existe é o erro que
# este projeto já cometeu de outras formas.
# ---------------------------------------------------------------------------
set -uo pipefail

CONF="$(dirname "$0")/../config/pipeline.conf"
[[ -f "$CONF" ]] && source "$CONF"
WORK="${WORK:-$HOME/PRosettaC_runs/vhl_crbn_pcsk9_protac}"
ENV_BOLTZ="${ENV_BOLTZ:-gustavo_boltz-2}"

CAND="${1:?uso: run_boltz_ternary.sh <candidate_id> [n_amostras]}"
N_AMOSTRAS="${2:-5}"
DIR="$WORK/boltz_ternario/$CAND"
YAML="$DIR/ternary.yaml"

[[ -s "$YAML" ]] || {
  echo "não achei $YAML"
  echo "rode antes: python scripts/boltz_ternary.py --candidato $CAND"; exit 1; }

echo "=============================================================="
echo "Boltz-2 — complexo ternário — $CAND"
echo "=============================================================="
grep -c "protein:" "$YAML" | xargs echo "  cadeias de proteína:"
echo "  amostras pedidas: $N_AMOSTRAS"
echo

# --- conferir as flags contra a versão instalada --------------------------
AJUDA=$(conda run -n "$ENV_BOLTZ" boltz predict --help 2>&1)
if [[ -z "$AJUDA" ]] || grep -qi "not found\|No such" <<< "$AJUDA"; then
  echo "*** não consegui rodar 'boltz predict --help' no env $ENV_BOLTZ"
  echo "$AJUDA" | tail -5
  exit 1
fi

OPCOES=()
for flag in --out_dir --use_msa_server --output_format --diffusion_samples; do
  if grep -q -- "$flag" <<< "$AJUDA"; then
    OPCOES+=("$flag")
  else
    echo "  [ausente nesta versão] $flag — não será usada"
  fi
done

CMD=(boltz predict "$YAML")
for flag in "${OPCOES[@]}"; do
  case "$flag" in
    --out_dir)            CMD+=(--out_dir "$DIR") ;;
    --use_msa_server)     CMD+=(--use_msa_server) ;;
    --output_format)      CMD+=(--output_format pdb) ;;
    --diffusion_samples)  CMD+=(--diffusion_samples "$N_AMOSTRAS") ;;
  esac
done

echo "  comando: conda run -n $ENV_BOLTZ ${CMD[*]}"
echo
if grep -q -- "--use_msa_server" <<< "$AJUDA"; then
  echo "  NOTA: --use_msa_server envia as sequências ao servidor de MSA do"
  echo "  ColabFold. PCSK9 e CRBN são proteínas humanas públicas, então não há"
  echo "  questão de confidencialidade — mas a máquina precisa de internet."
  echo
fi

setsid conda run -n "$ENV_BOLTZ" "${CMD[@]}" > "$DIR/boltz.log" 2>&1 &
disown -a
sleep 8
if pgrep -f "boltz predict $YAML" > /dev/null; then
  echo "  rodando — pode desligar o notebook"
  echo "  log:    tail -f $DIR/boltz.log"
  echo "  estado: ls $DIR"
else
  echo "*** o Boltz não ficou rodando. Primeiras linhas do log:"
  head -25 "$DIR/boltz.log" 2>/dev/null
  echo
  echo "*** Se for memória de GPU, tente menos amostras:"
  echo "***   bash $0 $CAND 1"
  exit 1
fi
