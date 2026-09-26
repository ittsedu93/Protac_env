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

# --- a GPU responde? -------------------------------------------------------
# O nvidia-smi usa a NVML, que é a biblioteca de GERENCIAMENTO. Ela pode estar
# desalinhada com o módulo do kernel (depois de um apt upgrade sem reboot)
# enquanto a API do CUDA continua funcionando. Quem decide é um contexto CUDA
# de verdade, não o nvidia-smi — e descobrir isso em 10 segundos é melhor que
# descobrir depois de 40 minutos de fila.
CUDA=$(conda run -n "$ENV_BOLTZ" python -c \
  "import torch; print(torch.cuda.is_available())" 2>/dev/null | tail -1)
if [[ "$CUDA" != "True" ]]; then
  echo "*** o CUDA não está disponível no env $ENV_BOLTZ (torch: '$CUDA')."
  echo "*** Se o nvidia-smi também falha com 'Driver/library version"
  echo "*** mismatch', o módulo do kernel e as bibliotecas estão em versões"
  echo "*** diferentes — só um reboot alinha, e ele precisa de root."
  echo "***"
  echo "*** Peça a quem administra a máquina. Enquanto isso, o PRosettaC não"
  echo "*** usa GPU e pode rodar:"
  echo "***   bash $(dirname "$0")/run_prosettac.sh $CAND"
  exit 1
fi
echo "  CUDA disponível"

# --- conferir as flags contra a versão instalada --------------------------
AJUDA=$(conda run -n "$ENV_BOLTZ" boltz predict --help 2>&1)
if [[ -z "$AJUDA" ]] || grep -qi "not found\|No such" <<< "$AJUDA"; then
  echo "*** não consegui rodar 'boltz predict --help' no env $ENV_BOLTZ"
  echo "$AJUDA" | tail -5
  exit 1
fi

OPCOES=()
for flag in --out_dir --use_msa_server --output_format --diffusion_samples \
            --no_kernels --override; do
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
    --no_kernels)         TEM_NO_KERNELS=1 ;;
    --override)           TEM_OVERRIDE=1 ;;
  esac
done
TEM_NO_KERNELS="${TEM_NO_KERNELS:-0}"
TEM_OVERRIDE="${TEM_OVERRIDE:-0}"

# Havendo predição anterior, o Boltz PULA em vez de refazer:
#     Found some existing predictions (1), skipping and running only the
#     missing ones, if any.
# Numa repetição isso é o contrário do que se quer — a repetição existe
# justamente para ver se a geometria se reproduz.
if compgen -G "$DIR/boltz_results_*/predictions/*/confidence_*.json" > /dev/null; then
  n_ant=$(ls "$DIR"/boltz_results_*/predictions/*/confidence_*.json 2>/dev/null | wc -l)
  if [[ "$TEM_OVERRIDE" == "1" ]]; then
    echo "  já há $n_ant predição(ões); --override vai refazer"
    echo "  (melhor_ternario.pdb e ranking_poses.csv ficam fora dessa pasta"
    echo "   e não se perdem)"
    CMD+=(--override)
  else
    echo "*** já há $n_ant predição(ões) e esta versão não tem --override."
    echo "*** Apague a pasta de resultados antes de repetir:"
    echo "***   rm -rf $DIR/boltz_results_*"
    exit 1
  fi
fi

echo "  comando: conda run -n $ENV_BOLTZ ${CMD[*]}"
echo
if grep -q -- "--use_msa_server" <<< "$AJUDA"; then
  echo "  NOTA: --use_msa_server envia as sequências ao servidor de MSA do"
  echo "  ColabFold. PCSK9 e CRBN são proteínas humanas públicas, então não há"
  echo "  questão de confidencialidade — mas a máquina precisa de internet."
  echo
fi

# Escada de kernels. O Boltz tenta os kernels otimizados da NVIDIA
# (cuequivariance) para o triangle multiplicative update; não estando
# instalados ou não casando com o CUDA, ele morre com
#     ImportError: Error importing triangle_multiplicative_update
#                  from cuequivariance_ops_torch
# `--no_kernels` usa a implementação nativa: mais lenta, mesmo resultado.
# Tentar e cair é melhor que decidir de antemão — quando os kernels existem,
# eles valem a pena.
RUNNER="$DIR/.boltz_run.sh"
{
  echo '#!/usr/bin/env bash'
  echo 'set -uo pipefail'
  printf 'conda run -n %q ' "$ENV_BOLTZ"
  printf '%q ' "${CMD[@]}"
  echo
  if [[ "$TEM_NO_KERNELS" == "1" ]]; then
    echo 'codigo=$?'
    echo 'if [[ $codigo -ne 0 ]] && grep -qi "cuequivariance\|triangle_multiplicative" "'"$DIR"'/boltz.log"; then'
    echo '  echo ""'
    echo '  echo "=== os kernels do cuequivariance não carregaram; repetindo com --no_kernels ==="'
    printf '  conda run -n %q ' "$ENV_BOLTZ"
    printf '%q ' "${CMD[@]}"
    printf -- '--no_kernels
'
    echo 'fi'
  fi
} > "$RUNNER"
chmod +x "$RUNNER"

setsid bash "$RUNNER" > "$DIR/boltz.log" 2>&1 &
disown -a
sleep 10
if pgrep -f "boltz predict $YAML" > /dev/null; then
  echo "  rodando — pode desligar o notebook"
  echo "  log:    tail -f $DIR/boltz.log"
  echo "  estado: ls $DIR"
else
  echo "*** o Boltz não ficou rodando. Primeiras linhas do log:"
  head -25 "$DIR/boltz.log" 2>/dev/null
  echo
  if [[ "$TEM_NO_KERNELS" != "1" ]]; then
    echo "*** Esta versão do boltz não tem --no_kernels. Se o erro for do"
    echo "*** cuequivariance, o caminho é instalar o pacote:"
    echo "***   conda run -n $ENV_BOLTZ pip install cuequivariance-ops-torch"
  fi
  echo "*** Se for memória de GPU, tente menos amostras:"
  echo "***   bash $0 $CAND 1"
  exit 1
fi
