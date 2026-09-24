#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Onde o pipeline está, agora. Só leitura — seguro rodar a qualquer momento.
#
#   bash scripts/status.sh
#
# Não existe fila (SLURM) aqui: o pipeline roda como processo solto, lançado
# com setsid. Quem responde "está rodando?" é o pgrep, não o squeue.
# ---------------------------------------------------------------------------
set -uo pipefail
CONF="$(dirname "$0")/../config/pipeline.conf"
[[ -f "$CONF" ]] && source "$CONF"
OUT="${PIPELINE_OUT:-$HOME/PRosettaC_runs/vhl_crbn_pcsk9_protac/pipeline}"
MD="$OUT/md"
LOG_LAB="$HOME/pipeline.log"

echo "=============================================================="
echo "PIPELINE PROTAC PCSK9 — $(date '+%d/%m %H:%M')"
echo "=============================================================="

# --- 1. está rodando? ------------------------------------------------------
proc=$(pgrep -af "run_pipeline.sh|gmx mdrun|acpype" | grep -v pgrep || true)
if [[ -n "$proc" ]]; then
  echo -e "\nRODANDO:"
  echo "$proc" | sed 's/^/  /'
else
  echo -e "\nNENHUM PROCESSO ATIVO — ou terminou, ou caiu (veja o log adiante)"
fi

# --- 2. fases concluídas ---------------------------------------------------
echo -e "\nFASES:"
nomes=("" "warheads" "docking PCSK9" "análise" "WP1 recrutador" \
       "WP2 linkers" "WP3 montagem" "ranking" "MD")
for n in 1 2 3 4 5 6 7 8; do
  if [[ -f "$OUT/.done_$n" ]]; then
    printf "  [x] %d %s  (%s)\n" "$n" "${nomes[$n]}" "$(cat "$OUT/.done_$n")"
  else
    printf "  [ ] %d %s\n" "$n" "${nomes[$n]}"
  fi
done

# --- 3. dentro da MD -------------------------------------------------------
echo -e "\nMD, etapa por etapa:"
for f in protac.pdb complexo.pdb PTC.acpype/PTC_GMX.itp topol.top complexo.gro \
         neutro.gro grupos.ndx em.gro nvt.gro npt.gro; do
  printf "  %s %s\n" "$([[ -s "$MD/$f" ]] && echo '[x]' || echo '[ ]')" "$f"
done

# --- 4. réplicas de produção ----------------------------------------------
if compgen -G "$MD/rep*" > /dev/null; then
  echo -e "\nPRODUÇÃO (200 ns por réplica):"
  for rep in "$MD"/rep*/; do
    nome=$(basename "$rep")
    if [[ -s "$rep/prod.gro" ]]; then
      echo "  $nome: CONCLUÍDA"
      continue
    fi
    if [[ -f "$rep/prod.log" ]]; then
      # o cabeçalho "Step Time" vem numa linha e os valores na seguinte
      linha=$(grep -A1 "^ *Step  *Time" "$rep/prod.log" | tail -1)
      passo=$(echo "$linha" | awk '{print $1}')
      tempo=$(echo "$linha" | awk '{print $2}')
      resta=$(awk -v t="${tempo:-0}" 'BEGIN{printf "%.1f", 200 - t/1000}')
      pct=$(awk -v t="${tempo:-0}" 'BEGIN{printf "%.1f", t/2000}')
      echo "  $nome: ${pct}% — passo ${passo:-?}, ${tempo:-?} ps (faltam ~${resta} ns)"
      # desempenho, se o mdrun já tiver escrito
      perf=$(grep -E "^Performance:" "$rep/prod.log" | tail -1)
      [[ -n "$perf" ]] && echo "      $perf"
    else
      echo "  $nome: iniciando"
    fi
  done
fi

# --- 5. GPU ----------------------------------------------------------------
if command -v nvidia-smi > /dev/null; then
  echo -e "\nGPU:"
  nvidia-smi --query-gpu=name,utilization.gpu,memory.used,temperature.gpu \
      --format=csv,noheader | sed 's/^/  /'
fi

# --- 6. veredito ou últimas linhas do log ---------------------------------
if [[ -s "$MD/md_veredito.json" ]]; then
  echo -e "\nVEREDITO DA MD:"
  sed 's/^/  /' "$MD/md_veredito.json"
elif [[ -f "$LOG_LAB" ]]; then
  echo -e "\nÚLTIMAS LINHAS DO LOG:"
  grep -vE "^\s*$" "$LOG_LAB" | tail -6 | sed 's/^/  /'
  falha=$(grep -n "^\*\*\*" "$LOG_LAB" | tail -1)
  [[ -n "$falha" ]] && echo -e "\n  ÚLTIMA FALHA REGISTRADA: $falha"
fi
echo
