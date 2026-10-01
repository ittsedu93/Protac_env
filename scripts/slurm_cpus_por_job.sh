#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Ajusta o pedido de CPU dos jobs PENDENTES do PRosettaC.
#
#   bash scripts/slurm_cpus_por_job.sh            # mostra o que faria
#   bash scripts/slurm_cpus_por_job.sh 2 --aplicar
#   bash scripts/slurm_cpus_por_job.sh 2 --aplicar --vigiar 20
#
# O problema que isto resolve
# ---------------------------
# O PRosettaC gera os scripts de lote com
#
#     #SBATCH --cpus-per-task=12
#
# e os processos que eles lançam usam UM núcleo cada:
#
#     859316  99.9% CPU  12 threads      <- 99,9% = um núcleo
#     859307  99.9% CPU  12 threads
#
# Com 12 reservados por job e 32 núcleos na máquina, cabem dois jobs por vez.
# Dois jobs de dois núcleos úteis cada: 4 dos 32 trabalham, e 2553 esperam.
# Medido no track da VHL: 125 jobs em 18 h, ou ~15 dias para os 2682.
#
# O pedido de job PENDENTE é alterável em voo, sem reiniciar nada e sem ser
# administrador. Com 2 por job cabem 16 simultâneos, e os 15 dias viram ~2.
#
# Por que não 1
# -------------
# Cada job roda DOIS processos de um núcleo (4 processos para 2 jobs, na
# medição). Pedir 1 faria dois processos disputarem um núcleo: a fila anda mais
# rápido e cada job anda mais devagar, sem ganho líquido. O número certo é o que
# o job usa, e ele se mede com `ps -o %cpu` — não se adivinha.
#
# Por que repetir (--vigiar)
# --------------------------
# O PRosettaC emite os jobs em LOTES, uma vez por etapa: docking local, depois
# conformações do linker, depois o agrupamento. Os lotes seguintes nascem com os
# 12 de novo. Em modo vigia, isto reaplica a cada N minutos.
# ---------------------------------------------------------------------------
set -uo pipefail

CPUS="${1:-2}"
[[ "$CPUS" =~ ^[0-9]+$ ]] || { echo "uso: $0 <cpus_por_job> [--aplicar] [--vigiar <min>]"; exit 1; }
shift || true

APLICAR=0
VIGIAR=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --aplicar) APLICAR=1; shift;;
    --vigiar) VIGIAR="${2:-20}"; shift 2;;
    *) echo "opção desconhecida: $1"; exit 1;;
  esac
done

U="${USER:-$(id -un)}"
command -v squeue >/dev/null || { echo "sem squeue nesta máquina"; exit 1; }

uma_passada() {
  local pendentes rodando total_cpu ajustados=0
  pendentes=$(squeue -h -u "$U" -t PENDING -o %i | wc -l)
  rodando=$(squeue -h -u "$U" -t RUNNING -o %i | wc -l)
  if [[ "$pendentes" -eq 0 ]]; then
    echo "  [$(date +%H:%M)] nenhum job pendente"
    return 0
  fi

  # quanto cada job pede hoje, pelo primeiro pendente
  local j atual
  j=$(squeue -h -u "$U" -t PENDING -o %i | head -1)
  atual=$(scontrol show job "$j" 2>/dev/null \
            | grep -o "NumCPUs=[0-9]*" | head -1 | cut -d= -f2)
  echo "  [$(date +%H:%M)] rodando: $rodando | pendentes: $pendentes"
  echo "              pedido atual: ${atual:-?} CPUs por job -> pedido novo: $CPUS"

  if [[ "${atual:-0}" == "$CPUS" ]]; then
    echo "              já estão em $CPUS — nada a fazer"
    return 0
  fi
  if [[ $APLICAR -eq 0 ]]; then
    echo "              (simulação; para aplicar: $0 $CPUS --aplicar)"
    return 0
  fi

  # Um job primeiro: se o SLURM recusar a alteração, não adianta tentar 2553.
  if ! scontrol update job "$j" NumCPUs="$CPUS" MinCPUsNode="$CPUS" 2>/dev/null; then
    echo "              [FALHOU] o SLURM recusou a alteração em $j."
    echo "              Alguns clusters proíbem reduzir CPU de job pendente;"
    echo "              nesse caso o caminho é o administrador, ou relançar com"
    echo "              o script de lote corrigido."
    return 1
  fi
  local conferido
  conferido=$(scontrol show job "$j" | grep -o "NumCPUs=[0-9]*" | head -1 | cut -d= -f2)
  if [[ "$conferido" != "$CPUS" ]]; then
    echo "              [NÃO PEGOU] $j continua com $conferido CPUs — parei aqui"
    return 1
  fi

  echo "              funcionou em $j; aplicando nos demais..."
  squeue -h -u "$U" -t PENDING -o %i \
    | xargs -P 4 -I{} scontrol update job {} NumCPUs="$CPUS" MinCPUsNode="$CPUS" \
      2>/dev/null
  ajustados=$(squeue -h -u "$U" -t PENDING -o "%i %C" | awk -v c="$CPUS" '$2==c' | wc -l)
  echo "              $ajustados de $pendentes pendentes agora com $CPUS CPUs"
  local nucleos antes depois
  nucleos=$(nproc)
  depois=$(( nucleos / CPUS )); [[ $depois -lt 1 ]] && depois=1
  antes=$(( nucleos / ${atual:-12} )); [[ $antes -lt 1 ]] && antes=1
  echo "              simultâneos esperados: ~$depois (eram ~$antes)"
}

echo "=============================================================="
echo "CPUs por job do PRosettaC — $U"
echo "=============================================================="
uma_passada
if [[ "$VIGIAR" -gt 0 ]]; then
  echo
  echo "  vigiando a cada $VIGIAR min (os lotes seguintes nascem com o pedido"
  echo "  antigo; Ctrl+C para sair, ou deixe destacado com setsid)"
  while true; do
    sleep $((VIGIAR * 60))
    squeue -h -u "$U" 2>/dev/null | grep -q . || {
      echo "  [$(date +%H:%M)] fila vazia — encerrando a vigia"; break; }
    uma_passada
  done
fi
