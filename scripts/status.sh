#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Onde o pipeline está, agora. Só leitura — seguro rodar a qualquer momento.
#
#   bash scripts/status.sh                                      # track da CRBN
#   PIPELINE_CONF=config/pipeline_vhl.conf bash scripts/status.sh   # VHL
#
# A config decide QUAL track ele olha, e isto é o ponto. A versão anterior
# fixava `config/pipeline.conf` no código e ignorava a variável PIPELINE_CONF
# que todos os outros scripts respeitam: pedir o status da VHL devolvia o
# estado da CRBN, com a mesma cara de certo. Um monitor que responde sobre o
# track errado é pior que nenhum monitor.
#
# Há DUAS formas de execução neste projeto, e as duas são olhadas aqui:
#   - o driver e a MD rodam como processo solto (lançado com setsid) -> pgrep
#   - o PRosettaC emite milhares de jobs no SLURM                   -> squeue
# ---------------------------------------------------------------------------
set -uo pipefail

AQUI="$(cd "$(dirname "$0")" && pwd)"
CONF="${PIPELINE_CONF:-$AQUI/../config/pipeline.conf}"
# Caminho relativo passado na variável (o uso natural: PIPELINE_CONF=config/x)
[[ -f "$CONF" ]] || [[ "$CONF" = /* ]] || CONF="$AQUI/../$CONF"
if [[ -f "$CONF" ]]; then
  # shellcheck disable=SC1090
  source "$CONF"
else
  echo "[aviso] config não encontrada: $CONF — usando os valores default"
fi
OUT="${PIPELINE_OUT:-$HOME/PRosettaC_runs/vhl_crbn_pcsk9_protac/pipeline}"
NS_PROD="${MD_NS_PROD:-200}"
U="${USER:-$(id -un)}"
# shellcheck disable=SC1091
source "$AQUI/proc_prosettac.sh"

echo "=============================================================="
echo "PIPELINE PROTAC PCSK9 — $(date '+%d/%m %H:%M')"
echo "=============================================================="
echo "  config: $(basename "$CONF")"
echo "  saídas: $OUT"
echo "  E3:     ${WP1_E3_LIST:-?}"

# --- 1. está rodando? ------------------------------------------------------
proc=$(pgrep -af "run_pipeline.sh|gmx mdrun|acpype|PRosettaC" | grep -v pgrep || true)
if [[ -n "$proc" ]]; then
  echo -e "\nPROCESSOS ATIVOS:"
  echo "$proc" | cut -c1-110 | sed 's/^/  /'
else
  echo -e "\nNENHUM PROCESSO ATIVO — ou terminou, ou caiu (veja o log adiante)"
fi

# --- 1b. a fila do SLURM, que é onde o PRosettaC vive ---------------------
if command -v squeue >/dev/null 2>&1; then
  rod=$(squeue -h -u "$U" -t RUNNING -o %i 2>/dev/null | wc -l)
  pen=$(squeue -h -u "$U" -t PENDING -o %i 2>/dev/null | wc -l)
  if [[ $((rod + pen)) -gt 0 ]]; then
    echo -e "\nFILA (SLURM):"
    echo "  rodando: $rod   |   pendentes: $pen"
    j=$(squeue -h -u "$U" -o %i 2>/dev/null | head -1)
    cpt=$(scontrol show job "$j" 2>/dev/null \
            | grep -o "CPUs/Task=[0-9]*" | head -1 | cut -d= -f2)
    nuc=$(nproc 2>/dev/null || echo "?")
    if [[ -n "${cpt:-}" && "$nuc" != "?" ]]; then
      sim=$(( nuc / cpt )); [[ $sim -lt 1 ]] && sim=1
      echo "  $cpt CPUs por job em $nuc núcleos -> ~$sim simultâneos"
      # 12 é o que o PRosettaC pede no script de lote, e cada processo usa UM
      # núcleo: com 12 reservados cabem 2 jobs em 32 núcleos.
      [[ "$cpt" -ge 12 ]] && \
        echo "  [ATENÇÃO] 12 CPUs por job é o pedido original e desperdiça a máquina:
            bash scripts/slurm_cpus_por_job.sh 4 --aplicar"

      # Quantos PROCESSOS cada job de fato cria. O pedido de CPU do SLURM é uma
      # promessa; isto é a medida. Se o job cria mais processos do que pediu
      # CPUs, a máquina fica com mais trabalho do que núcleos — e numa máquina
      # COMPARTILHADA essa diferença sai do bolso de outras pessoas, sem
      # acelerar a nossa fila (o sistema só reparte os mesmos núcleos).
      # Não `pgrep -c -f PRosettaC`: esse padrão conta o PRÓPRIO status.sh,
      # porque o caminho do track contém "PRosettaC_runs". Era por isso que ele
      # dizia 19 onde havia 18 — e um número que alimenta recomendação de CPU
      # não pode ter o medidor dentro da medida.
      nproc_pr=$(prosettac_n_trabalhadores)
      carga=$(awk '{printf "%.0f", $1}' /proc/loadavg 2>/dev/null)
      if [[ "$rod" -gt 0 && "$nproc_pr" -gt 0 ]]; then
        por_job=$(awk -v p="$nproc_pr" -v j="$rod" 'BEGIN{printf "%.1f", p/j}')
        echo "  processos do PRosettaC: $nproc_pr  (~$por_job por job)"
        echo "  carga do sistema: ${carga:-?}  em $nuc núcleos"
        if [[ -n "${carga:-}" ]] && [[ "$carga" -gt $((nuc + nuc / 4)) ]]; then
          justo=$(awk -v n="$nuc" -v pj="$por_job" 'BEGIN{v=int(n/pj); print (v<1?1:v)}')
          cpu_justo=$(awk -v n="$nuc" -v j="$justo" 'BEGIN{v=int(n/j); print (v<1?1:v)}')
          echo "  [ATENÇÃO] carga ${carga} em ${nuc} núcleos: a máquina está"
          echo "            sobrecarregada ~$(awk -v c="$carga" -v n="$nuc" 'BEGIN{printf "%.1fx", c/n}')."
          echo "            Cada job cria ~${por_job} processos, não ${cpt}. Para caber"
          echo "            nos núcleos (~${justo} jobs simultâneos), sem perder"
          echo "            vazão nossa e devolvendo a máquina a quem mais a usa:"
          echo "              bash scripts/slurm_cpus_por_job.sh ${cpu_justo} --aplicar --vigiar 20"
        fi
      fi
    fi
  fi
fi

# --- 2. fases concluídas ---------------------------------------------------
# A lista é a MESMA do run_pipeline.sh, e tem de ser: um status que para na
# fase 8 diz "acabou" de um pipeline que vai até a 10.
echo -e "\nFASES:"
fases=(1 2 3 4 5 6 6a 6b 7 8 9 10)
declare -A nome=(
  [1]="warheads PCSK9"            [2]="docking PCSK9"
  [3]="análise da triagem"        [4]="WP1 recrutador"
  [5]="WP2 linkers"               [6]="WP3 montagem + jobs"
  [6a]="vão exigido (PatchDock)"  [6b]="portão de alcance"
  [7]="ranking dos PROTACs"       [8]="MD nível (ii)"
  [9]="ternário (PRosettaC)"      [10]="MD nível (iii)"
)
for n in "${fases[@]}"; do
  if [[ -f "$OUT/.done_$n" ]]; then
    printf "  [x] %-3s %-26s (%s)\n" "$n" "${nome[$n]}" "$(cat "$OUT/.done_$n")"
  else
    printf "  [ ] %-3s %s\n" "$n" "${nome[$n]}"
  fi
done

# --- 2b. a fase 9, que é onde a VHL está ----------------------------------
# O diretório do job vem do prosetta_config.txt e não do nome da pasta: o
# `-name "$CAND"` casa wp3/protacs/<cand> (só o protac.smi) e também
# wp3/prosettac/<cand>, e pegar o primeiro pegava o errado.
# QUAL candidato, e isto importa mais do que parece. A primeira versão disto
# pegava o primeiro prosetta_config.txt em ordem alfabética, e a fase 6 emite um
# config para CADA candidato montado: com SC0006__WH014 e SC0013__WH023 no
# disco, o status anunciava o SC0006 enquanto o que rodava era o SC0013 — e o
# comando de taxa que ele imprimia media o candidato errado. Dois nomes com a
# mesma cara de certo.
#
# A fonte mais verdadeira é o processo VIVO: o cwd do main.py do PRosettaC é o
# diretório do job que está rodando, agora. Sem processo vivo, cai no ranking
# (a mesma fonte que o driver usa) e, por último, no log.txt mais recente.
DIRC=""; FONTE=""
# `pgrep -f` casa QUALQUER processo cuja linha de comando contenha o padrão —
# inclusive o shell que roda este script, um `tail` num caminho com esse nome,
# ou um editor com o arquivo aberto. No teste ele devolveu 4 PIDs e o `head -1`
# pegou o do próprio shell, cujo cwd não é diretório de job nenhum.
# Então não se escolhe por POSIÇÃO: percorre-se os candidatos e aceita-se o
# primeiro cujo cwd realmente contenha um prosetta_config.txt. É verificação
# de conteúdo, não de nome.
for pid_ps in $(pgrep -u "$U" -f "PRosettaC/main\.py" 2>/dev/null); do
  [[ -r "/proc/$pid_ps/cwd" ]] || continue
  cand_dir=$(readlink -f "/proc/$pid_ps/cwd" 2>/dev/null) || continue
  if [[ -n "${cand_dir:-}" && -f "$cand_dir/prosetta_config.txt" ]]; then
    DIRC="$cand_dir"; FONTE="processo vivo (pid $pid_ps)"; break
  fi
done
if [[ -z "$DIRC" && -s "$OUT/protac_ranking.csv" ]]; then
  cand=$(python3 -c "
import csv
r=list(csv.DictReader(open('$OUT/protac_ranking.csv')))
i=min(max(int('${MD_RANK:-1}')-1,0),len(r)-1) if r else 0
print(r[i]['candidate_id'] if r else '')" 2>/dev/null)
  if [[ -n "${cand:-}" ]]; then
    DIRC=$(dirname "$(find "$OUT" -maxdepth 5 -type f -name "prosetta_config.txt" \
             -path "*/$cand/*" 2>/dev/null | sort | head -1)" 2>/dev/null)
    [[ "$DIRC" == "." ]] && DIRC=""
    [[ -n "$DIRC" ]] && FONTE="ranking, rank ${MD_RANK:-1}"
  fi
fi
if [[ -z "$DIRC" ]]; then
  ultimo=$(find "$OUT" -maxdepth 5 -type f -name "log.txt" -path "*/prosettac/*" \
             -printf '%T@ %p\n' 2>/dev/null | sort -rn | head -1 | cut -d' ' -f2-)
  [[ -n "${ultimo:-}" ]] && { DIRC=$(dirname "$ultimo"); FONTE="log.txt mais recente"; }
fi
if [[ -n "$DIRC" ]]; then
  n_cfg=$(find "$OUT" -maxdepth 5 -type f -name "prosetta_config.txt" 2>/dev/null | wc -l)
  echo -e "\nTERNÁRIO — PRosettaC ($(basename "$DIRC")):"
  echo "  candidato identificado por: $FONTE"
  [[ "$n_cfg" -gt 1 ]] && \
    echo "  (há $n_cfg candidatos com job emitido; este é o que está em curso)"
  if [[ -d "$DIRC/Results" ]]; then
    nclu=$(find "$DIRC/Results" -maxdepth 1 -type d ! -path "$DIRC/Results" \
             2>/dev/null | wc -l)
    echo "  [x] Results existe — $nclu cluster(s)"
    [[ -s "$DIRC/result_summary.txt" ]] && \
      sed 's/^/      /' "$DIRC/result_summary.txt"
    [[ -s "$DIRC/concordancia_clusters.csv" ]] \
      && echo "      concordância medida (concordancia_clusters.csv)" \
      || echo "      falta o portão de concordância: --from 9"
  else
    etapa=$(tail -1 "$DIRC/log.txt" 2>/dev/null | cut -c1-70)
    [[ -n "$etapa" ]] && echo "  etapa: $etapa"
    npd=$(find "$DIRC/Patchdock_Results" -maxdepth 1 -type f 2>/dev/null | wc -l)
    [[ "$npd" -gt 0 ]] && echo "  arquivos em Patchdock_Results: $npd"
    echo "  Results ainda não existe — a fase 9 está em andamento"
    echo "  taxa medida:  python scripts/prosettac_progress.py --candidato $(basename "$DIRC")"
  fi
fi

# --- 3. as MDs: nível (ii) em md/, nível (iii) em md_ternario/ -------------
relatorio_md() {
  local MD="$1" rotulo="$2"
  [[ -d "$MD" ]] || return 0
  echo -e "\n$rotulo  ($MD)"

  echo "  PREPARO:"
  local RESNAME; RESNAME=$(python3 -c "
import json;print(json.load(open('$MD/md_sistema.json')).get('resname','PTC'))" \
    2>/dev/null || echo PTC)
  for f in md_sistema.json protac.sdf "${RESNAME}.acpype/${RESNAME}_GMX.itp" \
           topol.top complexo.gro neutro.gro grupos.ndx; do
    printf "    %s %s\n" "$([[ -s "$MD/$f" ]] && echo '[x]' || echo '[ ]')" "$f"
  done

  # O equilíbrio é uma escada: cada degrau alimenta o próximo, e depois do
  # npt.gro nenhum dos anteriores é necessário. Listar arquivo por arquivo
  # responde "o que existe no disco" quando a pergunta é "em que estado o
  # pipeline está" — e um intermediário ausente parece falta quando é sobra.
  echo "  EQUILÍBRIO:"
  if [[ -s "$MD/npt.gro" ]]; then
    echo "    [x] CONCLUÍDO — npt.gro é a estrutura de partida da produção"
    local faltando=()
    for f in em.gro em2.gro warm.gro nvt.gro; do
      [[ -s "$MD/$f" ]] || faltando+=("$f")
    done
    (( ${#faltando[@]} )) && {
      echo "        intermediários ausentes: ${faltando[*]}"
      echo "        (consumidos pelo degrau seguinte; não fazem falta)"; }
  else
    for f in em.gro em2.gro warm.gro nvt.gro npt.gro; do
      printf "    %s %s\n" "$([[ -s "$MD/$f" ]] && echo '[x]' || echo '[ ]')" "$f"
    done
  fi

  # O mdrun escreve "remaining wall clock time" enquanto roda. É a única
  # estimativa honesta de quanto falta: vem do desempenho medido.
  local f linha nomef dir
  for f in "$MD"/mdrun_*.out "$MD"/rep*/mdrun_prod.out; do
    [[ -f "$f" ]] || continue
    linha=$(grep -a "remaining wall clock" "$f" | tail -1)
    [[ -z "$linha" ]] && continue
    nomef=$(basename "$f" .out); nomef=${nomef#mdrun_}
    dir=$(basename "$(dirname "$f")")
    [[ "$dir" == rep* ]] && nomef="$dir"
    echo "    ${nomef}: ${linha}"
  done

  if compgen -G "$MD/rep*" > /dev/null; then
    echo "  PRODUÇÃO (${NS_PROD} ns por réplica):"
    local rep n perf horas passo tempo resta pct
    for rep in "$MD"/rep*/; do
      n=$(basename "$rep")
      if [[ -s "$rep/prod.gro" ]]; then
        perf=$(awk '/^Performance:/ {print $2}' "$rep/prod.log" 2>/dev/null | tail -1)
        if [[ -n "$perf" ]]; then
          horas=$(awk -v p="$perf" -v ns="$NS_PROD" 'BEGIN{printf "%.1f", ns/p*24}')
          echo "    $n: CONCLUÍDA — ${perf} ns/dia (${horas} h)"
        else
          echo "    $n: CONCLUÍDA"
        fi
        continue
      fi
      if [[ -f "$rep/prod.log" ]]; then
        linha=$(grep -A1 "^ *Step  *Time" "$rep/prod.log" | tail -1)
        passo=$(echo "$linha" | awk '{print $1}')
        tempo=$(echo "$linha" | awk '{print $2}')
        resta=$(awk -v t="${tempo:-0}" -v ns="$NS_PROD" 'BEGIN{printf "%.1f", ns - t/1000}')
        pct=$(awk -v t="${tempo:-0}" -v ns="$NS_PROD" 'BEGIN{printf "%.1f", t/(ns*10)}')
        echo "    $n: ${pct}% — passo ${passo:-?}, ${tempo:-?} ps (faltam ~${resta} ns)"
        perf=$(grep -E "^Performance:" "$rep/prod.log" | tail -1)
        [[ -n "$perf" ]] && echo "        $perf"
      else
        echo "    $n: iniciando"
      fi
    done
  fi

  for v in md_veredito.json ternario_metricas.csv; do
    [[ -s "$MD/$v" ]] && { echo "  VEREDITO ($v):";
                           head -12 "$MD/$v" | sed 's/^/    /'; }
  done
}

relatorio_md "$OUT/md"           "MD NÍVEL (ii) — E3 + recrutador-linker-warhead"
relatorio_md "$OUT/md_ternario"  "MD NÍVEL (iii) — ternário completo com PCSK9"

# --- 5. GPU ----------------------------------------------------------------
if command -v nvidia-smi > /dev/null; then
  echo -e "\nGPU:"
  nvidia-smi --query-gpu=name,utilization.gpu,memory.used,temperature.gpu \
      --format=csv,noheader | sed 's/^/  /'
fi

# --- 6. o log do driver ----------------------------------------------------
# O log mais recente DESTE track, não um caminho fixo no $HOME: com dois
# tracks o `~/pipeline.log` é de um deles, e não se sabe de qual.
LOG=$(ls -t "$OUT"/pipeline_*.log 2>/dev/null | head -1)
[[ -z "$LOG" ]] && LOG=$(ls -t "$HOME"/pipeline*.log 2>/dev/null | head -1)
if [[ -n "${LOG:-}" && -f "$LOG" ]]; then
  echo -e "\nÚLTIMAS LINHAS DE $(basename "$LOG"):"
  grep -vE "^\s*$" "$LOG" | tail -6 | cut -c1-100 | sed 's/^/  /'
  falha=$(grep -n "^\*\*\*" "$LOG" | tail -1)
  [[ -n "$falha" ]] && echo -e "\n  ÚLTIMA FALHA REGISTRADA: $(echo "$falha" | cut -c1-100)"
fi
echo
