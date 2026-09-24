#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Driver do pipeline PROTAC PCSK9 — roda do começo ao fim sem supervisão.
#
#   bash scripts/run_pipeline.sh                # tudo que estiver pendente
#   bash scripts/run_pipeline.sh --from 4       # retoma da fase 4
#   bash scripts/run_pipeline.sh --only 3       # só a fase 3
#   bash scripts/run_pipeline.sh --list         # o que já está feito
#   bash scripts/run_pipeline.sh --dry-run      # imprime sem executar
#
# Para deixar rodando e desligar o notebook:
#   tmux new -s protac
#   bash ~/Protac_env/scripts/run_pipeline.sh 2>&1 | tee ~/pipeline.log
#   (Ctrl+B, depois D para soltar a sessão)
#   tmux attach -t protac        # para voltar
#
# Cada fase grava um marcador em $PIPELINE_OUT/.done_<n>. Fases já concluídas
# são puladas, então relançar depois de uma queda custa só a fase interrompida.
# Os scripts pesados são retomáveis por conta própria, então mesmo uma fase
# interrompida no meio retoma de onde parou.
# ---------------------------------------------------------------------------
set -uo pipefail

CONF="${PIPELINE_CONF:-$(dirname "$0")/../config/pipeline.conf}"
[[ -f "$CONF" ]] || { echo "config não encontrada: $CONF"; exit 1; }
# shellcheck disable=SC1090
source "$CONF"
OBABEL_EXE="${OBABEL_EXE:-/home/soberano/miniconda3/envs/obabel_env/bin/obabel}"

FROM=1; ONLY=""; DRY=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --from) FROM="$2"; shift 2;;
    --only) ONLY="$2"; shift 2;;
    --dry-run) DRY=1; shift;;
    --list) LIST=1; shift;;
    *) echo "opção desconhecida: $1"; exit 1;;
  esac
done

mkdir -p "$PIPELINE_OUT"
LOG="$PIPELINE_OUT/pipeline_$(date +%Y%m%d_%H%M%S).log"

log()  { echo -e "$*" | tee -a "$LOG"; }
head_() { log "\n$(printf '=%.0s' {1..70})"; log "FASE $1 — $2"; log "$(printf '=%.0s' {1..70})"; }
done_marker() { echo "$PIPELINE_OUT/.done_$1"; }
is_done() { [[ -f "$(done_marker "$1")" ]]; }
mark_done() { date -Is > "$(done_marker "$1")"; }

# roda um comando num env conda, com log e parada em erro
run_in() {
  local env="$1"; shift
  log "\n[\$ $*]  (env $env)"
  if [[ $DRY -eq 1 ]]; then return 0; fi
  if ! conda run --no-capture-output -n "$env" "$@" 2>&1 | tee -a "$LOG"; then
    log "\n*** FALHOU na fase atual. O log completo está em $LOG"
    log "*** Corrija e relance com --from <n>; o que já rodou é aproveitado."
    exit 1
  fi
}

exige() {  # exige <variável> <mensagem>
  local nome="$1"; local val="${!1:-}"
  if [[ -z "$val" ]]; then
    log "\n*** Falta preencher '$nome' em $CONF"
    log "*** $2"
    exit 2
  fi
}

if [[ "${LIST:-0}" == "1" ]]; then
  echo "Fases concluídas:"
  for n in 1 2 3 4 5 6 7; do
    if is_done "$n"; then echo "  [x] fase $n  ($(cat "$(done_marker "$n")"))"
    else echo "  [ ] fase $n"; fi
  done
  exit 0
fi

quer() {  # quer <n> -> deve rodar esta fase?
  local n="$1"
  [[ -n "$ONLY" ]] && { [[ "$ONLY" == "$n" ]]; return; }
  (( n >= FROM )) && ! is_done "$n"
}

log "pipeline iniciado $(date -Is)"
log "config: $CONF"
log "log:    $LOG"

# ===========================================================================
# FASE 1 — warheads PCSK9
# ===========================================================================
if quer 1; then
  head_ 1 "Enumerar os warheads PCSK9 (~20 min)"
  run_in "$ENV_MDTOOLS" python "$REPO/scripts/generate_pcsk9_warheads.py" \
      --outdir "$WARHEADS_DIR" --n-confs "$WARHEADS_NCONFS"
  mark_done 1
else
  log "\n[fase 1 pulada]"
fi

# ===========================================================================
# FASE 2 — PDBQT + receptor PCSK9 + validação + triagem
# ===========================================================================
if quer 2; then
  head_ 2 "Ancorar os warheads na PCSK9 (horas, GPU)"

  [[ -f "$PCSK9_PDB" ]] || {
    log "baixando o cristal da PCSK9"
    mkdir -p "$STRUCTURES"
    wget -q -O "$PCSK9_PDB" "https://files.rcsb.org/download/$(basename "${PCSK9_PDB%.pdb}").pdb" \
      || { log "*** falha ao baixar $PCSK9_PDB"; exit 1; }
  }

  run_in "$ENV_PFVS" bash "$WARHEADS_DIR/prepare_pdbqt.sh"

  run_in "$ENV_MDTOOLS" python "$REPO/scripts/prep_pcsk9_receptor.py" \
      --pdb-file "$PCSK9_PDB" --outdir "$DOCKING_DIR" \
      --target-chain "$PCSK9_CHAIN" --keep-chains $PCSK9_KEEP_CHAINS \
      --site-mode ligand --ref-ligand-resname "$PCSK9_REF_LIGAND" \
      --burial-threshold "$PCSK9_BURIAL"

  # o portão: se o redocking reprovar, o script sai != 0 e o driver para aqui
  run_in "$ENV_PFVS" python -u "$REPO/scripts/dock_warheads_pcsk9.py" \
      --site "$DOCKING_DIR/pcsk9_site.json" \
      --warheads-pdbqt "$WARHEADS_DIR/pdbqt" \
      --warheads-csv "$WARHEADS_DIR/pcsk9_warheads.csv" \
      --series "$DOCK_SERIES" --batch-size "$DOCK_BATCH" \
      --round1-cutoff "$DOCK_ROUND1_CUTOFF" \
      --outdir "$DOCKING_DIR/docking"
  mark_done 2
else
  log "\n[fase 2 pulada]"
fi

# ===========================================================================
# FASE 3 — análise e escolha dos warheads
# ===========================================================================
if quer 3; then
  head_ 3 "Analisar a triagem e selecionar os warheads (segundos)"
  run_in "$ENV_MDTOOLS" python "$REPO/scripts/analyze_warhead_docking.py" \
      --docking "$DOCKING_DIR/docking" \
      --site "$DOCKING_DIR/pcsk9_site.json" \
      --warheads-csv "$WARHEADS_DIR/pcsk9_warheads.csv" \
      --sd-max "$ANALISE_SD_MAX" \
      --out "$PIPELINE_OUT/warhead_ranking.csv"
  mark_done 3
else
  log "\n[fase 3 pulada]"
fi

# ===========================================================================
# FASE 4 — WP1: revalidar o portão e escolher o recrutador
# ===========================================================================
if quer 4; then
  head_ 4 "WP1 — revalidar o redocking e escolher o recrutador (minutos)"

  # 4a. revalida o portão de redocking de cada E3 com o RMSD correto.
  #     Não interrompe o pipeline: o resultado é informativo e fica no log e
  #     no CSV, para você decidir o que fazer com a triagem do WP1.
  for E3 in $WP1_E3_LIST; do
    D=$(find "$WP1_PREP" -maxdepth 1 -type d -name "${E3}*" | head -1)
    [[ -z "$D" ]] && { log "  [${E3}] pasta não encontrada em $WP1_PREP"; continue; }
    REF=$(find "$D" -name "*ref_ligand*.sdf" | head -1)
    REC=$(find "$D" -name "*_receptor.pdb" | head -1)
    if [[ -z "$REF" ]]; then
      REFPDB=$(find "$D" -name "*ref_ligand*.pdb" | head -1)
      if [[ -n "$REFPDB" ]]; then
        log "  [${E3}] convertendo ligante de referência para SDF"
        [[ $DRY -eq 0 ]] && "$OBABEL_EXE" "$REFPDB" -O "${REFPDB%.pdb}.sdf" -h \
          >/dev/null 2>&1 && REF="${REFPDB%.pdb}.sdf"
      fi
    fi
    if [[ -n "$REF" && -n "$REC" ]]; then
      run_in "$ENV_MDTOOLS" python "$REPO/scripts/revalidate_redocking.py" \
          --receptor-pdb "$REC" --ref-ligand "$REF" \
          --poses "$D/redock*.pdbqt" --label "$E3" \
          --out "$PIPELINE_OUT/wp1_revalidation_${E3}.csv" || \
        log "  [${E3}] revalidação não concluiu — segue, o log tem o motivo"
    else
      log "  [${E3}] sem ligante de referência ou receptor; revalidação pulada"
    fi
  done

  # 4b. escolhe o recrutador e extrai o exit vector da pose
  run_in "$ENV_MDTOOLS" python "$REPO/scripts/wp1_select_recruiter.py" \
      --screening "$WP1_SCREENING" --prep "$WP1_PREP" \
      --e3 $WP1_E3_LIST --burial "$WP1_BURIAL" \
      --out "$PIPELINE_OUT/wp1_recruiter.json"
  mark_done 4
else
  log "\n[fase 4 pulada]"
fi

# --- lê o recrutador escolhido, a menos que você tenha forçado na config ---
RJ="$PIPELINE_OUT/wp1_recruiter.json"
if [[ -f "$RJ" ]]; then
  ler_json() { python3 -c "import json,sys; d=json.load(open('$RJ'))['escolhido']; v=d.get('$1'); print(','.join(map(str,v)) if isinstance(v,list) else (v or ''))"; }
  : "${E3_NAME:=$(ler_json e3)}"
  : "${E3_RECRUITER_SDF:=$(ler_json recruiter_sdf)}"
  : "${E3_EXIT_POINT:=$(ler_json exit_point)}"
  : "${E3_EXIT_DIRECTION:=$(ler_json exit_direction)}"
  : "${E3_RECEPTOR_PDB:=$(ler_json receptor_pdb)}"
  log "\nrecrutador: $E3_NAME | exit point [$E3_EXIT_POINT] | direção [$E3_EXIT_DIRECTION]"
fi

# ===========================================================================
# FASE 5 — WP2: linkers, geometria no exit vector, sub-complexos
# ===========================================================================
if quer 5; then
  head_ 5 "WP2 — linkers e sub-complexos recrutador-linker"
  exige E3_RECRUITER_SDF "Recrutador escolhido no WP1, átomo 0 = ponto de conjugação."
  exige E3_EXIT_POINT "Coordenadas do exit vector, formato x,y,z"
  exige E3_EXIT_DIRECTION "Direção do exit vector, formato dx,dy,dz"
  exige E3_RECEPTOR_PDB "Receptor da E3 ligase, para detectar clashes."

  run_in "$ENV_MDTOOLS" python "$REPO/scripts/wp2_build_subcomplexes.py" \
      --linkers "$LINKERS_SDF" ${LINKERS_EXTRA:+--linkers-extra "$LINKERS_EXTRA"} \
      --recruiter "$E3_RECRUITER_SDF" \
      --receptor-pdb "$E3_RECEPTOR_PDB" \
      --exit-point "$E3_EXIT_POINT" \
      --exit-direction "$E3_EXIT_DIRECTION" \
      --n-confs "$LINKER_NCONFS" --n-spins "$LINKER_NSPINS" \
      --atom-range "$LINKER_ATOM_MIN" "$LINKER_ATOM_MAX" \
      --rotb-max "$LINKER_ROTB_MAX" \
      --top-n "$N_LINKERS_WP3" \
      --outdir "$PIPELINE_OUT/wp2"
  mark_done 5
else
  log "\n[fase 5 pulada]"
fi

# ===========================================================================
# FASE 6 — WP3: montar os PROTACs e emitir os jobs do PRosettaC
# ===========================================================================
if quer 6; then
  head_ 6 "WP3 — montar PROTACs e emitir jobs"
  run_in "$ENV_MDTOOLS" python "$REPO/scripts/wp3_assemble_protacs.py" \
      --subcomplexes "$PIPELINE_OUT/wp2/subcomplexes" \
      --warhead-ranking "$PIPELINE_OUT/warhead_ranking.csv" \
      --warheads-sdf "$WARHEADS_DIR/sdf" \
      --heads-pcsk9 "$DOCKING_DIR/docking/heads_pcsk9" \
      --top-warheads "$N_WARHEADS_WP2" \
      --cos-min "$ANALISE_COS_MIN" \
      --e3-structure "$E3_RECEPTOR_PDB" \
      --pcsk9-structure "$DOCKING_DIR/receptor/$(basename "${PCSK9_PDB%.pdb}")_receptor.pdb" \
      --e3-head "$E3_RECRUITER_SDF" \
      --prosettac-dir "$PROSETTAC_DIR" \
      --anchor-serial "$WARHEAD_ANCHOR_SERIAL" \
      --outdir "$PIPELINE_OUT/wp3"
  mark_done 6
else
  log "\n[fase 6 pulada]"
fi

# ===========================================================================
log "\n$(printf '=%.0s' {1..70})"
log "PIPELINE CONCLUÍDO  $(date -Is)"
log "$(printf '=%.0s' {1..70})"
log "\nSaídas em $PIPELINE_OUT"
log "\nO que exige julgamento humano e NÃO foi automatizado:"
log "  - disparar os jobs do PRosettaC: rode UM antes do lote e confira o"
log "    'Anchor atoms' (0-based vs 1-based varia por build)"
log "  - submeter os JSONs do AlphaFold 3 em alphafoldserver.com"
log "  - inspecionar as poses no ChimeraX antes de comprometer dias de MD"
log "$LOG"
