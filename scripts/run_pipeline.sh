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
# Para deixar rodando e desligar o notebook, SEM instalar nada:
#   setsid bash ~/Protac_env/scripts/run_pipeline.sh > ~/pipeline.log 2>&1 &
#   disown -a
#   tail -f ~/pipeline.log       # Ctrl+C sai do tail, não do pipeline
#
# `setsid` cria uma sessão nova: o processo deixa de ser filho do terminal e
# sobrevive ao fim da conexão SSH e ao fechamento do VS Code.
#
# Com tmux instalado (snap install tmux), a alternativa é:
#   tmux new -s protac
#   bash ~/Protac_env/scripts/run_pipeline.sh 2>&1 | tee ~/pipeline.log
#   (Ctrl+B, depois D)  /  tmux attach -t protac para voltar
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
mark_done() { [[ $DRY -eq 1 ]] && return 0; date -Is > "$(done_marker "$1")"; }

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

# A ORDEM das fases vive aqui, e é ela que o --from respeita. É uma LISTA e não
# uma contagem porque 4b e 6b entraram depois, entre fases que já tinham
# marcador no disco: renumerar invalidaria os `.done_N` de quem já rodou, e um
# pipeline que se esquece do que fez refaz 46 h de MD.
#
# 6a e 6b são os portões que este projeto aprendeu a ter, e a POSIÇÃO deles é o
# ponto: medir o vão que o par E3/alvo exige, e reprovar quem não alcança, ANTES
# da MD. Com a CRBN os dois vieram depois, e o candidato que consumiu 46 h de MD
# era geometricamente impossível desde o começo.
#
# A 6a fica depois da montagem porque quem produz o Init0/Init1 que a medição
# usa é a preparação do PRosettaC, e ela precisa de um candidato. Não dá para
# medir antes de existir um PROTAC — mas dá, e é o que importa, para medir antes
# de gastar MD.
FASES=(1 2 3 4 5 6 6a 6b 7 8 9)

indice_da_fase() {
  local alvo="$1" i=0
  for f in "${FASES[@]}"; do
    [[ "$f" == "$alvo" ]] && { echo "$i"; return 0; }
    i=$((i+1))
  done
  echo "-1"
}

if [[ "${LIST:-0}" == "1" ]]; then
  echo "Fases concluídas:"
  for n in "${FASES[@]}"; do
    if is_done "$n"; then echo "  [x] fase $n  ($(cat "$(done_marker "$n")"))"
    else echo "  [ ] fase $n"; fi
  done
  exit 0
fi

quer() {  # quer <n> -> deve rodar esta fase?
  local n="$1"
  [[ -n "$ONLY" ]] && { [[ "$ONLY" == "$n" ]]; return; }
  local i_n i_from
  i_n=$(indice_da_fase "$n")
  i_from=$(indice_da_fase "$FROM")
  if [[ "$i_from" == "-1" ]]; then
    echo "*** --from $FROM não é uma fase. As fases são: ${FASES[*]}" >&2
    exit 1
  fi
  (( i_n >= i_from )) && ! is_done "$n"
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

  if [[ ! -f "$PCSK9_PDB" ]]; then
    if [[ $DRY -eq 1 ]]; then
      log "[dry-run] baixaria $PCSK9_PDB do RCSB"
    else
      log "baixando o cristal da PCSK9"
      mkdir -p "$STRUCTURES"
      wget -q -O "$PCSK9_PDB" \
        "https://files.rcsb.org/download/$(basename "${PCSK9_PDB%.pdb}").pdb" \
        || { log "*** falha ao baixar $PCSK9_PDB"; exit 1; }
    fi
  fi

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
      ${ANCHORS_SDF:+--anchors-sdf "$ANCHORS_SDF"} \
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
      "--exit-point=$E3_EXIT_POINT" \
      "--exit-direction=$E3_EXIT_DIRECTION" \
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

  # A fase 3 pode ter sido rodada à mão antes do pipeline existir, gravando o
  # ranking no diretório do docking em vez de $PIPELINE_OUT. Procura nos dois,
  # e se não achar em nenhum, refaz a análise (custa segundos).
  RANKING="$PIPELINE_OUT/warhead_ranking.csv"
  if [[ ! -f "$RANKING" ]]; then
    ALT="$DOCKING_DIR/docking/warhead_ranking.csv"
    if [[ -f "$ALT" ]]; then
      log "  ranking encontrado em $ALT"
      RANKING="$ALT"
    else
      log "  ranking ausente: refazendo a análise da fase 3"
      run_in "$ENV_MDTOOLS" python "$REPO/scripts/analyze_warhead_docking.py" \
          --docking "$DOCKING_DIR/docking" \
          --site "$DOCKING_DIR/pcsk9_site.json" \
          --warheads-csv "$WARHEADS_DIR/pcsk9_warheads.csv" \
          --sd-max "$ANALISE_SD_MAX" --out "$RANKING"
    fi
  fi
  run_in "$ENV_MDTOOLS" python "$REPO/scripts/wp3_assemble_protacs.py" \
      --subcomplexes "$PIPELINE_OUT/wp2/subcomplexes" \
      --warhead-ranking "$RANKING" \
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
# FASE 6a — quanto o par E3/alvo EXIGE de alcance (minutos)
# ===========================================================================
# Este é o portão que faltava. Ele mede, com o PatchDock sobre as duas poses
# validadas, a distância mínima entre os átomos de conjugação para existir
# alguma colocação das duas proteínas — e a distância a partir da qual existem
# soluções SUFICIENTES. Com a CRBN: 14 Å davam UMA transformada, 18 Å davam 40.
if quer 6a; then
  head_ 6a "Vão exigido pelo par E3/alvo — PatchDock (~30 min)"

  CANDS="$PIPELINE_OUT/wp3/protac_candidates.csv"
  [[ -s "$CANDS" ]] || { log "*** $CANDS não existe: rode a fase 6"; exit 1; }
  # A sonda é um candidato qualquer: o que se mede aqui é geometria das duas
  # PROTEÍNAS, e o PROTAC da sonda só define a restrição de distância — que a
  # varredura sobrescreve ponto a ponto.
  SONDA=$(python3 -c "
import csv,sys
r=list(csv.DictReader(open('$CANDS')))
print(r[0]['candidate_id'] if r else '')")
  [[ -n "$SONDA" ]] || { log "*** nenhum candidato em $CANDS"; exit 1; }
  log "  sonda: $SONDA  (mede as proteínas, não o PROTAC dela)"

  if [[ $DRY -eq 0 ]]; then
    PIPELINE_OUT="$PIPELINE_OUT" PROSETTAC_DIR="$PROSETTAC_DIR" \
      bash "$REPO/scripts/prosettac_prepare.sh" "$SONDA" 2>&1 | tee -a "$LOG" \
      || { log "*** a preparação do PRosettaC não produziu Init0/Init1."; exit 1; }
    PIPELINE_OUT="$PIPELINE_OUT" \
      bash "$REPO/scripts/patchdock_span_scan.sh" "$SONDA" \
        ${SPAN_DISTANCIAS:-10 12 14 16 18 20 22 25 28 32 100} 2>&1 | tee -a "$LOG" \
      || { log "*** a varredura do PatchDock falhou."; exit 1; }
  fi

  SCAN="$PIPELINE_OUT/wp3/span_scan_${SONDA}.csv"
  run_in "$ENV_MDTOOLS" python "$REPO/scripts/span_requirement.py" \
      --scan "$SCAN" \
      --transformadas-uteis "${SPAN_TRANSFORMADAS_UTEIS:-20}" \
      --razao "${SPAN_RAZAO_MEDIANA_TETO:-0.64}" \
      --out "$PIPELINE_OUT/span_requirement.json"
  mark_done 6a
else
  log "\n[fase 6a pulada]"
fi

# --- lê o requisito medido; daqui para baixo ele é critério, não palpite ----
REQ_JSON="$PIPELINE_OUT/span_requirement.json"
SPAN_MIN=""
if [[ -f "$REQ_JSON" ]]; then
  SPAN_MIN=$(python3 -c "import json;print(json.load(open('$REQ_JSON'))['alcance_exigido_A'])")
  SPAN_LIG=$(python3 -c "import json;print(json.load(open('$REQ_JSON'))['ligacoes_minimas'])")
  log "\nalcance exigido: $SPAN_MIN Å  (~$SPAN_LIG ligações de cadeia)"
fi

# ===========================================================================
# FASE 6b — portão de alcance: quem não vence o vão não vai para a MD
# ===========================================================================
if quer 6b; then
  head_ 6b "Portão de alcance dos PROTACs montados (segundos)"
  exige SPAN_MIN "Alcance exigido — sai da fase 6a, em span_requirement.json"

  run_in "$ENV_MDTOOLS" python "$REPO/scripts/linker_span.py" \
      --protacs "$PIPELINE_OUT/wp3/protac_candidates.csv" \
      --e3-head "$E3_RECRUITER_SDF" \
      --warheads-sdf "$DOCKING_DIR/docking/heads_pcsk9" \
      --requisito "$SPAN_MIN" \
      --minimo-absoluto "$(python3 -c "import json;print(json.load(open('$REQ_JSON'))['vao_util_A'])")" \
      --out "$PIPELINE_OUT/wp3/linker_span.csv"

  # O portão reprova, e reprovar aqui é o serviço dele. Zero aprovados NÃO é
  # falha do pipeline: é a resposta de que o pool de linkers é curto para este
  # par de sítios — e agora com um número para refiltrar, em vez de intuição.
  N_OK=$(python3 -c "
import csv
print(sum(1 for r in csv.DictReader(open('$PIPELINE_OUT/wp3/linker_span.csv'))
          if r.get('situacao') == 'alcança'))" 2>/dev/null || echo 0)
  log "  candidatos que alcançam $SPAN_MIN Å: $N_OK"
  if [[ "${N_OK:-0}" -eq 0 && $DRY -eq 0 ]]; then
    log ""
    log "*** NENHUM candidato alcança o vão exigido. O pipeline para AQUI, e"
    log "*** isto é o portão funcionando: com a CRBN, um candidato nesta mesma"
    log "*** situação consumiu 46 h de MD antes de o PatchDock devolver zero."
    log "***"
    log "*** A decisão é da fase 5: refiltre o catálogo Chemspace exigindo"
    log "*** cadeia de ~$SPAN_LIG ligações ou mais entre os pontos de"
    log "*** conjugação, e relance com --from 5."
    log "***   LINKER_ATOM_MIN=$SPAN_LIG em $CONF"
    exit 3
  fi
  mark_done 6b
else
  log "\n[fase 6b pulada]"
fi

# ===========================================================================
# FASE 7 — ranquear os PROTACs e escolher o candidato da MD
# ===========================================================================
if quer 7; then
  head_ 7 "Ranquear os PROTACs (segundos)"
  # O --span ELIMINA quem não alcança, em vez de dar nota baixa: um PROTAC que
  # não vence o vão não é pior, é impossível. Sem o --span o rank_protacs avisa,
  # alto, que não está verificando isso.
  run_in "$ENV_MDTOOLS" python "$REPO/scripts/rank_protacs.py" \
      --protacs "$PIPELINE_OUT/wp3/protac_candidates.csv" \
      --warhead-ranking "${RANKING:-$PIPELINE_OUT/warhead_ranking.csv}" \
      --subcomplexes "$PIPELINE_OUT/wp2/subcomplexes_manifest.json" \
      ${SPAN_MIN:+--span "$PIPELINE_OUT/wp3/linker_span.csv" --span-min "$SPAN_MIN"} \
      --out "$PIPELINE_OUT/protac_ranking.csv"
  mark_done 7
else
  log "\n[fase 7 pulada]"
fi

# ===========================================================================
# FASE 8 — MD do nível (ii) para o candidato escolhido
# ===========================================================================
if quer 8 && [[ "${MD_AUTO:-1}" == "1" ]]; then
  head_ 8 "MD nível (ii) — E3 + recrutador-linker-warhead (horas)"

  run_in "$ENV_MDTOOLS" python "$REPO/scripts/md_prepare.py" \
      --candidato "$PIPELINE_OUT/md_candidato.json" \
      --recruiter "$E3_RECRUITER_SDF" \
      --receptor "$E3_RECEPTOR_PDB" \
      --rank "${MD_RANK:-1}" \
      --outdir "$PIPELINE_OUT/md"

  log "\n[\$ md_run.sh $PIPELINE_OUT/md ${MD_N_REPLICAS:-3}]"
  if [[ $DRY -eq 0 ]]; then
    NS_PROD="${MD_NS_PROD:-200}" NS_NPT="${MD_NS_NPT:-5}" \
      MD_TRUNCAR_PERTO="${MD_TRUNCAR_PERTO:-}" \
      bash "$REPO/scripts/md_run.sh" "$PIPELINE_OUT/md" "${MD_N_REPLICAS:-3}" \
      2>&1 | tee -a "$LOG" || {
        log "*** a MD falhou; veja $LOG"; exit 1; }
  fi

  run_in "$ENV_MDTOOLS" python "$REPO/scripts/md_analyze.py" \
      --md-dir "$PIPELINE_OUT/md"
  mark_done 8
else
  log "\n[fase 8 pulada]"
fi

# ===========================================================================
# FASE 9 — complexo ternário: PRosettaC e o portão de concordância (horas)
# ===========================================================================
# Com a CRBN esta fase era manual, e foi onde tudo parou. Agora ela roda sozinha
# e, o que importa mais, ela tem VEREDITO: sem um grupo de modelos concordando,
# não existe pose, e não há MD nível (iii) a fazer.
if quer 9; then
  head_ 9 "Complexo ternário — PRosettaC + concordância (horas)"

  RANK="$PIPELINE_OUT/protac_ranking.csv"
  [[ -s "$RANK" ]] || { log "*** $RANK não existe: rode a fase 7"; exit 1; }
  CAND=$(python3 -c "
import csv
r=list(csv.DictReader(open('$RANK')))
i=min(max(int('${MD_RANK:-1}')-1,0),len(r)-1) if r else 0
print(r[i]['candidate_id'] if r else '')")
  [[ -n "$CAND" ]] || { log "*** ranking vazio"; exit 1; }
  log "  candidato: $CAND"

  DIRC=$(find "$PIPELINE_OUT" -maxdepth 5 -type d -name "$CAND" 2>/dev/null | head -1)
  if [[ -z "$DIRC" ]]; then
    log "*** não achei o diretório do job de $CAND — a fase 6 emite os configs"
    exit 1
  fi

  if [[ -d "$DIRC/Results" ]]; then
    log "  já há Results — pulando direto para a concordância"
  elif [[ $DRY -eq 0 ]]; then
    # O run_prosettac.sh valida, conserta head e âncora, e lança DESTACADO —
    # então aqui é preciso esperar o resultado, e não o lançamento.
    PIPELINE_OUT="$PIPELINE_OUT" PROSETTAC_DIR="$PROSETTAC_DIR" \
      bash "$REPO/scripts/run_prosettac.sh" "$CAND" 2>&1 | tee -a "$LOG" \
      || { log "*** o PRosettaC não subiu; o log acima diz por quê"; exit 1; }

    LIMITE_H="${TERNARIO_LIMITE_H:-24}"
    log "\n  esperando o Results (limite de ${LIMITE_H} h)..."
    pronto=0
    for _ in $(seq 1 $((LIMITE_H * 12))); do
      sleep 300
      [[ -d "$DIRC/Results" ]] && { pronto=1; break; }
      # morreu no meio: o log.txt do PRosettaC diz a última etapa alcançada
      if ! pgrep -u "${USER:-$(id -un)}" -f "PRosettaC.*main\.py" >/dev/null \
         && [[ -f "$DIRC/log.txt" ]] \
         && grep -q "run has finished" "$DIRC/log.txt"; then
        break
      fi
    done
    if [[ $pronto -eq 0 ]]; then
      log "\n  o PRosettaC terminou sem Results. Últimas linhas do log.txt:"
      tail -6 "$DIRC/log.txt" 2>/dev/null | sed 's/^/      /' | tee -a "$LOG"
      log ""
      log "*** Sem modelos ternários não há o que levar à MD nível (iii)."
      log "*** O resumo em result_summary.txt diz quantos modelos passaram o"
      log "*** limiar de energia — 3% significa ponte tensa (volte à fase 5 com"
      log "*** mais folga de alcance), e 10%+ com dispersão significa que a"
      log "*** interface não é definida por este método (troque de E3 ou sítio)."
      exit 4
    fi
  fi

  [[ -s "$DIRC/result_summary.txt" ]] && sed 's/^/  /' "$DIRC/result_summary.txt" \
    | tee -a "$LOG"

  run_in "$ENV_MDTOOLS" python "$REPO/scripts/prosettac_agreement.py" \
      --candidato "$CAND" --work "$(dirname "$DIRC")" \
      --corte "${TERNARIO_CORTE_A:-5.0}"

  # O veredito: um grupo de modelos concordando é pose; dois modelos não são.
  MAIOR=$(python3 -c "
import csv, itertools, math
import numpy as np
f='$DIRC/concordancia_clusters.csv'
rows=list(csv.reader(open(f)))
nomes=rows[0][1:]
M=np.full((len(nomes),len(nomes)), np.nan)
for i,r in enumerate(rows[1:]):
    for j,v in enumerate(r[1:]):
        if v: M[i,j]=float(v)
A=(M<=${TERNARIO_CORTE_A:-5.0}) & ~np.eye(len(nomes),dtype=bool) & ~np.isnan(M)
melhor=[]
def exp(at,ca):
    global melhor
    if len(at)>len(melhor): melhor=list(at)
    for k,v in enumerate(ca):
        if len(at)+len(ca)-k<=len(melhor): return
        exp(at+[v],[u for u in ca[k+1:] if A[v,u]])
exp([],list(range(len(nomes))))
print(len(melhor))" 2>/dev/null || echo 0)
  log "\n  maior grupo concordando a ${TERNARIO_CORTE_A:-5.0} Å: ${MAIOR:-?} modelos"
  if [[ "${MAIOR:-0}" -lt "${TERNARIO_GRUPO_MIN:-5}" && $DRY -eq 0 ]]; then
    log ""
    log "*** Os modelos NÃO convergem numa interface (${MAIOR:-0} < ${TERNARIO_GRUPO_MIN:-5})."
    log "*** Isto é resultado, não falha: nenhuma pose ternária testável saiu"
    log "*** daqui, e levar uma pose arbitrária à MD nível (iii) seria simular"
    log "*** uma hipótese escolhida por acaso."
    log "*** Próximo candidato: MD_RANK=$(( ${MD_RANK:-1} + 1 )) em $CONF, --from 9."
    exit 5
  fi
  mark_done 9
else
  log "\n[fase 9 pulada]"
fi

# ===========================================================================
log "\n$(printf '=%.0s' {1..70})"
log "PIPELINE CONCLUÍDO  $(date -Is)"
log "$(printf '=%.0s' {1..70})"
log "\nSaídas em $PIPELINE_OUT"
log "\nO que exige julgamento humano e NÃO foi automatizado:"
log "  - submeter os JSONs do AlphaFold 3 em alphafoldserver.com (o Boltz-2"
log "    local cobre a via ortogonal; o AF3 é confirmação independente)"
log "  - inspecionar no ChimeraX a pose ternária aprovada, antes da MD (iii)"
log "  - decidir trocar de E3 ou de sítio quando a fase 9 reprovar: o número"
log "    está medido, a decisão de projeto é sua"
log "$LOG"
