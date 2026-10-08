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
# RELANCE O DRIVER DEPOIS DE UM `git pull` QUE MEXA NESTE ARQUIVO.
# O bash lê o script do disco conforme executa, guardando a posição em bytes.
# Editar o arquivo embaixo de um processo em andamento faz ele retomar na
# posição antiga de um conteúdo novo — e o que ele executa a partir dali é
# lixo. Um laço já parseado (como a espera da fase 9) termina inteiro, então
# há tempo: `pkill -f run_pipeline.sh` e relançar é seguro, porque a guarda de
# execução em andamento do run_prosettac.sh impede relançamento duplicado.
#
# Cada fase grava um marcador em $PIPELINE_OUT/.done_<n>. Fases já concluídas
# são puladas, então relançar depois de uma queda custa só a fase interrompida.
# Os scripts pesados são retomáveis por conta própria, então mesmo uma fase
# interrompida no meio retoma de onde parou.
# ---------------------------------------------------------------------------
set -uo pipefail

# PIPELINE_CONF manda. Quatro scripts deste repositório fixavam
# `config/pipeline.conf` — a config da CRBN — e ignoravam a variável, então nem
# passá-la ajudava: o comando rodava no track abandonado devolvendo números com
# cara de certo.
CONF="${PIPELINE_CONF:-$(dirname "$0")/../config/pipeline.conf}"
[[ -f "$CONF" ]] || [[ "$CONF" = /* ]] || CONF="$(dirname "$0")/../$CONF"
[[ -f "$CONF" ]] || { echo "config não encontrada: $CONF"; exit 1; }
# shellcheck disable=SC1090
source "$CONF"
# EXPORTADA para os filhos. O driver chama run_prosettac.sh, prosettac_prepare.sh
# e patchdock_span_scan.sh, e cada um deles carrega a config por conta própria —
# sem exportar, eles cairiam no default (a CRBN) mesmo com o driver na VHL.
export PIPELINE_CONF="$CONF"
# A CRBN foi abandonada; a guarda recusa rodar nela sem PROTAC_TRACK_OK=1.
# shellcheck disable=SC1090
source "$(dirname "$CONF")/guarda_de_track.sh"

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

# O candidato do ranking e o diretório do job dele. As fases 9 e 10 precisam
# dos MESMOS dois valores, e este projeto já pagou por ter duas buscas do mesmo
# diretório: elas divergiram, o driver vigiou a pasta errada, e declararia
# falha em 168 h sobre uma execução que deu certo. Uma função, um lugar.
#
# A busca é pelo prosetta_config.txt e não pelo nome do diretório, porque
# `-type d -name "$CAND"` casa tanto wp3/protacs/<cand> (onde só vive o
# protac.smi) quanto wp3/prosettac/<cand> — e o `head -1` pegava o primeiro.
achar_candidato() {  # define CAND e DIRC, ou sai != 0
  local rank="${1:-1}"
  local rankcsv="$PIPELINE_OUT/protac_ranking.csv"
  [[ -s "$rankcsv" ]] || { log "*** $rankcsv não existe: rode a fase 7"; return 1; }
  CAND=$(python3 -c "
import csv
r=list(csv.DictReader(open('$rankcsv')))
i=min(max(int('$rank')-1,0),len(r)-1) if r else 0
print(r[i]['candidate_id'] if r else '')")
  [[ -n "$CAND" ]] || { log "*** ranking vazio"; return 1; }
  DIRC=$(dirname "$(find "$PIPELINE_OUT" -maxdepth 5 -type f \
                         -name "prosetta_config.txt" -path "*/$CAND/*" \
                         2>/dev/null | sort | head -1)" 2>/dev/null)
  [[ "$DIRC" == "." ]] && DIRC=""
  if [[ -z "$DIRC" ]]; then
    log "*** não achei o diretório do job de $CAND — a fase 6 emite os configs"
    return 1
  fi
  return 0
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
FASES=(1 2 3 4 5 6 6a 6b 7 8 9 10)

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
      ${E3_CHAIN:+--e3-chain "$E3_CHAIN"} \
      ${PCSK9_CHAIN:+--pcsk9-chain "$PCSK9_CHAIN"} \
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

  achar_candidato "${MD_RANK:-1}" || exit 1
  log "  candidato: $CAND"
  log "  job em:    $DIRC"

  if [[ -d "$DIRC/Results" ]]; then
    log "  já há Results — pulando direto para a concordância"
  elif [[ $DRY -eq 0 ]]; then
    # O run_prosettac.sh valida, conserta head e âncora, e lança DESTACADO —
    # então aqui é preciso esperar o resultado, e não o lançamento.
    PIPELINE_OUT="$PIPELINE_OUT" PROSETTAC_DIR="$PROSETTAC_DIR" \
      bash "$REPO/scripts/run_prosettac.sh" "$CAND" 2>&1 | tee -a "$LOG"
    # PIPESTATUS e não $?: com o `| tee` o status que chega é o do tee, e um
    # lançamento que falhou passaria por bem-sucedido.
    cod_lanc=${PIPESTATUS[0]}
    if [[ $cod_lanc -eq 7 ]]; then
      # 7 = já havia execução em andamento. Não é falha: é o motivo de existir
      # a guarda. Seguir para a espera é o certo — relançar embaralharia as duas.
      log "  (já estava rodando; sigo para a espera sem relançar)"
    elif [[ $cod_lanc -ne 0 ]]; then
      log "*** o PRosettaC não subiu; o log acima diz por quê"
      exit 1
    fi

    # A ESPERA vive em scripts/esperar_results.sh, e vive lá por um motivo:
    # ela já errou duas vezes, e das duas o erro só apareceria depois de horas
    # ou dias — um contador amarrado ao nome de arquivo de UMA etapa, e um
    # limite de RELÓGIO que declararia falha em 08/10 de uma execução medida
    # para terminar em 10/10. Um laço que só se exercita em 12 h não é um laço
    # testado. Como script próprio, os cinco desfechos dele se verificam em
    # segundos com TERNARIO_ESPERA_S=1.
    #
    #   0 = Results apareceu
    #   1 = terminou SEM Results  -> conclusão sobre o CANDIDATO
    #   2 = travamento            -> conclusão sobre a MÁQUINA
    # Os valores vão EXPLÍCITOS na linha. `source` da config define variáveis
    # de shell, não de ambiente: sem isto o esperar_results.sh não vê nada do
    # que a config declara e cai nos defaults DELE — a config diria 12 h e o
    # script usaria 12 h por coincidência, e no dia em que os dois divergissem
    # ninguém saberia qual valeu. Foi o teste que pegou isto, e é a mesma falha
    # de fiação que já custou duas corridas a este projeto: arquivo presente,
    # config correta, ninguém passando o valor.
    TERNARIO_SEM_PROGRESSO_H="${TERNARIO_SEM_PROGRESSO_H:-12}" \
    TERNARIO_LIMITE_H="${TERNARIO_LIMITE_H:-336}" \
    TERNARIO_ESPERA_S="${TERNARIO_ESPERA_S:-300}" \
    TERNARIO_BATIDAS="${TERNARIO_BATIDAS:-6}" \
      bash "$REPO/scripts/esperar_results.sh" "$DIRC" 2>&1 | tee -a "$LOG"
    cod_espera=${PIPESTATUS[0]}
    pronto=0; travado=0
    case "$cod_espera" in
      0) pronto=1;;
      2) travado=1;;
    esac
    if [[ $travado -eq 1 ]]; then
      # Travado NÃO é "sem geometria". Dizer a segunda coisa quando aconteceu a
      # primeira manda o projeto refazer a fase 5 por causa de uma fila presa.
      log "\n  TRAVAMENTO, não ausência de resultado. Últimas linhas do log.txt:"
      tail -8 "$DIRC/log.txt" 2>/dev/null | sed 's/^/      /' | tee -a "$LOG"
      log ""
      log "*** A execução parou de avançar: fila, arquivos e log.txt sem mudar"
      log "*** por ${TERNARIO_SEM_PROGRESSO_H:-12} h. Isto NÃO é 'não existe"
      log "*** geometria ternária' —"
      log "*** é a execução presa, e a conclusão sobre o candidato continua"
      log "*** em aberto."
      log "***"
      log "*** O que olhar, nesta ordem:"
      log "***   squeue -u ${USER:-$(id -un)} | head   jobs presos em PENDING?"
      log "***   sinfo -o '%P %a %l %D %t %N'      a partição está drenada?"
      log "***   df -h $(dirname "$PIPELINE_OUT")  o disco encheu?"
      log "***   tail -40 $DIRC/log.txt"
      log "***"
      log "***   bash $REPO/scripts/esperar_results.sh $DIRC   (só reatacha)"
      log "***"
      log "*** Resolvido o motivo, retome com --from 9: a guarda de execução"
      log "*** em andamento impede relançamento duplicado."
      exit 9
    fi
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

  # O VEREDITO, corrigido. A versão anterior exigia 5 REPRESENTANTES de
  # cluster mutuamente a menos de 5 Å — e isso era quase logicamente
  # impossível, porque o clustering.py do PRosettaC roda DBSCAN com eps = 4 Å
  # sobre a MESMA medida: representantes de clusters distintos estão a mais de
  # 4 Å por construção. O portão pedia que o DBSCAN tivesse produzido clusters
  # que ele mesmo teria fundido, e o "não convergiu" que saía dali era artefato
  # do critério. Isso se demonstra lendo os dois códigos, sem olhar resultado.
  #
  # A régua agora é a do PRÓPRIO PRosettaC, que imprime no result_summary.txt
  # "Out of them N have at least 5 members": 5 membros é o que os autores da
  # ferramenta escolheram como cluster com massa. Adotá-la não é mover o alvo —
  # ela existia antes deste resultado e não foi calibrada por ele.
  MEMBROS=$(python3 -c "
from pathlib import Path
res = Path('$DIRC/Results')
tam = [len(list(d.rglob('*.pdb'))) for d in res.iterdir() if d.is_dir()]
print(max(tam) if tam else 0)" 2>/dev/null || echo 0)
  log "\n  maior cluster: ${MEMBROS:-?} membros (régua do PRosettaC: 5+)"
  if [[ "${MEMBROS:-0}" -lt "${TERNARIO_MEMBROS_MIN:-5}" ]]; then
    log ""
    log "*** NENHUMA pose dominante: o maior cluster tem ${MEMBROS:-0} membros,"
    log "*** abaixo da régua de ${TERNARIO_MEMBROS_MIN:-5} da própria ferramenta."
    log "*** Isto é resultado, não falha: não há o que levar ao método"
    log "*** ortogonal nem à MD nível (iii)."
    log "*** Próximo candidato: MD_RANK=$(( ${MD_RANK:-1} + 1 )) em $CONF, --from 9."
    exit 5
  fi
  log "  há pose dominante a testar — mas ela é CANDIDATA, não confirmada."
  log "  Quem decide é o método ORTOGONAL (Boltz-2/AlphaFold 3), como o WP3 pede."
  mark_done 9
else
  log "\n[fase 9 pulada]"
fi

# ===========================================================================
# FASE 10 — MD do nível (iii): o complexo ternário completo (DIAS)
# ===========================================================================
# É a pergunta final do WP3, e ela só existe quando a fase 9 aprovou: a
# interface ternária que o PRosettaC modelou sobrevive ao solvente, ou ela
# existia só porque o docking a construiu?
#
# Esta fase NÃO roda sozinha por default (MD_TERNARIO_AUTO=0). A workstation é
# compartilhada, a produção são 3 x 200 ns, e o md_run.sh toma a GPU com
# `-update gpu` por dias seguidos. Tomar a placa de outras pessoas sem aviso
# não é uma decisão do driver — é uma combinação entre quem usa a máquina.
# Para ligar: MD_TERNARIO_AUTO=1 na config, ou rode os três passos à mão.
if quer 10; then
  if [[ "${MD_TERNARIO_AUTO:-0}" != "1" ]]; then
    log ""
    head_ 10 "MD nível (iii) — ternário completo (NÃO automática)"
    log "  A fase 9 aprovou a pose. O nível (iii) está pronto para rodar, e"
    log "  pede a GPU por dias — então ele espera a sua decisão."
    log ""
    log "  Os três comandos, em ordem (cada um destacável):"
    if achar_candidato "${MD_RANK:-1}"; then
      log "    MDT=$PIPELINE_OUT/md_ternario"
      log ""
      log "    # 1) preparar (segundos) — separa proteínas e ligante, confere"
      log "    #    o SMILES contra as coordenadas e imprime o custo do sistema"
      log "    conda run -n $ENV_MDTOOLS python $REPO/scripts/md_prepare_ternario.py \\"
      log "        --candidato $CAND --prosettac-dir $(dirname "$DIRC") \\"
      log "        --outdir \$MDT --ram-limite-gb ${MD_RAM_LIMITE_GB:-8}"
      log ""
      log "    # 2) rodar (dias; GPU dedicada) — destacado, sobrevive ao SSH"
      log "    setsid bash $REPO/scripts/md_run.sh \$MDT ${MD_N_REPLICAS:-3} \\"
      log "        > ~/md_ternario.log 2>&1 & disown -a"
      log ""
      log "    # 3) analisar (minutos) — PBC primeiro, senão o RMSD mede a caixa"
      log "    bash $REPO/scripts/md_fix_pbc.sh \$MDT"
      log "    conda run -n $ENV_MDTOOLS python $REPO/scripts/md_analyze_ternario.py \\"
      log "        --md-dir \$MDT"
    fi
    log ""
    log "  Para o driver fazer isso sozinho: MD_TERNARIO_AUTO=1 em $CONF"
  else
    head_ 10 "MD nível (iii) — ternário completo (dias, GPU dedicada)"
    achar_candidato "${MD_RANK:-1}" || exit 1
    log "  candidato: $CAND"

    # A pose que sai da fase 9 é CANDIDATA. Levar uma pose não confirmada a
    # 3 x 200 ns é gastar dias de GPU simulando uma hipótese que um segundo
    # método poderia derrubar em horas — e a metodologia do WP3 pede os dois:
    # "PRosettaC and AlphaFold 3 as orthogonal, complementary tools".
    ORTO="$PIPELINE_OUT/ortogonal_confirmado.json"
    if [[ ! -s "$ORTO" && "${TERNARIO_ACEITAR_SEM_ORTOGONAL:-0}" != "1" ]]; then
      log ""
      log "*** Falta a confirmação ORTOGONAL da pose: $ORTO"
      log "*** A fase 9 entrega um cluster dominante, não uma pose confirmada."
      log "*** Rode o Boltz-2 (ou o AF3) sobre este ternário e compare com o"
      log "*** representante do maior cluster antes de gastar a GPU por dias."
      log "***"
      log "*** Para seguir sem isso, de propósito e declarando na tese:"
      log "***   TERNARIO_ACEITAR_SEM_ORTOGONAL=1 em $CONF"
      exit 10
    fi
    [[ -d "$DIRC/Results" ]] || {
      log "*** $DIRC/Results não existe: a fase 9 não terminou"; exit 1; }

    MDT="$PIPELINE_OUT/md_ternario"

    # 10a. preparar. Daqui sai o md_sistema.json com as cinco chaves que o
    #      md_run.sh lê pelo nome, e o portão de custo do sistema solvatado.
    run_in "$ENV_MDTOOLS" python "$REPO/scripts/md_prepare_ternario.py" \
        --candidato "$CAND" --prosettac-dir "$(dirname "$DIRC")" \
        --outdir "$MDT" --ram-limite-gb "${MD_RAM_LIMITE_GB:-8}"

    # 10b. rodar. Mesmo script do nível (ii): a metodologia do WP3 é a mesma
    #      (ff14SB + GAFF2/AM1-BCC, TIP3P, caixa de 1,3 nm, corte 0,9 nm com
    #      PME, EM -> 5 ns NPT -> 3 x 200 ns a 310 K/1 atm, quadro a cada
    #      200 ps). O que muda é o conteúdo da caixa, não os parâmetros.
    log "\n[\$ md_run.sh $MDT ${MD_N_REPLICAS:-3}]"
    if [[ $DRY -eq 0 ]]; then
      NS_PROD="${MD_NS_PROD:-200}" NS_NPT="${MD_NS_NPT:-5}" \
        MD_TRUNCAR_PERTO="${MD_TRUNCAR_PERTO:-}" \
        bash "$REPO/scripts/md_run.sh" "$MDT" "${MD_N_REPLICAS:-3}" \
        2>&1 | tee -a "$LOG" || {
          log "*** a MD do ternário falhou; veja $LOG"; exit 1; }

      # PBC antes da análise, sempre. Os cortes nas lacunas do cristal fizeram
      # da proteína várias moléculas, e o GROMACS envolve cada uma por conta —
      # segmentos da mesma cadeia aparecem em lados opostos da caixa e o RMSD
      # passa a medir a aresta. Com DUAS proteínas o artefato é pior: E3 e alvo
      # podem ser envolvidos separadamente, e a interface "desaparece".
      bash "$REPO/scripts/md_fix_pbc.sh" "$MDT" 2>&1 | tee -a "$LOG" \
        || log "  [ATENÇÃO] a correção de PBC falhou; a análise vai avisar"
    fi

    # 10c. analisar. Sai 6 quando REPROVA, e reprovar não é falha do pipeline:
    #      é a resposta de que esta pose não sobrevive ao solvente.
    log "\n[\$ md_analyze_ternario.py --md-dir $MDT]  (env $ENV_MDTOOLS)"
    if [[ $DRY -eq 0 ]]; then
      conda run --no-capture-output -n "$ENV_MDTOOLS" python -u \
          "$REPO/scripts/md_analyze_ternario.py" --md-dir "$MDT" \
          2>&1 | tee -a "$LOG"
      cod=${PIPESTATUS[0]}
      if [[ $cod -eq 6 ]]; then
        log ""
        log "*** O complexo ternário REPROVOU nos critérios do nível (iii)."
        log "*** Isto é resultado, não falha: a interface modelada não se"
        log "*** sustentou em solvente, e essa é uma conclusão sobre o"
        log "*** candidato — para a tese, com os números no CSV."
        log "*** Próximo candidato: MD_RANK=$(( ${MD_RANK:-1} + 1 )) em $CONF,"
        log "*** depois --from 9."
        exit 8
      elif [[ $cod -ne 0 ]]; then
        log "*** a análise do ternário falhou (código $cod); veja $LOG"
        exit 1
      fi
    fi
    mark_done 10
  fi
else
  log "\n[fase 10 pulada]"
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
log "  - COMBINAR o horário da MD (iii) com quem mais usa a workstation: são"
log "    3 x 200 ns com a GPU dedicada, e por isso a fase 10 não roda sozinha"
log "  - decidir trocar de E3 ou de sítio quando a fase 9 reprovar: o número"
log "    está medido, a decisão de projeto é sua"
log "$LOG"
