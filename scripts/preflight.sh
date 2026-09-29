#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Confere, ANTES de começar, tudo que o pipeline inteiro vai precisar.
#
#   PIPELINE_CONF=config/pipeline_vhl.conf bash scripts/preflight.sh
#
# Por que este arquivo existe
# ---------------------------
# O track da CRBN parou quatro vezes por coisa que já dava para saber no minuto
# zero: a cadeia do receptor que o config pedia e o PDB não tinha, o
# `Full: False` herdado sem razão, o head com a contagem de hidrogênio presa, o
# alcance do linker incompatível com o vão. Cada uma custou horas — e nenhuma
# precisava de horas para ser descoberta.
#
# Este script não conserta nada e não roda nada pesado: ele LISTA o estado de
# cada dependência de cada fase, em ordem, e termina com a contagem de
# bloqueios. Passa em segundos.
#
# Convenção da saída:
#   [ok]      existe e serve
#   [aviso]   ausente, mas a fase que o produz vem antes — será criado
#   [BLOQUEIO] ausente e ninguém o produz: o pipeline vai parar nisso
# ---------------------------------------------------------------------------
set -uo pipefail

AQUI="$(cd "$(dirname "$0")" && pwd)"
CONF="${PIPELINE_CONF:-$AQUI/../config/pipeline.conf}"
[[ -f "$CONF" ]] || { echo "config não encontrada: $CONF"; exit 1; }
# shellcheck disable=SC1090
source "$CONF"

bloqueios=0
avisos=0
fase_atual=""

secao() { fase_atual="$1"; echo; echo "--- $1"; }
ok()    { printf '  [ok]       %s\n' "$1"; }
aviso() { printf '  [aviso]    %s\n' "$1"; avisos=$((avisos+1)); }
bloq()  { printf '  [BLOQUEIO] %s\n' "$1"; bloqueios=$((bloqueios+1)); }

# tem_arquivo <caminho> <descrição> <quem_produz>
# `quem_produz` vazio significa: é entrada externa, ninguém a produz aqui.
tem_arquivo() {
  local p="$1" desc="$2" produz="${3:-}"
  if [[ -s "$p" ]]; then ok "$desc"
  elif [[ -n "$produz" ]]; then aviso "$desc — ausente; a $produz cria"
  else bloq "$desc — ausente, e nada no pipeline o cria: $p"; fi
}

tem_dir() {
  local p="$1" desc="$2" produz="${3:-}"
  if [[ -d "$p" ]] && [[ -n "$(ls -A "$p" 2>/dev/null)" ]]; then ok "$desc"
  elif [[ -n "$produz" ]]; then aviso "$desc — vazio/ausente; a $produz preenche"
  else bloq "$desc — vazio ou ausente: $p"; fi
}

tem_env() {
  local env="$1"
  if conda env list 2>/dev/null | awk '{print $1}' | grep -qx "$env"; then
    ok "env conda '$env'"
  else bloq "env conda '$env' não existe (conda env list)"; fi
}

tem_exe() {
  local exe="$1" desc="$2"
  if [[ -x "$exe" ]]; then ok "$desc"
  else bloq "$desc não é executável: $exe"; fi
}

echo "=============================================================="
echo "Preflight — $CONF"
echo "=============================================================="
echo "  saída do pipeline: ${PIPELINE_OUT:-<não definida>}"
echo "  E3 do track:       ${WP1_E3_LIST:-<não definida>}"

secao "ferramentas"
for e in "${ENV_MDTOOLS:-mdtools}" "${ENV_PFVS:-pf_vs}"; do tem_env "$e"; done
tem_exe "${GMX_EXE:-/usr/local/gromacs/bin/gmx}" "GROMACS (gmx)"
tem_exe "${OBABEL_EXE:-/home/soberano/miniconda3/envs/obabel_env/bin/obabel}" \
        "Open Babel"
tem_exe "${PROSETTAC_DIR:-/mnt/hd2tb/Documentos/PRosettaC}/run_prosettac.sh" \
        "PRosettaC (run_prosettac.sh)"
# O PatchDock é deduzido do protLib nos params, mas na fase 4b precisa existir
PD_DIR="$(dirname "${PROSETTAC_DIR:-/mnt/hd2tb/Documentos/PRosettaC}")/PatchDock"
if [[ -d /mnt/hd2tb/Documentos/PatchDock ]]; then PD_DIR=/mnt/hd2tb/Documentos/PatchDock; fi
if ls "$PD_DIR"/patch_dock* >/dev/null 2>&1; then ok "PatchDock em $PD_DIR"
else bloq "PatchDock não encontrado (procurei em $PD_DIR) — a fase 4b depende dele"; fi

secao "fases 1-3 — warheads e ancoragem na PCSK9 (independentes da E3)"
tem_dir "${WARHEADS_DIR:-}" "warheads enumeradas ($WARHEADS_DIR)" "fase 1"
tem_arquivo "${PCSK9_PDB:-}" "cristal da PCSK9 ($PCSK9_PDB)"
tem_arquivo "${DOCKING_DIR:-}/docking/warhead_ranking.csv" \
            "ranking das warheads" "fase 2/3"
tem_dir "${DOCKING_DIR:-}/docking/heads_pcsk9" \
        "poses das warheads na PCSK9" "fase 2"
tem_arquivo "${DOCKING_DIR:-}/receptor/6U26_receptor.pdb" \
            "receptor da PCSK9 preparado" "fase 2"

secao "fase 4 — WP1: recrutador da E3"
for E3 in ${WP1_E3_LIST:-}; do
  D=$(find "${WP1_PREP:-/dev/null}" -maxdepth 1 -type d -name "${E3}*" 2>/dev/null | head -1)
  if [[ -z "$D" ]]; then
    bloq "preparo da $E3 não existe em ${WP1_PREP:-} — sem ele não há"
    printf '             recrutador, e a fase 4 não tem o que escolher.\n'
    printf '             É a etapa manual do WP1: baixar o PDB, limpar em\n'
    printf '             ChimeraX, extrair o ligante de referência, docar.\n'
    continue
  fi
  ok "preparo da $E3 em $(basename "$D")"
  REC=$(find "$D" -name "*_receptor.pdb" 2>/dev/null | head -1)
  REF=$(find "$D" -name "*ref_ligand*" 2>/dev/null | head -1)
  [[ -n "$REC" ]] && ok "  receptor da $E3: $(basename "$REC")" \
                  || bloq "  receptor da $E3 (*_receptor.pdb) ausente em $D"
  [[ -n "$REF" ]] && ok "  ligante de referência: $(basename "$REF")" \
                  || bloq "  ligante de referência da $E3 ausente — sem ele não"
  # a cadeia declarada existe mesmo no arquivo? foi isto que quebrou 104 configs
  if [[ -n "$REC" ]]; then
    CADS=$(awk '/^ATOM/{print substr($0,22,1)}' "$REC" | sort -u | tr -d '\n ')
    ok "  cadeias presentes no receptor: '$CADS'"
  fi
  SCR=$(find "${WP1_SCREENING:-/dev/null}" -maxdepth 2 -type d -name "${E3}*" 2>/dev/null | head -1)
  [[ -n "$SCR" ]] && ok "  triagem da $E3 em $(basename "$SCR")" \
                  || bloq "  triagem (docking) da $E3 ausente em ${WP1_SCREENING:-}"
done

secao "fase 4b — vão exigido pelo par E3/alvo (o portão novo)"
if [[ -s "${PIPELINE_OUT:-}/span_requirement.json" ]]; then
  ok "já medido: $(python3 -c "
import json;d=json.load(open('${PIPELINE_OUT}/span_requirement.json'))
print(f\"vão útil {d['vao_util_A']} Å -> alcance exigido {d['alcance_exigido_A']} Å ({d['ligacoes_minimas']} ligações)\")" 2>/dev/null)"
else
  aviso "ainda não medido; a fase 4b mede (precisa de PatchDock e de uma"
  printf '             preparação do PRosettaC, ~20 min)\n'
fi

secao "fase 5 — WP2: linkers"
tem_arquivo "${LINKERS_SDF:-}" "biblioteca de linkers Chemspace"
[[ -n "${ANCHORS_SDF:-}" ]] && tem_arquivo "${ANCHORS_SDF}" "âncoras Chemspace"
if [[ -s "${LINKERS_SDF:-}" ]]; then
  N=$(grep -c '^\$\$\$\$' "$LINKERS_SDF" 2>/dev/null || echo 0)
  ok "  $N linkers no arquivo"
fi

secao "fase 6-7 — montagem e ranking"
tem_arquivo "${PIPELINE_OUT:-}/wp2/subcomplexes_manifest.json" \
            "sub-complexos do WP2" "fase 5"

secao "fase 8-9 — MD e ternário"
tem_exe "${GMX_EXE:-/usr/local/gromacs/bin/gmx}" "GROMACS para a MD"
if [[ "${PROSETTAC_FULL:-False}" == "True" ]]; then
  ok "PROSETTAC_FULL=True (amostragem completa)"
else
  bloq "PROSETTAC_FULL='${PROSETTAC_FULL:-}' — com False são 500 soluções"
  printf '             globais em vez de 1000 e 10 refinamentos em vez de 50.\n'
  printf '             Concluir ausência de geometria com isso é concluir errado.\n'
fi

echo
echo "=============================================================="
if [[ $bloqueios -gt 0 ]]; then
  echo "  $bloqueios BLOQUEIO(S) e $avisos aviso(s)."
  echo "  Resolva os bloqueios antes de lançar: eles não somem sozinhos, e"
  echo "  cada um deles pararia o pipeline no meio."
  exit 1
fi
echo "  nenhum bloqueio; $avisos aviso(s) — são fases que ainda vão produzir"
echo "  o que falta, na ordem certa."
echo
echo "  Para lançar, destacado do terminal:"
echo "    setsid env PIPELINE_CONF=$CONF bash $AQUI/run_pipeline.sh \\"
echo "        --from 4 > ~/pipeline_vhl.log 2>&1 & disown -a"
