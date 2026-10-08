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
# PIPELINE_CONF manda. Quatro scripts deste repositório fixavam
# `config/pipeline.conf` — a config da CRBN — e ignoravam a variável, então nem
# passá-la ajudava: o comando rodava no track abandonado devolvendo números com
# cara de certo.
CONF="${PIPELINE_CONF:-$AQUI/../config/pipeline.conf}"
[[ -f "$CONF" ]] || [[ "$CONF" = /* ]] || CONF="$AQUI/../$CONF"
[[ -f "$CONF" ]] || { echo "config não encontrada: $CONF"; exit 1; }
# shellcheck disable=SC1090
source "$CONF"
# A CRBN foi abandonada; a guarda recusa rodar nela sem PROTAC_TRACK_OK=1.
# shellcheck disable=SC1090
source "$(dirname "$CONF")/guarda_de_track.sh"


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

# O detalhe de quem produz o quê é do check_wp1.py, que importa as funções do
# wp1_select_recruiter e por isso responde a pergunta certa: não "existe uma
# pasta VHL", mas "o seletor encontraria o que precisa dentro dela".
echo
echo "  detalhe do WP1 (poses, colunas do CSV, cadeias do receptor):"
echo "    conda run -n ${ENV_MDTOOLS:-mdtools} python $AQUI/check_wp1.py \\"
echo "        --screening ${WP1_SCREENING:-} --prep ${WP1_PREP:-} \\"
echo "        --e3 ${WP1_E3_LIST:-} ${ANCHORS_SDF:+--anchors-sdf ${ANCHORS_SDF}}"

# O preflight conferia ARQUIVOS e não VALORES de config — e foi por aí que o
# track da VHL parou na fase 6: o emissor exige a cadeia quando o PDB tem mais
# de uma, o config declarava PCSK9_CHAIN="B", e o driver não passava o valor
# adiante. Arquivo presente, config correto, fiação faltando. Conferir a cadeia
# CONTRA o arquivo pega as três coisas de uma vez.
secao "cadeias declaradas x cadeias que os arquivos têm"
confere_cadeia() {   # confere_cadeia <pdb> <valor_do_config> <nome_da_variavel>
  local pdb="$1" valor="$2" nome="$3"
  if [[ ! -s "$pdb" ]]; then
    aviso "$nome: $pdb ainda não existe — a fase que o cria vem antes"
    return
  fi
  local cads
  cads=$(awk '/^ATOM/{print substr($0,22,1)}' "$pdb" | sort -u | tr -d '\n ')
  if [[ -z "$valor" ]]; then
    if [[ ${#cads} -eq 1 ]]; then
      ok "$nome vazia, e $(basename "$pdb") tem só '$cads' — detectada sem risco"
    else
      bloq "$nome está VAZIA e $(basename "$pdb") tem as cadeias '$cads'."
      printf '             Nada aqui pode escolher por você: preencha %s.\n' "$nome"
    fi
  elif [[ "$cads" == *"$valor"* ]]; then
    ok "$nome='$valor' existe em $(basename "$pdb") (cadeias: '$cads')"
  else
    bloq "$nome='$valor' NÃO existe em $(basename "$pdb") — ele tem '$cads'"
  fi
}

PCSK9_REC="${DOCKING_DIR:-}/receptor/$(basename "${PCSK9_PDB:-x.pdb}" .pdb)_receptor.pdb"
confere_cadeia "$PCSK9_REC" "${PCSK9_CHAIN:-}" "PCSK9_CHAIN"
for E3 in ${WP1_E3_LIST:-}; do
  D=$(find "${WP1_PREP:-/dev/null}" -maxdepth 1 -type d -name "${E3}*" 2>/dev/null | head -1)
  [[ -n "$D" ]] || continue
  R=$(find "$D" -name "*_receptor.pdb" 2>/dev/null | head -1)
  [[ -n "$R" ]] && confere_cadeia "$R" "${E3_CHAIN:-}" "E3_CHAIN ($E3)"
done

secao "fase 6a — vão exigido pelo par E3/alvo (o portão novo)"
if [[ -s "${PIPELINE_OUT:-}/span_requirement.json" ]]; then
  ok "já medido: $(python3 -c "
import json;d=json.load(open('${PIPELINE_OUT}/span_requirement.json'))
print(f\"vão útil {d['vao_util_A']} Å -> alcance exigido {d['alcance_exigido_A']} Å ({d['ligacoes_minimas']} ligações)\")" 2>/dev/null)"
else
  aviso "ainda não medido; a fase 6a mede (precisa de PatchDock e de uma"
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

secao "fase 10 — MD do ternário (nível iii)"
# Estas conferências existem porque cada uma delas falha TARDE. O acpype é
# chamado depois de o PRosettaC ter rodado dias; o MDAnalysis só no fim da
# produção; o disco enche no meio da terceira réplica. Um `import` que custa
# um segundo agora vale horas depois.
tem_exe "${ACPYPE_EXE:-$HOME/miniconda3/envs/mdtools/bin/acpype}" \
        "acpype (GAFF2/AM1-BCC do PROTAC)"

if conda run -n "${ENV_MDTOOLS:-mdtools}" python -c \
     "import MDAnalysis; from MDAnalysis.lib.distances import capped_distance" \
     >/dev/null 2>&1; then
  ok "MDAnalysis com capped_distance (contagem de contatos por células)"
else
  bloq "MDAnalysis ausente ou sem capped_distance no env '${ENV_MDTOOLS:-mdtools}'"
  printf '             A análise do nível (iii) conta contatos por lista de\n'
  printf '             células; sem isso são 17 milhões de distâncias por quadro.\n'
fi

# md_analyze_ternario.py IMPORTA do md_analyze.py (carregar, descobrir_resname).
# Um erro de import ali só apareceria depois de dias de produção.
if conda run -n "${ENV_MDTOOLS:-mdtools}" python -c \
     "import sys; sys.path.insert(0, '$AQUI'); import md_analyze_ternario" \
     >/dev/null 2>&1; then
  ok "md_analyze_ternario importa (inclusive o md_analyze de onde ele reusa)"
else
  bloq "md_analyze_ternario.py não importa no env '${ENV_MDTOOLS:-mdtools}'"
  printf '             Rode para ver o motivo:\n'
  printf '               conda run -n %s python %s/md_analyze_ternario.py --help\n' \
         "${ENV_MDTOOLS:-mdtools}" "$AQUI"
fi

# Espaço em disco. Três réplicas de 200 ns de um sistema de ~120 mil átomos,
# com quadro a cada 200 ps, dão ~1 GB de .xtc cada; somados aos .trr, .edr,
# checkpoints e ao sistema solvatado, a ordem é 15 GB. Encher o disco na
# terceira réplica perde as três, porque o mdrun morre sem fechar o .xtc.
if [[ -n "${PIPELINE_OUT:-}" ]]; then
  livre_gb=$(df -BG --output=avail "$(dirname "${PIPELINE_OUT}")" 2>/dev/null \
               | tail -1 | tr -dc '0-9')
  if [[ -z "${livre_gb:-}" ]]; then
    aviso "não consegui medir o espaço livre em $(dirname "$PIPELINE_OUT")"
  elif [[ "$livre_gb" -lt 15 ]]; then
    bloq "só ${livre_gb} GB livres; a MD (iii) pede ~15 GB"
    printf '             O mdrun morre sem fechar o .xtc, então encher o disco\n'
    printf '             na terceira réplica perde as três.\n'
  else
    ok "espaço livre: ${livre_gb} GB (a MD (iii) pede ~15 GB)"
  fi
fi

if [[ "${MD_TERNARIO_AUTO:-0}" == "1" ]]; then
  aviso "MD_TERNARIO_AUTO=1 — a fase 10 vai TOMAR A GPU por dias"
  printf '             São 3 x %s ns com `-update gpu`. A workstation é\n' \
         "${MD_NS_PROD:-200}"
  printf '             compartilhada: isto precisa estar combinado com quem\n'
  printf '             mais a usa, e não é o preflight que combina.\n'
else
  ok "MD_TERNARIO_AUTO=0 — a fase 10 imprime os comandos e espera você"
fi

secao "WP1 — o recrutador recruta a E3 que dizemos?"
# Esta seção existe porque a ausência dela custou sete dias de PRosettaC. O WP1
# escolhe o recrutador por SCORE DE DOCKING e enterramento, e nenhuma das duas
# medidas sabe que E3 a molécula de fato recruta. O recrutador escolhido para o
# "track da VHL" era um análogo de talidomida — ligante canônico de CRBN.
RECR=$(find "${WP1_PREP:-}" "${PIPELINE_OUT:-}" -maxdepth 4 \
            -name "recruiter_*_H.sdf" -o -name "recruiter_*.sdf" 2>/dev/null \
         | grep -v "_H.sdf" | sort | head -1)
if [[ -z "${RECR:-}" ]]; then
  aviso "não achei o SDF do recrutador para checar o quimiotipo"
elif conda run -n "${ENV_MDTOOLS:-mdtools}" python "$AQUI/checar_recrutador.py" \
       --e3 "${WP1_E3_LIST:-VHL}" --sdf "$RECR" >/tmp/.recr.$$ 2>&1; then
  ok "recrutador $(basename "$RECR") tem assinatura de ${WP1_E3_LIST:-VHL}"
  rm -f /tmp/.recr.$$
else
  cod=$?
  if [[ $cod -eq 2 ]]; then
    bloq "o recrutador recruta OUTRA E3, não a ${WP1_E3_LIST:-VHL}"
  else
    aviso "o recrutador não tem assinatura canônica de ${WP1_E3_LIST:-VHL}"
  fi
  sed 's/^/             /' /tmp/.recr.$$ | head -14
  rm -f /tmp/.recr.$$
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
