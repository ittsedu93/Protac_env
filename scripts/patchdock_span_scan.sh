#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Quanto o linker PRECISA alcançar para existir um ternário?
#
#   bash scripts/patchdock_span_scan.sh <candidate_id> [D1 D2 ...]
#
# O PRosettaC parou aqui:
#
#   distanceConstraints 1 3 10
#   Transforms remain after interface clustering with RMSD 2 : 0
#   INFO: PatchDock did not find any global docking solution within the
#         geometrical constraints
#
# Exigindo que os dois átomos de conjugação — um dentro do bolsão da CRBN, o
# outro dentro do sítio da PCSK9 — fiquem a no máximo 10 Å, nenhuma colocação
# das duas proteínas sobreviveu. Isso tem duas leituras possíveis, e elas levam
# a decisões opostas:
#
#   (a) o par CRBN/PCSK9 não admite ternário nenhum        -> trocar de E3/sítio
#   (b) admite, mas exige um vão maior que 10 Å            -> trocar de LINKER
#
# Distinguir é relaxar a restrição e ver a partir de quanto aparecem soluções.
# O menor valor que produz transformadas é o ALCANCE MÍNIMO que o par exige, e
# isso é um requisito de projeto — entra no WP2 como filtro, com número.
#
# Este script NÃO refaz a preparação: reaproveita o Init0.pdb, o Init1.pdb e o
# Patchdock_params.txt que o PRosettaC já produziu, e roda só o PatchDock, uma
# vez por distância. Minutos por ponto, não horas.
#
# A distância 100 Å é o CONTROLE: sem restrição efetiva, o docking tem de
# devolver milhares de transformadas. Devolvendo zero também aí, o problema não
# é o linker — é a preparação, e o resto da tabela não significa nada.
# ---------------------------------------------------------------------------
set -uo pipefail

AQUI="$(cd "$(dirname "$0")" && pwd)"
_OUT_ENV="${PIPELINE_OUT:-}"
CONF="$AQUI/../config/pipeline.conf"
[[ -f "$CONF" ]] && source "$CONF"
OUT="${_OUT_ENV:-${PIPELINE_OUT:-$HOME/PRosettaC_runs/vhl_crbn_pcsk9_protac/pipeline}}"
CFG_NOME="prosetta_config.txt"

CAND="${1:?uso: patchdock_span_scan.sh <candidate_id> [distâncias...]}"
shift
DISTANCIAS=("$@")
[[ ${#DISTANCIAS[@]} -eq 0 ]] && DISTANCIAS=(10 12 14 16 18 20 25 30 100)

echo "=============================================================="
echo "Alcance mínimo do linker — $CAND"
echo "=============================================================="

# --- achar o diretório, do mesmo jeito que o run_prosettac.sh --------------
DIR=""
while IFS= read -r achado; do DIR="$(dirname "$achado")"; break
done < <(find "$OUT" -maxdepth 5 -type f -name "$CFG_NOME" -path "*/$CAND/*" \
              2>/dev/null | sort)
[[ -n "$DIR" ]] || { echo "não achei $CFG_NOME para $CAND em $OUT"; exit 1; }
cd "$DIR" || exit 1
echo "  diretório: $DIR"

for f in Init0.pdb Init1.pdb Patchdock_params.txt; do
  [[ -s "$f" ]] || {
    echo "  [FALTA] $f — este script reaproveita a preparação do PRosettaC."
    echo "          Rode primeiro:  bash scripts/run_prosettac.sh $CAND"
    exit 1; }
done

# --- o binário do PatchDock, deduzido do próprio params -------------------
# O params aponta o chem.lib; o binário está ao lado dele. Deduzir é melhor que
# pedir o caminho: o arquivo que a ferramenta usa é a fonte mais confiável.
PD_DIR="$(dirname "$(awk '/^protLib/{print $2; exit}' Patchdock_params.txt)")"
PD_BIN=""
for c in patch_dock.Linux patch_dock patchdock PatchDock; do
  [[ -x "$PD_DIR/$c" ]] && { PD_BIN="$PD_DIR/$c"; break; }
done
if [[ -z "$PD_BIN" ]]; then
  echo "  não achei o executável do PatchDock em $PD_DIR"
  echo "  o que há lá:"
  ls "$PD_DIR" 2>/dev/null | head -20 | sed 's/^/      /'
  exit 1
fi
echo "  PatchDock: $PD_BIN"

RESTRICAO="$(awk '/^distanceConstraints/{print; exit}' Patchdock_params.txt)"
A_REC="$(echo "$RESTRICAO" | awk '{print $2}')"
A_LIG="$(echo "$RESTRICAO" | awk '{print $3}')"
D_ORIG="$(echo "$RESTRICAO" | awk '{print $4}')"
echo "  restrição original: átomo $A_REC do receptor, $A_LIG do ligante,"
echo "                      no máximo $D_ORIG Å  -> 0 transformadas"
echo

# --- a varredura ----------------------------------------------------------
# Em subdiretório próprio: mexer no Patchdock_params.txt do PRosettaC faria a
# próxima execução dele rodar com uma restrição que não é a dele.
ESCANEIO=span_scan
mkdir -p "$ESCANEIO"
ln -sf ../Init0.pdb ../Init1.pdb "$ESCANEIO/"
cd "$ESCANEIO" || exit 1

CSV="$OUT/wp3/span_scan_${CAND}.csv"     # FORA do diretório de trabalho: a
mkdir -p "$(dirname "$CSV")"             # limpeza do run_prosettac.sh o apagaria
echo "dist_thr_A,transforms,solucoes_no_res" > "$CSV"

printf "  %8s  %13s  %10s\n" "máx (Å)" "transformadas" "soluções"
printf "  %8s  %13s  %10s\n" "--------" "-------------" "----------"

primeiro_ok=""
for D in "${DISTANCIAS[@]}"; do
  P="params_${D}.txt"
  sed "s|^distanceConstraints .*|distanceConstraints $A_REC $A_LIG $D|" \
      ../Patchdock_params.txt > "$P"
  # o params fala de Init0.pdb/Init1.pdb por nome relativo, e os links estão aqui
  "$PD_BIN" "$P" "out_${D}.res" > "pd_${D}.log" 2>&1

  # o número que interessa é o que o PRosettaC leu: o que sobra depois do
  # agrupamento de interface. `tail -1` porque o PatchDock imprime por patch.
  T="$(grep "Transforms remain after interface clustering" "pd_${D}.log" \
        | tail -1 | awk -F': ' '{print $2+0}')"
  [[ -z "$T" ]] && T="?"
  # `grep -c` já imprime 0 quando não casa nada, e AINDA sai com status 1 — um
  # `|| echo 0` aqui fazia a variável virar "0\n0", que suja a tabela e quebra
  # uma linha do CSV em duas.
  S="$(grep -cE "^ *[0-9]+ \|" "out_${D}.res" 2>/dev/null)" || true
  [[ -n "$S" ]] || S=0

  printf "  %8s  %13s  %10s\n" "$D" "$T" "$S"
  echo "$D,$T,$S" >> "$CSV"
  [[ -z "$primeiro_ok" && "$T" != "?" && "$T" -gt 0 ]] && primeiro_ok="$D"
done

echo
echo "  dados: $CSV"
echo

# --- ler a tabela ---------------------------------------------------------
CONTROLE="$(awk -F, '$1==100{print $2}' "$CSV")"
if [[ -n "$CONTROLE" && "$CONTROLE" == "0" ]]; then
  echo "  [ATENÇÃO] o CONTROLE de 100 Å também deu zero. Sem restrição efetiva"
  echo "  o docking tinha de devolver milhares de transformadas, então o"
  echo "  problema NÃO é o alcance do linker — é a preparação (Init0/Init1) ou"
  echo "  os parâmetros de superfície. O resto da tabela não significa nada."
  exit 1
fi

if [[ -z "$primeiro_ok" ]]; then
  echo "  Nenhuma distância produziu transformadas, mas o controle produziu."
  echo "  Leia isso como: a restrição de distância, e não o docking, é o que"
  echo "  elimina tudo — e nem relaxá-la até $( echo "${DISTANCIAS[*]}" | awk '{print $(NF-1)}' ) Å resolve."
  exit 0
fi

echo "  ALCANCE MÍNIMO: $primeiro_ok Å entre os dois átomos de conjugação."
echo
echo "  O linker deste candidato alcança até 11 Å (o initial_distances.hist)."
if [[ "$primeiro_ok" -le 11 ]]; then
  echo "  Isso CABE no alcance dele — então o zero de $D_ORIG Å não é falta de"
  echo "  comprimento, é a margem: a restrição do PRosettaC usa o topo da"
  echo "  distribuição, e sobra pouco. Vale relançar com Full: True e um"
  echo "  histograma mais amostrado."
else
  echo "  Isso NÃO cabe: faltam $(( primeiro_ok - 11 )) Å. É achado, não falha —"
  echo "  o par CRBN/PCSK9 exige um linker mais longo que o deste candidato, e"
  echo "  agora isso é um número que filtra o catálogo no WP2 em vez de uma"
  echo "  suspeita. Some ~1,3 Å por átomo de cadeia para estimar quantos"
  echo "  átomos faltam."
fi
