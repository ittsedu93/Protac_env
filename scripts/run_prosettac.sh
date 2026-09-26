#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Dispara o PRosettaC para um candidato.
#
#   bash scripts/run_prosettac.sh <candidate_id>
#
# A ORDEM aqui é a parte que importa, e cada passo depende do anterior:
#
#   1. achar o diretório do candidato   (procurando o config, não adivinhando)
#   2. ENTRAR nele                      (o PRosettaC só trabalha com nomes
#                                        relativos; tudo depois disto assume
#                                        que o diretório atual é o de trabalho)
#   3. trazer as entradas para dentro   (copiar + reescrever o config)
#   4. limpar restos de tentativa morta (senão o clean_pdb pula a limpeza)
#   5. validar                          (arquivos, cadeias, âncoras)
#   6. lançar
#   7. conferir que subiu               (três sinais, não um)
#
# Roda UM job e para. A convenção de `Anchor atoms` (0-based vs 1-based) varia
# por build, e um lote inteiro com o átomo errado é um lote perdido.
# ---------------------------------------------------------------------------
set -uo pipefail

AQUI="$(cd "$(dirname "$0")" && pwd)"
# O que veio do ambiente ganha do arquivo de configuração: quem passa a
# variável na linha de comando está dizendo explicitamente o que quer, e é
# também o que torna este script testável fora da máquina do laboratório.
_OUT_ENV="${PIPELINE_OUT:-}"
_PROS_ENV="${PROSETTAC_DIR:-}"
CONF="$AQUI/../config/pipeline.conf"
[[ -f "$CONF" ]] && source "$CONF"
OUT="${_OUT_ENV:-${PIPELINE_OUT:-$HOME/PRosettaC_runs/vhl_crbn_pcsk9_protac/pipeline}}"
PROSETTAC="${_PROS_ENV:-${PROSETTAC_DIR:-/mnt/hd2tb/Documentos/PRosettaC}}"
CFG_NOME="prosetta_config.txt"

CAND="${1:?uso: run_prosettac.sh <candidate_id>}"

echo "=============================================================="
echo "PRosettaC — $CAND"
echo "=============================================================="

# --- 1. achar o diretório do candidato ------------------------------------
DIR=""
while IFS= read -r achado; do DIR="$(dirname "$achado")"; break
done < <(find "$OUT" -maxdepth 5 -type f -name "$CFG_NOME" -path "*/$CAND/*" \
              2>/dev/null | sort)
if [[ -z "$DIR" ]]; then
  echo "não achei $CFG_NOME para $CAND em $OUT"
  echo
  echo "Candidatos com job emitido:"
  find "$OUT" -maxdepth 5 -type f -name "$CFG_NOME" 2>/dev/null \
    | sed "s|/$CFG_NOME$||" | xargs -r -n1 basename | sort -u \
    | sed 's/^/    /' | head -20
  exit 1
fi

# --- 2. entrar nele: daqui para baixo, tudo é relativo --------------------
cd "$DIR" || { echo "não consegui entrar em $DIR"; exit 1; }
echo "  diretório: $DIR"

[[ -x "$PROSETTAC/run_prosettac.sh" ]] || {
  echo "não achei $PROSETTAC/run_prosettac.sh (executável)"; exit 1; }

# --- 3. trazer as entradas para dentro ------------------------------------
python "$AQUI/prosettac_localize.py" --dir . --config "$CFG_NOME" || exit 1

# --- 4. limpar restos de uma tentativa que morreu no meio -----------------
# O clean_pdb vê um <struct>_<cadeia>.pdb existente, pula a limpeza, e o que
# segue opera sobre um arquivo pela metade.
rm -f ./*_[A-Z].fasta ./*_[A-Z].pdb ./log.txt

echo
echo "config:"
sed 's/^/    /' "$CFG_NOME"
echo

# --- 5. validar ------------------------------------------------------------
falhou=0

ESTRUTURAS=$(awk -F': ' '/^Structures:/{print $2}' "$CFG_NOME")
CADEIAS=$(awk -F': ' '/^Chains:/{print $2}' "$CFG_NOME")
HEADS=$(awk -F': ' '/^Heads:/{print $2}' "$CFG_NOME")
PROTAC=$(awk -F': ' '/^Protac:/{print $2}' "$CFG_NOME")
ANCHORS=$(awk -F': ' '/^Anchor atoms:/{print $2}' "$CFG_NOME")

for f in $ESTRUTURAS $HEADS $PROTAC; do
  [[ -s "$f" ]] || { echo "  [FALTA] $f"; falhou=1; }
done

i=1
for est in $ESTRUTURAS; do
  cad=$(echo "$CADEIAS" | cut -d' ' -f$i)
  if [[ -s "$est" ]]; then
    presentes=$(awk '/^ATOM/{print substr($0,22,1)}' "$est" | sort -u | tr -d '\n ')
    n=$(awk -v c="$cad" '/^ATOM/ && substr($0,22,1)==c{print substr($0,23,5)}' \
          "$est" | sort -u | wc -l)
    if [[ "$n" -eq 0 ]]; then
      echo "  [CADEIA ERRADA] $est: pediu '$cad', tem '$presentes'"
      falhou=1
    else
      echo "  cadeia $cad de $est: $n resíduos"
    fi
  fi
  i=$((i+1))
done

# os âncoras têm de existir dentro de cada head
i=1
for head in $HEADS; do
  anc=$(echo "$ANCHORS" | cut -d' ' -f$i)
  if [[ -s "$head" ]]; then
    n_at=$(sed -n '4p' "$head" | awk '{print $1+0}')
    if [[ -n "$anc" && "$n_at" -gt 0 && "$anc" -gt "$n_at" ]]; then
      echo "  [ÂNCORA FORA] $head tem $n_at átomos, âncora pedida: $anc"
      falhou=1
    else
      echo "  âncora $anc em $head ($n_at átomos)"
    fi
  fi
  i=$((i+1))
done

if [[ $falhou -ne 0 ]]; then
  echo
  echo "*** corrija o $CFG_NOME antes de lançar:"
  echo "***   nano $DIR/$CFG_NOME"
  exit 1
fi

echo
echo "  Anchor atoms = $ANCHORS"
echo "  >>> CONFIRA no resultado deste job que os âncoras são os átomos que"
echo "  >>> ligam cada head ao linker. 0-based vs 1-based varia por build."
echo

# --- 6. já há resultado? ---------------------------------------------------
if [[ -d Results || -d results ]]; then
  echo "  já há resultado aqui — nada a fazer"
  exit 0
fi

# --- 7. lançar -------------------------------------------------------------
# O run_prosettac.sh do PRosettaC roda com `set -euo pipefail` e faz
# `conda activate`. Isso dispara o hook de DESATIVAÇÃO do env ativo, e o
# mpivars.deactivate.sh do mdtools referencia SETVARS_CALL sem defini-la —
# sob `set -eu` isso mata o script antes da primeira linha de trabalho.
export SETVARS_CALL="${SETVARS_CALL:-}"

echo "[$(date -Is)] lançando..."
setsid "$PROSETTAC/run_prosettac.sh" "$DIR" "$CFG_NOME" \
    > run_prosettac.log 2>&1 &
disown -a
sleep 8

# --- 8. conferir que subiu -------------------------------------------------
vivo=0
pgrep -f "$CFG_NOME" > /dev/null && vivo=1
squeue -h -u "${USER:-$(id -un)}" 2>/dev/null | grep -q . && vivo=1
[[ -d Patchdock_Results ]] && vivo=1

if [[ $vivo -eq 1 ]]; then
  echo "  RODANDO — pode desligar o notebook"
  squeue -h -u "${USER:-$(id -un)}" 2>/dev/null | head -5 | sed 's/^/      /'
  echo "  log: tail -f $DIR/run_prosettac.log"
else
  echo "  [NÃO SUBIU] log completo:"
  sed 's/^/      /' run_prosettac.log 2>/dev/null | head -40
  exit 1
fi
