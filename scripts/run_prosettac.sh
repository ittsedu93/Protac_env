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
#   4. limpar restos de tentativa morta (o PRosettaC pula etapa cujo arquivo
#                                        já existe, mesmo com zero byte)
#   5. validar E CONSERTAR              (arquivos, cadeias, heads, âncoras —
#                                        head e âncora são calculáveis, então
#                                        são corrigidos aqui, não reportados
#                                        para o usuário digitar outro comando)
#   6. lançar
#   7. conferir que subiu               (cinco sinais e 60 s de paciência; e
#                                        não subindo, DESCOBRIR por quê aqui)
#
# Roda UM job e para: um lote inteiro com o átomo errado é um lote perdido.
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
# ANTES de limpar: um resultado pronto não é resto.
if [[ -d Results || -d results ]]; then
  echo "  já há resultado aqui — nada a fazer"
  exit 0
fi

# O PRosettaC decide o que fazer pela EXISTÊNCIA do arquivo: o clean_pdb vê um
# <struct>_<cadeia>.pdb e pula a limpeza, e um `random_sampling.sdf` de zero
# byte faz pular a amostragem que o produziria. A primeira versão disto apagava
# só os produtos que eu tinha visto naquele dia — e nomear produtos um a um é
# uma lista que sempre está incompleta. A regra agora é a inversa: fica o que o
# config declara como entrada, sai todo o resto.
python "$AQUI/prosettac_clean.py" --dir . --config "$CFG_NOME" --aplicar \
  || exit 1

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

# Os heads precisam ser legíveis pelo RDKit do jeito que o PRosettaC os lê:
# ele protona o head (`<head>_H.sdf`) e remapeia a âncora casando o original
# contra a versão protonada. Um .sdf que fixa 0 hidrogênio nos carbonos casa
# zero átomos, o translate_anchors devolve -1, e o -1 vira OverflowError
# cinquenta linhas adiante. Corrigir é reescrever o .sdf — coordenadas
# intactas, com .bak — então isto roda com --aplicar: conferir e não consertar
# só devolveria o erro ao usuário para ele digitar o comando seguinte.
if ! python "$AQUI/prosettac_fix_heads.py" --dir . --config "$CFG_NOME" \
       --aplicar > .heads.out 2>&1; then
  sed 's/^/  /' .heads.out
  rm -f .heads.out
  echo "*** os heads não descrevem moléculas válidas — veja acima."
  exit 1
fi
grep -E "OK —|CORRIGIDO|casa [0-9]+ átomos|movido para" .heads.out \
  | sed 's/^/  /'
rm -f .heads.out

# A âncora é CALCULÁVEL a partir do head e do PROTAC: é o átomo do head que
# se liga ao linker. Sendo calculável, ela é aplicada e não sugerida — quando
# está errada o PRosettaC não diz "âncora errada", o translate_anchors devolve
# -1 e o erro aparece trinta linhas adiante como índice negativo.
if python "$AQUI/prosettac_anchors.py" --dir . --config "$CFG_NOME" \
     --aplicar > .ancoras.out 2>&1; then
  grep -E "âncora =|linha correta|config atualizado" .ancoras.out \
    | sed 's/^/  /'
  ANCHORS=$(awk -F': ' '/^Anchor atoms:/{print $2}' "$CFG_NOME")
else
  sed 's/^/  /' .ancoras.out
  rm -f .ancoras.out
  echo "*** não consegui determinar as âncoras — veja acima."
  exit 1
fi
rm -f .ancoras.out

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
echo "  Anchor atoms = $ANCHORS   (calculado a partir do head e do PROTAC,"
echo "                             conferido contra a versão protonada)"
echo

# --- 6. lançar -------------------------------------------------------------
# O run_prosettac.sh do PRosettaC roda com `set -euo pipefail` e faz
# `conda activate`. Isso dispara o hook de DESATIVAÇÃO do env ativo, e o
# mpivars.deactivate.sh do mdtools referencia SETVARS_CALL sem defini-la —
# sob `set -eu` isso mata o script antes da primeira linha de trabalho.
export SETVARS_CALL="${SETVARS_CALL:-}"

LOG=run_prosettac.log
echo "[$(date -Is)] lançando..."
setsid "$PROSETTAC/run_prosettac.sh" "$DIR" "$CFG_NOME" > "$LOG" 2>&1 &
disown -a

# --- 8. conferir que subiu -------------------------------------------------
# Cinco sinais, e até 60 s de paciência. Oito segundos não bastavam: o
# `conda activate` mais a importação do RDKit e do Rosetta consomem isso antes
# de a ferramenta escrever a primeira linha, e um `pgrep` que não casa nesse
# instante não quer dizer que nada subiu.
usuario="${USER:-$(id -un)}"
esta_vivo() {
  pgrep -f "$CFG_NOME"    > /dev/null && return 0
  pgrep -f "$PROSETTAC"   > /dev/null && return 0
  # o main.py do PRosettaC, e não qualquer main.py do usuário — um `pgrep`
  # frouxo aqui daria "RODANDO" por causa de outro job na mesma máquina
  pgrep -u "$usuario" -f "PRosettaC.*main\.py" > /dev/null && return 0
  squeue -h -u "$usuario" 2>/dev/null | grep -q . && return 0
  [[ -d Patchdock_Results ]] && return 0
  return 1
}

vivo=0
for _ in $(seq 1 12); do
  sleep 5
  if esta_vivo; then vivo=1; break; fi
done

if [[ $vivo -eq 1 ]]; then
  echo "  RODANDO — pode desligar o notebook"
  squeue -h -u "$usuario" 2>/dev/null | head -5 | sed 's/^/      /'
  echo "  log: tail -f $DIR/$LOG"
  exit 0
fi

# --- 9. não subiu: descobrir POR QUÊ, aqui, agora -------------------------
# Dizer "[NÃO SUBIU]" e imprimir um log vazio devolve o problema sem nenhuma
# informação — foi o que aconteceu, e custou um ciclo. Quando o processo morre
# antes de escrever, a única fonte é a execução em primeiro plano com rastro.
echo "  [NÃO SUBIU]"
if [[ -s "$LOG" ]]; then
  echo "  log ($(wc -l < "$LOG") linhas), últimas 40:"
  tail -40 "$LOG" | sed 's/^/      /'
  exit 1
fi

echo "  o log está VAZIO: nada foi escrito, nem por bash nem por python."
echo
echo "  O script que estou chamando — $PROSETTAC/run_prosettac.sh — é onde a"
echo "  causa mora, então ele vai impresso junto com o rastro:"
echo "  --------------------------------------------------------------------"
sed -n '1,40p' "$PROSETTAC/run_prosettac.sh" | sed 's/^/      /'
echo "  --------------------------------------------------------------------"
echo
echo "  o processo morreu antes de escrever a primeira"
echo "  linha. Repetindo em primeiro plano, com rastro, por até 120 s — o que"
echo "  aparecer abaixo é a causa. (Conseguindo trabalhar, ele é interrompido"
echo "  no fim do tempo; os restos são limpos no próximo lançamento.)"
echo "  --------------------------------------------------------------------"
timeout --foreground 120 bash -x "$PROSETTAC/run_prosettac.sh" \
        "$DIR" "$CFG_NOME" 2>&1 | tail -60 | sed 's/^/      /'
codigo=${PIPESTATUS[0]}
echo "  --------------------------------------------------------------------"
if [[ $codigo -eq 124 ]]; then
  echo "  [ATENÇÃO] ele rodou os 120 s sem morrer. Então ele NÃO morre sozinho:"
  echo "  o que falhou foi o desligamento do terminal (setsid/disown). Relance"
  echo "  por conta própria, sem intermediário:"
  echo
  echo "    cd $DIR"
  echo "    nohup setsid $PROSETTAC/run_prosettac.sh $DIR $CFG_NOME \\"
  echo "        > $LOG 2>&1 < /dev/null &"
else
  echo "  saiu com código $codigo — a causa está nas linhas acima."
fi
exit 1
