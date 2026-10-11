#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Achar o diretório de job de um candidato. UMA função, UM lugar.
#
#   source scripts/achar_job.sh
#   DIR="$(achar_job_dir "$PIPELINE_OUT" "$CAND")" || { ... }
#   achar_job_cfg "$PIPELINE_OUT"            # todos os configs vivos
#
# Por que isto existe
# -------------------
# A mesma busca estava copiada em cinco scripts:
#
#   find "$OUT" -maxdepth 5 -type f -name prosetta_config.txt \
#        -path "*/$CAND/*" | sort        # e depois head -1
#
# Eu arquivei o material do recrutador antigo em
# $PIPELINE_OUT/_recrutador_antigo_<data>/, dentro da própria saída, e nenhuma
# dessas cinco buscas exclui nada. Pior que azar: `_` ordena antes de `w`, então
# `sort | head -1` escolhe o arquivo morto SEMPRE, e não às vezes. A fase 6a
# mediu o vão sobre os Init0/Init1 da talidomida e anunciou "preparação já
# existe — nada a fazer", exatamente o que o arquivamento devia ter evitado.
#
# Duas defesas, porque uma só já falhou:
#   1. o arquivamento passou a guardar FORA de $PIPELINE_OUT
#   2. esta busca ignora os diretórios de arquivo morto — inclusive os que já
#      estão no disco de quem arquivou antes deste conserto
#
# E uma função em vez de seis cópias, porque a sexta cópia escrita amanhã
# repetiria o defeito. Este projeto já pagou três vezes pelo mesmo `find`
# duplicado (o status.sh nomeando o candidato errado, o preflight bloqueando
# pelo recrutador velho, e agora isto).
# ---------------------------------------------------------------------------

PROSETTAC_CFG_NOME="${PROSETTAC_CFG_NOME:-prosetta_config.txt}"

# Diretórios que estão DENTRO da saída mas não são material vivo. Em array
# para entrar no find como argumentos separados — numa string única o find
# receberia tudo como um só termo e casaria nada.
ARQUIVO_MORTO_FIND=(
  -not -path "*/_recrutador_antigo_*"
  -not -path "*/_arquivo_morto_*"
)

achar_job_cfg() {  # <OUT> [CAND] -> um caminho por linha, ordenado
  local out="$1" cand="${2:-}"
  if [[ -n "$cand" ]]; then
    find "$out" -maxdepth 5 -type f -name "$PROSETTAC_CFG_NOME" \
         -path "*/$cand/*" "${ARQUIVO_MORTO_FIND[@]}" 2>/dev/null | sort
  else
    find "$out" -maxdepth 5 -type f -name "$PROSETTAC_CFG_NOME" \
         "${ARQUIVO_MORTO_FIND[@]}" 2>/dev/null | sort
  fi
}

achar_job_dir() {  # <OUT> <CAND> -> o diretório, ou nada e status 1
  local primeiro
  primeiro="$(achar_job_cfg "$1" "$2" | head -1)"
  [[ -n "$primeiro" ]] || return 1
  dirname "$primeiro"
}
