#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Roda um comando e MATA se o disco livre cair abaixo de um piso.
#
#   bash scripts/rodar_com_teto_de_disco.sh --piso-gb 30 --onde /home -- \
#        python algum_script.py args
#
# Por que isto existe
# -------------------
# A workstation é compartilhada, e encher o disco dela não atrasa só o nosso
# trabalho: trava o de todo mundo, e um `mdrun` ou um `rosetta` que perde a
# escrita no meio costuma corromper o arquivo que estava abrindo. O operador
# deste projeto foi explícito — se houver risco de exceder o disco, o projeto
# para ali.
#
# "Vou acompanhar e mato se precisar" não é uma trava, é uma intenção: ninguém
# olha o `df` às 3 da manhã. Isto transforma o limite numa regra que a máquina
# aplica, com um piso que fica REGISTRADO no comando — então a decisão de
# quanto é demais é tomada antes, com a cabeça fria, e não durante.
#
# O que ele faz
# -------------
#   mede o livre ANTES, e recusa se já começar abaixo do piso
#   lança o comando, mede a cada 20 s
#   cruzando o piso: mata o grupo de processos e sai com 9
#   no fim: diz quanto o comando consumiu, que é o número que faltava para
#           planejar a etapa seguinte
#
# Mata o GRUPO (setsid + kill -PID negativo) e não o processo: um script que
# lançou filhos deixaria os filhos escrevendo depois de o pai morrer, e aí a
# trava não teria travado nada.
# ---------------------------------------------------------------------------
set -uo pipefail

PISO_GB=30
ONDE=""
INTERVALO=20
while [[ $# -gt 0 ]]; do
  case "$1" in
    --piso-gb)   PISO_GB="$2"; shift 2;;
    --onde)      ONDE="$2"; shift 2;;
    --intervalo) INTERVALO="$2"; shift 2;;
    --) shift; break;;
    *) echo "opção desconhecida: $1"; exit 1;;
  esac
done
[[ $# -gt 0 ]] || { echo "uso: $0 [--piso-gb N] [--onde CAMINHO] -- <comando>"; exit 1; }
ONDE="${ONDE:-$PWD}"

livre_gb() { df -BG --output=avail "$ONDE" 2>/dev/null | tail -1 | tr -dc '0-9'; }

L0=$(livre_gb)
if [[ -z "$L0" ]]; then
  echo "*** não consegui medir o espaço livre em $ONDE — não vou rodar às cegas"
  exit 1
fi
echo "=============================================================="
echo "teto de disco: piso de ${PISO_GB} GB em $ONDE"
echo "  livre agora: ${L0} GB"
echo "  comando:     $*"
echo "=============================================================="
if [[ "$L0" -le "$PISO_GB" ]]; then
  echo "*** JÁ está em ${L0} GB, no piso ou abaixo. Não vou começar."
  echo "*** Libere espaço ou baixe o piso de propósito com --piso-gb."
  exit 2
fi

setsid "$@" &
PG=$!
echo "  rodando (grupo $PG); medindo a cada ${INTERVALO}s"

matou=0
while kill -0 "$PG" 2>/dev/null; do
  sleep "$INTERVALO"
  L=$(livre_gb); [[ -z "$L" ]] && continue
  if [[ "$L" -le "$PISO_GB" ]]; then
    echo ""
    echo "*** PISO ATINGIDO: ${L} GB livres (piso ${PISO_GB} GB). Matando o grupo."
    kill -TERM -"$PG" 2>/dev/null
    sleep 5
    kill -KILL -"$PG" 2>/dev/null
    matou=1
    break
  fi
done
wait "$PG" 2>/dev/null
cod=$?

L1=$(livre_gb)
gasto=$(( L0 - ${L1:-$L0} ))
echo ""
echo "  livre antes: ${L0} GB | depois: ${L1:-?} GB | consumido: ${gasto} GB"
if [[ $matou -eq 1 ]]; then
  echo "*** INTERROMPIDO pela trava de disco, não por erro do comando."
  echo "*** O que estava sendo escrito pode estar pela metade — confira antes"
  echo "*** de reusar a saída."
  exit 9
fi
echo "  comando terminou com código $cod"
exit "$cod"
