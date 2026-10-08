#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Roda um comando com TETO DE CONSUMO de disco, e mata ao passar dele.
#
#   bash scripts/rodar_com_teto_de_disco.sh --max-consumo-gb 15 --medir <dir> \
#        -- python algum_script.py args
#
# O limite é QUANTO O COMANDO PODE ESCREVER, não quanto sobra no disco. As duas
# leituras dão resultados opostos — com 77 GB livres, "consumir no máximo 15" e
# "deixar 15 livres" diferem por 47 GB — então o número é nomeado pelo que ele
# é, e os dois aparecem na tela antes de começar.
#
# Dois limites independentes, e o mais apertado vence
# ---------------------------------------------------
#   --max-consumo-gb   o que ESTE comando escreveu. É o limite do operador.
#   --piso-gb          o livre absoluto do disco. Protege a máquina mesmo que
#                      quem esteja enchendo seja outra pessoa.
#
# Dois instrumentos, porque a máquina é COMPARTILHADA
# ---------------------------------------------------
# Medir consumo pelo espaço livre do `df` está errado numa máquina com outros
# usuários: o `df` cai quando qualquer um escreve, e sobe quando qualquer um
# apaga. Isso culparia o nosso comando pelo que não é dele, e — pior —
# ESCONDERIA o consumo dele se outra pessoa apagasse arquivos no mesmo
# intervalo. Então:
#
#   du do diretório medido   atribui o consumo a NÓS, com precisão.
#                            Custa stat() em muitos arquivos, então é o
#                            instrumento lento (--intervalo-du, 30 s)
#   df do sistema de arquivos  barato, pega qualquer escrita de qualquer um,
#                            então é o instrumento rápido (--intervalo, 2 s)
#
# A queda do `df` também conta contra o teto de consumo, de propósito: ela pode
# ser culpa de outra pessoa, e nesse caso paramos à toa. Parar à toa custa um
# relançamento; não parar custa o disco da máquina do laboratório.
#
# O QUE ESTA TRAVA GARANTE, E O QUE NÃO GARANTE
# ---------------------------------------------
# Ela NÃO é um limite absoluto do kernel. Ela mede a cada poucos segundos e
# mata; num NVMe que escreve ~2 GB/s, a ultrapassagem teórica entre duas
# medições é de alguns GB, não zero. Dizer "absolutamente não vai exceder"
# seria mentira sobre o mecanismo.
#
# O que ela faz, e que é honesto:
#   - `ulimit -f` no filho, que é limite DE KERNEL por arquivo: nenhum arquivo
#     individual passa de --max-arquivo-gb, e o processo recebe SIGXFSZ
#   - medição a cada 2 s, e morte do GRUPO de processos (não do processo: um
#     script que lançou filhos deixaria os filhos escrevendo)
#   - recusa de começar se o teto já não cabe no livre
#
# Para um limite absoluto, de kernel, sem janela nenhuma, o caminho é quota de
# disco (`setquota`) ou um sistema de arquivos próprio montado com tamanho
# fixo — os dois exigem root, ou seja, o administrador da máquina.
#
# A outra metade da resposta não está aqui: é escolher um teto PEQUENO. Rodar
# com teto de 2 GB numa etapa cuja estimativa é 0,5 GB arrisca 2 GB. Diminuir
# o raio do estrago vale mais que confiar no freio.
# ---------------------------------------------------------------------------
set -uo pipefail

MAX_GB=""
PISO_GB=10
MEDIR=""
ONDE=""
INTERVALO=2
INTERVALO_DU=30
MAX_ARQ_GB=20
while [[ $# -gt 0 ]]; do
  case "$1" in
    --max-consumo-gb) MAX_GB="$2"; shift 2;;
    --piso-gb)        PISO_GB="$2"; shift 2;;
    --medir)          MEDIR="$2"; shift 2;;
    --onde)           ONDE="$2"; shift 2;;
    --intervalo)      INTERVALO="$2"; shift 2;;
    --intervalo-du)   INTERVALO_DU="$2"; shift 2;;
    --max-arquivo-gb) MAX_ARQ_GB="$2"; shift 2;;
    --) shift; break;;
    *) echo "opção desconhecida: $1"; exit 1;;
  esac
done
[[ $# -gt 0 ]] || {
  echo "uso: $0 --max-consumo-gb N [--medir DIR] [--piso-gb N] -- <comando>"; exit 1; }
[[ -n "$MAX_GB" ]] || { echo "*** --max-consumo-gb é obrigatório: sem teto eu não rodo"; exit 1; }
MEDIR="${MEDIR:-$PWD}"
ONDE="${ONDE:-$MEDIR}"

# EM MEGABYTES, e isto não é detalhe. Com `df -BG` o livre vem arredondado
# para GB inteiros, e no teste a trava disparou com 298 MB escritos porque o
# arredondamento fez o livre "cair 1 GB". Errar para o lado de parar é o certo,
# mas parar à toa é relançamento perdido — e um relatório final dizendo
# "consumido: 0 GB" de uma etapa que escreveu 300 MB não informa nada.
livre_mb() { df -BM --output=avail "$ONDE" 2>/dev/null | tail -1 | tr -dc '0-9'; }
usado_mb() { du -sBM "$MEDIR" 2>/dev/null | cut -f1 | tr -dc '0-9'; }
gb() { awk -v m="${1:-0}" 'BEGIN{printf "%.1f", m/1024}'; }

MAX_MB=$(( MAX_GB * 1024 ))
PISO_MB=$(( PISO_GB * 1024 ))
L0=$(livre_mb)
U0=$(usado_mb)
[[ -n "$L0" && -n "$U0" ]] || {
  echo "*** não consegui medir disco em $ONDE / $MEDIR — não rodo às cegas"; exit 1; }

PISO_DERIVADO_MB=$(( L0 - MAX_MB ))
echo "=============================================================="
echo "TETO DE CONSUMO: ${MAX_GB} GB"
echo "=============================================================="
echo "  livre agora no disco de $ONDE : $(gb "$L0") GB"
echo "  já ocupado por $MEDIR : $(gb "$U0") GB"
echo ""
echo "  o comando pode escrever no máximo ..... ${MAX_GB} GB"
echo "  deixando livre, no pior caso .......... $(gb "$PISO_DERIVADO_MB") GB"
echo "  e eu mato antes disso se o livre tocar  ${PISO_GB} GB"
echo "  nenhum arquivo individual passa de .... ${MAX_ARQ_GB} GB (limite de kernel)"
echo "  medindo: df a cada ${INTERVALO}s, du a cada ${INTERVALO_DU}s"
echo "=============================================================="
if [[ "$MAX_MB" -ge $(( L0 - PISO_MB )) ]]; then
  echo "*** O teto de ${MAX_GB} GB não cabe: com $(gb "$L0") GB livres e piso"
  echo "*** de ${PISO_GB} GB, o máximo que cabe é $(gb $(( L0 - PISO_MB ))) GB."
  exit 2
fi

# ulimit -f é em blocos de 1 KiB. É limite de KERNEL, não de medição: o
# processo recebe SIGXFSZ na escrita que passaria do tamanho, sem janela.
BLOCOS=$(( MAX_ARQ_GB * 1024 * 1024 ))
( ulimit -f "$BLOCOS" 2>/dev/null; exec setsid "$@" ) &
PG=$!
echo "  rodando (grupo $PG)"

matou=""
voltas=0
while kill -0 "$PG" 2>/dev/null; do
  sleep "$INTERVALO"
  voltas=$((voltas + 1))

  L=$(livre_mb)
  if [[ -n "$L" ]]; then
    if [[ "$L" -le "$PISO_MB" ]]; then
      matou="piso absoluto: $(gb "$L") GB livres (piso ${PISO_GB} GB)"
    elif [[ $(( L0 - L )) -ge "$MAX_MB" ]]; then
      # Pode ser escrita de outra pessoa. Paramos de qualquer forma: parar à
      # toa custa um relançamento, não parar custa o disco do laboratório.
      matou="o livre caiu $(gb $(( L0 - L ))) GB, teto de consumo ${MAX_GB} GB"
    fi
  fi

  if [[ -z "$matou" ]] && (( voltas % (INTERVALO_DU / INTERVALO + 1) == 0 )); then
    U=$(usado_mb)
    if [[ -n "$U" ]] && [[ $(( U - U0 )) -ge "$MAX_MB" ]]; then
      matou="o diretório cresceu $(gb $(( U - U0 ))) GB, teto ${MAX_GB} GB"
    fi
  fi

  if [[ -n "$matou" ]]; then
    echo ""
    echo "*** TETO ATINGIDO: $matou"
    echo "*** Matando o grupo $PG."
    kill -TERM -"$PG" 2>/dev/null
    sleep 5
    kill -KILL -"$PG" 2>/dev/null
    break
  fi
done
wait "$PG" 2>/dev/null
cod=$?

L1=$(livre_mb); U1=$(usado_mb)
echo ""
echo "  livre: $(gb "$L0") -> $(gb "${L1:-$L0}") GB"
echo "  $MEDIR: $(gb "$U0") -> $(gb "${U1:-$U0}") GB"
echo "  CONSUMIDO pelo comando: $(( ${U1:-$U0} - U0 )) MB"\
     "($(gb $(( ${U1:-$U0} - U0 ))) GB de um teto de ${MAX_GB} GB)"
if [[ -n "$matou" ]]; then
  echo "*** INTERROMPIDO pela trava de disco, não por erro do comando."
  echo "*** O que estava sendo escrito pode estar pela metade — confira antes"
  echo "*** de reusar a saída."
  exit 9
fi
echo "  comando terminou com código $cod"
exit "$cod"
