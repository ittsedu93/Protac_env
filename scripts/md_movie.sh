#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Filme .mp4 da trajetória, no ChimeraX headless.
#
#   bash scripts/md_movie.sh <dir_md> [replica] [segundos]
#
# Mostra o que o número diz: o bolso do recrutador parado (a câmera está
# alinhada nele) enquanto o domínio distal gira. O PROTAC em bastões, a E3 em
# fita, e uma faixa de cor marcando o sítio.
#
# Precisa de movie_<rep>.pdb, escrito por md_export_structures.py.
# ---------------------------------------------------------------------------
set -uo pipefail

MD_DIR="${1:?uso: md_movie.sh <dir_md> [replica] [segundos]}"
MD_DIR="$(cd "$MD_DIR" && pwd)"
REP="${2:-rep1}"
SEGUNDOS="${3:-12}"
CHIMERAX="${CHIMERAX_EXE:-/usr/bin/chimerax}"
EST="$MD_DIR/structures"
PDB="$EST/movie_${REP}.pdb"

[[ -x "$CHIMERAX" ]] || { echo "não achei o ChimeraX em $CHIMERAX"; exit 1; }
[[ -s "$PDB" ]] || {
  echo "não achei $PDB"
  echo "rode antes: python scripts/md_export_structures.py --md-dir $MD_DIR"
  exit 1; }

N=$(grep -c "^MODEL" "$PDB" || echo 0)
[[ "$N" -gt 1 ]] || { echo "$PDB não tem múltiplos modelos"; exit 1; }
FPS=$(python3 -c "print(max(5, round($N / $SEGUNDOS)))")
SAIDA="$EST/${REP}_trajectory.mp4"

echo "Filme: $N quadros, ${SEGUNDOS}s a ${FPS} fps -> $(basename "$SAIDA")"

CXC="$EST/.movie_${REP}.cxc"
cat > "$CXC" <<CX
open $PDB
set bgColor white
hide atoms
show cartoon
color #d8d8d4 target c
# o ligante: o PROTAC inteiro, que no .gro é o que não é proteína
select ~protein
show sel atoms
style sel stick
color sel #2a78d6 target a
# o sítio do recrutador, para o olho saber onde olhar
select zone sel 12 residues true
color sel #eb6834 target c
~select
lighting soft
graphics silhouettes true width 1.5
view
movie record supersample 3
coordset #1 1,$N
wait $N
movie encode "$SAIDA" framerate $FPS quality high
exit
CX

"$CHIMERAX" --nogui --exit "$CXC" > "$EST/movie_${REP}.log" 2>&1
if [[ -s "$SAIDA" ]]; then
  echo "OK — $SAIDA ($(du -h "$SAIDA" | cut -f1))"
  rm -f "$CXC"
else
  echo "*** o ChimeraX não produziu o mp4; veja $EST/movie_${REP}.log"
  echo "*** (falta de codec é o motivo mais comum; o log diz)"
  exit 1
fi
