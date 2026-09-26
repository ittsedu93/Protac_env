#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Desfaz a condição periódica de contorno antes da análise.
#
#   bash scripts/md_fix_pbc.sh <dir_md>
#
# O problema que isto resolve
# ---------------------------
# Ao cortar as lacunas do cristal (split_chain_gaps.py), a proteína deixou de
# ser UMA molécula e virou cinco. O GROMACS envolve cada molécula
# independentemente na caixa periódica, então segmentos da mesma proteína
# aparecem em lados opostos da caixa. Qualquer RMSD calculado sobre isso mede a
# aresta da caixa, não o movimento:
#
#     RMSD da proteína  30.00 Å      <- impossível para uma proteína enovelada
#     frações ancoradas  1.00        <- e incompatível com a linha de cima
#
# A contradição entre as duas linhas é a assinatura do artefato.
#
# A correção é a receita padrão: agrupar os solutos na mesma imagem periódica
# (-pbc cluster), depois impedir saltos entre quadros (-pbc nojump). A saída
# guarda só proteína + ligante, que é o que a análise usa — 5.934 átomos em vez
# de 111.809, então os arquivos ficam ~20x menores que a trajetória original.
# ---------------------------------------------------------------------------
set -uo pipefail

MD_DIR="${1:?uso: md_fix_pbc.sh <dir_md>}"
MD_DIR="$(cd "$MD_DIR" && pwd)"
GMX="${GMX_EXE:-/usr/local/gromacs/bin/gmx}"
GRP="${GRP_SOLUTOS:-Protein_LIG}"

[[ -s "$MD_DIR/grupos.ndx" ]] || { echo "não achei grupos.ndx em $MD_DIR"; exit 1; }

echo "=============================================================="
echo "Correção de PBC — $MD_DIR"
echo "  grupo de solutos: $GRP"
echo "=============================================================="

for rep in "$MD_DIR"/rep*/; do
  nome=$(basename "$rep")
  tpr="$rep/prod.tpr"; xtc="$rep/prod.xtc"
  [[ -s "$tpr" && -s "$xtc" ]] || { echo "  $nome: sem prod.tpr/prod.xtc, pulando"; continue; }
  if [[ -s "$rep/solutos.xtc" && -s "$rep/solutos.gro" ]]; then
    echo "  $nome: já corrigida, pulando"; continue
  fi

  echo -e "\n  $nome: agrupando os solutos na mesma imagem (-pbc cluster)"
  tmp="$rep/tmp_cluster.xtc"
  printf "%s\nSystem\n" "$GRP" | "$GMX" trjconv -s "$tpr" -f "$xtc" \
      -n "$MD_DIR/grupos.ndx" -o "$tmp" -pbc cluster > "$rep/pbc1.log" 2>&1
  [[ -s "$tmp" ]] || { echo "  *** -pbc cluster falhou; veja $rep/pbc1.log"; exit 1; }

  echo "  $nome: impedindo saltos entre quadros (-pbc nojump), só os solutos"
  printf "%s\n" "$GRP" | "$GMX" trjconv -s "$tpr" -f "$tmp" \
      -n "$MD_DIR/grupos.ndx" -o "$rep/solutos.xtc" -pbc nojump \
      > "$rep/pbc2.log" 2>&1
  printf "%s\n" "$GRP" | "$GMX" trjconv -s "$tpr" -f "$tmp" \
      -n "$MD_DIR/grupos.ndx" -o "$rep/solutos.gro" -dump 0 \
      > "$rep/pbc3.log" 2>&1

  if [[ -s "$rep/solutos.xtc" && -s "$rep/solutos.gro" ]]; then
    rm -f "$tmp"
    n=$(sed -n '2p' "$rep/solutos.gro" | tr -d ' ')
    echo "  $nome: OK — solutos.xtc + solutos.gro ($n átomos)"
  else
    echo "  *** $nome falhou; veja $rep/pbc2.log e $rep/pbc3.log"; exit 1
  fi
done

echo -e "\nPronto. Agora rode a análise — ela prefere solutos.xtc quando existe:"
echo "  python $(dirname "$0")/md_analyze.py --md-dir $MD_DIR"
