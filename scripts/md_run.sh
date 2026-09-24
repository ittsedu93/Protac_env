#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Fase 8b — MD do nível (ii): E3 ligase + recrutador-linker-warhead.
#
#   bash scripts/md_run.sh <dir_md> [n_replicas]
#
# Lê md_sistema.json (escrito por md_prepare.py), parametriza o PROTAC com
# GAFF2/AM1-BCC via acpype, monta o sistema em AMBER ff14SB + TIP3P, e roda
# minimização, NVT, NPT e produção no GROMACS.
#
# Retomável: cada etapa pula se a saída já existe. Uma queda custa a etapa em
# andamento, não a réplica inteira.
#
# A GPU faz nonbonded, PME, bonded e update. Numa Blackwell de 32 GB o sistema
# inteiro cabe na placa e a CPU fica só coordenando.
# ---------------------------------------------------------------------------
set -uo pipefail

MD_DIR="${1:?uso: md_run.sh <dir_md> [n_replicas]}"
N_REP="${2:-3}"
MD_DIR="$(cd "$MD_DIR" && pwd)"
INFO="$MD_DIR/md_sistema.json"
[[ -f "$INFO" ]] || { echo "não achei $INFO — rode md_prepare.py antes"; exit 1; }

GMX="${GMX_EXE:-/usr/local/gromacs/bin/gmx}"
ACPYPE="${ACPYPE_EXE:-/home/soberano/miniconda3/envs/mdtools/bin/acpype}"

jq_() { python3 -c "import json;print(json.load(open('$INFO'))['$1'])"; }
CAND=$(jq_ candidate_id); CARGA=$(jq_ carga_formal)
RESNAME=$(jq_ resname); LIG_PDB=$(jq_ lig_pdb); RECEPTOR=$(jq_ receptor_pdb)

# parâmetros da metodologia
TEMP=310; PRESSAO=1.0; DT=0.002
NS_PROD="${NS_PROD:-200}"; NS_NPT="${NS_NPT:-5}"; NS_NVT="${NS_NVT:-0.1}"
BOX_NM=1.3; RC=0.9
SAVE_PS=200

echo "=============================================================="
echo "MD nível (ii) — E3 + recrutador-linker-warhead"
echo "  candidato : $CAND"
echo "  carga     : $CARGA"
echo "  produção  : ${NS_PROD} ns x ${N_REP} réplicas"
echo "=============================================================="

cd "$MD_DIR"

# --- 1. parametrização do ligante (GAFF2 + AM1-BCC) -----------------------
if [[ ! -f "${RESNAME}.acpype/${RESNAME}_GMX.itp" ]]; then
  echo -e "\n[1/6] acpype — GAFF2/AM1-BCC (passo longo, ~10-40 min)"
  "$ACPYPE" -i "$LIG_PDB" -b "$RESNAME" -n "$CARGA" -a gaff2 -c bcc -o gmx \
    || { echo "*** acpype falhou. Carga formal ($CARGA) está correta?"; exit 1; }
else
  echo -e "\n[1/6] acpype — já feito, pulando"
fi
ITP="$MD_DIR/${RESNAME}.acpype/${RESNAME}_GMX.itp"
LIG_GRO="$MD_DIR/${RESNAME}.acpype/${RESNAME}_GMX.gro"

# --- 2. topologia da proteína ---------------------------------------------
if [[ ! -f topol.top ]]; then
  echo -e "\n[2/6] pdb2gmx — AMBER ff14SB + TIP3P"
  "$GMX" pdb2gmx -f "$RECEPTOR" -ff amber14sb -water tip3p \
      -o proteina.gro -p topol.top -i posre_prot.itp -ignh \
    || { echo "*** pdb2gmx falhou (resíduo não reconhecido?)"; exit 1; }

  # insere o ligante na topologia: #include depois dos includes do campo,
  # e a linha da molécula no fim de [ molecules ]
  python3 - "$ITP" <<'PY'
import sys, re
from pathlib import Path
itp = sys.argv[1]
top = Path("topol.top"); t = top.read_text()
if "PTC_GMX.itp" not in t:
    linhas = t.splitlines()
    # depois do último #include de campo de força
    ult = max(i for i, l in enumerate(linhas) if l.strip().startswith('#include'))
    linhas.insert(ult + 1, f'#include "{itp}"')
    t = "\n".join(linhas)
    t = t.rstrip() + "\nPTC                 1\n"
    top.write_text(t)
    print("      ligante inserido em topol.top")
PY

  echo "      juntando proteína e ligante"
  python3 - "$LIG_GRO" <<'PY'
import sys
from pathlib import Path
prot = Path("proteina.gro").read_text().splitlines()
lig  = Path(sys.argv[1]).read_text().splitlines()
n_p, n_l = int(prot[1]), int(lig[1])
saida = [prot[0], f"{n_p + n_l:5d}"] + prot[2:2+n_p] + lig[2:2+n_l] + [prot[-1]]
Path("complexo.gro").write_text("\n".join(saida) + "\n")
print(f"      complexo.gro: {n_p} + {n_l} = {n_p+n_l} átomos")
PY
else
  echo -e "\n[2/6] topologia — já feita, pulando"
fi

# --- 3. caixa, solvente e íons --------------------------------------------
if [[ ! -f neutro.gro ]]; then
  echo -e "\n[3/6] caixa cúbica (${BOX_NM} nm), TIP3P, neutralização"
  "$GMX" editconf -f complexo.gro -o caixa.gro -c -d "$BOX_NM" -bt cubic
  "$GMX" solvate -cp caixa.gro -cs spc216.gro -p topol.top -o solvatado.gro
  printf "integrator = steep\nnsteps = 1\ncutoff-scheme = Verlet\ncoulombtype = PME\nrvdw = %s\nrcoulomb = %s\n" "$RC" "$RC" > ions.mdp
  "$GMX" grompp -f ions.mdp -c solvatado.gro -p topol.top -o ions.tpr -maxwarn 2
  echo SOL | "$GMX" genion -s ions.tpr -o neutro.gro -p topol.top \
      -pname NA -nname CL -neutral
else
  echo -e "\n[3/6] solvatação — já feita, pulando"
fi

# --- 4. .mdp da metodologia ------------------------------------------------
echo -e "\n[4/6] escrevendo .mdp"
comum="cutoff-scheme = Verlet
coulombtype = PME
rvdw = $RC
rcoulomb = $RC
constraints = h-bonds
constraint-algorithm = LINCS"
acopl="tcoupl = Nose-Hoover
tc-grps = Protein_${RESNAME} Water_and_ions
tau-t = 1.0 1.0
ref-t = $TEMP $TEMP"
pacopl="pcoupl = Parrinello-Rahman
pcoupltype = isotropic
tau-p = 2.0
ref-p = $PRESSAO
compressibility = 4.5e-5"

printf "integrator = steep\nemtol = 1000.0\nemstep = 0.01\nnsteps = 50000\n%s\n" "$comum" > em.mdp
printf "integrator = md\ndt = %s\nnsteps = %d\n%s\n%s\npcoupl = no\ngen-vel = yes\ngen-temp = %s\ndefine = -DPOSRES\n" \
    "$DT" "$(python3 -c "print(int($NS_NVT*1000/$DT))")" "$comum" "$acopl" "$TEMP" > nvt.mdp
printf "integrator = md\ndt = %s\nnsteps = %d\n%s\n%s\n%s\ngen-vel = no\ndefine = -DPOSRES\n" \
    "$DT" "$(python3 -c "print(int($NS_NPT*1000/$DT))")" "$comum" "$acopl" "$pacopl" > npt.mdp
printf "integrator = md\ndt = %s\nnsteps = %d\n%s\n%s\n%s\ngen-vel = no\nnstxout-compressed = %d\nnstenergy = %d\n" \
    "$DT" "$(python3 -c "print(int($NS_PROD*1000/$DT))")" "$comum" "$acopl" "$pacopl" \
    "$(python3 -c "print(int($SAVE_PS/$DT))")" "$(python3 -c "print(int($SAVE_PS/$DT))")" > prod.mdp

# grupos de acoplamento: proteína+ligante contra água+íons
if [[ ! -f grupos.ndx ]]; then
  printf "1 | 13\nname 22 Protein_%s\n14 | 15\nname 23 Water_and_ions\nq\n" "$RESNAME" \
    | "$GMX" make_ndx -f neutro.gro -o grupos.ndx > make_ndx.log 2>&1 || true
  grep -q "Protein_${RESNAME}" grupos.ndx 2>/dev/null || {
    echo "      [ATENÇÃO] os índices dos grupos variam com o sistema."
    echo "      Confira grupos.ndx: os tc-grps dos .mdp precisam existir nele."
    echo "      Abra com: $GMX make_ndx -f neutro.gro"
  }
fi

# --- 5. minimização e equilíbrio ------------------------------------------
GPU="-nb gpu -pme gpu -bonded gpu -update gpu"
etapa() {  # etapa <nome> <mdp> <gro_entrada> [extra_grompp]
  local nome="$1" mdp="$2" entrada="$3"; shift 3
  if [[ -f "${nome}.gro" ]]; then echo "      $nome — já feito"; return 0; fi
  "$GMX" grompp -f "$mdp" -c "$entrada" -p topol.top -n grupos.ndx \
      -o "${nome}.tpr" -maxwarn 5 "$@" || return 1
  "$GMX" mdrun -deffnm "$nome" $GPU || return 1
}

echo -e "\n[5/6] minimização e equilíbrio"
etapa em  em.mdp  neutro.gro || { echo "*** minimização falhou"; exit 1; }
etapa nvt nvt.mdp em.gro -r em.gro || { echo "*** NVT falhou"; exit 1; }
etapa npt npt.mdp nvt.gro -r nvt.gro -t nvt.cpt || { echo "*** NPT falhou"; exit 1; }

# --- 6. produção, N réplicas com sementes diferentes ----------------------
echo -e "\n[6/6] produção — ${NS_PROD} ns x ${N_REP} réplicas"
for i in $(seq 1 "$N_REP"); do
  rep="rep${i}"
  if [[ -f "${rep}/prod.gro" ]]; then echo "      $rep — já concluída"; continue; fi
  mkdir -p "$rep"
  seed=$((1000 + (i - 1) * 97))
  sed "s/^gen-vel = no/gen-vel = yes\ngen-seed = $seed\ngen-temp = $TEMP/" prod.mdp > "$rep/prod.mdp"
  "$GMX" grompp -f "$rep/prod.mdp" -c npt.gro -t npt.cpt -p topol.top \
      -n grupos.ndx -o "$rep/prod.tpr" -maxwarn 5 || exit 1
  ( cd "$rep" && "$GMX" mdrun -deffnm prod $GPU ) || exit 1
  echo "      $rep concluída"
done

echo -e "\n=============================================================="
echo "MD CONCLUÍDA — $CAND"
echo "=============================================================="
echo "Trajetórias: $MD_DIR/rep*/prod.xtc"
echo
echo "Análise (env mdtools):"
echo "  python $(dirname "$0")/md_analyze.py --md-dir $MD_DIR"
