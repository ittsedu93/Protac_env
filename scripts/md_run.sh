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

# Roda um comando e confere que o arquivo esperado nasceu. Sem isto, um gmx
# que falha deixa o script seguir e todos os passos seguintes reclamam de
# "file does not exist" — o erro que aparece é o último, não o primeiro.
gmx_ok() {
  local saida="$1"; shift
  "$@" || { echo "*** falhou: $*"; exit 1; }
  [[ -s "$saida" ]] || { echo "*** $* não produziu $saida"; exit 1; }
}
CAND=$(jq_ candidate_id); CARGA=$(jq_ carga_formal)
RESNAME=$(jq_ resname); LIG_PDB=$(jq_ lig_pdb); RECEPTOR=$(jq_ receptor_pdb)

# parâmetros da metodologia
TEMP=310; PRESSAO=1.0; DT=0.002
# dt cresce ao longo do equilíbrio. Partir de 2 fs num sistema recém-montado é
# pedir estouro: o solvente ainda não acomodou e o PROTAC tem 19 torções.
DT_WARM=0.0005; DT_NVT=0.001
NS_WARM=0.02
# Nomes dos grupos de acoplamento térmico. Fixos de propósito: NÃO derivam do
# nome do resíduo, porque o resíduo escrito nas coordenadas ($RESNAME para o
# acpype, mas UNL no .gro) não é confiável. O build_index.py garante que
# grupos com estes nomes existam, qualquer que seja o nome do ligante.
GRP_SOLUTOS="Protein_LIG"; GRP_SOLVENTE="Water_and_ions"
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
if [[ ! -s "${RESNAME}.acpype/${RESNAME}_GMX.itp" ]] || \
   [[ ! -s "${RESNAME}.acpype/${RESNAME}_GMX.gro" ]]; then
  echo -e "\n[1/6] acpype — GAFF2/AM1-BCC (passo longo, ~10-40 min)"
  "$ACPYPE" -i "$LIG_PDB" -b "$RESNAME" -n "$CARGA" -a gaff2 -c bcc -o gmx \
    || { echo "*** acpype falhou. Carga formal ($CARGA) está correta?"; exit 1; }
  [[ -s "${RESNAME}.acpype/${RESNAME}_GMX.itp" ]] \
    || { echo "*** acpype terminou sem gerar o .itp"; exit 1; }
else
  echo -e "\n[1/6] acpype — já feito, pulando"
fi
ITP="$MD_DIR/${RESNAME}.acpype/${RESNAME}_GMX.itp"
LIG_GRO="$MD_DIR/${RESNAME}.acpype/${RESNAME}_GMX.gro"

# --- 2. topologia da proteína ---------------------------------------------
# A guarda olha complexo.gro, a saída FINAL desta etapa. Guardar por
# topol.top é errado: o pdb2gmx cria o arquivo antes de terminar, então uma
# execução que falhou no meio deixa topol.top no disco e a retomada pula a
# etapa inteira achando que deu certo.
if [[ ! -s complexo.gro ]]; then
  rm -f topol.top proteina.gro posre_prot.itp
  # Cristais têm cadeias laterais parcialmente resolvidas. O pdb2gmx recusa
  # a estrutura ("Incomplete ring in HIS68"), então completa antes.
  REC_FIX="$MD_DIR/receptor_fixed.pdb"
  if [[ ! -f "$REC_FIX" ]]; then
    echo -e "\n[2/6] completando resíduos incompletos do receptor"
    python "$(dirname "$0")/fix_receptor_for_md.py" \
        --receptor "$RECEPTOR" --ligante "$MD_DIR/protac.sdf" \
        ${MD_TRUNCAR_PERTO:+--truncar-perto} \
        --out "$REC_FIX" || {
      echo "*** não foi possível completar o receptor"; exit 1; }
  fi
  RECEPTOR="$REC_FIX"

  # Cristais têm alças não resolvidas. O pdb2gmx, vendo uma cadeia contínua,
  # liga os resíduos que LADEIAM a lacuna — uma ligação peptídica de 1,4 nm
  # onde cabem 0,133. A mola desmonta o sistema nos primeiros passos da
  # dinâmica, e o sintoma não parece com a causa: LINCS WARNING seguido de
  # cudaErrorIllegalAddress, porque a CUDA tocou em coordenadas que viraram
  # lixo. Cortar a cadeia nas lacunas resolve na raiz.
  REC_CAP="$MD_DIR/receptor_capped.pdb"
  if [[ ! -f "$REC_CAP" ]]; then
    echo -e "\n[2/6] cortando a cadeia nas lacunas do cristal"
    python "$(dirname "$0")/split_chain_gaps.py" --receptor "$RECEPTOR" \
        --ligante "$MD_DIR/protac.sdf" --out "$REC_CAP" || {
      echo "*** não consegui cortar as lacunas da cadeia"; exit 1; }
  fi
  RECEPTOR="$REC_CAP"

  echo -e "\n[2/6] pdb2gmx — AMBER ff14SB + TIP3P"
  "$GMX" pdb2gmx -f "$RECEPTOR" -ff amber14sb -water tip3p \
      -o proteina.gro -p topol.top -i posre_prot.itp -ignh \
      > pdb2gmx.log 2>&1
  cat pdb2gmx.log
  [[ -s proteina.gro ]] && [[ -s topol.top ]] || {
    echo "*** pdb2gmx não produziu proteina.gro/topol.top"; exit 1; }

  # Verificação, não confiança: uma ligação longa que sobre AQUI vira estouro
  # dez minutos depois, e o erro que aparece então não aponta para cá.
  if grep -q "Long Bond" pdb2gmx.log; then
    echo "*** o pdb2gmx ainda criou ligações longas:"
    grep "Long Bond" pdb2gmx.log | sed 's/^/      /'
    echo "*** Uma ligação peptídica mede 0,133 nm. Qualquer coisa acima de"
    echo "*** 0,25 nm é uma ligação que não existe, e ela VAI estourar a MD."
    if [[ "${MD_ACEITA_LIGACOES_LONGAS:-0}" == "1" ]]; then
      echo "*** MD_ACEITA_LIGACOES_LONGAS=1 — seguindo mesmo assim, por sua"
      echo "*** conta. Isto fica registrado no log para constar na tese."
    else
      echo "*** Para seguir de propósito: MD_ACEITA_LIGACOES_LONGAS=1"
      exit 1
    fi
  else
    echo "      nenhuma ligação longa — a cadeia está sã"
  fi

  # O .itp do acpype traz [ atomtypes ] e [ moleculetype ] juntos, e o
  # GROMACS exige todo atomtypes antes da primeira molécula. build_topology.py
  # separa os dois e insere cada um no lugar certo, verificando a ordem final.
  python "$(dirname "$0")/build_topology.py" \
      --topol topol.top --ligand-itp "$ITP" --ligand-gro "$LIG_GRO" \
      --resname "$RESNAME" || {
    echo "*** não consegui montar a topologia do complexo"; exit 1; }

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
  gmx_ok caixa.gro "$GMX" editconf -f complexo.gro -o caixa.gro -c \
      -d "$BOX_NM" -bt cubic
  gmx_ok solvatado.gro "$GMX" solvate -cp caixa.gro -cs spc216.gro \
      -p topol.top -o solvatado.gro
  printf "integrator = steep\nnsteps = 1\ncutoff-scheme = Verlet\ncoulombtype = PME\nrvdw = %s\nrcoulomb = %s\n" "$RC" "$RC" > ions.mdp
  gmx_ok ions.tpr "$GMX" grompp -f ions.mdp -c solvatado.gro -p topol.top \
      -o ions.tpr -maxwarn 2
  echo SOL | "$GMX" genion -s ions.tpr -o neutro.gro -p topol.top \
      -pname NA -nname CL -neutral
  [[ -s neutro.gro ]] || { echo "*** genion não produziu neutro.gro"; exit 1; }
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
constraint-algorithm = LINCS
lincs-order = 8
lincs-iter = 2"

# A MINIMIZAÇÃO é outra história. Vínculos ligados e água rígida fazem o LINCS
# tentar satisfazer restrições enquanto a geometria ainda está tensa — fonte
# clássica de estouro. Minimiza-se com tudo flexível; os vínculos entram depois.
comum_em="cutoff-scheme = Verlet
coulombtype = PME
rvdw = $RC
rcoulomb = $RC
constraints = none
define = -DFLEXIBLE"

# Dois termostatos, de propósito. O Nose-Hoover é o da metodologia e fica na
# PRODUÇÃO, que é onde ele importa: ele reproduz o ensemble canônico
# corretamente, mas não é robusto longe do equilíbrio — oscila, e num sistema
# recém-solvatado a oscilação vira estouro. O equilíbrio usa V-rescale, que é
# dissipativo e perdoa geometria ruim.
acopl_eq="tcoupl = V-rescale
tc-grps = ${GRP_SOLUTOS} ${GRP_SOLVENTE}
tau-t = 0.1 0.1
ref-t = $TEMP $TEMP"
acopl="tcoupl = Nose-Hoover
tc-grps = ${GRP_SOLUTOS} ${GRP_SOLVENTE}
tau-t = 1.0 1.0
ref-t = $TEMP $TEMP"
pacopl="pcoupl = Parrinello-Rahman
pcoupltype = isotropic
tau-p = 2.0
ref-p = $PRESSAO
compressibility = 4.5e-5"

# ATENÇÃO às unidades: `passos` recebe NANOssegundos, `passos_ps` recebe
# PICOssegundos. Confundir as duas escreve nstxout-compressed = 1e8 num run de
# 1e8 passos — um único frame em 200 ns, e a análise sem nada para medir.
passos() { python3 -c "print(int($1*1000/$2))"; }
passos_ps() { python3 -c "print(int($1/$2))"; }

# 1) minimização grosseira, e 2) uma fina que tira a tensão que sobrou
printf "integrator = steep\nemtol = 1000.0\nemstep = 0.01\nnsteps = 50000\n%s\n" \
    "$comum_em" > em.mdp
printf "integrator = steep\nemtol = 100.0\nemstep = 0.001\nnsteps = 50000\n%s\n" \
    "$comum_em" > em2.mdp

# 3) aquecimento: dt de 0,5 fs partindo de 100 K, com o soluto preso
printf "integrator = md\ndt = %s\nnsteps = %d\n%s\n%s\npcoupl = no\ngen-vel = yes\ngen-temp = 100\ndefine = -DPOSRES\n" \
    "$DT_WARM" "$(passos "$NS_WARM" "$DT_WARM")" "$comum" "$acopl_eq" > warm.mdp

# 4) NVT a 1 fs, continuando o aquecimento
printf "integrator = md\ndt = %s\nnsteps = %d\n%s\n%s\npcoupl = no\ngen-vel = no\ncontinuation = yes\ndefine = -DPOSRES\n" \
    "$DT_NVT" "$(passos "$NS_NVT" "$DT_NVT")" "$comum" "$acopl_eq" > nvt.mdp

# 5) NPT a 2 fs, ainda com V-rescale e o soluto preso
printf "integrator = md\ndt = %s\nnsteps = %d\n%s\n%s\n%s\ngen-vel = no\ncontinuation = yes\ndefine = -DPOSRES\n" \
    "$DT" "$(passos "$NS_NPT" "$DT")" "$comum" "$acopl_eq" "$pacopl" > npt.mdp

# 6) produção: aqui sim Nose-Hoover + Parrinello-Rahman, sem restrições
printf "integrator = md\ndt = %s\nnsteps = %d\n%s\n%s\n%s\ngen-vel = no\nnstxout-compressed = %d\nnstenergy = %d\n" \
    "$DT" "$(passos "$NS_PROD" "$DT")" "$comum" "$acopl" "$pacopl" \
    "$(passos_ps "$SAVE_PS" "$DT")" "$(passos_ps "$SAVE_PS" "$DT")" > prod.mdp

# Grupos de acoplamento: proteína+ligante contra água+íons.
#
# Não se escreve isto à mão. Os números dos grupos mudam com a composição do
# sistema e o nome do resíduo do ligante não é o que se pediu ao acpype — era
# assim que a versão anterior falhava, procurando um `Protein_PTC` que o
# make_ndx nunca criaria porque o grupo se chamava `UNL`. O build_index.py
# descobre os números pelos nomes, monta as fusões e confere o resultado.
# Um grupos.ndx que ficou pela metade numa tentativa anterior é pior que
# nenhum: ele existe, a retomada o aceita, e o grompp morre lá na frente
# reclamando de tc-grps. A guarda confere o conteúdo, não só a existência.
if [[ -s grupos.ndx ]] && \
   ! { grep -q "\[ ${GRP_SOLUTOS} \]" grupos.ndx && \
       grep -q "\[ ${GRP_SOLVENTE} \]" grupos.ndx; }; then
  echo "      grupos.ndx incompleto (sem ${GRP_SOLUTOS}/${GRP_SOLVENTE}) — refazendo"
  rm -f grupos.ndx
fi
if [[ ! -s grupos.ndx ]]; then
  echo "      montando grupos de acoplamento"
  python "$(dirname "$0")/build_index.py" --gro neutro.gro --out grupos.ndx \
      --gmx "$GMX" --nome-solutos "$GRP_SOLUTOS" \
      --nome-solvente "$GRP_SOLVENTE" || {
    echo "*** não consegui montar os grupos de acoplamento."
    echo "*** Rode à mão, crie ${GRP_SOLUTOS} e ${GRP_SOLVENTE}, salve como grupos.ndx:"
    echo "***   $GMX make_ndx -f neutro.gro -o grupos.ndx"
    exit 1; }
else
  echo "      grupos.ndx — já feito, pulando"
fi

# --- 5. minimização e equilíbrio ------------------------------------------
# Escada de offload. O ideal é a GPU fazer tudo, mas parte das combinações da
# metodologia o GROMACS recusa em `-update gpu` (o termostato Nose-Hoover e as
# restrições de posição do equilíbrio são os casos conhecidos), e a recusa é um
# erro de configuração que mata a etapa em segundos. Em vez de escolher de
# antemão, tenta-se do mais rápido para o mais conservador — e só se a mensagem
# do próprio GROMACS for sobre offload. Um erro de física (LINCS, explosão)
# falha de uma vez, sem repetir horas de simulação.
GPU="-nb gpu -pme gpu -bonded gpu -update gpu"
GPU_MENOS="-nb gpu -pme gpu"

# A escada usada em cada etapa. A minimização tem a sua porque `steep` NÃO é um
# integrador dinâmico: não integra equações de movimento, e por isso o PME e o
# update na GPU são inválidos ali por definição, não por combinação infeliz
# ("PME GPU does not support: Non-dynamical integrator"). Só o cálculo de
# não-ligadas aproveita a placa na minimização.
ESCADA_MD=("$GPU" "$GPU_MENOS" " ")
ESCADA_EM=("-nb gpu" " ")
ESCADA=("${ESCADA_MD[@]}")

mdrun_ok() {  # mdrun_ok <deffnm>
  local nome="$1"; shift
  local log="mdrun_${nome}.out" extra
  local flags
  for flags in "${ESCADA[@]}"; do
    # retomada: havendo checkpoint, continua de onde parou em vez de reiniciar
    extra=""
    [[ -f "${nome}.cpt" ]] && extra="-cpi ${nome}.cpt -append"
    echo "      mdrun $nome ${extra:+(retomando) }[offload:${flags# }]"
    if "$GMX" mdrun -deffnm "$nome" $extra $flags "$@" 2>&1 | tee "$log"; then
      return 0
    fi
    # As recusas de offload não têm uma redação só. Já vistas:
    #   "...following condition(s) were not satisfied"   (update na GPU)
    #   "PME GPU does not support: Non-dynamical integrator"
    #   "Inconsistency in user input"
    # Procurar "not supported" não casa com "does not support" — é assim que a
    # escada deixou de descer numa recusa que era, sim, de offload.
    if grep -qi "gpu" "$log" && grep -qiE \
         "not satisfied|not support|cannot be used|cannot compute|incompatible|inconsistency in user input" \
         "$log"; then
      echo "      o GROMACS recusou este offload — tentando com menos GPU"
      continue
    fi
    # O erro do GROMACS já está no log, mas dezenas de linhas acima do fim —
    # quem olha `tail` vê o rodapé e não a causa. Repete a causa AQUI, no fim.
    echo "      mdrun falhou por motivo que não é offload. O erro do GROMACS:"
    sed -n "/^-\{20,\}$/,\$p" "$log" | head -30 | sed 's/^/      | /'
    if grep -qiE "LINCS WARNING|can not be settled|Too many LINCS" "$log"; then
      echo "      >>> Isto é INSTABILIDADE NUMÉRICA, não configuração: alguma"
      echo "      >>> geometria do sistema está tensa."
      # "átomo 539" não diz nada; o resíduo diz tudo. Os índices são do
      # sistema SOLVATADO, então o mapeamento sai do neutro.gro.
      python "$(dirname "$0")/explain_lincs.py" --gro "$MD_DIR/neutro.gro" \
          --log "$log" 2>/dev/null || true
    fi
    echo "      (completo em $MD_DIR/$log)"
    return 1
  done
  return 1
}

etapa() {  # etapa <nome> <mdp> <gro_entrada> [extra_grompp]
  local nome="$1" mdp="$2" entrada="$3"; shift 3
  if [[ -f "${nome}.gro" ]]; then echo "      $nome — já feito"; return 0; fi
  "$GMX" grompp -f "$mdp" -c "$entrada" -p topol.top -n grupos.ndx \
      -o "${nome}.tpr" -maxwarn 5 "$@" || return 1
  mdrun_ok "$nome" || return 1
}

echo -e "\n[5/6] minimização e equilíbrio"

# A produção já começada é prova de que o equilíbrio foi feito. Se o npt.gro
# tiver sido perdido (apagado por engano, disco cheio), ele é recuperável de
# qualquer prod.tpr: o .tpr guarda o estado completo com que a réplica começou.
# Recuperar é melhor que refazer — um equilíbrio novo geraria uma estrutura de
# partida DIFERENTE da que as réplicas já rodadas usaram, e as réplicas
# deixariam de ser o que a metodologia diz que são: mesmas coordenadas,
# velocidades independentes.
if [[ ! -s npt.gro ]] && compgen -G "rep*/prod.tpr" > /dev/null; then
  TPR_REF=$(ls -1 rep*/prod.tpr | head -1)
  echo "      npt.gro ausente, mas $TPR_REF existe — recuperando dele"
  "$GMX" editconf -f "$TPR_REF" -o npt.gro || {
    echo "*** não consegui recuperar npt.gro de $TPR_REF"; exit 1; }
fi

if [[ -s npt.gro ]]; then
  echo "      equilíbrio — já feito, pulando"
else
ESCADA=("${ESCADA_EM[@]}")
etapa em  em.mdp  neutro.gro || { echo "*** minimização falhou"; exit 1; }
etapa em2 em2.mdp em.gro    || { echo "*** minimização fina falhou"; exit 1; }

# A minimização "termina" mesmo num sistema tenso. O Fmax final é o que
# distingue relaxou de desistiu — entrar no NVT com Fmax na casa de 1e4
# kJ/mol/nm é entrar num estouro, e o estouro custa minutos de GPU para dizer
# o que esta linha diz em um segundo.
if [[ -f em2.log ]]; then
  fmax=$(awk '/Maximum force/ {print $4; exit}' em2.log)
  conv=""; grep -q "converged to Fmax" em2.log && conv=" (convergiu)"
  echo "      minimização: Fmax = ${fmax:-?} kJ/mol/nm${conv}"
  if [[ -n "${fmax:-}" ]] && awk -v f="$fmax" 'BEGIN{exit !(f>1000)}'; then
    echo "      [ATENÇÃO] Fmax alto: o sistema segue tenso. Se o NVT estourar,"
    echo "      [ATENÇÃO] o problema é a geometria de partida, não o NVT."
  fi
fi

ESCADA=("${ESCADA_MD[@]}")
# A referência das restrições de posição é sempre a estrutura minimizada, não a
# etapa anterior: encadear referências deixa o soluto derivar de degrau em
# degrau e chegar na produção longe de onde o docking o colocou.
etapa warm warm.mdp em2.gro -r em2.gro || { echo "*** aquecimento falhou"; exit 1; }
etapa nvt nvt.mdp warm.gro -r em2.gro -t warm.cpt || { echo "*** NVT falhou"; exit 1; }
etapa npt npt.mdp nvt.gro  -r em2.gro -t nvt.cpt  || { echo "*** NPT falhou"; exit 1; }
fi

# --- 6. produção, N réplicas com sementes diferentes ----------------------
echo -e "\n[6/6] produção — ${NS_PROD} ns x ${N_REP} réplicas"
for i in $(seq 1 "$N_REP"); do
  rep="rep${i}"
  if [[ -f "${rep}/prod.gro" ]]; then echo "      $rep — já concluída"; continue; fi
  mkdir -p "$rep"
  seed=$((1000 + (i - 1) * 97))
  sed "s/^gen-vel = no/gen-vel = yes\ngen-seed = $seed\ngen-temp = $TEMP/" prod.mdp > "$rep/prod.mdp"
  # O -t é dispensável: com gen-vel = yes e semente própria, as velocidades
  # são geradas de novo. Exigir o checkpoint faria uma réplica falhar por falta
  # de um arquivo cujo conteúdo ia ser descartado.
  T_NPT=(); [[ -s npt.cpt ]] && T_NPT=(-t npt.cpt)
  "$GMX" grompp -f "$rep/prod.mdp" -c npt.gro "${T_NPT[@]}" -p topol.top \
      -n grupos.ndx -o "$rep/prod.tpr" -maxwarn 5 || exit 1
  ( cd "$rep" && mdrun_ok prod ) || exit 1
  echo "      $rep concluída"
done

echo -e "\n=============================================================="
echo "MD CONCLUÍDA — $CAND"
echo "=============================================================="
echo "Trajetórias: $MD_DIR/rep*/prod.xtc"
echo
echo "Análise (env mdtools):"
echo "  python $(dirname "$0")/md_analyze.py --md-dir $MD_DIR"
