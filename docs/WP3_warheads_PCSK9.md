# Warheads PCSK9 (feniletilamina) na Fase 3 do `protac_wp1_wp3_pipeline.ipynb`

Assume WP1 e WP2 concluídos (recrutadores VHL/CRBN com exit vector definido;
sub-complexos recrutador–linker validados, com cap dummy `[*]` no terminal).

---

## O que o co-cristal 6U26 mudou

`6U26` = PCSK9 + "composto 16" (HET `063`), raios-X 1,53 Å. Cadeia A =
pró-domínio (61–152), cadeia B = catalítico + CHRD (153–682), ligante na B.

Análise da estrutura (`scripts/prep_pcsk9_receptor.py`):

| medida | valor |
|---|---|
| centro do sítio | `[38.59, 25.82, 26.44]` |
| box (núcleo enterrado + 4 Å) | `20,5 × 18,5 × 18,2 Å` |
| resíduos de contato (≤ 4,5 Å) | 25, todos na cadeia B |
| vetor de saída | `[0.222, 0.816, 0.533]`, braço exposto de 15,7 Å |
| distância mínima até a interface EGF(A)/LDLR | **23,4 Å** |

**O sítio não é a interface clássica com o LDLR.** Os contatos são Tyr293,
Cys323–Thr335, Arg357–Asp360, Arg412, Arg458–Ala463, Ile474–Ala478, Lys506 —
um bolsão a 23 Å da região Ser153/Arg194/Asp374 que a literatura de inibição
PCSK9–LDLR costuma citar. Se você tinha em mente ancorar a série na interface
EGF(A), este cristal aponta para outro lugar, e a decisão precisa ser
consciente. A lista completa de contatos fica em `pcsk9_site.json`.

Três coisas que o cristal resolve de graça:

1. **Onde dockar** — o bolsão é medido da estrutura, não estimado por detecção
   de cavidade.
2. **Como validar** — há ligante co-cristalizado, então o redocking com corte
   RMSD ≤ 2,0 Å do WP1 passa a ser aplicável. Sem ele, a etapa 3.0 seria docking
   sem controle.
3. **Para onde o linker sai** — o braço PEG/guanidina do `063` tem 17 átomos com
   *zero* vizinhos proteicos a 6 Å. É um vetor de saída medido
   experimentalmente, não inferido: exatamente onde o linker do PROTAC deve
   projetar.

---

## Dois bugs no notebook, um deles sério

**1. `redocking_rmsd()` (célula 1.6) não valida nada.** Ela usa
`rdMolAlign.GetBestRMS`, que **superpõe as moléculas antes de medir**. Teste
reproduzível em `scripts/rmsd_inplace.py`:

```
pose deslocada 25 Å do cristal:
  rmsd_inplace : 25.00   <- correto
  GetBestRMS   :  0.00   <- "VALIDADO" no bolsão errado
  CalcRMS      : 25.00   <- correto
```

Uma pose ancorada no bolsão errado passa pelo corte de 2,0 Å. **Isso afeta o
WP1 também**, onde a validação do protocolo de docking é o que autoriza toda a
triagem. Troque por `rdMolAlign.CalcRMS` ou pelo `rmsd_inplace()`.

**2. `GetBestRMS` altera a molécula probe**, deixando as coordenadas superpostas
à referência. Quem chamar e depois reaproveitar o objeto da pose (salvar em
disco, medir contatos, alimentar o PRosettaC) trabalha com a pose deslocada. Se
precisar chamá-la, passe `Chem.Mol(pose)`.

O `rmsd_inplace()` existe pelo que a `CalcRMS` não faz: restringir o RMSD a um
subconjunto de átomos. No `063` isso é necessário — exigir que o docking
reproduza um braço de 15 Å solto no solvente reprovaria um protocolo que acerta
o núcleo ancorado. A validação mede o núcleo e reporta o ligante inteiro à parte.

---

# Ordem de execução

Numerada, com env e tempo. Os passos 5 e 7 são os longos.

### 1. Gerar a biblioteca de warheads — `mdtools`, ~20 min

```bash
conda activate mdtools
cd ~/Protac_env
python scripts/generate_pcsk9_warheads.py --outdir ~/PCSK9_warheads --n-confs 50
```

180 moléculas únicas: R1 (8) × R2 (6) × R3 (4), deduplicadas por SMILES
canônico. `--no-3d` confere a enumeração em segundos; `--only-figure` restringe
aos 60 que usam só os substituintes nomeados na figura.

**Confira antes de seguir:** abra `pcsk9_warheads.csv` e olhe `filter_flags`
(139/180 na faixa típica) e `net_charge_ph74` (90/180 com carga ≠ 0).

### 2. PDBQT dos warheads — `pf_vs`, minutos

```bash
conda activate pf_vs
bash ~/PCSK9_warheads/prepare_pdbqt.sh
```

### 3. Inspecionar o cristal — `mdtools`, segundos

```bash
conda activate mdtools
python scripts/prep_pcsk9_receptor.py \
    --pdb-file ~/structures/6U26_1.pdb \
    --outdir ~/PCSK9_docking --inspect-only
```

Confirma cadeias e ligante. Só siga se bater com A=61–152, B=153–682, HET `063`.

### 4. Preparar receptor e sítio — `mdtools`, minutos

```bash
python scripts/prep_pcsk9_receptor.py \
    --pdb-file ~/structures/6U26_1.pdb \
    --outdir ~/PCSK9_docking \
    --target-chain B --keep-chains A B \
    --site-mode ligand --ref-ligand-resname 063 \
    --burial-threshold 20
```

`--keep-chains A B` mantém o pró-domínio: ele fica a 9,1 Å do ligante, perto o
bastante para compor a borda do bolsão. Removê-lo abriria uma face que na
proteína real não está aberta.

`--burial-threshold 20` não é cosmético. Com o ligante inteiro o box sai com
~28 Å de aresta e o docking vira busca cega; restrito ao núcleo ancorado fica em
~20 Å, tamanho certo para 230–420 Da.

**Inspeção obrigatória** — o box sobre a estrutura, antes de gastar horas:

```bash
chimerax ~/PCSK9_docking/receptor/*_receptor.pdb ~/PCSK9_docking/receptor/ref_ligand.sdf
```

### 5. Validar o protocolo — `pf_vs`, ~30–60 min

```bash
conda activate pf_vs
python scripts/dock_warheads_pcsk9.py \
    --site ~/PCSK9_docking/pcsk9_site.json \
    --warheads-pdbqt ~/PCSK9_warheads/pdbqt \
    --outdir ~/PCSK9_docking/docking --validate-only
```

Redocking do `063` em 5 seeds. **Portão de decisão**: se o melhor RMSD de núcleo
não fechar em ≤ 2,0 Å, pare. O script se recusa a seguir para a triagem
(`--skip-validation` existe, mas use só se for registrar a decisão por escrito).

Se reprovar, mexa nesta ordem: (a) protonação do receptor, (b) tamanho do box,
(c) `--exhaustiveness`. O `063` é grande e flexível — se o núcleo não reproduzir,
é sinal de problema no receptor ou no box, não de azar.

### 6. Triagem em dois estágios — `pf_vs`, horas, `nohup`

```bash
nohup python scripts/dock_warheads_pcsk9.py \
    --site ~/PCSK9_docking/pcsk9_site.json \
    --warheads-pdbqt ~/PCSK9_warheads/pdbqt \
    --warheads-csv ~/PCSK9_warheads/pcsk9_warheads.csv --series A \
    --outdir ~/PCSK9_docking/docking \
    > ~/PCSK9_docking/docking.log 2>&1 &
```

30 runs por warhead → corte por score → 100 runs nos sobreviventes.
**Retomável**: poses já gravadas são puladas, uma queda não custa o lote.

`--series A` são os 45 de R3 = H (N primário, conjugação limpa). A Série B
(`--series B`) só depois que a pose mostrar que o N-H do R3 não faz ligação de
hidrogênio essencial — se fizer, conjugar ali destrói a interação que o R3 foi
proposto para criar.

Saída: `docking/heads_pcsk9/<ID>_in_pcsk9.sdf` = os `HeadB` posicionados.

### 7. Conferir as poses — notebook, minutos

Célula 3.0d. Duas checagens que custam minutos e evitam dias perdidos:

- `check_pose_consistency()` — `std_score` alto entre runs independentes é pose
  instável. Score bom com dispersão alta é candidato a falso positivo.
- `check_exit_vector_accessible()` — o N de conjugação está exposto ao solvente
  e alinhado com o braço do `063`? Se ficar enterrado, não há por onde o linker
  sair e o candidato morre no WP2. Melhor descobrir agora.

### 8. Montar os PROTACs — notebook, minutos

Célula 3.2. `assemble_full_protac()` funciona **sem alteração**: o gerador
renumera cada warhead para que o átomo 0 seja o N de conjugação, que é o que
`ReplaceSubstructs(replacementConnectionPoint=0)` espera. Testado nas quatro
variantes de R3.

Vigie o tamanho: 10 sub-complexos × 30 warheads = 300 jobs PRosettaC de horas
cada. Corte na priorização, e pela geometria do linker do WP2, não só por
descritores.

### 9. PRosettaC e AlphaFold 3 — workstation, dias

`emit_prosettac_jobs()` gera um config por candidato. **Rode um sozinho
primeiro** e confira o `Anchor atoms`: nos SDFs deste pipeline o átomo de
conjugação é o primeiro, mas sua build pode contar de 0 ou de 1. Ajuste
`WARHEAD_ANCHOR_SERIAL` e só então libere o lote.

### 10. MD e análise — 3.4 a 3.8, sem alteração

- **`acpype -n`**: use `net_charge_ph74` somado à carga do recrutador–linker
  (`protac_net_charge()`). Carga errada contamina a MD silenciosamente.
- **Nível (ii)** (`E3_recruiter_linker_warhead`): detecta interação espúria
  warhead–E3 sem a PCSK9. Risco maior na série de naftaleno (R2 = `1-Naph` /
  `2-Naph`, vários marcados por cLogP). Rode este nível **antes** de gastar
  tempo com o ternário completo.

---

## Resumo dos arquivos

| arquivo | env | papel |
|---|---|---|
| `scripts/generate_pcsk9_warheads.py` | mdtools | enumera a série, SDFs com átomo 0 = N |
| `scripts/prep_pcsk9_receptor.py` | mdtools | receptor, sítio, vetor de saída, ligante de referência |
| `scripts/dock_warheads_pcsk9.py` | pf_vs | validação + triagem + `Heads` posicionados |
| `scripts/rmsd_inplace.py` | qualquer | RMSD correto para redocking (substitui `GetBestRMS`) |
| `notebooks/wp3_pcsk9_warheads_cells.py` | mdtools | células 3.0–3.3 para o notebook |

Tudo determinístico (`--seed`, default `0xC0FFEE`). Para a tese, registre as
linhas de comando exatas e a versão do RDKit.
