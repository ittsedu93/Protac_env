# Warheads PCSK9 (feniletilamina) na Fase 3 do `protac_wp1_wp3_pipeline.ipynb`

Guia de execução. Assume WP1 e WP2 concluídos (recrutadores VHL/CRBN com exit
vector definido; sub-complexos recrutador–linker validados, com cap dummy `[*]`
no terminal livre).

---

## Visão geral do que muda na Seção 3

A metodologia original entrava no WP3 com os **ligantes GMF**, que já tinham
pose conhecida na PCSK9. Os warheads da Figura 3 são moléculas novas, então a
Fase 3 ganha uma etapa que antes não existia:

| Etapa | Original (GMF) | Com warheads da Figura 3 |
|---|---|---|
| 3.1 | apontar `GMF_LIGANDS_DIR` | enumerar a série e priorizar (`scripts/generate_pcsk9_warheads.py`) |
| **3.0** | — (pose já existia) | **docking dos warheads na PCSK9** para produzir o `HeadB.sdf` posicionado |
| 3.2 | conjugar GMF ao linker | idêntico — `assemble_full_protac()` funciona sem alteração |
| 3.3–3.8 | inalterado | inalterado |

A etapa 3.0 é obrigatória: o PRosettaC exige `Heads` **já posicionados** dentro
das respectivas estruturas. Um confôrmero isolado do warhead não serve.

---

## Passo 1 — Gerar a biblioteca de warheads (roda em segundos/minutos)

No env `mdtools`:

```bash
conda activate mdtools
cd ~/Protac_env
python scripts/generate_pcsk9_warheads.py --outdir ~/PCSK9_warheads --n-confs 50
```

Enumera **180 moléculas únicas**: R1 (8) × R2 (6) × R3 (4), deduplicadas por
SMILES canônico (R1 e R2 são posições meta equivalentes, então pares simétricos
colapsam). Com `--only-figure` ficam só os 60 que usam exclusivamente os
substituintes nomeados na figura, sem as variações F/Cl/OH dos quadros.

Tempo: ~20 min para os 180 com 50 confôrmeros. Use `--no-3d` para conferir só a
enumeração (segundos).

### Grupos enumerados

| Posição | Grupos |
|---|---|
| **R1** | `OPh`, `O-4F-Ph`, `O-4Cl-Ph` (-OAr), `cPr`, `cBu`, `F`, `Cl`, `OH` |
| **R2** | `CF3`, `2-Naph`, `1-Naph`, `F`, `Cl`, `OH` |
| **R3** | `H` (amina livre), `Ms` (metilsulfonil), `Tet` (1H-tetrazol-5-il), `CH2Tet` |

Para ampliar o -OAr, edite `R1_GROUPS` no topo do script — cada entrada é só um
par `rótulo -> (SMILES do substituinte, origem)`.

### O que sai em `~/PCSK9_warheads/`

| Arquivo | Uso |
|---|---|
| `pcsk9_warheads.csv` | manifesto: SMILES, R1/R2/R3, descritores, vetor de saída, carga em pH 7,4 |
| `pcsk9_warheads.smi` | entrada rápida para outras ferramentas |
| `pcsk9_warheads_3d.sdf` | todos os confôrmeros de menor energia num arquivo |
| `sdf/<id>.sdf` | **um por warhead, átomo 0 = N de conjugação** — base do `HeadB.sdf` |
| `sdf_dummy/<id>_dummy.sdf` | mesma molécula com `[*]` explícito no ponto de conjugação (conferência visual) |
| `prepare_pdbqt.sh` | gera os PDBQT no env `pf_vs` |

**A renumeração para átomo 0 é o detalhe que faz o resto encaixar**:
`assemble_full_protac()` no notebook chama
`ReplaceSubstructs(..., replacementConnectionPoint=0)`, que liga o linker ao
átomo de índice 0 da molécula de substituição. Testado com as quatro variantes
de R3 — todas montam e sanitizam.

### Colunas do manifesto que você vai usar

- `filter_flags` — `ok`, ou quais descritores saem da faixa típica de warhead
  (MW ≤ 400, HBD ≤ 3, HBA ≤ 6, RotB ≤ 8, cLogP ≤ 4,5). **Nada é descartado**;
  a priorização é explícita na célula 3.1a, para poder ser justificada.
- `net_charge_ph74` — estimativa: amina alifática protona (+1), tetrazol
  desprotona (−1), sulfonamida N-H neutra. Alimenta `acpype -n` e o `genion`.
  Os combos `R3=CH2Tet` saem **zwitteriônicos** (carga líquida 0, mas com dois
  centros carregados) — trate com atenção no preparo do sistema.
- `exit_vector_note` — distingue a série de N primário da de N já substituído.

---

## Passo 2 — PDBQT dos warheads (env `pf_vs`)

```bash
conda activate pf_vs
bash ~/PCSK9_warheads/prepare_pdbqt.sh
```

Se o `mk_prepare_ligand.py` não estiver no PATH, use
`python -m meeko.cli.mk_prepare_ligand` ou o caminho completo do env.

---

## Passo 3 — Substituir a célula 3.1 do notebook

Abra `notebooks/wp3_pcsk9_warheads_cells.py` no VS Code (formato `# %%`, cada
bloco vira célula) e copie:

1. **3.1** — substitui a célula de configuração original: troca
   `GMF_LIGANDS_DIR` por `WARHEADS_DIR` e acrescenta os PDB IDs da PCSK9.
2. **3.1a** — carrega o manifesto e define **Série A** (R3 = H, N primário,
   conjugação limpa) e **Série B** (R3 farmacofórico mantido, o linker ocupa o
   H restante e o N vira terciário).

Recomendação: leve a **Série A** como principal. A Série B só depois que a pose
em PCSK9 mostrar que o N-H não faz uma ligação de hidrogênio essencial — se
fizer, conjugar ali destrói a interação que você está tentando explorar.

---

## Passo 4 (pesado, workstation) — Etapa 3.0: docking dos warheads na PCSK9

Esta é a etapa nova. Antes de rodar:

1. **Escolher a estrutura e o sítio.** `2P4E` (PCSK9 madura) ou `3BPS`
   (PCSK9:LDLR-EGF(A)), conforme o bolsão que a série deve ocupar. Confira as
   cadeias no PDB baixado.
2. **Preencher `PCSK9_SITE_CENTER` e `PCSK9_BOX_SIZE`.** Sem ligante
   co-cristalizado não há redocking para validar o protocolo como no WP1 —
   declare o sítio por análise de cavidade (ChimeraX / fpocket) e **registre
   essa escolha por escrito**, porque ela é o ponto mais frágil da cadeia.
3. Rodar `staged_screening()` (30 runs → 100 runs) com o receptor PCSK9.
4. Converter a melhor pose de cada warhead priorizado em SDF:
   `obabel best_pose.pdbqt -O <id>_in_pcsk9.sdf`. **Este** é o `HeadB.sdf`.

Rode fora do kernel Jupyter (`nohup`/SLURM) — são horas.

---

## Passo 5 — 3.2: montar os PROTACs completos

A célula 3.2 do arquivo faz o produto cartesiano sub-complexo × warhead,
escreve um `protac.smi` por candidato e devolve a tabela `protac_df`.

Vigie o tamanho: 10 sub-complexos × 30 warheads = 300 candidatos, e cada um vira
um job PRosettaC de horas. Corte antes, na priorização do 3.1a — e priorize
também pela geometria do linker validada no WP2, não só pelos descritores.

---

## Passo 6 (pesado) — PRosettaC + AlphaFold 3

`emit_prosettac_jobs()` gera um `prosetta_config.txt` por candidato e um
`launch_all_prosettac.sh`.

**Rode um candidato sozinho primeiro** e confira o config gerado. Em especial o
`Anchor atoms`: nos SDFs deste pipeline o átomo de conjugação é o primeiro, mas
sua build do PRosettaC pode contar a partir de 0 ou de 1. Ajuste
`WARHEAD_ANCHOR_SERIAL` conforme o teste, aí libere o lote.

`emit_af3_jobs()` escreve os JSONs para submissão manual em
alphafoldserver.com (AF3 não roda local nesta máquina).

---

## Passo 7 — 3.4 a 3.8 sem alteração

Filtragem por confiança, checagens estruturais, MD nos três níveis e análise de
trajetória seguem como no notebook. Dois pontos de atenção:

- **`acpype -n`**: use `net_charge_ph74` do manifesto somado à carga do
  recrutador–linker (`protac_net_charge()` na célula 3.2b). Carga errada aqui
  contamina toda a MD silenciosamente.
- **Nível (ii) da MD** (`E3_recruiter_linker_warhead`): é onde se detecta
  interação espúria warhead–E3 na ausência da PCSK9. Com a série de naftaleno
  (R2 = `1-Naph`/`2-Naph`, os mais lipofílicos do conjunto — vários marcados em
  `filter_flags` por cLogP) esse risco é maior. Vale rodar esse nível antes de
  comprometer tempo com o ternário completo.

---

## Reprodutibilidade

O script é determinístico (`--seed`, default `0xC0FFEE`). Para a tese, registre
a linha de comando exata e a versão do RDKit — `--n-confs` muda qual
confôrmero é escolhido como representante.
