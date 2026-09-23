# Runbook — WP1 → WP2 → WP3 do início ao fim

Ordem de execução com env, tempo e critério de parada de cada etapa. Os passos
marcados **[PORTÃO]** decidem se faz sentido seguir; os marcados **[PESADO]**
vão para a workstation com `nohup`/SLURM, nunca no kernel do Jupyter.

## Por onde entrar

**`notebooks/protac_wp1_wp3_pipeline_v2.ipynb`** é o notebook integrado: já traz
as células do original que continuam válidas, as correções do WP1/WP2 e a etapa
3.0, na ordem certa. Abra ele com o kernel `mdtools` e rode a célula de Setup —
ela detecta o repositório sozinha, não importa onde você clonou.

Os arquivos `notebooks/wp1_wp2_revised_cells.py` e
`notebooks/wp3_pcsk9_warheads_cells.py` continuam disponíveis para quem preferir
enxertar as células no notebook antigo, mas não são necessários se usar a v2.

Para reconstruir a v2 depois de editar: `python scripts/build_notebook.py`.

## Pré-requisitos do ambiente

### Não instale o Vina — use o Uni-Dock

O cabeçalho do notebook original já registrava: **Uni-Dock está instalado no
`pf_vs`**, é "motor compatível c/ scoring Vina, GPU", mesma `scoring=vina`.

Isso resolve de vez o problema de permissão: nesta workstation os envs vivem em
`/home/soberano/miniconda3/envs/`, que pertence ao usuário `soberano`. Rodando
como outro usuário você ativa os envs mas não instala neles
(`EnvironmentNotWritableError`), e o `pf_vs` nem pip tem. Com o Uni-Dock não há
o que instalar.

Confirme os motores disponíveis antes de qualquer coisa:

```bash
conda activate pf_vs
python ~/Protac_env/scripts/docking_engines.py
```

Deve imprimir `unidock OK — unidock: Uni-Dock vX.Y`. O `--engine unidock` é o
padrão do `dock_warheads_pcsk9.py`, então os comandos do runbook já usam ele.

Além de dispensar instalação, o Uni-Dock docka em **lote na GPU**: o script
inverte os laços e faz uma chamada por seed com todos os warheads, em vez de um
processo por molécula. Para 45 warheads × 30 runs, são 30 chamadas em vez de
1350.

**Nada disso roda em célula de notebook.** `conda activate` não muda o
interpretador de um kernel já no ar, e o IPython converte `conda ...` em
`%conda`, devolvendo `CondaError: Run 'conda init' before 'conda activate'`.
Terminal do VS Code, sempre.

### O cristal da PCSK9

Baixe direto, sem depender de onde você salvou:

```bash
mkdir -p ~/structures && cd ~/structures
wget https://files.rcsb.org/download/6U26.pdb
```

Confira que bate com o que o pipeline espera — cadeia A 61-152, cadeia B
153-682, um HET `063`:

```bash
conda activate mdtools
python ~/Protac_env/scripts/prep_pcsk9_receptor.py \
    --pdb-file ~/structures/6U26.pdb --outdir ~/PCSK9_docking --inspect-only
```

Se aparecer mais de uma cópia do `063`, você pegou a unidade assimétrica com
várias moléculas no cristal; use a assembly biológica (`6U26.pdb1`) ou mantenha
só um par de cadeias.

---

## Correções aplicadas ao pipeline

Quatro defeitos encontrados ao revisar o notebook. Os dois primeiros
invalidavam resultados sem dar erro.

| # | onde | efeito |
|---|---|---|
| 1 | `cap_free_terminus` + `assemble_full_protac` (2.5 / 3.2) | o cap consome o `[#0]`; a montagem não acha o ponto, o RDKit devolve a molécula **inalterada** sem exceção, e o lote vira recrutador-linker **sem warhead** |
| 2 | `redocking_rmsd` (1.6) | `GetBestRMS` **superpõe antes de medir**: pose a 25 Å do sítio correto passa como RMSD 0,00. Também **altera as coordenadas do probe** |
| 3 | `filter_linkers` (2.2) | aceita linker **monofuncional**; `[NH2]` casa amina não terminal |
| 4 | `exit_vector_score` (2.3) | mede clash e ângulo **sem nunca posicionar** o linker no exit vector |

Substituições: `scripts/rmsd_inplace.py`, `scripts/wp2_linker_tools.py`,
`scripts/revalidate_redocking.py`, e as células em
`notebooks/wp1_wp2_revised_cells.py`.

---

# Fase 0 — Reconferir o WP1

## 0.1 [PORTÃO] Revalidar o redocking de VHL e CRBN — `mdtools`, minutos

Para cada alvo que já passou pelo portão original:

```bash
conda activate mdtools
cd ~/Protac_env
python scripts/revalidate_redocking.py \
    --receptor-pdb ~/PRosettaC_runs/vhl_crbn_pcsk9_protac/prep/VHL_6GFZ/6GFZ_receptor.pdb \
    --ref-ligand   ~/PRosettaC_runs/vhl_crbn_pcsk9_protac/prep/VHL_6GFZ/6GFZ_ref_ligand.sdf \
    --poses        '~/PRosettaC_runs/.../redock*.pdbqt' \
    --label        VHL_6GFZ \
    --out          ~/PRosettaC_runs/.../wp1_revalidation_VHL_6GFZ.csv
```

Repita para CRBN. Consolide com `consolidar_revalidacao()`
(`notebooks/wp1_wp2_revised_cells.py`).

**Critério:** se nenhum alvo tiver pose ≤ 2,0 Å pelo RMSD correto, a triagem do
WP1 não tem protocolo validado e os recrutadores priorizados precisam ser
refeitos antes de qualquer coisa. Se passar, siga — os recrutadores continuam
válidos e você tem o número certo para a tese.

## 0.2 Trocar a função no notebook

Substitua `redocking_rmsd()` da célula 1.6 por `redocking_rmsd_correto()`.
Se algum código reaproveita o objeto da pose depois de medir RMSD, passe
`Chem.Mol(pose)` para a medida.

---

# Fase 1 — Warheads PCSK9

## 1.1 Enumerar a série — `mdtools`, ~20 min

```bash
conda activate mdtools
python scripts/generate_pcsk9_warheads.py --outdir ~/PCSK9_warheads --n-confs 50
```

180 moléculas únicas (R1 8 × R2 6 × R3 4, deduplicadas por SMILES canônico).
`--no-3d` confere a enumeração em segundos; `--only-figure` restringe aos 60
que usam só os substituintes nomeados na figura.

**Confira:** `pcsk9_warheads.csv`, colunas `filter_flags` (139/180 na faixa
típica) e `net_charge_ph74` (90/180 com carga ≠ 0).

## 1.2 PDBQT — `pf_vs`, minutos

```bash
conda activate pf_vs
bash ~/PCSK9_warheads/prepare_pdbqt.sh
```

---

# Fase 2 — Ancorar os warheads na PCSK9 (etapa 3.0, nova)

O PRosettaC exige `Heads` já posicionados. Ao contrário dos ligantes GMF, estes
warheads não têm pose conhecida — é o que esta fase produz.

## 2.1 Inspecionar o cristal — `mdtools`, segundos

```bash
python scripts/prep_pcsk9_receptor.py \
    --pdb-file ~/structures/6U26.pdb \
    --outdir ~/PCSK9_docking --inspect-only
```

Confirme: cadeia A = 61–152 (pró-domínio), cadeia B = 153–682, HET `063`.

## 2.2 Receptor e sítio — `mdtools`, minutos

```bash
python scripts/prep_pcsk9_receptor.py \
    --pdb-file ~/structures/6U26.pdb \
    --outdir ~/PCSK9_docking \
    --target-chain B --keep-chains A B \
    --site-mode ligand --ref-ligand-resname 063 \
    --burial-threshold 20
```

Valores esperados: centro `[38.59, 25.82, 26.44]`, box `20,5 × 18,5 × 18,2 Å`,
vetor de saída `[0.222, 0.816, 0.533]` (braço exposto de 15,7 Å), 25 resíduos de
contato.

`--keep-chains A B` mantém o pró-domínio: ele fica a 9,1 Å do ligante e compõe
a borda do bolsão. `--burial-threshold 20` restringe o box ao núcleo ancorado —
com o ligante inteiro a aresta vai a ~28 Å e o docking vira busca cega.

**Inspeção obrigatória:**

```bash
chimerax ~/PCSK9_docking/receptor/*_receptor.pdb ~/PCSK9_docking/receptor/ref_ligand.sdf
```

## 2.3 [PORTÃO] Validar o protocolo — `pf_vs`, ~30–60 min

```bash
conda activate pf_vs
python scripts/dock_warheads_pcsk9.py \
    --site ~/PCSK9_docking/pcsk9_site.json \
    --warheads-pdbqt ~/PCSK9_warheads/pdbqt \
    --outdir ~/PCSK9_docking/docking --validate-only
```

**Critério:** melhor RMSD de núcleo ≤ 2,0 Å. O script recusa seguir se reprovar.
Reprovando, mexa nesta ordem: (a) protonação do receptor, (b) tamanho do box,
(c) `--exhaustiveness`. O `063` é grande e flexível — se o núcleo não
reproduz, o problema é receptor ou box, não azar.

## 2.4 [PESADO] Triagem em dois estágios — `pf_vs`, horas

```bash
nohup python -u scripts/dock_warheads_pcsk9.py \
    --site ~/PCSK9_docking/pcsk9_site.json \
    --warheads-pdbqt ~/PCSK9_warheads/pdbqt \
    --warheads-csv ~/PCSK9_warheads/pcsk9_warheads.csv --series A \
    --outdir ~/PCSK9_docking/docking \
    > ~/PCSK9_docking/docking.log 2>&1 &
```

30 runs → corte → 100 runs nos sobreviventes. **Retomável**: poses já gravadas
são puladas.

`--series A` = os 45 de R3 = H (N primário, conjugação limpa). A Série B só
depois que a pose mostrar que o N-H do R3 não faz ligação de hidrogênio
essencial — se fizer, conjugar ali destrói a interação que o R3 cria.

Saída: `docking/heads_pcsk9/<ID>_in_pcsk9.sdf` = os `HeadB` posicionados.

## 2.5 Conferir as poses — notebook, minutos

Célula 3.0d: `check_pose_consistency()` (desvio alto entre runs = pose
instável) e `check_exit_vector_accessible()` (o N de conjugação está exposto e
alinhado com o braço do `063`? se estiver enterrado, o candidato morre no WP2).

---

# Fase 3 — WP2 revisado

## 3.1 Filtrar linkers — `mdtools`, minutos

Células 2.1–2.2 de `notebooks/wp1_wp2_revised_cells.py`:

```python
linkers = load_library(ENAMINE_LINKER_DIR)
linkers_ok = filter_linkers(linkers, **LINKER_FILTERS)
```

Agora exige **bifuncionalidade** com separação topológica mínima. O relatório
diz quantos caíram por monofuncionalidade — se for a maioria, o subconjunto
exportado do catálogo Enamine provavelmente veio errado.

## 3.2 Avaliar a geometria no exit vector — `mdtools`, minutos a ~1 h

```python
df = avaliar_linkers(linkers_ok, exit_point, exit_direction, receptor_pdb)
```

`exit_point` e `exit_direction` vêm do WP1 (átomo do exit vector do recrutador
na pose validada). Cada confôrmero é **posicionado** no exit vector e a rotação
em torno do eixo é amostrada antes de medir clash.

**Critério:** priorize por `frac_sem_clash`, não pelo melhor confôrmero isolado.
Um linker que encaixa em 1 de 50 confôrmeros está geometricamente forçado,
mesmo que esse confôrmero pontue bem.

## 3.3 [PESADO] Docking restrito — `pf_vs`, horas (opcional)

Célula 2.4 original, box pequeno centrado no exit vector. A rationale já diz
que aqui o critério é adequação espacial, não afinidade — trate o score como
secundário ao `frac_sem_clash`.

## 3.4 Montar os sub-complexos — `mdtools`, minutos

```python
info = preparar_subcomplexo(linker, recruiter_mol, WP2_DIR, sc_id)
```

Rotula as pontas (`[1*]` recrutador, `[2*]` warhead), conjuga o recrutador e
grava **duas** versões:

- `<id>.sdf` → para o WP3, com o `[2*]` **intacto**
- `<id>_capped_md.sdf` → cópia capeada, só para a MD de nível (i)

Nunca capeie o objeto que segue para o WP3.

## 3.5 [PORTÃO] Validar os sub-complexos — segundos

```python
validar_subcomplexos()
```

Confere que cada SDF tem `[2*]` e não tem `[1*]` sobrando. **É o que impede o
lote de PROTACs sem warhead.** Não siga com nenhum que não passe.

## 3.6 [PESADO] Boltz-2 — horas

Use `filter_boltz_models_estrito()`, que exige o nome da coluna de confiança em
vez de cair num fallback silencioso para a última coluna.

---

# Fase 4 — WP3

## 4.1 Montar os PROTACs — notebook, minutos

Célula 3.2. `build_protac_matrix()` usa `assemble_protac()`, que **levanta
exceção** se o ponto `[2*]` não existir. Carregue os sub-complexos com
`rglob("*.sdf")` filtrando `"_capped_md" not in p.name`.

Vigie o tamanho: 10 sub-complexos × 30 warheads = 300 jobs PRosettaC de horas
cada. Corte na priorização.

## 4.2 [PESADO] PRosettaC — dias

`emit_prosettac_jobs()` gera um config por candidato. **Rode um sozinho
primeiro** e confira o `Anchor atoms`: nos SDFs deste pipeline o átomo de
conjugação é o primeiro, mas sua build pode contar de 0 ou de 1. Ajuste
`WARHEAD_ANCHOR_SERIAL` e só então libere o lote.

## 4.3 AlphaFold 3 — submissão manual

`emit_af3_jobs()` escreve os JSONs para alphafoldserver.com.

## 4.4 [PESADO] MD — dias

- **`acpype -n`**: use `net_charge_ph74` somado à carga do recrutador–linker
  (`protac_net_charge()`). Carga errada contamina a MD silenciosamente.
- **Nível (i)** usa os `*_capped_md.sdf`; níveis (ii) e (iii) usam o PROTAC
  montado.
- **Rode o nível (ii) antes do ternário completo**: é onde aparece interação
  espúria warhead–E3 sem a PCSK9. Risco maior na série de naftaleno (R2 =
  `1-Naph`/`2-Naph`, vários marcados por cLogP).

---

## Mapa dos arquivos

| arquivo | env | papel |
|---|---|---|
| `scripts/rmsd_inplace.py` | qualquer | RMSD correto para redocking |
| `scripts/revalidate_redocking.py` | mdtools | reconfere o portão do WP1 |
| `scripts/generate_pcsk9_warheads.py` | mdtools | enumera a série de warheads |
| `scripts/prep_pcsk9_receptor.py` | mdtools | receptor, sítio, vetor de saída |
| `scripts/dock_warheads_pcsk9.py` | pf_vs | validação + triagem + `Heads` |
| `scripts/wp2_linker_tools.py` | mdtools | linkers, colocação, montagem verificada |
| `notebooks/wp1_wp2_revised_cells.py` | mdtools | células do WP1/WP2 |
| `notebooks/wp3_pcsk9_warheads_cells.py` | mdtools | células 3.0–3.3 |

Todos determinísticos (`--seed`, default `0xC0FFEE`). Para a tese, registre as
linhas de comando exatas e a versão do RDKit.

## Autotestes

Cada módulo roda sozinho e demonstra o defeito que corrige:

```bash
python scripts/rmsd_inplace.py       # GetBestRMS dá 0,00 para pose a 25 Å
python scripts/wp2_linker_tools.py   # filtro, montagem verificada, colocação
```
