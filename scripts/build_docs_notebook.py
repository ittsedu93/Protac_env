#!/usr/bin/env python
"""
Monta notebooks/protac_pipeline_documentado.ipynb — o registro do pipeline.

Este notebook não *executa* o pipeline: ele conta o que cada fase faz, por que
faz assim, o que já se sabe dos resultados, e traz células que **leem** as
saídas do disco. Nenhuma célula escreve, roda docking ou toca na MD, então é
seguro abrir e executar enquanto a simulação roda na mesma máquina.

Gerado por script para ser reproduzível: para mudar o notebook, edite aqui e
rode de novo.

    python scripts/build_docs_notebook.py
"""

import json
from pathlib import Path

CELLS = []


def md(src: str):
    CELLS.append({"cell_type": "markdown", "metadata": {},
                  "source": src.strip("\n").splitlines(keepends=True)})


def code(src: str):
    CELLS.append({"cell_type": "code", "execution_count": None, "metadata": {},
                  "outputs": [],
                  "source": src.strip("\n").splitlines(keepends=True)})


# ===========================================================================
# Capa
# ===========================================================================
md(r"""
# PROTAC PCSK9 — registro do pipeline, fase por fase

Este notebook é o **caderno de laboratório computacional** do projeto: o que
foi feito em cada fase, por quê, com que código, o que os resultados mostraram,
e onde o caminho teve de ser corrigido.

> **Nenhuma célula aqui escreve arquivos, roda docking ou mexe na MD.**
> Todas só leem as saídas que já estão no disco. É seguro executar este
> notebook com a simulação rodando na mesma máquina.

## O alvo e a ideia

A PCSK9 se liga ao receptor de LDL (LDLR) e o manda para degradação lisossomal;
menos LDLR na superfície do hepatócito significa mais LDL circulante. Os
anticorpos aprovados (evolocumabe, alirocumabe) bloqueiam essa interação, mas
são injetáveis e caros. Um **PROTAC** ataca por outra via: em vez de bloquear a
PCSK9, recruta uma ligase E3 para ubiquitiná-la e destruí-la.

Um PROTAC tem três partes:

```
   [ recrutador E3 ] — [ linker ] — [ warhead ]
     VHL ou CRBN                     liga na PCSK9
```

A molécula é uma ponte: ela precisa segurar as duas proteínas ao mesmo tempo,
na geometria certa, por tempo suficiente para a transferência de ubiquitina. É
por isso que o pipeline não termina no docking — termina numa MD que pergunta
se a ponte **aguenta**.

## O mapa das fases

| Fase | O que faz | Custo típico | Script |
|---|---|---|---|
| 1 | enumera as warheads fenetilamina | segundos | `generate_pcsk9_warheads.py` |
| 2 | sítio na PCSK9, portão de validação, triagem por docking | horas (GPU) | `prep_pcsk9_receptor.py`, `dock_warheads_pcsk9.py` |
| 3 | análise da triagem e escolha das warheads | segundos | `analyze_warhead_docking.py` |
| 4 | WP1: revalida o redocking e escolhe o recrutador E3 | minutos | `wp1_select_recruiter.py` |
| 5 | WP2: linkers do Chemspace, geometria no exit vector | ~1 h | `wp2_linker_tools.py`, `wp2_build_subcomplexes.py` |
| 6 | WP3: monta os PROTACs completos | minutos | `wp3_assemble_protacs.py` |
| 7 | ranqueia os PROTACs e escolhe o candidato da MD | segundos | `rank_protacs.py` |
| 8 | MD do nível (ii): E3 + recrutador-linker-warhead | dias (GPU) | `md_prepare.py`, `md_run.sh`, `md_analyze.py` |

Tudo isso é encadeado por `scripts/run_pipeline.sh`, que grava um marcador
`.done_<n>` por fase concluída — relançar depois de uma queda custa só a fase
interrompida.

## Duas ideias que valem para todo o pipeline

**1. Portões, não filtros.** Em três pontos o pipeline *para* se um critério
não é atingido, em vez de seguir com números ruins: o redocking do WP1, a
validação do protocolo na PCSK9 e o veredito da MD. Um filtro descarta
moléculas; um portão descarta o **método**.

**2. Falha silenciosa é o inimigo.** Vários erros encontrados aqui não
levantavam exceção — devolviam um número plausível e errado. Um RMSD de 0,00 Å
para uma pose a 25 Å do lugar. Um lote inteiro de PROTACs montado sem warhead.
Por isso quase todo script verifica a própria saída antes de devolvê-la, e o
apêndice deste notebook lista cada caso.
""")

code(r"""
# Configuração: aponta para as saídas no disco. Só leitura.
from pathlib import Path
import json, os
import pandas as pd

HOME = Path.home()
WORK = HOME / "PRosettaC_runs" / "vhl_crbn_pcsk9_protac"
OUT = WORK / "pipeline"                 # $PIPELINE_OUT
WARHEADS = HOME / "PCSK9_warheads"
DOCKING = HOME / "PCSK9_docking"
MD = OUT / "md"

pd.set_option("display.width", 160, "display.max_columns", 40)

def existe(p):
    p = Path(p)
    return f"{'OK ' if p.exists() else '-- '} {p}"

for p in [WARHEADS / "pcsk9_warheads.csv",
          DOCKING / "validation.json",
          DOCKING / "round2_scores.csv",
          DOCKING / "warhead_ranking.csv",
          OUT / "wp1_recruiter.json",
          OUT / "linker_ranking.csv",
          OUT / "protac_candidates.csv",
          OUT / "protac_ranking.csv",
          OUT / "md_candidato.json",
          MD / "md_sistema.json",
          MD / "md_resultados.csv",
          MD / "md_veredito.json"]:
    print(existe(p))
""")

# ===========================================================================
md(r"""
---
# Fase 1 — as warheads

**Script:** `scripts/generate_pcsk9_warheads.py` · **saída:**
`~/PCSK9_warheads/pcsk9_warheads.csv` + um `.sdf` por warhead

## O que foi feito

O ponto de partida é o esqueleto **fenetilamina** dos inibidores de PCSK9 da
literatura, com três posições variáveis:

| Posição | Variações |
|---|---|
| **R1** | `-OAr` (fenoxi e análogos), ciclopropil, ciclobutil |
| **R2** | `-CF3`, naftaleno, F, Cl, OH |
| **R3** | metilsulfonil (`-Ms`), tetrazol (`-Tet`), H |

O produto cartesiano dessas listas, depois de remover duplicatas por SMILES
canônico, dá **180 warheads únicas**. Cada uma sai como `.sdf` com confôrmeros
gerados por ETKDGv3 e minimizados em MMFF94s.

## Duas decisões de implementação que importam

**A montagem é por `Chem.molzip`, não por concatenação de texto.** A primeira
versão montava as moléculas colando pedaços de SMILES. Isso funcionou para os
primeiros casos e quebrou nos anéis: os dígitos de fechamento de anel (`c1ccc1`)
de dois fragmentos diferentes colidem quando as strings são emendadas, e o
resultado é uma molécula com a conectividade errada — que o RDKit aceita sem
reclamar. Com `molzip` e mapas de átomos (`[*:1]`), a ligação é feita no grafo,
não no texto, e não existe colisão possível.

**O átomo 0 de cada warhead é o nitrogênio de conjugação.** A função
`reorder_attachment_first` renumera os átomos para que o N que vai receber o
linker seja sempre o índice 0. Isso é o que permite, lá na fase 6, montar o
PROTAC sem procurar o ponto de ligação por SMARTS em cada molécula.
""")

code(r"""
# As warheads enumeradas, com os grupos R decodificados
p = WARHEADS / "pcsk9_warheads.csv"
if p.exists():
    w = pd.read_csv(p)
    print(f"{len(w)} warheads\n")
    print(w.head(8).to_string(index=False))
    for col in ("R1", "R2", "R3"):
        if col in w:
            print(f"\n{col}:")
            print(w[col].value_counts().to_string())
else:
    print("fase 1 ainda não rodada")
""")

# ===========================================================================
md(r"""
---
# Fase 2 — o sítio na PCSK9, o portão de validação, e a triagem

**Scripts:** `prep_pcsk9_receptor.py`, `docking_engines.py`,
`dock_warheads_pcsk9.py` · **saída:** `~/PCSK9_docking/`

## 2a. O sítio

O receptor vem do co-cristal **6U26**, que traz o ligante de referência `063`
num bolso da PCSK9. O sítio de docking é definido pelo **núcleo enterrado**
desse ligante (`--site-mode ligand`, corte de enterramento em 20):

```
centro   = [38.589, 25.818, 26.443]
caixa    = [20.47, 18.46, 18.16] Å
exit vec = [0.222, 0.8162, 0.5333]
```

O *exit vector* é a direção em que o solvente está acessível a partir do sítio.
É por ali que o linker tem de sair para alcançar a E3 — e é o critério
geométrico central da fase 5.

> **Uma correção importante para a tese:** eu havia dito antes que este sítio
> ficava na interface com o LDLR. Não fica. Medindo, ele está a **23,4 Å da
> interface EGF(A) do LDLR**. A warhead aqui, portanto, **não** é um
> antagonista competitivo da ligação ao LDLR — ela é um *ponto de apoio* para
> a degradação. Isso é uma vantagem conceitual do PROTAC (não precisa bloquear
> nada, só segurar), mas precisa estar escrito assim, não como bloqueio.

## 2b. O portão: o protocolo sabe reproduzir o cristal?

Antes de docar 180 moléculas novas, o protocolo tem de reencontrar a pose
cristalográfica do próprio `063`. Critério: **RMSD ≤ 2,0 Å** em 5 sementes
diferentes.

A primeira tentativa deu **9,92 Å** e levou ~25 min por semente. O diagnóstico:
o script redocava o ligante **inteiro** (76 átomos) numa caixa desenhada para o
**núcleo enterrado**. 19 dos 76 átomos ficavam fora da caixa — a pose
cristalográfica era literalmente inalcançável, e o docking estava sendo
reprovado por uma pergunta impossível. Redocando o núcleo truncado, que é o que
a caixa comporta:

```
RMSD = 0.21 Å  (5/5 sementes)
```

O portão passa, e agora ele mede o que se propôs a medir.

## 2c. A triagem, em duas rodadas, na GPU

O motor é o **Uni-Dock** (já instalado no env `pf_vs` — descobri isso tarde,
depois de tentar clonar envs sem permissão de escrita). Ele docka em lote na
GPU via `--gpu_batch`, o que muda a escala do problema: o que o Vina faria em
dias sai em horas.

- **Rodada 1**, exaustividade baixa, todas as 180 warheads. Sobrevive quem
  pontua abaixo de `-7.0` kcal/mol.
- **Rodada 2**, exaustividade alta, 5 sementes, só nos sobreviventes.

Duas coisas foram ajustadas aqui por experiência, não por teoria:

1. **A semente é o laço externo**, não interno (`runs_by_seed`). Assim cada
   semente é *um* lote na GPU, em vez de um lote por ligante. É a diferença
   entre usar a placa e enfileirar nela.
2. **Progresso a cada semente, com `flush=True`**, e `python -u` no driver. Sem
   isso o buffer de blocos do Python segura a saída, `tail -f` mostra um arquivo
   parado, e a pessoa conclui que o job travou quando ele está rodando.
""")

code(r"""
# O portão de validação e o resultado da triagem
p = DOCKING / "validation.json"
if p.exists():
    v = json.loads(p.read_text())
    print("PORTÃO DE VALIDAÇÃO DO PROTOCOLO")
    for k, val in v.items():
        print(f"  {k}: {val}")
else:
    print("validation.json ausente")

p = DOCKING / "round2_scores.csv"
if p.exists():
    r2 = pd.read_csv(p)
    print(f"\nRodada 2: {len(r2)} warheads")
    print(r2.sort_values("best_score").head(10).to_string(index=False))
""")

# ===========================================================================
md(r"""
---
# Fase 3 — a análise que mudou a estratégia

**Script:** `scripts/analyze_warhead_docking.py` · **saída:**
`warhead_ranking.csv`

Esta fase é a mais importante do pipeline inteiro, e não porque calcula algo
difícil: porque ela **desmonta a leitura ingênua do docking**. Quatro achados:

## 1. O score bruto é uma medida de massa (r = −0,92)

A correlação entre score de docking e número de átomos pesados é **−0,92**.
Isto é, 85% da variação do score se explica só pelo tamanho da molécula.
Ranquear por score bruto é ranquear por peso molecular com passos extras.

O contorno é a **eficiência de ligação** (*ligand efficiency*):

$$\mathrm{LE} = \frac{-\,\mathrm{score}}{N_{\text{átomos pesados}}}$$

E as duas listas não se parecem:

| | top-20 por score | top-20 por LE |
|---|---|---|
| naftil (R2) | **100%** | 10% |

O naftil domina o ranking por score porque é grande, não porque encaixa bem. Se
o critério fosse o score bruto, o pipeline inteiro teria seguido com naftalenos
por um artefato de normalização.

## 2. R3 = H é inviável — 0%

Isto contradiz o que **eu mesmo havia recomendado** ("comece pela Série A").
O teste é geométrico: a warhead só serve se o **nitrogênio de conjugação ficar
exposto ao solvente** na pose docada, porque é ali que o linker se liga. Se o N
fica enterrado, a warhead pode ter score ótimo e ser inútil.

| R3 | fração com o N exposto |
|---|---|
| metilsulfonil (`-Ms`) | **48%** |
| tetrazol (`-Tet`) | intermediária |
| H | **0%** (só 14% têm o N acessível, nenhuma passa o conjunto de critérios) |

A explicação é simples depois de vista: sem o grupo volumoso em R3, a molécula
acomoda-se mais fundo no bolso e enterra justamente a amina.

## 3. R2 = CF3 não aparece em nenhuma das duas listas

Nem por score, nem por LE. É uma ausência informativa: vale registrar na tese
como variação testada e descartada por evidência, não por omissão.

## 4. Estabilidade de pose entre sementes

Uma warhead cuja pose muda de lugar entre as 5 sementes não tem um modo de
ligação — tem ruído. O desvio entre sementes entra como critério
(`ANALISE_SD_MAX=0.5`) junto ao alinhamento com o exit vector
(`ANALISE_COS_MIN=0.3`).

> **Honestidade sobre o corte:** `COS_MIN=0.3` foi escolhido sem base empírica;
> ele dá 30 warheads de 146. Afrouxar para 0.2 aumenta o conjunto. A escolha
> precisa aparecer na tese **como escolha**, não como constante natural.
""")

code(r"""
# Os quatro achados, reproduzidos a partir do ranking gravado
p = DOCKING / "warhead_ranking.csv"
if p.exists():
    r = pd.read_csv(p)
    print(f"{len(r)} warheads analisadas\n")

    if {"best_score", "heavy_atoms"} <= set(r.columns):
        print(f"1) correlação score x átomos pesados: "
              f"{r['best_score'].corr(r['heavy_atoms']):.3f}")

    if {"best_score", "ligand_efficiency"} <= set(r.columns):
        top_s = set(r.nsmallest(20, "best_score").index)
        top_le = set(r.nlargest(20, "ligand_efficiency").index)
        print(f"2) top-20 por score e por LE compartilham "
              f"{len(top_s & top_le)} moléculas de 20")

    col_ok = next((c for c in ("viavel", "passa", "aprovada") if c in r), None)
    if col_ok and "R3" in r:
        print("\n3) viabilidade por R3:")
        print((r.groupby("R3")[col_ok].mean() * 100).round(1)
              .to_string() + "  (% viáveis)")

    print("\nTop 10 por eficiência de ligação:")
    cols = [c for c in ("warhead_id", "R1", "R2", "R3", "best_score",
                        "heavy_atoms", "ligand_efficiency") if c in r]
    print(r.nlargest(10, "ligand_efficiency")[cols].to_string(index=False))
else:
    print("warhead_ranking.csv ausente")
""")

# ===========================================================================
md(r"""
---
# Fase 4 — WP1: o recrutador E3, e o portão que não media nada

**Script:** `scripts/wp1_select_recruiter.py` · **saída:**
`$PIPELINE_OUT/wp1_recruiter.json`

## O que faz

Para cada E3 candidata (**VHL** e **CRBN**), o script pega o recrutador do
co-cristal, revalida o protocolo por redocking, e escolhe o par
(E3, recrutador) que passa o portão com a melhor geometria. O JSON de saída
carrega o recrutador escolhido, o **exit point** (onde o linker sai do
recrutador) e a **direção de saída**.

## O bug que invalidava o portão

O portão do WP1 usava `GetBestRMS` do RDKit. Essa função **superpõe** as duas
moléculas antes de medir — e, de quebra, **modifica a molécula sonda** no
processo. O efeito: uma pose a 25 Å do lugar certo dava RMSD **0,00 Å**. O
portão aprovava qualquer coisa, sem nunca dar erro.

A função correta é `CalcRMS`, que mede *in place*. Mas `CalcRMS` não permite
restringir a um subconjunto de átomos, e era disso que eu precisava para medir
só o núcleo enterrado. Está resolvido em `scripts/rmsd_inplace.py`, que
implementa `rmsd_inplace()` e traz a armadilha documentada no cabeçalho para
não ser reintroduzida.

> Registro de honestidade: ao diagnosticar isso eu também acusei `CalcRMS` de
> estar errada. Não estava — meu teste tinha sido contaminado pela mutação que
> o `GetBestRMS` faz na sonda. O erro era um, não dois.

## A guarda de farmacóforo

O ponto de conjugação do recrutador tem de ser escolhido, e a escolha automática
inicial pegou o **NH da glutarimida** da CRBN. Esse NH é justamente o que faz as
três ligações de hidrogênio com o trio de triptofanos da CRBN: usá-lo como ponto
de ligação do linker destrói o reconhecimento que o recrutador existe para ter.

A correção é uma lista de SMARTS de farmacóforo (`FARMACOFOROS`) com exclusão
ativa: átomos que participam do reconhecimento da E3 não podem ser pontos de
conjugação, ponto.
""")

code(r"""
p = OUT / "wp1_recruiter.json"
if p.exists():
    w1 = json.loads(p.read_text())
    e = w1.get("escolhido", w1)
    print("RECRUTADOR ESCOLHIDO")
    for k in ("e3", "nome", "pdb", "recruiter_sdf", "rmsd_redocking_A",
              "exit_point", "exit_direction", "receptor_pdb"):
        if k in e:
            print(f"  {k}: {e[k]}")
    if "por_e3" in w1:
        print("\nComparação entre E3:")
        for nome, d in w1["por_e3"].items():
            rm = d.get("rmsd_redocking_A", d.get("rmsd"))
            print(f"  {nome}: RMSD do redocking = {rm}")
else:
    print("wp1_recruiter.json ausente")
""")

# ===========================================================================
md(r"""
---
# Fase 5 — WP2: os linkers e a geometria da ponte

**Scripts:** `wp2_linker_tools.py`, `wp2_build_subcomplexes.py` · **saída:**
`linker_ranking.csv`, `subcomplexes_manifest.json`

## A biblioteca

**Só Chemspace** — `Chemspace_PROTACs_linkers_SDF.sdf` e
`Chemspace_PROTACs_anchors_SDF.sdf`. O Enamine ficou fora por decisão sua, e o
`config/pipeline.conf` tem `LINKERS_EXTRA=""` reservado caso volte.

Filtros: 5–40 átomos pesados, no máximo 15 ligações rotáveis, e — o critério
que separa linker de reagente — **bifuncionalidade obrigatória**. Um linker com
um só ponto de ligação não é linker.

## Como a geometria é avaliada

Para cada linker, 50 confôrmeros × 12 rotações em torno do eixo de saída. Cada
confôrmero é **colocado no exit vector do recrutador** e pontuado por: alcance
(consegue atravessar a distância até a warhead?), clash com o receptor, e
alinhamento com a direção de saída.

O detalhe: a versão original de `exit_vector_score` media clash **sem colocar o
confôrmero no exit vector**. Media a molécula onde ela estava, na origem do seu
próprio sistema de coordenadas. Os números saíam, eram plausíveis, e não
correspondiam a nada. Corrigido em `place_conformer_at_exit` +
`score_conformer_at_exit`, que são funções separadas de propósito: colocar e
pontuar são passos distintos e testáveis.

## Três armadilhas de SMARTS

Encontrar amina terminal por SMARTS é mais difícil do que parece.

| Tentativa | O que quebra |
|---|---|
| `[NH2]` | casa aminas **não terminais**, aceitando linkers que não podem ser conjugados |
| `!$(N[!#6])` | exclui qualquer N com hidrogênio **explícito** — rejeitaria silenciosamente **todos** os linkers amino-terminados, que são a química dominante do catálogo Chemspace |
| **`[NX3;H2;!$(N[!#6;!#1])][CX4]`** | o que ficou: N trivalente, dois H, vizinhos só C ou H, ligado a carbono sp³ |

A segunda linha é o tipo de erro mais perigoso deste pipeline: ela não falharia,
ela devolveria zero linkers válidos e a conclusão seria "o catálogo não serve".

## Rótulos, não posições

Os pontos de ligação do linker são marcados como `[1*]` (lado do recrutador) e
`[2*]` (lado da warhead) por `label_attachment_points`. A montagem da fase 6
**exige** os dois rótulos e levanta exceção se faltar um. O motivo está na
próxima seção.
""")

code(r"""
p = OUT / "linker_ranking.csv"
if p.exists():
    lk = pd.read_csv(p)
    print(f"{len(lk)} linkers avaliados")
    cols = [c for c in ("linker_id", "smiles", "heavy_atoms", "n_rotb",
                        "reach_A", "clash", "cos_exit", "score") if c in lk]
    ordem = "score" if "score" in lk else cols[-1]
    print(lk.nlargest(10, ordem)[cols].to_string(index=False))
else:
    print("linker_ranking.csv ausente")

p = OUT / "subcomplexes_manifest.json"
if p.exists():
    sc = json.loads(p.read_text())
    print(f"\n{len(sc)} sub-complexos recrutador-linker montados")
""")

# ===========================================================================
md(r"""
---
# Fase 6 — WP3: montar os PROTACs

**Script:** `scripts/wp3_assemble_protacs.py` · **saída:**
`protac_candidates.csv`, um diretório por candidato

## O bug que teria estragado o lote inteiro, em silêncio

A primeira versão montava em duas etapas: `cap_free_terminus` e depois
`assemble_full_protac`, ambas com `ReplaceSubstructs`. E `ReplaceSubstructs`,
quando o padrão **não casa**, não levanta erro: devolve a molécula **inalterada**
dentro de uma tupla de um elemento. O código seguia adiante feliz.

O resultado seria um lote completo de "PROTACs" que são, de fato, apenas
recrutador + linker, **sem warhead nenhuma** — moléculas válidas, com SMILES
bonitos, que passariam pelo ranking e chegariam à MD. Dias de GPU para simular
a molécula errada.

A correção tem duas partes: os rótulos `[1*]`/`[2*]` da fase 5, e um
`assemble_protac` que **levanta exceção** se o `[2*]` não estiver presente. A
lição geral: quando uma API devolve "não fiz nada" com o mesmo tipo de
"fiz", a verificação tem de ser explícita.

## O outro problema: o PDBQT perde informação

A pose que sai do Uni-Dock está em **PDBQT**, um formato que descarta **ordens
de ligação e hidrogênios**. Reconstruir a molécula a partir dele produziu:

```
33 elétrons radicalares
Explicit valence for atom N, 6 ...        (75 falhas de montagem)
antechamber: number of electrons is odd (429)
```

A solução não é reconstruir: é **transplantar as coordenadas da pose para a
molécula preparada**, que já tem ordens de ligação e hidrogênios corretos.
Casar os átomos dos dois lados exigiu quatro estratégias em cascata; o que
funcionou nos dados reais foi `AdjustQueryProperties(makeBondsGeneric)` seguido
de alinhamento por MCS (`rdFMCS` com `CompareAny` e `ringMatchesRingOnly`).

## E a identificação do catálogo

Identificar o composto de catálogo **pelo índice** estava errado: o índice 87
tinha 22 átomos onde a pose tinha 37 — os arquivos não estavam na mesma ordem.
Passou a ser por **InChIKey do esqueleto**, que independe de ordem: índice 209.
""")

code(r"""
p = OUT / "protac_candidates.csv"
if p.exists():
    pc = pd.read_csv(p)
    print(f"{len(pc)} PROTACs montados")
    cols = [c for c in ("candidate_id", "warhead_id", "linker_id", "protac_mw",
                        "n_rotb", "formal_charge", "protac_smiles") if c in pc]
    print(pc[cols].head(10).to_string(index=False))
    if "protac_mw" in pc:
        print(f"\nMassa molar: {pc['protac_mw'].min():.0f} – "
              f"{pc['protac_mw'].max():.0f} Da "
              f"(mediana {pc['protac_mw'].median():.0f})")
else:
    print("protac_candidates.csv ausente")
""")

# ===========================================================================
md(r"""
---
# Fase 7 — ranquear e escolher **um** candidato

**Script:** `scripts/rank_protacs.py` · **saída:** `protac_ranking.csv`,
`md_candidato.json`

A MD de nível (ii) custa dias de GPU por candidato. Simular todos é inviável,
então a fase 7 existe para transformar a lista em uma **fila**: um escore
composto que combina a qualidade da warhead (eficiência de ligação, não score
bruto), a geometria do linker no exit vector, o desvio do recrutador em relação
à pose docada, e penalidades de flexibilidade e massa.

O `md_candidato.json` guarda o candidato de posto 1 — e a fila continua ali. Se
a MD reprovar o primeiro, basta mudar `MD_RANK=2` no `config/pipeline.conf` e
relançar: nada precisa ser recalculado.

## O candidato que está rodando

| | |
|---|---|
| **ID** | `SC0006__WH023` |
| warhead | **WH023** — R1 = OPh, R2 = OH, R3 = tetrazol |
| linker | `NCCOCCOCCO` (tri-etilenoglicol amino-terminado) |
| massa molar | 802,9 Da |
| ligações rotáveis | 19 |
| carga formal | 0 |
| escore composto | 0,814 |
| desvio do recrutador em relação ao docking | **0,000 Å** |

Vale notar o que esse último número diz: o recrutador na geometria de partida da
MD está exatamente onde o docking o colocou. A conformação que a MD vai testar
não foi inventada pelo gerador de confôrmeros — ela preserva a pose que passou
pelo portão de validação.

E vale notar o que o candidato **não** tem: 802,9 Da e 19 rotáveis estão bem
fora de Lipinski. Isso é normal e esperado em PROTACs, que ocupam o espaço
"beyond rule of five" — mas é um ponto a tratar explicitamente na tese, não a
esconder.
""")

code(r"""
p = OUT / "protac_ranking.csv"
if p.exists():
    rk = pd.read_csv(p)
    cols = [c for c in ("rank", "candidate_id", "warhead_id", "linker_id",
                        "score_composto", "protac_mw", "n_rotb") if c in rk]
    print("FILA DA MD (os 10 primeiros)")
    print(rk.head(10)[cols].to_string(index=False))

p = OUT / "md_candidato.json"
if p.exists():
    print("\nCANDIDATO ESCOLHIDO")
    for k, v in json.loads(p.read_text()).items():
        print(f"  {k}: {v}")
""")

# ===========================================================================
md(r"""
---
# Fase 8 — a MD, e a cadeia de sete obstáculos

**Scripts:** `md_prepare.py`, `fix_receptor_for_md.py`, `build_topology.py`,
`build_index.py`, `md_run.sh`, `md_analyze.py`

## O que a MD pergunta

O nível (ii) simula **E3 ligase + recrutador-linker-warhead** — o lado da ponte
que é possível montar sem depender de uma predição de complexo ternário. A
pergunta é concreta: ao longo de 200 ns, o recrutador **fica** no bolso da E3
enquanto o resto da molécula se mexe, ou o linker arrasta o recrutador para
fora?

**Protocolo:** AMBER ff14SB + TIP3P, GAFF2/AM1-BCC para o PROTAC (acpype),
caixa cúbica de 1,3 nm, neutralização, minimização (steep), NVT 0,1 ns,
NPT 5 ns com restrições de posição, e **200 ns × 3 réplicas** com sementes
diferentes. Termostato Nose-Hoover, barostato Parrinello-Rahman, cortes em
0,9 nm, LINCS em h-bonds.

Três réplicas não é luxo: é a única forma de distinguir "o complexo é estável"
de "aquela trajetória em particular foi sortuda".

## Os obstáculos, em ordem, e por que cada um era o que era

Esta é a parte mais útil deste registro, porque cada obstáculo aqui é um caso em
que o pipeline **parecia** funcionar.

### 1. `Incomplete ring in HIS68`

Cristais têm cadeias laterais parcialmente resolvidas — 15 resíduos, neste caso.
O `pdb2gmx` recusa a estrutura e está certo em recusar. `fix_receptor_for_md.py`
reconstrói as cadeias laterais no ChimeraX headless via `swapaa`, e verifica o
resultado: **15 → 0 incompletos**.

Pelo caminho, três detalhes de sintaxe do ChimeraX que custaram tempo:
`delete ~(/A,/B)` é inválido (é `~/A,B`); `:063` é ambíguo (é `::name="063"`);
e sem `--exit` o ChimeraX fica esperando em `cmd>` para sempre, o que num
pipeline em background parece um travamento.

Uma nota sobre `swapaa ... criteria highest`: não existe. Os critérios são
letras (`d`/`c`/`h`/`p`), e passar `highest` dá `Unknown criteria: 'i'` — a
mensagem cita a segunda letra, o que torna o erro mais difícil de ler do que
precisava.

### 2. Ligações peptídicas que não existem

O `pdb2gmx` avisou, e o aviso passou batido:

```
Long Bond (1247-1249 = 0.558481 nm)
Long Bond (3358-3360 = 1.39835 nm)
```

Uma ligação peptídica mede **0,133 nm**. Essas medem de 0,56 a 1,40 nm. O
cristal 4TZ4 da CRBN tem alças não resolvidas, e o `pdb2gmx`, vendo uma cadeia
contínua no arquivo, criou ligações entre os resíduos que **ladeiam** cada
lacuna — resíduos que estão a nanômetros de distância. Os índices diferem de 2,
que é o padrão de `C(i)–N(i+1)` com o `O` no meio: são ligações peptídicas
falsas.

Uma mola harmônica esticada dez vezes tem energia para desmontar o sistema. E o
sintoma **não se parece com a causa**:

```
Step 5  LINCS WARNING ... 539 540   0.9983  863453.1250   0.1090
CUDA error #700 (cudaErrorIllegalAddress): an illegal memory access
```

O erro de CUDA é consequência de a placa tocar coordenadas que viraram lixo —
na CPU o mesmo sistema estouraria com outra mensagem. Quem persegue o erro de
CUDA procura driver, versão, memória da GPU; nada disso tem a ver.

`split_chain_gaps.py` mede o `C–N` de cada par consecutivo e insere `TER` onde
a distância passa de 2,5 Å. O `pdb2gmx` então trata cada segmento como cadeia
própria, e a carga total não muda: cada corte soma um NH3+ (+1) e um COO- (−1).
O script também avisa quando uma lacuna cai a menos de 12 Å do ligante — ali os
términos carregados ficam perto do sítio, e o certo passa a ser modelar a alça
em vez de cortá-la.

Depois disso, o `md_run.sh` **verifica**: sobrando qualquer `Long Bond` no log
do `pdb2gmx`, ele para, porque uma ligação longa aqui vira estouro dez minutos
adiante e o erro de então não aponta para cá.

### 3. `[ atomtypes ]` depois de `[ moleculetype ]`

O `.itp` do acpype traz as duas seções no mesmo arquivo. O GROMACS exige que
**todo** `[ atomtypes ]` venha antes da primeira molécula, e o `topol.top` do
`pdb2gmx` já define a proteína logo após o campo de força. Um `#include` único,
em qualquer posição, viola a ordem:

```
Fatal error: Invalid order for directive atomtypes
```

`build_topology.py` separa o arquivo em dois — tipos de átomo e molécula — e
insere cada um no lugar certo. Depois **verifica** a ordem final expandindo os
includes.

Meta-erro divertido: meu próprio verificador reprovou meu próprio comentário,
porque ele continha o texto `[ moleculetype ]`. O GROMACS ignora tudo depois de
`;`, e o verificador passou a ignorar também.

### 4. A retomada guardava o arquivo errado

`md_run.sh` pulava o `pdb2gmx` se `topol.top` existisse. Mas o `pdb2gmx` **cria**
o `topol.top` antes de terminar: uma execução que morreu no meio (a do HIS68)
deixou o arquivo no disco, e a retomada concluiu que a etapa tinha dado certo.
Todos os passos seguintes reclamaram de "file does not exist", e o erro que
aparecia no log era o **último**, não o primeiro.

Duas correções: a guarda passou a olhar a saída **final** de cada etapa
(`complexo.gro`), e todo comando passa por `gmx_ok`, que confere o código de
saída **e** a existência do arquivo esperado.

### 5. O `make_ndx` que nunca poderia funcionar

Os `tc-grps` dos `.mdp` precisam nomear grupos que existam no índice. A versão
escrita à mão adivinhava os números (`1 | 13`, `name 22 Protein_PTC`) e supunha
que o resíduo do ligante se chamasse `PTC`.

Só que o acpype batiza a *moleculetype* como se pede em `-b`, mas propaga para o
`.gro` o resíduo que veio do PDB do RDKit — e o default do RDKit é **`UNL`**:

```
Group    12 (          Other) has   105 elements
Group    13 (            UNL) has   105 elements
```

Um grupo `Protein_PTC` nunca apareceria, por mais que o número estivesse certo.

`build_index.py` não adivinha nada: roda `make_ndx` só para **listar**, lê os
números pelos nomes, identifica o ligante como o menor grupo que não é proteína,
água nem íon (assim `UNL`, `PTC`, `LIG` ou `MOL` são achados igualmente), monta
as fusões com os números reais, e **verifica** o índice antes de aceitá-lo: os
dois grupos têm de existir e a soma deles tem de cobrir o sistema inteiro — um
átomo sem acoplamento térmico faz o `grompp` recusar a simulação. Se a
verificação falha, o índice pela metade é **apagado**, não deixado no disco para
a próxima retomada aceitar.

Dois cuidados a mais: os nomes dos grupos passaram a ser fixos
(`Protein_LIG`/`Water_and_ions`), então os `.mdp` não dependem mais do nome do
resíduo; e um `Water_and_ions` que o GROMACS já tenha criado é **reaproveitado**
em vez de duplicado, para o `grompp` não ter de escolher entre homônimos.

### 6. O nome do resíduo, resolvido nas duas pontas

`md_prepare.py` agora batiza o resíduo como `PTC` na única vez em que as
coordenadas são escritas, e `md_analyze.py` **descobre** o nome na trajetória
quando o pedido não está lá (tirando proteína, água e íons, o que sobra é o
PROTAC). A primeira correção conserta o futuro; a segunda faz a análise
funcionar no sistema que já está parametrizado como `UNL`, sem refazer os 27
minutos de acpype.

### 7. O offload total da GPU nem sempre é aceito

`-nb gpu -pme gpu -bonded gpu -update gpu` é o ideal numa Blackwell de 32 GB,
onde o sistema inteiro (111.809 átomos) cabe na placa. Mas o GROMACS **recusa**
`-update gpu` em algumas combinações desta metodologia — o termostato
Nose-Hoover e as restrições de posição do equilíbrio são os casos conhecidos — e
a recusa é um erro de configuração que mata a etapa em segundos.

Em vez de escolher de antemão, `mdrun_ok` desce uma escada: offload total →
`-nb gpu -pme gpu` → CPU. Mas **só** quando a mensagem do próprio GROMACS é
sobre offload. Um erro de física (LINCS, explosão) falha de uma vez, porque
repetir três vezes uma simulação que vai estourar de novo custaria horas por
nada.

### 8. Retomada de verdade na produção

200 ns × 3 réplicas é medido em dias. Havendo checkpoint, o `mdrun` retoma com
`-cpi` em vez de recomeçar. Sem isso, uma queda na 190ª nanosegundo custaria a
réplica inteira.

## Os critérios do veredito

`md_analyze.py` aplica critérios explícitos, em `CRITERIOS`:

| Critério | Corte | Por quê |
|---|---|---|
| RMSD da proteína (CA) | ≤ 3,5 Å | acima disso o sistema não equilibrou, e nada mais no resultado significa algo |
| RMSD do PROTAC | ≤ 5,0 Å | generoso de propósito: o PROTAC **é** flexível; o que não pode é o recrutador sair |
| contatos do recrutador | mantidos | é o que a MD existe para medir |

Reprovando, a fila da fase 7 já tem o próximo: `MD_RANK=2`.
""")

code(r"""
# Estado do sistema montado para a MD
p = MD / "md_sistema.json"
if p.exists():
    for k, v in json.loads(p.read_text()).items():
        print(f"  {k}: {v}")
else:
    print("md_sistema.json ausente")
""")

code(r"""
# Progresso da MD, sem tocar em nada: tamanho das trajetórias e último log
import subprocess
for rep in sorted(MD.glob("rep*")):
    xtc, log = rep / "prod.xtc", rep / "prod.log"
    if xtc.exists():
        print(f"{rep.name}: prod.xtc = {xtc.stat().st_size/1e6:.1f} MB"
              f"{'  [CONCLUÍDA]' if (rep/'prod.gro').exists() else ''}")
    if log.exists():
        # no .log do GROMACS o cabeçalho "Step  Time" vem numa linha e os
        # valores na seguinte; imprimir só o cabeçalho não diz nada
        ls = log.read_text(errors="ignore").splitlines()
        idx = [i for i, l in enumerate(ls) if l.split()[:1] == ["Step"]]
        if idx and idx[-1] + 1 < len(ls):
            passo, tempo = (ls[idx[-1] + 1].split() + ["?", "?"])[:2]
            print(f"   passo {passo} = {tempo} ps")
if not any(MD.glob("rep*")):
    print("nenhuma réplica iniciada ainda — a MD está no preparo/equilíbrio")
    for f in ("proteina.gro", "complexo.gro", "neutro.gro", "grupos.ndx",
              "em.gro", "nvt.gro", "npt.gro"):
        print(f"  {'OK' if (MD/f).exists() else '--'} {f}")
""")

code(r"""
# O veredito, quando a MD terminar
p = MD / "md_resultados.csv"
if p.exists():
    print(pd.read_csv(p).to_string(index=False))
p = MD / "md_veredito.json"
if p.exists():
    print("\nVEREDITO")
    print(json.dumps(json.loads(p.read_text()), indent=2, ensure_ascii=False))
else:
    print("\n(veredito ainda não emitido — a MD está rodando)")
""")

# ===========================================================================
md(r"""
---
# Apêndice A — catálogo dos erros

Cada linha aqui é um erro que **não levantava exceção**, ou que levantava a
exceção errada. Estão juntos de propósito: a lista é a parte transferível deste
trabalho.

| # | Onde | Sintoma | Causa | Correção |
|---|---|---|---|---|
| 1 | fase 1 | conectividade errada, sem erro | dígitos de fechamento de anel colidiam ao emendar SMILES | `Chem.molzip` com mapas de átomos |
| 2 | WP1 | RMSD 0,00 Å para pose a 25 Å | `GetBestRMS` **superpõe** e muta a sonda | `rmsd_inplace.py` |
| 3 | WP3 | lote sem warhead nenhuma | `ReplaceSubstructs` sem casamento devolve a molécula intacta | rótulos `[1*]`/`[2*]` + exceção |
| 4 | WP2 | linkers monofuncionais aceitos | bifuncionalidade não era exigida | `filter_linkers` |
| 5 | WP2 | `[NH2]` casava aminas não terminais | SMARTS frouxo | `[NX3;H2;!$(N[!#6;!#1])][CX4]` |
| 6 | WP2 | **zero** linkers válidos | `!$(N[!#6])` exclui N-H explícito | SMARTS acima |
| 7 | WP2 | clash medido fora do sítio | confôrmero não era colocado no exit vector | `place_conformer_at_exit` |
| 8 | fase 2 | RMSD 9,92 Å, ~25 min/semente | ligante inteiro numa caixa do núcleo | redocar o núcleo truncado → 0,21 Å |
| 9 | fase 2 | `tail -f` parecia travado | buffer de blocos do Python | progresso por semente + `flush=True` + `python -u` |
| 10 | WP3 | 33 radicais, valência 6 no N | PDBQT perde ordens de ligação e H | transplante de coordenadas + MCS |
| 11 | WP1 | conjugação no NH da glutarimida | ponto escolhido sem olhar o farmacóforo | lista `FARMACOFOROS` com exclusão |
| 12 | WP3 | composto errado do catálogo | identificação por índice | InChIKey do esqueleto |
| 13 | driver | `--dry-run` baixava arquivos e marcava fases | `DRY` não era respeitado | guardas em `mark_done` e nos downloads |
| 14 | MD | "file does not exist" em cascata | guarda de retomada em `topol.top`, criado antes do fim | guardar pela saída final + `gmx_ok` |
| 15 | MD | `Incomplete ring in HIS68` | cadeias laterais parciais no cristal | `fix_receptor_for_md.py` (15 → 0) |
| 16 | MD | `Invalid order for directive atomtypes` | acpype junta as duas seções | `build_topology.py` separa e verifica |
| 17 | MD | `make_ndx` falhava sempre | grupo é `UNL`, não `PTC`; números adivinhados | `build_index.py` descobre e verifica |
| 18 | MD | `-update gpu` recusado | Nose-Hoover / posres não suportados | escada de offload em `mdrun_ok` |
| 19 | MD | `make_ndx` "não listou grupos" | parser lido do log do genion; o make_ndx usa outro formato | os dois formatos |
| 20 | MD | offload recusado tratado como erro de física | classificador procurava `not supported`, mensagem dizia `does not support` | redações extras + escada própria da minimização |
| 21 | MD | `cudaErrorIllegalAddress` no NVT | ligações peptídicas falsas de até 1,4 nm nas lacunas do cristal | `split_chain_gaps.py` + verificação de `Long Bond` |

O padrão que atravessa a lista: **API que devolve "não fiz nada" com o mesmo
tipo de "fiz"**. `ReplaceSubstructs`, `GetBestRMS`, `MolFromMolFile` (que
levanta `OSError` em vez de devolver `None`), o `pdb2gmx` que cria o arquivo
antes de terminar. Em todos, a defesa foi a mesma: **verificar a saída, não
confiar na chamada**.

## Apêndice B — o que continua sendo julgamento humano

O pipeline é automático de ponta a ponta, mas quatro coisas não foram — e não
devem ser — automatizadas:

1. **Disparar os jobs do PRosettaC.** Rode **um** antes do lote e confira o
   `Anchor atoms`: a convenção 0-based *vs* 1-based varia por build, e um lote
   inteiro com o átomo errado é um lote perdido.
2. **Submeter os JSONs do AlphaFold 3** em `alphafoldserver.com`.
3. **Inspecionar as poses no ChimeraX** antes de comprometer dias de GPU.
4. **Escolher os cortes.** `ANALISE_COS_MIN=0.3`, `DOCK_ROUND1_CUTOFF=-7.0`,
   `ANALISE_SD_MAX=0.5` são decisões suas, não constantes da natureza. A tese
   fica mais forte se aparecerem como decisões, com o efeito de afrouxá-las
   registrado.

## Apêndice C — mapa dos arquivos

| Arquivo | Papel |
|---|---|
| `config/pipeline.conf` | **tudo** que muda de máquina ou de projeto |
| `scripts/run_pipeline.sh` | driver das 8 fases, marcadores `.done_<n>`, `--from/--only/--list/--dry-run` |
| `scripts/rmsd_inplace.py` | RMSD que mede sem superpor (a armadilha nº 2) |
| `scripts/docking_engines.py` | Vina e Uni-Dock atrás da mesma assinatura |
| `scripts/wp2_linker_tools.py` | filtros, geometria no exit vector, rótulos, montagem |
| `scripts/build_topology.py` | ordem das diretivas do GROMACS |
| `scripts/build_index.py` | grupos de acoplamento descobertos e verificados |
| `docs/RUNBOOK.md` | como rodar, fase por fase |
| `docs/WP3_warheads_PCSK9.md` | a química das warheads |
""")

# ===========================================================================
nb = {
    "cells": CELLS,
    "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python",
                       "name": "python3"},
        "language_info": {"name": "python", "version": "3.11"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}

out = Path(__file__).resolve().parent.parent / "notebooks" / \
    "protac_pipeline_documentado.ipynb"
out.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n")
n_md = sum(1 for c in CELLS if c["cell_type"] == "markdown")
print(f"{out}  ({len(CELLS)} células: {n_md} markdown, "
      f"{len(CELLS) - n_md} código)")
