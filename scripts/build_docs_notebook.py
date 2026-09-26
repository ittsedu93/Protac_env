#!/usr/bin/env python
"""
Monta notebooks/protac_pipeline_documentado.ipynb — o registro do pipeline.

Este notebook não *executa* o pipeline: ele explica o que cada fase faz, por
que faz assim, o que os resultados mostraram, e traz células que **leem** as
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
# CAPA
# ===========================================================================
md(r"""
# PROTAC contra a PCSK9 — registro completo do pipeline

**Caderno de laboratório computacional.** Oito fases, do esqueleto químico da
warhead até a dinâmica molecular que testa se a molécula aguenta. Cada seção
diz o que a fase faz, por que faz assim, o que entra, o que sai, e o que os
resultados mostraram.

> **Nenhuma célula deste caderno escreve arquivos, roda docking ou mexe na MD.**
> Todas só leem o que já está no disco. É seguro abrir e executar com a
> simulação rodando na mesma máquina.

---

## 1. O problema

A **PCSK9** é uma protease que se liga ao receptor de LDL (LDLR) na superfície
do hepatócito e o arrasta para degradação lisossomal. Menos LDLR na superfície
significa menos captação de LDL, e portanto mais colesterol LDL circulante.
Mutações que aumentam a atividade da PCSK9 causam hipercolesterolemia familiar;
mutações que a inativam dão LDL baixo por toda a vida, sem prejuízo aparente —
o que faz dela um alvo com validação genética humana, não apenas farmacológica.

Os fármacos aprovados são **anticorpos monoclonais** (evolocumabe,
alirocumabe): eficazes, mas injetáveis, caros, e limitados a bloquear a PCSK9
**extracelular**. Um inibidor de pequena molécula tradicional enfrenta um
problema de superfície: a interface PCSK9–LDLR é plana e extensa, o tipo de
sítio que não acomoda bem uma molécula pequena.

## 2. A ideia do PROTAC

Um **PROTAC** (*PROteolysis TArgeting Chimera*) não bloqueia — **destrói**. É
uma molécula com três partes:

```
      [ recrutador E3 ] ─── [ linker ] ─── [ warhead ]
        liga na ligase           ponte       liga na PCSK9
        (VHL ou CRBN)
```

Ela aproxima a proteína-alvo de uma **ligase E3 de ubiquitina**. A E3
ubiquitina a PCSK9, o proteassomo a degrada, e o PROTAC se solta para repetir o
ciclo. Três consequências mudam o jogo:

1. **Não precisa de sítio ativo.** Basta um ponto de apoio qualquer na
   superfície, porque a warhead não tem de inibir nada — só segurar.
2. **É catalítico.** Uma molécula degrada muitas cópias da proteína, então a
   ocupação necessária é muito menor que a de um inibidor competitivo.
3. **O efeito é a ausência da proteína**, não a redução de sua atividade.

O preço é geométrico e é o centro deste trabalho: o PROTAC precisa segurar
**as duas proteínas ao mesmo tempo**, na orientação certa, por tempo
suficiente para a transferência de ubiquitina. Uma molécula que liga bem nas
duas pontas e não consegue formar o complexo ternário produtivo é inútil. É
por isso que o pipeline não termina no docking — termina numa MD que pergunta
se a ponte **aguenta**.

## 3. Como o pipeline responde a isso

| Fase | Pergunta que responde | Custo | Script |
|---|---|---|---|
| **1** | Que warheads existem no espaço químico proposto? | segundos | `generate_pcsk9_warheads.py` |
| **2** | Alguma delas liga na PCSK9 — e o protocolo sabe medir isso? | horas (GPU) | `prep_pcsk9_receptor.py`, `dock_warheads_pcsk9.py` |
| **3** | Quais ligam **e** têm o ponto de conjugação exposto? | segundos | `analyze_warhead_docking.py` |
| **4** | Qual E3 usar, e por onde o linker sai do recrutador? | minutos | `wp1_select_recruiter.py` |
| **5** | Que linkers do catálogo alcançam a warhead a partir dali? | ~1 h | `wp2_linker_tools.py`, `wp2_build_subcomplexes.py` |
| **6** | As moléculas completas se montam corretamente? | minutos | `wp3_assemble_protacs.py` |
| **7** | Qual candidato merece os dias de GPU da MD? | segundos | `rank_protacs.py` |
| **8** | O complexo aguenta 200 ns? | dias (GPU) | `md_prepare.py`, `md_run.sh`, `md_analyze.py` |

O fluxo de dados, de ponta a ponta:

```
  esqueleto fenetilamina
        │  R1 × R2 × R3
        ▼
  [1] 180 warheads (.sdf)
        │  docking na PCSK9 (6U26)
        ▼
  [2] scores + poses           ──┐
        │                        │ portão: redocking do 063 ≤ 2,0 Å
  [3] warheads viáveis  ◄────────┘
        │  (liga E tem o N exposto)
        │
        │      cristal da E3 (4TZ4)
        │            │
        │      [4] recrutador + exit vector  ◄── portão: redocking
        │            │
        │      [5] linkers do Chemspace com geometria compatível
        │            │
        └──────► [6] PROTACs completos
                     │  escore composto
                     ▼
                [7] fila de candidatos ──► [8] MD do melhor ──► veredito
                                                  │                 │
                                             reprovou? ─────────────┘
                                             MD_RANK=2, sem recalcular nada
```

## 4. Duas ideias que valem para todo o pipeline

**Portões, não filtros.** Em três pontos o pipeline **para** se um critério não
é atingido, em vez de seguir com números ruins: o redocking do WP1, a validação
do protocolo na PCSK9, e o veredito da MD. Um filtro descarta moléculas; um
portão descarta o **método**. Se o protocolo não reproduz a pose
cristalográfica que ele deveria reproduzir, nenhum resultado dele significa
nada, e é melhor descobrir isso antes de gastar dias de GPU.

**Falha silenciosa é o inimigo.** Vários erros encontrados aqui não levantavam
exceção — devolviam um número plausível e errado. Um RMSD de 0,00 Å para uma
pose a 25 Å do lugar. Um lote inteiro de PROTACs montado sem warhead nenhuma.
Um índice de acoplamento procurando um grupo que nunca poderia existir. Por
isso quase todo script deste repositório **verifica a própria saída** antes de
devolvê-la, e o Apêndice A cataloga cada caso com o sintoma, a causa e a
correção.
""")

md(r"""
---
## Como usar este caderno

Execute a célula abaixo primeiro: ela aponta para as saídas no disco e mostra o
que já existe. Todas as outras células são independentes — rode a da fase que
quiser inspecionar.

Os caminhos vêm de `config/pipeline.conf`, que é o único lugar onde mora
qualquer coisa que muda de máquina ou de projeto.
""")

code(r"""
# Configuração e inventário. Só leitura.
from pathlib import Path
import json
import pandas as pd

HOME = Path.home()
WORK = HOME / "PRosettaC_runs" / "vhl_crbn_pcsk9_protac"
OUT = WORK / "pipeline"                 # $PIPELINE_OUT
WARHEADS = HOME / "PCSK9_warheads"
DOCKING = HOME / "PCSK9_docking"
MD = OUT / "md"

pd.set_option("display.width", 170, "display.max_columns", 40)

ESPERADOS = [
    ("1", WARHEADS / "pcsk9_warheads.csv",      "warheads enumeradas"),
    ("2", DOCKING / "validation.json",          "portão do protocolo"),
    ("2", DOCKING / "round2_scores.csv",         "scores da rodada 2"),
    ("3", DOCKING / "warhead_ranking.csv",       "análise e viabilidade"),
    ("4", OUT / "wp1_recruiter.json",            "recrutador E3 escolhido"),
    ("5", OUT / "linker_ranking.csv",            "linkers avaliados"),
    ("5", OUT / "subcomplexes_manifest.json",    "sub-complexos montados"),
    ("6", OUT / "protac_candidates.csv",         "PROTACs completos"),
    ("7", OUT / "protac_ranking.csv",            "fila da MD"),
    ("7", OUT / "md_candidato.json",             "candidato escolhido"),
    ("8", MD / "md_sistema.json",                "sistema montado"),
    ("8", MD / "md_resultados.csv",              "métricas por réplica"),
    ("8", MD / "md_veredito.json",               "veredito final"),
]

print(f"{'fase':>4}  {'':3}  arquivo")
print("-" * 78)
for fase, p, desc in ESPERADOS:
    marca = "OK " if p.exists() else "-- "
    print(f"{fase:>4}  {marca}  {desc:<28s} {p.name}")
print("-" * 78)
print("OK = existe no disco   -- = fase não rodada ou ainda em andamento")
""")

# ===========================================================================
# FASE 1
# ===========================================================================
md(r"""
---
# Fase 1 — Enumerar as warheads

**Script:** `scripts/generate_pcsk9_warheads.py` · **env:** `mdtools` ·
**custo:** segundos
**Saída:** `~/PCSK9_warheads/pcsk9_warheads.csv` + um `.sdf` por warhead

## O que esta fase faz

Transforma uma proposta de química em um conjunto concreto de moléculas 3D
prontas para docking. O ponto de partida é o esqueleto **fenetilamina** dos
inibidores de PCSK9 descritos na literatura, com três posições variáveis:

| Posição | Variações exploradas | Papel estrutural |
|---|---|---|
| **R1** | `-OAr` (fenoxi e análogos), ciclopropil, ciclobutil | preenche o bolso hidrofóbico |
| **R2** | `-CF3`, naftaleno, F, Cl, OH | modula eletrônica e volume |
| **R3** | metilsulfonil (`-Ms`), tetrazol (`-Tet`), H | grupo polar de ancoragem |

O produto cartesiano dessas listas, após remover duplicatas por SMILES
canônico, dá **180 warheads únicas**. Cada uma sai como `.sdf` com confôrmeros
gerados por **ETKDGv3** e minimizados em **MMFF94s**.

## Por que assim

**A montagem é por `Chem.molzip`, não por concatenação de texto.** A primeira
versão colava pedaços de SMILES. Funcionou nos casos simples e quebrou nos
anéis: os dígitos de fechamento de anel (`c1ccc1`) de dois fragmentos
diferentes **colidem** quando as strings são emendadas, e o resultado é uma
molécula com conectividade errada — que o RDKit aceita sem reclamar. Com
`molzip` e mapas de átomos (`[*:1]`), a ligação é feita no grafo, não no texto,
e não existe colisão possível.

**O átomo 0 de cada warhead é o nitrogênio de conjugação.** A função
`reorder_attachment_first` renumera os átomos para que o N que vai receber o
linker seja sempre o índice 0. Isso é o que permite, na fase 6, montar o PROTAC
sem ter de procurar o ponto de ligação por SMARTS em cada molécula — o ponto é
conhecido por construção.

## Como rodar

```bash
python scripts/generate_pcsk9_warheads.py --outdir ~/PCSK9_warheads --nconfs 50
```

## O que conferir

- **180 warheads** no CSV (menos que 8 × 6 × 4 = 192, porque duplicatas saem)
- a distribuição de R1/R2/R3 deve ser aproximadamente uniforme
- cada `.sdf` tem confôrmeros, não uma estrutura só
""")

code(r"""
# Fase 1 — as warheads enumeradas, com os grupos R decodificados
p = WARHEADS / "pcsk9_warheads.csv"
if p.exists():
    w = pd.read_csv(p)
    print(f"{len(w)} warheads únicas\n")
    print(w.head(8).to_string(index=False))
    for col in ("R1", "R2", "R3"):
        if col in w:
            print(f"\n{col}:")
            print(w[col].value_counts().to_string())
else:
    print("fase 1 ainda não rodada")
""")

# ===========================================================================
# FASE 2
# ===========================================================================
md(r"""
---
# Fase 2 — O sítio na PCSK9, o portão de validação, e a triagem

**Scripts:** `prep_pcsk9_receptor.py`, `docking_engines.py`,
`dock_warheads_pcsk9.py` · **envs:** `mdtools` → `pf_vs` · **custo:** horas (GPU)
**Saída:** `~/PCSK9_docking/` — `validation.json`, `round1_scores.csv`,
`round2_scores.csv`, `heads_manifest.json`

## O que esta fase faz

Três coisas, em ordem obrigatória: define **onde** docar, prova que o protocolo
**sabe** docar ali, e só então docar as 180 warheads.

### 2a. O sítio

O receptor vem do co-cristal **6U26**, que traz o ligante de referência `063`
num bolso da PCSK9. O sítio é definido pelo **núcleo enterrado** desse ligante
(`--site-mode ligand`, corte de enterramento 20):

```
centro       = [38.589, 25.818, 26.443]
caixa        = [20.47, 18.46, 18.16] Å
exit vector  = [0.222, 0.8162, 0.5333]
```

O **exit vector** é a direção em que o solvente está acessível a partir do
sítio. É por ali que o linker tem de sair para alcançar a E3, e é o critério
geométrico central da fase 5. Sem ele, o docking mede afinidade e ignora a
única coisa que importa para um PROTAC: se há saída.

> ### Uma correção importante para a tese
> Eu havia afirmado antes que este sítio ficava na interface com o LDLR.
> **Não fica.** Medido, ele está a **23,4 Å da interface EGF(A) do LDLR**.
> A warhead aqui, portanto, **não** é um antagonista competitivo da ligação
> ao LDLR — ela é um **ponto de apoio** para a degradação. Isso é uma
> vantagem conceitual do PROTAC (não precisa bloquear nada, só segurar), mas
> precisa estar escrito assim, e não como bloqueio.

### 2b. O portão: o protocolo sabe reproduzir o cristal?

Antes de docar 180 moléculas novas, o protocolo tem de reencontrar a pose
cristalográfica do próprio `063`. Critério: **RMSD ≤ 2,0 Å** em 5 sementes.

A primeira tentativa deu **9,92 Å** e levou ~25 min por semente. O diagnóstico:
o script redocava o ligante **inteiro** (76 átomos) numa caixa desenhada para o
**núcleo enterrado**. 19 dos 76 átomos ficavam fora da caixa — a pose
cristalográfica era literalmente inalcançável, e o protocolo estava sendo
reprovado por uma pergunta impossível. Redocando o núcleo truncado, que é o que
a caixa comporta:

```
RMSD = 0,21 Å  (5/5 sementes)
```

O portão passa, e agora ele mede o que se propôs a medir. A lição é geral:
**um portão que reprova pode estar errado sobre o que pergunta**, e vale checar
a pergunta antes de afrouxar o critério.

### 2c. A triagem, em duas rodadas, na GPU

O motor é o **Uni-Dock** (já instalado no env `pf_vs` — descobri isso tarde,
depois de tentar clonar envs sem permissão de escrita). Ele docka em lote na
GPU via `--gpu_batch`, o que muda a escala: o que o Vina faria em dias sai em
horas. `docking_engines.py` mantém Vina e Uni-Dock atrás da **mesma
assinatura**, então trocar de motor é um argumento, não uma reescrita.

- **Rodada 1**: exaustividade baixa, todas as 180. Sobrevive quem pontua abaixo
  de `-7.0` kcal/mol.
- **Rodada 2**: exaustividade alta, **5 sementes**, só nos sobreviventes.

Duas coisas ajustadas por experiência, não por teoria:

1. **A semente é o laço externo**, não interno (`runs_by_seed`). Assim cada
   semente é *um* lote na GPU, em vez de um lote por ligante. É a diferença
   entre usar a placa e enfileirar nela.
2. **Progresso a cada semente, com `flush=True`**, e `python -u` no driver. Sem
   isso o buffer de blocos do Python segura a saída, `tail -f` mostra um arquivo
   parado, e a pessoa conclui que o job travou quando ele está rodando.

## Como rodar

```bash
# no env mdtools — receptor e sítio
python scripts/prep_pcsk9_receptor.py --pdb ~/structures/6U26.pdb \
    --chain B --keep-chains "A B" --ref-ligand 063 --site-mode ligand

# no env pf_vs — portão e triagem
python -u scripts/dock_warheads_pcsk9.py --engine unidock --validate-protocol
python -u scripts/dock_warheads_pcsk9.py --engine unidock --series all
```

## O que conferir

- `validation.json` com RMSD ≤ 2,0 Å nas 5 sementes — **se falhar, pare aqui**
- quantas warheads passaram a rodada 1 (se passarem quase todas, o corte está
  frouxo; se passar quase nenhuma, reveja o sítio antes de afrouxar)
- desvio entre sementes na rodada 2: pose que muda de lugar não é pose
""")

code(r"""
# Fase 2 — o portão e o resultado da triagem
p = DOCKING / "validation.json"
if p.exists():
    v = json.loads(p.read_text())
    print("PORTÃO DE VALIDAÇÃO DO PROTOCOLO")
    for k, val in v.items():
        print(f"  {k}: {val}")
    passou = v.get("passou", v.get("ok"))
    if passou is not None:
        print(f"\n  -> {'PASSOU' if passou else 'REPROVOU — nada depois disto vale'}")
else:
    print("validation.json ausente")

for nome, rotulo in (("round1_scores.csv", "Rodada 1 (triagem grossa)"),
                     ("round2_scores.csv", "Rodada 2 (5 sementes)")):
    p = DOCKING / nome
    if p.exists():
        r = pd.read_csv(p)
        print(f"\n{rotulo}: {len(r)} warheads")
        if "best_score" in r:
            print(f"  melhor {r['best_score'].min():.2f} | "
                  f"mediana {r['best_score'].median():.2f} | "
                  f"pior {r['best_score'].max():.2f} kcal/mol")
            print(r.nsmallest(8, "best_score").to_string(index=False))
""")

# ===========================================================================
# FASE 3
# ===========================================================================
md(r"""
---
# Fase 3 — A análise que mudou a estratégia

**Script:** `scripts/analyze_warhead_docking.py` · **env:** `mdtools` ·
**custo:** segundos
**Saída:** `warhead_ranking.csv`

## O que esta fase faz

Esta é a fase mais importante do pipeline, e não porque calcula algo difícil:
porque ela **desmonta a leitura ingênua do docking**. Ela responde "quais
warheads seguem para virar PROTAC" com quatro critérios, e três deles não têm
nada a ver com afinidade.

## Achado 1 — o score bruto é uma medida de massa (r = −0,92)

A correlação entre score de docking e número de átomos pesados é **−0,92**.
Isto é, ~85% da variação do score se explica só pelo tamanho da molécula.
Ranquear por score bruto é ranquear por peso molecular com passos extras.

O contorno padrão é a **eficiência de ligação**:

$$\mathrm{LE} = \frac{-\,\mathrm{score}}{N_{\text{átomos pesados}}}$$

E as duas listas quase não se tocam:

| | top-20 por score | top-20 por LE |
|---|---|---|
| naftil (R2) | **100%** | 10% |

O naftil domina o ranking por score porque é grande, não porque encaixa bem. Se
o critério fosse o score bruto, o pipeline inteiro teria seguido com naftalenos
**por um artefato de normalização** — e essa é exatamente a série que o
`docs/RUNBOOK.md` marca como risco de interação espúria com a E3.

## Achado 2 — R3 = H é inviável: 0%

Isto **contradiz uma recomendação que eu mesmo havia dado** ("comece pela Série
A"). O teste é geométrico, não energético: a warhead só serve se o **nitrogênio
de conjugação ficar exposto ao solvente** na pose docada, porque é ali que o
linker se liga. N enterrado = warhead com score ótimo e utilidade zero.

| R3 | fração viável |
|---|---|
| metilsulfonil (`-Ms`) | **48%** |
| tetrazol (`-Tet`) | intermediária |
| **H** | **0%** — só 14% têm o N acessível, e nenhuma passa o conjunto de critérios |

A explicação é simples depois de vista: sem o grupo volumoso em R3, a molécula
acomoda-se mais fundo no bolso e enterra justamente a amina que precisa ficar
livre.

## Achado 3 — R2 = CF3 não aparece em nenhuma das duas listas

Nem por score, nem por LE. É uma **ausência informativa**: vale registrar na
tese como variação testada e descartada por evidência, não por omissão.

## Achado 4 — estabilidade de pose entre sementes

Uma warhead cuja pose muda de lugar entre as 5 sementes não tem um modo de
ligação: tem ruído. O desvio entre sementes entra como critério
(`ANALISE_SD_MAX=0.5`) junto ao alinhamento com o exit vector
(`ANALISE_COS_MIN=0.3`).

> **Honestidade sobre os cortes.** `COS_MIN=0.3` foi escolhido **sem base
> empírica**; ele dá 30 warheads de 146. Afrouxar para 0.2 aumenta o conjunto.
> Isso precisa aparecer na tese **como escolha**, com o efeito de afrouxá-la
> registrado — não como constante natural.

## Como rodar

```bash
python scripts/analyze_warhead_docking.py --docking ~/PCSK9_docking
```

## O que conferir

- a correlação score × átomos pesados (se não for fortemente negativa, algo
  está diferente do esperado e o raciocínio acima precisa ser revisto)
- a sobreposição entre os dois top-20
- as taxas de viabilidade por R3
""")

code(r"""
# Fase 3 — os quatro achados, reproduzidos a partir do ranking gravado
p = DOCKING / "warhead_ranking.csv"
if p.exists():
    r = pd.read_csv(p)
    print(f"{len(r)} warheads analisadas\n")

    if {"best_score", "heavy_atoms"} <= set(r.columns):
        rho = r["best_score"].corr(r["heavy_atoms"])
        print(f"1) correlação score x átomos pesados: {rho:.3f}"
              f"   (r² = {rho**2:.2f} da variação é só tamanho)")

    if {"best_score", "ligand_efficiency"} <= set(r.columns):
        top_s = set(r.nsmallest(20, "best_score").index)
        top_le = set(r.nlargest(20, "ligand_efficiency").index)
        print(f"2) os dois top-20 compartilham {len(top_s & top_le)} de 20 moléculas")
        if "R2" in r:
            for rot, grupo in (("por score", top_s), ("por LE", top_le)):
                sub = r.loc[sorted(grupo), "R2"].value_counts(normalize=True) * 100
                print(f"   R2 no top-20 {rot}: "
                      + ", ".join(f"{k} {v:.0f}%" for k, v in sub.items()))

    col_ok = next((c for c in ("viavel", "passa", "aprovada", "ok") if c in r), None)
    if col_ok and "R3" in r:
        print("\n3) viabilidade por R3 (% que passam todos os critérios):")
        print((r.groupby("R3")[col_ok].mean() * 100).round(1).to_string())
        print(f"   total viável: {int(r[col_ok].sum())} de {len(r)}")

    print("\nTop 10 por eficiência de ligação (o critério que vale):")
    cols = [c for c in ("warhead_id", "R1", "R2", "R3", "best_score",
                        "heavy_atoms", "ligand_efficiency") if c in r]
    print(r.nlargest(10, "ligand_efficiency")[cols].to_string(index=False))
else:
    print("warhead_ranking.csv ausente")
""")

# ===========================================================================
# FASE 4
# ===========================================================================
md(r"""
---
# Fase 4 — WP1: a ligase E3 e o portão que não media nada

**Script:** `scripts/wp1_select_recruiter.py` · **env:** `mdtools` ·
**custo:** minutos
**Saída:** `$PIPELINE_OUT/wp1_recruiter.json`

## O que esta fase faz

Escolhe **qual ligase E3 recrutar** e **por onde o linker sai do recrutador**.
Para cada E3 candidata — **VHL** e **CRBN** — o script pega o recrutador do
co-cristal, revalida o protocolo por redocking, identifica o ponto de
conjugação e a direção de saída, e escolhe o par (E3, recrutador) que passa o
portão com a melhor geometria.

**Resultado desta execução: CRBN**, do cristal **4TZ4**, com

```
exit point  = [-41.719, 60.152, -86.692]
exit vector = [ 0.3075,  0.8534,   0.4209]
```

A escolha importa muito mais que parece. VHL e CRBN têm geometrias de
apresentação completamente diferentes, e o comprimento de linker que funciona
para uma raramente funciona para a outra.

## O bug que invalidava o portão

O portão do WP1 usava `GetBestRMS` do RDKit. Essa função **superpõe** as duas
moléculas antes de medir — e, de quebra, **modifica a molécula sonda** no
processo. O efeito: uma pose a **25 Å** do lugar certo dava RMSD **0,00 Å**.
O portão aprovava qualquer coisa, sem nunca dar erro.

A função correta é `CalcRMS`, que mede *in place*. Mas `CalcRMS` não permite
restringir a um subconjunto de átomos, e era disso que eu precisava para medir
só o núcleo enterrado. Está resolvido em `scripts/rmsd_inplace.py`, que
implementa `rmsd_inplace()` e traz a armadilha documentada no cabeçalho para
que não seja reintroduzida. O módulo roda sozinho e **demonstra** o defeito:

```bash
python scripts/rmsd_inplace.py    # GetBestRMS dá 0,00 para pose a 25 Å
```

> Registro de honestidade: ao diagnosticar isso eu também acusei `CalcRMS` de
> estar errada. Não estava — meu teste tinha sido contaminado pela mutação que
> o `GetBestRMS` faz na sonda. O erro era um, não dois.

## A guarda de farmacóforo

O ponto de conjugação do recrutador tem de ser escolhido, e a escolha
automática inicial pegou o **NH da glutarimida** da CRBN. Esse NH é justamente
o que faz as três ligações de hidrogênio com o trio de triptofanos da CRBN:
usá-lo como ponto de ligação do linker **destrói o reconhecimento que o
recrutador existe para ter**. A molécula ficaria perfeita no papel e cega na
célula.

A correção é uma lista de SMARTS de farmacóforo (`FARMACOFOROS`) com exclusão
ativa: átomos que participam do reconhecimento da E3 não podem ser pontos de
conjugação, ponto.

## Como rodar

```bash
python scripts/wp1_select_recruiter.py --e3 "VHL CRBN" --burial 20 \
    --out $PIPELINE_OUT/wp1_recruiter.json
```

## O que conferir

- RMSD do redocking de **cada** E3, não só da escolhida
- que o exit point **não** seja um átomo de farmacóforo
- a direção de saída aponta para o solvente, não para dentro da proteína
""")

code(r"""
# Fase 4 — o recrutador escolhido e a comparação entre E3
p = OUT / "wp1_recruiter.json"
if p.exists():
    w1 = json.loads(p.read_text())
    e = w1.get("escolhido", w1)
    print("RECRUTADOR ESCOLHIDO")
    for k in ("e3", "nome", "pdb", "recruiter_sdf", "receptor_pdb",
              "rmsd_redocking_A", "exit_point", "exit_direction",
              "exit_atom_idx"):
        if k in e:
            print(f"  {k}: {e[k]}")
    if "por_e3" in w1:
        print("\nComparação entre E3 (o portão vale para as duas):")
        for nome, d in w1["por_e3"].items():
            rm = d.get("rmsd_redocking_A", d.get("rmsd", "?"))
            print(f"  {nome:<6s} RMSD do redocking = {rm}")
else:
    print("wp1_recruiter.json ausente")
""")

# ===========================================================================
# FASE 5
# ===========================================================================
md(r"""
---
# Fase 5 — WP2: os linkers e a geometria da ponte

**Scripts:** `wp2_linker_tools.py`, `wp2_build_subcomplexes.py` ·
**env:** `mdtools` · **custo:** ~1 h
**Saída:** `linker_ranking.csv`, `subcomplexes_manifest.json`

## O que esta fase faz

Encontra, num catálogo comercial, os linkers que **conseguem atravessar** a
distância entre o exit point do recrutador e a warhead, na direção certa, sem
colidir com a E3. É aqui que o PROTAC deixa de ser duas moléculas e passa a ser
uma ponte com geometria.

## A biblioteca

**Só Chemspace** — `Chemspace_PROTACs_linkers_SDF.sdf` e
`Chemspace_PROTACs_anchors_SDF.sdf`. O Enamine ficou fora por decisão sua, e
`config/pipeline.conf` tem `LINKERS_EXTRA=""` reservado caso volte.

Filtros aplicados:

| Critério | Valor | Por quê |
|---|---|---|
| átomos pesados | 5 – 40 | abaixo não alcança, acima a entropia mata |
| ligações rotáveis | ≤ 15 | flexibilidade excessiva custa entropia de ligação |
| **bifuncionalidade** | **obrigatória** | linker com um só ponto de ligação não é linker |

## Como a geometria é avaliada

Para cada linker: **50 confôrmeros × 12 rotações** em torno do eixo de saída.
Cada confôrmero é **colocado no exit vector do recrutador** e pontuado por
alcance (atravessa a distância até a warhead?), clash com o receptor, e
alinhamento com a direção de saída.

O detalhe que importa: a versão original de `exit_vector_score` media clash
**sem colocar o confôrmero no exit vector**. Media a molécula onde ela estava,
na origem do seu próprio sistema de coordenadas. Os números saíam, eram
plausíveis, e não correspondiam a nada. Corrigido em `place_conformer_at_exit`
+ `score_conformer_at_exit` — funções separadas de propósito: **colocar** e
**pontuar** são passos distintos e testáveis isoladamente.

## Três armadilhas de SMARTS

Encontrar amina terminal por SMARTS é mais difícil do que parece.

| Tentativa | O que quebra |
|---|---|
| `[NH2]` | casa aminas **não terminais** — aceita linkers que não podem ser conjugados ali |
| `!$(N[!#6])` | exclui qualquer N com hidrogênio **explícito**: rejeitaria silenciosamente **todos** os linkers amino-terminados, que são a química dominante do catálogo Chemspace |
| **`[NX3;H2;!$(N[!#6;!#1])][CX4]`** | o que ficou: N trivalente, dois H, vizinhos só C ou H, ligado a carbono sp³ |

A segunda linha é o tipo de erro mais perigoso deste pipeline: ela não
falharia. Devolveria **zero** linkers válidos, e a conclusão seria "o catálogo
não serve" — uma conclusão sobre o catálogo tirada de um bug no filtro.

## Rótulos, não posições

Os pontos de ligação do linker são marcados como `[1*]` (lado do recrutador) e
`[2*]` (lado da warhead) por `label_attachment_points`. A montagem da fase 6
**exige** os dois rótulos e levanta exceção se faltar um. O motivo está na
próxima seção.

## Como rodar

```bash
python scripts/wp2_build_subcomplexes.py \
    --linkers $WORK/chemspace_linkers/Chemspace_PROTACs_linkers_SDF.sdf \
    --anchors $WORK/chemspace_anchors/Chemspace_PROTACs_anchors_SDF.sdf \
    --recruiter $PIPELINE_OUT/wp1_recruiter.json --nconfs 50 --nspins 12
```

Autoteste que demonstra os defeitos corrigidos:

```bash
python scripts/wp2_linker_tools.py
```

## O que conferir

- quantos linkers sobraram depois do filtro (zero = suspeite do filtro, não do
  catálogo)
- a distribuição de alcance: se nenhum atravessa, o problema é o exit vector
- os sub-complexos montados têm o `[2*]` livre para a warhead
""")

code(r"""
# Fase 5 — linkers avaliados e sub-complexos montados
p = OUT / "linker_ranking.csv"
if p.exists():
    lk = pd.read_csv(p)
    print(f"{len(lk)} linkers avaliados")
    num = [c for c in ("heavy_atoms", "n_rotb", "reach_A", "clash",
                       "cos_exit", "score") if c in lk]
    if num:
        print("\nfaixas:")
        print(lk[num].describe().loc[["min", "50%", "max"]].round(3).to_string())
    cols = [c for c in ("linker_id", "smiles", "heavy_atoms", "n_rotb",
                        "reach_A", "clash", "cos_exit", "score") if c in lk]
    ordem = "score" if "score" in lk else cols[-1]
    print(f"\nTop 10 por {ordem}:")
    print(lk.nlargest(10, ordem)[cols].to_string(index=False))
else:
    print("linker_ranking.csv ausente")

p = OUT / "subcomplexes_manifest.json"
if p.exists():
    sc = json.loads(p.read_text())
    print(f"\n{len(sc)} sub-complexos recrutador-linker montados")
""")

# ===========================================================================
# FASE 6
# ===========================================================================
md(r"""
---
# Fase 6 — WP3: montar os PROTACs completos

**Script:** `scripts/wp3_assemble_protacs.py` · **env:** `mdtools` ·
**custo:** minutos
**Saída:** `protac_candidates.csv`, um diretório por candidato com
`protac.smi`, `prosetta_config.txt`, e os lançadores `launch_all.sh` /
`launch_one.sh`

## O que esta fase faz

Costura recrutador + linker + warhead numa única molécula, com coordenadas 3D
consistentes com as poses que passaram pelos portões, e emite os arquivos de
entrada para a predição de complexo ternário (PRosettaC) e para a MD.

## O bug que teria estragado o lote inteiro, em silêncio

A primeira versão montava em duas etapas: `cap_free_terminus` e depois
`assemble_full_protac`, ambas com `ReplaceSubstructs`. E `ReplaceSubstructs`,
quando o padrão **não casa**, não levanta erro: devolve a molécula
**inalterada** dentro de uma tupla de um elemento. O código seguia adiante
feliz.

O resultado seria um lote completo de "PROTACs" que são, de fato, apenas
**recrutador + linker, sem warhead nenhuma** — moléculas válidas, com SMILES
bonitos, que passariam pelo ranking e chegariam à MD. Dias de GPU para simular
a molécula errada, e um resultado publicável e falso.

A correção tem duas partes: os rótulos `[1*]`/`[2*]` da fase 5, e um
`assemble_protac` que **levanta exceção** se o `[2*]` não estiver presente. A
lição geral: **quando uma API devolve "não fiz nada" com o mesmo tipo de "fiz",
a verificação tem de ser explícita.**

## O outro problema: o PDBQT perde informação

A pose que sai do Uni-Dock está em **PDBQT**, formato que descarta **ordens de
ligação e hidrogênios**. Reconstruir a molécula a partir dele produziu:

```
33 elétrons radicalares
Explicit valence for atom N, 6 ...          (75 falhas de montagem)
antechamber: number of electrons is odd (429)
```

A solução não é reconstruir: é **transplantar as coordenadas da pose para a
molécula preparada**, que já tem ordens de ligação e hidrogênios corretos.
Casar os átomos dos dois lados exigiu quatro estratégias em cascata; o que
funcionou nos dados reais foi `AdjustQueryProperties(makeBondsGeneric)` seguido
de alinhamento por MCS (`rdFMCS` com `CompareAny` e `ringMatchesRingOnly`).

## E a identificação no catálogo

Identificar o composto de catálogo **pelo índice** estava errado: o índice 87
tinha 22 átomos onde a pose tinha 37 — os arquivos não estavam na mesma ordem.
Passou a ser por **InChIKey do esqueleto**, que independe de ordem: índice 209.

## Como rodar

```bash
python scripts/wp3_assemble_protacs.py --subcomplexes $PIPELINE_OUT \
    --warheads ~/PCSK9_docking --n-linkers 10 --n-warheads 15
```

## O que conferir

- **cada candidato tem warhead** (é literalmente o bug acima; confira a massa
  molar: recrutador + linker sozinhos ficam bem abaixo de 700 Da)
- zero radicais e valências válidas
- o `Anchor atoms` do `prosetta_config.txt` — a convenção 0-based *vs* 1-based
  varia por build, e um lote inteiro com o átomo errado é um lote perdido
""")

code(r"""
# Fase 6 — os PROTACs montados
p = OUT / "protac_candidates.csv"
if p.exists():
    pc = pd.read_csv(p)
    print(f"{len(pc)} PROTACs montados")
    cols = [c for c in ("candidate_id", "warhead_id", "linker_id", "protac_mw",
                        "n_rotb", "formal_charge") if c in pc]
    print(pc[cols].head(10).to_string(index=False))
    if "protac_mw" in pc:
        print(f"\nmassa molar: {pc['protac_mw'].min():.0f} – "
              f"{pc['protac_mw'].max():.0f} Da  "
              f"(mediana {pc['protac_mw'].median():.0f})")
        baixas = pc[pc["protac_mw"] < 600]
        if len(baixas):
            print(f"  [ATENÇÃO] {len(baixas)} candidato(s) abaixo de 600 Da — "
                  f"leves demais para ter as três partes; confira se a warhead "
                  f"entrou")
    if "n_rotb" in pc:
        print(f"ligações rotáveis: {pc['n_rotb'].min()} – {pc['n_rotb'].max()}")
else:
    print("protac_candidates.csv ausente")
""")

# ===========================================================================
# FASE 7
# ===========================================================================
md(r"""
---
# Fase 7 — Ranquear, e escolher **um**

**Script:** `scripts/rank_protacs.py` · **env:** `mdtools` · **custo:** segundos
**Saída:** `protac_ranking.csv`, `md_candidato.json`

## O que esta fase faz

A MD de nível (ii) custa **dias de GPU por candidato**. Simular todos é
inviável, então esta fase transforma a lista em uma **fila**: um escore
composto que combina

- a qualidade da warhead (**eficiência de ligação**, não score bruto)
- a geometria do linker no exit vector
- o desvio do recrutador em relação à pose docada
- penalidades de flexibilidade e massa

O `md_candidato.json` guarda o candidato de posto 1 — e a fila **continua ali**.
Se a MD reprovar o primeiro, basta mudar `MD_RANK=2` em
`config/pipeline.conf` e relançar: nada é recalculado.

## O candidato que está na MD

| | |
|---|---|
| **ID** | `SC0006__WH023` |
| warhead | **WH023** — R1 = OPh, R2 = OH, R3 = tetrazol |
| linker | `NCCOCCOCCO` (tri-etilenoglicol amino-terminado) |
| E3 | **CRBN** (cristal 4TZ4) |
| massa molar | **802,9 Da** |
| ligações rotáveis | 19 |
| carga formal | 0 |
| escore composto | 0,814 |
| átomos do recrutador travados | 29 |
| **desvio do recrutador vs. docking** | **0,000 Å** |

Aquele último número é o que mais importa: o recrutador, na geometria de
partida da MD, está **exatamente** onde o docking o colocou. A conformação que
a MD vai testar não foi inventada pelo gerador de confôrmeros — ela preserva a
pose que passou pelo portão de validação.

E vale registrar o que o candidato **não** tem: 802,9 Da e 19 rotáveis estão
bem fora de Lipinski. Isso é **normal e esperado** em PROTACs, que ocupam o
espaço *beyond rule of five* — mas é um ponto a tratar explicitamente na tese,
não a esconder. A permeabilidade celular é o gargalo conhecido da classe.

## Como rodar

```bash
python scripts/rank_protacs.py --candidates $PIPELINE_OUT/protac_candidates.csv \
    --out $PIPELINE_OUT/protac_ranking.csv
```

## O que conferir

- o desvio do recrutador do posto 1 (acima de 0,5 Å já merece atenção: a MD
  partiria de uma geometria que a triagem não validou)
- se os primeiros postos são todos do mesmo linker ou da mesma warhead, a
  diversidade da fila é baixa e o plano B se parece demais com o plano A
""")

code(r"""
# Fase 7 — a fila e o candidato escolhido
p = OUT / "protac_ranking.csv"
if p.exists():
    rk = pd.read_csv(p)
    cols = [c for c in ("rank", "candidate_id", "warhead_id", "linker_id",
                        "score_composto", "protac_mw", "n_rotb",
                        "desvio_recrutador_A") if c in rk]
    print("FILA DA MD — os 10 primeiros")
    print(rk.head(10)[cols].to_string(index=False))
    for col, rot in (("warhead_id", "warheads"), ("linker_id", "linkers")):
        if col in rk:
            n = rk.head(10)[col].nunique()
            print(f"  diversidade no top-10: {n} {rot} distintos")
else:
    print("protac_ranking.csv ausente")

p = OUT / "md_candidato.json"
if p.exists():
    print("\nCANDIDATO ESCOLHIDO (MD_RANK do pipeline.conf)")
    for k, v in json.loads(p.read_text()).items():
        print(f"  {k}: {v}")
""")

# ===========================================================================
# FASE 8
# ===========================================================================
md(r"""
---
# Fase 8 — A dinâmica molecular

**Scripts:** `md_prepare.py`, `fix_receptor_for_md.py`, `split_chain_gaps.py`,
`build_topology.py`, `build_index.py`, `md_run.sh`, `explain_lincs.py`,
`md_analyze.py` · **envs:** `mdtools` + GROMACS 2026.2 · **custo:** dias (GPU)

## O que esta fase pergunta

O **nível (ii)** simula **E3 ligase + recrutador-linker-warhead** — o lado da
ponte que é possível montar sem depender de uma predição de complexo ternário.
A pergunta é concreta:

> Ao longo de 200 ns, o recrutador **fica** no bolso da CRBN enquanto o resto
> da molécula se mexe, ou o linker arrasta o recrutador para fora?

É a pergunta que o docking não pode responder, porque o docking é estático e
o problema do PROTAC é dinâmico. Um recrutador que escapa em 50 ns não vai
sustentar a transferência de ubiquitina, por bom que seja o score.

**Três réplicas não são luxo:** é a única forma de distinguir "o complexo é
estável" de "aquela trajetória em particular foi sortuda".

## O protocolo

| | |
|---|---|
| campo de força da proteína | **AMBER ff14SB** |
| água | **TIP3P**, caixa cúbica com 1,3 nm de folga |
| ligante | **GAFF2 + AM1-BCC** via acpype (~27 min) |
| eletrostática | PME, cortes de 0,9 nm |
| vínculos | LINCS em h-bonds, `lincs-order 8`, `lincs-iter 2` |
| termostato (produção) | **Nose-Hoover**, τ = 1,0 ps, 310 K |
| barostato | **Parrinello-Rahman**, τ = 2,0 ps, 1 bar |
| produção | **200 ns × 3 réplicas**, dt = 2 fs |
| gravação | 1000 frames por réplica, um a cada 200 ps |

### O sistema, em números

```
CRBN 4TZ4, cadeia C     361 resíduos, 2908 átomos do cristal
depois do pdb2gmx       5829 átomos (com hidrogênios), carga +5
PROTAC                  105 átomos
complexo                5934 átomos
solvatado               111.819 átomos (35.295 águas)
neutralizado            111.809 átomos (5 CL substituindo águas)
grupos de acoplamento   Protein_LIG 5934  |  Water_and_ions 105.875
```

### O equilíbrio, em degraus

Isto **não** é o protocolo que eu escrevi primeiro. O primeiro estourava; o que
está abaixo é o que funciona, e a razão de cada degrau está no Apêndice A.

| etapa | dt | duração | termostato | restrições |
|---|---|---|---|---|
| `em` | — | `emtol 1000`, tudo flexível | — | — |
| `em2` | — | `emtol 100`, `emstep 0.001` | — | — |
| `warm` | 0,5 fs | 20 ps, partindo de 100 K | V-rescale | soluto preso |
| `nvt` | 1 fs | 100 ps | V-rescale | soluto preso |
| `npt` | 2 fs | 5 ns | V-rescale + P-R | soluto preso |
| `prod` | 2 fs | **200 ns × 3** | **Nose-Hoover** + P-R | livre |

O **Nose-Hoover** da metodologia fica onde importa: na produção. Ele reproduz o
ensemble canônico corretamente, mas **não é robusto longe do equilíbrio** —
oscila, e num sistema recém-solvatado a oscilação vira estouro. O equilíbrio
usa **V-rescale**, que é dissipativo e perdoa geometria ruim. A escolha do
termostato de equilíbrio é detalhe prático e vale declará-la como tal.

A referência das restrições de posição é **sempre `em2.gro`**, nunca a etapa
anterior: encadear referências deixaria o soluto derivar degrau a degrau e
chegar na produção longe da pose do docking.

### O receptor teve de ser consertado duas vezes

**Cadeias laterais incompletas.** Cristais têm resíduos parcialmente
resolvidos — 15 neste caso. O `pdb2gmx` recusa a estrutura
(`Incomplete ring in HIS68`) e está certo em recusar.
`fix_receptor_for_md.py` reconstrói no ChimeraX headless via `swapaa` e
**verifica**: 15 → 0 incompletos.

**Lacunas na cadeia.** O 4TZ4 tem alças não resolvidas, e o `pdb2gmx`, vendo
uma cadeia contínua no arquivo, criou **ligações peptídicas entre resíduos que
estão a nanômetros de distância**. `split_chain_gaps.py` mede o C–N de cada par
consecutivo e insere `TER` onde passa de 2,5 Å. As quatro lacunas encontradas
batem uma a uma com os quatro avisos `Long Bond` do GROMACS:

| lacuna | C–N medido | `Long Bond` do pdb2gmx | distância ao ligante |
|---|---|---|---|
| SER126 \| GLU132 | 5,58 Å | 0,558481 nm | 17,1 Å |
| THR172 \| GLY176 | 7,00 Å | 0,699799 nm | 16,3 Å |
| VAL213 \| SER220 | 8,15 Å | 0,814510 nm | 42,7 Å |
| ASP265 \| SER272 | 13,98 Å | 1,398350 nm | 29,0 Å |

**Para a tese:** 4 quebras de cadeia introduzidas nas lacunas do 4TZ4, todas a
**mais de 16 Å do sítio**, carga total preservada em **+5** (cada corte soma um
NH3+ e um COO−, que se cancelam). Se alguma lacuna caísse dentro do bolso, o
correto seria **modelar a alça**, não cortá-la — e o script avisa quando isso
acontece.

## Os critérios do veredito

`md_analyze.py` aplica critérios explícitos, em `CRITERIOS`:

### Antes de qualquer critério: desfazer a PBC

Ao cortar as 4 lacunas do cristal, a proteína deixou de ser **uma** molécula e
virou **cinco**. O GROMACS envolve cada molécula independentemente na caixa
periódica, então segmentos da mesma proteína aparecem em lados opostos da
caixa. Qualquer RMSD calculado sobre isso mede a **aresta da caixa**:

```
RMSD da proteína  30.00 Å      <- impossível para uma proteína enovelada
frações ancoradas  1.00        <- e incompatível com a linha de cima
```

A contradição entre as duas linhas é a assinatura do artefato, e é por isso que
`md_analyze.py` hoje **se recusa a emitir veredito** com RMSD acima de 10 Å.
`md_fix_pbc.sh` aplica a receita padrão — `-pbc cluster` para juntar os solutos
na mesma imagem, `-pbc nojump` para impedir saltos entre quadros — e guarda só
proteína + ligante (5.934 átomos em vez de 111.809).

### Os critérios do veredito

Os contatos são contados como **pares de átomos** a menos de 4,5 Å. Um fragmento
ancorado na superfície de uma proteína tem **centenas** desses pares, então
cortes absolutos não significam nada — a primeira versão usava 20 e 15, e 15
pares seria exigir que a warhead flutuasse no vácuo. Os critérios são
**relativos**: ao frame inicial, que é a geometria validada pelo docking, e ao
próprio recrutador, que dá a escala do que é "muito contato com a E3" nesta
molécula e neste sistema.

| Critério | Corte | Por quê |
|---|---|---|
| **RMSD do sítio** (CA a ≤12 Å do recrutador) | ≤ 3,5 Å | o bolso que segura o recrutador é o que o nível (ii) pergunta |
| RMSD do PROTAC | ≤ 5,0 Å | generoso de propósito: o PROTAC **é** flexível; o que não pode é o recrutador sair |
| ancoragem do recrutador | ≥ 0,5 × inicial | metade dos contatos de partida mantida |
| warhead / recrutador | ≤ 1,0 | a warhead não pode engajar a E3 tanto quanto o recrutador, que foi **desenhado** para isso |
| crescimento warhead–E3 | ≤ 2,0 × inicial | ela engajou mais do que engajava no início? |
| frações ancoradas | ≥ 0,7 | não escapou e voltou |

O quarto critério é o que o nível (ii) existe para detectar: se a warhead gruda
na **própria** E3 na ausência da PCSK9, a orientação do linker é improdutiva —
no ternário ela competiria com a ligação ao alvo.

### Por que o portão é o sítio e não a proteína inteira

O RMSD global **continua sendo medido e impresso**, e um valor alto tem de ser
declarado na tese. Ele só deixou de reprovar sozinho, e a razão saiu dos dados
deste projeto:

| | rep1 | rep2 | rep3 | média |
|---|---|---|---|---|
| RMSD global | 4,71 | 2,88 | 5,08 | **4,22 Å** |
| RMSD do núcleo | 4,19 | 2,71 | 4,59 | 3,83 Å |
| **RMSD do sítio** | 1,55 | 1,63 | 1,91 | **1,70 Å** |

O bolso fica em 1,5–2,2 Å **do primeiro bloco de 25 ns ao último, sem
tendência**, enquanto o global vai a 5,9 e 5,6 Å em duas réplicas. Os 4 Å de
diferença estão inteiramente fora do sítio. A CRBN é multidomínio e tem
dobradiça conhecida entre o domínio tipo-Lon e o de ligação à talidomida; um
RMSD global alto com o sítio rígido é giro de domínio, não desenovelamento. E a
dobradiça foi provavelmente amplificada por nós: as 4 alças cortadas são o tipo
de conexão que limita movimento interdomínio.

Reprovando, a fila da fase 7 já tem o próximo: `MD_RANK=2` e relançar. Mas
**atenção ao modo de falha**: se o que falha é a proteína, e não o PROTAC, o
próximo candidato usa o mesmo receptor e vai falhar igual — trocar de candidato
seria gastar 45 h para reencontrar o mesmo número.

## O resultado: SC0006__WH023

| Critério | Valor | Corte | |
|---|---|---|---|
| RMSD do sítio | **1,70 Å** | ≤ 3,5 | OK |
| RMSD do PROTAC | 2,87 Å | ≤ 5,0 | OK |
| ancoragem mantida | 0,77 × inicial | ≥ 0,5 | OK |
| warhead / recrutador | 0,65 | ≤ 1,0 | OK |
| crescimento warhead–E3 | 0,99 × inicial | ≤ 2,0 | OK |
| frações ancoradas | 1,00 | ≥ 0,7 | OK |

**Aprovado no nível (ii)**, com uma observação a declarar: movimento de domínio
do receptor, com RMSD global de 4,22 Å, presente em 2 das 3 réplicas.

O número mais informativo é o **0,99×**: a warhead termina os 200 ns com
exatamente o mesmo engajamento com a CRBN que tinha no início. O risco central
do nível (ii) — a warhead migrar para a própria E3 e competir com o alvo no
ternário — **não se materializou**.

Desempenho: 310, 329 e 326 ns/dia; ~15 h por réplica, ~46 h no total.

## Como rodar

O pipeline faz tudo sozinho. Para deixar rodando e desligar o notebook:

```bash
rm -f $PIPELINE_OUT/.done_8
setsid bash ~/Protac_env/scripts/run_pipeline.sh > ~/pipeline.log 2>&1 &
disown -a
```

Acompanhar:

```bash
bash ~/Protac_env/scripts/status.sh
```

> **Não apague arquivos para relançar.** Cada etapa se guarda pela própria
> saída, então o script sabe o que pular — mas ele **não** sabe distinguir um
> arquivo que você apagou de propósito de um que nunca existiu. Apagar as
> saídas do equilíbrio com a produção em andamento não derruba a réplica que
> já roda, mas derrubaria a próxima (hoje o `npt.gro` é recuperado de um
> `prod.tpr`, justamente por isso).
""")

code(r"""
# Fase 8 — o sistema montado para a MD
p = MD / "md_sistema.json"
if p.exists():
    info = json.loads(p.read_text())
    print("SISTEMA DA MD")
    for k, v in info.items():
        print(f"  {k}: {v}")
    d = info.get("desvio_do_docking_A")
    if d is not None:
        print(f"\n  -> recrutador a {d} Å da pose do docking"
              + ("  (a MD parte da geometria validada)" if d <= 0.5 else
                 "  [ATENÇÃO] acima de 0,5 Å"))
else:
    print("md_sistema.json ausente — a fase 8 ainda não preparou o sistema")
""")

code(r"""
# Fase 8 — progresso, sem tocar em nada
import re

etapas = ["protac.pdb", "complexo.pdb", "topol.top", "complexo.gro",
          "neutro.gro", "grupos.ndx", "em.gro", "em2.gro", "warm.gro",
          "nvt.gro", "npt.gro"]
print("PREPARO E EQUILÍBRIO")
for f in etapas:
    print(f"  {'[x]' if (MD / f).exists() else '[ ]'} {f}")

reps = sorted(MD.glob("rep*"))
if reps:
    print("\nPRODUÇÃO (200 ns por réplica)")
    for rep in reps:
        if (rep / "prod.gro").exists():
            print(f"  {rep.name}: CONCLUÍDA")
            continue
        log = rep / "prod.log"
        if log.exists():
            ls = log.read_text(errors="ignore").splitlines()
            idx = [i for i, l in enumerate(ls) if l.split()[:1] == ["Step"]]
            if idx and idx[-1] + 1 < len(ls):
                campos = ls[idx[-1] + 1].split()
                ps = float(campos[1]) if len(campos) > 1 else 0.0
                print(f"  {rep.name}: {ps/2000:.1f}% — {ps:.0f} ps de 200000 "
                      f"(faltam ~{200 - ps/1000:.1f} ns)")
            perf = [l for l in ls if l.startswith("Performance:")]
            if perf:
                print(f"      {perf[-1].strip()}   (o 1o numero e ns/dia)")
        else:
            print(f"  {rep.name}: iniciando")
else:
    print("\nnenhuma réplica iniciada — a MD está no preparo ou equilíbrio")
""")

code(r"""
# Fase 8 — métricas e veredito, quando a MD terminar
p = MD / "md_resultados.csv"
if p.exists():
    res = pd.read_csv(p)
    print("MÉTRICAS POR RÉPLICA")
    print(res.to_string(index=False))
    for col, corte, rot in (("rmsd_proteina_media", 3.5, "RMSD da proteína"),
                            ("rmsd_protac_media", 5.0, "RMSD do PROTAC")):
        if col in res:
            m = res[col].mean()
            print(f"  {rot}: média {m:.2f} Å  (corte {corte}) -> "
                  f"{'OK' if m <= corte else 'FORA'}")
else:
    print("md_resultados.csv ausente — a MD ainda não terminou")

p = MD / "md_veredito.json"
if p.exists():
    print("\nVEREDITO")
    print(json.dumps(json.loads(p.read_text()), indent=2, ensure_ascii=False))
else:
    print("\n(veredito ainda não emitido)")
    print("Reprovando: MD_RANK=2 em config/pipeline.conf e relançar —")
    print("a fila da fase 7 já tem o próximo candidato, nada é recalculado.")
""")

# ===========================================================================
# APÊNDICES
# ===========================================================================
md(r"""
---
# Apêndice A — Catálogo dos erros

Cada linha é um erro que **não levantava exceção**, ou que levantava a exceção
errada, ou cuja mensagem apontava para o lugar errado. Estão juntos de
propósito: esta lista é a parte transferível do trabalho.

## Nas fases 1–7

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

## Na fase 8 — a cadeia de dez obstáculos

| # | Sintoma | Causa | Correção |
|---|---|---|---|
| 14 | "file does not exist" em cascata | guarda de retomada em `topol.top`, que o `pdb2gmx` cria **antes** de terminar | guardar pela saída **final** + `gmx_ok` |
| 15 | `Incomplete ring in HIS68` | 15 cadeias laterais parciais no cristal | `fix_receptor_for_md.py` (15 → 0) |
| 16 | `Invalid order for directive atomtypes` | acpype junta `[atomtypes]` e `[moleculetype]` no mesmo arquivo | `build_topology.py` separa e **verifica** |
| 17 | `make_ndx` falhava sempre | grupo é `UNL` (default do RDKit), não `PTC`; números adivinhados | `build_index.py` descobre e verifica |
| 18 | `make_ndx` "não listou grupos" | parser aprendido no log do **genion**; o `make_ndx` usa outro formato | os dois formatos |
| 19 | offload recusado tratado como erro de física | classificador procurava `not supported`, mensagem dizia `does not support` | redações extras + escada própria para a minimização |
| 20 | `cudaErrorIllegalAddress` no NVT | **ligações peptídicas falsas** de até 1,4 nm nas lacunas do cristal | `split_chain_gaps.py` + verificação de `Long Bond` |
| 21 | NVT estourava mesmo com a cadeia sã | minimização com vínculos e água rígida; NVT partindo a 310 K com Nose-Hoover e dt = 2 fs | `-DFLEXIBLE` + `em2` + equilíbrio em degraus com V-rescale |
| 22 | `nstxout-compressed = 1e8` | reusei a função que converte **ns** em passos para um valor em **ps** | duas funções, com a unidade no nome |
| 23 | réplica 2 falharia por `npt.cpt` ausente | `-t` exigido para um conteúdo que ia ser descartado | `-t` condicional + `npt.gro` recuperável de um `prod.tpr` |
| 24 | `Your tpx version is 138` | o `.tpr` do GROMACS 2026 é novo demais para o `TPRParser` do MDAnalysis | cascata de topologia: `.tpr` → `prod.gro` → `npt.gro` → conversão pelo `gmx` |
| 25 | **RMSD da proteína = 30 Å** com "frações ancoradas = 1,00" | os cortes nas lacunas fizeram da proteína 5 moléculas, envolvidas separadamente pela caixa periódica | `md_fix_pbc.sh` + recusa de emitir veredito acima de 10 Å |
| 26 | "interação espúria" onde não havia | corte absoluto de 15 **pares de átomos** para contatos warhead–E3, nunca calibrado | critérios **relativos** ao frame inicial e ao próprio recrutador |
| 27 | reprovação por RMSD global de 4,22 Å | o critério media a proteína inteira, e o receptor é multidomínio com dobradiça | portão no **sítio** (1,70 Å); o global vira contexto a declarar |

## O padrão

Três coisas atravessam a lista.

**API que devolve "não fiz nada" com o mesmo tipo de "fiz".**
`ReplaceSubstructs`, `GetBestRMS`, `MolFromMolFile` (que levanta `OSError` em
vez de devolver `None`), o `pdb2gmx` que cria o arquivo antes de terminar. Em
todos, a defesa foi a mesma: **verificar a saída, não confiar na chamada.**

**O sintoma raramente fica perto da causa.** `cudaErrorIllegalAddress` levaria
a investigar driver, versão de CUDA e memória da placa. A causa era uma
ligação peptídica de 1,4 nm criada quinze minutos antes, cujo aviso
(`Long Bond`) estava na tela e foi lido como ruído. Hoje `explain_lincs.py`
traduz os índices em nomes de resíduo, e o `md_run.sh` repete o bloco de erro
no fim do log, onde quem roda `tail` o vê.

**Número impossível é bug, não descoberta.** Uma proteína enovelada não se
move 30 Å em 200 ns, e "ancorado em 100% dos frames" não convive com isso. Duas
saídas que se contradizem são sinal de artefato, e o pipeline hoje **para**
antes de emitir veredito em vez de reprovar um candidato por causa da aresta da
caixa periódica. O mesmo vale para limiares: um corte que nunca foi comparado
com dado real não é critério, é chute com aparência de critério.

**Escrever código contra a mensagem que se tem à mão, não contra a que a
ferramenta produz.** Os obstáculos 18, 19 e 22 são meus, e todos dessa forma.
Os stubs de teste do repositório hoje reproduzem as mensagens **literais** do
GROMACS 2026.2 por isso.
""")

md(r"""
---
# Apêndice B — O que continua sendo julgamento humano

O pipeline é automático de ponta a ponta, mas quatro coisas não foram — e não
devem ser — automatizadas.

**1. Disparar os jobs do PRosettaC.** Rode **um** antes do lote e confira o
`Anchor atoms`: a convenção 0-based *vs* 1-based varia por build, e um lote
inteiro com o átomo errado é um lote perdido.

**2. Submeter os JSONs do AlphaFold 3** em `alphafoldserver.com`.

**3. Inspecionar as poses no ChimeraX** antes de comprometer dias de GPU. O
pipeline verifica geometria e números; ele não olha.

**4. Escolher os cortes.** `ANALISE_COS_MIN=0.3`, `DOCK_ROUND1_CUTOFF=-7.0`,
`ANALISE_SD_MAX=0.5`, `MD_TRUNCAR_PERTO=0` são **decisões suas**, não
constantes da natureza. A tese fica mais forte se aparecerem como decisões, com
o efeito de afrouxá-las registrado.

Some-se a isso o que é limitação **do método**, não do código, e que precisa
constar:

- O nível (ii) simula E3 + PROTAC. Ele **não** testa o complexo ternário com a
  PCSK9 — essa é a etapa do PRosettaC/AF3, e ela depende de predição.
- 200 ns é curto para reorganização de interface. Estabilidade em 200 ns é
  evidência **necessária, não suficiente**.
- GAFF2/AM1-BCC para um PROTAC de 19 torções é o padrão da área, mas é
  parametrização genérica; a barreira torsional do linker é o ponto fraco.
- 802,9 Da e 19 rotáveis: permeabilidade celular é o gargalo conhecido da
  classe, e nada neste pipeline a mede.
""")

md(r"""
---
# Apêndice C — Mapa dos arquivos

## Configuração e driver

| Arquivo | Papel |
|---|---|
| `config/pipeline.conf` | **tudo** que muda de máquina ou de projeto |
| `scripts/run_pipeline.sh` | driver das 8 fases, marcadores `.done_<n>`, `--from/--only/--list/--dry-run` |
| `scripts/status.sh` | onde o pipeline está, num comando |
| `docs/RUNBOOK.md` | como rodar, fase por fase |
| `docs/WP3_warheads_PCSK9.md` | a química das warheads |

## Fases 1–7

| Arquivo | Papel |
|---|---|
| `generate_pcsk9_warheads.py` | enumera as 180 warheads (`molzip`, átomo 0 = N de conjugação) |
| `prep_pcsk9_receptor.py` | receptor, sítio, exit vector |
| `docking_engines.py` | Vina e Uni-Dock atrás da mesma assinatura |
| `dock_warheads_pcsk9.py` | portão de validação + triagem em duas rodadas |
| `analyze_warhead_docking.py` | LE, viabilidade geométrica, estabilidade de pose |
| `rmsd_inplace.py` | RMSD que mede **sem** superpor (a armadilha nº 2) |
| `revalidate_redocking.py` | reconfere o portão do WP1 |
| `wp1_select_recruiter.py` | E3, recrutador, exit point, guarda de farmacóforo |
| `wp2_linker_tools.py` | filtros, geometria no exit vector, rótulos, montagem verificada |
| `wp2_build_subcomplexes.py` | sub-complexos recrutador-linker |
| `wp3_assemble_protacs.py` | PROTACs completos + entradas do PRosettaC |
| `rank_protacs.py` | escore composto, fila da MD |

## Fase 8

| Arquivo | Papel |
|---|---|
| `md_prepare.py` | embebe o PROTAC com o recrutador travado na pose do docking |
| `fix_receptor_for_md.py` | reconstrói cadeias laterais incompletas (ChimeraX) |
| `split_chain_gaps.py` | corta a cadeia nas lacunas do cristal |
| `build_topology.py` | ordem das diretivas que o GROMACS exige |
| `build_index.py` | grupos de acoplamento descobertos e **verificados** |
| `md_run.sh` | acpype → topologia → solvatação → equilíbrio em degraus → produção |
| `explain_lincs.py` | traduz índices do LINCS em nomes de resíduo |
| `md_fix_pbc.sh` | desfaz a condição periódica antes da análise |
| `md_analyze.py` | RMSD, contatos, veredito contra critérios explícitos |

## Autotestes

Módulos que rodam sozinhos e **demonstram** o defeito que corrigem:

```bash
python scripts/rmsd_inplace.py       # GetBestRMS dá 0,00 para pose a 25 Å
python scripts/wp2_linker_tools.py   # filtro, montagem verificada, colocação
```

Todos os scripts são determinísticos (`--seed`, default `0xC0FFEE`). Para a
tese, registre as linhas de comando exatas e as versões de RDKit e GROMACS.
""")

md(r"""
---
# Apêndice D — Parâmetros da metodologia, para a tese

Tabela única com tudo que precisa ser declarado. Os valores vêm de
`config/pipeline.conf` e dos scripts; a célula abaixo relê o `.conf` do disco
para você conferir que nada divergiu.

## Estruturas

| Item | Valor |
|---|---|
| PCSK9 | **6U26**, cadeias A + B, ligante de referência `063` |
| sítio | núcleo enterrado do `063`, centro `[38.589, 25.818, 26.443]`, caixa `[20.47, 18.46, 18.16]` Å |
| exit vector (PCSK9) | `[0.222, 0.8162, 0.5333]` |
| distância do sítio à interface EGF(A) do LDLR | **23,4 Å** |
| E3 | **CRBN**, cristal **4TZ4**, cadeia C |
| exit point (CRBN) | `[-41.719, 60.152, -86.692]`, direção `[0.3075, 0.8534, 0.4209]` |
| bibliotecas | **Chemspace** linkers + anchors (Enamine excluído por decisão) |

## Docking

| Item | Valor |
|---|---|
| motor | Uni-Dock v1.2.0 (GPU, `--gpu_batch`), função de escore do AutoDock Vina |
| portão de validação | redocking do núcleo do `063`, **RMSD 0,21 Å em 5/5 sementes** (critério ≤ 2,0 Å) |
| rodada 1 | exaustividade baixa, 180 warheads, corte `-7,0` kcal/mol |
| rodada 2 | exaustividade alta, 5 sementes |
| critérios de viabilidade | `SD_MAX = 0,5` Å · `COS_MIN = 0,3` (escolhido sem base empírica) |

## Dinâmica molecular

| Item | Valor |
|---|---|
| campos de força | AMBER **ff14SB** (proteína) · **GAFF2/AM1-BCC** (PROTAC, acpype) |
| água / caixa | **TIP3P**, cúbica, folga 1,3 nm |
| neutralização | 5 CL⁻ (carga do sistema +5) |
| eletrostática | PME, `rvdw = rcoulomb = 0,9` nm |
| vínculos | LINCS em h-bonds, `order 8`, `iter 2` |
| equilíbrio | `em` (emtol 1000, `-DFLEXIBLE`) → `em2` (emtol 100) → `warm` 20 ps @ 0,5 fs → `nvt` 100 ps @ 1 fs → `npt` 5 ns @ 2 fs, **V-rescale**, soluto restrito |
| produção | **200 ns × 3 réplicas**, dt 2 fs, **Nose-Hoover** τ 1,0 ps @ 310 K, **Parrinello-Rahman** τ 2,0 ps @ 1 bar |
| amostragem | 1000 frames por réplica (um a cada 200 ps) |
| tratamento de PBC | `-pbc cluster` + `-pbc nojump` sobre proteína + ligante, obrigatório por a proteína ser 5 moléculas |
| critérios de aprovação | **RMSD do sítio** (CA a ≤12 Å do recrutador) ≤ 3,5 Å · RMSD PROTAC ≤ 5,0 Å · ancoragem ≥ 0,5× inicial · warhead/recrutador ≤ 1,0 · crescimento warhead ≤ 2,0× inicial · frações ancoradas ≥ 0,7 |
| declarado como contexto | RMSD global e do núcleo, sempre impressos; acima de 3,5 Å exigem declaração |
| reparos no receptor | 15 cadeias laterais reconstruídas (→ 0 incompletas) · 4 quebras de cadeia nas lacunas, todas a > 16 Å do sítio |
""")

code(r"""
# Apêndice D — o pipeline.conf como está no disco, para conferência
conf = Path.home() / "Protac_env" / "config" / "pipeline.conf"
if conf.exists():
    interessa = (
        "PCSK9_PDB", "PCSK9_CHAIN", "PCSK9_REF_LIGAND", "PCSK9_BURIAL",
        "DOCK_ROUND1_CUTOFF", "DOCK_BATCH", "ANALISE_SD_MAX", "ANALISE_COS_MIN",
        "N_WARHEADS_WP2", "LINKERS_SDF", "ANCHORS_SDF", "LINKERS_EXTRA",
        "LINKER_NCONFS", "LINKER_NSPINS", "LINKER_ATOM_MIN", "LINKER_ATOM_MAX",
        "LINKER_ROTB_MAX", "N_LINKERS_WP3", "MD_N_REPLICAS", "MD_NS_PROD",
        "MD_NS_NPT", "MD_RANK", "MD_AUTO", "MD_TRUNCAR_PERTO",
    )
    for linha in conf.read_text().splitlines():
        nome = linha.split("=")[0].strip()
        if nome in interessa:
            print("  " + linha.strip())
else:
    print(f"não achei {conf}")
""")

md(r"""
---
# Apêndice E — Glossário

| Termo | O que é |
|---|---|
| **PROTAC** | molécula bifuncional que recruta uma ligase E3 para degradar uma proteína-alvo |
| **warhead** | a ponta do PROTAC que liga na proteína-alvo (aqui, na PCSK9) |
| **recrutador** | a ponta que liga na ligase E3 (aqui, na CRBN) |
| **linker** | a ponte entre as duas, bifuncional por definição |
| **complexo ternário** | alvo + PROTAC + E3 juntos; é ele que tem de ser produtivo |
| **exit vector** | direção em que o solvente é acessível a partir de um sítio; por onde o linker sai |
| **LE** (*ligand efficiency*) | −score dividido pelo número de átomos pesados; corrige o viés de tamanho |
| **portão** | critério que interrompe o pipeline, em vez de filtrar moléculas |
| **nível (ii)** | a MD de E3 + recrutador-linker-warhead, sem a proteína-alvo |
| **ff14SB / GAFF2** | campos de força para proteína e para molécula pequena |
| **AM1-BCC** | método de cargas parciais usado com o GAFF2 |
| **PME** | *Particle Mesh Ewald*, tratamento da eletrostática de longo alcance |
| **LINCS** | algoritmo que mantém os vínculos (aqui, as ligações com H) |
| **V-rescale / Nose-Hoover** | termostatos; o primeiro é robusto, o segundo é correto no ensemble |
| **Parrinello-Rahman** | barostato usado na produção |
| **offload** | passar tarefas (não-ligadas, PME, ligadas, integração) para a GPU |
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
n_linhas = sum(len(c["source"]) for c in CELLS)
print(f"{out}")
print(f"  {len(CELLS)} células: {n_md} markdown, {len(CELLS) - n_md} código")
print(f"  {n_linhas} linhas de conteúdo")
