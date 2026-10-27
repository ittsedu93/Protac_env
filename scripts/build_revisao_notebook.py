#!/usr/bin/env python
"""
Monta notebooks/protac_revisao_parametros.ipynb — a coleta dos parâmetros da
aba PAPERS nas quatro bases de PROTAC.

O notebook é a interface: você preenche UMA célula (o alvo) e roda o resto.
A lógica toda está em scripts/protac_fontes.py, que é testado em
tests/test_protac_fontes.py — o notebook não repete nada dela, só chama.

Gerado por script para ser reproduzível: para mudar o notebook, edite aqui e
rode de novo.

    python scripts/build_revisao_notebook.py
"""

import json
from pathlib import Path

CELLS = []


def _id() -> str:
    """nbformat 4.5+ exige id por célula; sem ele o Jupyter avisa a cada leitura."""
    return f"c{len(CELLS):03d}"


def md(src: str):
    CELLS.append({"cell_type": "markdown", "id": _id(), "metadata": {},
                  "source": src.strip("\n").splitlines(keepends=True)})


def code(src: str):
    CELLS.append({"cell_type": "code", "id": _id(), "execution_count": None,
                  "metadata": {}, "outputs": [],
                  "source": src.strip("\n").splitlines(keepends=True)})


# ===========================================================================
# CAPA
# ===========================================================================
md(r"""
# Revisão PROTAC — coleta automatizada dos parâmetros da aba PAPERS

Preenche as 28 colunas da aba **PAPERS** da sua planilha a partir de quatro
bases: **PROTAC-DB**, **PROTACpedia**, **PROTAC-Bench** e **PROTAC-PatentDB**.

**Você preenche uma célula só** — a §1, com o alvo. O resto roda em sequência.

## Como usar

1. Rode a §0 (preparação).
2. Preencha a **§1** com o seu alvo e rode.
3. Rode a §2 → §6 na ordem. A §6 grava o `.xlsx`.

## O que esperar de cada base

Nenhuma base preenche as 28 colunas, e isso é esperado — você já disse que não
há problema. Vale saber de onde vem o quê, porque muda o que você deve
procurar na mão depois:

| Base | O que preenche bem | O que não tem |
|---|---|---|
| **PROTAC-DB** | DC₅₀, Dmax, IC₅₀, MW, DOI, linhagem, E3, tipo de linker; na 3.0 também PK (Tmax, T½, Cmax, AUC, CL, Vss) e permeabilidade | — é a mais completa para esta planilha |
| **PROTACpedia** | tipo de linker, anotação de ensaio, DC₅₀/Dmax, DOI | PK; **o download exige cadastro gratuito** |
| **PROTAC-Bench** | SMILES, E3, rótulo de atividade, DC₅₀ quando existe | bibliografia, PK, permeabilidade; **identifica o alvo só por UniProt** |
| **PROTAC-PatentDB** | MW e permeabilidade (**previstas**, não medidas), patente, data | potência medida, DOI de artigo |

Duas consequências práticas, que o notebook trata mas vale você saber:

- **O PROTAC-Bench só conhece o alvo por acessão UniProt.** Filtrar por `BRD4`
  zera essa base. Por isso a §1 pede o UniProt junto — e a §3 avisa quando
  isso aconteceu, em vez de deixar passar como "não há dado".
- **O PROTAC-PatentDB traz valor previsto por ADMETlab, não medido.** Entra na
  planilha, mas a aba `PROVENIENCIA` registra a origem de cada linha. Não
  misture com medida experimental sem marcar — numa revisão isso é erro grave.

## Honestidade sobre as colunas

O notebook **não depende do nome exato das colunas de cada base**. Cada
parâmetro tem uma lista de apelidos (`APELIDOS`, em `scripts/protac_fontes.py`)
comparada contra o nome normalizado da coluna de origem — `DC50 (nM)`, `DC₅₀`
e `dc50_nm` são reconhecidos como a mesma coisa.

Isso é deliberado. Se os nomes estivessem fixos no código, a próxima versão de
qualquer uma das bases quebraria a coleta **em silêncio**: a coluna viria vazia
e pareceria "a base não tem esse dado", que é o erro mais difícil de notar numa
revisão. Em vez disso, a §3 imprime as **colunas de origem que ficaram sem
destino** — é ali que aparece um parâmetro que a base tem e o mapa não
reconheceu. Para corrigir, você acrescenta um apelido e roda de novo.

Pelo mesmo motivo, **nenhum caminho de arquivo está chutado no código**: o
módulo lê a página de download de cada base e pega os links que encontrar (ou
usa a API pública do Figshare e do Hugging Face).
""")


# ===========================================================================
# 0. PREPARAÇÃO
# ===========================================================================
md(r"""
---
## §0 — Preparação

Precisa de `pandas` e `openpyxl`, que já estão nos envs do projeto. `rdkit` é
**opcional**: com ele a deduplicação usa SMILES canônico e o comprimento do
linker pode ser estimado; sem ele a comparação é por texto e o notebook avisa.
""")

code(r"""
import sys
from pathlib import Path

# A lógica vive em scripts/protac_fontes.py. Se você rodar este notebook de
# fora do repositório, ajuste RAIZ para a pasta do projeto.
RAIZ = Path.cwd()
if not (RAIZ / "scripts" / "protac_fontes.py").exists():
    RAIZ = RAIZ.parent          # típico: o notebook está em notebooks/
if not (RAIZ / "scripts" / "protac_fontes.py").exists():
    raise SystemExit(
        "não achei scripts/protac_fontes.py — ajuste RAIZ para a pasta do projeto"
    )

sys.path.insert(0, str(RAIZ / "scripts"))
import protac_fontes as pf

print("módulo:", pf.__file__)
print("parâmetros da aba PAPERS:", len(pf.CHAVES))
print("bases:", ", ".join(pf.FONTES))
print("RDKit:", "sim" if pf.rdkit_disponivel() else
      "não (deduplicação por texto; sem estimativa de comprimento de linker)")
""")


# ===========================================================================
# 1. O INPUT
# ===========================================================================
md(r"""
---
## §1 — O input (a única célula que você precisa editar)

`ALVO` costuma bastar. Os outros campos existem por motivos concretos:

- **`UNIPROT`** — o PROTAC-Bench não traz nome de gene. Sem a acessão, aquela
  base sai zerada. Pegue em [uniprot.org](https://www.uniprot.org/).
- **`SINONIMOS`** — as bases escrevem a mesma proteína de formas diferentes
  (`PCSK9` vs. `proprotein convertase subtilisin/kexin type 9`). Sem o
  sinônimo, linhas do seu alvo são descartadas.
- **`E3`**, **`LINHAGEM`** — só se você quiser restringir. Deixe `None` para não.
- **`SO_COM_DEGRADACAO`** — `True` descarta o que não tem DC₅₀ nem Dmax. Útil
  para a revisão de atividade; atrapalha se você quer cobertura química.
- **`UNIDADES_FIXAS`** — `True` converte potências para nM como número puro
  (bom para ordenar e fazer gráfico). `False` mantém o texto com unidade e
  qualificador (`>10 μM`), que é o que faz sentido ler na planilha.

> Os critérios combinam em **E** (alvo **e** E3 **e** linhagem). O alvo, os
> sinônimos e o UniProt combinam entre si em **OU** — são nomes do mesmo alvo.
""")

code(r"""
# =====================  PREENCHA AQUI  =====================================

ALVO = "BRD4"                  # proteína-alvo
UNIPROT = "O60885"             # acessão; necessária para o PROTAC-Bench
SINONIMOS = ["BRD-4", "bromodomain containing 4", "bromodomain 4"]

E3 = None                      # "VHL", "CRBN", ... ou None
LINHAGEM = None                # "MV4-11", ... ou None
TEXTO_LIVRE = None             # busca em qualquer coluna, ou None

BASES = ["PROTAC-DB", "PROTACpedia", "PROTAC-Bench", "PROTAC-PatentDB"]

SO_COM_DEGRADACAO = False      # True = exige DC₅₀ ou Dmax
UNIDADES_FIXAS = False         # True = potências em nM, número puro
CONSOLIDAR = True              # True = junta o mesmo composto vindo de 2 bases

DIR_DADOS = RAIZ / "dados" / "bases"   # onde os arquivos das bases ficam
BAIXAR = True                          # False = usa só o que já está em disco

# Sua planilha, para a saída sair com a formatação dela. None = cria do zero.
# O arquivo apontado aqui NUNCA é sobrescrito.
MODELO = None                  # p.ex. RAIZ / "dados" / "Revisao_PROTAC.xlsx"

SAIDA = RAIZ / "resultados" / f"revisao_{(ALVO or 'tudo').replace('/', '_')}.xlsx"
SAIDA_CSV = SAIDA.with_suffix(".csv")

# ===========================================================================

# A ORDEM DE `BASES` É ORDEM DE PREFERÊNCIA: quando o mesmo composto aparece em
# duas bases, o valor da primeira ganha. Com PROTAC-DB antes de
# PROTAC-PatentDB, a medida experimental ganha da prevista — que é o que se
# quer. Se você inverter, inverte isso também.
pedido = pf.Pedido(
    alvo=ALVO,
    uniprot=UNIPROT,
    sinonimos=SINONIMOS,
    e3=E3,
    linhagem=LINHAGEM,
    texto_livre=TEXTO_LIVRE,
    bases=BASES,
    so_com_degradacao=SO_COM_DEGRADACAO,
    unidades_fixas=UNIDADES_FIXAS,
    consolidar_duplicatas=CONSOLIDAR,
)

print("pedido montado:")
print(f"  alvo .......... {ALVO}  (uniprot {UNIPROT})")
print(f"  sinônimos ..... {', '.join(SINONIMOS) if SINONIMOS else '—'}")
print(f"  E3 / linhagem . {E3 or '—'} / {LINHAGEM or '—'}")
print(f"  bases ......... {', '.join(BASES)}")
print(f"  dados em ...... {DIR_DADOS}")
print(f"  saída ......... {SAIDA}")
""")


# ===========================================================================
# 2. OS ARQUIVOS DAS BASES
# ===========================================================================
md(r"""
---
## §2 — Conseguir os arquivos

A célula abaixo, para cada base: procura em `DIR_DADOS` o que já existe; se não
houver nada e `BAIXAR=True`, descobre os links e baixa.

**O que pode dar errado, e o que fazer:**

- **PROTACpedia** exige cadastro (gratuito) para o download completo. O módulo
  não automatiza o login — seria guardar sua senha num notebook. Faça o login
  em [protacpedia.weizmann.ac.il](https://protacpedia.weizmann.ac.il/ptcb/main),
  baixe o CSV/SDF e salve em `DIR_DADOS/PROTACpedia/`.
- **Rede fechada ou base fora do ar**: a mensagem diz em qual página baixar na
  mão e em qual pasta salvar. Uma base que falha **não derruba as outras três**.
- **Formato**: `.csv`, `.tsv`, `.xlsx`, `.json`, `.parquet`, `.sdf`, `.zip` e
  `.gz` são lidos. Para SDF só o bloco de propriedades é lido (sem RDKit), então
  se o SDF não tiver o SMILES como propriedade, prefira o CSV da mesma base.

Para só inspecionar sem baixar nada, ponha `BAIXAR = False` na §1.
""")

code(r"""
DIR_DADOS.mkdir(parents=True, exist_ok=True)

for nome_base in BASES:
    fonte = pf.FONTES[nome_base]
    arquivos, registro = pf.obter_arquivos(fonte, DIR_DADOS, baixar=BAIXAR)
    print(f"\n{'=' * 70}\n{nome_base}\n{'=' * 70}")
    print(f"  {fonte.descricao}")
    print(f"  página: {fonte.pagina}")
    if fonte.observacao:
        print(f"  nota:   {fonte.observacao}")
    for linha in registro:
        print(f"  · {linha}")
    if arquivos:
        for caminho in arquivos:
            tamanho = caminho.stat().st_size / 1024
            print(f"  → {caminho.name}  ({tamanho:,.0f} kB)")
    else:
        print(f"  → nada em disco. Salve o arquivo em: {DIR_DADOS / nome_base}/")
""")


# ===========================================================================
# 3. A COLETA
# ===========================================================================
md(r"""
---
## §3 — Coletar, mapear e filtrar

Lê cada arquivo, traduz as colunas para o vocabulário da planilha, aplica o
filtro da §1, padroniza unidades e junta duplicatas entre bases.

**Leia o relatório inteiro**, não só a contagem final. Três partes importam:

1. **Preenchimento por parâmetro** — quantas linhas têm cada dado. É o mapa do
   que vai sobrar para você buscar na mão.
2. **Colunas de origem sem destino** — colunas que as bases trazem e o mapa não
   reconheceu. **É aqui que aparece dado perdido.** Se algo útil estiver nessa
   lista, acrescente um apelido em `pf.APELIDOS` (§5) e rode de novo.
3. **Avisos** — em especial o de base que não identifica o alvo por nome.
""")

code(r"""
resultado = pf.coletar(pedido, dir_dados=DIR_DADOS, baixar=BAIXAR)
print(resultado.relatorio())
""")

md(r"""
### A tabela

As colunas aparecem com os nomes internos (ASCII); os rótulos da planilha
entram na exportação. Vazio (`<NA>`) significa "a base não traz", e é esperado.
""")

code(r"""
import pandas as pd

pd.set_option("display.max_columns", 50)
pd.set_option("display.width", 220)
pd.set_option("display.max_colwidth", 28)

tabela = resultado.tabela
print(f"{len(tabela)} linha(s)\n")

if tabela.empty:
    print("Nada encontrado. Antes de concluir que não há dado, verifique:")
    print("  · a §2 conseguiu algum arquivo? (sem arquivo não há o que filtrar)")
    print("  · o alvo está escrito como as bases escrevem? acrescente SINONIMOS")
    print("  · o UNIPROT está preenchido? sem ele o PROTAC-Bench sai zerado")
    print("  · SO_COM_DEGRADACAO=True pode estar descartando tudo")
else:
    colunas_vista = [c for c in ["idx", "fonte", "alvo", "e3", "dc50", "dmax",
                                 "ic50", "mw", "linhagem", "doi", "smiles"]
                     if c in tabela.columns]
    display(tabela[colunas_vista].head(30))
""")


# ===========================================================================
# 4. CONFERÊNCIA
# ===========================================================================
md(r"""
---
## §4 — Conferência antes de exportar

Três perguntas que vale responder antes de levar isso para a revisão.
""")

code(r"""
if not tabela.empty:
    print("1. Quantas linhas cada base contribuiu (após consolidação)")
    print("   Uma linha pode citar duas bases: é o mesmo composto confirmado")
    print("   nas duas, não contagem dobrada.")
    for fonte_nome, n in tabela["fonte"].value_counts().items():
        print(f"     {fonte_nome}: {n}")

    print("\n2. Parâmetros que nenhuma linha preencheu")
    print("   Candidatos a busca manual no artigo original.")
    vazios = [pf.ROTULOS[c] for c in pf.CHAVES
              if c != "idx" and c in tabela.columns and tabela[c].isna().all()]
    print("     " + (", ".join(vazios) if vazios else "nenhum — todas têm algum dado"))

    print("\n3. Linhas sem DOI")
    print("   Sem referência não dá para citar na revisão; as de patente têm")
    print("   número de patente na aba PROVENIENCIA.")
    if "doi" in tabela.columns:
        sem_doi = int(tabela["doi"].isna().sum())
        print(f"     {sem_doi} de {len(tabela)}")
""")


# ===========================================================================
# 5. AJUSTAR O MAPA
# ===========================================================================
md(r"""
---
## §5 — Quando uma coluna ficou sem destino

Se a §3 listou uma coluna de origem útil em **"COLUNAS DE ORIGEM SEM DESTINO"**,
é porque o mapa não a reconheceu. Acrescente o apelido e rode a §3 de novo.

O padrão é comparado contra o nome **normalizado**: minúsculas, sem acento, sem
pontuação, e com `₅₀` virando `50`. Confira com `pf.normalizar_nome(...)` antes
de escrever o padrão — é mais rápido que adivinhar.

A célula abaixo está comentada de propósito: descomente e adapte.
""")

code(r"""
# Como o nome chega ao mapa:
for exemplo in ["DC₅₀ (nM)", "Degradation DC50", "dc50_nm", "T1/2 (h)"]:
    print(f"  {exemplo!r:26} -> {pf.normalizar_nome(exemplo)!r}")

# --- descomente e adapte -------------------------------------------------
# pf.APELIDOS["dc50"].append(r"\bdegradation potency\b")
# pf.APELIDOS["t_meia"].append(r"\bterminal half life\b")
#
# Para um parâmetro que a planilha não tem, acrescente antes a coluna em
# pf.PARAMETROS — senão o valor é lido e descartado na exportação.

# Confira o efeito numa lista de colunas de verdade, sem reprocessar tudo:
# mapa, auditoria = pf.mapear_colunas(["Smiles", "Degradation Potency (nM)"])
# print(mapa); print(auditoria["sem_destino"])
""")


# ===========================================================================
# 6. EXPORTAR
# ===========================================================================
md(r"""
---
## §6 — Exportar

Grava a aba **PAPERS** com as 28 colunas na ordem exata da planilha: cabeçalho
na linha 3, dados a partir da linha 4.

Com `MODELO` apontando para a sua planilha, a formatação dela é preservada.
**O modelo nunca é sobrescrito** — a saída é sempre `SAIDA`, e o módulo recusa
gravar por cima do modelo.

Vão duas abas de auditoria, e elas são o que torna a revisão defensável:

- **`PROVENIENCIA`** — de qual base e de qual arquivo veio cada linha, com o id
  do composto na base de origem. É aqui que fica registrado que o valor do
  PROTAC-PatentDB é **previsto**, não medido.
- **`MAPEAMENTO`** — qual coluna de qual base alimentou cada parâmetro, mais o
  que ficou sem destino e o que a base não traz.
""")

code(r"""
if tabela.empty:
    print("nada para exportar — reveja a §3")
else:
    caminho = pf.exportar_planilha(tabela, SAIDA, modelo=MODELO,
                                   resultado=resultado)
    print(f"planilha: {caminho}")
    print(f"csv:      {pf.exportar_csv(tabela, SAIDA_CSV)}")
    print(f"\n{len(tabela)} linha(s), {len(pf.CHAVES)} colunas na aba PAPERS")
    print("abas: PAPERS (dados) + PROVENIENCIA (origem) + MAPEAMENTO (auditoria)")
""")

md(r"""
### Conferência da saída

Lê de volta o arquivo gravado e compara o cabeçalho com o que a planilha espera.
Rápido, e pega o caso em que o modelo tem a aba com outro nome ou outro layout.
""")

code(r"""
if not tabela.empty:
    from openpyxl import load_workbook

    ws = load_workbook(caminho)["PAPERS"]
    cabecalho = [ws.cell(row=3, column=c).value for c in range(1, len(pf.CHAVES) + 1)]
    esperado = [pf.ROTULOS[k] for k in pf.CHAVES]

    if cabecalho == esperado:
        print(f"cabeçalho confere ({len(esperado)} colunas)")
    else:
        print("CABEÇALHO DIVERGENTE — compare:")
        for i, (saiu, devia) in enumerate(zip(cabecalho, esperado), start=1):
            if saiu != devia:
                print(f"  col {i}: saiu {saiu!r}, devia ser {devia!r}")

    preenchidas = sum(1 for linha in ws.iter_rows(min_row=4, max_row=3 + len(tabela))
                      if any(c.value is not None for c in linha))
    print(f"linhas com dado a partir da linha 4: {preenchidas}")
""")


# ===========================================================================
# 7. REPETIR PARA VÁRIOS ALVOS
# ===========================================================================
md(r"""
---
## §7 — Vários alvos de uma vez

Para uma revisão com vários alvos, o laço abaixo evita reabrir o notebook a
cada um. Um arquivo por alvo. Os arquivos das bases já estão em disco depois da
§2, então aqui vai `baixar=False` — não baixa de novo a cada alvo.
""")

code(r"""
# Preencha com os seus alvos. O UniProt de cada um importa pelo PROTAC-Bench.
ALVOS = [
    # ("BRD4",  "O60885", ["bromodomain 4"]),
    # ("PCSK9", "Q8NBP7", ["proprotein convertase subtilisin/kexin type 9"]),
]

for nome_alvo, uniprot_alvo, sinonimos_alvo in ALVOS:
    ped = pf.Pedido(alvo=nome_alvo, uniprot=uniprot_alvo,
                    sinonimos=sinonimos_alvo, bases=BASES,
                    unidades_fixas=UNIDADES_FIXAS,
                    consolidar_duplicatas=CONSOLIDAR)
    res = pf.coletar(ped, dir_dados=DIR_DADOS, baixar=False)
    if res.tabela.empty:
        print(f"{nome_alvo}: nada encontrado")
        continue
    destino = RAIZ / "resultados" / f"revisao_{nome_alvo}.xlsx"
    pf.exportar_planilha(res.tabela, destino, modelo=MODELO, resultado=res)
    com_dc50 = int(res.tabela["dc50"].notna().sum()) if "dc50" in res.tabela else 0
    print(f"{nome_alvo}: {len(res.tabela)} linha(s) "
          f"({com_dc50} com DC₅₀) -> {destino.name}")
""")

md(r"""
---
## Fora do notebook

O mesmo trabalho pela linha de comando, para entrar em script ou SLURM:

```bash
python scripts/protac_fontes.py \
    --alvo BRD4 \
    --sinonimo "bromodomain 4" \
    --dir-dados dados/bases \
    --modelo dados/Revisao_PROTAC.xlsx \
    --saida resultados/revisao_BRD4.xlsx \
    --saida-csv resultados/revisao_BRD4.csv
```

`python scripts/protac_fontes.py --help` lista o resto. Os testes da lógica:

```bash
python -m unittest discover -s tests -v
```

## Para citar

As quatro bases devem ser citadas na revisão; `pf.FONTES[...].pagina` tem o
endereço de cada uma. O que a aba `PROVENIENCIA` registra por linha é o que
permite dizer, na defesa, de onde veio cada número — e separar o medido
(PROTAC-DB, PROTACpedia) do previsto (PROTAC-PatentDB).
""")


# ===========================================================================
NOTEBOOK = {
    "cells": CELLS,
    "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python",
                       "name": "python3"},
        "language_info": {"name": "python", "version": "3.11"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}

destino = Path(__file__).resolve().parent.parent / "notebooks" / "protac_revisao_parametros.ipynb"
destino.parent.mkdir(parents=True, exist_ok=True)
destino.write_text(json.dumps(NOTEBOOK, ensure_ascii=False, indent=1) + "\n",
                   encoding="utf-8")
print(f"{destino}  ({len(CELLS)} células)")
