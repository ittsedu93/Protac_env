#!/usr/bin/env python
"""
Monta notebooks/protac_revisao_parametros.ipynb.

O notebook extrai TPD-DB + PROTAC-DB + PROTACpedia para a aba PAPERS da
RevisãoPROTAC.xlsx, remove duplicatas e pinta cada linha pela base de origem.

É AUTOCONTIDO de propósito: roda na workstation do lab sem importar nada deste
repositório, só pandas + openpyxl (RDKit é opcional).

    python scripts/build_revisao_notebook.py
"""

import json
from pathlib import Path

CELLS = []


def _id() -> str:
    return f"c{len(CELLS):03d}"


def md(src: str):
    CELLS.append({"cell_type": "markdown", "id": _id(), "metadata": {},
                  "source": src.strip("\n").splitlines(keepends=True)})


def code(src: str):
    CELLS.append({"cell_type": "code", "id": _id(), "execution_count": None,
                  "metadata": {}, "outputs": [],
                  "source": src.strip("\n").splitlines(keepends=True)})


# ===========================================================================
md(r'''
# Revisão PROTAC — TPD-DB + PROTAC-DB + PROTACpedia → RevisãoPROTAC.xlsx

Preenche as 28 colunas da aba **PAPERS**, remove PROTACs duplicados e pinta
cada linha com a cor da base de origem.

**Edite só a §1** (os caminhos). Depois `Run All`.

## Entrada

| Base | Arquivo(s) | Linhas |
|---|---|---|
| TPD-DB | `PROTAC_main_table.txt` + `PROTAC_activity.txt` | ~21 400 + ~23 300 |
| PROTAC-DB | `protacdb_20220210.csv` | ~1 200 |
| PROTACpedia | `protac.xlsx` | ~15 500 |

## Cores da aba PAPERS

| Cor | Origem |
|---|---|
| 🟦 azul | só TPD-DB |
| 🟩 verde | só PROTAC-DB |
| 🟨 amarelo | só PROTACpedia |
| 🟧 laranja | **mesma molécula em mais de uma base** (linha consolidada) |

## Duas coisas que você precisa saber

**1. Os rótulos das bases podem estar trocados.** Pelo conteúdo,
`protacdb_20220210.csv` tem campos de curadoria colaborativa (`Curator`,
`Status`, `Best PROTAC`) que são a assinatura da **PROTACpedia**, enquanto
`protac.xlsx` tem as 89 colunas com PAMPA/Caco-2 e anotação de ensaio que são a
assinatura do **PROTAC-DB 3.0**. Mantive os nomes que você usou, mas eles ficam
em `ROTULOS` na §1 — se concordar que estão invertidos, troque as duas strings
ali e o resto (cores, legenda, proveniência) acompanha.

**2. Oito colunas vão sair vazias, e não é defeito.** Nenhuma das três bases
traz farmacocinética: `Tmax`, `T 1/2`, `Cmax`, `AUC`, `Vz`, `Vss`, `CL`, `MRT`.
Também não há `Título Paper`, `Data` nem `Classe do Alvo` em nenhuma delas. A
§5 lista exatamente o que ficou vazio.
''')

# ===========================================================================
md(r'''
---
## §1 — Entrada

Os arquivos são localizados **pelo nome**, dentro de `PASTA`. Se eles estão na
raiz do repositório (ao lado de `notebooks/`), não há nada para editar: rode e
siga. A célula imprime o que encontrou antes de usar.
''')

code(r'''
from pathlib import Path

# ====================  PREENCHA AQUI  ======================================

# Pasta com os arquivos das bases. "." = a raiz do repositório.
# Rodando o notebook de dentro de notebooks/, a raiz é a pasta de cima.
PASTA = Path("..") if Path("../scripts").is_dir() else Path(".")

SAIDA = PASTA / "resultados" / "RevisaoPROTAC_preenchida.xlsx"

# Filtro opcional. None = todas as moléculas das três bases (~35 mil linhas).
# Para restringir, p.ex.: ALVO = "BRD4"   /   ALVO = "AR"
ALVO = None

# Se os rótulos estiverem trocados (ver a nota na capa), inverta as duas
# strings abaixo — cores, legenda e proveniência seguem automaticamente.
ROTULOS = {
    "csv":   "PROTAC-DB",     # o que está em protacdb_20220210.csv
    "xlsx":  "PROTACpedia",   # o que está em protac.xlsx
    "tpd":   "TPD-DB",
}

CORES = {
    ROTULOS["tpd"]:  "DDEBF7",   # azul
    ROTULOS["csv"]:  "E2EFDA",   # verde
    ROTULOS["xlsx"]: "FFF2CC",   # amarelo
    "múltiplas":     "FCE4D6",   # laranja
}

# ===========================================================================

SAIDA.parent.mkdir(parents=True, exist_ok=True)


def achar(descricao, padroes, excluir=(), obrigatorio=True):
    """Localiza um arquivo por padrão de nome, sem depender da grafia exata.

    Os nomes variam (espaço em 'Revisão PROTAC.xlsx', data no
    'protacdb_20220210.csv'), e um caminho fixo falharia por uma letra. Procura
    em PASTA e nas subpastas, ignora maiúsculas/minúsculas e, quando acha mais
    de um, fica com o mais recente e avisa.
    """
    vistos = {}
    for padrao in padroes:
        for caminho in PASTA.rglob("*"):
            if not caminho.is_file():
                continue
            nome = caminho.name.lower()
            if re.fullmatch(padrao.lower(), nome) and not any(
                    e.lower() in nome for e in excluir):
                vistos[caminho.resolve()] = caminho
    achados = sorted(vistos.values(), key=lambda c: c.stat().st_mtime, reverse=True)
    if not achados:
        if obrigatorio:
            raise FileNotFoundError(
                f"não encontrei {descricao} em {PASTA.resolve()}\n"
                f"  procurei por: {', '.join(padroes)}\n"
                f"  ajuste PASTA, ou atribua o caminho na mão a esta variável."
            )
        return None
    if len(achados) > 1:
        print(f"  ! vários candidatos a {descricao}; usando o mais recente: "
              f"{achados[0].name}  (outros: {', '.join(c.name for c in achados[1:])})")
    return achados[0]


import re  # usado por achar()

TPD_MAIN     = achar("TPD-DB tabela principal", [r".*main.?table.*\.txt"])
TPD_ACTIVITY = achar("TPD-DB atividade",        [r".*activity.*\.txt"])
PROTACDB_CSV = achar("PROTAC-DB (csv)",         [r".*protacdb.*\.csv", r".*protac.?db.*\.csv"])
PROTACPEDIA  = achar("PROTACpedia (xlsx)",      [r"protac\.xlsx", r".*protacpedia.*\.xlsx"],
                     excluir=("revis",))
# A planilha-molde. Se não aparecer, o notebook cria o arquivo do zero com os
# mesmos cabeçalhos. ESTE ARQUIVO NUNCA É SOBRESCRITO.
MODELO = achar("planilha RevisãoPROTAC", [r".*revis.*protac?.*\.xlsx"],
               obrigatorio=False)

print("entrada encontrada:")
for rotulo, caminho in [("TPD-DB main", TPD_MAIN), ("TPD-DB activity", TPD_ACTIVITY),
                        ("PROTAC-DB", PROTACDB_CSV), ("PROTACpedia", PROTACPEDIA)]:
    print(f"  {rotulo:<16} {caminho}")
print(f"  {'molde':<16} {MODELO if MODELO else '(nenhum — criarei a planilha do zero)'}")
print(f"\nsaída:  {SAIDA}")
print(f"filtro: {ALVO or 'nenhum (todas as moléculas)'}")
''')

# ===========================================================================
md(r'''
---
## §2 — As 28 colunas e as funções de leitura

Os rótulos são copiados da planilha **com os defeitos de digitação** (`Gl₅₀`
com L, `Comprimeto`, os espaços sobrando). Se eu os corrigisse, o arquivo
deixaria de bater coluna a coluna com a sua planilha.
''')

code(r'''
import re
import warnings
import pandas as pd

warnings.simplefilter("ignore")

# (chave interna, rótulo exato da planilha)
PARAMETROS = [
    ("idx", "-"), ("titulo", "Título Paper "), ("doi", "DOI"),
    ("smiles", "SMILES"), ("mw", "MW"), ("data", "Data"),
    ("ic50", "IC₅₀"), ("ec50", "EC₅₀"), ("gi50", "Gl₅₀"), ("ed50", "ED₅₀"),
    ("gr50", "GR₅₀"), ("pampa", "PAMPA "), ("caco2", "Caco-2 "),
    ("tmax", "Tmax"), ("t_meia", "T 1/2"), ("cmax", "Cmax"), ("auc", "AUC"),
    ("vz", "Vz"), ("vss", "Vss"), ("cl", "CL"), ("mrt", "MRT"),
    ("dc50", "DC₅₀"), ("dmax", "Dmax"), ("e3", "E3 Ligase Usada"),
    ("linker_tipo", "Tipo de linker"), ("linker_comprimento", "Comprimeto do Linker"),
    ("classe_alvo", "Classe do Alvo "), ("linhagem", "Linhagem Celular"),
]
CHAVES = [k for k, _ in PARAMETROS]
ROTULOS_COL = dict(PARAMETROS)

# colunas extras, que saem na aba PROVENIENCIA (não na PAPERS)
EXTRAS = ["fonte", "id_fonte", "nome", "alvo", "uniprot", "referencia"]


def limpar(v):
    """Normaliza célula e transforma os vários 'vazios' das bases em NA."""
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return pd.NA
    s = str(v).strip()
    if s in ("", ".", "-", "NA", "N/A", "nan", "None", "null", "ND", "n.d."):
        return pd.NA
    return re.sub(r"\s+", " ", s)


_RE_NUM = re.compile(r"^\s*(?P<q>[<>≤≥~≈]{0,2})\s*(?P<v>[-+]?\d+(?:[.,]\d+)?(?:[eE][-+]?\d+)?)")


def valor_num(v):
    """Número de uma célula tipo '417nM', '>10000.0nM', '95%'. None se não houver.

    Serve só para ESCOLHER entre medições; o que vai para a planilha é sempre
    o texto original, com unidade e qualificador.
    """
    s = limpar(v)
    if s is pd.NA:
        return None
    m = _RE_NUM.match(str(s))
    if not m:
        return None
    try:
        return float(m.group("v").replace(",", "."))
    except ValueError:
        return None


# Massa atômica média, para calcular MW a partir da fórmula do TPD-DB (que não
# traz MW). Cobre os elementos que aparecem em PROTAC.
_MASSAS = {
    "C": 12.011, "H": 1.008, "N": 14.007, "O": 15.999, "S": 32.06,
    "P": 30.974, "F": 18.998, "Cl": 35.45, "Br": 79.904, "I": 126.904,
    "Si": 28.085, "B": 10.81, "Se": 78.971, "Na": 22.990, "K": 39.098,
    "Li": 6.94, "D": 2.014,
}
_RE_FORMULA = re.compile(r"([A-Z][a-z]?)(\d*)")


def mw_da_formula(formula):
    """MW médio a partir de 'C34H31F3N8O4'. None se a fórmula tiver algo fora da tabela."""
    s = limpar(formula)
    if s is pd.NA:
        return None
    s = str(s).replace(" ", "")
    if not s or not re.fullmatch(r"(?:[A-Z][a-z]?\d*)+", s):
        return None
    total = 0.0
    for elem, qtd in _RE_FORMULA.findall(s):
        if elem not in _MASSAS:
            return None
        total += _MASSAS[elem] * (int(qtd) if qtd else 1)
    return round(total, 2) if total else None


def quadro_vazio(n):
    """DataFrame com todas as colunas canônicas, para cada base preencher."""
    return pd.DataFrame({c: pd.Series([pd.NA] * n, dtype="object")
                         for c in CHAVES + EXTRAS})


print(f"{len(CHAVES)} colunas da aba PAPERS definidas")
''')

# ===========================================================================
md(r'''
---
## §3 — Ler as três bases

### TPD-DB

A atividade vem em **formato longo** (uma linha por medição), então precisa ser
pivotada: ~23 300 medições viram uma linha por molécula. Quando há várias
medições do mesmo tipo para a mesma molécula, a escolha é:

- **DC₅₀, IC₅₀, EC₅₀, GI₅₀, ED₅₀, GR₅₀** → a **menor** (mais potente)
- **Dmax** → o **maior** (maior degradação)

O que vai para a célula é sempre o **texto original** da medição escolhida
(`417nM`, `>10000.0nM`), não um número reconvertido — e a linhagem celular é a
daquela medição. Medições qualitativas (`+++`) entram como texto quando não há
nenhuma numérica.

`Activity Type` tem sufixos (`DC50_1`, `DC50 (SU-DHL-2)`, `IC50 (EGFR WT)`), então
o tipo é reconhecido pelo **início** da string.
''')

code(r'''
# --- TPD-DB: tabela principal -------------------------------------------
tpd = pd.read_csv(TPD_MAIN, sep="\t", dtype=str, keep_default_na=False)
print(f"TPD-DB main:     {len(tpd):>6} linhas, colunas: {list(tpd.columns)}")

# --- TPD-DB: atividade (formato longo) ----------------------------------
act = pd.read_csv(TPD_ACTIVITY, sep="\t", dtype=str, keep_default_na=False)
print(f"TPD-DB activity: {len(act):>6} linhas, colunas: {list(act.columns)}")

# Qual chave da planilha cada Activity Type alimenta. A ordem importa: 'DIC50'
# começa com 'D', então 'DC50' tem de ser testado antes de 'IC50' não casar
# errado — por isso comparo com startswith sobre a lista nesta ordem.
TIPOS = [("DC50", "dc50"), ("DMAX", "dmax"), ("IC50", "ic50"),
         ("EC50", "ec50"), ("GI50", "gi50"), ("ED50", "ed50"), ("GR50", "gr50")]
MAIOR_MELHOR = {"dmax"}          # Dmax: quanto maior, melhor


def tipo_da_medicao(texto):
    t = str(texto).strip().upper()
    for prefixo, chave in TIPOS:
        if t.startswith(prefixo):
            return chave
    return None


act["_chave"] = act["Activity Type"].map(tipo_da_medicao)
act["_num"] = act["Activity"].map(valor_num)
reconhecidas = act["_chave"].notna().sum()
print(f"  medições reconhecidas: {reconhecidas} de {len(act)} "
      f"({len(act) - reconhecidas} são de tipos que a planilha não tem, "
      f"p.ex. 'Residual rate', 'Cell viability')")

# pivot: uma linha por TPD ID, uma coluna por parâmetro
melhor = {}
for (tpd_id, chave), grupo in act[act["_chave"].notna()].groupby(
        ["TPD ID", "_chave"], sort=False):
    com_num = grupo[grupo["_num"].notna()]
    if len(com_num):
        linha = (com_num.loc[com_num["_num"].idxmax()] if chave in MAIOR_MELHOR
                 else com_num.loc[com_num["_num"].idxmin()])
    else:
        linha = grupo.iloc[0]
    melhor.setdefault(tpd_id, {})[chave] = (limpar(linha["Activity"]),
                                            limpar(linha["Cell Line"]))

print(f"  moléculas com alguma atividade: {len(melhor)}")

# --- TPD-DB: montar no formato da planilha ------------------------------
t = quadro_vazio(len(tpd))
t["smiles"] = tpd["SMILES"].map(limpar)
t["e3"] = tpd["Ligase"].map(limpar)
t["nome"] = tpd["TPD NAME"].map(limpar)
t["alvo"] = tpd["Target Symbol"].map(limpar)
t["uniprot"] = tpd["Target ID"].map(limpar)
t["id_fonte"] = tpd["TPD ID"].map(limpar)
# O TPD-DB não traz MW; a fórmula permite calcular o valor médio exato.
t["mw"] = tpd["Fomula"].map(lambda f: mw_da_formula(f) or pd.NA)
# 'Source' é a referência (patente ou artigo), e vai para a coluna DOI porque
# a planilha não tem coluna de patente. Fica registrado em PROVENIENCIA.
t["doi"] = tpd["Source"].map(limpar)
t["referencia"] = tpd["Source"].map(limpar)
t["fonte"] = ROTULOS["tpd"]

for i, tpd_id in enumerate(tpd["TPD ID"]):
    dados = melhor.get(tpd_id)
    if not dados:
        continue
    for chave, (valor, celula) in dados.items():
        t.iat[i, t.columns.get_loc(chave)] = valor
        if celula is not pd.NA and t.iat[i, t.columns.get_loc("linhagem")] is pd.NA:
            t.iat[i, t.columns.get_loc("linhagem")] = celula

print(f"  -> {len(t)} linhas, {int(t['mw'].notna().sum())} com MW calculado da fórmula")
''')

md(r'''
### PROTAC-DB (`protacdb_20220210.csv`)

Esta base não traz DOI, só `Pubmed`. O identificador entra na coluna DOI como
`PMID:xxxxxxxx` — é a referência que existe, e fica marcada como PMID para não
se confundir com um DOI. `linker_ha` (átomos pesados do linker) preenche
*Comprimeto do Linker*.
''')

code(r'''
db = pd.read_csv(PROTACDB_CSV, dtype=str, keep_default_na=False, encoding="utf-8-sig")
print(f"PROTAC-DB: {len(db)} linhas, {len(db.columns)} colunas")

d = quadro_vazio(len(db))
d["smiles"] = db["PROTAC SMILES"].map(limpar)
d["mw"] = db["MW"].map(limpar)
d["dc50"] = db["Dc50"].map(limpar)
d["dmax"] = db["Dmax"].map(limpar)
d["ic50"] = db["Ic50 of PROTAC"].map(limpar)
d["ec50"] = db["Ec50 of PROTAC Cells"].map(limpar)
d["e3"] = db["E3 Ligase"].map(limpar)
d["linker_tipo"] = db["Linker Type"].map(limpar)
d["linker_comprimento"] = db["linker_ha"].map(limpar)
d["linhagem"] = db["Cells"].map(limpar)
d["uniprot"] = db["Target"].map(limpar)          # nesta base, Target é acessão
d["alvo"] = db["Target"].map(limpar)
d["nome"] = db["PROTAC Name"].map(limpar)
d["id_fonte"] = db["PROTACDB ID"].map(limpar)
pmid = db["Pubmed"].map(limpar)
d["doi"] = pmid.map(lambda p: f"PMID:{p}" if p is not pd.NA else pd.NA)
d["referencia"] = d["doi"]
d["fonte"] = ROTULOS["csv"]

print(f"  -> {len(d)} linhas, {int(d['dc50'].notna().sum())} com DC₅₀")
''')

md(r'''
### PROTACpedia (`protac.xlsx`)

A base mais rica: 89 colunas, incluindo PAMPA e Caco-2. Duas escolhas:

- **IC₅₀/EC₅₀**: há três conjuntos (PROTAC→alvo, PROTAC→E3, atividade celular).
  Uso o **celular** e caio para o PROTAC→alvo quando o celular está vazio — é o
  que as outras duas bases reportam, então fica comparável.
- **Linhagem celular**: não existe como coluna, mas está escrita no texto do
  ensaio (`"Degradation of BRD9 in HeLa cells after 4 h"`). Extraio com regex o
  termo antes de *cells*. Onde o padrão não casa, a célula fica vazia em vez de
  receber o texto inteiro do ensaio.
''')

code(r'''
pedia = pd.read_excel(PROTACPEDIA, dtype=str)
print(f"PROTACpedia: {len(pedia)} linhas, {len(pedia.columns)} colunas")

# "Degradation of BRD9 in HeLa cells after 4 h" -> "HeLa"
_RE_CELULA = re.compile(r"\bin\s+([A-Za-z0-9][A-Za-z0-9\-\./+]*)\s+cells?\b", re.I)


def celula_do_ensaio(texto):
    s = limpar(texto)
    if s is pd.NA:
        return pd.NA
    m = _RE_CELULA.search(str(s))
    return m.group(1) if m else pd.NA


def coalesce(*series):
    """Primeira série com valor, célula a célula."""
    out = series[0].copy()
    for s in series[1:]:
        out = out.where(out.notna(), s)
    return out


p = quadro_vazio(len(pedia))
p["smiles"] = pedia["Smiles"].map(limpar)
p["mw"] = pedia["Molecular Weight"].map(limpar)
p["doi"] = pedia["Article DOI"].map(limpar)
p["referencia"] = p["doi"]
p["dc50"] = pedia["DC50 (nM)"].map(limpar)
p["dmax"] = pedia["Dmax (%)"].map(limpar)
p["ic50"] = coalesce(pedia["IC50 (nM, Cellular activities)"].map(limpar),
                     pedia["IC50 (nM, Protac to Target)"].map(limpar))
p["ec50"] = coalesce(pedia["EC50 (nM, Cellular activities)"].map(limpar),
                     pedia["EC50 (nM, Protac to Target)"].map(limpar))
p["gi50"] = pedia["GI50 (nM, Cellular activities)"].map(limpar)
p["ed50"] = pedia["ED50 (nM, Cellular activities)"].map(limpar)
p["gr50"] = pedia["GR50 (nM, Cellular activities)"].map(limpar)
p["pampa"] = pedia["PAMPA Papp (nm/s, Permeability)"].map(limpar)
p["caco2"] = pedia["Caco-2 A2B Papp (nm/s, Permeability)"].map(limpar)
p["e3"] = pedia["E3 ligase"].map(limpar)
p["alvo"] = pedia["Target"].map(limpar)
p["uniprot"] = pedia["Uniprot"].map(limpar)
p["nome"] = pedia["Name"].map(limpar)
p["id_fonte"] = pedia["Compound ID"].map(limpar)
p["linhagem"] = coalesce(pedia["Assay (DC50/Dmax)"].map(celula_do_ensaio),
                         pedia["Assay (Percent degradation)"].map(celula_do_ensaio))
p["fonte"] = ROTULOS["xlsx"]

print(f"  -> {len(p)} linhas, {int(p['linhagem'].notna().sum())} com linhagem "
      f"extraída do texto do ensaio")
''')

# ===========================================================================
md(r'''
---
## §4 — Juntar e remover duplicatas

A chave é o **SMILES canônico**. Com RDKit, duas grafias da mesma molécula
colapsam corretamente; **sem RDKit**, a comparação é textual e duplicatas entre
bases escapam (cada base escreve o SMILES do seu jeito). O notebook roda nos
dois casos e avisa qual está em uso — se for rodar no `mdtools`, o RDKit está lá.

Ao consolidar, para cada parâmetro fica o **primeiro valor não-vazio** na ordem
PROTACpedia → PROTAC-DB → TPD-DB (a mais anotada primeiro). A coluna `fonte`
registra **todas** as bases em que a molécula apareceu, e a linha ganha a cor
laranja.
''')

code(r'''
try:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    TEM_RDKIT = True
except ImportError:
    TEM_RDKIT = False

print("RDKit:", "sim — deduplicação por SMILES canônico" if TEM_RDKIT else
      "NÃO — deduplicação textual; duplicatas entre bases vão escapar")

_cache = {}


def chave_molecula(smiles):
    s = limpar(smiles)
    if s is pd.NA:
        return pd.NA
    s = str(s)
    if s in _cache:
        return _cache[s]
    chave = s
    if TEM_RDKIT:
        mol = Chem.MolFromSmiles(s)
        if mol is not None:
            chave = Chem.MolToSmiles(mol)
    _cache[s] = chave
    return chave


# Ordem de preferência: a base mais anotada ganha o valor em caso de conflito.
ORDEM = [ROTULOS["xlsx"], ROTULOS["csv"], ROTULOS["tpd"]]
por_fonte = {ROTULOS["xlsx"]: p, ROTULOS["csv"]: d, ROTULOS["tpd"]: t}
bruto = pd.concat([por_fonte[f] for f in ORDEM], ignore_index=True)
print(f"\nlinhas somadas: {len(bruto)}")

if ALVO:
    alvo_norm = str(ALVO).strip().lower()
    mascara = False
    for col in ("alvo", "uniprot", "nome"):
        mascara = mascara | bruto[col].astype("string").str.lower().str.contains(
            re.escape(alvo_norm), na=False)
    bruto = bruto[mascara].copy()
    print(f"após filtro ALVO={ALVO!r}: {len(bruto)} linhas")

bruto["_chave"] = bruto["smiles"].map(chave_molecula)
sem_smiles = bruto[bruto["_chave"].isna()]
com_smiles = bruto[bruto["_chave"].notna()]
print(f"sem SMILES (não dá para deduplicar, passam intactas): {len(sem_smiles)}")

COLS_JUNTAR = {"fonte", "id_fonte", "referencia"}

linhas = []
for _, grupo in com_smiles.groupby("_chave", sort=False):
    if len(grupo) == 1:
        linhas.append(grupo.iloc[0].to_dict())
        continue
    nova = {}
    for col in grupo.columns:
        if col == "_chave":
            continue
        if col in COLS_JUNTAR:
            vistos, unicos = set(), []
            for v in grupo[col]:
                v = limpar(v)
                if v is not pd.NA and str(v) not in vistos:
                    vistos.add(str(v))
                    unicos.append(str(v))
            nova[col] = "; ".join(unicos) if unicos else pd.NA
        else:
            validos = grupo[col].dropna()
            nova[col] = validos.iloc[0] if len(validos) else pd.NA
    nova["_chave"] = grupo.iloc[0]["_chave"]
    linhas.append(nova)

final = pd.concat(
    [pd.DataFrame(linhas), sem_smiles.to_frame().T if isinstance(sem_smiles, pd.Series)
     else sem_smiles],
    ignore_index=True,
).drop(columns=["_chave"], errors="ignore")

# Útil para revisar: quem tem DC₅₀ primeiro, que é o dado mais caro de obter.
final["_sem_dc50"] = final["dc50"].isna()
final = final.sort_values("_sem_dc50", kind="stable").drop(columns=["_sem_dc50"])
final = final.reset_index(drop=True)
final["idx"] = range(1, len(final) + 1)

n_multi = int(final["fonte"].astype("string").str.contains(";", na=False).sum())
print(f"\nlinhas após remover duplicatas: {len(final)}  "
      f"({len(bruto) - len(final)} duplicatas eliminadas)")
print(f"moléculas presentes em mais de uma base: {n_multi}")
print("\ncontribuição por base (uma linha pode citar várias):")
for fonte in ORDEM:
    n = int(final["fonte"].astype("string").str.contains(re.escape(fonte), na=False).sum())
    print(f"  {fonte:<14} {n:>6}")
''')

# ===========================================================================
md(r'''
---
## §5 — O que ficou preenchido

Antes de exportar: quantas linhas têm cada parâmetro. As oito colunas de
farmacocinética saem em 0% porque nenhuma das três bases as traz.
''')

code(r'''
print(f"{len(final)} linhas\n")
print(f"{'parâmetro':<24} {'linhas':>8} {'%':>7}")
print("-" * 42)
for chave in CHAVES:
    if chave == "idx":
        continue
    n = int(final[chave].notna().sum())
    marca = "   " if n else " ··"
    print(f"{marca}{ROTULOS_COL[chave]:<21} {n:>8} {100 * n / max(len(final), 1):>6.1f}%")

vazias = [ROTULOS_COL[c] for c in CHAVES if c != "idx" and final[c].isna().all()]
print(f"\ncolunas sem nenhum dado ({len(vazias)}): {', '.join(vazias) if vazias else '—'}")

display(final[["idx", "fonte", "alvo", "e3", "dc50", "dmax", "ic50", "mw",
               "linhagem", "doi"]].head(15))
''')

# ===========================================================================
md(r'''
---
## §6 — Exportar com as linhas coloridas

Escreve a aba **PAPERS** (cabeçalho na linha 3, dados da linha 4) usando a sua
planilha como molde quando ela existe — **sem nunca sobrescrevê-la**. Saem mais
duas abas: `LEGENDA` (as cores) e `PROVENIENCIA` (de qual base e com qual id
veio cada linha).
''')

code(r'''
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Font, PatternFill

LINHA_INICIAL = 4

if MODELO is not None and MODELO.exists():
    if SAIDA.resolve() == MODELO.resolve():
        raise ValueError("SAIDA igual ao MODELO — a planilha original não deve "
                         "ser sobrescrita. Troque SAIDA na §1.")
    wb = load_workbook(MODELO)
    if "PAPERS" not in wb.sheetnames:
        raise KeyError(f"o modelo não tem a aba 'PAPERS' (abas: {wb.sheetnames})")
    ws = wb["PAPERS"]
    print(f"usando {MODELO.name} como molde")
else:
    wb = Workbook()
    ws = wb.active
    ws.title = "PAPERS"
    for c, chave in enumerate(CHAVES, start=1):
        ws.cell(row=3, column=c, value=ROTULOS_COL[chave]).font = Font(bold=True)
    print("modelo não encontrado — planilha criada do zero")


def cor_da_linha(fonte):
    """Laranja quando a molécula veio de mais de uma base."""
    s = "" if fonte is pd.NA or pd.isna(fonte) else str(fonte)
    return CORES["múltiplas"] if ";" in s else CORES.get(s, "FFFFFF")


preenchimentos = {cor: PatternFill("solid", fgColor=cor)
                  for cor in set(CORES.values())}

for desloc, registro in enumerate(final.to_dict("records")):
    linha = LINHA_INICIAL + desloc
    fill = preenchimentos.get(cor_da_linha(registro.get("fonte")))
    for col, chave in enumerate(CHAVES, start=1):
        valor = registro.get(chave, pd.NA)
        if valor is pd.NA or (not isinstance(valor, str) and pd.isna(valor)):
            valor = None
        elif not isinstance(valor, (int, float, str)):
            valor = str(valor)
        celula = ws.cell(row=linha, column=col, value=valor)
        if fill is not None:
            celula.fill = fill

print(f"aba PAPERS: {len(final)} linhas escritas (linha {LINHA_INICIAL} a "
      f"{LINHA_INICIAL + len(final) - 1})")

# --- LEGENDA -------------------------------------------------------------
if "LEGENDA" in wb.sheetnames:
    del wb["LEGENDA"]
wsl = wb.create_sheet("LEGENDA")
wsl.cell(row=1, column=1, value="Cor").font = Font(bold=True)
wsl.cell(row=1, column=2, value="Origem da molécula").font = Font(bold=True)
wsl.cell(row=1, column=3, value="Linhas").font = Font(bold=True)
descricao = {
    ROTULOS["tpd"]: f"só {ROTULOS['tpd']}",
    ROTULOS["csv"]: f"só {ROTULOS['csv']}",
    ROTULOS["xlsx"]: f"só {ROTULOS['xlsx']}",
    "múltiplas": "mesma molécula em mais de uma base (linha consolidada)",
}
coluna_fonte = final["fonte"].astype("string")
contagem = {
    ROTULOS["tpd"]: int((coluna_fonte == ROTULOS["tpd"]).sum()),
    ROTULOS["csv"]: int((coluna_fonte == ROTULOS["csv"]).sum()),
    ROTULOS["xlsx"]: int((coluna_fonte == ROTULOS["xlsx"]).sum()),
    "múltiplas": int(coluna_fonte.str.contains(";", na=False).sum()),
}
for r, (nome, cor) in enumerate(CORES.items(), start=2):
    wsl.cell(row=r, column=1).fill = preenchimentos[cor]
    wsl.cell(row=r, column=2, value=descricao.get(nome, nome))
    wsl.cell(row=r, column=3, value=contagem.get(nome, 0))
wsl.column_dimensions["A"].width = 6
wsl.column_dimensions["B"].width = 52
wsl.column_dimensions["C"].width = 10

# --- PROVENIENCIA --------------------------------------------------------
if "PROVENIENCIA" in wb.sheetnames:
    del wb["PROVENIENCIA"]
wsp = wb.create_sheet("PROVENIENCIA")
cols_prov = ["idx"] + EXTRAS
for c, nome in enumerate(cols_prov, start=1):
    wsp.cell(row=1, column=c, value=nome).font = Font(bold=True)
for desloc, registro in enumerate(final.to_dict("records")):
    for c, chave in enumerate(cols_prov, start=1):
        valor = registro.get(chave, pd.NA)
        if valor is pd.NA or (not isinstance(valor, str) and pd.isna(valor)):
            valor = None
        elif not isinstance(valor, (int, float, str)):
            valor = str(valor)
        wsp.cell(row=2 + desloc, column=c, value=valor)

wb.save(SAIDA)
print(f"\nsalvo: {SAIDA.resolve()}")
print(f"abas: {', '.join(wb.sheetnames)}")
''')

md(r'''
### Conferência da saída

Relê o arquivo gravado e confere o cabeçalho, a contagem de linhas e as cores.
''')

code(r'''
conf = load_workbook(SAIDA)["PAPERS"]
cabecalho = [conf.cell(row=3, column=c).value for c in range(1, len(CHAVES) + 1)]
esperado = [ROTULOS_COL[k] for k in CHAVES]

if cabecalho == esperado:
    print(f"cabeçalho confere ({len(esperado)} colunas)")
else:
    for i, (saiu, devia) in enumerate(zip(cabecalho, esperado), start=1):
        if saiu != devia:
            print(f"  col {i}: saiu {saiu!r}, esperado {devia!r}")

pintadas = {}
for linha in conf.iter_rows(min_row=LINHA_INICIAL,
                            max_row=LINHA_INICIAL + len(final) - 1,
                            max_col=1):
    cor = linha[0].fill.fgColor.rgb
    pintadas[cor] = pintadas.get(cor, 0) + 1

print(f"linhas de dados: {LINHA_INICIAL} a {LINHA_INICIAL + len(final) - 1}")
print("cores aplicadas (contagem por cor):")
for cor, n in sorted(pintadas.items(), key=lambda kv: -kv[1]):
    nome = next((k for k, v in CORES.items() if str(cor).endswith(v)), str(cor))
    print(f"  {nome:<14} {n:>6}")
''')

md(r'''
---
## Rodar na workstation

Com o VSCode conectado via Remote-SSH, selecione o kernel do env que tem
RDKit (`mdtools`) — sem ele a deduplicação entre bases fica fraca. Depois
ajuste `PASTA` na §1 e `Run All`.

Dependências: `pandas`, `openpyxl` e, recomendado, `rdkit`.

```bash
conda activate mdtools
python -c "import pandas, openpyxl, rdkit; print('ok')"
```
''')

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

destino = (Path(__file__).resolve().parent.parent / "notebooks"
           / "protac_revisao_parametros.ipynb")
destino.parent.mkdir(parents=True, exist_ok=True)
destino.write_text(json.dumps(NOTEBOOK, ensure_ascii=False, indent=1) + "\n",
                   encoding="utf-8")
print(f"{destino}  ({len(CELLS)} células)")
