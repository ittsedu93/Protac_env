#!/usr/bin/env python
"""
Coleta dos parâmetros da aba PAPERS a partir das quatro bases de PROTAC.

Motor por trás de `notebooks/protac_revisao_parametros.ipynb`. Fica num
módulo, e não dentro do notebook, por dois motivos: dá para testar sem abrir
o Jupyter (`python -m unittest discover -s tests`) e dá para rodar na mão
(`python scripts/protac_fontes.py --alvo BRD4 --saida revisao.xlsx`).

Bases cobertas
--------------
- **PROTAC-DB** (ZJU) — literatura curada. É a única que traz, de fábrica,
  quase todas as colunas da planilha: potência, degradação, e na 3.0 também
  os parâmetros farmacocinéticos e de permeabilidade.
- **PROTACpedia** (Weizmann) — curadoria colaborativa. O download completo
  exige cadastro, então esta base entra pelo caminho manual (ver `requer_cadastro`).
- **PROTAC-Bench** — conjunto de benchmark (SMILES + alvo + E3 + rótulo de
  atividade + DC₅₀ quando existe). Preenche poucas colunas, mas é bom para
  conferir cobertura.
- **PROTAC-PatentDB** — compostos extraídos de patentes, com propriedades
  previstas por ADMETlab. Traz MW e permeabilidade para muita molécula que a
  literatura não descreve, e quase nunca traz potência medida.

O princípio de projeto: NÃO dependo do nome exato das colunas
-------------------------------------------------------------
Cada base escreve o mesmo parâmetro com um nome diferente, e muda esses nomes
entre versões — "DC50 (nM)", "DC₅₀", "dc50_nm" e "Degradation DC50" são a
mesma coisa. Se eu fixasse os nomes no código, o notebook quebraria calado na
próxima versão de qualquer uma das bases: a coluna viria vazia e pareceria
"a base não tem esse dado", que é exatamente o erro mais difícil de notar numa
revisão.

Então o mapeamento é por APELIDOS (`APELIDOS`, abaixo): para cada parâmetro da
planilha há uma lista de padrões, comparados contra o nome da coluna de origem
depois de normalizado (minúscula, sem acento, sem pontuação, com os subscritos
₅₀ virando 50). O resultado é auditado: `mapear_colunas` devolve o que casou,
o que casou em mais de um lugar e, principalmente, **as colunas de origem que
ficaram sem destino** — é ali que aparece um parâmetro que a base tem e eu não
reconheci. Para corrigir, acrescente um padrão em `APELIDOS` e rode de novo.

Pelo mesmo motivo os arquivos para baixar são DESCOBERTOS, não fixados: o
módulo lê a página de download de cada base e pega os links .csv/.sdf/.zip que
encontrar (ou usa a API pública, no caso do Figshare e do Hugging Face). Não
há caminho de arquivo chutado no código.

Uso típico
----------
    import protac_fontes as pf

    pedido = pf.Pedido(alvo="BRD4", sinonimos=["BRD-4", "bromodomain 4"])
    bruto  = pf.coletar(pedido, dir_dados="dados/bases")
    pf.exportar_planilha(bruto.tabela, "revisao_BRD4.xlsx")
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import os
import re
import sys
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

import pandas as pd

# ===========================================================================
# 1. O DESTINO: as colunas da aba PAPERS, exatamente como estão na planilha
# ===========================================================================

# Os rótulos abaixo foram copiados da planilha COM OS DEFEITOS: o espaço
# sobrando em "Título Paper ", o "Gl₅₀" com L minúsculo (é GI₅₀, de *growth
# inhibition*), o "Comprimeto" sem o N, o "PAMPA " com espaço no fim. Isso é
# deliberado. Se eu "consertasse" os rótulos aqui, o arquivo exportado
# deixaria de bater coluna a coluna com a planilha que você já usa — e a
# conferência manual, que é o ponto da revisão, ficaria impossível.
#
# A CHAVE (primeiro item) é o nome interno, ASCII e estável, usado em todo o
# resto do código. O RÓTULO (segundo) só aparece no cabeçalho da saída.
PARAMETROS: list[tuple[str, str]] = [
    ("idx",                 "-"),
    ("titulo",              "Título Paper "),
    ("doi",                 "DOI"),
    ("smiles",              "SMILES"),
    ("mw",                  "MW"),
    ("data",                "Data"),
    ("ic50",                "IC₅₀"),
    ("ec50",                "EC₅₀"),
    ("gi50",                "Gl₅₀"),
    ("ed50",                "ED₅₀"),
    ("gr50",                "GR₅₀"),
    ("pampa",               "PAMPA "),
    ("caco2",               "Caco-2 "),
    ("tmax",                "Tmax"),
    ("t_meia",              "T 1/2"),
    ("cmax",                "Cmax"),
    ("auc",                 "AUC"),
    ("vz",                  "Vz"),
    ("vss",                 "Vss"),
    ("cl",                  "CL"),
    ("mrt",                 "MRT"),
    ("dc50",                "DC₅₀"),
    ("dmax",                "Dmax"),
    ("e3",                  "E3 Ligase Usada"),
    ("linker_tipo",         "Tipo de linker"),
    ("linker_comprimento",  "Comprimeto do Linker"),
    ("classe_alvo",         "Classe do Alvo "),
    ("linhagem",            "Linhagem Celular"),
]

CHAVES: list[str] = [k for k, _ in PARAMETROS]
ROTULOS: dict[str, str] = dict(PARAMETROS)

# Colunas que NÃO vão para a aba PAPERS, mas que eu preciso carregar para
# poder filtrar, deduplicar e depois justificar cada linha. Saem na aba
# PROVENIENCIA do arquivo exportado.
CHAVES_PROVENIENCIA: list[str] = [
    "fonte",          # qual base trouxe a linha
    "id_fonte",       # identificador do composto naquela base
    "nome",           # nome/apelido do composto, quando existe
    "alvo",           # proteína-alvo, como a base escreve
    "uniprot",
    "linker_smiles",  # usado para estimar o comprimento do linker
    "patente",        # só no PROTAC-PatentDB
    "atividade",      # rótulo ativo/inativo, só no PROTAC-Bench
    "arquivo_origem",
]

TODAS_AS_CHAVES: list[str] = CHAVES + CHAVES_PROVENIENCIA


# ===========================================================================
# 2. NORMALIZAÇÃO DE NOMES DE COLUNA
# ===========================================================================

_SUBSCRITOS = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")
_SOBRESCRITOS = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹", "0123456789")


def normalizar_nome(nome: Any) -> str:
    """Reduz um nome de coluna à forma em que dá para comparar.

    'DC₅₀ (nM)', 'DC50 [nM]' e 'dc50_nm' viram os três 'dc50 nm'. Sem isso,
    cada base precisaria de uma lista de apelidos própria.
    """
    if nome is None:
        return ""
    s = str(nome).translate(_SUBSCRITOS).translate(_SOBRESCRITOS)
    # NFKD separa o acento da letra; descarto as combinantes e fico com ASCII
    s = unicodedata.normalize("NFKD", s)
    s = "".join(c for c in s if not unicodedata.combining(c))
    s = re.sub(r"[^A-Za-z0-9]+", " ", s.lower())
    return re.sub(r"\s+", " ", s).strip()


# ===========================================================================
# 3. OS APELIDOS: como cada base pode ter chamado cada parâmetro
# ===========================================================================

# Padrões de regex, aplicados com re.search contra o nome JÁ NORMALIZADO.
# Acrescentar um apelido aqui é a única coisa que precisa ser feita quando uma
# base renomeia uma coluna, ou quando a auditoria mostra uma coluna de origem
# que ficou sem destino.
APELIDOS: dict[str, list[str]] = {
    # --- identificação e referência bibliográfica ---------------------------
    "doi": [r"\bdoi\b", r"article doi", r"reference doi", r"publication doi"],
    "titulo": [r"\btitle\b", r"article title", r"paper title", r"titulo",
               r"publication title", r"\bpaper\b", r"citation"],
    "data": [r"^date", r"\bdate\b", r"\byear\b", r"publication year",
             r"pub year", r"article date", r"patent date", r"publication date",
             r"priority date", r"\bdata\b"],

    # --- estrutura e propriedade física -----------------------------------
    # 'smiles' tem de vir antes de qualquer padrão genérico de estrutura, e o
    # do linker é resolvido antes (ver ORDEM_RESOLUCAO) para não roubar esta.
    "smiles": [r"\bsmiles\b", r"canonical smiles", r"isomeric smiles",
               r"structure smiles", r"protac smiles", r"compound smiles"],
    "mw": [r"^mw$", r"^mw\b", r"molecular weight", r"\bmol wt\b", r"molwt",
           r"\bmw g mol\b", r"exact mass", r"average mass", r"massa molar"],

    # --- potência ---------------------------------------------------------
    # Note que todos aceitam o subscrito, porque normalizar_nome já o converteu.
    "ic50": [r"\bic50\b", r"\bic 50\b", r"inhibition ic50", r"\bpic50\b"],
    "ec50": [r"\bec50\b", r"\bec 50\b", r"\bpec50\b"],
    "gi50": [r"\bgi50\b", r"\bgl50\b", r"\bgi 50\b", r"growth inhibition"],
    "ed50": [r"\bed50\b", r"\bed 50\b", r"effective dose"],
    "gr50": [r"\bgr50\b", r"\bgr 50\b", r"growth rate inhibition"],

    # --- degradação (o que é específico de PROTAC) ------------------------
    "dc50": [r"\bdc50\b", r"\bdc 50\b", r"degradation dc50", r"\bdc50 nm\b",
             r"half maximal degradation"],
    "dmax": [r"\bdmax\b", r"\bd max\b", r"maximum degradation",
             r"max degradation", r"degradation max"],

    # --- permeabilidade ---------------------------------------------------
    "pampa": [r"\bpampa\b", r"\bpapp\b.*pampa", r"pampa permeability"],
    "caco2": [r"\bcaco 2\b", r"\bcaco2\b", r"caco 2 permeability",
              r"\bpapp\b.*caco"],

    # --- farmacocinética --------------------------------------------------
    # Os limites de palavra importam muito aqui: sem eles, 'cl' casaria dentro
    # de 'cell line' e 'auc' dentro de 'aucs'.
    "tmax": [r"^t ?max\b", r"\btmax\b", r"time to cmax"],
    "t_meia": [r"\bt1 2\b", r"\bt 1 2\b", r"\bt half\b", r"half life",
               r"\bhalflife\b", r"\bthalf\b", r"\bt12\b"],
    "cmax": [r"^c ?max\b", r"\bcmax\b", r"peak concentration"],
    "auc": [r"^auc\b", r"\bauc\b", r"\bauc0\b", r"area under the curve"],
    "vz": [r"^vz\b", r"\bvz\b", r"\bv z\b", r"volume of distribution vz"],
    "vss": [r"^vss\b", r"\bvss\b", r"\bv ss\b", r"steady state volume"],
    "cl": [r"^cl$", r"^cl\b", r"\bclearance\b", r"\bcl ml min kg\b",
           r"plasma clearance", r"total clearance"],
    "mrt": [r"^mrt\b", r"\bmrt\b", r"mean residence time"],

    # --- anatomia do PROTAC ----------------------------------------------
    "e3": [r"e3 ligase", r"\be3\b", r"ligase", r"e3 ligand name",
           r"\brecruiter\b", r"\bcrbn vhl\b", r"e3 type"],
    "linker_tipo": [r"linker type", r"tipo de linker", r"linker class",
                    r"linker category", r"linker chemistry", r"tipo linker"],
    "linker_comprimento": [r"linker length", r"comprimento do linker",
                           r"linker n atoms", r"linker atoms",
                           r"linker heavy atoms", r"n linker"],
    "linker_smiles": [r"linker smiles", r"smiles linker", r"linker structure"],

    # --- alvo e ensaio ----------------------------------------------------
    # 'classe_alvo' resolve antes de 'alvo' para não perder 'Target Class'
    # para o padrão genérico de 'target'.
    "classe_alvo": [r"target class", r"target family", r"protein class",
                    r"protein family", r"classe do alvo", r"target type",
                    r"\bfamily\b"],
    "linhagem": [r"cell line", r"cell lines", r"\bcells\b", r"cell type",
                 r"linhagem", r"\bcell\b"],
    "alvo": [r"^target$", r"^target\b", r"\btarget name\b", r"\btarget protein\b",
             r"\bprotein\b", r"\bgene\b", r"\balvo\b", r"\bpoi\b"],
    "uniprot": [r"\buniprot\b", r"uniprot id", r"uniprot accession"],

    # --- identificação interna da base -----------------------------------
    "id_fonte": [r"^id$", r"compound id", r"protac id", r"\bcid\b",
                 r"^index$", r"entry id", r"\bidentifier\b"],
    "nome": [r"^name$", r"compound name", r"protac name", r"\balias\b",
             r"common name"],
    "patente": [r"\bpatent\b", r"patent number", r"patent family",
                r"publication number", r"\bassignee\b"],
    "atividade": [r"^active$", r"\bactivity\b", r"activity label",
                  r"\blabel\b", r"is active"],
}

# A ordem em que os parâmetros disputam as colunas de origem. Uma coluna de
# origem é reivindicada por UM parâmetro só (o primeiro que casar), então o
# específico tem de vir antes do genérico. Os três pares que realmente
# importam: linker_smiles antes de smiles, classe_alvo antes de alvo,
# linhagem antes de alvo (por causa de 'target cell line').
ORDEM_RESOLUCAO: list[str] = [
    "doi", "titulo", "data",
    "linker_smiles", "linker_tipo", "linker_comprimento",
    "smiles", "mw",
    "dc50", "dmax",
    "ic50", "ec50", "gi50", "ed50", "gr50",
    "pampa", "caco2",
    "tmax", "t_meia", "cmax", "auc", "vz", "vss", "cl", "mrt",
    "e3", "classe_alvo", "linhagem", "uniprot", "alvo",
    "id_fonte", "nome", "patente", "atividade",
]


def mapear_colunas(colunas: Iterable[Any]) -> tuple[dict[str, str], dict]:
    """Decide qual coluna de origem alimenta cada parâmetro da planilha.

    Devolve `(mapa, auditoria)`:
      - `mapa[chave] = nome_da_coluna_de_origem`
      - `auditoria["alternativas"][chave]` = outras colunas que também casaram
        (ficaram de fora porque só uma pode alimentar a célula)
      - `auditoria["sem_destino"]` = colunas de origem que nenhum parâmetro
        reivindicou. **É a lista que importa revisar**: um parâmetro que a base
        tem e eu não reconheci aparece aqui, e não como uma coluna vazia.
      - `auditoria["sem_origem"]` = parâmetros da planilha que esta base não tem

    Uma coluna de origem só é usada uma vez; a ordem da disputa é
    ORDEM_RESOLUCAO.
    """
    originais = [c for c in colunas]
    normalizados = {c: normalizar_nome(c) for c in originais}

    mapa: dict[str, str] = {}
    alternativas: dict[str, list[str]] = {}
    usadas: set[Any] = set()

    for chave in ORDEM_RESOLUCAO:
        padroes = APELIDOS.get(chave, [])
        casaram: list[Any] = []
        for coluna in originais:
            if coluna in usadas:
                continue
            alvo_norm = normalizados[coluna]
            if not alvo_norm:
                continue
            if any(re.search(p, alvo_norm) for p in padroes):
                casaram.append(coluna)
        if casaram:
            # O primeiro na ordem da base vira o escolhido. Em geral a coluna
            # mais "crua" vem antes das derivadas, e é a que eu quero.
            mapa[chave] = casaram[0]
            usadas.add(casaram[0])
            if len(casaram) > 1:
                alternativas[chave] = [str(c) for c in casaram[1:]]

    auditoria = {
        "alternativas": alternativas,
        "sem_destino": [str(c) for c in originais if c not in usadas],
        "sem_origem": [k for k in CHAVES if k not in mapa and k != "idx"],
        "mapa": {k: str(v) for k, v in mapa.items()},
    }
    return mapa, auditoria


def aplicar_mapa(df: pd.DataFrame, mapa: dict[str, str]) -> pd.DataFrame:
    """Reescreve um DataFrame qualquer no vocabulário da planilha."""
    saida = pd.DataFrame(index=df.index)
    for chave in TODAS_AS_CHAVES:
        if chave in mapa and mapa[chave] in df.columns:
            saida[chave] = df[mapa[chave]]
        else:
            saida[chave] = pd.NA
    return saida


# ===========================================================================
# 4. VALORES: unidades, qualificadores e faixas
# ===========================================================================

# As bases escrevem potência de todo jeito: '13 nM', '0.013 uM', '>10 μM',
# '~5 nM', '1.2-3.4 nM', '13' (sem unidade, com a unidade no cabeçalho).
# Guardo o número E o texto: a planilha é lida por gente, e '>10 μM' diz mais
# do que '10000'.

_FATOR_PARA_NM = {
    "m": 1e9, "mol/l": 1e9, "molar": 1e9,
    "mm": 1e6, "mmol/l": 1e6,
    "um": 1e3, "umol/l": 1e3,
    "nm": 1.0, "nmol/l": 1.0,
    "pm": 1e-3, "pmol/l": 1e-3,
}

# O micro aparece como U+00B5 (MICRO SIGN) ou U+03BC (mu grego), e os dois
# somem se eu filtrar só [A-Za-z0-9] — aí 'μM' vira 'M' e o valor sai um
# milhão de vezes maior, calado. Traduzo antes de filtrar.
_MICRO = str.maketrans({"\u00b5": "u", "\u03bc": "u"})


def normalizar_unidade(unidade: Any) -> str:
    """'μM', 'uM', 'µmol/L' -> 'um', 'um', 'umol/l'."""
    if unidade is None:
        return ""
    s = str(unidade).strip().translate(_MICRO)
    s = unicodedata.normalize("NFKD", s)
    s = "".join(c for c in s if not unicodedata.combining(c))
    return re.sub(r"[^a-z0-9/]+", "", s.lower())

_RE_NUMERO = r"[-+]?\d+(?:[.,]\d+)?(?:[eE][-+]?\d+)?"
_RE_QUANTIDADE = re.compile(
    r"^\s*(?P<qual>[<>≤≥~≈]{0,2})\s*"
    r"(?P<v1>" + _RE_NUMERO + r")"
    r"(?:\s*(?:-|–|—|a|to|±)\s*(?P<v2>" + _RE_NUMERO + r"))?"
    r"\s*(?P<unidade>[^\s,;()]*)",
)


def _para_float(texto: str | None) -> float | None:
    if texto is None:
        return None
    try:
        return float(str(texto).replace(",", "."))
    except (TypeError, ValueError):
        return None


@dataclass
class Quantidade:
    """Um valor medido, decomposto: '>10 μM' -> qual='>', v1=10, unidade='μM'."""
    bruto: str
    qualificador: str = ""
    v1: float | None = None
    v2: float | None = None
    unidade: str = ""

    @property
    def vazia(self) -> bool:
        return self.v1 is None

    def em_nm(self) -> float | None:
        """Converte para nM quando a unidade é de concentração; senão, None."""
        if self.v1 is None:
            return None
        fator = _FATOR_PARA_NM.get(normalizar_unidade(self.unidade))
        if fator is None:
            return None
        return self.v1 * fator

    def texto(self) -> str:
        """Forma canônica para escrever na célula, preservando o sentido."""
        if self.v1 is None:
            return self.bruto.strip()
        def fmt(x: float) -> str:
            if x == int(x) and abs(x) < 1e15:
                return str(int(x))
            return f"{x:g}"
        nucleo = fmt(self.v1)
        if self.v2 is not None:
            nucleo += f"–{fmt(self.v2)}"
        return f"{self.qualificador}{nucleo}{(' ' + self.unidade) if self.unidade else ''}".strip()


def ler_quantidade(valor: Any) -> Quantidade:
    """Decompõe uma célula de potência/PK em número, qualificador e unidade."""
    if valor is None or (isinstance(valor, float) and pd.isna(valor)):
        return Quantidade(bruto="")
    if isinstance(valor, (int, float)) and not isinstance(valor, bool):
        return Quantidade(bruto=str(valor), v1=float(valor))
    bruto = str(valor).strip()
    if not bruto or normalizar_nome(bruto) in {"na", "nan", "nd", "none", "null", ""}:
        return Quantidade(bruto=bruto)
    m = _RE_QUANTIDADE.match(bruto)
    if not m:
        return Quantidade(bruto=bruto)
    qual = m.group("qual") or ""
    qual = qual.replace("≤", "<=").replace("≥", ">=").replace("≈", "~")
    unidade = (m.group("unidade") or "").strip()
    # '13 nM in MV4-11' -> a unidade é 'nM'; '13 %' -> '%'. Se o que sobrou não
    # parece unidade, descarto em vez de inventar.
    if unidade and not re.fullmatch(r"[%µμA-Za-z/·\.]+", unidade):
        unidade = ""
    return Quantidade(
        bruto=bruto,
        qualificador=qual,
        v1=_para_float(m.group("v1")),
        v2=_para_float(m.group("v2")),
        unidade=unidade,
    )


def ler_dmax(valor: Any) -> Quantidade:
    """Dmax é porcentagem, e as bases discordam na escala.

    '95%' e '95' querem dizer a mesma coisa; '0.95' quer dizer o mesmo em
    fração. Converto fração para porcentagem quando o número é <= 1 e não há
    sinal de '%' — é heurística, e está documentada aqui porque muda o número
    que vai para a planilha.
    """
    q = ler_quantidade(valor)
    if q.v1 is None:
        return q
    tem_pct = "%" in q.bruto
    if not tem_pct and 0 < q.v1 <= 1.0:
        q.v1 *= 100.0
        if q.v2 is not None:
            q.v2 *= 100.0
    q.unidade = "%"
    return q


def limpar_texto(valor: Any) -> Any:
    """Normaliza células de texto e transforma os vários 'vazios' em NA."""
    if valor is None:
        return pd.NA
    if isinstance(valor, float) and pd.isna(valor):
        return pd.NA
    s = str(valor).strip()
    if not s or normalizar_nome(s) in {"na", "nan", "nd", "none", "null", "n a", "unknown"}:
        return pd.NA
    return re.sub(r"\s+", " ", s)


# Parâmetros tratados como concentração/medida (passam por ler_quantidade) e
# os tratados como texto livre.
CHAVES_QUANTIDADE = ["ic50", "ec50", "gi50", "ed50", "gr50", "dc50",
                     "pampa", "caco2", "tmax", "t_meia", "cmax", "auc",
                     "vz", "vss", "cl", "mrt", "mw"]


def padronizar_valores(df: pd.DataFrame, unidades_fixas: bool = False) -> pd.DataFrame:
    """Passa cada coluna pelo tratamento que o seu tipo pede.

    Com `unidades_fixas=True`, as potências e o DC₅₀ saem convertidos para nM
    (número puro), o que serve para ordenar e comparar. Por padrão sai o texto
    canônico com unidade, que é o que faz sentido ler na planilha.
    """
    out = df.copy()
    for chave in CHAVES_QUANTIDADE:
        if chave not in out.columns:
            continue
        quantidades = [ler_quantidade(v) for v in out[chave]]
        if unidades_fixas and chave in {"ic50", "ec50", "gi50", "ed50", "gr50", "dc50"}:
            out[chave] = [q.em_nm() if not q.vazia else pd.NA for q in quantidades]
        else:
            out[chave] = [q.texto() if not q.vazia else pd.NA for q in quantidades]
    if "dmax" in out.columns:
        out["dmax"] = [ (q.texto() if not q.vazia else pd.NA)
                        for q in (ler_dmax(v) for v in out["dmax"]) ]
    for chave in out.columns:
        if chave not in CHAVES_QUANTIDADE and chave != "dmax":
            out[chave] = [limpar_texto(v) for v in out[chave]]
    return out


# ===========================================================================
# 5. LEITURA DE ARQUIVOS (csv, tsv, xlsx, json, parquet, sdf, zip, gz)
# ===========================================================================

def _abrir_bytes(caminho: Path) -> bytes:
    dados = caminho.read_bytes()
    if caminho.suffix.lower() == ".gz":
        return gzip.decompress(dados)
    return dados


def ler_sdf(texto: str) -> pd.DataFrame:
    """Lê só o BLOCO DE PROPRIEDADES de um SDF, sem RDKit.

    O molblock (as coordenadas) é ignorado de propósito: para esta revisão o
    que interessa são as propriedades anotadas (`> <Smiles>`, `> <DC50>` ...),
    e assim o módulo não passa a depender do RDKit. A consequência está
    documentada: se o SDF não tiver uma propriedade com o SMILES, a coluna
    SMILES vem vazia — nesse caso prefira o CSV da mesma base.
    """
    registros: list[dict[str, str]] = []
    for bloco in texto.split("$$$$"):
        if not bloco.strip():
            continue
        registro: dict[str, str] = {}
        # A primeira linha do molblock é o nome do composto.
        linhas = bloco.splitlines()
        if linhas and linhas[0].strip():
            registro["_sdf_name"] = linhas[0].strip()
        for m in re.finditer(r"^>\s*<([^>]+)>[^\n]*\n(.*?)(?=\n>|\n\s*\n|\Z)",
                             bloco, flags=re.MULTILINE | re.DOTALL):
            registro[m.group(1).strip()] = m.group(2).strip()
        if len(registro) > 1:
            registros.append(registro)
    return pd.DataFrame(registros)


def ler_tabela(caminho: str | Path) -> pd.DataFrame:
    """Lê uma tabela de qualquer um dos formatos que estas bases distribuem."""
    caminho = Path(caminho)
    nome = caminho.name.lower()
    sufixos = [s.lower() for s in caminho.suffixes]

    if ".zip" in sufixos:
        quadros = []
        with zipfile.ZipFile(caminho) as z:
            for membro in z.namelist():
                if membro.endswith("/"):
                    continue
                if not re.search(r"\.(csv|tsv|txt|sdf|json|xlsx)$", membro, re.I):
                    continue
                with z.open(membro) as fh:
                    bruto = fh.read()
                quadro = _ler_bytes(bruto, membro)
                if quadro is not None and not quadro.empty:
                    quadro["arquivo_origem"] = f"{caminho.name}:{membro}"
                    quadros.append(quadro)
        if not quadros:
            return pd.DataFrame()
        return pd.concat(quadros, ignore_index=True) if len(quadros) > 1 else quadros[0]

    if ".parquet" in sufixos:
        return pd.read_parquet(caminho)
    if ".xlsx" in sufixos or ".xls" in sufixos:
        return pd.read_excel(caminho)

    quadro = _ler_bytes(_abrir_bytes(caminho), nome)
    return quadro if quadro is not None else pd.DataFrame()


def _ler_bytes(bruto: bytes, nome: str) -> pd.DataFrame | None:
    """Despacha pelo nome do arquivo, já com os bytes em mão."""
    nome = nome.lower()
    if nome.endswith(".gz"):
        bruto = gzip.decompress(bruto)
        nome = nome[:-3]

    if nome.endswith(".sdf"):
        return ler_sdf(bruto.decode("utf-8", errors="replace"))
    if nome.endswith(".parquet"):
        return pd.read_parquet(io.BytesIO(bruto))
    if nome.endswith((".xlsx", ".xls")):
        return pd.read_excel(io.BytesIO(bruto))
    if nome.endswith(".jsonl") or nome.endswith(".ndjson"):
        return pd.read_json(io.BytesIO(bruto), lines=True)
    if nome.endswith(".json"):
        dados = json.loads(bruto.decode("utf-8", errors="replace"))
        if isinstance(dados, dict):
            # Formato comum: {"data": [...]} ou {"results": [...]}
            for chave in ("data", "results", "records", "entries", "rows"):
                if isinstance(dados.get(chave), list):
                    dados = dados[chave]
                    break
        return pd.json_normalize(dados)

    # csv/tsv/txt: deixo o pandas farejar o separador, que varia entre bases
    texto = bruto.decode("utf-8", errors="replace")
    separador = "\t" if nome.endswith((".tsv", ".tab")) else None
    try:
        return pd.read_csv(io.StringIO(texto), sep=separador, engine="python",
                           on_bad_lines="skip")
    except Exception:
        return pd.read_csv(io.StringIO(texto), sep=None, engine="python",
                           on_bad_lines="skip")


# ===========================================================================
# 6. AS BASES: descoberta dos arquivos e download
# ===========================================================================

@dataclass
class Fonte:
    """Uma base de dados e como chegar nos arquivos dela."""
    nome: str
    descricao: str
    pagina: str                        # página humana, para citar e para abrir
    paginas_download: tuple[str, ...] = ()   # onde varrer por links de arquivo
    api_figshare: str | None = None    # id do artigo no Figshare
    repo_hf: str | None = None         # 'owner/nome' de dataset no Hugging Face
    extensoes: tuple[str, ...] = (".csv", ".tsv", ".sdf", ".zip", ".xlsx",
                                  ".json", ".parquet")
    padroes_locais: tuple[str, ...] = ()
    requer_cadastro: bool = False
    observacao: str = ""


# Nenhuma URL de arquivo está fixada aqui de propósito: os caminhos diretos
# mudam a cada versão das bases, e um caminho errado no código viraria uma
# coluna vazia silenciosa. O que está fixado é o ENDEREÇO ESTÁVEL de cada base
# (a página de download, o id do Figshare, o repo do Hugging Face) e a varredura
# de links faz o resto.
FONTES: dict[str, Fonte] = {
    "PROTAC-DB": Fonte(
        nome="PROTAC-DB",
        descricao="Literatura curada (ZJU). A base mais completa para esta "
                  "planilha: a versão 3.0 acrescentou os parâmetros "
                  "farmacocinéticos e de permeabilidade.",
        pagina="http://cadd.zju.edu.cn/protacdb/",
        paginas_download=("http://cadd.zju.edu.cn/protacdb/downloads",
                          "http://cadd.zju.edu.cn/protacdb/download",
                          "http://cadd.zju.edu.cn/protacdb/"),
        padroes_locais=("*protac*db*.csv", "*protac*db*.sdf", "protac.csv",
                        "protac*.csv", "protac*.sdf", "*protacdb*.zip"),
        observacao="Baixe 'protac' (a tabela dos PROTACs completos). As "
                   "tabelas warhead/linker/e3_ligand são dos fragmentos e não "
                   "preenchem as colunas de atividade.",
    ),
    "PROTACpedia": Fonte(
        nome="PROTACpedia",
        descricao="Curadoria colaborativa (Weizmann). Boa em 'Linker Type' e "
                  "em anotação de ensaio, que as outras quase não têm.",
        pagina="https://protacpedia.weizmann.ac.il/ptcb/main",
        paginas_download=("https://protacpedia.weizmann.ac.il/ptcb/main",),
        padroes_locais=("*protacpedia*.csv", "*protacpedia*.sdf",
                        "*protacpedia*.xlsx", "*ptcb*.csv"),
        requer_cadastro=True,
        observacao="O download completo exige cadastro (gratuito) no site. "
                   "Faça o login, baixe o CSV/SDF e salve em dir_dados — o "
                   "notebook o encontra pelo nome. Não automatizo o login: "
                   "seria guardar a sua senha num notebook.",
    ),
    "PROTAC-Bench": Fonte(
        nome="PROTAC-Bench",
        descricao="Conjunto de benchmark de degradação (SMILES, alvo, E3, "
                  "rótulo de atividade, DC₅₀ quando existe).",
        pagina="https://huggingface.co/datasets/ThorKl/protac-bench",
        repo_hf="ThorKl/protac-bench",
        padroes_locais=("*protac*bench*.csv", "*protac*bench*.parquet",
                        "*protac_bench*.json*"),
        observacao="Atenção ao nome: 'PROTAC-Bench' também é usado para "
                   "benchmarks de DOCKING de complexo ternário, que são "
                   "estruturas (PDB) e não preenchem esta planilha. Se for "
                   "esse que você quer, troque o repo em FONTES.",
    ),
    "PROTAC-PatentDB": Fonte(
        nome="PROTAC-PatentDB",
        descricao="Compostos extraídos de patentes, com propriedades previstas "
                  "por ADMETlab. Preenche MW e permeabilidade em larga escala; "
                  "potência medida, quase nunca.",
        pagina="https://www.nature.com/articles/s41597-025-06136-9",
        paginas_download=("https://protacpatentdb.com/",),
        api_figshare="29351321",
        padroes_locais=("*patent*db*.csv", "*protac*patent*.csv",
                        "*patent*.xlsx", "*patentdb*.zip"),
        observacao="Os valores de permeabilidade/PK desta base são PREVISTOS, "
                   "não medidos. A coluna 'fonte' na aba PROVENIENCIA registra "
                   "isso — não misture com medida experimental sem marcar.",
    ),
}


def _requisicao(url: str, aceita: str = "*/*", tempo: float = 60.0) -> bytes:
    """GET simples, respeitando o proxy do ambiente (urllib lê HTTPS_PROXY)."""
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "protac-fontes/1.0 (revisao academica)",
            "Accept": aceita,
        },
    )
    with urllib.request.urlopen(req, timeout=tempo) as resp:
        return resp.read()


def descobrir_links(pagina: str, extensoes: Sequence[str]) -> list[str]:
    """Varre uma página de download por links de arquivo de dados.

    É isto que evita ter caminho de arquivo chutado no código: em vez de supor
    '/static/protac.csv', eu leio a página e pego o que ela oferece.
    """
    try:
        html = _requisicao(pagina, aceita="text/html").decode("utf-8", errors="replace")
    except Exception:
        return []
    achados: list[str] = []
    for m in re.finditer(r"""(?:href|src)\s*=\s*["']([^"']+)["']""", html, re.I):
        alvo = m.group(1)
        if any(alvo.lower().split("?")[0].endswith(e) for e in extensoes):
            achados.append(urllib.parse.urljoin(pagina, alvo))
    # sem perder a ordem, remove repetido
    vistos: set[str] = set()
    return [u for u in achados if not (u in vistos or vistos.add(u))]


def descobrir_figshare(id_artigo: str) -> list[tuple[str, str]]:
    """Lista os arquivos de um artigo do Figshare pela API pública.

    Devolve [(nome_do_arquivo, url_de_download)].
    """
    url = f"https://api.figshare.com/v2/articles/{id_artigo}/files"
    try:
        dados = json.loads(_requisicao(url, aceita="application/json"))
    except Exception:
        return []
    saida = []
    for item in dados if isinstance(dados, list) else []:
        nome = item.get("name")
        link = item.get("download_url")
        if nome and link:
            saida.append((nome, link))
    return saida


def descobrir_hf(repo: str) -> list[tuple[str, str]]:
    """Lista os arquivos de um dataset público do Hugging Face."""
    url = f"https://huggingface.co/api/datasets/{repo}/tree/main?recursive=1"
    try:
        dados = json.loads(_requisicao(url, aceita="application/json"))
    except Exception:
        return []
    saida = []
    for item in dados if isinstance(dados, list) else []:
        if item.get("type") != "file":
            continue
        caminho = item.get("path", "")
        if re.search(r"\.(csv|tsv|parquet|json|jsonl)$", caminho, re.I):
            saida.append((caminho.replace("/", "__"),
                          f"https://huggingface.co/datasets/{repo}/resolve/main/{caminho}"))
    return saida


def arquivos_locais(fonte: Fonte, dir_dados: Path) -> list[Path]:
    """Procura, em dir_dados, arquivos que pareçam ser desta base.

    O cuidado com a pasta de outra base não é teórico: o padrão 'protac*.csv'
    do PROTAC-DB casa com 'protac_bench.csv'. Sem o descarte abaixo, uma base
    comeria o arquivo da outra e a coluna 'fonte' mentiria.
    """
    if not dir_dados.exists():
        return []
    outras = {n for n in FONTES if n != fonte.nome}

    def de_outra_base(caminho: Path) -> bool:
        return any(parte in outras for parte in caminho.parts)

    achados: list[Path] = []
    for padrao in fonte.padroes_locais:
        achados.extend(p for p in sorted(dir_dados.rglob(padrao))
                       if not de_outra_base(p))
    # a subpasta com o nome da base também vale, qualquer que seja o arquivo
    pasta = dir_dados / fonte.nome
    if pasta.exists():
        for p in sorted(pasta.rglob("*")):
            if p.is_file() and re.search(r"\.(csv|tsv|sdf|xlsx|json|jsonl|parquet|zip|gz)$",
                                         p.name, re.I):
                achados.append(p)
    vistos: set[Path] = set()
    return [p for p in achados if not (p in vistos or vistos.add(p))]


def obter_arquivos(fonte: Fonte, dir_dados: str | Path,
                   baixar: bool = True, forcar: bool = False) -> tuple[list[Path], list[str]]:
    """Garante os arquivos da base em disco. Devolve (arquivos, registro).

    Ordem: (1) o que já está em dir_dados; (2) se `baixar`, descobre e baixa.
    Nunca estoura exceção por falta de rede — devolve lista vazia e o registro
    explica o que faltou, porque uma base indisponível não deve derrubar a
    coleta das outras três.
    """
    dir_dados = Path(dir_dados)
    destino = dir_dados / fonte.nome
    registro: list[str] = []

    if not forcar:
        locais = arquivos_locais(fonte, dir_dados)
        if locais:
            registro.append(f"{fonte.nome}: {len(locais)} arquivo(s) já em disco "
                            f"({', '.join(p.name for p in locais[:4])}"
                            f"{', ...' if len(locais) > 4 else ''})")
            return locais, registro

    if not baixar:
        registro.append(f"{fonte.nome}: nada em disco e download desligado.")
        return [], registro

    if fonte.requer_cadastro:
        registro.append(
            f"{fonte.nome}: exige cadastro, não baixo automaticamente. "
            f"Baixe em {fonte.pagina} e salve em {destino}/"
        )
        return [], registro

    candidatos: list[tuple[str, str]] = []
    if fonte.api_figshare:
        candidatos.extend(descobrir_figshare(fonte.api_figshare))
    if fonte.repo_hf:
        candidatos.extend(descobrir_hf(fonte.repo_hf))
    for pagina in fonte.paginas_download:
        for url in descobrir_links(pagina, fonte.extensoes):
            candidatos.append((Path(urllib.parse.urlparse(url).path).name or "dados.csv", url))

    if not candidatos:
        registro.append(
            f"{fonte.nome}: não consegui descobrir nenhum arquivo para baixar "
            f"(rede bloqueada, ou a página mudou). Baixe na mão em "
            f"{fonte.pagina} e salve em {destino}/"
        )
        return [], registro

    destino.mkdir(parents=True, exist_ok=True)
    baixados: list[Path] = []
    for nome, url in candidatos:
        caminho = destino / nome
        if caminho.exists() and not forcar:
            baixados.append(caminho)
            continue
        try:
            caminho.write_bytes(_requisicao(url))
            baixados.append(caminho)
            registro.append(f"{fonte.nome}: baixado {nome}")
        except Exception as erro:
            registro.append(f"{fonte.nome}: falhou {url} ({type(erro).__name__}: {erro})")
    if not baixados:
        registro.append(f"{fonte.nome}: descobri {len(candidatos)} link(s), "
                        f"mas nenhum download completou.")
    return baixados, registro


# ===========================================================================
# 7. O PEDIDO (o "input" que você preenche) E O FILTRO
# ===========================================================================

@dataclass
class Pedido:
    """Tudo o que você precisa preencher para rodar a coleta.

    `alvo` é o único campo que costuma bastar. Os sinônimos existem porque as
    bases escrevem a mesma proteína de formas diferentes ('PCSK9', 'proprotein
    convertase subtilisin/kexin type 9'): sem eles, um filtro por 'PCSK9'
    descarta linhas que são do seu alvo.
    """
    alvo: str | None = None
    sinonimos: list[str] = field(default_factory=list)
    e3: str | None = None                      # 'VHL', 'CRBN', ...
    uniprot: str | None = None
    linhagem: str | None = None
    texto_livre: str | None = None             # busca em qualquer coluna
    bases: list[str] = field(default_factory=lambda: list(FONTES))
    so_com_degradacao: bool = False            # exige DC₅₀ ou Dmax preenchido
    unidades_fixas: bool = False               # potências em nM, número puro
    consolidar_duplicatas: bool = True

    def termos_alvo(self) -> list[str]:
        """Alvo + sinônimos + o UniProt, que para algumas bases é o único jeito.

        O PROTAC-Bench, por exemplo, não traz o nome do gene: só a acessão.
        Incluir o UniProt aqui faz um pedido só servir para as quatro bases.
        """
        termos = [t for t in ([self.alvo] + list(self.sinonimos) + [self.uniprot]) if t]
        return [normalizar_nome(t) for t in termos if normalizar_nome(t)]


def _contem(serie: pd.Series, termos: Sequence[str]) -> pd.Series:
    """Série booleana: a célula contém algum dos termos (já normalizados)?"""
    normalizada = serie.astype("string").fillna("").map(normalizar_nome)
    resultado = pd.Series(False, index=serie.index)
    for termo in termos:
        resultado |= normalizada.str.contains(re.escape(termo), regex=True, na=False)
    return resultado


def filtrar(df: pd.DataFrame, pedido: Pedido) -> pd.DataFrame:
    """Aplica o pedido. Cada critério é um E; cada critério casa por 'contém'.

    O alvo é procurado em `alvo`, `uniprot`, `nome` e `titulo` — e não só em
    `alvo` — porque em várias bases o alvo só aparece no nome do composto ou no
    título do artigo.

    `Pedido.uniprot` entra como MAIS UM NOME do alvo (em OU com o nome do
    gene), não como critério separado. A diferença importa: como critério em E,
    informar o UniProt para resgatar o PROTAC-Bench zerava as outras três
    bases, que não trazem a coluna — e zeravam caladas, parecendo ausência de
    dado. Para filtrar só por acessão, deixe `alvo=None` e preencha `uniprot`.
    """
    if df.empty:
        return df
    manter = pd.Series(True, index=df.index)

    termos = pedido.termos_alvo()
    if termos:
        em_algum = pd.Series(False, index=df.index)
        for coluna in ("alvo", "uniprot", "nome", "titulo", "classe_alvo"):
            if coluna in df.columns:
                em_algum |= _contem(df[coluna], termos)
        manter &= em_algum

    if pedido.e3:
        termos_e3 = [normalizar_nome(pedido.e3)]
        em_e3 = pd.Series(False, index=df.index)
        for coluna in ("e3", "nome", "titulo"):
            if coluna in df.columns:
                em_e3 |= _contem(df[coluna], termos_e3)
        manter &= em_e3

    if pedido.linhagem and "linhagem" in df.columns:
        manter &= _contem(df["linhagem"], [normalizar_nome(pedido.linhagem)])

    if pedido.texto_livre:
        termo = [normalizar_nome(pedido.texto_livre)]
        em_qualquer = pd.Series(False, index=df.index)
        for coluna in df.columns:
            em_qualquer |= _contem(df[coluna], termo)
        manter &= em_qualquer

    if pedido.so_com_degradacao:
        tem = pd.Series(False, index=df.index)
        for coluna in ("dc50", "dmax"):
            if coluna in df.columns:
                tem |= df[coluna].notna()
        manter &= tem

    return df[manter].copy()


# ===========================================================================
# 8. DEDUPLICAÇÃO E CONSOLIDAÇÃO ENTRE BASES
# ===========================================================================

def canonizar_smiles(smiles: Any) -> Any:
    """SMILES canônico via RDKit se houver; senão, o texto limpo.

    O RDKit é opcional de propósito — sem ele a deduplicação fica menos
    precisa (duas grafias do mesmo composto viram duas linhas), e o relatório
    avisa. Nenhuma outra parte do módulo precisa dele.
    """
    texto = limpar_texto(smiles)
    if texto is pd.NA:
        return pd.NA
    try:
        from rdkit import Chem                      # noqa: PLC0415
        from rdkit import RDLogger                  # noqa: PLC0415
        RDLogger.DisableLog("rdApp.*")
        mol = Chem.MolFromSmiles(str(texto))
        if mol is not None:
            return Chem.MolToSmiles(mol)
    except Exception:
        pass
    return str(texto).strip()


def rdkit_disponivel() -> bool:
    try:
        import rdkit  # noqa: F401,PLC0415
        return True
    except Exception:
        return False


def derivar_comprimento_linker(df: pd.DataFrame) -> pd.DataFrame:
    """Preenche 'Comprimeto do Linker' contando átomos pesados do linker.

    Só age onde a coluna está vazia E existe um SMILES de linker, e só com
    RDKit. O valor sai marcado com '(calc.)' para não se confundir com o que a
    base afirmou: a contagem de átomos pesados não é a mesma definição de
    "comprimento" que todo artigo usa.
    """
    if "linker_smiles" not in df.columns or not rdkit_disponivel():
        return df
    from rdkit import Chem                          # noqa: PLC0415
    from rdkit import RDLogger                      # noqa: PLC0415
    RDLogger.DisableLog("rdApp.*")
    out = df.copy()
    novos = []
    for atual, linker in zip(out.get("linker_comprimento", pd.Series(dtype="object")),
                             out["linker_smiles"]):
        if not (atual is pd.NA or pd.isna(atual)):
            novos.append(atual)
            continue
        texto = limpar_texto(linker)
        if texto is pd.NA:
            novos.append(pd.NA)
            continue
        mol = Chem.MolFromSmiles(str(texto))
        novos.append(f"{mol.GetNumHeavyAtoms()} (calc.)" if mol is not None else pd.NA)
    out["linker_comprimento"] = novos
    return out


def consolidar(df: pd.DataFrame) -> pd.DataFrame:
    """Junta linhas do mesmo composto vindas de bases diferentes.

    A chave é o SMILES canônico. Para cada parâmetro fico com o primeiro valor
    não-vazio, na ordem em que as bases foram pedidas — então a ordem de
    `Pedido.bases` é uma ordem de PREFERÊNCIA, e vale pensar nela: PROTAC-DB
    antes de PROTAC-PatentDB faz a medida experimental ganhar da prevista.

    Linhas sem SMILES não são consolidadas (não há como saber se são o mesmo
    composto), e passam intactas.
    """
    if df.empty or "smiles" not in df.columns:
        return df
    out = df.copy()
    out["_chave"] = [canonizar_smiles(s) for s in out["smiles"]]

    sem_chave = out[out["_chave"].isna()].drop(columns=["_chave"])
    com_chave = out[out["_chave"].notna()]
    if com_chave.empty:
        return sem_chave

    linhas = []
    for _, grupo in com_chave.groupby("_chave", sort=False):
        linha: dict[str, Any] = {}
        for coluna in grupo.columns:
            if coluna == "_chave":
                continue
            if coluna in ("fonte", "id_fonte", "arquivo_origem"):
                valores = [str(v) for v in grupo[coluna] if not pd.isna(v)]
                # sem perder a ordem, remove repetido
                vistos: set[str] = set()
                unicos = [v for v in valores if not (v in vistos or vistos.add(v))]
                linha[coluna] = "; ".join(unicos) if unicos else pd.NA
            else:
                nao_vazios = grupo[coluna].dropna()
                linha[coluna] = nao_vazios.iloc[0] if len(nao_vazios) else pd.NA
        linhas.append(linha)

    consolidado = pd.DataFrame(linhas)
    return pd.concat([consolidado, sem_chave], ignore_index=True)


# ===========================================================================
# 9. A COLETA
# ===========================================================================

@dataclass
class Resultado:
    """O que a coleta devolve: a tabela e tudo que justifica a tabela."""
    tabela: pd.DataFrame
    auditorias: dict[str, dict] = field(default_factory=dict)
    registro: list[str] = field(default_factory=list)
    contagens: dict[str, int] = field(default_factory=dict)

    def relatorio(self) -> str:
        """Texto para imprimir no notebook: o que entrou, de onde, e o que falta."""
        linhas = ["=" * 72, "COLETA — RESUMO", "=" * 72]
        for item in self.registro:
            linhas.append(f"  {item}")
        linhas.append("")
        linhas.append(f"  linhas na tabela final: {len(self.tabela)}")
        for base, n in self.contagens.items():
            linhas.append(f"    {base}: {n}")
        if not self.tabela.empty:
            linhas.append("")
            linhas.append("  preenchimento por parâmetro (quantas linhas têm o dado):")
            for chave in CHAVES:
                if chave == "idx" or chave not in self.tabela.columns:
                    continue
                n = int(self.tabela[chave].notna().sum())
                pct = 100.0 * n / len(self.tabela)
                marca = "  " if n else "··"
                linhas.append(f"    {marca} {ROTULOS[chave]:<24} {n:>6}  ({pct:5.1f}%)")
        faltando = {}
        for base, aud in self.auditorias.items():
            if aud.get("sem_destino"):
                faltando[base] = aud["sem_destino"]
        if faltando:
            linhas.append("")
            linhas.append("  COLUNAS DE ORIGEM SEM DESTINO — revise, pode ter")
            linhas.append("  parâmetro que a base traz e o mapa não reconheceu:")
            for base, colunas in faltando.items():
                mostra = ", ".join(colunas[:12])
                resto = f" (+{len(colunas) - 12})" if len(colunas) > 12 else ""
                linhas.append(f"    {base}: {mostra}{resto}")
            linhas.append("  Para corrigir: acrescente o padrão em APELIDOS e rode de novo.")
        return "\n".join(linhas)


# O PROTAC-DB distribui, ao lado da tabela dos PROTACs, três tabelas de
# FRAGMENTO: warhead, linker e e3_ligand. Elas têm SMILES, então "tem SMILES?"
# não serve para distinguir — e um fragmento que entrasse na planilha pareceria
# um PROTAC de 200 Da sem atividade. Separo pelo nome do arquivo e pelo nome da
# coluna de identificação, que é onde a diferença está escrita.
_MARCAS_DE_FRAGMENTO = (r"\bwarhead\b", r"\be3 ligand\b", r"\blinker id\b",
                        r"\bwarhead id\b", r"\be3 ligand id\b")
_ARQUIVOS_DE_FRAGMENTO = (r"warhead", r"^linker", r"e3[_ -]?ligand")

# Para valer como tabela de composto, preciso de atividade OU de um SMILES
# acompanhado de contexto de PROTAC (massa, bibliografia, E3, linhagem...).
_CHAVES_DE_ATIVIDADE = {"dc50", "dmax", "ic50", "ec50", "gi50", "ed50", "gr50"}
_CHAVES_DE_CONTEXTO = {"mw", "doi", "e3", "linhagem", "linker_tipo", "titulo",
                       "data", "uniprot", "patente", "atividade", "classe_alvo"}


def motivo_para_ignorar(caminho: Path, colunas: Iterable[Any],
                        mapa: dict[str, str]) -> str | None:
    """Devolve o motivo para ignorar o arquivo, ou None para processá-lo."""
    nome = normalizar_nome(caminho.stem)
    if any(re.search(p, nome) for p in _ARQUIVOS_DE_FRAGMENTO):
        return "o nome do arquivo é de tabela de fragmento"
    normalizadas = [normalizar_nome(c) for c in colunas]
    marcadas = [n for n in normalizadas
                if any(re.search(p, n) for p in _MARCAS_DE_FRAGMENTO)]
    if marcadas and not (_CHAVES_DE_ATIVIDADE & set(mapa)):
        return f"colunas de fragmento ({', '.join(marcadas[:3])}) e nenhuma atividade"
    if _CHAVES_DE_ATIVIDADE & set(mapa):
        return None
    if "smiles" in mapa and len(_CHAVES_DE_CONTEXTO & set(mapa)) >= 2:
        return None
    if "smiles" not in mapa:
        return "nenhum SMILES reconhecido"
    return "só SMILES, sem atividade nem contexto suficiente"


def carregar_base(nome_base: str, dir_dados: str | Path, baixar: bool = True,
                  forcar: bool = False) -> tuple[pd.DataFrame, dict, list[str]]:
    """Lê uma base inteira e devolve (tabela no vocabulário da planilha, auditoria, registro)."""
    fonte = FONTES[nome_base]
    arquivos, registro = obter_arquivos(fonte, dir_dados, baixar=baixar, forcar=forcar)
    if not arquivos:
        return pd.DataFrame(columns=TODAS_AS_CHAVES), {"erro": "sem arquivos"}, registro

    quadros: list[pd.DataFrame] = []
    auditoria_geral: dict[str, Any] = {"arquivos": {}}
    for caminho in arquivos:
        try:
            bruto = ler_tabela(caminho)
        except Exception as erro:
            registro.append(f"{nome_base}: não li {caminho.name} "
                            f"({type(erro).__name__}: {erro})")
            continue
        if bruto is None or bruto.empty:
            continue
        mapa, auditoria = mapear_colunas(bruto.columns)
        motivo = motivo_para_ignorar(caminho, bruto.columns, mapa)
        if motivo:
            registro.append(f"{nome_base}: {caminho.name} não parece tabela de "
                            f"composto ({motivo}) — ignorado")
            continue
        posto = aplicar_mapa(bruto, mapa)
        posto["fonte"] = nome_base
        posto["arquivo_origem"] = caminho.name
        quadros.append(posto)
        auditoria_geral["arquivos"][caminho.name] = auditoria
        registro.append(f"{nome_base}: {caminho.name} -> {len(posto)} linhas, "
                        f"{len(mapa)} coluna(s) reconhecida(s)")

    if not quadros:
        return pd.DataFrame(columns=TODAS_AS_CHAVES), auditoria_geral, registro

    junto = pd.concat(quadros, ignore_index=True)
    # A auditoria "achatada" junta o que cada arquivo relatou, para o relatório.
    sem_destino: list[str] = []
    sem_origem: set[str] = set(CHAVES) - {"idx"}
    for aud in auditoria_geral["arquivos"].values():
        sem_destino.extend(aud.get("sem_destino", []))
        sem_origem &= set(aud.get("sem_origem", []))
    vistos: set[str] = set()
    auditoria_geral["sem_destino"] = [c for c in sem_destino
                                      if not (c in vistos or vistos.add(c))]
    auditoria_geral["sem_origem"] = sorted(sem_origem)
    return junto, auditoria_geral, registro


def coletar(pedido: Pedido, dir_dados: str | Path = "dados/bases",
            baixar: bool = True, forcar: bool = False) -> Resultado:
    """Percorre as bases pedidas, filtra, padroniza e consolida.

    Uma base que falhe (sem rede, exige cadastro, formato novo) NÃO interrompe
    as outras: o `registro` do resultado diz o que aconteceu com cada uma.
    """
    quadros: list[pd.DataFrame] = []
    auditorias: dict[str, dict] = {}
    registro: list[str] = []
    contagens: dict[str, int] = {}

    for nome_base in pedido.bases:
        if nome_base not in FONTES:
            registro.append(f"{nome_base}: base desconhecida, ignorada "
                            f"(conhecidas: {', '.join(FONTES)})")
            continue
        bruto, auditoria, reg = carregar_base(nome_base, dir_dados,
                                              baixar=baixar, forcar=forcar)
        registro.extend(reg)
        auditorias[nome_base] = auditoria
        if bruto.empty:
            contagens[nome_base] = 0
            continue
        padronizado = padronizar_valores(bruto, unidades_fixas=pedido.unidades_fixas)
        selecionado = filtrar(padronizado, pedido)
        contagens[nome_base] = len(selecionado)
        registro.append(f"{nome_base}: {len(selecionado)} de {len(padronizado)} "
                        f"linha(s) passaram no filtro")
        if selecionado.empty and not padronizado.empty and pedido.termos_alvo():
            tem_nome_de_alvo = any(
                coluna in padronizado.columns and padronizado[coluna].notna().any()
                for coluna in ("alvo", "nome", "titulo")
            )
            if not tem_nome_de_alvo:
                registro.append(
                    f"    ATENÇÃO: {nome_base} não identifica o alvo por nome "
                    f"(só por UniProt). Filtrar por '{pedido.alvo}' zera esta "
                    f"base — preencha também uniprot= no Pedido."
                )
        if not selecionado.empty:
            quadros.append(selecionado)

    if not quadros:
        return Resultado(pd.DataFrame(columns=TODAS_AS_CHAVES), auditorias,
                         registro, contagens)

    tabela = pd.concat(quadros, ignore_index=True)
    tabela = derivar_comprimento_linker(tabela)
    if pedido.consolidar_duplicatas:
        antes = len(tabela)
        tabela = consolidar(tabela)
        registro.append(f"consolidação por SMILES: {antes} -> {len(tabela)} linha(s)"
                        + ("" if rdkit_disponivel() else
                           "  [sem RDKit: comparação por texto, menos precisa]"))

    # Ordem útil para revisar: quem tem DC₅₀ primeiro, que é o dado mais caro.
    if "dc50" in tabela.columns:
        tabela = tabela.sort_values(
            by=["dc50", "dmax"], key=lambda s: s.isna(), kind="stable"
        ).reset_index(drop=True)

    tabela["idx"] = range(1, len(tabela) + 1)
    colunas = [c for c in TODAS_AS_CHAVES if c in tabela.columns]
    tabela = tabela[colunas]
    return Resultado(tabela, auditorias, registro, contagens)


# ===========================================================================
# 10. EXPORTAÇÃO PARA A PLANILHA
# ===========================================================================

def _escrever_aba_papers(ws, tabela: pd.DataFrame, linha_inicial: int = 4) -> None:
    """Escreve os dados nas 28 colunas da aba PAPERS, na ordem da planilha."""
    for deslocamento, (_, registro) in enumerate(tabela.iterrows()):
        linha = linha_inicial + deslocamento
        for coluna_idx, chave in enumerate(CHAVES, start=1):
            valor = registro.get(chave, pd.NA)
            if valor is pd.NA or (not isinstance(valor, str) and pd.isna(valor)):
                valor = None
            elif not isinstance(valor, (int, float, str)):
                valor = str(valor)
            ws.cell(row=linha, column=coluna_idx, value=valor)


def exportar_planilha(tabela: pd.DataFrame, destino: str | Path,
                      modelo: str | Path | None = None,
                      aba: str = "PAPERS", linha_inicial: int = 4,
                      resultado: "Resultado | None" = None) -> Path:
    """Grava o arquivo final com as colunas da aba PAPERS.

    Com `modelo` apontando para a sua planilha, a formatação dela é
    preservada e os dados entram a partir de `linha_inicial` (4, logo abaixo do
    cabeçalho da linha 3). Sem modelo, o arquivo é criado do zero com o mesmo
    desenho. **O modelo nunca é sobrescrito** — a saída é sempre `destino`.

    As abas PROVENIENCIA (de onde veio cada linha) e MAPEAMENTO (qual coluna da
    base alimentou qual parâmetro) vão junto. São elas que tornam a revisão
    defensável: sem elas, a planilha é um monte de número sem procedência.
    """
    from openpyxl import Workbook, load_workbook              # noqa: PLC0415
    from openpyxl.styles import Font                           # noqa: PLC0415

    destino = Path(destino)
    destino.parent.mkdir(parents=True, exist_ok=True)

    if modelo is not None:
        modelo = Path(modelo)
        if not modelo.exists():
            raise FileNotFoundError(f"modelo não encontrado: {modelo}")
        if destino.resolve() == modelo.resolve():
            raise ValueError("destino igual ao modelo: a planilha original não "
                             "deve ser sobrescrita. Escolha outro nome.")
        wb = load_workbook(modelo)
        if aba not in wb.sheetnames:
            raise KeyError(f"a aba {aba!r} não existe no modelo "
                           f"(abas: {', '.join(wb.sheetnames)})")
        ws = wb[aba]
    else:
        wb = Workbook()
        ws = wb.active
        ws.title = aba
        for coluna_idx, chave in enumerate(CHAVES, start=1):
            celula = ws.cell(row=3, column=coluna_idx, value=ROTULOS[chave])
            celula.font = Font(bold=True)

    _escrever_aba_papers(ws, tabela, linha_inicial=linha_inicial)

    # --- PROVENIENCIA: de onde veio cada linha ----------------------------
    if "PROVENIENCIA" in wb.sheetnames:
        del wb["PROVENIENCIA"]
    wsp = wb.create_sheet("PROVENIENCIA")
    colunas_prov = ["idx"] + [c for c in CHAVES_PROVENIENCIA if c in tabela.columns]
    for coluna_idx, chave in enumerate(colunas_prov, start=1):
        wsp.cell(row=1, column=coluna_idx, value=chave).font = Font(bold=True)
    for deslocamento, (_, registro) in enumerate(tabela.iterrows()):
        for coluna_idx, chave in enumerate(colunas_prov, start=1):
            valor = registro.get(chave, pd.NA)
            if valor is pd.NA or (not isinstance(valor, str) and pd.isna(valor)):
                valor = None
            wsp.cell(row=2 + deslocamento, column=coluna_idx,
                     value=valor if isinstance(valor, (int, float, str, type(None)))
                     else str(valor))

    # --- MAPEAMENTO: qual coluna da base alimentou qual parâmetro ---------
    if resultado is not None:
        if "MAPEAMENTO" in wb.sheetnames:
            del wb["MAPEAMENTO"]
        wsm = wb.create_sheet("MAPEAMENTO")
        for coluna_idx, titulo in enumerate(
                ["base", "parâmetro (planilha)", "coluna de origem", "situação"], start=1):
            wsm.cell(row=1, column=coluna_idx, value=titulo).font = Font(bold=True)
        linha = 2
        for base, auditoria in resultado.auditorias.items():
            for arquivo_aud in auditoria.get("arquivos", {}).values():
                for chave, origem in arquivo_aud.get("mapa", {}).items():
                    wsm.cell(row=linha, column=1, value=base)
                    wsm.cell(row=linha, column=2, value=ROTULOS.get(chave, chave))
                    wsm.cell(row=linha, column=3, value=str(origem))
                    wsm.cell(row=linha, column=4, value="mapeado")
                    linha += 1
            for coluna in auditoria.get("sem_destino", []):
                wsm.cell(row=linha, column=1, value=base)
                wsm.cell(row=linha, column=2, value="(nenhum)")
                wsm.cell(row=linha, column=3, value=str(coluna))
                wsm.cell(row=linha, column=4, value="coluna da base sem destino")
                linha += 1
            for chave in auditoria.get("sem_origem", []):
                wsm.cell(row=linha, column=1, value=base)
                wsm.cell(row=linha, column=2, value=ROTULOS.get(chave, chave))
                wsm.cell(row=linha, column=3, value="(nenhuma)")
                wsm.cell(row=linha, column=4, value="base não traz este parâmetro")
                linha += 1

    wb.save(destino)
    return destino


def exportar_csv(tabela: pd.DataFrame, destino: str | Path) -> Path:
    """CSV com os rótulos da planilha, para quem preferir conferir fora do Excel."""
    destino = Path(destino)
    destino.parent.mkdir(parents=True, exist_ok=True)
    saida = tabela[[c for c in CHAVES if c in tabela.columns]].rename(columns=ROTULOS)
    saida.to_csv(destino, index=False)
    return destino


# ===========================================================================
# 11. LINHA DE COMANDO
# ===========================================================================

def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Coleta os parâmetros da aba PAPERS nas bases de PROTAC.")
    ap.add_argument("--alvo", help="proteína-alvo, p.ex. BRD4")
    ap.add_argument("--sinonimo", action="append", default=[],
                    help="outra grafia do alvo (pode repetir)")
    ap.add_argument("--e3", help="E3 ligase, p.ex. VHL ou CRBN")
    ap.add_argument("--linhagem", help="linhagem celular")
    ap.add_argument("--texto", help="busca livre em qualquer coluna")
    ap.add_argument("--bases", default=",".join(FONTES),
                    help="bases a usar, separadas por vírgula")
    ap.add_argument("--dir-dados", default="dados/bases")
    ap.add_argument("--sem-download", action="store_true",
                    help="usa só o que já está em --dir-dados")
    ap.add_argument("--forcar-download", action="store_true")
    ap.add_argument("--so-com-degradacao", action="store_true")
    ap.add_argument("--unidades-fixas", action="store_true",
                    help="potências em nM, número puro (bom para ordenar)")
    ap.add_argument("--sem-consolidar", action="store_true")
    ap.add_argument("--modelo", help="planilha a usar como molde (não é sobrescrita)")
    ap.add_argument("--saida", default="revisao_protac.xlsx")
    ap.add_argument("--saida-csv")
    args = ap.parse_args(argv)

    pedido = Pedido(
        alvo=args.alvo,
        sinonimos=args.sinonimo,
        e3=args.e3,
        linhagem=args.linhagem,
        texto_livre=args.texto,
        bases=[b.strip() for b in args.bases.split(",") if b.strip()],
        so_com_degradacao=args.so_com_degradacao,
        unidades_fixas=args.unidades_fixas,
        consolidar_duplicatas=not args.sem_consolidar,
    )
    resultado = coletar(pedido, dir_dados=args.dir_dados,
                        baixar=not args.sem_download, forcar=args.forcar_download)
    print(resultado.relatorio())
    if resultado.tabela.empty:
        print("\nNada para exportar.")
        return 1
    caminho = exportar_planilha(resultado.tabela, args.saida, modelo=args.modelo,
                                resultado=resultado)
    print(f"\nplanilha: {caminho}")
    if args.saida_csv:
        print(f"csv:      {exportar_csv(resultado.tabela, args.saida_csv)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
