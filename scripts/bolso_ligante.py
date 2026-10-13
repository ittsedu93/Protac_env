#!/usr/bin/env python
"""
Tudo que está a menos de 5 Å do ligante: resíduos, outras cadeias, cofatores, íons.

Só biblioteca padrão. Lê PDB e mmCIF. Instantâneo.

    python bolso_ligante.py --estrutura modelo.cif --ligante-auto
    python bolso_ligante.py --estrutura modelo.pdb --ligante-resname LIG --raio 5.0
    python bolso_ligante.py --estrutura modelo.cif --ligante-cadeia C --ligante-resnum 1

Por que não basta "os resíduos do bolso"
----------------------------------------
A pergunta usual ("quais resíduos formam o bolso") esconde três coisas que
decidem o resultado num ternário de PROTAC:

  1. QUANTAS CADEIAS o ligante toca. Se o PROTAC está a 5 Å só da cadeia do
     alvo, o modelo não montou o ternário — montou um binário com a E3
     pendurada. O bolso de interface (>= 2 cadeias poliméricas) é a assinatura
     geométrica do que se quer.
  2. COFATORES e ÍONS. Um Zn²⁺ ou um Mg²⁺ a 3 Å do ligante muda a química do
     sítio e some de qualquer lista que filtre HETATM. Em co-cristais reais
     metade da vizinhança a 5 Å costuma ser glicerol, EDO e sulfato da
     cristalização — e isso é informação, não ruído, desde que etiquetada.
  3. ÁGUAS. Nos modelos de cofold não existem; num co-cristal existem e
     ocupam o bolso. Entram na lista, etiquetadas, para poder somar ou tirar.

Então aqui a saída é uma lista COMPLETA de tudo dentro do raio, cada item com
categoria (proteina / acido_nucleico / ion / cofator / aditivo / glicano /
lipideo / agua / hetero_desconhecido / mesma_cadeia_do_ligante), a distância
mínima e o par de átomos que a realiza.

Convenções
----------
- Distância entre ÁTOMOS PESADOS (hidrogênio entra só com --com-hidrogenios):
  é o que a literatura de contato usa, e os modelos de cofold não têm H.
- altLoc: fica a principal (' ' ou 'A'); as alternativas seriam contagem dupla.
- Só o primeiro MODEL, salvo --modelo N.
- mmCIF: usa auth_asym_id/auth_seq_id quando existem (é a numeração que o
  usuário vê no PyMOL); cai para label_* quando não.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

# ---------------------------------------------------------------------------
# dicionários de resíduo -> o que é

AA3_PARA_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
    # protonação/variantes que pdb2gmx, Amber e os modelos escrevem
    "HID": "H", "HIE": "H", "HIP": "H", "HISD": "H", "HISE": "H", "HISP": "H",
    "CYX": "C", "CYM": "C", "ASH": "D", "GLH": "E", "LYN": "K", "ARN": "R",
    "TYM": "Y", "CYT": "C",
    # modificados frequentes (contam como proteína: estão na cadeia)
    "MSE": "M", "SEC": "C", "PYL": "K", "SEP": "S", "TPO": "T", "PTR": "Y",
    "CSO": "C", "CME": "C", "CSD": "C", "OCS": "C", "CAS": "C", "CSS": "C",
    "KCX": "K", "MLY": "K", "M3L": "K", "ALY": "K", "MLZ": "K",
    "HYP": "P", "PCA": "Q", "FME": "M", "SAC": "S", "AIB": "A",
    "ABA": "A", "ORN": "K", "DAL": "A", "DLE": "L", "DVA": "V",
    "NLE": "L", "PHI": "F", "TYS": "Y", "TRQ": "W", "LLP": "K",
}
TAPAS = {"ACE", "NME", "NH2", "FOR", "NHE"}          # caps N/C-terminais

NUCLEICOS = {
    "DA", "DC", "DG", "DT", "DU", "DI",
    "A", "C", "G", "U", "I", "N",
    "5MC", "7MG", "1MA", "2MG", "M2G", "H2U", "PSU", "5MU", "OMC", "OMG",
    "4SU", "G7M", "YG", "QUO", "UR3", "RIA",
}

AGUAS = {"HOH", "WAT", "DOD", "H2O", "SOL", "TIP", "TIP3", "TIP4", "SPC", "OH2"}

IONS_MONOATOMICOS = {
    "NA", "K", "LI", "RB", "CS", "MG", "CA", "SR", "BA", "ZN", "MN", "FE",
    "FE2", "FE3", "CO", "3CO", "NI", "CU", "CU1", "CU2", "CD", "HG", "AU",
    "AU3", "AG", "PT", "PT4", "PD", "PB", "TL", "AL", "GA", "IN", "SN",
    "V", "CR", "MO", "W", "YB", "SM", "EU", "GD", "LU", "CE", "TB", "HO",
    "ER", "DY", "NA1", "CL", "BR", "IOD", "F", "SE", "SE4", "YT3", "4MO",
    "6MO", "ZN2", "NH4", "OH", "O",
}
IONS_POLIATOMICOS = {
    "SO4", "SUL", "PO4", "2HP", "PI", "NO3", "NO2", "CO3", "SCN", "AZI",
    "N3", "CN", "BF4", "PF6", "CLO", "CLO3", "CLO4", "BR3", "I3M", "TCN",
    "WO4", "MOO", "MO4", "VO4", "ASO", "AST", "SEO", "SO3", "PO3", "BO3",
    "BO4", "CAC", "ACT", "ACY", "FLC", "OXL", "FMT", "PER", "PEO", "IUM",
}
IONS = IONS_MONOATOMICOS | IONS_POLIATOMICOS

COFATORES = {
    # hemes e porfirinas
    "HEM", "HEA", "HEB", "HEC", "HED", "HEV", "HAS", "DHE", "HNI", "VER",
    "COH", "CLA", "CHL", "BCL", "BCB", "PHO", "SRM", "MQ7", "MQ8", "MQ9",
    # nicotinamidas / flavinas
    "NAD", "NAI", "NAJ", "NAP", "NDP", "NAX", "NAQ", "DND", "NDC", "ZID",
    "FAD", "FDA", "FAS", "FMN", "FNS", "RBF", "FFO",
    # nucleotídeos e análogos (sinalização, ATPases, quinases)
    "ATP", "ADP", "AMP", "ANP", "ACP", "AGS", "APC", "ADX", "A2P", "AP5",
    "GTP", "GDP", "GMP", "GNP", "GSP", "GCP", "G2P", "5GP",
    "CTP", "CDP", "CMP", "UTP", "UDP", "UMP", "TTP", "TDP", "TMP", "DTP",
    "PPV", "POP", "PPS", "PPF", "P3G",
    # coenzimas de transferência
    "COA", "ACO", "COO", "CAA", "SCA", "MCA", "HAX", "3HC",
    "SAM", "SAH", "SFG", "MTA", "DCS",
    "PLP", "PMP", "PXL", "P1T", "TPP", "THD", "TDL", "TPO1",
    "BTN", "BTI", "B12", "COB", "CNC", "B1Z", "BYA",
    "THF", "THG", "FOL", "5FU0", "MTF", "5MF",
    "GSH", "GDS", "GTT", "GSF", "GSB",
    "PQQ", "LPA", "LPB", "TRQ1", "RET", "A1R", "SIN1",
    # clusters metálicos e metalocentros
    "SF4", "FES", "FS4", "F3S", "CLF", "CFM", "ICS", "HCA", "HC0", "HC1",
    "MGD", "MSS", "MOS", "MOO1", "WCC", "NFU", "SRI", "CUA", "CUZ", "OEC",
    "CL1", "NI1", "NIK", "MNH", "CO5", "COY", "ZNO",
}

ADITIVOS = {                      # crioprotetor, tampão, precipitante
    "GOL", "EDO", "PEG", "PG4", "PGE", "1PE", "2PE", "7PE", "12P", "P6G",
    "PE3", "PE4", "PE5", "PE8", "XPE", "P33", "P4C", "MPD", "MRD", "DIO",
    "DOX", "TRS", "MES", "EPE", "BTB", "TAM", "IMD", "IMZ", "HEZ", "BU1",
    "BU3", "BUD", "PDO", "PGO", "PGR", "ETE", "DMS", "DMF", "DMA", "NMO",
    "ACN", "ACE1", "MOH", "EOH", "IPA", "IPH", "TFA", "TFE", "URE", "GAI",
    "FLU1", "BME", "DTT", "DTU", "DTV", "TCE", "MBO", "MYR1", "SIN",
    "CIT", "FLC1", "TAR", "TLA", "MLA", "MLI", "MAE", "SUC", "GLV",
    "AKG", "PYR", "LAC", "ACM", "ACA", "GLY1", "BEZ", "PHQ", "SPD", "SPM",
    "PUT", "CAD", "EDT", "EDA", "TBU", "NHE1", "C2E1", "AZE",
}

GLICANOS = {
    "NAG", "NDG", "NGA", "A2G", "BMA", "MAN", "BGC", "GLC", "GAL", "GLA",
    "GALA", "FUC", "FUL", "XYL", "XYS", "XYP", "SIA", "NAN", "NGC", "RAM",
    "RIB", "ARA", "ARB", "LXC", "GCU", "ADA", "IDS", "SGN", "GCS", "BDP",
    "MMA", "MBG", "TRE", "SUC1", "LAT", "LBT", "MAL1", "CBI", "GLP",
}

LIPIDEOS = {
    "LDA", "LMT", "DDQ", "BOG", "BNG", "OCT", "C8E", "CE9", "HTG", "HTO",
    "MYR", "PLM", "STE", "OLA", "OLB", "OLC", "ELA", "LAP", "DAO", "D10",
    "D12", "DEP", "HEX", "HP6", "UND", "DKA", "PEE", "PEF", "PEH", "PEV",
    "PC1", "PCW", "PCF", "POV", "LHG", "CDL", "SQD", "DGA", "MGE", "CHD",
    "CHS", "CLR", "CHO", "Y01", "SOG", "F09", "LPP", "PX4", "TRD", "TWT",
}

METAIS = {
    "NA", "K", "LI", "RB", "CS", "MG", "CA", "SR", "BA", "ZN", "MN", "FE",
    "CO", "NI", "CU", "CD", "HG", "AU", "AG", "PT", "PD", "PB", "TL", "AL",
    "GA", "IN", "SN", "V", "CR", "MO", "W", "YB", "SM", "EU", "GD", "LU",
    "CE", "TB", "HO", "ER", "DY", "BE", "TI", "ZR", "NB", "RU", "RH", "OS",
    "IR", "RE", "TA", "HF", "SC", "Y", "BI", "AS", "SB",
}


def categorizar(resname: str, n_atomos: int = 0, elementos=()) -> tuple[str, str]:
    """(categoria, subtipo). A ordem das perguntas é a que importa.

    Água antes de tudo (HOH é HETATM e seria 'hetero'); proteína/nucleico pela
    tabela de resíduos; íon monoatômico também se for um resíduo de 1 átomo
    cujo elemento é metal/halogênio, porque modelos e campos de força inventam
    nomes de resíduo para íons.
    """
    r = (resname or "").strip().upper()
    els = {str(e).strip().upper() for e in elementos if e}
    if r in AGUAS:
        return "agua", "agua"
    if r in AA3_PARA_1:
        return "proteina", "residuo_padrao" if r in (
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
            "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
            "TYR", "VAL",
        ) else "residuo_modificado"
    if r in TAPAS:
        return "proteina", "cap_terminal"
    if r in NUCLEICOS:
        return "acido_nucleico", "dna" if r.startswith("D") else "rna"
    if r in IONS_MONOATOMICOS:
        return "ion", "monoatomico"
    if r in IONS_POLIATOMICOS:
        return "ion", "poliatomico"
    if r in COFATORES:
        return "cofator", "cofator"
    if r in GLICANOS:
        return "glicano", "glicano"
    if r in LIPIDEOS:
        return "lipideo", "lipideo_ou_detergente"
    if r in ADITIVOS:
        return "aditivo", "cristalizacao"
    # só DEPOIS das tabelas por nome: um resíduo de um átomo cujo elemento é
    # metal/halogênio é íon mesmo com nome exótico (campos de força e modelos
    # inventam nomes). Antes das tabelas, isto chamaria de íon um HEM do qual
    # só sobrou o Fe.
    if n_atomos == 1 and els and els <= METAIS | {"CL", "BR", "I", "F", "SE"}:
        return "ion", "monoatomico_por_elemento"
    return "hetero_desconhecido", "hetero"


POLIMEROS = {"proteina", "acido_nucleico"}

# ---------------------------------------------------------------------------
# leitura de estrutura


@dataclass(frozen=True)
class Atomo:
    serial: int
    nome: str
    altloc: str
    resname: str
    cadeia: str
    resnum: int
    icode: str
    x: float
    y: float
    z: float
    elemento: str
    hetero: bool
    modelo: int = 1

    @property
    def chave_residuo(self):
        return (self.cadeia, self.resnum, self.icode, self.resname)

    @property
    def rotulo_residuo(self):
        return f"{self.cadeia}/{self.resname}{self.resnum}{self.icode.strip()}"


def _elemento_do_nome(nome: str, resname: str = "") -> str:
    """Quando a coluna de elemento falta (PDB velho, saída de ferramenta)."""
    n = nome.strip()
    if not n:
        return ""
    r = (resname or "").strip().upper()
    if r in IONS_MONOATOMICOS or r in METAIS:
        # 'ZN', 'NA', 'FE2' -> o próprio nome do resíduo diz o elemento
        base = "".join(c for c in r if c.isalpha())
        if base in METAIS or base in ("CL", "BR", "I", "F", "SE"):
            return base
    if n[0].isdigit():
        n = n[1:]
    if len(n) >= 2 and n[:2].upper() in ("CL", "BR", "ZN", "MG", "NA", "FE",
                                         "MN", "CU", "CA", "NI", "CO", "SE"):
        if r not in AA3_PARA_1:              # em proteína, CA é carbono alfa
            return n[:2].upper()
    return n[0].upper()


def ler_pdb(texto: str, modelo: int = 1) -> list[Atomo]:
    atomos: list[Atomo] = []
    mod_atual = 1
    for linha in texto.splitlines():
        reg = linha[:6]
        if reg == "MODEL ":
            try:
                mod_atual = int(linha[10:14])
            except ValueError:
                mod_atual += 1
            continue
        if reg == "ENDMDL":
            mod_atual += 1
            continue
        if reg not in ("ATOM  ", "HETATM"):
            continue
        if mod_atual != modelo:
            continue
        if len(linha) < 54:
            continue
        altloc = linha[16]
        if altloc not in (" ", "A", "1"):
            continue
        try:
            x, y, z = float(linha[30:38]), float(linha[38:46]), float(linha[46:54])
        except ValueError:
            continue
        try:
            resnum = int(linha[22:26])
        except ValueError:
            resnum = 0
        nome = linha[12:16].strip()
        resname = linha[17:20].strip().upper()
        el = linha[76:78].strip().upper() if len(linha) >= 78 else ""
        try:
            serial = int(linha[6:11])
        except ValueError:
            serial = len(atomos) + 1
        atomos.append(Atomo(
            serial=serial, nome=nome, altloc=altloc, resname=resname,
            cadeia=linha[21].strip() or "_", resnum=resnum, icode=linha[26],
            x=x, y=y, z=z,
            elemento=el or _elemento_do_nome(nome, resname),
            hetero=(reg == "HETATM"), modelo=mod_atual,
        ))
    return atomos


def _tokens_cif(linha: str) -> list[str]:
    """Tokeniza uma linha de mmCIF respeitando aspas simples e duplas."""
    toks, atual, quote = [], [], None
    for ch in linha:
        if quote:
            if ch == quote:
                quote = None
                toks.append("".join(atual))
                atual = []
            else:
                atual.append(ch)
        elif ch in ("'", '"'):
            quote = ch
        elif ch.isspace():
            if atual:
                toks.append("".join(atual))
                atual = []
        else:
            atual.append(ch)
    if atual:
        toks.append("".join(atual))
    return toks


def ler_cif(texto: str, modelo: int = 1) -> list[Atomo]:
    """Parser mínimo do loop `_atom_site` do mmCIF/PDBx.

    Mínimo de propósito: Boltz-2, Chai-1/2 e OpenFold escrevem um mmCIF
    regular, e depender de gemmi/biopython só para ler coordenadas colocaria
    uma dependência pesada no caminho de uma análise que é aritmética.
    """
    linhas = texto.splitlines()
    atomos: list[Atomo] = []
    i = 0
    while i < len(linhas):
        if linhas[i].strip() != "loop_":
            i += 1
            continue
        j = i + 1
        cabecalhos = []
        while j < len(linhas) and linhas[j].strip().startswith("_"):
            cabecalhos.append(linhas[j].strip().split()[0])
            j += 1
        if not cabecalhos or not cabecalhos[0].startswith("_atom_site."):
            i = j
            continue
        idx = {c.split(".", 1)[1]: k for k, c in enumerate(cabecalhos)}

        def pega(toks, *nomes, padrao=""):
            for n in nomes:
                k = idx.get(n)
                if k is not None and k < len(toks):
                    v = toks[k]
                    if v not in (".", "?"):
                        return v
            return padrao

        while j < len(linhas):
            s = linhas[j].strip()
            if not s or s.startswith("#"):
                j += 1
                continue
            if s.startswith("_") or s == "loop_" or s.startswith("data_"):
                break
            toks = _tokens_cif(linhas[j])
            if len(toks) < len(cabecalhos):
                # linha continuada (raro em atom_site): junta a seguinte
                k = j + 1
                while k < len(linhas) and len(toks) < len(cabecalhos):
                    toks += _tokens_cif(linhas[k])
                    k += 1
                j = k - 1
            if len(toks) < len(cabecalhos):
                j += 1
                continue
            mod = pega(toks, "pdbx_PDB_model_num", padrao="1")
            try:
                mod_i = int(float(mod))
            except ValueError:
                mod_i = 1
            if mod_i != modelo:
                j += 1
                continue
            altloc = pega(toks, "label_alt_id", padrao=" ") or " "
            if altloc not in (" ", "A", "1", ".", ""):
                j += 1
                continue
            try:
                x = float(pega(toks, "Cartn_x", padrao="nan"))
                y = float(pega(toks, "Cartn_y", padrao="nan"))
                z = float(pega(toks, "Cartn_z", padrao="nan"))
            except ValueError:
                j += 1
                continue
            if math.isnan(x) or math.isnan(y) or math.isnan(z):
                j += 1
                continue
            resnum_txt = pega(toks, "auth_seq_id", "label_seq_id", padrao="0")
            try:
                resnum = int(float(resnum_txt))
            except ValueError:
                resnum = 0
            nome = pega(toks, "auth_atom_id", "label_atom_id")
            resname = pega(toks, "auth_comp_id", "label_comp_id").upper()
            cadeia = pega(toks, "auth_asym_id", "label_asym_id", padrao="_") or "_"
            grupo = pega(toks, "group_PDB", padrao="ATOM").upper()
            el = pega(toks, "type_symbol").upper()
            try:
                serial = int(float(pega(toks, "id", padrao=str(len(atomos) + 1))))
            except ValueError:
                serial = len(atomos) + 1
            icode = pega(toks, "pdbx_PDB_ins_code", padrao=" ") or " "
            atomos.append(Atomo(
                serial=serial, nome=nome,
                altloc=" " if altloc in (".", "") else altloc,
                resname=resname, cadeia=cadeia, resnum=resnum, icode=icode,
                x=x, y=y, z=z,
                elemento=el or _elemento_do_nome(nome, resname),
                hetero=(grupo == "HETATM"), modelo=mod_i,
            ))
            j += 1
        i = j
    return atomos


def ler_estrutura(caminho, modelo: int = 1) -> list[Atomo]:
    """PDB ou mmCIF, decidido pela extensão e confirmado pelo conteúdo."""
    caminho = Path(caminho)
    if not caminho.exists():
        raise ValueError(f"estrutura não encontrada: {caminho}")
    if caminho.is_dir():
        raise ValueError(f"{caminho} é um diretório, não um arquivo de estrutura")
    if caminho.stat().st_size == 0:
        raise ValueError(f"estrutura vazia (0 bytes): {caminho}")
    texto = caminho.read_text(encoding="utf-8", errors="replace")
    ext = caminho.suffix.lower()
    if ext in (".cif", ".mmcif", ".bcif"):
        atomos = ler_cif(texto, modelo)
    elif ext in (".pdb", ".ent", ".pdb1", ".brk"):
        atomos = ler_pdb(texto, modelo)
    else:
        atomos = ler_cif(texto, modelo) if "_atom_site." in texto else ler_pdb(texto, modelo)
    if not atomos:
        # o motivo mais comum é pedir um MODEL que não existe
        raise ValueError(
            f"{caminho}: nenhum átomo lido (modelo {modelo}). Confira o formato "
            "e, em NMR/ensembles, o número do modelo."
        )
    return atomos


def sem_hidrogenios(atomos: list[Atomo]) -> list[Atomo]:
    return [a for a in atomos if a.elemento.upper() not in ("H", "D")]


# ---------------------------------------------------------------------------
# resíduos e cadeias


def agrupar_residuos(atomos: list[Atomo]) -> dict:
    por_res: dict[tuple, list[Atomo]] = defaultdict(list)
    for a in atomos:
        por_res[a.chave_residuo].append(a)
    return por_res


def inventario(atomos: list[Atomo]) -> dict:
    """O que tem nesta estrutura, por cadeia e por categoria."""
    por_res = agrupar_residuos(atomos)
    cadeias: dict[str, dict] = {}
    categorias: dict[str, int] = defaultdict(int)
    for (cad, num, ic, rn), ats in por_res.items():
        cat, _ = categorizar(rn, len(ats), [a.elemento for a in ats])
        categorias[cat] += 1
        d = cadeias.setdefault(cad, {"n_residuos": 0, "n_atomos": 0,
                                     "categorias": defaultdict(int)})
        d["n_residuos"] += 1
        d["n_atomos"] += len(ats)
        d["categorias"][cat] += 1
    for d in cadeias.values():
        d["categorias"] = dict(d["categorias"])
        cats = d["categorias"]
        polim = {c: n for c, n in cats.items() if c in POLIMEROS}
        d["tipo"] = max(polim, key=polim.get) if polim else (
            max(cats, key=cats.get) if cats else "vazia"
        )
    return {"n_atomos": len(atomos), "n_residuos": len(por_res),
            "cadeias": cadeias, "categorias": dict(categorias)}


def cadeias_polimericas(atomos: list[Atomo], minimo_residuos: int = 2) -> list[str]:
    inv = inventario(atomos)
    return sorted(
        c for c, d in inv["cadeias"].items()
        if d["tipo"] in POLIMEROS and d["categorias"].get(d["tipo"], 0) >= minimo_residuos
    )


def sequencia_cadeia(atomos: list[Atomo], cadeia: str):
    """Sequência de 1 letra da cadeia + a chave de cada resíduo, na ordem.

    É a ponte para o MSA: `chaves[k]` é o resíduo da estrutura que corresponde
    ao caractere `seq[k]`. Resíduo modificado entra com o aminoácido que ele é
    (MSE -> M), porque é assim que ele aparece no MSA.
    """
    vistos = set()
    seq, chaves = [], []
    for a in atomos:
        if a.cadeia != cadeia:
            continue
        ch = (a.resnum, a.icode, a.resname)
        if ch in vistos:
            continue
        um = AA3_PARA_1.get(a.resname)
        if um is None:
            continue
        vistos.add(ch)
        seq.append(um)
        chaves.append(ch)
    ordem = sorted(range(len(chaves)), key=lambda k: (chaves[k][0], chaves[k][1]))
    return ("".join(seq[k] for k in ordem), [chaves[k] for k in ordem])


# ---------------------------------------------------------------------------
# escolha do ligante


def candidatos_a_ligante(atomos: list[Atomo], min_atomos: int = 6) -> list[dict]:
    """Hetero-resíduos que podem ser 'o ligante', do maior para o menor."""
    out = []
    for (cad, num, ic, rn), ats in agrupar_residuos(atomos).items():
        cat, sub = categorizar(rn, len(ats), [a.elemento for a in ats])
        if cat in POLIMEROS or cat == "agua":
            continue
        pesados = [a for a in ats if a.elemento.upper() not in ("H", "D")]
        out.append({
            "cadeia": cad, "resnum": num, "icode": ic, "resname": rn,
            "categoria": cat, "subtipo": sub, "n_atomos": len(pesados),
            "plausivel": cat in ("hetero_desconhecido", "cofator")
                         and len(pesados) >= min_atomos,
        })
    out.sort(key=lambda d: (d["plausivel"], d["n_atomos"]), reverse=True)
    return out


def _peptideo_curto_como_ligante(atomos: list[Atomo], max_residuos: int = 30):
    """Cadeia polimérica curta (peptídeo/macrociclo) também pode ser o ligante."""
    inv = inventario(atomos)
    curtas = [
        (c, d) for c, d in inv["cadeias"].items()
        if d["tipo"] == "proteina" and d["n_residuos"] <= max_residuos
    ]
    curtas.sort(key=lambda cd: cd[1]["n_residuos"], reverse=True)
    return curtas


def selecionar_ligante(
    atomos: list[Atomo],
    resname: str | None = None,
    cadeia: str | None = None,
    resnum: int | None = None,
    auto: bool = False,
    copia: int | None = None,
):
    """Devolve (átomos do ligante, descrição, alternativas).

    Falha com mensagem útil em vez de escolher por acaso: 'achei 3 cópias de
    LIG' é um erro que o usuário resolve em 5 segundos; escolher a cópia errada
    em silêncio é um resultado errado que ninguém nota.
    """
    por_res = agrupar_residuos(atomos)
    chaves = list(por_res)

    def filtro(ch):
        cad, num, ic, rn = ch
        if resname and rn != resname.strip().upper():
            return False
        if cadeia and cad != cadeia:
            return False
        if resnum is not None and num != resnum:
            return False
        return True

    def categoria_de(ch):
        return categorizar(ch[3], len(por_res[ch]),
                           [a.elemento for a in por_res[ch]])[0]

    if resname or cadeia or resnum is not None:
        achados = [ch for ch in chaves if filtro(ch)]
        # Pedir só a cadeia do ligante não deve trazer a proteína inteira. Se a
        # cadeia tem hetero, o ligante é o hetero; se é toda polimérica e curta,
        # o "ligante" é o peptídeo inteiro (caso dos macrociclos e stapled).
        if cadeia and not resname and resnum is None and achados:
            hetero = [ch for ch in achados
                      if categoria_de(ch) not in POLIMEROS and categoria_de(ch) != "agua"]
            if hetero:
                achados = hetero
            elif len(achados) <= 30:
                ats = [a for a in atomos if a.cadeia == cadeia]
                desc = {"cadeia": cadeia, "resnum": -1, "icode": " ",
                        "resname": f"(peptídeo cadeia {cadeia})",
                        "categoria": "proteina", "subtipo": "peptideo_curto",
                        "n_atomos": len(ats), "como": "seletor:cadeia_inteira"}
                return ats, desc, candidatos_a_ligante(atomos)
            else:
                raise ValueError(
                    f"a cadeia {cadeia!r} tem {len(achados)} resíduos poliméricos e "
                    "nenhum hetero-resíduo: ela é proteína, não ligante. Passe "
                    "--ligante-resname (e/ou --ligante-resnum)."
                )
        if not achados:
            disp = sorted({ch[3] for ch in chaves})
            raise ValueError(
                "nenhum resíduo casa com o ligante pedido "
                f"(resname={resname}, cadeia={cadeia}, resnum={resnum}). "
                f"Resíduos presentes: {', '.join(disp[:40])}"
                + (" ..." if len(disp) > 40 else "")
            )
        if len(achados) > 1 and copia is None:
            # várias cópias do mesmo ligante: o usuário escolhe
            if len({ch[3] for ch in achados}) == 1 and len(achados) <= 50:
                lista = ", ".join(f"{c}/{n}{i.strip()}" for c, n, i, _ in sorted(achados))
                raise ValueError(
                    f"{len(achados)} cópias de {achados[0][3]} na estrutura: {lista}. "
                    "Escolha com --ligante-cadeia/--ligante-resnum, ou --ligante-copia N."
                )
            raise ValueError(
                f"{len(achados)} resíduos casam com o ligante pedido; refine o seletor"
            )
        alvo = sorted(achados)[(copia - 1) if copia else 0]
        ats = por_res[alvo]
        cat, sub = categorizar(alvo[3], len(ats), [a.elemento for a in ats])
        desc = {"cadeia": alvo[0], "resnum": alvo[1], "icode": alvo[2],
                "resname": alvo[3], "categoria": cat, "subtipo": sub,
                "n_atomos": len(ats), "como": "seletor"}
        return ats, desc, candidatos_a_ligante(atomos)

    if not auto:
        raise ValueError(
            "diga qual é o ligante: --ligante-resname (e/ou --ligante-cadeia / "
            "--ligante-resnum), ou --ligante-auto para pegar o maior hetero"
        )

    # Ordem de preferência, e o motivo de cada degrau:
    #   1. hetero grande (>= 6 átomos pesados): é o ligante em 99% dos casos
    #   2. qualquer hetero que não seja íon: ligante pequeno, fragmento
    #   3. cadeia polimérica curta: peptídeo/macrociclo como ligante
    #   4. íon: se não sobrou mais nada, o "ligante" é o íon mesmo
    cands = candidatos_a_ligante(atomos)
    plaus = [c for c in cands if c["plausivel"]]
    nao_ion = [c for c in cands if c["categoria"] != "ion"]
    escolha = plaus[0] if plaus else (nao_ion[0] if nao_ion else None)
    if escolha is None:
        curtas = _peptideo_curto_como_ligante(atomos)
        if curtas:
            cad = curtas[0][0]
            ats = [a for a in atomos if a.cadeia == cad]
            desc = {"cadeia": cad, "resnum": -1, "icode": " ",
                    "resname": f"(peptídeo cadeia {cad})", "categoria": "proteina",
                    "subtipo": "peptideo_curto", "n_atomos": len(ats),
                    "como": "auto:peptideo_curto"}
            return ats, desc, cands
        if cands:
            escolha = cands[0]
    if escolha is None:
        raise ValueError(
            "--ligante-auto não achou ligante: a estrutura não tem hetero-resíduo "
            "(fora água) nem cadeia curta. É um modelo só de proteína?"
        )
    alvo = (escolha["cadeia"], escolha["resnum"], escolha["icode"], escolha["resname"])
    desc = dict(escolha, como="auto")
    return por_res[alvo], desc, cands


# ---------------------------------------------------------------------------
# vizinhança


def _celulas(atomos, lado):
    grade = defaultdict(list)
    for a in atomos:
        grade[(int(math.floor(a.x / lado)), int(math.floor(a.y / lado)),
               int(math.floor(a.z / lado)))].append(a)
    return grade


def contatos(
    atomos: list[Atomo],
    ligante: list[Atomo],
    raio: float = 5.0,
    com_hidrogenios: bool = False,
    incluir_agua: bool = True,
):
    """Tudo a <= raio de qualquer átomo do ligante, agrupado por resíduo.

    Grade de células do tamanho do raio: com um ligante de 50 átomos e uma
    estrutura de 20 mil, a conta é instantânea e não depende de numpy.
    """
    if raio <= 0:
        raise ValueError("raio deve ser > 0")
    alvo = atomos if com_hidrogenios else sem_hidrogenios(atomos)
    lig = ligante if com_hidrogenios else sem_hidrogenios(ligante)
    if not lig:
        raise ValueError("o ligante não tem átomos pesados (só hidrogênio?)")
    chaves_lig = {a.chave_residuo for a in ligante}
    seriais_lig = {id(a) for a in ligante}
    cadeias_lig = {a.cadeia for a in ligante}
    por_res_alvo = agrupar_residuos(alvo)

    grade = _celulas(alvo, raio)
    r2 = raio * raio
    por_res: dict[tuple, dict] = {}
    pares: list[dict] = []

    for la in lig:
        cx, cy, cz = (int(math.floor(la.x / raio)), int(math.floor(la.y / raio)),
                      int(math.floor(la.z / raio)))
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for a in grade.get((cx + dx, cy + dy, cz + dz), ()):
                        if id(a) in seriais_lig or a.chave_residuo in chaves_lig:
                            continue
                        d2 = (a.x - la.x) ** 2 + (a.y - la.y) ** 2 + (a.z - la.z) ** 2
                        if d2 > r2:
                            continue
                        d = math.sqrt(d2)
                        ch = a.chave_residuo
                        pares.append({
                            "atomo_ligante": la.nome, "elemento_ligante": la.elemento,
                            "cadeia": a.cadeia, "resname": a.resname,
                            "resnum": a.resnum, "icode": a.icode.strip(),
                            "atomo": a.nome, "elemento": a.elemento,
                            "distancia": d,
                        })
                        reg = por_res.get(ch)
                        if reg is None:
                            ats_res = por_res_alvo.get(ch, [a])
                            # a categoria é a VERDADEIRA (íon é íon, água é
                            # água); estar na cadeia do ligante é um campo
                            # separado, para não apagar a química do item
                            cat, sub = categorizar(
                                a.resname, len(ats_res),
                                [b.elemento for b in ats_res],
                            )
                            reg = por_res[ch] = {
                                "cadeia": a.cadeia, "resname": a.resname,
                                "resnum": a.resnum, "icode": a.icode.strip(),
                                "residuo": a.rotulo_residuo,
                                "categoria": cat, "subtipo": sub,
                                "aminoacido": AA3_PARA_1.get(a.resname, ""),
                                "mesma_cadeia_que_o_ligante":
                                    int(a.cadeia in cadeias_lig),
                                "n_atomos_no_residuo": len(ats_res),
                                "n_atomos_em_contato": 0,
                                "distancia_min": float("inf"),
                                "atomo_mais_proximo": "", "atomo_ligante_mais_proximo": "",
                                "atomos": set(),
                            }
                        reg["atomos"].add(a.nome)
                        reg["n_atomos_em_contato"] = len(reg["atomos"])
                        if d < reg["distancia_min"]:
                            reg["distancia_min"] = d
                            reg["atomo_mais_proximo"] = a.nome
                            reg["atomo_ligante_mais_proximo"] = la.nome

    lista, descartadas = [], set()
    for reg in por_res.values():
        if not incluir_agua and reg["categoria"] == "agua":
            descartadas.add((reg["cadeia"], reg["resnum"], reg["resname"]))
            continue
        reg = dict(reg)
        reg["atomos"] = ";".join(sorted(reg.pop("atomos")))
        lista.append(reg)
    if descartadas:          # o CSV por átomo segue a mesma regra da lista
        pares = [d for d in pares
                 if (d["cadeia"], d["resnum"], d["resname"]) not in descartadas]
    lista.sort(key=lambda d: (d["distancia_min"], d["cadeia"], d["resnum"]))
    pares.sort(key=lambda d: d["distancia"])
    return lista, pares


def resumir_bolso(lista: list[dict], descricao: dict, raio: float, inv: dict) -> dict:
    """O que o relatório precisa dizer em uma tela."""
    por_cat: dict[str, int] = defaultdict(int)
    por_cadeia: dict[str, int] = defaultdict(int)
    for c in lista:
        por_cat[c["categoria"]] += 1
        por_cadeia[c["cadeia"]] += 1
    cad_polim = {
        c["cadeia"] for c in lista if c["categoria"] in POLIMEROS
    }
    nao_proteicos = [
        {k: c[k] for k in ("residuo", "categoria", "subtipo", "distancia_min",
                           "n_atomos_em_contato", "mesma_cadeia_que_o_ligante")}
        for c in lista if c["categoria"] not in POLIMEROS
    ]
    return {
        "ligante": descricao,
        "raio_angstrom": raio,
        "n_itens_no_raio": len(lista),
        "n_residuos_proteicos": sum(1 for c in lista if c["categoria"] == "proteina"),
        "por_categoria": dict(sorted(por_cat.items())),
        "por_cadeia": dict(sorted(por_cadeia.items())),
        "cadeias_polimericas_em_contato": sorted(cad_polim),
        "bolso_de_interface": len(cad_polim) >= 2,
        "nao_proteicos_no_raio": nao_proteicos,
        "distancia_minima": min((c["distancia_min"] for c in lista), default=None),
        "inventario_estrutura": inv,
    }


def selecao_pymol(lista: list[dict], descricao: dict, raio: float) -> str:
    """Um .pml para olhar o bolso — conferir com o olho pega o que o CSV não mostra."""
    sel = " or ".join(
        f"(chain {c['cadeia']} and resi {c['resnum']}{c['icode']})" for c in lista
    ) or "none"
    lig = (f"(chain {descricao['cadeia']} and resi {descricao['resnum']})"
           if descricao.get("resnum", -1) >= 0 else f"(chain {descricao['cadeia']})")
    return (
        f"# bolso a {raio} A de {descricao.get('resname')} "
        f"({descricao.get('cadeia')}/{descricao.get('resnum')})\n"
        f"select ligante, {lig}\n"
        f"select bolso, {sel}\n"
        "show sticks, ligante\nshow sticks, bolso\nshow spheres, bolso and inorganic\n"
        "color yellow, ligante\ncolor cyan, bolso\norient ligante\n"
    )


def analisar_bolso(
    estrutura,
    resname: str | None = None,
    cadeia: str | None = None,
    resnum: int | None = None,
    auto: bool = False,
    copia: int | None = None,
    raio: float = 5.0,
    modelo: int = 1,
    com_hidrogenios: bool = False,
    incluir_agua: bool = True,
):
    """Pipeline do bolso: ler -> escolher ligante -> listar vizinhança."""
    atomos = ler_estrutura(estrutura, modelo=modelo)
    inv = inventario(atomos)
    lig, desc, cands = selecionar_ligante(
        atomos, resname=resname, cadeia=cadeia, resnum=resnum, auto=auto, copia=copia
    )
    lista, pares = contatos(
        atomos, lig, raio=raio, com_hidrogenios=com_hidrogenios,
        incluir_agua=incluir_agua,
    )
    resumo = resumir_bolso(lista, desc, raio, inv)
    resumo["candidatos_a_ligante"] = cands[:10]
    return {"atomos": atomos, "ligante": lig, "descricao": desc,
            "contatos": lista, "pares": pares, "resumo": resumo,
            "inventario": inv}


# ---------------------------------------------------------------------------


def escrever_csv(caminho: Path, linhas: list[dict]):
    caminho.parent.mkdir(parents=True, exist_ok=True)
    if not linhas:
        caminho.write_text("")
        return caminho
    with caminho.open("w", newline="", encoding="utf-8") as fh:
        esc = csv.DictWriter(fh, fieldnames=list(linhas[0].keys()))
        esc.writeheader()
        for l in linhas:
            esc.writerow({k: (round(v, 3) if isinstance(v, float) else v)
                          for k, v in l.items()})
    return caminho


_escrever_csv = escrever_csv          # nome antigo, mantido


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Lista tudo a menos de N Å do ligante: resíduos, cadeias, "
                    "cofatores, íons, águas.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--estrutura", required=True, type=Path, help="PDB ou mmCIF")
    g = ap.add_argument_group("qual é o ligante")
    g.add_argument("--ligante-resname", default=None, help="ex.: LIG, UNK, 063")
    g.add_argument("--ligante-cadeia", default=None)
    g.add_argument("--ligante-resnum", type=int, default=None)
    g.add_argument("--ligante-copia", type=int, default=None,
                   help="qual cópia usar quando há várias (1 = a primeira)")
    g.add_argument("--ligante-auto", action="store_true",
                   help="pegar o maior hetero-resíduo que não é água/íon")
    ap.add_argument("--raio", type=float, default=5.0)
    ap.add_argument("--modelo", type=int, default=1)
    ap.add_argument("--com-hidrogenios", action="store_true")
    ap.add_argument("--sem-agua", action="store_true", help="não listar águas")
    ap.add_argument("--saida", type=Path, default=Path("bolso"))
    args = ap.parse_args(argv)

    try:
        r = analisar_bolso(
            args.estrutura, resname=args.ligante_resname, cadeia=args.ligante_cadeia,
            resnum=args.ligante_resnum, auto=args.ligante_auto, copia=args.ligante_copia,
            raio=args.raio, modelo=args.modelo, com_hidrogenios=args.com_hidrogenios,
            incluir_agua=not args.sem_agua,
        )
    except ValueError as erro:
        print(f"ERRO: {erro}", file=sys.stderr)
        return 2

    saida = args.saida
    saida.mkdir(parents=True, exist_ok=True)
    escrever_csv(saida / "contatos_residuos.csv", r["contatos"])
    escrever_csv(saida / "contatos_atomos.csv", r["pares"])
    (saida / "bolso_resumo.json").write_text(
        json.dumps(r["resumo"], indent=2, ensure_ascii=False, default=str)
    )
    (saida / "bolso.pml").write_text(
        selecao_pymol(r["contatos"], r["descricao"], args.raio)
    )

    d, s = r["descricao"], r["resumo"]
    print(f"ligante: {d['resname']} cadeia {d['cadeia']} resnum {d['resnum']} "
          f"({d['n_atomos']} átomos, escolha por {d['como']})")
    print(f"a <= {args.raio} Å: {s['n_itens_no_raio']} itens "
          f"({s['n_residuos_proteicos']} resíduos proteicos)")
    for cat, n in s["por_categoria"].items():
        print(f"  {cat:28s} {n}")
    print(f"cadeias poliméricas em contato: "
          f"{', '.join(s['cadeias_polimericas_em_contato']) or '(nenhuma)'}"
          + ("  <- bolso de INTERFACE" if s["bolso_de_interface"] else ""))
    if s["nao_proteicos_no_raio"]:
        print("não-proteicos no raio:")
        for c in s["nao_proteicos_no_raio"]:
            print(f"  {c['residuo']:20s} {c['categoria']:24s} "
                  f"{c['distancia_min']:.2f} Å")
    print(f"-> {saida}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
