#!/usr/bin/env python
"""
Conservação por coluna de um MSA, com pesos de Henikoff.

Só biblioteca padrão (usa numpy se existir, mas não exige). Segundos.

    python msa_henikoff.py --msa alinhamento.fasta --saida conservacao/
    python msa_henikoff.py --msa msa_boltz.a3m --rotulo boltz2 --saida conservacao/

Por que pesar as sequências
---------------------------
Os MSAs que saem dos servidores de cofold (ColabFold/MMseqs2 para o Boltz-2 e
o Chai-1/2, jackhmmer+hhblits para o OpenFold) são buscas em bancos
REDUNDANTES. Dez mil sequências de SARS-CoV-2 ou dez mil ortólogos de primata
não são dez mil observações independentes: são uma observação repetida dez mil
vezes. Contar frequência crua nesses MSAs faz a coluna parecer conservada
porque o banco está enviesado, não porque a evolução a protegeu.

Os pesos de Henikoff & Henikoff (1994), "position-based sequence weights",
corrigem isso sem precisar de árvore filogenética nem de limiar de identidade:
para cada coluna j, um resíduo do tipo `a` que aparece `s` vezes entre `r`
tipos distintos contribui com

    1 / (r_j * s_ja)

para o peso da sequência que o carrega. Uma sequência que só repete o que todo
mundo já diz acumula contribuições pequenas; uma sequência que traz um tipo
raro acumula contribuições grandes. Somando sobre as colunas sai o peso da
sequência.

O efeito prático é medível e este módulo o relata: `N/Neff` é o fator de
redundância do MSA. Um MSA com 8000 sequências e Neff 120 tem a informação de
~120, e qualquer conservação lida dele sem pesos está inflada.

Decisões que não são óbvias, e por que foram tomadas assim
----------------------------------------------------------
1. GAP não é um tipo de resíduo (padrão). Henikoff original conta tipos de
   aminoácido; tratar "-" como 21º tipo faz um bloco de sequências truncadas
   virar "tipo raro" e ganhar peso alto por não ter sido sequenciado. Quem
   quiser o contrário tem `--gap-como-tipo`.

2. A soma das contribuições é dividida pelo NÚMERO DE COLUNAS QUE A SEQUÊNCIA
   DE FATO COBRE (`--normalizacao media`, padrão). Os MSAs desses servidores
   são cheios de hits parciais; na soma crua um fragmento de 40 resíduos num
   alinhamento de 400 fica com peso ~10x menor só por ser curto, e não é isso
   que Henikoff mede. Com `--normalizacao soma` sai a soma clássica.

3. Os pesos são escalados para somar N (o número de sequências). Assim o peso
   uniforme vale 1,0 e o peso de cada sequência se lê direto: 0,2 = "esta
   sequência vale 1/5 de uma observação independente".

4. O escore principal relatado é a JSD ao fundo de Robinson & Robinson
   (Capra & Singh 2007, o método que ganhou a comparação deles), não a entropia
   de Shannon pura: a JSD sabe que uma coluna toda de triptofano é mais
   informativa que uma coluna toda de leucina, porque W é raro e L é comum.
   A entropia, a conservação normalizada e a KL saem também, por coluna.

5. Coluna com pouca gente é marcada, não silenciada: `ocupacao_pond` e
   `neff_coluna` entram no CSV e `confiavel` resume os dois. Uma coluna com
   duas sequências não-gap é "100% conservada" por acidente.

Formatos de entrada aceitos (o `--formato auto` decide pela extensão/conteúdo)
-----------------------------------------------------------------------------
  .fasta/.fa/.afa/.aln/.mfa  FASTA ALINHADO (todas as sequências do mesmo
                             comprimento). É o formato recomendado.
  .a3m                       o que o ColabFold/MMseqs2 entrega ao Boltz e ao
                             Chai; minúsculas são inserções e são removidas
                             (vira FASTA alinhado nas colunas da query).
  .a2m                       idem .a3m
  .sto/.stockholm            o que o jackhmmer/hmmer entrega ao OpenFold/AF
  .csv/.tsv                  o CSV de MSA do Boltz (colunas key,sequence)
  .pqt/.parquet              o `.aligned.pqt` do Chai-1/2 (exige pandas)

A PRIMEIRA sequência é tratada como a query/referência — é a convenção de
todos esses formatos, e é dela que vem o mapeamento para a estrutura.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

try:                                     # opcional: só acelera
    import numpy as _np
except Exception:                        # pragma: no cover - ambiente sem numpy
    _np = None

ALFABETO = "ACDEFGHIKLMNPQRSTVWY"        # os 20 canônicos, e só eles contam
GAP = "-"
DESCONHECIDO = "X"
LOG2_20 = math.log2(20)

# Frequências de fundo de Robinson & Robinson (1991), as mesmas que o
# score_conservation de Capra & Singh usa. Normalizadas em tempo de import
# porque a tabela publicada soma 1,003.
_FUNDO_CRU = {
    "A": 0.0780, "R": 0.0512, "N": 0.0448, "D": 0.0536, "C": 0.0193,
    "Q": 0.0408, "E": 0.0632, "G": 0.0732, "H": 0.0229, "I": 0.0534,
    "L": 0.0905, "K": 0.0594, "M": 0.0220, "F": 0.0387, "P": 0.0520,
    "S": 0.0683, "T": 0.0585, "W": 0.0136, "Y": 0.0323, "V": 0.0673,
}
_S = sum(_FUNDO_CRU.values())
FUNDO = {a: _FUNDO_CRU[a] / _S for a in ALFABETO}

# ---------------------------------------------------------------------------
# leitura


@dataclass
class MSA:
    """Um alinhamento já normalizado: maiúsculas, gap `-`, resto `X`."""

    ids: list[str]
    seqs: list[str]
    caminho: Path | None = None
    formato: str = "?"
    rotulo: str = "msa"
    avisos: list[str] = field(default_factory=list)
    descartadas: list[str] = field(default_factory=list)

    @property
    def n(self) -> int:
        return len(self.seqs)

    @property
    def largura(self) -> int:
        return len(self.seqs[0]) if self.seqs else 0

    @property
    def query(self) -> str:
        return self.seqs[0] if self.seqs else ""

    @property
    def query_sem_gap(self) -> str:
        return self.query.replace(GAP, "")

    def colunas(self) -> list[str]:
        return ["".join(c) for c in zip(*self.seqs)] if self.seqs else []


def _normalizar(seq: str) -> str:
    """Maiúsculas; '.'/'~' -> '-'; qualquer outra coisa -> 'X'.

    Depois disso o alfabeto é exatamente ALFABETO + '-' + 'X', o que garante
    ASCII puro (o caminho rápido com numpy depende disso) e tira do caminho
    as letras ambíguas (B, Z, J, U, O, *) sem fingir que são resíduos.
    """
    saida = []
    for ch in seq.upper():
        if ch in (".", "~", GAP):
            saida.append(GAP)
        elif ch in ALFABETO:
            saida.append(ch)
        elif ch.isspace():
            continue
        else:
            saida.append(DESCONHECIDO)
    return "".join(saida)


def _ler_fasta_cru(texto: str) -> tuple[list[str], list[str]]:
    """FASTA -> (ids, sequências) SEM normalizar (o caso .a3m precisa do caso)."""
    ids: list[str] = []
    seqs: list[str] = []
    atual: list[str] = []
    for linha in texto.splitlines():
        if linha.startswith(">"):
            if ids:
                seqs.append("".join(atual))
            ids.append(linha[1:].strip() or f"seq{len(ids) + 1}")
            atual = []
        elif linha.startswith((";", "#")):      # comentário de FASTA antigo
            continue
        elif ids:
            atual.append(linha.strip())
        elif linha.strip():
            raise ValueError(
                "o arquivo tem sequência antes de qualquer cabeçalho '>' — "
                "não é FASTA"
            )
    if ids:
        seqs.append("".join(atual))
    return ids, seqs


def _remover_insercoes(ids: list[str], seqs: list[str]) -> tuple[list[str], list[str], list[str]]:
    """a3m/a2m -> FASTA alinhado: minúscula é inserção, cai fora.

    Em a3m as inserções em relação à query não ocupam coluna: a query não tem
    nenhuma minúscula e, tirando as minúsculas de cada sequência, todas voltam
    ao comprimento da query. Se não voltam, o arquivo não é a3m — e é melhor
    dizer isso do que analisar um alinhamento torto.
    """
    limpas = ["".join(c for c in s if not c.islower()) for s in seqs]
    avisos = []
    n_ins = sum(1 for a, b in zip(seqs, limpas) if len(a) != len(b))
    if n_ins:
        avisos.append(f"{n_ins} sequências tinham inserções (minúsculas) removidas")
    return ids, limpas, avisos


def _ler_stockholm(texto: str) -> tuple[list[str], list[str]]:
    """Stockholm do hmmer/jackhmmer, inclusive quebrado em blocos."""
    ordem: list[str] = []
    partes: dict[str, list[str]] = {}
    for linha in texto.splitlines():
        s = linha.strip()
        if not s or s.startswith("#") or s == "//":
            continue
        campos = s.split(None, 1)
        if len(campos) != 2:
            continue
        nome, pedaco = campos[0], campos[1].replace(" ", "")
        if nome not in partes:
            partes[nome] = []
            ordem.append(nome)
        partes[nome].append(pedaco)
    return ordem, ["".join(partes[n]) for n in ordem]


def _ler_tabela(caminho: Path, col_id: str | None, col_seq: str | None):
    """CSV/TSV de MSA (o formato que o Boltz aceita: key,sequence)."""
    delim = "\t" if caminho.suffix.lower() in (".tsv", ".tab") else ","
    with caminho.open(newline="", encoding="utf-8", errors="replace") as fh:
        leitor = csv.DictReader(fh, delimiter=delim)
        campos = [c.strip() for c in (leitor.fieldnames or [])]
        if not campos:
            raise ValueError(f"{caminho}: tabela sem cabeçalho")
        alvo_seq = col_seq or next(
            (c for c in campos if c.lower() in ("sequence", "seq", "sequencia", "aligned_sequence")),
            None,
        )
        if alvo_seq is None:
            raise ValueError(
                f"{caminho}: não achei a coluna de sequência (cabeçalho: {campos}). "
                "Passe --coluna-seq."
            )
        alvo_id = col_id or next(
            (c for c in campos if c.lower() in ("key", "id", "name", "header", "description")),
            None,
        )
        ids, seqs = [], []
        for i, linha in enumerate(leitor):
            s = (linha.get(alvo_seq) or "").strip()
            if not s:
                continue
            seqs.append(s)
            ids.append(str(linha.get(alvo_id) or f"seq{i + 1}") if alvo_id else f"seq{i + 1}")
    return ids, seqs


def _ler_parquet(caminho: Path, col_id: str | None, col_seq: str | None):
    """`.aligned.pqt` do Chai-1/2. Exige pandas + pyarrow (ou fastparquet).

    Os dois erros possíveis (pandas ausente, engine de parquet ausente) viram
    ValueError com a saída: todo o resto do script trata ValueError como "erro
    de entrada do usuário" e imprime a mensagem em vez de um traceback.
    """
    try:
        import pandas as pd
    except Exception as erro:                   # pragma: no cover
        raise ValueError(
            f"{caminho}: ler parquet exige pandas neste env ({erro}). "
            "Alternativa: converta o .pqt para FASTA alinhado antes."
        ) from erro
    try:
        df = pd.read_parquet(caminho)
    except Exception as erro:
        raise ValueError(
            f"{caminho}: não deu para ler o parquet ({type(erro).__name__}: "
            f"{str(erro).splitlines()[0]}). Instale o engine "
            "(`pip install pyarrow`) ou exporte o MSA como FASTA alinhado."
        ) from erro
    cols = [str(c) for c in df.columns]
    alvo_seq = col_seq or next(
        (c for c in cols if c.lower() in ("sequence", "seq", "aligned_sequence")), None
    )
    if alvo_seq is None:
        raise ValueError(f"{caminho}: não achei a coluna de sequência (colunas: {cols})")
    alvo_id = col_id or next(
        (c for c in cols if c.lower() in ("key", "id", "name", "source_database", "header")), None
    )
    seqs = [str(x) for x in df[alvo_seq].tolist()]
    ids = [str(x) for x in df[alvo_id].tolist()] if alvo_id else [
        f"seq{i + 1}" for i in range(len(seqs))
    ]
    return ids, seqs


def detectar_formato(caminho: Path) -> str:
    """Extensão primeiro; se ela não disser nada, o conteúdo diz."""
    ext = caminho.suffix.lower()
    mapa = {
        ".a3m": "a3m", ".a2m": "a3m",
        ".fasta": "fasta", ".fa": "fasta", ".afa": "fasta", ".aln": "fasta",
        ".mfa": "fasta", ".fas": "fasta", ".fst": "fasta", ".seq": "fasta",
        ".sto": "stockholm", ".stockholm": "stockholm", ".sth": "stockholm",
        ".csv": "csv", ".tsv": "csv", ".tab": "csv",
        ".pqt": "parquet", ".parquet": "parquet",
    }
    if ext in mapa:
        return mapa[ext]
    try:
        with caminho.open("rb") as fh:
            cabeca = fh.read(4096).decode("utf-8", errors="replace")
    except OSError as erro:
        raise ValueError(f"{caminho}: não deu para ler ({erro})") from erro
    if cabeca.startswith("# STOCKHOLM"):
        return "stockholm"
    if cabeca.lstrip().startswith(">"):
        return "a3m" if any(c.islower() for c in cabeca) else "fasta"
    primeira = next(iter(cabeca.splitlines()), "")
    if "," in primeira or "\t" in primeira:
        return "csv"
    raise ValueError(
        f"{caminho}: não reconheci o formato. Use --formato "
        "{fasta,a3m,stockholm,csv,parquet}."
    )


def ler_msa(
    caminho,
    formato: str = "auto",
    rotulo: str | None = None,
    max_seqs: int | None = None,
    semente: int = 0,
    coluna_id: str | None = None,
    coluna_seq: str | None = None,
) -> MSA:
    """Lê, normaliza e devolve um MSA pronto para pesar.

    `max_seqs` amostra (com semente fixa, e a query SEMPRE fica) para MSAs
    gigantes. Amostrar é honesto aqui porque os pesos de Henikoff já são uma
    média sobre sequências — mas o resumo registra que houve amostragem.
    """
    caminho = Path(caminho)
    if not caminho.exists():
        raise ValueError(f"MSA não encontrado: {caminho}")
    if caminho.is_dir():
        raise ValueError(f"{caminho} é um diretório, não um arquivo de MSA")
    if caminho.stat().st_size == 0:
        raise ValueError(f"MSA vazio (0 bytes): {caminho}")

    fmt = detectar_formato(caminho) if formato in (None, "auto") else formato
    avisos: list[str] = []

    if fmt == "parquet":
        ids, cruas = _ler_parquet(caminho, coluna_id, coluna_seq)
    elif fmt == "csv":
        ids, cruas = _ler_tabela(caminho, coluna_id, coluna_seq)
    else:
        texto = caminho.read_text(encoding="utf-8", errors="replace")
        if fmt == "stockholm":
            ids, cruas = _ler_stockholm(texto)
        else:
            ids, cruas = _ler_fasta_cru(texto)
    if not ids:
        raise ValueError(f"{caminho}: nenhuma sequência lida (formato {fmt})")

    if fmt in ("a3m", "csv", "parquet") and any(
        c.islower() for s in cruas[: min(len(cruas), 200)] for c in s
    ):
        ids, cruas, av = _remover_insercoes(ids, cruas)
        avisos.extend(av)
        if fmt != "a3m":
            avisos.append(f"havia minúsculas em {fmt}: tratadas como inserções (estilo a3m)")

    seqs = [_normalizar(s) for s in cruas]

    # comprimentos: é AQUI que um "FASTA alinhado" que não está alinhado morre
    comps = Counter(len(s) for s in seqs)
    if len(comps) > 1:
        esperado = comps.most_common(1)[0][0]
        culpados = [
            f"{ids[i]!r} ({len(s)})" for i, s in enumerate(seqs) if len(s) != esperado
        ][:5]
        fora = len(seqs) - comps[esperado]       # sum(), não Counter.total(): 3.9
        raise ValueError(
            f"{caminho}: o alinhamento não é retangular — a maioria tem "
            f"{esperado} colunas, mas {fora} sequências "
            f"não: {', '.join(culpados)}"
            + (" ..." if fora > 5 else "")
            + ". FASTA ALINHADO significa todas do mesmo comprimento; "
            "se este arquivo é a3m, passe --formato a3m."
        )

    # sequências sem nenhum resíduo não carregam informação e quebrariam a
    # normalização por média; saem com registro, não em silêncio
    guardadas_ids, guardadas, descartadas = [], [], []
    for i, s in enumerate(seqs):
        if any(c in ALFABETO for c in s):
            guardadas_ids.append(ids[i])
            guardadas.append(s)
        else:
            descartadas.append(ids[i])
    if descartadas:
        avisos.append(
            f"{len(descartadas)} sequências sem nenhum resíduo canônico foram "
            f"descartadas (ex.: {descartadas[0]!r})"
        )
    if not guardadas:
        raise ValueError(f"{caminho}: nenhuma sequência com resíduo canônico")

    # ids repetidos: comuns em a3m de servidor, e atrapalham o CSV de pesos
    vistos: dict[str, int] = {}
    unicos = []
    for i in guardadas_ids:
        if i in vistos:
            vistos[i] += 1
            unicos.append(f"{i}__{vistos[i]}")
        else:
            vistos[i] = 0
            unicos.append(i)
    n_rep = sum(1 for v in vistos.values() if v)
    if n_rep:
        avisos.append(f"{n_rep} identificadores repetidos ganharam sufixo __n")

    if max_seqs and len(guardadas) > max_seqs:
        rng = random.Random(semente)
        resto = list(range(1, len(guardadas)))
        rng.shuffle(resto)
        indices = [0] + sorted(resto[: max_seqs - 1])
        avisos.append(
            f"amostradas {max_seqs} de {len(guardadas)} sequências "
            f"(semente {semente}; a query ficou)"
        )
        unicos = [unicos[i] for i in indices]
        guardadas = [guardadas[i] for i in indices]

    return MSA(
        ids=unicos,
        seqs=guardadas,
        caminho=caminho,
        formato=fmt,
        rotulo=rotulo or caminho.stem,
        avisos=avisos,
        descartadas=descartadas,
    )


def validar_msa(msa: MSA, minimo_sequencias: int = 10, permitir_unica: bool = False):
    """(erros, avisos). Erro = não dá para medir conservação; aviso = cuidado.

    O caso da sequência única não é um detalhe: o ESMFold é single-sequence por
    construção. Pedir "conservação do MSA do ESMFold" é uma pergunta sem
    objeto, e o lugar de dizer isso é aqui, não no meio do relatório.
    """
    erros, avisos = [], list(msa.avisos)
    if msa.n == 0:
        erros.append("MSA sem sequências")
        return erros, avisos
    if msa.largura == 0:
        erros.append("MSA com largura 0")
        return erros, avisos
    if msa.n == 1:
        msg = (
            "MSA com UMA sequência: não existe conservação para medir. "
            "Modelos single-sequence (ESMFold, e Boltz/Chai rodados com "
            "--use_msa_server desligado) não produzem MSA — use o MSA do "
            "OpenFold/ColabFold para a mesma cadeia se quiser o perfil."
        )
        (avisos if permitir_unica else erros).append(msg)
        if not permitir_unica:
            return erros, avisos
    elif msa.n < minimo_sequencias:
        avisos.append(
            f"MSA raso: {msa.n} sequências. Conservação por coluna com menos de "
            f"{minimo_sequencias} sequências é anedota, não estatística."
        )

    gaps_query = msa.query.count(GAP)
    if gaps_query:
        frac = gaps_query / msa.largura
        avisos.append(
            f"a query (1ª sequência, {msa.ids[0]!r}) tem {gaps_query} gaps "
            f"({frac:.1%} das colunas)"
            + (" — confira se a 1ª sequência é mesmo a query" if frac > 0.05 else "")
        )
    n_x = sum(s.count(DESCONHECIDO) for s in msa.seqs)
    if n_x:
        avisos.append(
            f"{n_x} posições com caractere fora dos 20 canônicos (B/Z/J/X/U/O/*): "
            "contam como 'X' e ficam fora das frequências"
        )
    todas_gap = sum(1 for c in msa.colunas() if set(c) <= {GAP, DESCONHECIDO})
    if todas_gap:
        avisos.append(f"{todas_gap} colunas sem nenhum resíduo (só gap/X)")
    return erros, avisos


# ---------------------------------------------------------------------------
# pesos de Henikoff


def _matriz_bytes(seqs: list[str]):
    return _np.frombuffer("".join(seqs).encode("ascii"), dtype=_np.uint8).reshape(
        len(seqs), len(seqs[0])
    )


def _pesos_cru_python(seqs, validos, colunas_usar):
    n = len(seqs)
    soma = [0.0] * n
    contadas = [0] * n
    for j in colunas_usar:
        col = [s[j] for s in seqs]
        cnt = Counter(col)
        presentes = {c: k for c, k in cnt.items() if c in validos}
        r = len(presentes)
        if r == 0:
            continue
        contrib = {c: 1.0 / (r * k) for c, k in presentes.items()}
        for i, ch in enumerate(col):
            v = contrib.get(ch)
            if v is not None:
                soma[i] += v
                contadas[i] += 1
    return soma, contadas


def _pesos_cru_numpy(seqs, validos, colunas_usar):
    m = _matriz_bytes(seqs)
    n = m.shape[0]
    soma = _np.zeros(n, dtype=_np.float64)
    contadas = _np.zeros(n, dtype=_np.int64)
    codigos = [ord(c) for c in validos]
    for j in colunas_usar:
        col = m[:, j]
        cnt = _np.bincount(col, minlength=256)
        presentes = [c for c in codigos if cnt[c]]
        r = len(presentes)
        if r == 0:
            continue
        lut = _np.zeros(256, dtype=_np.float64)
        marca = _np.zeros(256, dtype=_np.int64)
        for c in presentes:
            lut[c] = 1.0 / (r * int(cnt[c]))
            marca[c] = 1
        soma += lut[col]
        contadas += marca[col]
    return soma.tolist(), contadas.tolist()


def pesos_henikoff(
    msa: MSA,
    gap_como_tipo: bool = False,
    normalizacao: str = "media",
    so_colunas_query: bool = False,
    usar_numpy: bool = True,
):
    """Pesos posicionais de Henikoff & Henikoff (1994), somando N.

    Devolve (pesos, info). `info` guarda tudo que o relatório precisa citar:
    Neff de Kish, redundância, quantas colunas entraram, e os extremos.
    """
    if normalizacao not in ("media", "soma"):
        raise ValueError("normalizacao deve ser 'media' ou 'soma'")
    n = msa.n
    if n == 0:
        raise ValueError("MSA sem sequências")
    validos = set(ALFABETO) | ({GAP} if gap_como_tipo else set())
    colunas_usar = [
        j for j in range(msa.largura)
        if not so_colunas_query or msa.query[j] != GAP
    ]
    if n == 1:
        return [1.0], {
            "n_seqs": 1, "neff": 1.0, "redundancia": 1.0,
            "colunas_usadas": len(colunas_usar), "normalizacao": normalizacao,
            "gap_como_tipo": gap_como_tipo, "so_colunas_query": so_colunas_query,
            "motor": "-", "uniforme": True,
            "peso_min": 1.0, "peso_max": 1.0, "peso_query": 1.0,
            "motivo": "sequência única: peso 1,0 e nada a corrigir",
        }

    motor = "numpy" if (usar_numpy and _np is not None) else "python"
    calc = _pesos_cru_numpy if motor == "numpy" else _pesos_cru_python
    soma, contadas = calc(msa.seqs, validos, colunas_usar)

    if normalizacao == "media":
        cru = [s / c if c else 0.0 for s, c in zip(soma, contadas)]
    else:
        cru = list(soma)

    total = sum(cru)
    uniforme = False
    motivo = ""
    if total <= 0:
        cru = [1.0] * n
        total = float(n)
        uniforme = True
        motivo = "nenhuma coluna informativa: caí para pesos uniformes"
    pesos = [w * n / total for w in cru]

    sq = sum(w * w for w in pesos)
    neff = (sum(pesos) ** 2) / sq if sq > 0 else float(n)   # Kish
    info = {
        "n_seqs": n,
        "neff": neff,
        "redundancia": n / neff if neff > 0 else float("nan"),
        "colunas_usadas": len(colunas_usar),
        "normalizacao": normalizacao,
        "gap_como_tipo": gap_como_tipo,
        "so_colunas_query": so_colunas_query,
        "motor": motor,
        "uniforme": uniforme,
        "peso_min": min(pesos),
        "peso_max": max(pesos),
        "peso_query": pesos[0],
        "motivo": motivo,
    }
    return pesos, info


def pesos_uniformes(msa: MSA):
    """Para comparar: é o que a conservação 'crua' usa, sem saber que usa."""
    pesos = [1.0] * msa.n
    return pesos, {
        "n_seqs": msa.n, "neff": float(msa.n), "redundancia": 1.0,
        "colunas_usadas": msa.largura, "normalizacao": "uniforme",
        "gap_como_tipo": False, "motor": "-", "uniforme": True,
        "peso_min": 1.0, "peso_max": 1.0, "peso_query": 1.0, "motivo": "",
    }


# ---------------------------------------------------------------------------
# escores por coluna


def entropia_bits(p: dict[str, float]) -> float:
    return -sum(v * math.log2(v) for v in p.values() if v > 0)


def kl_bits(p: dict[str, float], q: dict[str, float] = FUNDO) -> float:
    return sum(v * math.log2(v / q[a]) for a, v in p.items() if v > 0)


def jsd(p: dict[str, float], q: dict[str, float] = FUNDO) -> float:
    """Jensen-Shannon divergence em base 2: 0 (igual ao fundo) a 1 (máxima).

    É a medida de Capra & Singh (2007). Simétrica, limitada, e penaliza menos
    a coluna conservada num resíduo COMUM do que a entropia penalizaria.
    """
    total = 0.0
    for a in ALFABETO:
        pa, qa = p.get(a, 0.0), q[a]
        m = 0.5 * (pa + qa)
        if m <= 0:
            continue
        if pa > 0:
            total += 0.5 * pa * math.log2(pa / m)
        if qa > 0:
            total += 0.5 * qa * math.log2(qa / m)
    return max(0.0, min(1.0, total))


def _contagens_coluna_python(seqs, pesos, j):
    """(freq ponderada por resíduo, peso_gap, peso_X, n_bruto, soma dos w²).

    A soma dos quadrados sai daqui porque é o que dá o Neff da coluna, e
    percorrer a coluna uma segunda vez só para isso dobraria o custo.
    """
    f: dict[str, float] = {}
    peso_gap = peso_x = soma_q = 0.0
    bruta = 0
    for s, w in zip(seqs, pesos):
        ch = s[j]
        if ch == GAP:
            peso_gap += w
        elif ch == DESCONHECIDO:
            peso_x += w
        else:
            f[ch] = f.get(ch, 0.0) + w
            soma_q += w * w
            bruta += 1
    return f, peso_gap, peso_x, bruta, soma_q


def _contagens_coluna_numpy(m, pesos_np, pesos2_np, j):
    col = m[:, j]
    acc = _np.bincount(col, weights=pesos_np, minlength=256)
    acc2 = _np.bincount(col, weights=pesos2_np, minlength=256)
    bruto = _np.bincount(col, minlength=256)
    f, bruta, soma_q = {}, 0, 0.0
    for a in ALFABETO:
        c = ord(a)
        if bruto[c]:
            f[a] = float(acc[c])
            bruta += int(bruto[c])
            soma_q += float(acc2[c])
    return (f, float(acc[ord(GAP)]), float(acc[ord(DESCONHECIDO)]), bruta, soma_q)


def perfil_colunas(
    msa: MSA,
    pesos: list[float],
    meia_janela: int = 3,
    peso_janela: float = 0.5,
    min_ocupacao: float = 0.5,
    min_neff_coluna: float = 3.0,
    usar_numpy: bool = True,
) -> list[dict]:
    """Uma linha por coluna do alinhamento, com todos os escores.

    `jsd_janela` é a média ponderada com a vizinhança (Capra & Singh mostram
    que isso melhora a detecção de sítio funcional: o sinal de um bolso é
    local, não pontual). `confiavel` cruza ocupação e Neff da coluna.
    """
    n = msa.n
    soma_pesos = sum(pesos) or 1.0
    usar_np = usar_numpy and _np is not None
    if usar_np:
        m = _matriz_bytes(msa.seqs)
        pesos_np = _np.asarray(pesos, dtype=_np.float64)
        pesos2_np = pesos_np * pesos_np

    linhas: list[dict] = []
    pos_query = 0
    for j in range(msa.largura):
        if usar_np:
            f, peso_gap, peso_x, bruta, soma_q = _contagens_coluna_numpy(
                m, pesos_np, pesos2_np, j
            )
        else:
            f, peso_gap, peso_x, bruta, soma_q = _contagens_coluna_python(
                msa.seqs, pesos, j
            )

        peso_res = sum(f.values())
        p = {a: v / peso_res for a, v in f.items()} if peso_res > 0 else {}
        res_q = msa.query[j]
        if res_q != GAP:
            pos_query += 1

        if p:
            h = entropia_bits(p)
            cons = max(0.0, min(1.0, 1.0 - h / LOG2_20))
            d_jsd = jsd(p)
            d_kl = kl_bits(p)
            mais, frac_mais = max(p.items(), key=lambda kv: kv[1])
        else:
            h, cons, d_jsd, d_kl, mais, frac_mais = 0.0, 0.0, 0.0, 0.0, "", 0.0

        # Neff da coluna: tamanho de amostra efetivo entre os NÃO-gap (Kish).
        # É o número que diz se "100% conservada" vale algo.
        neff_col = (peso_res * peso_res / soma_q) if soma_q > 0 else 0.0

        ocup_pond = peso_res / soma_pesos
        linhas.append({
            "coluna": j + 1,
            "pos_query": pos_query if res_q != GAP else "",
            "res_query": res_q,
            "profundidade": bruta,
            "ocupacao": bruta / n if n else 0.0,
            "ocupacao_pond": ocup_pond,
            "frac_gap_pond": peso_gap / soma_pesos,
            "frac_x_pond": peso_x / soma_pesos,
            "n_tipos": len(f),
            "res_mais_frequente": mais,
            "frac_mais_frequente": frac_mais,
            "entropia_bits": h,
            "conservacao_shannon": cons,
            "kl_bits": d_kl,
            "jsd": d_jsd,
            "neff_coluna": neff_col,
            "confiavel": int(ocup_pond >= min_ocupacao and neff_col >= min_neff_coluna),
        })

    suavizado = janela_movel([l["jsd"] for l in linhas], meia_janela, peso_janela)
    for l, s in zip(linhas, suavizado):
        l["jsd_janela"] = s
        # penalizar gap é o que evita que uma coluna de 3 sequências lidere o
        # ranking de conservação do bolso
        l["jsd_penalizado"] = s * l["ocupacao_pond"]
    return linhas


def janela_movel(valores: list[float], meia_janela: int = 3, peso_vizinhos: float = 0.5):
    """(1-λ)·x_i + λ·média(vizinhos), como em Capra & Singh (2007)."""
    if meia_janela <= 0 or peso_vizinhos <= 0 or not valores:
        return list(valores)
    n = len(valores)
    saida = []
    for i in range(n):
        ini, fim = max(0, i - meia_janela), min(n, i + meia_janela + 1)
        viz = [valores[k] for k in range(ini, fim) if k != i]
        media = sum(viz) / len(viz) if viz else valores[i]
        saida.append((1 - peso_vizinhos) * valores[i] + peso_vizinhos * media)
    return saida


# ---------------------------------------------------------------------------
# resumo / veredito


def resumo_msa(msa: MSA, info_pesos: dict, linhas: list[dict]) -> dict:
    """O cartão de confiabilidade do MSA — o que decide se vale confiar nele."""
    l_query = len(msa.query_sem_gap)
    neff = info_pesos["neff"]
    neff_por_res = neff / l_query if l_query else 0.0
    confiaveis = sum(l["confiavel"] for l in linhas)
    colunas_query = [l for l in linhas if l["res_query"] != GAP]
    conf_query = sum(l["confiavel"] for l in colunas_query)
    ocup_media = (
        sum(l["ocupacao_pond"] for l in colunas_query) / len(colunas_query)
        if colunas_query else 0.0
    )

    if msa.n == 1:
        veredito = "SEM MSA (single-sequence): conservação não se aplica"
    elif neff < 10:
        veredito = "RASO: Neff < 10 — conservação por coluna não é confiável"
    elif neff_por_res < 1:
        veredito = ("POBRE: Neff/resíduo < 1 — é a faixa em que o cofold também "
                    "erra; trate a conservação como indicativa")
    elif neff_por_res < 10:
        veredito = "RAZOÁVEL: Neff/resíduo entre 1 e 10"
    else:
        veredito = "PROFUNDO: Neff/resíduo >= 10"

    return {
        "rotulo": msa.rotulo,
        "arquivo": str(msa.caminho) if msa.caminho else "",
        "formato": msa.formato,
        "n_sequencias": msa.n,
        "n_descartadas": len(msa.descartadas),
        "largura_alinhamento": msa.largura,
        "comprimento_query": l_query,
        "id_query": msa.ids[0] if msa.ids else "",
        "neff_henikoff": neff,
        "neff_por_residuo": neff_por_res,
        "redundancia_n_sobre_neff": info_pesos["redundancia"],
        "peso_min": info_pesos["peso_min"],
        "peso_max": info_pesos["peso_max"],
        "peso_query": info_pesos["peso_query"],
        "colunas_confiaveis": confiaveis,
        "colunas_query_confiaveis": conf_query,
        "ocupacao_media_colunas_query": ocup_media,
        "jsd_media_colunas_query": (
            sum(l["jsd"] for l in colunas_query) / len(colunas_query)
            if colunas_query else 0.0
        ),
        "veredito": veredito,
        "avisos": msa.avisos,
    }


def analisar_msa(
    caminho,
    rotulo: str | None = None,
    formato: str = "auto",
    gap_como_tipo: bool = False,
    normalizacao: str = "media",
    so_colunas_query: bool = False,
    meia_janela: int = 3,
    peso_janela: float = 0.5,
    min_ocupacao: float = 0.5,
    min_neff_coluna: float = 3.0,
    max_seqs: int | None = None,
    semente: int = 0,
    permitir_unica: bool = False,
    comparar_sem_peso: bool = True,
    coluna_id: str | None = None,
    coluna_seq: str | None = None,
):
    """Pipeline de um MSA: ler -> validar -> pesar -> perfil -> resumo.

    Devolve um dict com tudo. Erro de validação virá como ValueError, com a
    mensagem que o usuário precisa para arrumar o arquivo.
    """
    msa = ler_msa(
        caminho, formato=formato, rotulo=rotulo, max_seqs=max_seqs,
        semente=semente, coluna_id=coluna_id, coluna_seq=coluna_seq,
    )
    erros, avisos = validar_msa(msa, permitir_unica=permitir_unica)
    msa.avisos = avisos
    if erros:
        raise ValueError(f"{caminho}: " + " | ".join(erros))

    pesos, info = pesos_henikoff(
        msa, gap_como_tipo=gap_como_tipo, normalizacao=normalizacao,
        so_colunas_query=so_colunas_query,
    )
    linhas = perfil_colunas(
        msa, pesos, meia_janela=meia_janela, peso_janela=peso_janela,
        min_ocupacao=min_ocupacao, min_neff_coluna=min_neff_coluna,
    )
    resumo = resumo_msa(msa, info, linhas)

    if comparar_sem_peso and msa.n > 1:
        p_u, i_u = pesos_uniformes(msa)
        linhas_u = perfil_colunas(
            msa, p_u, meia_janela=meia_janela, peso_janela=peso_janela,
            min_ocupacao=min_ocupacao, min_neff_coluna=min_neff_coluna,
        )
        for l, lu in zip(linhas, linhas_u):
            l["jsd_sem_peso"] = lu["jsd"]
            l["conservacao_shannon_sem_peso"] = lu["conservacao_shannon"]
        cq = [l for l in linhas if l["res_query"] != GAP]
        resumo["jsd_media_sem_peso"] = (
            sum(l["jsd_sem_peso"] for l in cq) / len(cq) if cq else 0.0
        )
        resumo["efeito_dos_pesos_jsd"] = (
            resumo["jsd_media_colunas_query"] - resumo["jsd_media_sem_peso"]
        )

    return {"msa": msa, "pesos": pesos, "info_pesos": info,
            "colunas": linhas, "resumo": resumo}


# ---------------------------------------------------------------------------
# escrita


def escrever_csv(caminho: Path, linhas: list[dict], colunas: list[str] | None = None):
    caminho.parent.mkdir(parents=True, exist_ok=True)
    if not linhas:
        caminho.write_text("")
        return caminho
    campos = colunas or list(linhas[0].keys())
    with caminho.open("w", newline="", encoding="utf-8") as fh:
        esc = csv.DictWriter(fh, fieldnames=campos, extrasaction="ignore")
        esc.writeheader()
        for l in linhas:
            esc.writerow({
                c: (round(v, 6) if isinstance(v, float) else v)
                for c, v in l.items() if c in campos
            })
    return caminho


def escrever_pesos(caminho: Path, msa: MSA, pesos: list[float]):
    linhas = [
        {"i": i + 1, "id": ident, "peso": w,
         "residuos": sum(1 for c in seq if c in ALFABETO),
         "frac_gap": seq.count(GAP) / len(seq) if seq else 0.0,
         "identidade_com_query": identidade(msa.query, seq)}
        for i, (ident, seq, w) in enumerate(zip(msa.ids, msa.seqs, pesos))
    ]
    return escrever_csv(caminho, linhas)


def identidade(a: str, b: str) -> float:
    """Identidade sobre as colunas em que AMBAS têm resíduo."""
    pares = [(x, y) for x, y in zip(a, b) if x in ALFABETO and y in ALFABETO]
    return sum(1 for x, y in pares if x == y) / len(pares) if pares else 0.0


# ---------------------------------------------------------------------------


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Conservação por coluna de um MSA, com pesos de Henikoff.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--msa", required=True, type=Path, help="FASTA alinhado, a3m, sto, csv ou pqt")
    ap.add_argument("--rotulo", default=None, help="nome curto (ex.: boltz2, chai1, openfold)")
    ap.add_argument("--formato", default="auto",
                    choices=["auto", "fasta", "a3m", "stockholm", "csv", "parquet"])
    ap.add_argument("--saida", type=Path, default=Path("conservacao_msa"))
    ap.add_argument("--normalizacao", default="media", choices=["media", "soma"],
                    help="como fechar o peso de Henikoff (ver docstring)")
    ap.add_argument("--gap-como-tipo", action="store_true",
                    help="contar '-' como 21º tipo de resíduo nos pesos")
    ap.add_argument("--pesos-so-colunas-query", action="store_true",
                    help="pesar usando só as colunas em que a query tem resíduo")
    ap.add_argument("--janela", type=int, default=3, help="meia-janela da suavização da JSD")
    ap.add_argument("--peso-janela", type=float, default=0.5)
    ap.add_argument("--min-ocupacao", type=float, default=0.5,
                    help="ocupação ponderada mínima para a coluna ser 'confiavel'")
    ap.add_argument("--min-neff-coluna", type=float, default=3.0)
    ap.add_argument("--max-seqs", type=int, default=None, help="amostrar MSAs gigantes")
    ap.add_argument("--semente", type=int, default=0)
    ap.add_argument("--permitir-sequencia-unica", action="store_true",
                    help="não falhar quando o 'MSA' tem 1 sequência (ESMFold)")
    ap.add_argument("--coluna-seq", default=None, help="nome da coluna de sequência (csv/pqt)")
    ap.add_argument("--coluna-id", default=None)
    args = ap.parse_args(argv)

    try:
        r = analisar_msa(
            args.msa, rotulo=args.rotulo, formato=args.formato,
            gap_como_tipo=args.gap_como_tipo, normalizacao=args.normalizacao,
            so_colunas_query=args.pesos_so_colunas_query,
            meia_janela=args.janela, peso_janela=args.peso_janela,
            min_ocupacao=args.min_ocupacao, min_neff_coluna=args.min_neff_coluna,
            max_seqs=args.max_seqs, semente=args.semente,
            permitir_unica=args.permitir_sequencia_unica,
            coluna_id=args.coluna_id, coluna_seq=args.coluna_seq,
        )
    except ValueError as erro:
        print(f"ERRO: {erro}", file=sys.stderr)
        return 2

    saida = args.saida
    saida.mkdir(parents=True, exist_ok=True)
    rot = r["resumo"]["rotulo"]
    escrever_csv(saida / f"conservacao_colunas__{rot}.csv", r["colunas"])
    escrever_pesos(saida / f"pesos_henikoff__{rot}.csv", r["msa"], r["pesos"])
    (saida / f"resumo__{rot}.json").write_text(
        json.dumps(r["resumo"], indent=2, ensure_ascii=False)
    )

    s = r["resumo"]
    print(f"{s['rotulo']}: {s['n_sequencias']} sequências x {s['largura_alinhamento']} colunas "
          f"(formato {s['formato']})")
    print(f"  Neff (Henikoff/Kish) = {s['neff_henikoff']:.1f}  "
          f"Neff/resíduo = {s['neff_por_residuo']:.2f}  "
          f"redundância N/Neff = {s['redundancia_n_sobre_neff']:.1f}x")
    print(f"  colunas da query confiáveis: {s['colunas_query_confiaveis']}/"
          f"{s['comprimento_query']}")
    print(f"  {s['veredito']}")
    for a in s["avisos"]:
        print(f"  aviso: {a}")
    print(f"  -> {saida}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
