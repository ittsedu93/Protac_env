#!/usr/bin/env python
"""
Conservação (Henikoff) dos MSAs de cofold, cruzada com o bolso a 5 Å do ligante.

Só biblioteca padrão. Segundos.

    # 1) ver se tudo funciona, sem ter dado nenhum (gera dado sintético)
    python conservacao_bolso.py --demo --saida /tmp/demo_conservacao

    # 2) escrever o arquivo de configuração em branco, para preencher depois
    python conservacao_bolso.py --escrever-config config_conservacao.json

    # 3) rodar de verdade
    python conservacao_bolso.py --config config_conservacao.json --saida resultados/

    # ou direto na linha de comando
    python conservacao_bolso.py \
        --estrutura modelo_boltz2.cif --ligante-resname LIG --raio 5.0 \
        --msa boltz2:A:msa_alvo.fasta --msa boltz2:B:msa_e3.fasta \
        --saida resultados/

O que este script responde
--------------------------
1. O MSA que alimentou o cofold é confiável? (Neff de Henikoff, Neff/resíduo,
   redundância N/Neff, ocupação por coluna — ver msa_henikoff.py)
2. O bolso do ligante cai em região conservada ou em região variável? Com
   teste estatístico (Mann-Whitney + permutação) contra o resto da MESMA
   cadeia, não contra uma expectativa inventada.
3. O que exatamente está a menos de 5 Å do ligante — incluindo outras cadeias,
   cofatores, íons, glicanos, aditivos de cristalização e águas, cada um
   etiquetado (ver bolso_ligante.py).

Por que 2 e 3 andam juntas
--------------------------
Um escore de conservação alto no bolso é a evidência barata de que o cofold
colocou o ligante onde a evolução guarda sítio funcional. Um escore baixo tem
duas leituras possíveis, e elas se distinguem olhando o item 3: se o bolso é
de INTERFACE (duas cadeias poliméricas a 5 Å), conservação baixa é esperada —
interfaces induzidas de PROTAC são neo-interfaces, não sítios ancestrais. Se o
bolso é de uma cadeia só e a conservação é baixa, o modelo provavelmente
colocou o ligante num bolso de superfície qualquer.

Mapeamento estrutura <-> MSA
----------------------------
A query do MSA (1ª sequência) é alinhada à sequência da CADEIA da estrutura
com Needleman-Wunsch de pontas livres. Isso cobre o caso normal (modelo com
loops faltando, numeração deslocada) e, principalmente, PEGA O ERRO de passar
o MSA da cadeia errada: a identidade do alinhamento entra no relatório e
abaixo de 90% vira aviso.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import sys
from pathlib import Path

_AQUI = Path(__file__).resolve().parent
if str(_AQUI) not in sys.path:
    sys.path.insert(0, str(_AQUI))

import bolso_ligante as bl                                   # noqa: E402
import msa_henikoff as mh                                    # noqa: E402

# ---------------------------------------------------------------------------
# alinhamento query do MSA <-> sequência da cadeia


MAX_CELULAS_ALINHAMENTO = 25_000_000      # ~5000 x 5000; acima disso, 75 MB


def alinhar(a: str, b: str, acerto=2.0, erro=-1.0, abrir=-4.0, estender=-0.3,
            pontas_livres=True):
    """Needleman-Wunsch com gap AFIM (Gotoh) e pontas livres.

    Devolve (pares, identidade), em que `pares` são índices (em `a`, em `b`) e
    identidade é a fração de resíduos do menor dos dois que casaram.

    Pontas livres porque o modelo quase sempre é um pedaço da query (ou a query
    é um pedaço do construto): penalizar as pontas empurraria o alinhamento
    para fora de registro só para economizar gap terminal.

    O gap tem de ser AFIM, e isso não é refinamento: com gap linear de -1,5,
    pular um loop de 10 resíduos que falta no modelo custa -15, enquanto
    emendar a sequência fora de registro custa ~-4 em desencontros. O
    alinhamento linear escolhe o registro errado e todo o mapeamento
    estrutura->coluna sai deslocado, SEM dar erro. Com abertura -4 e extensão
    -0,3 o buraco custa -6,7 e o registro certo ganha.
    """
    n, m = len(a), len(b)
    if n == 0 or m == 0:
        return [], 0.0
    # atalho honesto: contido exatamente, que é o caso mais comum
    pos = a.find(b)
    if pos >= 0:
        return [(pos + k, k) for k in range(m)], 1.0
    if (n + 1) * (m + 1) > MAX_CELULAS_ALINHAMENTO:
        raise ValueError(
            f"alinhamento {n}x{m} é grande demais para a memória deste script. "
            "Recorte a query do MSA para a região da cadeia, ou rode com "
            "--sem-estrutura."
        )

    MENOS = float("-inf")
    # M: par alinhado; X: gap em b (consome a); Y: gap em a (consome b)
    m_ant = [MENOS] * (m + 1)
    x_ant = [MENOS] * (m + 1)
    y_ant = [MENOS] * (m + 1)
    m_ant[0] = 0.0
    for j in range(1, m + 1):
        y_ant[j] = 0.0 if pontas_livres else abrir + (j - 1) * estender
    # traceback: 1 = veio de M, 2 = de X, 3 = de Y
    tm = [bytearray(m + 1) for _ in range(n + 1)]
    tx = [bytearray(m + 1) for _ in range(n + 1)]
    ty = [bytearray(m + 1) for _ in range(n + 1)]
    melhor_i, melhor_j, melhor_v = 0, 0, MENOS

    for i in range(1, n + 1):
        mc = [MENOS] * (m + 1)
        xc = [MENOS] * (m + 1)
        yc = [MENOS] * (m + 1)
        xc[0] = 0.0 if pontas_livres else abrir + (i - 1) * estender
        tx[i][0] = 2 if i > 1 else 1
        ai = a[i - 1]
        for j in range(1, m + 1):
            s = acerto if ai == b[j - 1] else erro
            # M
            melhor, t = m_ant[j - 1], 1
            if x_ant[j - 1] > melhor:
                melhor, t = x_ant[j - 1], 2
            if y_ant[j - 1] > melhor:
                melhor, t = y_ant[j - 1], 3
            mc[j] = melhor + s if melhor > MENOS else MENOS
            tm[i][j] = t
            # X (gap em b): abre a partir de M, ou estende X
            ab = m_ant[j] + abrir
            ex = x_ant[j] + estender
            if ex > ab:
                xc[j], tx[i][j] = ex, 2
            else:
                xc[j], tx[i][j] = ab, 1
            # Y (gap em a): abre a partir de M, ou estende Y
            ab = mc[j - 1] + abrir
            ey = yc[j - 1] + estender
            if ey > ab:
                yc[j], ty[i][j] = ey, 3
            else:
                yc[j], ty[i][j] = ab, 1
        if pontas_livres and mc[m] > melhor_v:      # sobra de `a` pendurada
            melhor_i, melhor_j, melhor_v = i, m, mc[m]
        m_ant, x_ant, y_ant = mc, xc, yc

    if pontas_livres:
        for j in range(1, m + 1):                   # sobra de `b` pendurada
            if m_ant[j] > melhor_v:
                melhor_i, melhor_j, melhor_v = n, j, m_ant[j]
        i, j = melhor_i, melhor_j
    else:
        i, j = n, m
    if melhor_v == MENOS and pontas_livres:
        return [], 0.0

    pares = []
    estado = 1                                      # termina sempre num par
    while i > 0 and j > 0:
        if estado == 1:
            pares.append((i - 1, j - 1))
            prox = tm[i][j]
            i, j, estado = i - 1, j - 1, prox
        elif estado == 2:
            prox = tx[i][j]
            i, estado = i - 1, prox
        else:
            prox = ty[i][j]
            j, estado = j - 1, prox
    pares.reverse()
    if not pares:
        return [], 0.0
    iguais = sum(1 for ia, ib in pares if a[ia] == b[ib])
    # a identidade é sobre o MENOR dos dois, não sobre os pares alinhados:
    # senão um alinhamento de 2 resíduos casados marcaria 100% e o aviso de
    # "MSA da cadeia errada" nunca dispararia
    return pares, iguais / min(n, m)


def mapa_query_para_coluna(query_alinhada: str) -> list[int]:
    """índice (0-based) na query sem gap -> índice da coluna no alinhamento."""
    return [j for j, c in enumerate(query_alinhada) if c != mh.GAP]


# ---------------------------------------------------------------------------
# estatística sem scipy


def _ranks(valores):
    ordem = sorted(range(len(valores)), key=lambda k: valores[k])
    ranks = [0.0] * len(valores)
    i = 0
    empates = []
    while i < len(ordem):
        j = i
        while j + 1 < len(ordem) and valores[ordem[j + 1]] == valores[ordem[i]]:
            j += 1
        r = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[ordem[k]] = r
        if j > i:
            empates.append(j - i + 1)
        i = j + 1
    return ranks, empates


def mann_whitney(x, y):
    """U de Mann-Whitney, bilateral, aproximação normal com correção de empates.

    Com n pequeno (bolso de 8 resíduos) a aproximação normal é grosseira — por
    isso o relatório sempre traz TAMBÉM o p de permutação, que é exato por
    construção. Os dois juntos evitam tanto o falso positivo quanto a
    sensação de rigor que um único p daria.
    """
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return {"U": float("nan"), "z": float("nan"), "p": float("nan"),
                "cliff_delta": float("nan"), "n1": n1, "n2": n2}
    todos = list(x) + list(y)
    ranks, empates = _ranks(todos)
    r1 = sum(ranks[:n1])
    u1 = r1 - n1 * (n1 + 1) / 2.0
    n = n1 + n2
    mu = n1 * n2 / 2.0
    corr = sum(t ** 3 - t for t in empates)
    var = (n1 * n2 / 12.0) * ((n + 1) - corr / (n * (n - 1))) if n > 1 else 0.0
    z = (u1 - mu) / math.sqrt(var) if var > 0 else 0.0
    p = math.erfc(abs(z) / math.sqrt(2.0))
    return {"U": u1, "z": z, "p": p,
            "cliff_delta": 2.0 * u1 / (n1 * n2) - 1.0, "n1": n1, "n2": n2}


def teste_permutacao(x, y, n_perm=10000, semente=0, estatistica="mediana"):
    """p bilateral por permutação da etiqueta bolso/resto. Determinístico."""
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return {"p": float("nan"), "n_perm": 0, "observado": float("nan")}
    f = (lambda v: _mediana(v)) if estatistica == "mediana" else (
        lambda v: sum(v) / len(v))
    todos = list(x) + list(y)
    obs = f(list(x)) - f(list(y))
    rng = random.Random(semente)
    extremos = 0
    for _ in range(n_perm):
        rng.shuffle(todos)
        d = f(todos[:n1]) - f(todos[n1:])
        if abs(d) >= abs(obs) - 1e-12:
            extremos += 1
    return {"p": (extremos + 1) / (n_perm + 1), "n_perm": n_perm, "observado": obs}


def _mediana(v):
    s = sorted(v)
    n = len(s)
    if n == 0:
        return float("nan")
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def spearman(x, y):
    """ρ de Spearman (Pearson sobre ranks). Para comparar ferramentas."""
    if len(x) != len(y) or len(x) < 3:
        return float("nan")
    rx, _ = _ranks(list(x))
    ry, _ = _ranks(list(y))
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx > 0 and dy > 0 else float("nan")


# ---------------------------------------------------------------------------
# cruzamento


def residuos_da_cadeia(
    analise_msa: dict,
    atomos,
    cadeia: str,
    contatos: list[dict],
):
    """Uma linha por resíduo da cadeia: conservação + se está no bolso.

    Devolve (linhas, diagnostico). O diagnóstico é o que impede a leitura
    errada: identidade do alinhamento query<->estrutura, resíduos sem coluna
    e discordâncias de identidade de resíduo.
    """
    msa = analise_msa["msa"]
    colunas = analise_msa["colunas"]
    seq_est, chaves = bl.sequencia_cadeia(atomos, cadeia)
    if not seq_est:
        raise ValueError(
            f"cadeia {cadeia!r} não tem resíduo de aminoácido na estrutura "
            f"(cadeias com proteína: {', '.join(bl.cadeias_polimericas(atomos)) or 'nenhuma'})"
        )
    query = msa.query_sem_gap
    pares, identidade = alinhar(query, seq_est)
    col_de_query = mapa_query_para_coluna(msa.query)
    mapa = {ib: ia for ia, ib in pares}              # índice estrutura -> query

    no_bolso = {}
    for c in contatos:
        if c["cadeia"] == cadeia and c["categoria"] == "proteina":
            no_bolso[(c["resnum"], (c["icode"] or " ").ljust(1))] = c

    linhas, discordancias, sem_coluna = [], 0, 0
    for k, (resnum, icode, resname) in enumerate(chaves):
        iq = mapa.get(k)
        col = colunas[col_de_query[iq]] if iq is not None else None
        chave_bolso = (resnum, icode if icode.strip() else " ")
        cont = no_bolso.get(chave_bolso) or no_bolso.get((resnum, " "))
        aa = seq_est[k]
        if col is not None and col["res_query"] not in (aa, mh.GAP):
            discordancias += 1
        if col is None:
            sem_coluna += 1
        linhas.append({
            "cadeia": cadeia,
            "resnum": resnum,
            "icode": icode.strip(),
            "resname": resname,
            "aa_estrutura": aa,
            "pos_query": (iq + 1) if iq is not None else "",
            "coluna_msa": col["coluna"] if col else "",
            "res_query_msa": col["res_query"] if col else "",
            "profundidade": col["profundidade"] if col else "",
            "ocupacao_pond": col["ocupacao_pond"] if col else "",
            "neff_coluna": col["neff_coluna"] if col else "",
            "conservacao_shannon": col["conservacao_shannon"] if col else "",
            "jsd": col["jsd"] if col else "",
            "jsd_janela": col["jsd_janela"] if col else "",
            "jsd_penalizado": col["jsd_penalizado"] if col else "",
            "jsd_sem_peso": col.get("jsd_sem_peso", "") if col else "",
            "coluna_confiavel": col["confiavel"] if col else "",
            "no_bolso": 1 if cont else 0,
            "dist_min_ligante": round(cont["distancia_min"], 3) if cont else "",
            "atomos_em_contato": cont["n_atomos_em_contato"] if cont else "",
        })
    diagnostico = {
        "cadeia": cadeia,
        "residuos_na_estrutura": len(chaves),
        "comprimento_query_msa": len(query),
        "pares_alinhados": len(pares),
        "identidade_query_estrutura": identidade,
        "residuos_sem_coluna_no_msa": sem_coluna,
        "discordancias_de_residuo": discordancias,
        "residuos_no_bolso": sum(1 for l in linhas if l["no_bolso"]),
        "avisos": [],
    }
    if identidade < 0.9:
        diagnostico["avisos"].append(
            f"identidade query do MSA x cadeia {cadeia} = {identidade:.1%}: "
            "provavelmente este MSA não é desta cadeia (ou é um ortólogo "
            "distante). Confira o par ferramenta:cadeia:arquivo."
        )
    if sem_coluna > 0.2 * len(chaves):
        diagnostico["avisos"].append(
            f"{sem_coluna}/{len(chaves)} resíduos da cadeia não têm coluna no MSA"
        )
    if discordancias:
        diagnostico["avisos"].append(
            f"{discordancias} resíduos com identidade diferente da query do MSA "
            "(mutação, construto diferente, ou desalinhamento)"
        )
    if not diagnostico["residuos_no_bolso"]:
        diagnostico["avisos"].append(
            f"nenhum resíduo da cadeia {cadeia} está no raio do ligante — "
            "a comparação bolso x resto não se aplica a esta cadeia"
        )
    return linhas, diagnostico


def estatistica_do_bolso(linhas: list[dict], escore: str = "jsd_janela",
                         so_confiaveis: bool = True, semente: int = 0,
                         n_perm: int = 10000):
    """Bolso x resto da cadeia, no escore pedido."""
    def val(l):
        v = l.get(escore, "")
        return v if isinstance(v, (int, float)) and v != "" else None

    bolso, resto = [], []
    for l in linhas:
        v = val(l)
        if v is None:
            continue
        if so_confiaveis and l.get("coluna_confiavel") in (0, "", None):
            continue
        (bolso if l["no_bolso"] else resto).append(v)

    pesos = [
        (val(l), l["atomos_em_contato"] or 1) for l in linhas
        if l["no_bolso"] and val(l) is not None
    ]
    media_pond = (
        sum(v * w for v, w in pesos) / sum(w for _, w in pesos) if pesos else float("nan")
    )
    mw = mann_whitney(bolso, resto)
    perm = teste_permutacao(bolso, resto, n_perm=n_perm, semente=semente)
    return {
        "escore": escore,
        "so_colunas_confiaveis": so_confiaveis,
        "n_bolso": len(bolso),
        "n_resto": len(resto),
        "media_bolso": sum(bolso) / len(bolso) if bolso else float("nan"),
        "mediana_bolso": _mediana(bolso),
        "media_resto": sum(resto) / len(resto) if resto else float("nan"),
        "mediana_resto": _mediana(resto),
        "media_bolso_ponderada_por_contato": media_pond,
        "mann_whitney": mw,
        "permutacao": perm,
        "veredito": _veredito_bolso(bolso, resto, mw, perm),
    }


def _veredito_bolso(bolso, resto, mw, perm):
    if len(bolso) < 3 or len(resto) < 3:
        return ("amostra pequena demais para teste "
                f"(bolso n={len(bolso)}, resto n={len(resto)})")
    d = mw["cliff_delta"]
    p = perm["p"]
    direcao = "MAIS conservado" if d > 0 else "MENOS conservado"
    if p > 0.05:
        return (f"bolso indistinguível do resto da cadeia (p={p:.3f}, "
                f"Cliff δ={d:+.2f})")
    força = "forte" if abs(d) >= 0.5 else ("moderado" if abs(d) >= 0.3 else "fraco")
    return (f"bolso {direcao} que o resto da cadeia: efeito {força} "
            f"(Cliff δ={d:+.2f}, p_perm={p:.4f})")


# ---------------------------------------------------------------------------
# configuração


MODELO_CONFIG = {
    "_comentario": [
        "Preencha os campos. 'null' = não informado.",
        "estrutura: o modelo de cofold (PDB ou mmCIF). Boltz-2 e Chai-1/2 saem",
        "  em .cif; ESMFold e OpenFold saem em .pdb.",
        "ligante: como achar o ligante no arquivo. resname 'LIG' cobre Boltz e",
        "  Chai; use 'auto': true para pegar o maior hetero-resíduo.",
        "msas: um por (ferramenta, cadeia). FASTA ALINHADO é o formato",
        "  recomendado; .a3m, .sto, .csv e .pqt também são lidos.",
        "  'cadeia' é a cadeia da ESTRUTURA que aquele MSA descreve.",
    ],
    "estrutura": "/caminho/para/modelo.cif",
    "modelo_nmr": 1,
    "ligante": {
        "resname": "LIG",
        "cadeia": None,
        "resnum": None,
        "copia": None,
        "auto": False,
    },
    "raio_angstrom": 5.0,
    "com_hidrogenios": False,
    "incluir_agua": True,
    "msas": [
        {"ferramenta": "boltz2", "cadeia": "A",
         "arquivo": "/caminho/para/msa_cadeia_A.fasta", "formato": "auto"},
        {"ferramenta": "boltz2", "cadeia": "B",
         "arquivo": "/caminho/para/msa_cadeia_B.fasta", "formato": "auto"},
        {"ferramenta": "chai1", "cadeia": "A",
         "arquivo": "/caminho/para/chai_cadeia_A.aligned.pqt", "formato": "auto"},
        {"ferramenta": "openfold", "cadeia": "A",
         "arquivo": "/caminho/para/uniref90_hits.sto", "formato": "auto"},
    ],
    "pesos": {
        "normalizacao": "media",
        "gap_como_tipo": False,
        "so_colunas_query": False,
        "max_seqs": None,
        "semente": 0,
    },
    "escores": {
        "escore_principal": "jsd_janela",
        "janela": 3,
        "peso_janela": 0.5,
        "min_ocupacao": 0.5,
        "min_neff_coluna": 3.0,
        "so_colunas_confiaveis": True,
        "n_permutacoes": 10000,
        "comparar_sem_peso": True,
    },
    "permitir_sequencia_unica": False,
}

ESCORES_VALIDOS = ("jsd", "jsd_janela", "jsd_penalizado", "conservacao_shannon",
                   "kl_bits", "jsd_sem_peso")


def mesclar_config(padrao, dado):
    """Config do usuário sobre o modelo, sem exigir que ele repita tudo."""
    if not isinstance(dado, dict):
        return dado
    out = dict(padrao)
    for k, v in dado.items():
        out[k] = mesclar_config(padrao[k], v) if isinstance(padrao.get(k), dict) else v
    return out


def ler_config(caminho) -> dict:
    caminho = Path(caminho)
    if not caminho.exists():
        raise ValueError(f"config não encontrada: {caminho}")
    try:
        dado = json.loads(caminho.read_text(encoding="utf-8"))
    except json.JSONDecodeError as erro:
        raise ValueError(f"{caminho}: JSON inválido — {erro}") from erro
    if not isinstance(dado, dict):
        raise ValueError(f"{caminho}: o JSON deve ser um objeto {{...}}")
    return mesclar_config(MODELO_CONFIG, dado)


def validar_config(cfg: dict, exigir_estrutura: bool = True) -> list[str]:
    """Todos os problemas de uma vez — não um por execução."""
    erros = []
    placeholder = ("/caminho/para", "CAMINHO", "<", "preencha", "PREENCHA")

    def em_branco(v):
        return v is None or (isinstance(v, str) and
                             (not v.strip() or any(p in v for p in placeholder)))

    if exigir_estrutura:
        est = cfg.get("estrutura")
        if em_branco(est):
            erros.append("campo 'estrutura' não preenchido")
        elif not Path(est).exists():
            erros.append(f"estrutura não encontrada: {est}")

    msas = cfg.get("msas") or []
    if not isinstance(msas, list) or not msas:
        erros.append("campo 'msas' vazio: informe pelo menos um MSA")
    else:
        vistos = set()
        for i, m in enumerate(msas, 1):
            if not isinstance(m, dict):
                erros.append(f"msas[{i}] não é um objeto")
                continue
            arq = m.get("arquivo")
            if em_branco(arq):
                erros.append(f"msas[{i}]: 'arquivo' não preenchido")
            elif not Path(arq).exists():
                erros.append(f"msas[{i}]: MSA não encontrado: {arq}")
            if exigir_estrutura and em_branco(m.get("cadeia")):
                erros.append(f"msas[{i}]: 'cadeia' não preenchida "
                             "(qual cadeia da estrutura este MSA descreve?)")
            chave = (m.get("ferramenta"), m.get("cadeia"))
            if chave in vistos:
                erros.append(f"msas[{i}]: par ferramenta/cadeia repetido: {chave}")
            vistos.add(chave)

    raio = cfg.get("raio_angstrom", 5.0)
    if not isinstance(raio, (int, float)) or raio <= 0 or raio > 50:
        erros.append(f"raio_angstrom inválido: {raio!r} (use algo entre 0 e 50)")

    esc_cfg = cfg.get("escores") or {}
    esc = esc_cfg.get("escore_principal", "jsd_janela")
    if esc not in ESCORES_VALIDOS:
        erros.append(f"escore_principal {esc!r} não existe. Use um de: "
                     f"{', '.join(ESCORES_VALIDOS)}")
    elif esc == "jsd_sem_peso" and not esc_cfg.get("comparar_sem_peso", True):
        erros.append(
            "escore_principal 'jsd_sem_peso' exige escores.comparar_sem_peso = true "
            "(é a passagem que calcula a conservação sem os pesos)"
        )
    norm = (cfg.get("pesos") or {}).get("normalizacao", "media")
    if norm not in ("media", "soma"):
        erros.append(f"pesos.normalizacao {norm!r}: use 'media' ou 'soma'")

    lig = cfg.get("ligante") or {}
    if exigir_estrutura and not lig.get("auto") and em_branco(lig.get("resname")) \
            and em_branco(lig.get("cadeia")) and lig.get("resnum") is None:
        erros.append("ligante: informe 'resname' (e/ou cadeia/resnum) ou use "
                     "\"auto\": true")
    return erros


def config_de_argumentos(args) -> dict:
    cfg = json.loads(json.dumps(MODELO_CONFIG))
    cfg.pop("_comentario", None)
    cfg["estrutura"] = str(args.estrutura) if args.estrutura else None
    cfg["modelo_nmr"] = args.modelo
    cfg["raio_angstrom"] = args.raio
    cfg["com_hidrogenios"] = args.com_hidrogenios
    cfg["incluir_agua"] = not args.sem_agua
    cfg["ligante"] = {
        "resname": args.ligante_resname, "cadeia": args.ligante_cadeia,
        "resnum": args.ligante_resnum, "copia": args.ligante_copia,
        "auto": args.ligante_auto,
    }
    cfg["msas"] = [_parse_msa_arg(s) for s in (args.msa or [])]
    cfg["pesos"].update({
        "normalizacao": args.normalizacao, "gap_como_tipo": args.gap_como_tipo,
        "so_colunas_query": args.pesos_so_colunas_query,
        "max_seqs": args.max_seqs, "semente": args.semente,
    })
    cfg["escores"].update({
        "escore_principal": args.escore, "janela": args.janela,
        "min_ocupacao": args.min_ocupacao, "min_neff_coluna": args.min_neff_coluna,
        "so_colunas_confiaveis": not args.todas_as_colunas,
        "n_permutacoes": args.n_permutacoes,
    })
    cfg["permitir_sequencia_unica"] = args.permitir_sequencia_unica
    return cfg


def _parse_msa_arg(s: str) -> dict:
    """'ferramenta:cadeia:caminho' (ou 'cadeia:caminho')."""
    partes = s.split(":", 2)
    if len(partes) == 3:
        ferr, cad, arq = partes
    elif len(partes) == 2:
        cad, arq = partes
        ferr = Path(arq).stem
    else:
        raise ValueError(
            f"--msa {s!r} malformado. Use ferramenta:cadeia:caminho, "
            "ex.: boltz2:A:/dados/msa_A.fasta"
        )
    return {"ferramenta": ferr.strip(), "cadeia": cad.strip(),
            "arquivo": arq.strip(), "formato": "auto"}


# ---------------------------------------------------------------------------
# execução


def executar(cfg: dict, saida: Path, sem_estrutura: bool = False) -> dict:
    """O pipeline inteiro. Devolve o dict que vai para resumo.json."""
    saida = Path(saida)
    saida.mkdir(parents=True, exist_ok=True)
    erros = validar_config(cfg, exigir_estrutura=not sem_estrutura)
    if erros:
        raise ValueError("configuração incompleta:\n  - " + "\n  - ".join(erros))

    pesos_cfg = cfg["pesos"]
    esc_cfg = cfg["escores"]
    escore = esc_cfg["escore_principal"]
    resultado = {
        "configuracao": {k: v for k, v in cfg.items() if k != "_comentario"},
        "msas": [], "cadeias": [], "avisos": [],
    }

    # --- bolso ------------------------------------------------------------
    atomos, contatos, resumo_bolso = [], [], None
    if not sem_estrutura:
        lig = cfg["ligante"]
        r = bl.analisar_bolso(
            cfg["estrutura"], resname=lig.get("resname"), cadeia=lig.get("cadeia"),
            resnum=lig.get("resnum"), auto=bool(lig.get("auto")),
            copia=lig.get("copia"), raio=cfg["raio_angstrom"],
            modelo=cfg.get("modelo_nmr", 1),
            com_hidrogenios=bool(cfg.get("com_hidrogenios")),
            incluir_agua=bool(cfg.get("incluir_agua", True)),
        )
        atomos, contatos = r["atomos"], r["contatos"]
        resumo_bolso = r["resumo"]
        bl.escrever_csv(saida / "contatos_5A_residuos.csv", r["contatos"])
        bl.escrever_csv(saida / "contatos_5A_atomos.csv", r["pares"])
        (saida / "bolso.pml").write_text(
            bl.selecao_pymol(r["contatos"], r["descricao"], cfg["raio_angstrom"])
        )
        resultado["bolso"] = resumo_bolso
        if not resumo_bolso["cadeias_polimericas_em_contato"]:
            resultado["avisos"].append(
                "o ligante não toca nenhuma cadeia polimérica no raio pedido — "
                "confira se o ligante selecionado é o certo"
            )

    # --- MSAs -------------------------------------------------------------
    analises = {}
    for m in cfg["msas"]:
        rotulo = f"{m.get('ferramenta', 'msa')}__{m.get('cadeia', '?')}"
        try:
            a = mh.analisar_msa(
                m["arquivo"], rotulo=rotulo, formato=m.get("formato", "auto"),
                gap_como_tipo=pesos_cfg["gap_como_tipo"],
                normalizacao=pesos_cfg["normalizacao"],
                so_colunas_query=pesos_cfg["so_colunas_query"],
                meia_janela=esc_cfg["janela"], peso_janela=esc_cfg["peso_janela"],
                min_ocupacao=esc_cfg["min_ocupacao"],
                min_neff_coluna=esc_cfg["min_neff_coluna"],
                max_seqs=pesos_cfg["max_seqs"], semente=pesos_cfg["semente"],
                permitir_unica=bool(cfg.get("permitir_sequencia_unica")),
                comparar_sem_peso=bool(esc_cfg.get("comparar_sem_peso", True)),
                coluna_id=m.get("coluna_id"), coluna_seq=m.get("coluna_seq"),
            )
        except ValueError as erro:
            resultado["msas"].append({
                "rotulo": rotulo, "ferramenta": m.get("ferramenta"),
                "cadeia": m.get("cadeia"), "arquivo": m.get("arquivo"),
                "erro": str(erro),
            })
            resultado["avisos"].append(f"{rotulo}: {erro}")
            continue
        analises[(m.get("ferramenta"), m.get("cadeia"))] = (m, a)
        mh.escrever_csv(saida / f"conservacao_colunas__{rotulo}.csv", a["colunas"])
        mh.escrever_pesos(saida / f"pesos_henikoff__{rotulo}.csv", a["msa"], a["pesos"])
        info = dict(a["resumo"])
        info.update({"ferramenta": m.get("ferramenta"), "cadeia": m.get("cadeia")})
        resultado["msas"].append(info)

    if not analises:
        raise ValueError(
            "nenhum MSA pôde ser analisado:\n  - " + "\n  - ".join(resultado["avisos"])
        )

    # --- cruzamento com o bolso -------------------------------------------
    if not sem_estrutura:
        for (ferr, cad), (m, a) in analises.items():
            try:
                linhas, diag = residuos_da_cadeia(a, atomos, cad, contatos)
            except ValueError as erro:
                resultado["avisos"].append(f"{ferr}/{cad}: {erro}")
                continue
            rotulo = f"{ferr}__{cad}"
            mh.escrever_csv(saida / f"residuos__{rotulo}.csv", linhas)
            est = estatistica_do_bolso(
                linhas, escore=escore,
                so_confiaveis=bool(esc_cfg["so_colunas_confiaveis"]),
                semente=pesos_cfg["semente"], n_perm=esc_cfg["n_permutacoes"],
            )
            top = sorted(
                (l for l in linhas if l["no_bolso"] and isinstance(l.get(escore), float)),
                key=lambda l: l[escore], reverse=True,
            )
            est["residuos_do_bolso"] = [
                {"residuo": f"{l['cadeia']}/{l['resname']}{l['resnum']}",
                 "dist_min": l["dist_min_ligante"], escore: round(l[escore], 4),
                 "profundidade": l["profundidade"],
                 "coluna_confiavel": l["coluna_confiavel"]}
                for l in top
            ]
            resultado["cadeias"].append(
                {"ferramenta": ferr, "cadeia": cad, "diagnostico": diag,
                 "estatistica": est}
            )
            resultado["avisos"].extend(f"{ferr}/{cad}: {x}" for x in diag["avisos"])

        # --- comparação entre ferramentas na MESMA cadeia -----------------
        comp = []
        por_cadeia = {}
        for (ferr, cad), (m, a) in analises.items():
            por_cadeia.setdefault(cad, []).append((ferr, a))
        for cad, itens in sorted(por_cadeia.items()):
            if len(itens) < 2:
                continue
            for i in range(len(itens)):
                for j in range(i + 1, len(itens)):
                    f1, a1 = itens[i]
                    f2, a2 = itens[j]
                    s1 = {l["pos_query"]: l for l in a1["colunas"] if l["pos_query"] != ""}
                    s2 = {l["pos_query"]: l for l in a2["colunas"] if l["pos_query"] != ""}
                    comuns = sorted(set(s1) & set(s2))
                    x = [s1[p][escore] for p in comuns]
                    y = [s2[p][escore] for p in comuns]
                    comp.append({
                        "cadeia": cad, "ferramenta_a": f1, "ferramenta_b": f2,
                        "posicoes_comparadas": len(comuns),
                        "neff_a": round(a1["resumo"]["neff_henikoff"], 2),
                        "neff_b": round(a2["resumo"]["neff_henikoff"], 2),
                        "media_a": round(sum(x) / len(x), 4) if x else "",
                        "media_b": round(sum(y) / len(y), 4) if y else "",
                        "spearman": round(spearman(x, y), 4) if len(comuns) >= 3 else "",
                    })
        if comp:
            mh.escrever_csv(saida / "comparacao_ferramentas.csv", comp)
            resultado["comparacao_ferramentas"] = comp

    (saida / "resumo.json").write_text(
        json.dumps(resultado, indent=2, ensure_ascii=False, default=str)
    )
    (saida / "relatorio.md").write_text(relatorio_markdown(resultado))
    return resultado


# ---------------------------------------------------------------------------
# relatório


def relatorio_markdown(r: dict) -> str:
    L = []
    cfg = r.get("configuracao", {})
    esc = (cfg.get("escores") or {}).get("escore_principal", "jsd_janela")
    L.append("# Conservação do MSA (pesos de Henikoff) e bolso do ligante\n")
    L.append(f"- estrutura: `{cfg.get('estrutura')}`")
    L.append(f"- raio de contato: {cfg.get('raio_angstrom')} Å")
    L.append(f"- escore principal: `{esc}`")
    pz = cfg.get("pesos", {})
    L.append(f"- pesos: Henikoff, normalização `{pz.get('normalizacao')}`, "
             f"gap como tipo = {pz.get('gap_como_tipo')}\n")

    b = r.get("bolso")
    if b:
        d = b["ligante"]
        L.append("## 1. O que está a menos de "
                 f"{b['raio_angstrom']} Å do ligante\n")
        L.append(f"Ligante: **{d.get('resname')}** (cadeia `{d.get('cadeia')}`, "
                 f"resnum {d.get('resnum')}, {d.get('n_atomos')} átomos, "
                 f"seleção por {d.get('como')}).\n")
        L.append(f"- itens no raio: **{b['n_itens_no_raio']}** "
                 f"({b['n_residuos_proteicos']} resíduos proteicos)")
        L.append(f"- cadeias poliméricas em contato: "
                 f"`{', '.join(b['cadeias_polimericas_em_contato']) or '—'}`"
                 + ("  → **bolso de interface**" if b["bolso_de_interface"] else ""))
        L.append(f"- distância mínima: {b['distancia_minima']:.2f} Å"
                 if b.get("distancia_minima") is not None else "")
        L.append("\n| categoria | itens |\n|---|---|")
        for k, v in b["por_categoria"].items():
            L.append(f"| {k} | {v} |")
        if b["nao_proteicos_no_raio"]:
            L.append("\n### Não-proteicos no raio "
                     "(cofatores, íons, glicanos, aditivos, águas)\n")
            L.append("| item | categoria | subtipo | dist. mín. (Å) | átomos |\n"
                     "|---|---|---|---|---|")
            for c in b["nao_proteicos_no_raio"]:
                L.append(f"| `{c['residuo']}` | {c['categoria']} | {c['subtipo']} "
                         f"| {c['distancia_min']:.2f} | {c['n_atomos_em_contato']} |")
        else:
            L.append("\nNenhum cofator, íon, glicano, aditivo ou água no raio.")
        L.append("\nListas completas: `contatos_5A_residuos.csv` (por resíduo), "
                 "`contatos_5A_atomos.csv` (par de átomos), `bolso.pml` (PyMOL).\n")

    L.append("## 2. Confiabilidade de cada MSA\n")
    L.append("| MSA | formato | N | Neff | Neff/res | N/Neff | colunas confiáveis | veredito |")
    L.append("|---|---|---|---|---|---|---|---|")
    for m in r["msas"]:
        if m.get("erro"):
            L.append(f"| {m['rotulo']} | — | — | — | — | — | — | **ERRO**: {m['erro']} |")
            continue
        L.append(
            f"| {m['rotulo']} | {m['formato']} | {m['n_sequencias']} | "
            f"{m['neff_henikoff']:.1f} | {m['neff_por_residuo']:.2f} | "
            f"{m['redundancia_n_sobre_neff']:.1f}x | "
            f"{m['colunas_query_confiaveis']}/{m['comprimento_query']} | "
            f"{m['veredito']} |"
        )
    for m in r["msas"]:
        if m.get("erro") or not m.get("avisos"):
            continue
        L.append(f"\n**{m['rotulo']}**")
        for a in m["avisos"]:
            L.append(f"- {a}")
    efeitos = [m for m in r["msas"] if m.get("efeito_dos_pesos_jsd") is not None
               and not m.get("erro")]
    if efeitos:
        L.append("\n### Efeito dos pesos de Henikoff\n")
        L.append("JSD média nas colunas da query, com e sem peso. A diferença é "
                 "a inflação (ou deflação) que a redundância do banco causava.\n")
        L.append("| MSA | com Henikoff | sem peso | diferença |\n|---|---|---|---|")
        for m in efeitos:
            L.append(f"| {m['rotulo']} | {m['jsd_media_colunas_query']:.4f} | "
                     f"{m['jsd_media_sem_peso']:.4f} | "
                     f"{m['efeito_dos_pesos_jsd']:+.4f} |")

    if r.get("cadeias"):
        L.append("\n## 3. O bolso é conservado?\n")
        L.append(f"Bolso x resto da mesma cadeia, escore `{esc}`. "
                 "O p de permutação é o que vale com n pequeno; "
                 "Cliff δ é o tamanho do efeito (−1 a +1).\n")
        L.append("| MSA/cadeia | n bolso | n resto | mediana bolso | mediana resto "
                 "| Cliff δ | p (MW) | p (perm) | veredito |")
        L.append("|---|---|---|---|---|---|---|---|---|")
        for c in r["cadeias"]:
            e = c["estatistica"]
            mw, pm = e["mann_whitney"], e["permutacao"]
            L.append(
                f"| {c['ferramenta']}/{c['cadeia']} | {e['n_bolso']} | {e['n_resto']} | "
                f"{e['mediana_bolso']:.4f} | {e['mediana_resto']:.4f} | "
                f"{mw['cliff_delta']:+.2f} | {mw['p']:.4f} | {pm['p']:.4f} | "
                f"{e['veredito']} |"
            )
        for c in r["cadeias"]:
            e = c["estatistica"]
            d = c["diagnostico"]
            L.append(f"\n### {c['ferramenta']} / cadeia {c['cadeia']}\n")
            L.append(f"- alinhamento query do MSA x estrutura: "
                     f"identidade {d['identidade_query_estrutura']:.1%}, "
                     f"{d['pares_alinhados']} pares, "
                     f"{d['residuos_sem_coluna_no_msa']} resíduos sem coluna")
            L.append(f"- resíduos do bolso nesta cadeia: {d['residuos_no_bolso']}")
            L.append(f"- média do bolso ponderada por átomos em contato: "
                     f"{e['media_bolso_ponderada_por_contato']:.4f}")
            for a in d["avisos"]:
                L.append(f"- ⚠ {a}")
            if e["residuos_do_bolso"]:
                L.append(f"\n| resíduo | dist. mín. (Å) | {esc} | profundidade "
                         "| coluna confiável |\n|---|---|---|---|---|")
                for x in e["residuos_do_bolso"]:
                    L.append(f"| `{x['residuo']}` | {x['dist_min']} | {x[esc]} | "
                             f"{x['profundidade']} | "
                             f"{'sim' if x['coluna_confiavel'] else 'não'} |")

    if r.get("comparacao_ferramentas"):
        L.append("\n## 4. As ferramentas concordam?\n")
        L.append("Mesma cadeia, MSAs de origens diferentes. ρ de Spearman alto com "
                 "Neff muito diferente significa: o perfil é o mesmo, mas uma das "
                 "buscas é mais rasa.\n")
        L.append("| cadeia | A | B | posições | Neff A | Neff B | média A | média B | ρ |")
        L.append("|---|---|---|---|---|---|---|---|---|")
        for c in r["comparacao_ferramentas"]:
            L.append(f"| {c['cadeia']} | {c['ferramenta_a']} | {c['ferramenta_b']} | "
                     f"{c['posicoes_comparadas']} | {c['neff_a']} | {c['neff_b']} | "
                     f"{c['media_a']} | {c['media_b']} | {c['spearman']} |")

    if r.get("avisos"):
        L.append("\n## Avisos\n")
        for a in dict.fromkeys(r["avisos"]):
            L.append(f"- {a}")

    L.append("\n---\n")
    L.append("## Como ler os números\n")
    L.append("- **Neff/resíduo < 1**: MSA raso. É a faixa em que o próprio cofold "
             "erra mais, e em que conservação por coluna não sustenta conclusão.")
    L.append("- **N/Neff grande** (ex.: 50x): o banco estava redundante; sem os "
             "pesos de Henikoff a conservação estaria inflada nessa proporção.")
    L.append("- **coluna confiável**: ocupação ponderada e Neff da coluna acima dos "
             "limiares. Coluna rasa pode marcar 100% de conservação por acidente.")
    L.append("- **bolso de interface**: o ligante toca ≥2 cadeias poliméricas. Aí "
             "conservação baixa no bolso é esperada — neo-interface de PROTAC não é "
             "sítio ancestral.")
    return "\n".join(x for x in L if x is not None) + "\n"


# ---------------------------------------------------------------------------
# demo sintética (serve de autoteste e de exemplo de formato)


def _linha_pdb(serial, nome, resname, cadeia, resnum, x, y, z, elemento, hetero=False):
    nome_fmt = f" {nome:<3s}" if len(nome) < 4 else nome[:4]
    reg = "HETATM" if hetero else "ATOM  "
    return (f"{reg}{serial:5d} {nome_fmt}{'':1s}{resname:>3s} {cadeia:1s}"
            f"{resnum:4d}{'':1s}   {x:8.3f}{y:8.3f}{z:8.3f}"
            f"{1.0:6.2f}{0.0:6.2f}          {elemento:>2s}")


def gerar_demo(destino) -> dict:
    """Escreve um caso sintético COMPLETO e devolve os caminhos.

    O caso é desenhado para ter resposta conhecida: a cadeia A tem 30 resíduos,
    os resíduos 11-13 estão a ~4 Å do ligante, a cadeia B toca o ligante com 2
    resíduos (bolso de interface), e há um Zn, uma água e um glicerol no raio.
    No MSA, as colunas do bolso são conservadas e o resto é variável — então o
    teste estatístico tem de dizer "bolso mais conservado".

    Também existe para servir de exemplo de FORMATO: o FASTA alinhado gerado
    aqui é exatamente o que o script espera receber.
    """
    destino = Path(destino)
    destino.mkdir(parents=True, exist_ok=True)
    rng = random.Random(20251013)

    n_a, n_b = 30, 20
    bolso_a = [10, 11, 12]                      # 0-based -> resnums 11,12,13
    aminos = mh.ALFABETO
    seq_a = "".join(rng.choice(aminos) for _ in range(n_a))
    seq_b = "".join(rng.choice(aminos) for _ in range(n_b))

    linhas, serial = [], 1
    for i, aa in enumerate(seq_a):
        resname = {v: k for k, v in list(bl.AA3_PARA_1.items())[:20]}.get(aa, "ALA")
        linhas.append(_linha_pdb(serial, "CA", resname, "A", i + 1,
                                 i * 6.0, 0.0, 0.0, "C"))
        serial += 1
    for j, aa in enumerate(seq_b):
        resname = {v: k for k, v in list(bl.AA3_PARA_1.items())[:20]}.get(aa, "ALA")
        linhas.append(_linha_pdb(serial, "CA", resname, "B", j + 1,
                                 60.0 + j * 6.0, 8.5, 0.0, "C"))
        serial += 1
    # ligante: 4 átomos sobre os resíduos 11-13 da cadeia A, a 4,0 Å
    for k, (x, y, z, el, nome) in enumerate([
        (60.0, 4.0, 0.0, "C", "C1"), (66.0, 4.0, 0.0, "C", "C2"),
        (72.0, 4.0, 0.0, "N", "N1"), (66.0, 4.0, 2.0, "O", "O1"),
    ]):
        linhas.append(_linha_pdb(serial, nome, "LIG", "L", 1, x, y, z, el, hetero=True))
        serial += 1
    # vizinhança não-proteica: Zn²⁺ na cadeia da proteína (como num
    # metaloenzima real), uma água em cadeia própria e um glicerol de
    # cristalização. Os três caem no raio e cada um numa categoria diferente.
    extras = [
        ("ZN", "ZN", "A", 201, 66.0, 4.0, -3.0, "ZN"),
        ("O", "HOH", "W", 1, 69.0, 4.0, 0.0, "O"),
        ("C1", "GOL", "L", 2, 72.0, 4.0, 4.0, "C"),
        ("C2", "GOL", "L", 2, 73.2, 4.0, 4.8, "C"),
        ("C3", "GOL", "L", 2, 74.4, 4.0, 5.6, "C"),
    ]
    for nome, resname, cad, num, x, y, z, el in extras:
        linhas.append(_linha_pdb(serial, nome, resname, cad, num, x, y, z, el,
                                 hetero=True))
        serial += 1
    linhas.append("END")
    pdb = destino / "modelo_demo.pdb"
    pdb.write_text("\n".join(linhas) + "\n")

    # MSA da cadeia A: 60 sequências. Bolso conservado, resto variável, e um
    # bloco de 20 quase-idênticas para o peso de Henikoff ter o que corrigir.
    def mutar(seq, p_bolso, p_resto):
        out = []
        for i, aa in enumerate(seq):
            p = p_bolso if i in bolso_a else p_resto
            out.append(rng.choice(aminos) if rng.random() < p else aa)
        return "".join(out)

    registros = [(f"query_cadeia_A", seq_a)]
    for k in range(20):                                   # bloco redundante
        registros.append((f"redundante_{k:02d}", mutar(seq_a, 0.0, 0.02)))
    for k in range(25):                                   # diversos
        registros.append((f"diverso_{k:02d}", mutar(seq_a, 0.03, 0.45)))
    for k in range(10):                                   # parciais (com gap)
        s = list(mutar(seq_a, 0.03, 0.40))
        for i in range(0, 8):
            s[i] = mh.GAP
        registros.append((f"parcial_{k:02d}", "".join(s)))
    fasta_a = destino / "msa_cadeia_A.fasta"
    fasta_a.write_text("".join(f">{i}\n{s}\n" for i, s in registros))

    registros_b = [("query_cadeia_B", seq_b)]
    for k in range(30):
        registros_b.append((f"homologo_{k:02d}",
                            "".join(rng.choice(aminos) if rng.random() < 0.3 else aa
                                    for aa in seq_b)))
    fasta_b = destino / "msa_cadeia_B.fasta"
    fasta_b.write_text("".join(f">{i}\n{s}\n" for i, s in registros_b))

    cfg = json.loads(json.dumps(MODELO_CONFIG))
    cfg["_comentario"] = ["Config da demo sintética gerada por gerar_demo()."]
    cfg["estrutura"] = str(pdb)
    cfg["ligante"] = {"resname": "LIG", "cadeia": None, "resnum": None,
                      "copia": None, "auto": False}
    cfg["msas"] = [
        {"ferramenta": "demo", "cadeia": "A", "arquivo": str(fasta_a), "formato": "auto"},
        {"ferramenta": "demo", "cadeia": "B", "arquivo": str(fasta_b), "formato": "auto"},
    ]
    cfg["escores"]["n_permutacoes"] = 2000
    config = destino / "config_demo.json"
    config.write_text(json.dumps(cfg, indent=2, ensure_ascii=False))
    return {"pdb": pdb, "msa_a": fasta_a, "msa_b": fasta_b, "config": config,
            "cfg": cfg, "seq_a": seq_a, "seq_b": seq_b,
            "bolso_esperado_a": [i + 1 for i in bolso_a]}


# ---------------------------------------------------------------------------


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Conservação (Henikoff) do MSA + tudo a 5 Å do ligante.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("O que este script responde")[0],
    )
    ap.add_argument("--config", type=Path, default=None,
                    help="JSON com todos os campos (ver --escrever-config)")
    ap.add_argument("--escrever-config", type=Path, default=None,
                    metavar="ARQUIVO.json",
                    help="escrever o modelo de config em branco e sair")
    ap.add_argument("--demo", action="store_true",
                    help="gerar dado sintético e rodar o pipeline inteiro nele")
    ap.add_argument("--saida", type=Path, default=Path("conservacao_bolso"))

    e = ap.add_argument_group("estrutura e ligante")
    e.add_argument("--estrutura", type=Path, default=None, help="PDB ou mmCIF")
    e.add_argument("--ligante-resname", default=None)
    e.add_argument("--ligante-cadeia", default=None)
    e.add_argument("--ligante-resnum", type=int, default=None)
    e.add_argument("--ligante-copia", type=int, default=None)
    e.add_argument("--ligante-auto", action="store_true")
    e.add_argument("--raio", type=float, default=5.0)
    e.add_argument("--modelo", type=int, default=1)
    e.add_argument("--com-hidrogenios", action="store_true")
    e.add_argument("--sem-agua", action="store_true")
    e.add_argument("--sem-estrutura", action="store_true",
                   help="só a análise dos MSAs, sem bolso")

    m = ap.add_argument_group("MSAs")
    m.add_argument("--msa", action="append", metavar="FERRAMENTA:CADEIA:CAMINHO",
                   help="pode repetir: --msa boltz2:A:msa_A.fasta")
    m.add_argument("--normalizacao", default="media", choices=["media", "soma"])
    m.add_argument("--gap-como-tipo", action="store_true")
    m.add_argument("--pesos-so-colunas-query", action="store_true")
    m.add_argument("--max-seqs", type=int, default=None)
    m.add_argument("--semente", type=int, default=0)
    m.add_argument("--permitir-sequencia-unica", action="store_true")

    s = ap.add_argument_group("escores")
    s.add_argument("--escore", default="jsd_janela", choices=list(ESCORES_VALIDOS))
    s.add_argument("--janela", type=int, default=3)
    s.add_argument("--min-ocupacao", type=float, default=0.5)
    s.add_argument("--min-neff-coluna", type=float, default=3.0)
    s.add_argument("--todas-as-colunas", action="store_true",
                   help="usar também as colunas rasas na estatística")
    s.add_argument("--n-permutacoes", type=int, default=10000)

    args = ap.parse_args(argv)

    if args.escrever_config:
        args.escrever_config.parent.mkdir(parents=True, exist_ok=True)
        args.escrever_config.write_text(
            json.dumps(MODELO_CONFIG, indent=2, ensure_ascii=False) + "\n"
        )
        print(f"modelo de config escrito em {args.escrever_config}")
        print("preencha 'estrutura', 'ligante' e a lista 'msas' e rode:")
        print(f"  python {Path(__file__).name} --config {args.escrever_config} "
              "--saida resultados/")
        return 0

    try:
        if args.demo:
            demo_dir = args.saida / "entrada_demo"
            d = gerar_demo(demo_dir)
            print(f"dado sintético em {demo_dir}")
            cfg = d["cfg"]
            cfg["raio_angstrom"] = args.raio
            cfg["escores"]["escore_principal"] = args.escore
        elif args.config:
            cfg = ler_config(args.config)
        else:
            cfg = config_de_argumentos(args)
        r = executar(cfg, args.saida, sem_estrutura=args.sem_estrutura)
    except ValueError as erro:
        print(f"ERRO: {erro}", file=sys.stderr)
        return 2

    b = r.get("bolso")
    if b:
        print(f"\nbolso ({b['raio_angstrom']} Å de {b['ligante'].get('resname')}): "
              f"{b['n_itens_no_raio']} itens, "
              f"cadeias {', '.join(b['cadeias_polimericas_em_contato']) or '—'}"
              + ("  [INTERFACE]" if b["bolso_de_interface"] else ""))
        for c in b["nao_proteicos_no_raio"]:
            print(f"  {c['residuo']:18s} {c['categoria']:24s} {c['distancia_min']:.2f} Å")
    print()
    for m_ in r["msas"]:
        if m_.get("erro"):
            print(f"{m_['rotulo']}: ERRO {m_['erro']}")
            continue
        print(f"{m_['rotulo']}: N={m_['n_sequencias']} "
              f"Neff={m_['neff_henikoff']:.1f} "
              f"Neff/res={m_['neff_por_residuo']:.2f} "
              f"({m_['veredito']})")
    for c in r.get("cadeias", []):
        print(f"{c['ferramenta']}/{c['cadeia']}: {c['estatistica']['veredito']}")
    for a in dict.fromkeys(r.get("avisos", [])):
        print(f"aviso: {a}")
    print(f"\n-> {args.saida}  (relatorio.md, resumo.json, *.csv)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
