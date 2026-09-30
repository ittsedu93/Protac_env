#!/usr/bin/env python
"""
Relatório do pipeline PROTAC/PCSK9 em PDF — tracks CRBN e VHL.

Roda no env `mdtools` (matplotlib, pandas). ~30 s.

    python relatorio_pdf.py                      # detecta tudo sozinho
    python relatorio_pdf.py --out ~/relatorio.pdf

O que ele faz
-------------
Monta um PDF A4 com o que o pipeline produziu até agora: o que foi medido, com
que número, e o que cada número decidiu. As figuras da MD (geradas pelo
`md_figures.py`) entram como imagens; as curvas comparativas entre as duas E3
são plotadas aqui, a partir dos CSVs no disco.

Sobre os números embutidos
--------------------------
Os resultados já medidos estão no código como literais, com a data. Isso é
deliberado: o relatório precisa ser completo mesmo rodando numa máquina que não
tenha mais os arquivos intermediários, e um relatório que se degrada em silêncio
quando falta um CSV é pior que um que declara a fonte de cada número.

Havendo o CSV no disco, ele é preferido ao literal — e a legenda da figura diz
qual dos dois foi usado.

Só matplotlib para o PDF: é a única dependência garantida no env, e um relatório
que não gera por falta de biblioteca não é um relatório.
"""

from __future__ import annotations

import argparse
import csv
import json
import textwrap
from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg          # noqa: E402
import matplotlib.pyplot as plt           # noqa: E402
import numpy as np                        # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages   # noqa: E402

# --------------------------------------------------------------------------
# Paleta e tipografia. As duas cores categóricas passaram os seis testes do
# validador (CVD ΔE 24,7 protan; normal 33,6; contraste ≥ 3:1 na superfície) —
# e são as mesmas das figuras da MD, para o relatório não trocar de linguagem
# visual no meio.
# --------------------------------------------------------------------------
COR = {
    "crbn": "#2a78d6",
    "vhl": "#eb6834",
    "terceira": "#1baf7a",
    "tinta": "#1a1a1a",
    "tinta2": "#4a4a4a",
    "tinta3": "#767676",
    "grade": "#e4e4e2",
    "superficie": "#fcfcfb",
}
A4 = (8.27, 11.69)
MARGEM_X, MARGEM_TOPO, MARGEM_BASE = 0.95, 0.85, 0.75

TAM = {"h1": 15, "h2": 11.5, "h3": 9.5, "p": 8.8, "cap": 7.2, "mono": 7.4}
ALTURA = {"h1": 0.42, "h2": 0.30, "h3": 0.24, "p": 0.165, "mono": 0.150}
LARGURA = {"p": 96, "mono": 104}     # caracteres por linha, medidos no DejaVu


# --------------------------------------------------------------------------
# Motor de páginas: uma lista de blocos que flui e pagina sozinha.
# --------------------------------------------------------------------------
class Relatorio:
    def __init__(self, pdf: PdfPages):
        self.pdf = pdf
        self.fig = None
        self.y = 0.0
        self.pagina = 0

    def nova_pagina(self):
        if self.fig is not None:
            self._rodape()
            self.pdf.savefig(self.fig, facecolor=COR["superficie"])
            plt.close(self.fig)
        self.fig = plt.figure(figsize=A4, facecolor=COR["superficie"])
        self.pagina += 1
        self.y = A4[1] - MARGEM_TOPO

    def _rodape(self):
        self.fig.text(1 - MARGEM_X / A4[0], MARGEM_BASE / 2 / A4[1],
                      f"{self.pagina}", ha="right", va="center",
                      size=TAM["cap"], color=COR["tinta3"])
        self.fig.text(MARGEM_X / A4[0], MARGEM_BASE / 2 / A4[1],
                      "PROTAC/PCSK9 — relatório do pipeline", ha="left",
                      va="center", size=TAM["cap"], color=COR["tinta3"])

    def espaco(self, pol: float):
        if self.y - pol < MARGEM_BASE:
            self.nova_pagina()
        else:
            self.y -= pol

    def _cabe(self, pol: float):
        if self.y - pol < MARGEM_BASE:
            self.nova_pagina()

    def _texto(self, s: str, size, cor, x_pol=None, peso="normal",
               familia="sans-serif"):
        x = (x_pol if x_pol is not None else MARGEM_X) / A4[0]
        self.fig.text(x, self.y / A4[1], s, ha="left", va="top", size=size,
                      color=cor, weight=peso, family=familia)

    def h1(self, s: str):
        # Capítulo novo só quebra a página quando não sobra espaço útil. Antes
        # cada h1 forçava página nova, e seções curtas deixavam três quartos de
        # folha em branco.
        usada = (A4[1] - MARGEM_TOPO) - self.y
        restante = self.y - MARGEM_BASE
        if usada > 0.3 and restante < 4.0:
            self.nova_pagina()
        elif usada > 0.3:
            self.y -= 0.45
        self._cabe(0.75)
        self._texto(s, TAM["h1"], COR["tinta"], peso="bold")
        self.y -= ALTURA["h1"]
        eixo = self.fig.add_axes([MARGEM_X / A4[0], self.y / A4[1],
                                  1 - 2 * MARGEM_X / A4[0], 0.0015])
        eixo.set_facecolor(COR["tinta"])
        eixo.set_xticks([]); eixo.set_yticks([])
        for s_ in eixo.spines.values():
            s_.set_visible(False)
        self.y -= 0.16

    def h2(self, s: str):
        self.espaco(0.14)
        self._cabe(1.15)          # título + ~5 linhas, para não deixar órfã
        self._texto(s, TAM["h2"], COR["tinta"], peso="bold")
        self.y -= ALTURA["h2"]

    def h3(self, s: str):
        self.espaco(0.08)
        self._cabe(0.80)          # idem, em escala menor
        self._texto(s, TAM["h3"], COR["tinta2"], peso="bold")
        self.y -= ALTURA["h3"]

    def p(self, s: str, recuo=0.0):
        for linha in textwrap.wrap(" ".join(s.split()),
                                   width=LARGURA["p"] - int(recuo * 11)) or [""]:
            self._cabe(ALTURA["p"])
            self._texto(linha, TAM["p"], COR["tinta2"],
                        x_pol=MARGEM_X + recuo)
            self.y -= ALTURA["p"]
        self.y -= 0.06

    def item(self, s: str, marca="•"):
        linhas = textwrap.wrap(" ".join(s.split()), width=LARGURA["p"] - 4)
        for i, linha in enumerate(linhas or [""]):
            self._cabe(ALTURA["p"])
            if i == 0:
                self._texto(marca, TAM["p"], COR["tinta3"],
                            x_pol=MARGEM_X + 0.06)
            self._texto(linha, TAM["p"], COR["tinta2"], x_pol=MARGEM_X + 0.26)
            self.y -= ALTURA["p"]

    def mono(self, s: str, cor=None):
        for linha in s.rstrip("\n").split("\n"):
            self._cabe(ALTURA["mono"])
            self._texto(linha[:LARGURA["mono"]], TAM["mono"],
                        cor or COR["tinta2"], familia="monospace")
            self.y -= ALTURA["mono"]
        self.y -= 0.06

    def tabela(self, cabecalho: list[str], linhas: list[list[str]],
               larguras: list[int] | None = None):
        """Tabela em monoespaçado, com o texto QUEBRADO dentro da célula.

        Sem a quebra, uma célula mais larga que a coluna invade a vizinha e a
        linha some pela margem — foi o que aconteceu com o catálogo de erros,
        cujas frases são longas por natureza. A largura total é limitada ao que
        cabe na página, e o que sobra é redistribuído da coluna mais larga.
        """
        n = len(cabecalho)
        larg = list(larguras) if larguras else [
            max(len(str(cabecalho[i])), *(len(str(l[i])) for l in linhas)) + 2
            for i in range(n)]
        # cabe na página? encolhe a maior coluna até caber
        while sum(larg) > LARGURA["mono"]:
            larg[larg.index(max(larg))] -= 1

        # Alinhamento por CONTEÚDO, não por posição: número à direita, frase à
        # esquerda. Alinhar frase à direita deixa a margem esquerda serrilhada
        # e o olho não acha o começo da linha seguinte.
        def e_numero(v):
            return bool(str(v).strip()) and all(
                c.isdigit() or c in " .,%×–-/ÅnsmhdkD" for c in str(v))
        a_direita = [i > 0 and all(e_numero(l[i]) for l in linhas)
                     for i in range(n)]

        def celulas(vals, cabeca=False):
            """Cada célula quebrada; devolve as linhas físicas da linha lógica."""
            colunas = [textwrap.wrap(str(v), width=max(larg[i] - 2, 6)) or [""]
                       for i, v in enumerate(vals)]
            altura = max(len(c) for c in colunas)
            fisicas = []
            for k in range(altura):
                partes = []
                for i, c in enumerate(colunas):
                    txt = c[k] if k < len(c) else ""
                    partes.append(txt.rjust(larg[i]) if a_direita[i]
                                  else txt.ljust(larg[i]))
                fisicas.append("".join(partes))
            return fisicas

        self._cabe(0.55)
        for linha in celulas(cabecalho, cabeca=True):
            self._texto(linha, TAM["mono"], COR["tinta"],
                        familia="monospace", peso="bold")
            self.y -= ALTURA["mono"]
        self._texto("─" * sum(larg), TAM["mono"], COR["grade"],
                    familia="monospace")
        self.y -= ALTURA["mono"] * 0.8
        for l in linhas:
            fisicas = celulas(l)
            # a linha lógica não se parte entre páginas: quebrar no meio de uma
            # célula de três linhas separa a causa do sintoma
            self._cabe(len(fisicas) * ALTURA["mono"])
            for linha in fisicas:
                self._texto(linha, TAM["mono"], COR["tinta2"],
                            familia="monospace")
                self.y -= ALTURA["mono"]
            # respiro depois de TODA linha lógica: sem ele, uma linha de
            # uma física colada numa de três parece uma linha só
            self.y -= ALTURA["mono"] * 0.32
        self.y -= 0.10

    def legenda(self, s: str):
        for linha in textwrap.wrap(" ".join(s.split()), width=118):
            self._cabe(0.14)
            self._texto(linha, TAM["cap"], COR["tinta3"])
            self.y -= 0.14
        self.y -= 0.06

    # ---- figuras ---------------------------------------------------------
    def eixos(self, altura_pol: float, pad_topo=0.40, pad_base=0.50,
              pad_esq=0.62):
        """Reserva espaço e devolve eixos INSET dentro dele.

        O título, os rótulos de eixo e a legenda são desenhados FORA da caixa
        dos eixos. Reservar só a altura da caixa fazia o gráfico invadir o
        parágrafo de cima e a legenda de baixo — que é o que acontecia até
        aqui. As folgas abaixo são o espaço desses elementos, medido nas
        figuras deste relatório.
        """
        total = pad_topo + altura_pol + pad_base
        self._cabe(total + 0.10)
        larg = 1 - (MARGEM_X + pad_esq) / A4[0] - MARGEM_X / A4[0]
        y0 = self.y - pad_topo - altura_pol
        eixo = self.fig.add_axes([(MARGEM_X + pad_esq) / A4[0], y0 / A4[1],
                                  larg, altura_pol / A4[1]])
        self.y -= total + 0.10
        return eixo

    def imagem(self, caminho: Path, altura_pol: float = 3.4):
        try:
            img = mpimg.imread(str(caminho))
        except Exception as e:
            self.p(f"[figura não lida: {caminho.name} — {str(e)[:60]}]")
            return False
        h, w = img.shape[0], img.shape[1]
        larg_pol = 1 * A4[0] - 2 * MARGEM_X
        alt = min(altura_pol, larg_pol * h / w)
        largura_usada = alt * w / h
        self._cabe(alt + 0.1)
        x = (A4[0] - largura_usada) / 2 / A4[0]
        eixo = self.fig.add_axes([x, (self.y - alt) / A4[1],
                                  largura_usada / A4[0], alt / A4[1]])
        eixo.imshow(img)
        eixo.axis("off")
        self.y -= alt + 0.10
        return True

    def fechar(self):
        if self.fig is not None:
            self._rodape()
            self.pdf.savefig(self.fig, facecolor=COR["superficie"])
            plt.close(self.fig)
            self.fig = None


def estiliza(eixo, titulo=None, xlabel=None, ylabel=None):
    """Grade e eixos recessivos; o dado é que tem de ficar visível."""
    eixo.set_facecolor(COR["superficie"])
    for lado in ("top", "right"):
        eixo.spines[lado].set_visible(False)
    for lado in ("left", "bottom"):
        eixo.spines[lado].set_color(COR["grade"])
        eixo.spines[lado].set_linewidth(0.8)
    eixo.tick_params(colors=COR["tinta3"], labelsize=TAM["cap"], length=3,
                     width=0.8)
    eixo.grid(True, color=COR["grade"], linewidth=0.6, alpha=0.9)
    eixo.set_axisbelow(True)
    if titulo:
        eixo.set_title(titulo, size=TAM["h3"], color=COR["tinta"],
                       weight="bold", loc="left", pad=8)
    if xlabel:
        eixo.set_xlabel(xlabel, size=TAM["cap"], color=COR["tinta2"])
    if ylabel:
        eixo.set_ylabel(ylabel, size=TAM["cap"], color=COR["tinta2"])


# --------------------------------------------------------------------------
# Dados medidos. Literais com data, e o disco tem preferência quando existe.
# --------------------------------------------------------------------------
VARREDURA_CRBN = [(10, 0), (12, 0), (14, 1), (15, 7), (16, 14), (17, 19),
                  (18, 40), (20, 77), (25, 295), (30, 639), (100, 8312)]
VARREDURA_VHL = [(10, 204), (12, 423), (14, 774), (16, 1433), (18, 2468),
                 (20, 3666), (22, 5160)]

GRUPO_POR_CORTE = {          # maior grupo mutuamente próximo (clique máxima)
    "CRBN · WH022": [(5, 1), (10, 3), (15, 4), (20, 6)],
    "CRBN · WH023": [(5, 2), (10, 3), (15, 4), (20, 6)],
}

ENERGIA = [               # (rótulo, soluções, abaixo do limiar)
    ("CRBN · WH022\nponte no limite", 658, 21),
    ("CRBN · WH023\nponte com folga", 1184, 134),
]


def ler_varredura(csv_path: Path):
    if not csv_path or not csv_path.exists():
        return None
    pts = []
    with open(csv_path) as fh:
        for r in csv.DictReader(fh):
            try:
                pts.append((float(r["dist_thr_A"]), int(float(r["transforms"]))))
            except (KeyError, ValueError):
                pass
    return sorted(pts) or None


def achar(raiz: Path, padrao: str):
    try:
        return sorted(raiz.rglob(padrao))
    except Exception:
        return []


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", type=Path,
                    default=Path.home() / "PRosettaC_runs"
                    / "vhl_crbn_pcsk9_protac")
    ap.add_argument("--out", type=Path,
                    default=Path.home() / "relatorio_protac_pcsk9.pdf")
    args = ap.parse_args()

    W = args.work.expanduser()
    P_CRBN = W / "pipeline"
    P_VHL = W / "pipeline_vhl"

    # --- o que existe no disco -------------------------------------------
    figs_md = [p for p in [P_CRBN / "md" / "figures" / f
                           for f in ("fig1_rmsd.png", "fig2_rmsf.png",
                                     "fig3_contacts.png", "fig4_protac.png")]
               if p.exists()]
    scan_crbn = next(iter(achar(P_CRBN, "span_scan_*.csv")), None)
    scan_vhl = next(iter(achar(P_VHL, "span_scan_*.csv")), None)
    conc = achar(P_CRBN, "concordancia_clusters.csv")
    span_csv = next(iter(achar(P_CRBN, "linker_span.csv")), None)
    req_vhl = P_VHL / "span_requirement.json"

    v_crbn = ler_varredura(scan_crbn) or VARREDURA_CRBN
    v_vhl = ler_varredura(scan_vhl) or VARREDURA_VHL
    de_disco = {"CRBN": scan_crbn is not None, "VHL": scan_vhl is not None}

    args.out.expanduser().parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(str(args.out.expanduser())) as pdf:
        R = Relatorio(pdf)

        # =================================================== CAPA
        R.nova_pagina()
        R.y -= 1.6
        R._texto("PROTAC contra a PCSK9", 26, COR["tinta"], peso="bold")
        R.y -= 0.62
        R._texto("Relatório do pipeline computacional", 14, COR["tinta2"])
        R.y -= 0.40
        R._texto("Tracks CRBN e VHL — do warhead ao complexo ternário",
                 11, COR["tinta3"])
        R.y -= 0.9
        eixo = R.eixos(0.02, pad_topo=0.0, pad_base=0.0, pad_esq=0.0)
        eixo.set_facecolor(COR["crbn"]); eixo.set_xticks([]); eixo.set_yticks([])
        for s_ in eixo.spines.values():
            s_.set_visible(False)
        R.y -= 0.3
        R.p(f"Gerado em {date.today().strftime('%d/%m/%Y')} a partir dos "
            f"resultados em {W}.")
        R.espaco(0.2)
        R.h3("O que este documento contém")
        R.item("as fases do pipeline na ordem em que rodam, e o que cada uma "
               "decide;")
        R.item("os números medidos de cada etapa, com a fonte de cada um;")
        R.item("as figuras da dinâmica molecular do nível (ii);")
        R.item("a comparação entre as duas E3 ligases, que é o resultado "
               "central até aqui;")
        R.item("o catálogo de erros e o que cada um virou no código — a parte "
               "transferível do trabalho.")
        R.espaco(0.35)
        R.h3("Aviso de escopo")
        R.p("Nada aqui é validação experimental. Todos os resultados são "
            "computacionais: docking, predição de estrutura e dinâmica "
            "molecular. Onde um método falha em convergir, isso é relatado "
            "como ausência de geometria determinada pelo método — não como "
            "ausência de complexo.")

        # =================================================== SUMÁRIO
        R.nova_pagina()
        R.h1("Sumário executivo")
        R.p("O pipeline vai de warheads feniletilamina contra a PCSK9 até a "
            "predição do complexo ternário e a dinâmica molecular. Ele foi "
            "executado por inteiro com a E3 ligase CRBN e está em execução "
            "com a VHL. O quadro abaixo é o estado em "
            f"{date.today().strftime('%d/%m/%Y')}.")
        R.tabela(
            ["Etapa", "CRBN", "VHL"],
            [["Warheads e ancoragem na PCSK9", "concluída", "reaproveitada"],
             ["Recrutador da E3 (WP1)", "ligand_116", "ligand_087"],
             ["Linkers e sub-complexos (WP2)", "concluída", "concluída"],
             ["PROTACs montados (WP3)", "105", "120"],
             ["Vão exigido pelo par E3/alvo", "14 Å mín · 18 Å útil", "< 10 Å"],
             ["Portão de alcance do linker", "30 de 105 passam", "aprovado"],
             ["MD nível (ii), 600 ns", "APROVADO", "pendente"],
             ["Ternário — Boltz-2", "negativo", "não aplicado"],
             ["Ternário — PRosettaC", "sem convergência", "em execução"],
             ["MD nível (iii)", "sem pose", "aguarda o ternário"]])
        R.espaco(0.1)
        R.h2("Os três achados que o trabalho já sustenta")
        R.h3("1. A validação binária não prevê a ternária")
        R.p("O candidato SC0006__WH023 passou 600 ns de dinâmica molecular com "
            "RMSD do sítio de 1,70 Å e engajamento da warhead preservado em "
            "0,99× o inicial — e era geometricamente impossível de formar "
            "ternário, o que o PatchDock mostrou depois em 30 minutos. A MD do "
            "nível (ii) não podia perceber: a PCSK9 não está na caixa de "
            "simulação.")
        R.h3("2. O requisito geométrico é medível, e é propriedade do par de "
             "proteínas")
        R.p("A distância mínima entre os dois átomos de conjugação para existir "
            "alguma colocação sem clash não depende do candidato: depende das "
            "duas poses. Mede-se uma vez, com o PatchDock, e vale para o "
            "catálogo inteiro. Para CRBN/PCSK9 são 14 Å (regime utilizável em "
            "18 Å); para VHL/PCSK9, menos de 10 Å.")
        R.h3("3. As duas E3 falham — ou decidem — por motivos diferentes")
        R.p("A CRBN não tem onde pôr a PCSK9: zero colocações a 12 Å. A VHL tem "
            "milhares: 204 já a 10 Å. Com a CRBN o gargalo é geometria; com a "
            "VHL a decisão passa para a energia do Rosetta. É a primeira vez "
            "que o pipeline chega nessa pergunta.")

        # =================================================== ORDEM
        R.h1("O pipeline, na ordem em que roda")
        R.p("Cada fase grava um marcador e é retomável. As fases 6a e 6b são "
            "portões que este projeto aprendeu a ter, e a posição delas é o "
            "ponto principal do relatório: medir o vão exigido e reprovar quem "
            "não alcança ANTES de gastar dinâmica molecular.")
        R.tabela(
            ["Fase", "O que faz", "Custo", "Portão"],
            [["1", "Enumerar warheads feniletilamina", "minutos", "—"],
             ["2", "Ancorar na PCSK9 + revalidar redocking", "horas", "RMSD ≤ 2 Å"],
             ["3", "Analisar e escolher warheads", "segundos", "SD, vetor de saída"],
             ["4", "WP1 — recrutador da E3 + exit vector", "minutos", "—"],
             ["5", "WP2 — linkers Chemspace + geometria", "~1 h", "clash no exit vector"],
             ["6", "WP3 — montar PROTACs, emitir jobs", "minutos", "—"],
             ["6a", "MEDIR o vão exigido (PatchDock)", "~30 min", "controle a 100 Å"],
             ["6b", "REPROVAR quem não alcança", "segundos", "alcance ≥ exigido"],
             ["7", "Ranquear", "segundos", "portão elimina"],
             ["8", "MD nível (ii), 3 × 200 ns", "46 h", "RMSD do sítio, ancoragem"],
             ["9", "Ternário (PRosettaC) + concordância", "horas a dias", "grupo ≥ 5 modelos"],
             ["10", "MD nível (iii) — ternário completo", "dias", "pendente"]])
        R.legenda("No track da VHL as fases 8 e 9 estão invertidas por "
                  "aritmética: o ternário custa horas e foi o que barrou a "
                  "CRBN; a MD do nível (ii) custa 46 h. Gastar 46 h antes de "
                  "saber se existe pose ternária é a troca que o projeto já "
                  "fez uma vez.")

        # =================================================== MD (ii)
        R.h1("Dinâmica molecular do nível (ii) — CRBN + PROTAC")
        R.p("Sistema: CRBN (cristal 4TZ4, cadeia C) com o PROTAC completo "
            "recrutador–linker–warhead. Campo de força AMBER ff14SB para a "
            "proteína e GAFF2/AM1-BCC para o ligante, água TIP3P, caixa cúbica "
            "com 13 Å de distância soluto–borda, PME, 3 réplicas independentes "
            "de 200 ns a 310 K e 1 atm.")
        R.tabela(
            ["Critério", "Valor", "Corte", "Veredito"],
            [["RMSD do sítio (115 resíduos)", "1,70 Å", "≤ 3,5", "OK"],
             ["RMSD do PROTAC", "2,87 Å", "≤ 5,0", "OK"],
             ["Ancoragem do recrutador", "0,77 × inicial", "≥ 0,5", "OK"],
             ["Razão warhead/recrutador", "0,65", "≤ 1,0", "OK"],
             ["Crescimento warhead–E3", "0,99 × inicial", "≤ 2,0", "OK"],
             ["Frações de frames ancorado", "1,00", "≥ 0,7", "OK"],
             ["RMSD global da proteína", "4,22 Å", "declarar", "movimento de domínio"]])
        R.p("O número mais informativo é o 0,99×: a warhead termina os 200 ns "
            "com o mesmo engajamento com a CRBN que tinha no início. O risco "
            "central do nível (ii) — a warhead migrar para a própria E3 e "
            "competir com o alvo no ternário — não se materializou. "
            "Desempenho: 310, 329 e 326 ns/dia; ~46 h no total.")
        R.p("O RMSD global de 4,22 Å, presente em 2 das 3 réplicas, é "
            "movimento de domínio do receptor (a dobradiça Lon/TBD da CRBN) e "
            "não perda de estrutura — por isso o portão é o RMSD do SÍTIO, e o "
            "global entra como contexto a declarar.")

        if figs_md:
            R.p("As figuras abaixo vêm do `md_figures.py`, em inglês e a 300 "
                "dpi, cada uma com o CSV dos dados-fonte ao lado do arquivo — "
                "uma figura cujos dados ninguém pode replotar é uma imagem, "
                "não um resultado.")
            TITULOS = {
                "fig1_rmsd": "RMSD vs tempo — sítio, núcleo e proteína inteira",
                "fig2_rmsf": "RMSF por resíduo, com o sítio sombreado",
                "fig3_contacts": "Contatos nas duas interfaces vs tempo",
                "fig4_protac": "RMSD do PROTAC vs tempo",
            }
            for f in figs_md:
                # título e figura juntos: um sozinho no rodapé separa a legenda
                # do que ela legenda
                R._cabe(3.9)
                R.h3(TITULOS.get(f.stem, f.stem.replace("_", " ")))
                R.imagem(f, altura_pol=3.4)
                R.legenda(f"Fonte: {f.name}")
        else:
            R.h2("Figuras da MD")
            R.p("As figuras não foram encontradas em "
                f"{P_CRBN / 'md' / 'figures'}. Gere-as com:")
            R.mono("python scripts/md_figures.py --md-dir $PIPELINE_OUT/md")
            R.p("e rode este relatório novamente — elas entram "
                "automaticamente.")

        # =================================================== TERNÁRIO CRBN
        R.h1("Complexo ternário — track CRBN")
        R.h2("Boltz-2: um negativo, e ele é informativo")
        R.p("25 amostras (5 lotes de 5; 20 de uma vez estouram os 32 GB da "
            "placa). O critério de aprovação não é o confidence_score global — "
            "esse número é dominado pelo enovelamento das proteínas, que é "
            "fácil e já conhecido. O que decide é o ipTM do pior par entre uma "
            "cadeia da E3 e uma do alvo: abaixo de 0,5 o modelo não sabe onde "
            "as duas proteínas se encontram.")
        R.tabela(
            ["Métrica", "Valor"],
            [["Amostras geradas", "25"],
             ["Aprovadas (ipTM do pior par ≥ 0,5)", "3 (12%)"],
             ["Distância entre as 3 aprovadas", "42,7 / 74,6 / 82,2 Å"],
             ["Pares aprovados a menos de 5 Å", "0 de 3"],
             ["Agrupamento reprodutível encontrado", "10 modelos a 2–16 Å"],
             ["…e aprovado pelo ipTM", "nenhum dos 10"]])
        R.p("As três aprovadas põem a PCSK9 em três lugares diferentes. Não é "
            "uma interface com 12% de taxa de acerto: são interfaces "
            "incompatíveis, e o corte foi satisfeito por acaso três vezes. "
            "Onde o método é reprodutível ele não tem confiança; onde tem "
            "confiança não é reprodutível.")

        R.h2("PRosettaC: dois candidatos, e o que a comparação entre eles "
             "revelou")
        R.p("O PRosettaC parte das duas poses experimentais — o recrutador no "
            "cristal da CRBN e a warhead no da PCSK9 — amostra o linker e "
            "agrupa as soluções. Foram executados dois candidatos com o MESMO "
            "linker (15 ligações, teto de alcance 18,03 Å) e poses de warhead "
            "diferentes, o que isolou a variável que importa.")
        R.tabela(
            ["", "WH022", "WH023"],
            [["Vão exigido pela pose", "~18 Å (no limite)", "14 Å (folga de 4 Å)"],
             ["Transformadas do PatchDock", "19", "34"],
             ["Soluções de docking local", "658", "1184"],
             ["Abaixo do limiar de energia", "21 (3,2%)", "134 (11,3%)"],
             ["Clusters gerados", "20", "105"],
             ["Clusters com ≥ 5 membros", "0", "1"],
             ["Maior grupo a 5 Å", "1 de 20", "2 de 30"]])

        eixo = R.eixos(2.6)
        rotulos = [e[0] for e in ENERGIA]
        frac = [100 * e[2] / e[1] for e in ENERGIA]
        barras = eixo.bar(range(len(frac)), frac, width=0.45,
                          color=[COR["crbn"], COR["crbn"]], linewidth=0)
        for i, (b, f, e) in enumerate(zip(barras, frac, ENERGIA)):
            eixo.text(b.get_x() + b.get_width() / 2, f + 0.4,
                      f"{f:.1f}%", ha="center", va="bottom",
                      size=TAM["cap"], color=COR["tinta"], weight="bold")
            eixo.text(b.get_x() + b.get_width() / 2, f - 0.6,
                      f"{e[2]} de {e[1]}", ha="center", va="top",
                      size=TAM["cap"] - 0.6, color=COR["superficie"])
        eixo.set_xticks(range(len(rotulos)))
        eixo.set_xticklabels(rotulos, size=TAM["cap"], color=COR["tinta2"])
        eixo.set_ylim(0, max(frac) * 1.25)
        estiliza(eixo, "Modelos abaixo do limiar de energia do Rosetta",
                 ylabel="% das soluções de docking local")
        eixo.grid(axis="x", visible=False)
        R.legenda("Mesmo linker nos dois; muda a folga entre o alcance do "
                  "linker e o vão que a pose exige. Fechar a ponte na "
                  "conformação estendida paga entropia e tensão, e a tensão "
                  "aparece como energia desfavorável: o filtro corta 97% dos "
                  "modelos no candidato sem folga e 89% no candidato com "
                  "folga. Fonte: result_summary.txt de cada execução.")

        R.nova_pagina()
        R.h2("Concordância: 105 clusters são 105 interfaces, ou uma região?")
        R.p("Contar clusters não responde a pergunta. O que responde é "
            "superpor cada modelo pela E3 — que é a âncora experimental, vinda "
            "do cristal — e medir o quanto o alvo se desloca entre modelos. O "
            "gráfico mostra o maior grupo MUTUAMENTE próximo (clique máxima) "
            "em função do corte de distância.")

        eixo = R.eixos(2.9)
        # As duas séries coincidem de 10 Å em diante: sem o traço distinto a de
        # baixo desaparece sob a de cima, e o leitor conclui que só há uma.
        estilos = [(COR["crbn"], (0, (5, 2)), "o"), (COR["vhl"], "-", "s")]
        finais = []
        for (rotulo, pts), (cor, traco, marca) in zip(GRUPO_POR_CORTE.items(),
                                                      estilos):
            x = [p[0] for p in pts]; y = [p[1] for p in pts]
            eixo.plot(x, y, linestyle=traco, marker=marca, color=cor,
                      linewidth=2, markersize=5,
                      markeredgecolor=COR["superficie"], markeredgewidth=1.5,
                      label=rotulo)
            finais.append((x[-1], y[-1], rotulo, cor))
        # Rótulo direto só quando os fins não se sobrepõem; sobrepostos, a
        # legenda já resolve a identidade e dois rótulos empilhados só sujam.
        if abs(finais[0][1] - finais[1][1]) > 0.5:
            for xf, yf, rot, cor in finais:
                eixo.annotate(rot.split(" · ")[1], (xf, yf),
                              textcoords="offset points", xytext=(8, -2),
                              size=TAM["cap"], color=cor, weight="bold")
        eixo.axhline(5, color=COR["tinta3"], linewidth=1, linestyle=(0, (4, 3)))
        eixo.annotate("portão: 5 modelos", (5.2, 5.15), size=TAM["cap"],
                      color=COR["tinta3"])
        eixo.set_xlim(3.5, 23)
        eixo.set_ylim(0, 8)
        estiliza(eixo, "Convergência dos modelos ternários (track CRBN)",
                 xlabel="corte de distância entre modelos (Å)",
                 ylabel="maior grupo mutuamente próximo")
        eixo.legend(frameon=False, fontsize=TAM["cap"], loc="upper left",
                    labelcolor=COR["tinta2"])
        R.legenda("Uma interface convergente daria um grupo com metade dos "
                  "modelos no corte de 5 Å. Os dois candidatos ficam em 1 e 2 "
                  "modelos, e mesmo a 20 Å — deslocamento que muda a interface "
                  "inteira de uma proteína de 479 resíduos — o grupo chega a 6. "
                  "Fonte: concordancia_clusters.csv, via prosettac_agreement.py.")

        if conc:
            m = []
            with open(conc[0]) as fh:
                rows = list(csv.reader(fh))
            for r in rows[1:]:
                for v in r[1:]:
                    if v:
                        try:
                            m.append(float(v))
                        except ValueError:
                            pass
            if len(m) > 10:
                eixo = R.eixos(2.4)
                eixo.hist(m, bins=24, color=COR["crbn"], linewidth=0)
                eixo.axvline(5, color=COR["tinta3"], linewidth=1,
                             linestyle=(0, (4, 3)))
                eixo.annotate("5 Å", (5.4, eixo.get_ylim()[1] * 0.9),
                              size=TAM["cap"], color=COR["tinta3"])
                estiliza(eixo,
                         "Distância par a par entre modelos, após superpor "
                         "pela E3",
                         xlabel="RMSD do alvo (Å)", ylabel="pares")
                eixo.grid(axis="x", visible=False)
                R.legenda(f"Fonte: {conc[0]}. A massa da distribuição está "
                          f"entre 20 e 45 Å: os modelos não discordam por "
                          f"pouco, discordam de lado da proteína.")

        R.h2("Veredito do track CRBN")
        R.p("Para o recrutador de CRBN na pose do cristal 4TZ4 e este sítio da "
            "PCSK9, nenhum dos dois métodos produz geometria ternária "
            "testável. Duas linhas independentes, com o controle que elimina "
            "'faltou amostragem' como explicação: quando a amostragem "
            "triplicou (3,2% → 11,3% de aproveitamento energético), a "
            "dispersão não mudou.")
        R.p("O que isso NÃO sustenta: que estes PROTACs não formem complexo "
            "ternário. Ausência de convergência computacional não é ausência "
            "de complexo, e os benchmarks do próprio PRosettaC têm taxa de "
            "acerto limitada. O experimento que separaria as duas afirmações é "
            "um controle positivo — rodar o pipeline num sistema com ternário "
            "cristalográfico resolvido — e ele está listado como próximo passo.")

        # =================================================== A MEDIÇÃO DO VÃO
        R.h1("A medição que reorganizou o pipeline")
        R.p("O PRosettaC restringe o docking global a uma distância máxima "
            "entre os dois átomos de conjugação. Essa distância é o alcance do "
            "linker. Variando a restrição e contando as colocações que "
            "sobrevivem, mede-se o que o par de proteínas exige — e isso é "
            "independente do candidato.")

        eixo = R.eixos(3.2)
        for (rot, pts, cor) in (("CRBN", v_crbn, COR["crbn"]),
                                ("VHL", v_vhl, COR["vhl"])):
            xs = [p[0] for p in pts if p[0] < 100 and p[1] > 0]
            ys = [p[1] for p in pts if p[0] < 100 and p[1] > 0]
            eixo.plot(xs, ys, "-o", color=cor, linewidth=2, markersize=5,
                      markeredgecolor=COR["superficie"], markeredgewidth=1.5,
                      label=rot)
            if xs:
                eixo.annotate(rot, (xs[-1], ys[-1]),
                              textcoords="offset points", xytext=(8, 0),
                              size=TAM["cap"], color=cor, weight="bold")
        zeros = [p[0] for p in v_crbn if p[1] == 0]
        if zeros:
            eixo.plot(zeros, [0.7] * len(zeros), "x", color=COR["crbn"],
                      markersize=6)
            eixo.annotate("CRBN = 0 (fora da escala log)",
                          (min(zeros), 1.0), size=TAM["cap"] - 0.4,
                          color=COR["crbn"])
        eixo.set_yscale("log")
        eixo.set_ylim(0.5, 1e4)
        estiliza(eixo, "Colocações das duas proteínas que sobrevivem à "
                       "restrição de distância",
                 xlabel="distância máxima entre os átomos de conjugação (Å)",
                 ylabel="transformadas do PatchDock (escala log)")
        eixo.legend(frameon=False, fontsize=TAM["cap"], loc="lower right",
                    labelcolor=COR["tinta2"])
        R.legenda(
            "Escala logarítmica porque os valores vão de 0 a 8312. Os pontos "
            "em × são zeros da CRBN, que não existem em escala log. "
            f"Fonte: {'CSV do disco' if de_disco['CRBN'] else 'medição de 28/09'}"
            f" (CRBN) e "
            f"{'CSV do disco' if de_disco['VHL'] else 'medição de 29/09'} (VHL). "
            "A distância de 100 Å é controle e não aparece: sem restrição "
            "efetiva o docking devolve milhares, e zero ali invalidaria a "
            "curva inteira.")

        R.tabela(
            ["Distância (Å)", "CRBN", "VHL", "Razão"],
            [[f"{d}", f"{c}",
              f"{next((v for dd, v in v_vhl if dd == d), '—')}",
              (f"{next((v for dd, v in v_vhl if dd == d), 0) / c:.0f}×"
               if c and next((v for dd, v in v_vhl if dd == d), 0) else
               ("∞" if next((v for dd, v in v_vhl if dd == d), 0) else "—"))]
             for d, c in v_crbn if d in (10, 12, 14, 16, 18, 20)])

        R.h2("O critério que saiu disso")
        R.p("O piso absoluto não serve de critério: com 14 Å o par CRBN/PCSK9 "
            "dava UMA transformada, e uma solução não sustenta cluster. O "
            "regime utilizável começa onde há transformadas suficientes para o "
            "Rosetta refinar — adotamos 20. E o teto de alcance do linker deve "
            "exceder esse vão com margem, para a ponte fechar perto do meio da "
            "distribuição conformacional e não no extremo dela.")
        R.mono("alcance exigido = vão útil / (mediana/teto)\n"
               "                = 18,0 Å / 0,64  =  28,1 Å   (CRBN)\n"
               "                = 10,0 Å / 0,64  =  15,6 Å   (VHL)")
        R.p("A razão mediana/teto de 0,64 foi medida em um candidato (11,6 / "
            "18,03). É um ponto experimental, não uma constante — está como "
            "parâmetro no código, com a origem registrada no JSON de saída.")

        R.h2("Como o alcance é calculado, e o estimador que estava errado")
        R.p("A primeira versão do portão media o alcance como o máximo sobre "
            "confôrmeros gerados por ETKDG. Esse estimador SATURA muito abaixo "
            "do real, porque confôrmero aleatório quase nunca cai na "
            "conformação estendida — e ele reprovou como 'curto' justamente o "
            "candidato que o PRosettaC aceitou.")
        R.tabela(
            ["Estimador", "SC0006__WH023", "SC0013__WH022", "Custo"],
            [["Limite de distance-geometry", "10,89 Å", "18,03 Å", "3 ms"],
             ["Histograma do PRosettaC", "11 Å", "18 Å", "—"],
             ["Máximo de 200 confôrmeros", "9,31 Å", "14,44 Å", "27 s"],
             ["Máximo de 800 confôrmeros", "9,31 Å", "14,44 Å", "108 s"]])
        R.legenda("Quadruplicar a amostra não move o máximo em 0,01 Å. O "
                  "limite de distance-geometry acerta o número do PRosettaC "
                  "nas duas calibrações independentes e é 9000× mais rápido — "
                  "e como o portão existe para prever o que o PRosettaC vai "
                  "fazer, o teto é a medida certa.")

        # =================================================== VHL
        R.h1("Track VHL — em execução")
        R.p("A troca de E3 é configuração, não reescrita: o driver já escolhia "
            "entre VHL e CRBN, e a CRBN havia vencido por score de docking. O "
            "que mudou de verdade foi a ORDEM — os portões de alcance entraram "
            "antes da dinâmica molecular, e o ternário antes da MD do nível "
            "(ii).")
        R.tabela(
            ["Item", "CRBN", "VHL"],
            [["Cristal do recrutador", "4TZ4", "6GFZ"],
             ["Cadeia do receptor preparado", "C (361 res)", "C (141 res)"],
             ["Ligantes na triagem round2", "97", "14"],
             ["Recrutador escolhido", "ligand_116", "ligand_087"],
             ["Átomo de conjugação (âncora)", "1", "15"],
             ["PROTACs montados", "105", "120"],
             ["Massa molecular", "802–953 Da", "911–1102 Da"],
             ["Ligações rotacionáveis", "19–27", "19–27"],
             ["Vão exigido (útil)", "18 Å", "< 10 Å"],
             ["Alcance exigido do linker", "28,1 Å", "15,6 Å (~13 ligações)"],
             ["Transformadas do PatchDock", "19–34", "2464"]])
        R.h2("O que já se pode afirmar")
        R.item("A restrição geométrica que barrou a CRBN não existe na VHL: "
               "204 colocações já a 10 Å, contra zero da CRBN.")
        R.item("A explicação provável é tamanho — o pVHL tem 141 resíduos "
               "contra 361 da CRBN preparada, portanto menos volume estérico "
               "entre os dois sítios. Isso é hipótese; o número é o que está "
               "medido.")
        R.item("O custo dessa permissividade é tempo de máquina: o Rosetta "
               "escala com as transformadas, e 2464 (limitadas a 1000 pelo "
               "PRosettaC) contra 34 significam dias em vez de horas.")
        R.item("A massa molecular subiu: recrutadores de VHL são maiores que a "
               "lenalidomida, e a faixa 911–1102 Da encosta no limite de 1000 "
               "Da que o próprio ranking usa como penalidade de "
               "permeabilidade.")
        R.item("A triagem da VHL tem 14 ligantes contra 97 da CRBN — o "
               "conjunto de recrutadores a escolher é 7× menor, e isso limita "
               "a diversidade do WP1 para esta E3. É limitação a declarar.")
        if req_vhl.exists():
            try:
                d = json.loads(req_vhl.read_text())
                R.h2("Requisito medido no disco")
                R.mono(json.dumps({k: v for k, v in d.items() if k != "curva"},
                                  indent=2, ensure_ascii=False))
            except Exception:
                pass

        # =================================================== ERROS
        R.h1("Catálogo dos erros, e o que cada um virou")
        R.p("Esta é a parte transferível do trabalho. Cada linha é um erro que "
            "não levantava exceção, ou levantava a exceção errada, ou cuja "
            "mensagem apontava para o lugar errado. A coluna da direita é o "
            "que ficou no código — guarda, não instrução: instrução se perde, "
            "guarda não.")
        R.h2("Erros de método")
        R.tabela(
            ["Sintoma", "Causa", "O que ficou"],
            [["teste dizia OK, ferramenta recusava",
              "comparava o arquivo com uma releitura dele mesmo",
              "alvo do teste = versão protonada"],
             ["portão reprovava o que funciona",
              "máximo amostrado satura abaixo do real",
              "limite de distance-geometry"],
             ["'faltou amostragem' sobre 20 clusters",
              "veredito saía da menor distância, não do tamanho do grupo",
              "clique máxima, varredura de cortes"],
             ["46 h de MD num candidato impossível",
              "portão geométrico no fim da fila",
              "fases 6a e 6b antes da MD"],
             ["conclusão de ausência com metade da amostragem",
              "`Full: False` herdado sem razão declarada",
              "default True + aviso no lançamento"]],
            larguras=[34, 40, 30])
        R.h2("Erros de acoplamento")
        R.tabela(
            ["Sintoma", "Causa", "O que ficou"],
            [["OverflowError 50 linhas depois",
              "head com 0 H fixado nos carbonos; match vazio → −1",
              "conserto verificado"],
             ["morria no `conda activate`",
              "SETVARS_CALL não definida sob `set -eu`",
              "export antes de lançar"],
             ["'rodando' com a fila vazia",
              "`grep` casava o cabeçalho do squeue",
              "`squeue -h`"],
             ["etapa pulada por arquivo de 0 byte",
              "limpeza nomeava produtos um a um",
              "regra inversa: fica o que é entrada"],
             ["config pedia cadeia A, PDB tinha C",
              "default plausível em 104 configs",
              "detectar no arquivo, ou exigir"],
             ["cadeia no config, não passada ao emissor",
              "fiação faltando entre driver e script",
              "preflight compara config × arquivo"],
             ["1490 arquivos apagados",
              "limpeza checava Results, não execução em andamento",
              "guarda com três sinais, código 7"],
             ["[NÃO SUBIU] com log vazio",
              "`sed` de arquivo vazio imprime nada",
              "rastro com `bash -x` + o wrapper"]],
            larguras=[34, 40, 30])
        R.h2("O padrão")
        R.p("Três coisas atravessam a lista. Primeiro: API que devolve 'não "
            "fiz nada' com o mesmo tipo de 'fiz' — a defesa é verificar a "
            "saída, não confiar na chamada. Segundo: o sintoma raramente fica "
            "perto da causa; um OverflowError num índice negativo tinha origem "
            "em contagem de hidrogênio três funções antes. Terceiro: um teste "
            "que compara uma coisa consigo mesma passa sempre e não testa "
            "nada.")

        # =================================================== LIMITAÇÕES
        R.h1("Limitações e próximos passos")
        R.h2("O que este trabalho não estabelece")
        R.item("Nenhuma validação experimental foi feita. Todos os resultados "
               "são computacionais.")
        R.item("Ausência de convergência ternária não é ausência de complexo. "
               "Sem um controle positivo, não é possível separar limitação do "
               "método de propriedade do sistema.")
        R.item("A razão mediana/teto (0,64) que converte vão exigido em "
               "alcance exigido vem de um único candidato medido.")
        R.item("O vão mínimo da VHL está abaixo da faixa varrida: os 15,6 Å de "
               "alcance exigido são limite superior, não medida.")
        R.item("O sítio na PCSK9 é um só, escolhido na fase 3. O requisito "
               "geométrico depende da pose da warhead, e 161 poses foram "
               "geradas — apenas duas tiveram o requisito medido.")
        R.h2("Próximos passos, em ordem de valor por hora de máquina")
        R.tabela(
            ["#", "O que", "Custo", "O que resolve"],
            [["1", "Varrer o vão da VHL abaixo de 10 Å", "minutos",
              "linker mais curto = menos massa"],
             ["2", "Concluir o ternário da VHL (em execução)", "~1 dia",
              "existe pose para o nível (iii)?"],
             ["3", "Controle positivo com ternário cristalográfico", "horas",
              "separa método de sistema"],
             ["4", "MD nível (ii) do candidato da VHL", "46 h",
              "interação espúria warhead–E3"],
             ["5", "MD nível (iii) sobre a pose aprovada", "dias",
              "estabilidade do ternário"],
             ["6", "Medir o requisito para outras poses de warhead", "horas",
              "sítio que dispense linker longo"]])
        R.h2("Reprodutibilidade")
        R.p("Todo o pipeline está em um repositório git, com cada fase "
            "retomável por marcador em disco e cada decisão registrada no "
            "commit que a implementou. O preflight confere, em segundos, cada "
            "dependência de cada fase antes de qualquer execução pesada — "
            "incluindo a comparação entre valores de configuração e o conteúdo "
            "real dos arquivos.")
        R.mono("PIPELINE_CONF=config/pipeline_vhl.conf bash scripts/preflight.sh\n"
               "PIPELINE_CONF=config/pipeline_vhl.conf bash scripts/run_pipeline.sh --from 4")

        R.fechar()

    print(f"  {args.out.expanduser()}")
    nfig = len(figs_md)
    print(f"  figuras da MD incluídas: {nfig}"
          + ("" if nfig else "  (rode o md_figures.py e refaça)"))
    print(f"  varredura CRBN: {'CSV do disco' if de_disco['CRBN'] else 'literal'}")
    print(f"  varredura VHL:  {'CSV do disco' if de_disco['VHL'] else 'literal'}")


if __name__ == "__main__":
    main()
