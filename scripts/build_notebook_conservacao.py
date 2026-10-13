#!/usr/bin/env python
"""
Monta o notebook notebooks/conservacao_msa_cofold.ipynb.

Gerado por script para ser reproduzível: para mudar o notebook, edite aqui e
rode de novo.

    python scripts/build_notebook_conservacao.py

O notebook é fino de propósito: a conta mora em scripts/msa_henikoff.py,
scripts/bolso_ligante.py e scripts/conservacao_bolso.py, que têm teste
(tests/test_conservacao_msa.py). O notebook configura, chama e desenha.
"""

import json
from pathlib import Path

CELLS = []


def md(src: str):
    CELLS.append({"cell_type": "markdown", "metadata": {},
                  "source": src.strip("\n").splitlines(keepends=True)})


def code(src: str):
    CELLS.append({"cell_type": "code", "execution_count": None, "metadata": {},
                  "outputs": [], "source": src.strip("\n").splitlines(keepends=True)})


# ===========================================================================
md(r"""
# Conservação dos MSAs de cofold (pesos de Henikoff) + bolso a 5 Å do ligante

Avalia **quanto se pode confiar** no MSA que alimentou um modelo de cofold
(Boltz-2, Chai-1/2, ESMFold, OpenFold) e **se o bolso do ligante cai em região
conservada**.

Três perguntas, nesta ordem:

| # | Pergunta | Onde |
|---|---|---|
| 1 | O MSA é profundo e não-redundante? | §3 — Neff de Henikoff, N/Neff, ocupação |
| 2 | O que está a menos de 5 Å do ligante? | §4 — resíduos, **outras cadeias, cofatores, íons**, glicanos, aditivos, águas |
| 3 | O bolso é mais conservado que o resto da cadeia? | §5 — Mann-Whitney + permutação |

## Por que pesos de Henikoff

Os MSAs desses servidores vêm de bancos **redundantes**. Dez mil ortólogos de
primata não são dez mil observações independentes, e contar frequência crua
faz a coluna parecer conservada porque o banco está enviesado — não porque a
evolução a protegeu.

Os pesos posicionais de Henikoff & Henikoff (1994) corrigem isso sem árvore
filogenética nem limiar de identidade: numa coluna com `r` tipos de resíduo,
um resíduo do tipo `a` que aparece `s` vezes contribui `1/(r·s)` para o peso da
sequência que o carrega. Quem só repete o que todos já dizem acumula pouco.

O notebook relata o efeito em número: **N/Neff** é o fator de redundância, e a
tabela "efeito dos pesos" mostra a conservação média com e sem peso. Sem os
pesos, a conservação vem inflada nessa proporção.

## Escore de conservação

O principal é a **JSD** (Jensen-Shannon divergence) contra o fundo de
Robinson & Robinson, suavizada em janela — o método que ganhou a comparação de
Capra & Singh (2007) para detecção de sítio funcional. Ela sabe que uma coluna
toda de triptofano informa mais que uma coluna toda de leucina, porque W é
raro. Entropia de Shannon, conservação normalizada e KL saem também, por
coluna, no CSV.

## Onde achar o MSA de cada ferramenta

| Ferramenta | O que ela usa | Formato |
|---|---|---|
| **Boltz-2** | o `.a3m`/`.csv` do YAML, ou o que o `--use_msa_server` baixou (procure `*.a3m` / `*.csv` no diretório de saída) | a3m, csv |
| **Chai-1 / Chai-2** | `*.aligned.pqt`, ou o a3m do ColabFold se você gerou por fora | parquet, a3m |
| **OpenFold / AlphaFold** | o diretório `msas/`: `*.sto` (jackhmmer, UniRef90/MGnify) e `*.a3m` (hhblits, BFD/Uniclust) | sto, a3m |
| **ESMFold** | **nenhum** — é single-sequence por construção | — |

O script lê os quatro formatos direto, além de **FASTA alinhado**, que é o
formato recomendado quando você exportar. Para o ESMFold não há conservação a
medir: use o MSA de outra ferramenta para a mesma cadeia se quiser o perfil
(§3 avisa em vez de inventar número).

---
**Pré-requisito:** nada além da biblioteca padrão. `numpy` acelera, `pandas` e
`matplotlib` deixam as tabelas e os gráficos bonitos — os três são opcionais e
o notebook roda sem eles.
""")

# ---------------------------------------------------------------------------
md(r"""
## 1. Caminhos do projeto

Só isto precisa existir para o notebook importar os scripts.
""")

code(r"""
from pathlib import Path
import json
import sys

# raiz do projeto: ajuste se o notebook não estiver em Protac_env/notebooks/
PROJETO = Path.home() / "Protac_env"
if not (PROJETO / "scripts").is_dir():
    PROJETO = Path.cwd().parent if (Path.cwd().parent / "scripts").is_dir() else Path.cwd()

sys.path.insert(0, str(PROJETO / "scripts"))

import msa_henikoff as mh
import bolso_ligante as bl
import conservacao_bolso as cb

# `display` existe no kernel do Jupyter; fora dele (ao rodar este notebook
# como script, ou com jupyter nbconvert em ambiente mínimo) cai para print
try:
    display
except NameError:
    display = print

print("projeto:", PROJETO)
print("numpy  :", "sim" if mh._np is not None else "não (usa python puro, mais lento)")
try:
    import pandas as pd
    TEM_PANDAS = True
except ImportError:
    TEM_PANDAS = False
try:
    import matplotlib
    if not hasattr(sys, "ps1") and "ipykernel" not in sys.modules:
        matplotlib.use("Agg")      # nó sem tela (nbconvert, SLURM)
    import matplotlib.pyplot as plt
    TEM_PLOT = True
except ImportError:
    TEM_PLOT = False
print(f"pandas : {'sim' if TEM_PANDAS else 'não'}   matplotlib: {'sim' if TEM_PLOT else 'não'}")
""")

# ---------------------------------------------------------------------------
md(r"""
## 2. ENTRADAS — preencha aqui

**Todos os campos abaixo são genéricos.** Enquanto os resultados não
chegarem, deixe `USAR_EXEMPLO_SINTETICO = True`: o notebook gera um caso de
teste com resposta conhecida e roda inteiro, provando que a instalação está
boa. Quando o cofold terminar, ponha `False` e preencha os caminhos.

### O que vai em cada campo

- `ESTRUTURA` — o modelo do cofold. **Boltz-2 e Chai** entregam `.cif`;
  **ESMFold e OpenFold**, `.pdb`. Os dois formatos são lidos.
- `LIGANTE` — como achar o ligante no arquivo. `resname="LIG"` cobre Boltz e
  Chai; num co-cristal é o código do HET (ex.: `"063"` no 6U26). Se houver
  várias cópias, o script **para e lista** as cópias em vez de escolher uma
  por acaso — aí você passa `cadeia`/`resnum`, ou `copia=1`.
  Com `auto=True` ele pega o maior hetero-resíduo que não é água nem íon.
- `MSAS` — **um por (ferramenta, cadeia)**. `cadeia` é a cadeia da ESTRUTURA
  que aquele MSA descreve; é por ela que o perfil é mapeado nos resíduos.
  Pode repetir a mesma cadeia com ferramentas diferentes: aí sai a tabela de
  concordância (§6).
- `RAIO` — 5,0 Å, como pedido. Mude se quiser outra vizinhança.
""")

code(r"""
# ===========================================================================
# ENTRADAS
# ===========================================================================
USAR_EXEMPLO_SINTETICO = True        # True = roda a demo; False = usa o que está abaixo

# --- estrutura do modelo de cofold ------------------------------------------
ESTRUTURA = None                     # ex.: Path.home() / "resultados/boltz/pred.model_idx_0.cif"
MODELO_NMR = 1                       # qual MODEL usar (ensembles/NMR)

# --- ligante ----------------------------------------------------------------
LIGANTE = {
    "resname": "LIG",                # ex.: "LIG" (Boltz/Chai), "063" (co-cristal 6U26)
    "cadeia": None,                   # ex.: "C" — só se precisar desambiguar
    "resnum": None,                   # ex.: 1
    "copia": None,                     # ex.: 1 — qual cópia, quando há várias
    "auto": False,                     # True = pega o maior hetero automaticamente
}

# --- raio do bolso ----------------------------------------------------------
RAIO = 5.0                            # Å
INCLUIR_AGUA = True                   # listar águas (modelo de cofold não tem; cristal tem)
COM_HIDROGENIOS = False               # distância entre átomos pesados é a convenção

# --- MSAs: um por (ferramenta, cadeia) --------------------------------------
# formato aceito: FASTA alinhado (.fasta/.afa), a3m, Stockholm (.sto),
#                 csv do Boltz, parquet do Chai (.aligned.pqt)
MSAS = [
    {"ferramenta": "boltz2",   "cadeia": "A", "arquivo": None},   # ex.: ".../msa_alvo.fasta"
    {"ferramenta": "boltz2",   "cadeia": "B", "arquivo": None},   # ex.: ".../msa_e3.fasta"
    # {"ferramenta": "chai1",    "cadeia": "A", "arquivo": None},
    # {"ferramenta": "openfold", "cadeia": "A", "arquivo": None},
]

# --- pesos de Henikoff ------------------------------------------------------
PESOS = {
    "normalizacao": "media",     # "media" (protege hits parciais) ou "soma" (clássico)
    "gap_como_tipo": False,      # gap NÃO é o 21º tipo de resíduo
    "so_colunas_query": False,   # pesar só nas colunas em que a query tem resíduo
    "max_seqs": None,            # ex.: 20000 para amostrar MSA gigante
    "semente": 0,                # amostragem e permutação determinísticas
}

# --- escores ----------------------------------------------------------------
ESCORES = {
    "escore_principal": "jsd_janela",   # jsd | jsd_janela | jsd_penalizado |
                                        # conservacao_shannon | kl_bits | jsd_sem_peso
    "janela": 3,                        # meia-janela da suavização
    "peso_janela": 0.5,
    "min_ocupacao": 0.5,                # coluna "confiável": ocupação ponderada mínima
    "min_neff_coluna": 3.0,             #                     e Neff mínimo da coluna
    "so_colunas_confiaveis": True,      # estatística só nas colunas confiáveis
    "n_permutacoes": 10000,
    "comparar_sem_peso": True,          # calcular também sem peso, para comparar
}

PERMITIR_SEQUENCIA_UNICA = False        # True só se você quiser forçar o caso ESMFold

# --- saída ------------------------------------------------------------------
SAIDA = Path.home() / "conservacao_cofold"
""")

# ---------------------------------------------------------------------------
md(r"""
### 2.1 Montar e validar a configuração

A validação roda **antes** de qualquer conta e lista **todos** os problemas de
uma vez: campo em branco, arquivo que não existe, par ferramenta/cadeia
repetido, escore que não existe. É de propósito: descobrir o terceiro erro
depois de três execuções é pior que ver os três agora.
""")

code(r"""
if USAR_EXEMPLO_SINTETICO:
    DEMO = cb.gerar_demo(SAIDA / "exemplo_sintetico")
    cfg = DEMO["cfg"]
    cfg["raio_angstrom"] = RAIO
    cfg["escores"].update(ESCORES)
    print("MODO EXEMPLO: dado sintético em", SAIDA / "exemplo_sintetico")
    print("Resposta conhecida: bolso = resíduos 11-13 da cadeia A (+3 da B),")
    print("Zn/água/glicerol no raio, e o bolso TEM de sair mais conservado.\n")
else:
    cfg = {
        "estrutura": str(ESTRUTURA) if ESTRUTURA else None,
        "modelo_nmr": MODELO_NMR,
        "ligante": LIGANTE,
        "raio_angstrom": RAIO,
        "com_hidrogenios": COM_HIDROGENIOS,
        "incluir_agua": INCLUIR_AGUA,
        "msas": [dict(m, arquivo=(str(m["arquivo"]) if m["arquivo"] else None))
                 for m in MSAS],
        "pesos": PESOS,
        "escores": ESCORES,
        "permitir_sequencia_unica": PERMITIR_SEQUENCIA_UNICA,
    }
    cfg = cb.mesclar_config(cb.MODELO_CONFIG, cfg)

erros = cb.validar_config(cfg)
if erros:
    print("CONFIGURAÇÃO INCOMPLETA — corrija na célula anterior:")
    for e in erros:
        print("  -", e)
else:
    print("configuração ok")
    print(f"  estrutura: {cfg['estrutura']}")
    print(f"  ligante:   {cfg['ligante']}")
    print(f"  raio:      {cfg['raio_angstrom']} Å")
    for m in cfg["msas"]:
        print(f"  MSA {m['ferramenta']:10s} cadeia {m['cadeia']}: {m['arquivo']}")
""")

# ---------------------------------------------------------------------------
md(r"""
### 2.2 Guardar a configuração em disco (opcional)

Útil para rodar a mesma análise fora do notebook, no terminal:

```bash
python scripts/conservacao_bolso.py --config config_conservacao.json --saida resultados/
```
""")

code(r"""
SAIDA.mkdir(parents=True, exist_ok=True)
caminho_cfg = SAIDA / "config_conservacao.json"
caminho_cfg.write_text(json.dumps(cfg, indent=2, ensure_ascii=False), encoding="utf-8")
print("configuração salva em", caminho_cfg)
print(f"\npara rodar no terminal:\n  python {PROJETO}/scripts/conservacao_bolso.py "
      f"--config {caminho_cfg} --saida {SAIDA}")
""")

# ---------------------------------------------------------------------------
md(r"""
## 3. O MSA é confiável?

O cartão de confiabilidade de cada MSA. Como ler:

- **Neff** — número efetivo de sequências pelos pesos de Henikoff (Kish).
  É quantas observações independentes o MSA realmente tem.
- **Neff/resíduo** — a medida que importa. `< 1` é raso (e é a faixa em que o
  próprio cofold erra mais); `1–10` razoável; `> 10` profundo.
- **N/Neff** — redundância. `50x` significa que, sem os pesos, a conservação
  estaria inflada nessa proporção.
- **colunas confiáveis** — quantas colunas da query passam de `min_ocupacao` e
  `min_neff_coluna`. Coluna rasa pode marcar 100% de conservação por acidente.
""")

code(r"""
analises = {}
falhas = []
for m in cfg["msas"]:
    rotulo = f"{m['ferramenta']}__{m['cadeia']}"
    try:
        analises[(m["ferramenta"], m["cadeia"])] = mh.analisar_msa(
            m["arquivo"], rotulo=rotulo, formato=m.get("formato", "auto"),
            gap_como_tipo=cfg["pesos"]["gap_como_tipo"],
            normalizacao=cfg["pesos"]["normalizacao"],
            so_colunas_query=cfg["pesos"]["so_colunas_query"],
            meia_janela=cfg["escores"]["janela"],
            peso_janela=cfg["escores"]["peso_janela"],
            min_ocupacao=cfg["escores"]["min_ocupacao"],
            min_neff_coluna=cfg["escores"]["min_neff_coluna"],
            max_seqs=cfg["pesos"]["max_seqs"], semente=cfg["pesos"]["semente"],
            permitir_unica=cfg.get("permitir_sequencia_unica", False),
            comparar_sem_peso=cfg["escores"].get("comparar_sem_peso", True),
        )
    except ValueError as erro:
        falhas.append((rotulo, str(erro)))

for rotulo, erro in falhas:
    print(f"[FALHOU] {rotulo}: {erro}\n")

tabela = []
for (ferr, cad), a in analises.items():
    s = a["resumo"]
    tabela.append({
        "ferramenta": ferr, "cadeia": cad, "formato": s["formato"],
        "N": s["n_sequencias"], "L_query": s["comprimento_query"],
        "Neff": round(s["neff_henikoff"], 1),
        "Neff/res": round(s["neff_por_residuo"], 2),
        "N/Neff": round(s["redundancia_n_sobre_neff"], 1),
        "ocup_media": round(s["ocupacao_media_colunas_query"], 2),
        "col_confiaveis": f"{s['colunas_query_confiaveis']}/{s['comprimento_query']}",
        "JSD_media": round(s["jsd_media_colunas_query"], 4),
        "JSD_sem_peso": round(s.get("jsd_media_sem_peso", float('nan')), 4),
        "efeito_pesos": round(s.get("efeito_dos_pesos_jsd", float('nan')), 4),
        "veredito": s["veredito"],
    })

if TEM_PANDAS and tabela:
    df_msa = pd.DataFrame(tabela)
    display(df_msa)
else:
    for t in tabela:
        print(t)

for (ferr, cad), a in analises.items():
    for aviso in a["resumo"]["avisos"]:
        print(f"aviso {ferr}/{cad}: {aviso}")
""")

# ---------------------------------------------------------------------------
md(r"""
### 3.1 Distribuição dos pesos de Henikoff

O histograma mostra o que a redundância faz. Peso 1,0 é o peso uniforme (eles
somam N). Uma cauda grossa perto de zero = um bloco de sequências quase
idênticas que, sem peso, votaria N vezes na mesma conservação.

`pesos_henikoff__*.csv` traz o peso de cada sequência com o id, o número de
resíduos e a identidade com a query — dá para ver **quem** foi rebaixado.
""")

code(r"""
if TEM_PLOT and analises:
    n = len(analises)
    fig, eixos = plt.subplots(1, n, figsize=(5 * n, 3.4), squeeze=False)
    for ax, ((ferr, cad), a) in zip(eixos[0], analises.items()):
        pesos = a["pesos"]
        ax.hist(pesos, bins=40, color="#4878a8", edgecolor="white")
        ax.axvline(1.0, color="crimson", ls="--", lw=1,
                   label="peso uniforme (1,0)")
        s = a["resumo"]
        ax.set_title(f"{ferr} / cadeia {cad}\nN={s['n_sequencias']}  "
                     f"Neff={s['neff_henikoff']:.0f}  "
                     f"N/Neff={s['redundancia_n_sobre_neff']:.1f}x", fontsize=9)
        ax.set_xlabel("peso de Henikoff")
        ax.set_ylabel("sequências")
        ax.legend(fontsize=7)
    fig.tight_layout()
    plt.show()
elif analises:
    for (ferr, cad), a in analises.items():
        p = sorted(a["pesos"])
        print(f"{ferr}/{cad}: peso min={p[0]:.3f} mediana={p[len(p)//2]:.3f} "
              f"max={p[-1]:.3f}")
""")

# ---------------------------------------------------------------------------
md(r"""
## 4. O bolso: tudo a menos de 5 Å do ligante

A lista é **completa e etiquetada**: resíduos da própria cadeia, resíduos de
**outras cadeias**, **cofatores**, **íons**, glicanos, lipídeos/detergentes,
aditivos de cristalização (glicerol, EDO, sulfato…) e **águas**. Cada item com
categoria, distância mínima e o par de átomos que a realiza.

Três coisas para olhar:

1. **Quantas cadeias poliméricas** o ligante toca. Num ternário de PROTAC,
   tocar só uma significa que o modelo não montou o ternário — montou um
   binário com a outra proteína pendurada.
2. **Íons e cofatores.** Um Zn²⁺ ou Mg²⁺ a 3 Å muda a química do sítio e
   desaparece de qualquer lista que filtre HETATM.
3. **Aditivos.** Em co-cristal, metade da vizinhança costuma ser glicerol e
   sulfato. É informação — desde que etiquetada, e não contada como "bolso".
""")

code(r"""
lig_cfg = cfg["ligante"]
bolso = bl.analisar_bolso(
    cfg["estrutura"], resname=lig_cfg.get("resname"), cadeia=lig_cfg.get("cadeia"),
    resnum=lig_cfg.get("resnum"), auto=bool(lig_cfg.get("auto")),
    copia=lig_cfg.get("copia"), raio=cfg["raio_angstrom"],
    modelo=cfg.get("modelo_nmr", 1),
    com_hidrogenios=bool(cfg.get("com_hidrogenios")),
    incluir_agua=bool(cfg.get("incluir_agua", True)),
)
d, s = bolso["descricao"], bolso["resumo"]
print(f"ligante escolhido: {d['resname']}  cadeia {d['cadeia']}  resnum {d['resnum']}  "
      f"({d['n_atomos']} átomos; seleção por {d['como']})")
print(f"\nestrutura: {s['inventario_estrutura']['n_atomos']} átomos, "
      f"{s['inventario_estrutura']['n_residuos']} resíduos, "
      f"cadeias: " + ", ".join(
          f"{c}({v['tipo']}, {v['n_residuos']} res)"
          for c, v in sorted(s["inventario_estrutura"]["cadeias"].items())))
print(f"\na <= {s['raio_angstrom']} Å do ligante: {s['n_itens_no_raio']} itens "
      f"({s['n_residuos_proteicos']} resíduos proteicos), "
      f"mínima {s['distancia_minima']:.2f} Å" if s["n_itens_no_raio"] else
      f"\nNADA a <= {s['raio_angstrom']} Å — confira se o ligante é o certo")
print("por categoria:", s["por_categoria"])
print("cadeias poliméricas em contato:",
      ", ".join(s["cadeias_polimericas_em_contato"]) or "(nenhuma)",
      "  <<< BOLSO DE INTERFACE" if s["bolso_de_interface"] else "")

if TEM_PANDAS and bolso["contatos"]:
    df_bolso = pd.DataFrame(bolso["contatos"])[[
        "residuo", "cadeia", "resname", "resnum", "categoria", "subtipo",
        "aminoacido", "mesma_cadeia_que_o_ligante", "n_atomos_em_contato",
        "distancia_min", "atomo_mais_proximo", "atomo_ligante_mais_proximo",
        "atomos"]].round({"distancia_min": 2})
    display(df_bolso)
else:
    for c in bolso["contatos"]:
        print(f"  {c['residuo']:16s} {c['categoria']:20s} "
              f"{c['distancia_min']:.2f} Å  {c['n_atomos_em_contato']} átomos")
""")

# ---------------------------------------------------------------------------
md(r"""
### 4.1 O que NÃO é proteína no raio

Separado porque é o que costuma escapar: cofatores, íons, glicanos, aditivos e
águas. Nenhum deles tem coluna de MSA (não são resíduos da cadeia), então eles
entram na descrição do bolso mas **ficam fora** da estatística de conservação.
""")

code(r"""
if s["nao_proteicos_no_raio"]:
    if TEM_PANDAS:
        display(pd.DataFrame(s["nao_proteicos_no_raio"]).round({"distancia_min": 2}))
    else:
        for c in s["nao_proteicos_no_raio"]:
            print(f"  {c['residuo']:16s} {c['categoria']:20s} {c['subtipo']:24s} "
                  f"{c['distancia_min']:.2f} Å")
else:
    print(f"nenhum cofator, íon, glicano, aditivo ou água a <= {cfg['raio_angstrom']} Å")

# contatos átomo-a-átomo (o detalhe fino, para conferir geometria)
if TEM_PANDAS and bolso["pares"]:
    display(pd.DataFrame(bolso["pares"]).head(25).round({"distancia": 2}))
""")

# ---------------------------------------------------------------------------
md(r"""
## 5. O bolso é conservado?

Para cada (ferramenta, cadeia): a query do MSA é alinhada à sequência daquela
cadeia da estrutura (Needleman-Wunsch com gap afim e pontas livres), o perfil
de conservação é mapeado resíduo a resíduo, e os resíduos do bolso são
comparados ao **resto da mesma cadeia**.

O diagnóstico do alinhamento é parte do resultado, não enfeite: **identidade
abaixo de 90% quer dizer que aquele MSA provavelmente não é dessa cadeia** —
o erro mais fácil de cometer quando se tem um MSA por entidade.

Dois testes, de propósito:

- **Mann-Whitney** (aproximação normal) — rápido, mas grosseiro com n pequeno.
- **Permutação** (determinística, com semente) — é o que vale quando o bolso
  tem 8 resíduos.
- **Cliff δ** — tamanho do efeito, de −1 a +1. É o número que diz se a
  diferença importa, e não só se ela é "significativa".
""")

code(r"""
escore = cfg["escores"]["escore_principal"]
resultados_cadeia = []
residuos_por_cadeia = {}

for (ferr, cad), a in analises.items():
    try:
        linhas, diag = cb.residuos_da_cadeia(a, bolso["atomos"], cad,
                                             bolso["contatos"])
    except ValueError as erro:
        print(f"[{ferr}/{cad}] {erro}\n")
        continue
    est = cb.estatistica_do_bolso(
        linhas, escore=escore,
        so_confiaveis=cfg["escores"]["so_colunas_confiaveis"],
        semente=cfg["pesos"]["semente"], n_perm=cfg["escores"]["n_permutacoes"],
    )
    residuos_por_cadeia[(ferr, cad)] = linhas
    resultados_cadeia.append({"ferramenta": ferr, "cadeia": cad,
                              "diagnostico": diag, "estatistica": est})

    print(f"=== {ferr} / cadeia {cad} ===")
    print(f"  alinhamento query x estrutura: identidade "
          f"{diag['identidade_query_estrutura']:.1%} em {diag['pares_alinhados']} pares; "
          f"{diag['residuos_sem_coluna_no_msa']} resíduos sem coluna")
    print(f"  resíduos no bolso: {diag['residuos_no_bolso']}")
    print(f"  {escore}: bolso mediana {est['mediana_bolso']:.4f} "
          f"(n={est['n_bolso']})  resto {est['mediana_resto']:.4f} (n={est['n_resto']})")
    print(f"  Cliff δ={est['mann_whitney']['cliff_delta']:+.2f}  "
          f"p(MW)={est['mann_whitney']['p']:.4f}  "
          f"p(perm)={est['permutacao']['p']:.4f}")
    print(f"  >>> {est['veredito']}")
    for aviso in diag["avisos"]:
        print(f"  (!) {aviso}")
    print()
""")

# ---------------------------------------------------------------------------
md(r"""
### 5.1 Resíduos do bolso, do mais conservado ao menos

`coluna_confiavel = não` quer dizer que aquela coluna do MSA é rasa demais
(pouca ocupação ou Neff baixo) — o número existe, mas não sustenta conclusão.
""")

code(r"""
for (ferr, cad), linhas in residuos_por_cadeia.items():
    do_bolso = [l for l in linhas if l["no_bolso"]]
    if not do_bolso:
        print(f"{ferr}/{cad}: nenhum resíduo desta cadeia no raio do ligante\n")
        continue
    do_bolso.sort(key=lambda l: (l[escore] if isinstance(l[escore], float) else -1),
                  reverse=True)
    print(f"=== {ferr} / cadeia {cad} — {len(do_bolso)} resíduos no bolso ===")
    if TEM_PANDAS:
        display(pd.DataFrame(do_bolso)[[
            "cadeia", "resnum", "resname", "aa_estrutura", "coluna_msa",
            "dist_min_ligante", "atomos_em_contato", "profundidade",
            "ocupacao_pond", "neff_coluna", escore, "jsd_sem_peso",
            "coluna_confiavel"]].round(4))
    else:
        for l in do_bolso:
            print(f"  {l['resname']}{l['resnum']:<5} {escore}={l[escore]:.4f} "
                  f"dist={l['dist_min_ligante']} confiável={l['coluna_confiavel']}")
    print()
""")

# ---------------------------------------------------------------------------
md(r"""
### 5.2 Perfil ao longo da cadeia

A linha é a conservação por resíduo; as barras verticais marcam os resíduos do
bolso. A linha pontilhada é a mesma conservação **sem** os pesos de Henikoff —
a distância entre as duas é a redundância do banco, vista de perto.
""")

code(r"""
if TEM_PLOT and residuos_por_cadeia:
    n = len(residuos_por_cadeia)
    fig, eixos = plt.subplots(n, 1, figsize=(11, 2.8 * n), squeeze=False)
    for ax, ((ferr, cad), linhas) in zip(eixos[:, 0], residuos_por_cadeia.items()):
        xs = [l["resnum"] for l in linhas if isinstance(l[escore], float)]
        ys = [l[escore] for l in linhas if isinstance(l[escore], float)]
        ax.plot(xs, ys, lw=1.2, color="#2b5d8a", label=f"{escore} (Henikoff)")
        sp = [(l["resnum"], l["jsd_sem_peso"]) for l in linhas
              if isinstance(l.get("jsd_sem_peso"), float)]
        if sp:
            ax.plot([p[0] for p in sp], [p[1] for p in sp], lw=0.9, ls=":",
                    color="gray", label="JSD sem peso")
        for l in linhas:
            if l["no_bolso"]:
                ax.axvspan(l["resnum"] - 0.45, l["resnum"] + 0.45,
                           color="orange", alpha=0.35)
        for l in linhas:
            if l["no_bolso"] and isinstance(l[escore], float):
                ax.annotate(f"{l['resname']}{l['resnum']}",
                            (l["resnum"], l[escore]), fontsize=7,
                            rotation=90, ha="center", va="bottom")
        ax.set_title(f"{ferr} / cadeia {cad} — laranja = a <= {cfg['raio_angstrom']} Å "
                     "do ligante", fontsize=10)
        ax.set_xlabel("resíduo (numeração da estrutura)")
        ax.set_ylabel(escore)
        ax.set_ylim(0, 1)
        ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    plt.show()
else:
    print("matplotlib ausente: o perfil está nos CSVs residuos__*.csv")
""")

# ---------------------------------------------------------------------------
md(r"""
### 5.3 Bolso x resto, em distribuição

O boxplot é o teste estatístico desenhado. Se as caixas se sobrepõem quase
inteiramente, o `p` pequeno (quando houver) está vindo de n grande, não de
efeito grande — e é por isso que o Cliff δ aparece no título.
""")

code(r"""
if TEM_PLOT and resultados_cadeia:
    n = len(resultados_cadeia)
    fig, eixos = plt.subplots(1, n, figsize=(3.6 * n, 3.6), squeeze=False)
    for ax, r in zip(eixos[0], resultados_cadeia):
        linhas = residuos_por_cadeia[(r["ferramenta"], r["cadeia"])]
        conf = cfg["escores"]["so_colunas_confiaveis"]
        def vals(no_bolso):
            return [l[escore] for l in linhas
                    if l["no_bolso"] == no_bolso and isinstance(l[escore], float)
                    and (not conf or l["coluna_confiavel"] == 1)]
        dados = [vals(1), vals(0)]
        try:                                  # matplotlib >= 3.9
            ax.boxplot(dados, tick_labels=["bolso", "resto"], showmeans=True)
        except TypeError:                     # matplotlib < 3.9
            ax.boxplot(dados, labels=["bolso", "resto"], showmeans=True)
        for k, v in enumerate(dados, start=1):
            if v:
                ax.scatter([k] * len(v), v, s=12, alpha=0.5,
                           color="orange" if k == 1 else "#2b5d8a")
        e = r["estatistica"]
        ax.set_title(f"{r['ferramenta']}/{r['cadeia']}\n"
                     f"δ={e['mann_whitney']['cliff_delta']:+.2f}  "
                     f"p_perm={e['permutacao']['p']:.3f}", fontsize=9)
        ax.set_ylabel(escore)
        ax.set_ylim(0, 1)
    fig.tight_layout()
    plt.show()
""")

# ---------------------------------------------------------------------------
md(r"""
## 6. As ferramentas concordam?

Com dois MSAs da **mesma cadeia** vindos de origens diferentes (Boltz via
ColabFold vs. OpenFold via jackhmmer, por exemplo), ρ de Spearman alto com
Neff muito diferente tem leitura clara: o perfil é o mesmo, mas uma das buscas
é mais rasa — e a mais profunda é a que vale.

ρ baixo é alarme: os dois MSAs não descrevem a mesma família, e pelo menos um
deles pegou sequências que não deveria.
""")

code(r"""
por_cadeia = {}
for (ferr, cad), a in analises.items():
    por_cadeia.setdefault(cad, []).append((ferr, a))

comparacoes = []
for cad, itens in sorted(por_cadeia.items()):
    for i in range(len(itens)):
        for j in range(i + 1, len(itens)):
            f1, a1 = itens[i]
            f2, a2 = itens[j]
            s1 = {l["pos_query"]: l for l in a1["colunas"] if l["pos_query"] != ""}
            s2 = {l["pos_query"]: l for l in a2["colunas"] if l["pos_query"] != ""}
            comuns = sorted(set(s1) & set(s2))
            x = [s1[p][escore] for p in comuns]
            y = [s2[p][escore] for p in comuns]
            comparacoes.append({
                "cadeia": cad, "A": f1, "B": f2, "posicoes": len(comuns),
                "Neff_A": round(a1["resumo"]["neff_henikoff"], 1),
                "Neff_B": round(a2["resumo"]["neff_henikoff"], 1),
                "media_A": round(sum(x) / len(x), 4) if x else None,
                "media_B": round(sum(y) / len(y), 4) if y else None,
                "spearman": round(cb.spearman(x, y), 4) if len(comuns) >= 3 else None,
            })

if not comparacoes:
    print("só um MSA por cadeia: nada a comparar.\n"
          "Para comparar, acrescente outra entrada em MSAS com a MESMA cadeia "
          "e outra 'ferramenta'.")
elif TEM_PANDAS:
    display(pd.DataFrame(comparacoes))
else:
    for c in comparacoes:
        print(c)
""")

# ---------------------------------------------------------------------------
md(r"""
## 7. Exportar tudo

Roda o pipeline completo pelo script (a mesma conta das células acima, com a
mesma configuração) e escreve os arquivos. Vale como conferência: se o script
e o notebook discordassem, apareceria aqui.

Sai no diretório de saída:

| arquivo | o que tem |
|---|---|
| `relatorio.md` | o relatório completo, legível, com os vereditos |
| `resumo.json` | tudo em forma de máquina |
| `contatos_5A_residuos.csv` | **a lista do bolso**: um item por resíduo/cofator/íon/água, com categoria e distância |
| `contatos_5A_atomos.csv` | par de átomos a par de átomos |
| `bolso.pml` | seleção para abrir no PyMOL e conferir com o olho |
| `conservacao_colunas__*.csv` | por coluna do MSA: entropia, JSD, KL, ocupação, Neff |
| `pesos_henikoff__*.csv` | o peso de cada sequência, com id e identidade |
| `residuos__*.csv` | por resíduo da cadeia: conservação + se está no bolso |
| `comparacao_ferramentas.csv` | concordância entre ferramentas (se houver 2+) |
""")

code(r"""
resultado = cb.executar(cfg, SAIDA)

print("arquivos escritos em", SAIDA)
for p in sorted(SAIDA.glob("*")):
    if p.is_file():
        print(f"  {p.name:42s} {p.stat().st_size / 1024:8.1f} kB")

print("\n--- relatorio.md ---\n")
print((SAIDA / "relatorio.md").read_text(encoding="utf-8"))
""")

# ---------------------------------------------------------------------------
md(r"""
## 8. Como ler o resultado — e onde ele engana

**Conservação alta no bolso** é a evidência barata de que o cofold pôs o
ligante onde a evolução guarda sítio funcional. Com `Neff/resíduo > 1` e as
colunas marcadas como confiáveis, é um argumento.

**Conservação baixa no bolso** tem duas leituras, e a §4 as distingue:

- **bolso de interface** (duas ou mais cadeias poliméricas a 5 Å): esperado.
  A interface induzida por um PROTAC é uma neo-interface — ela não existe sem
  o ligante, logo não há pressão evolutiva ancestral para conservá-la. Aqui
  conservação baixa **não** é defeito do modelo.
- **bolso de uma cadeia só, conservação baixa**: suspeito. O modelo
  provavelmente encaixou o ligante num bolso de superfície qualquer.

**Três armadilhas que o relatório marca, e que valem reler:**

1. `identidade query x estrutura < 90%` — o MSA é de outra cadeia. Todo o
   mapeamento abaixo disso é ficção.
2. `coluna_confiavel = 0` — a coluna é rasa. "100% conservada" com três
   sequências não-gap é acidente, não evidência.
3. `Neff/resíduo < 1` — o MSA é raso. É exatamente a faixa em que o cofold
   também erra mais: conservação e estrutura são, as duas, pouco confiáveis
   pelo mesmo motivo. Não trate um como validação do outro.

**O que este notebook não faz:** não mede a qualidade do modelo (isso é
pLDDT/ipTM/PAE — ver `select_ternary_pose.py`), e não decide se o ternário
existe. Conservação do bolso é uma evidência independente, para somar.

---

## Checklist

- [ ] `USAR_EXEMPLO_SINTETICO = True` rodou inteiro sem erro (instalação ok)
- [ ] `ESTRUTURA` aponta para o modelo certo, e o `MODELO_NMR` é o que se quer
- [ ] o ligante escolhido em §4 é o PROTAC/ligante, e não um aditivo
- [ ] um MSA por cadeia, e a identidade query x estrutura ≥ 90% em todos
- [ ] `Neff/resíduo` registrado para cada MSA (é o que calibra tudo)
- [ ] §4 conferida: cadeias, **cofatores e íons** no raio, e aditivos separados
- [ ] o bolso é de interface? Se sim, a leitura da §8 muda
- [ ] `relatorio.md` e os CSVs arquivados junto com o modelo de cofold
""")

# ===========================================================================
nb = {
    "cells": CELLS,
    "metadata": {
        "kernelspec": {"display_name": "mdtools", "language": "python",
                       "name": "python3"},
        "language_info": {"name": "python", "version": "3.11"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}

out = Path(__file__).resolve().parent.parent / "notebooks" / \
    "conservacao_msa_cofold.ipynb"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(nb, indent=1, ensure_ascii=False))
n_code = sum(1 for c in CELLS if c["cell_type"] == "code")
print(f"{out}\n{len(CELLS)} células ({n_code} de código, "
      f"{len(CELLS) - n_code} markdown)")
