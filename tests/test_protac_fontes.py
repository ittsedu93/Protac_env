#!/usr/bin/env python
"""
Testes de protac_fontes.py — a coleta dos parâmetros da aba PAPERS.

Só unittest da biblioteca padrão, como o resto do projeto:

    python -m unittest discover -s tests -v
    python tests/test_protac_fontes.py

O que está coberto, e por que estes casos
-----------------------------------------
Nenhuma das quatro bases é acessível daqui (rede fechada), então os testes NÃO
tocam a rede: eles cobrem o que é meu, que é o caminho do dado depois que o
arquivo chega. Os fixtures abaixo são tabelas sintéticas pequenas, escritas de
propósito com os nomes de coluna em grafias DIFERENTES entre si — 'DC50 (nM)',
'DC₅₀', 'dc50_nm' — porque o ponto do módulo é justamente não depender da
grafia. Se as bases mudarem os nomes, é este arquivo que diz se o mapa
aguentou.

- **Disputa de colunas**: os três pares que podem se roubar são testados um a
  um — 'Linker SMILES' não pode virar o SMILES do PROTAC, 'Target Class' não
  pode virar o alvo, e 'CL' (clearance) não pode casar dentro de 'Cell Line'.
  Esse último é o mais traiçoeiro: um regex sem limite de palavra erra, e o
  erro sai como um número de clearance que é, na verdade, nome de linhagem.
- **Colunas sem destino**: uma coluna de origem que o mapa não reconheceu tem
  de aparecer na auditoria. É o contrário do defeito que importa evitar: um
  parâmetro que a base TEM vindo vazio na planilha, parecendo ausência de dado.
- **Unidades**: '>10 μM', '0.013 uM' e '13 nM' têm de significar a mesma ordem
  de grandeza depois de convertidos, e o qualificador '>' não pode se perder.
- **Dmax**: '0.95' e '95%' têm de sair iguais, porque as bases discordam da
  escala e a planilha é lida por gente.
- **Consolidação**: o mesmo composto em duas bases tem de virar uma linha, com
  a medida experimental ganhando da prevista (ordem de preferência).
- **Contaminação entre bases**: o padrão 'protac*.csv' do PROTAC-DB casa com
  'protac_bench.csv'. O teste garante que não come o arquivo da outra base.
- **Exportação**: a aba PAPERS tem de sair com as 28 colunas na ordem exata da
  planilha, cabeçalho na linha 3 e dados da linha 4 — e o modelo não pode ser
  sobrescrito nunca.
"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ / "scripts"))

import pandas as pd                      # noqa: E402
import protac_fontes as pf               # noqa: E402


# ===========================================================================
# Fixtures: as quatro bases como elas PLAUSIVELMENTE se parecem, cada uma com
# a sua grafia. Não são cópias das bases reais (que não dá para acessar daqui);
# são o formato, que é o que o módulo precisa aguentar.
# ===========================================================================

PROTACDB_CSV = (
    "Compound ID,Smiles,Target,E3 ligase,DC50 (nM),Dmax (%),IC50 (nM),"
    "Molecular Weight,Cell Line,Article DOI,Linker Type,PAMPA,Caco-2,"
    "T1/2 (h),Cmax (ng/mL),AUC (h*ng/mL),CL (mL/min/kg),Vss (L/kg)\n"
    "1,CC(=O)Nc1ccccc1,BRD4,VHL,13,95,4.2,845.9,MV4-11,10.1021/aaa,PEG,"
    "1.2e-6,8.5,2.1,320,1450,22,1.8\n"
    "2,CCOc1ccccc1N,BRD4,CRBN,>10 uM,45,120,912.4,HeLa,10.1038/bbb,alkyl,,,,,,,\n"
    "3,CNc1ccc(F)cc1,PCSK9,VHL,0.013 uM,0.87,,1024.1,HepG2,10.1016/ccc,PEG,,,,,,,\n"
)

# PROTACpedia escreve com subscrito e usa 'Cells' em vez de 'Cell Line'.
PROTACPEDIA_CSV = (
    "Name,SMILES,Uniprot,Target,E3 Ligase,DC₅₀,Dmax,Cells,Article DOI,"
    "Linker Type,Target Class\n"
    "MZ1,CC(=O)Nc1ccccc1,O60885,BRD4,VHL,13 nM,95%,MV4-11,10.1021/aaa,PEG,"
    "Epigenetic reader\n"
    "ARV-110,c1ccccc1O,P10275,AR,CRBN,1 nM,98%,VCaP,10.1021/ddd,alkyl,"
    "Nuclear receptor\n"
)

# PROTAC-Bench: minúsculo, com underscore, e sem bibliografia.
PROTACBENCH_CSV = (
    "smiles,uniprot_id,e3_ligase,active,dc50_nm,mw\n"
    "CC(=O)Nc1ccccc1,O60885,VHL,1,13.0,845.9\n"
    "CCCCCCO,Q9NZJ5,Other,0,,500.2\n"
)

# PROTAC-PatentDB: propriedades PREVISTAS, com nomes do ADMETlab.
PATENTDB_CSV = (
    "compound_id,SMILES,Patent Number,Target,Publication Date,MW,"
    "Caco-2 Permeability,PAMPA,Target Family\n"
    "P-1,CC(=O)Nc1ccccc1,WO2019/123456,BRD4,2019-06-01,845.90,-5.2,-4.8,"
    "Bromodomain\n"
    "P-2,C1CCCCC1N,WO2020/999999,PCSK9,2020-01-15,760.30,-6.1,-5.5,Protease\n"
)

# Tabela de FRAGMENTO: o módulo tem de ignorar, não misturar com os PROTACs.
PROTACDB_WARHEAD_CSV = (
    "Warhead ID,Warhead Smiles,Target\n"
    "w1,c1ccccc1,BRD4\n"
)


def escrever(base: Path, nome: str, texto: str) -> Path:
    caminho = base / nome
    caminho.parent.mkdir(parents=True, exist_ok=True)
    caminho.write_text(texto, encoding="utf-8")
    return caminho


def montar_dir_dados(raiz: Path) -> Path:
    """Monta dados/bases com um arquivo por base, cada um na sua pasta."""
    dados = raiz / "bases"
    escrever(dados / "PROTAC-DB", "protac.csv", PROTACDB_CSV)
    escrever(dados / "PROTACpedia", "protacpedia_dump.csv", PROTACPEDIA_CSV)
    escrever(dados / "PROTAC-Bench", "protac_bench.csv", PROTACBENCH_CSV)
    escrever(dados / "PROTAC-PatentDB", "protac_patentdb.csv", PATENTDB_CSV)
    return dados


# ===========================================================================
class TestNormalizacaoDeNomes(unittest.TestCase):

    def test_subscrito_e_pontuacao_convergem(self):
        """As três grafias de DC50 têm de virar a mesma chave."""
        for grafia in ("DC₅₀ (nM)", "DC50 [nM]", "dc50_nm", "Dc50   nM"):
            self.assertEqual(pf.normalizar_nome(grafia), "dc50 nm", grafia)

    def test_acento_cai(self):
        self.assertEqual(pf.normalizar_nome("Título Paper "), "titulo paper")
        self.assertEqual(pf.normalizar_nome("Comprimeto do Linker"),
                         "comprimeto do linker")

    def test_vazio_e_none(self):
        self.assertEqual(pf.normalizar_nome(None), "")
        self.assertEqual(pf.normalizar_nome("   "), "")


# ===========================================================================
class TestMapeamentoDeColunas(unittest.TestCase):

    def test_protacdb(self):
        df = pd.read_csv(pd.io.common.StringIO(PROTACDB_CSV))
        mapa, aud = pf.mapear_colunas(df.columns)
        self.assertEqual(mapa["smiles"], "Smiles")
        self.assertEqual(mapa["dc50"], "DC50 (nM)")
        self.assertEqual(mapa["dmax"], "Dmax (%)")
        self.assertEqual(mapa["ic50"], "IC50 (nM)")
        self.assertEqual(mapa["mw"], "Molecular Weight")
        self.assertEqual(mapa["doi"], "Article DOI")
        self.assertEqual(mapa["e3"], "E3 ligase")
        self.assertEqual(mapa["linhagem"], "Cell Line")
        self.assertEqual(mapa["linker_tipo"], "Linker Type")
        self.assertEqual(mapa["alvo"], "Target")
        self.assertEqual(aud["sem_destino"], [],
                         "esta tabela foi escrita para casar inteira")

    def test_protacpedia_com_subscrito(self):
        df = pd.read_csv(pd.io.common.StringIO(PROTACPEDIA_CSV))
        mapa, _ = pf.mapear_colunas(df.columns)
        self.assertEqual(mapa["dc50"], "DC₅₀")
        self.assertEqual(mapa["linhagem"], "Cells")
        self.assertEqual(mapa["uniprot"], "Uniprot")
        self.assertEqual(mapa["classe_alvo"], "Target Class")
        self.assertEqual(mapa["alvo"], "Target",
                         "'Target Class' não pode ter roubado 'Target'")

    def test_cl_nao_casa_dentro_de_cell_line(self):
        """O caso traiçoeiro: clearance vs. linhagem celular.

        Um regex sem limite de palavra casaria 'cl' dentro de 'cell line', e o
        erro sairia como um nome de linhagem na coluna de clearance.
        """
        mapa, _ = pf.mapear_colunas(["Cell Line", "Cells", "Cell Type"])
        self.assertNotIn("cl", mapa)
        self.assertEqual(mapa["linhagem"], "Cell Line")

        mapa2, _ = pf.mapear_colunas(["Cell Line", "CL (mL/min/kg)"])
        self.assertEqual(mapa2["cl"], "CL (mL/min/kg)")
        self.assertEqual(mapa2["linhagem"], "Cell Line")

    def test_linker_smiles_nao_vira_smiles_do_protac(self):
        mapa, _ = pf.mapear_colunas(["Linker SMILES", "SMILES"])
        self.assertEqual(mapa["linker_smiles"], "Linker SMILES")
        self.assertEqual(mapa["smiles"], "SMILES")

    def test_linker_smiles_sozinho_nao_preenche_smiles(self):
        """Se a base só tem o SMILES do linker, a coluna SMILES fica vazia —
        e não com o linker no lugar do PROTAC."""
        mapa, _ = pf.mapear_colunas(["Linker SMILES", "Target"])
        self.assertEqual(mapa["linker_smiles"], "Linker SMILES")
        self.assertNotIn("smiles", mapa)

    def test_coluna_sem_destino_aparece_na_auditoria(self):
        """É o defeito que mais importa pegar: dado que existe e passou batido."""
        mapa, aud = pf.mapear_colunas(["SMILES", "Degradation Rate Constant",
                                       "Ternary Kd"])
        self.assertIn("Degradation Rate Constant", aud["sem_destino"])
        self.assertIn("Ternary Kd", aud["sem_destino"])

    def test_parametros_ausentes_sao_listados(self):
        _, aud = pf.mapear_colunas(["SMILES", "Target"])
        self.assertIn("dc50", aud["sem_origem"])
        self.assertIn("auc", aud["sem_origem"])

    def test_uma_coluna_alimenta_um_parametro_so(self):
        mapa, _ = pf.mapear_colunas(["Target", "Cell Line", "E3 Ligase"])
        self.assertEqual(len(set(mapa.values())), len(mapa))

    def test_tmax_e_cmax_nao_se_confundem(self):
        mapa, _ = pf.mapear_colunas(["Tmax (h)", "Cmax (ng/mL)", "T1/2 (h)",
                                     "MRT (h)", "Vz (L/kg)", "Vss (L/kg)",
                                     "AUC (h*ng/mL)"])
        self.assertEqual(mapa["tmax"], "Tmax (h)")
        self.assertEqual(mapa["cmax"], "Cmax (ng/mL)")
        self.assertEqual(mapa["t_meia"], "T1/2 (h)")
        self.assertEqual(mapa["mrt"], "MRT (h)")
        self.assertEqual(mapa["vz"], "Vz (L/kg)")
        self.assertEqual(mapa["vss"], "Vss (L/kg)")
        self.assertEqual(mapa["auc"], "AUC (h*ng/mL)")


# ===========================================================================
class TestQuantidades(unittest.TestCase):

    def test_unidades_viram_nm(self):
        self.assertAlmostEqual(pf.ler_quantidade("13 nM").em_nm(), 13.0)
        self.assertAlmostEqual(pf.ler_quantidade("0.013 uM").em_nm(), 13.0)
        self.assertAlmostEqual(pf.ler_quantidade("0,013 μM").em_nm(), 13.0)
        self.assertAlmostEqual(pf.ler_quantidade("13000 pM").em_nm(), 13.0)

    def test_micro_grego_e_micro_sinal_valem_o_mesmo(self):
        """'μM' não pode ser lido como 'M'.

        O mu some num filtro ASCII ingênuo, e aí 0,013 μM (13 nM) sai como
        13 000 000 nM — erro de um milhão de vezes, sem nenhuma mensagem.
        """
        for grafia in ("0.013 uM", "0.013 \u03bcM", "0.013 \u00b5M",
                       "0.013 umol/L"):
            self.assertAlmostEqual(pf.ler_quantidade(grafia).em_nm(), 13.0,
                                   msg=grafia)
        self.assertAlmostEqual(pf.ler_quantidade("1 M").em_nm(), 1e9)
        self.assertAlmostEqual(pf.ler_quantidade("1 mM").em_nm(), 1e6)

    def test_qualificador_nao_se_perde(self):
        q = pf.ler_quantidade(">10 uM")
        self.assertEqual(q.qualificador, ">")
        self.assertAlmostEqual(q.em_nm(), 10000.0)
        self.assertTrue(q.texto().startswith(">"))

    def test_qualificador_unicode(self):
        self.assertEqual(pf.ler_quantidade("≥5 nM").qualificador, ">=")
        self.assertEqual(pf.ler_quantidade("≈5 nM").qualificador, "~")

    def test_faixa(self):
        q = pf.ler_quantidade("1.2-3.4 nM")
        self.assertAlmostEqual(q.v1, 1.2)
        self.assertAlmostEqual(q.v2, 3.4)
        self.assertIn("–", q.texto())

    def test_sem_unidade_nao_inventa(self):
        q = pf.ler_quantidade("13")
        self.assertAlmostEqual(q.v1, 13.0)
        self.assertEqual(q.unidade, "")
        self.assertIsNone(q.em_nm(), "sem unidade não há conversão honesta")

    def test_vazios(self):
        for valor in (None, "", "   ", "NA", "n.d.", "ND", float("nan")):
            self.assertTrue(pf.ler_quantidade(valor).vazia, repr(valor))

    def test_texto_junto_do_numero(self):
        q = pf.ler_quantidade("13 nM in MV4-11")
        self.assertAlmostEqual(q.v1, 13.0)
        self.assertEqual(q.unidade, "nM")

    def test_dmax_escala(self):
        """'0.95' e '95%' são a mesma coisa; a saída tem de ser a mesma."""
        self.assertEqual(pf.ler_dmax("0.95").texto(), pf.ler_dmax("95%").texto())
        self.assertEqual(pf.ler_dmax("95").texto(), "95 %")
        self.assertAlmostEqual(pf.ler_dmax("0.87").v1, 87.0)
        self.assertAlmostEqual(pf.ler_dmax("98%").v1, 98.0)

    def test_dmax_1_0_e_fracao(self):
        self.assertAlmostEqual(pf.ler_dmax("1.0").v1, 100.0)

    def test_limpar_texto(self):
        self.assertIs(pf.limpar_texto("NA"), pd.NA)
        self.assertIs(pf.limpar_texto(None), pd.NA)
        self.assertEqual(pf.limpar_texto("  BRD4   humano "), "BRD4 humano")


# ===========================================================================
class TestLeituraDeArquivos(unittest.TestCase):

    def test_sdf_le_propriedades(self):
        sdf = (
            "MZ1\n  dummy\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n"
            "> <Smiles>\nCC(=O)Nc1ccccc1\n\n"
            "> <DC50 (nM)>\n13\n\n"
            "> <Target>\nBRD4\n\n"
            "$$$$\n"
        )
        df = pf.ler_sdf(sdf)
        self.assertEqual(len(df), 1)
        self.assertEqual(df.iloc[0]["Smiles"], "CC(=O)Nc1ccccc1")
        self.assertEqual(df.iloc[0]["DC50 (nM)"], "13")
        self.assertEqual(df.iloc[0]["_sdf_name"], "MZ1")

    def test_sdf_mapeia_como_csv(self):
        """Um SDF e um CSV com as mesmas propriedades têm de dar o mesmo mapa."""
        sdf = ("x\n\n\n> <Smiles>\nCCO\n\n> <DC50 (nM)>\n5\n\n$$$$\n")
        mapa_sdf, _ = pf.mapear_colunas(pf.ler_sdf(sdf).columns)
        mapa_csv, _ = pf.mapear_colunas(["Smiles", "DC50 (nM)"])
        self.assertEqual(mapa_sdf["smiles"], mapa_csv["smiles"])
        self.assertEqual(mapa_sdf["dc50"], mapa_csv["dc50"])

    def test_le_csv_tsv_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            raiz = Path(tmp)
            escrever(raiz, "a.csv", "SMILES,Target\nCCO,BRD4\n")
            escrever(raiz, "b.tsv", "SMILES\tTarget\nCCO\tBRD4\n")
            escrever(raiz, "c.json", '[{"SMILES": "CCO", "Target": "BRD4"}]')
            for nome in ("a.csv", "b.tsv", "c.json"):
                df = pf.ler_tabela(raiz / nome)
                self.assertEqual(list(df.columns), ["SMILES", "Target"], nome)
                self.assertEqual(df.iloc[0]["Target"], "BRD4", nome)


# ===========================================================================
class TestLocalizacaoDeArquivos(unittest.TestCase):

    def test_base_nao_come_arquivo_da_outra(self):
        """'protac*.csv' do PROTAC-DB casa com 'protac_bench.csv'.

        Sem o descarte por pasta, o PROTAC-DB ingeriria o arquivo do
        PROTAC-Bench e a coluna 'fonte' passaria a mentir.
        """
        with tempfile.TemporaryDirectory() as tmp:
            dados = montar_dir_dados(Path(tmp))
            for nome_base in pf.FONTES:
                achados = pf.arquivos_locais(pf.FONTES[nome_base], dados)
                self.assertTrue(achados, f"{nome_base}: nada encontrado")
                for caminho in achados:
                    self.assertIn(nome_base, caminho.parts,
                                  f"{nome_base} pegou {caminho}")

    def test_sem_download_nao_toca_a_rede(self):
        with tempfile.TemporaryDirectory() as tmp:
            arquivos, registro = pf.obter_arquivos(
                pf.FONTES["PROTAC-DB"], Path(tmp) / "vazio", baixar=False)
            self.assertEqual(arquivos, [])
            self.assertTrue(any("download desligado" in r for r in registro))

    def test_base_com_cadastro_nao_tenta_baixar(self):
        """PROTACpedia exige login: a mensagem tem de dizer onde pôr o arquivo."""
        with tempfile.TemporaryDirectory() as tmp:
            arquivos, registro = pf.obter_arquivos(
                pf.FONTES["PROTACpedia"], Path(tmp) / "vazio", baixar=True)
            self.assertEqual(arquivos, [])
            texto = " ".join(registro)
            self.assertIn("cadastro", texto)
            self.assertIn("PROTACpedia", texto)


# ===========================================================================
class TestDeteccaoDeFragmento(unittest.TestCase):
    """As tabelas de warhead/linker/e3_ligand têm SMILES, e não são PROTACs.

    Se entrassem, apareceriam na planilha como compostos de 200 Da sem
    atividade — e contariam como "achados" na revisão.
    """

    def test_por_nome_de_arquivo(self):
        for nome in ("warhead.csv", "linker.csv", "e3_ligand.csv",
                     "e3-ligand.csv"):
            mapa, _ = pf.mapear_colunas(["Smiles", "Target"])
            motivo = pf.motivo_para_ignorar(Path(nome), ["Smiles", "Target"], mapa)
            self.assertIsNotNone(motivo, nome)

    def test_por_nome_de_coluna(self):
        colunas = ["Warhead ID", "Warhead Smiles", "Target"]
        mapa, _ = pf.mapear_colunas(colunas)
        self.assertIsNotNone(pf.motivo_para_ignorar(Path("x.csv"), colunas, mapa))

    def test_tabela_de_protac_passa(self):
        colunas = ["Compound ID", "Smiles", "Target", "DC50 (nM)", "Dmax (%)"]
        mapa, _ = pf.mapear_colunas(colunas)
        self.assertIsNone(pf.motivo_para_ignorar(Path("protac.csv"), colunas, mapa))

    def test_smiles_com_contexto_passa_mesmo_sem_atividade(self):
        """O PROTAC-PatentDB quase nunca tem atividade medida, e tem de entrar."""
        colunas = ["compound_id", "SMILES", "Patent Number", "Target", "MW"]
        mapa, _ = pf.mapear_colunas(colunas)
        self.assertIsNone(pf.motivo_para_ignorar(Path("patent.csv"), colunas, mapa))

    def test_so_smiles_sem_nada_nao_passa(self):
        mapa, _ = pf.mapear_colunas(["SMILES"])
        self.assertIsNotNone(pf.motivo_para_ignorar(Path("x.csv"), ["SMILES"], mapa))


# ===========================================================================
class TestFiltro(unittest.TestCase):

    def montar(self) -> pd.DataFrame:
        return pd.DataFrame({
            "alvo": ["BRD4", "AR", "PCSK9", None],
            "uniprot": ["O60885", "P10275", "Q8NBP7", None],
            "nome": ["MZ1", "ARV-110", None, "composto do PCSK9"],
            "titulo": [None, None, None, None],
            "classe_alvo": [None, None, None, None],
            "e3": ["VHL", "CRBN", "VHL", "VHL"],
            "linhagem": ["MV4-11", "VCaP", "HepG2", None],
            "dc50": ["13 nM", "1 nM", None, None],
            "dmax": ["95 %", None, None, None],
        })

    def test_filtro_por_alvo(self):
        r = pf.filtrar(self.montar(), pf.Pedido(alvo="BRD4"))
        self.assertEqual(list(r["alvo"]), ["BRD4"])

    def test_alvo_tambem_no_nome(self):
        """A última linha só diz 'PCSK9' no nome do composto, e tem de entrar."""
        r = pf.filtrar(self.montar(), pf.Pedido(alvo="PCSK9"))
        self.assertEqual(len(r), 2)

    def test_sinonimo(self):
        r = pf.filtrar(self.montar(),
                       pf.Pedido(alvo="bromodomain 4", sinonimos=["BRD4"]))
        self.assertEqual(len(r), 1)

    def test_uniprot(self):
        r = pf.filtrar(self.montar(), pf.Pedido(uniprot="O60885"))
        self.assertEqual(list(r["alvo"]), ["BRD4"])

    def test_uniprot_nao_zera_base_sem_a_coluna(self):
        """Informar o UniProt não pode custar as bases que não o trazem.

        Como critério em E, o UniProt zerava PROTAC-DB, PROTACpedia e
        PROTAC-PatentDB sempre que a coluna faltasse — calado, parecendo que o
        alvo não estava nessas bases. Ele é mais um NOME do alvo, em OU.
        """
        df = self.montar().drop(columns=["uniprot"])
        r = pf.filtrar(df, pf.Pedido(alvo="BRD4", uniprot="O60885"))
        self.assertEqual(len(r), 1)

    def test_uniprot_em_ou_com_o_nome(self):
        """Alvo OU acessão: as duas linhas entram, cada uma por um caminho."""
        df = pd.DataFrame({
            "alvo": ["BRD4", None],
            "uniprot": [None, "O60885"],
            "nome": [None, None], "titulo": [None, None],
            "classe_alvo": [None, None],
        })
        r = pf.filtrar(df, pf.Pedido(alvo="BRD4", uniprot="O60885"))
        self.assertEqual(len(r), 2)

    def test_e3(self):
        r = pf.filtrar(self.montar(), pf.Pedido(e3="CRBN"))
        self.assertEqual(list(r["alvo"]), ["AR"])

    def test_alvo_e_e3_sao_E(self):
        r = pf.filtrar(self.montar(), pf.Pedido(alvo="BRD4", e3="CRBN"))
        self.assertEqual(len(r), 0)

    def test_so_com_degradacao(self):
        r = pf.filtrar(self.montar(), pf.Pedido(so_com_degradacao=True))
        self.assertEqual(len(r), 2)

    def test_sem_criterio_passa_tudo(self):
        r = pf.filtrar(self.montar(), pf.Pedido())
        self.assertEqual(len(r), 4)

    def test_maiuscula_e_acento_nao_importam(self):
        r = pf.filtrar(self.montar(), pf.Pedido(alvo="brd4"))
        self.assertEqual(len(r), 1)


# ===========================================================================
class TestConsolidacao(unittest.TestCase):

    def test_mesmo_composto_em_duas_bases_vira_uma_linha(self):
        df = pd.DataFrame({
            "smiles": ["CC(=O)Nc1ccccc1", "CC(=O)Nc1ccccc1"],
            "dc50": ["13 nM", None],
            "mw": [None, "845.9"],
            "fonte": ["PROTAC-DB", "PROTAC-PatentDB"],
            "id_fonte": ["1", "P-1"],
            "arquivo_origem": ["protac.csv", "patent.csv"],
        })
        r = pf.consolidar(df)
        self.assertEqual(len(r), 1)
        self.assertEqual(r.iloc[0]["dc50"], "13 nM")
        self.assertEqual(r.iloc[0]["mw"], "845.9")
        self.assertIn("PROTAC-DB", r.iloc[0]["fonte"])
        self.assertIn("PROTAC-PatentDB", r.iloc[0]["fonte"])

    def test_ordem_das_bases_e_preferencia(self):
        """A primeira base da lista ganha o valor — medido antes de previsto."""
        df = pd.DataFrame({
            "smiles": ["CCO", "CCO"],
            "mw": ["846 (medido)", "845 (previsto)"],
            "fonte": ["PROTAC-DB", "PROTAC-PatentDB"],
        })
        r = pf.consolidar(df)
        self.assertEqual(len(r), 1)
        self.assertEqual(r.iloc[0]["mw"], "846 (medido)")

    def test_compostos_diferentes_nao_se_juntam(self):
        df = pd.DataFrame({"smiles": ["CCO", "CCN"], "fonte": ["a", "b"]})
        self.assertEqual(len(pf.consolidar(df)), 2)

    def test_linha_sem_smiles_passa_intacta(self):
        """Sem SMILES não há como afirmar que é o mesmo composto."""
        df = pd.DataFrame({"smiles": [None, None], "dc50": ["13 nM", "9 nM"],
                           "fonte": ["a", "b"]})
        self.assertEqual(len(pf.consolidar(df)), 2)


# ===========================================================================
class TestColetaCompleta(unittest.TestCase):
    """Fim a fim, sem rede: arquivos em disco -> tabela -> planilha."""

    def test_coleta_brd4(self):
        with tempfile.TemporaryDirectory() as tmp:
            dados = montar_dir_dados(Path(tmp))
            # O uniprot não é decoração: o PROTAC-Bench não traz o nome do
            # gene, e sem a acessão aquela base sai zerada (ver o teste
            # test_aviso_quando_base_so_tem_uniprot).
            pedido = pf.Pedido(alvo="BRD4", uniprot="O60885",
                               bases=list(pf.FONTES))
            res = pf.coletar(pedido, dir_dados=dados, baixar=False)

            self.assertFalse(res.tabela.empty)
            # BRD4 aparece nas quatro bases; PROTACpedia entra porque o arquivo
            # já está em disco (o cadastro só vale para o download).
            for base in ("PROTAC-DB", "PROTACpedia", "PROTAC-Bench",
                         "PROTAC-PatentDB"):
                self.assertGreater(res.contagens.get(base, 0), 0, base)

            self.assertNotIn("AR", list(res.tabela.get("alvo", [])),
                             "o filtro deixou passar alvo de outra proteína")
            # O MZ1 está em três bases com o mesmo SMILES: tem de ser uma linha.
            mz1 = res.tabela[res.tabela["smiles"] == "CC(=O)Nc1ccccc1"]
            self.assertEqual(len(mz1), 1)
            self.assertIn("PROTAC-DB", str(mz1.iloc[0]["fonte"]))

            self.assertEqual(list(res.tabela["idx"]),
                             list(range(1, len(res.tabela) + 1)))

    def test_tabela_de_fragmento_e_ignorada(self):
        with tempfile.TemporaryDirectory() as tmp:
            dados = Path(tmp) / "bases"
            escrever(dados / "PROTAC-DB", "protac.csv", PROTACDB_CSV)
            escrever(dados / "PROTAC-DB", "warhead.csv", PROTACDB_WARHEAD_CSV)
            res = pf.coletar(pf.Pedido(bases=["PROTAC-DB"]),
                             dir_dados=dados, baixar=False)
            self.assertTrue(any("não parece tabela de composto" in r
                                for r in res.registro))
            self.assertEqual(len(res.tabela), 3)

    def test_base_ausente_nao_derruba_as_outras(self):
        with tempfile.TemporaryDirectory() as tmp:
            dados = Path(tmp) / "bases"
            escrever(dados / "PROTAC-DB", "protac.csv", PROTACDB_CSV)
            res = pf.coletar(pf.Pedido(bases=list(pf.FONTES)),
                             dir_dados=dados, baixar=False)
            self.assertEqual(len(res.tabela), 3)
            self.assertEqual(res.contagens["PROTACpedia"], 0)

    def test_unidades_fixas(self):
        with tempfile.TemporaryDirectory() as tmp:
            dados = montar_dir_dados(Path(tmp))
            res = pf.coletar(pf.Pedido(alvo="BRD4", unidades_fixas=True),
                             dir_dados=dados, baixar=False)
            valores = [v for v in res.tabela["dc50"] if pd.notna(v)]
            self.assertTrue(all(isinstance(v, float) for v in valores), valores)
            # '>10 uM' tem de ter virado 10000 nM
            self.assertIn(10000.0, valores)

    def test_aviso_quando_base_so_tem_uniprot(self):
        """Zero linhas por falta de COLUNA não é zero linhas por falta do alvo.

        Filtrar o PROTAC-Bench por 'BRD4' devolve nada, porque a base só tem a
        acessão. Sem o aviso, isso se leria como "não há PROTAC de BRD4 aqui".
        """
        with tempfile.TemporaryDirectory() as tmp:
            dados = Path(tmp) / "bases"
            escrever(dados / "PROTAC-Bench", "protac_bench.csv", PROTACBENCH_CSV)
            res = pf.coletar(pf.Pedido(alvo="BRD4", bases=["PROTAC-Bench"]),
                             dir_dados=dados, baixar=False)
            self.assertTrue(res.tabela.empty)
            self.assertTrue(any("não identifica o alvo por nome" in r
                                for r in res.registro), res.registro)

    def test_uniprot_resgata_a_base_sem_nome_de_alvo(self):
        with tempfile.TemporaryDirectory() as tmp:
            dados = Path(tmp) / "bases"
            escrever(dados / "PROTAC-Bench", "protac_bench.csv", PROTACBENCH_CSV)
            res = pf.coletar(
                pf.Pedido(alvo="BRD4", uniprot="O60885", bases=["PROTAC-Bench"]),
                dir_dados=dados, baixar=False)
            self.assertEqual(len(res.tabela), 1)

    def test_relatorio_menciona_coluna_sem_destino(self):
        with tempfile.TemporaryDirectory() as tmp:
            dados = Path(tmp) / "bases"
            escrever(dados / "PROTAC-DB", "protac.csv",
                     "SMILES,DC50 (nM),Ternary Complex Kd\nCCO,13,45\n")
            res = pf.coletar(pf.Pedido(bases=["PROTAC-DB"]),
                             dir_dados=dados, baixar=False)
            self.assertIn("Ternary Complex Kd", res.relatorio())


# ===========================================================================
class TestExportacao(unittest.TestCase):

    def tabela(self) -> pd.DataFrame:
        with tempfile.TemporaryDirectory() as tmp:
            dados = montar_dir_dados(Path(tmp))
            return pf.coletar(pf.Pedido(alvo="BRD4"), dir_dados=dados,
                              baixar=False)

    def test_aba_papers_bate_com_a_planilha(self):
        from openpyxl import load_workbook
        res = self.tabela()
        with tempfile.TemporaryDirectory() as tmp:
            destino = Path(tmp) / "saida.xlsx"
            pf.exportar_planilha(res.tabela, destino, resultado=res)
            ws = load_workbook(destino)["PAPERS"]

            cabecalho = [ws.cell(row=3, column=c).value
                         for c in range(1, len(pf.CHAVES) + 1)]
            self.assertEqual(cabecalho, [pf.ROTULOS[k] for k in pf.CHAVES])
            self.assertEqual(ws.cell(row=4, column=1).value, 1)
            self.assertIsNotNone(ws.cell(row=4, column=4).value)  # SMILES

    def test_abas_de_auditoria_existem(self):
        from openpyxl import load_workbook
        res = self.tabela()
        with tempfile.TemporaryDirectory() as tmp:
            destino = Path(tmp) / "saida.xlsx"
            pf.exportar_planilha(res.tabela, destino, resultado=res)
            wb = load_workbook(destino)
            self.assertIn("PROVENIENCIA", wb.sheetnames)
            self.assertIn("MAPEAMENTO", wb.sheetnames)
            self.assertGreater(wb["MAPEAMENTO"].max_row, 1)

    def test_modelo_e_preservado_e_nunca_sobrescrito(self):
        from openpyxl import Workbook, load_workbook
        res = self.tabela()
        with tempfile.TemporaryDirectory() as tmp:
            modelo = Path(tmp) / "modelo.xlsx"
            wb = Workbook()
            ws = wb.active
            ws.title = "PAPERS"
            for c, chave in enumerate(pf.CHAVES, start=1):
                ws.cell(row=3, column=c, value=pf.ROTULOS[chave])
            wb.create_sheet("CRONORAMA")
            wb.save(modelo)
            antes = modelo.read_bytes()

            destino = Path(tmp) / "saida.xlsx"
            pf.exportar_planilha(res.tabela, destino, modelo=modelo, resultado=res)

            self.assertEqual(modelo.read_bytes(), antes, "o modelo foi alterado")
            wb2 = load_workbook(destino)
            self.assertIn("CRONORAMA", wb2.sheetnames, "abas do modelo sumiram")
            self.assertEqual(wb2["PAPERS"].cell(row=4, column=1).value, 1)

    def test_recusa_escrever_por_cima_do_modelo(self):
        from openpyxl import Workbook
        res = self.tabela()
        with tempfile.TemporaryDirectory() as tmp:
            modelo = Path(tmp) / "modelo.xlsx"
            wb = Workbook()
            wb.active.title = "PAPERS"
            wb.save(modelo)
            with self.assertRaises(ValueError):
                pf.exportar_planilha(res.tabela, modelo, modelo=modelo)

    def test_csv_sai_com_os_rotulos_da_planilha(self):
        res = self.tabela()
        with tempfile.TemporaryDirectory() as tmp:
            destino = Path(tmp) / "saida.csv"
            pf.exportar_csv(res.tabela, destino)
            cabecalho = destino.read_text(encoding="utf-8").splitlines()[0]
            self.assertIn("DC₅₀", cabecalho)
            self.assertIn("Comprimeto do Linker", cabecalho)

    def test_celula_vazia_e_vazia_mesmo(self):
        """Nada de 'nan' escrito como texto na planilha."""
        from openpyxl import load_workbook
        res = self.tabela()
        with tempfile.TemporaryDirectory() as tmp:
            destino = Path(tmp) / "saida.xlsx"
            pf.exportar_planilha(res.tabela, destino)
            ws = load_workbook(destino)["PAPERS"]
            for linha in ws.iter_rows(min_row=4, max_row=4 + len(res.tabela)):
                for celula in linha:
                    if isinstance(celula.value, str):
                        self.assertNotIn(celula.value.strip().lower(),
                                         {"nan", "<na>", "none"})


if __name__ == "__main__":
    unittest.main(verbosity=2)
