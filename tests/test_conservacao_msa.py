#!/usr/bin/env python
"""
Testes de msa_henikoff.py, bolso_ligante.py e conservacao_bolso.py.

Só unittest da biblioteca padrão — não exige pytest, para rodar em qualquer
env do projeto (`mdtools`, `pf_vs`, ou o python do sistema).

    python -m unittest discover -s tests -v
    python tests/test_conservacao_msa.py

O que está coberto, e por que estes casos
-----------------------------------------
- Os pesos de Henikoff são checados contra o EXEMPLO PUBLICADO (Henikoff &
  Henikoff 1994): 4 sequências, 5 colunas, pesos 0,267/0,267/0,200/0,267. É o
  único jeito honesto de validar a implementação.
- O bolso é checado contra uma estrutura sintética com distâncias EXATAS
  construídas para ter resposta conhecida, inclusive um átomo a 5,1 Å que tem
  de ficar fora e um a 4,99 Å que tem de entrar.
- Os dois caminhos de cálculo (numpy e python puro) são comparados entre si:
  o script roda nos dois, então os dois precisam concordar.
- Os formatos (FASTA alinhado, a3m, Stockholm, CSV) são comparados ao MESMO
  alinhamento escrito de quatro formas: todos têm de dar o mesmo resultado.
- Os erros de entrada são testados pela MENSAGEM, porque a mensagem é o
  produto: um alinhamento não-retangular tem de dizer qual sequência está
  fora, não "IndexError".
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ / "scripts"))

import bolso_ligante as bl               # noqa: E402
import conservacao_bolso as cb           # noqa: E402
import msa_henikoff as mh                # noqa: E402

PY = sys.executable


def escrever(dirn: Path, nome: str, texto: str) -> Path:
    p = dirn / nome
    p.write_text(texto)
    return p


def fasta(registros) -> str:
    return "".join(f">{i}\n{s}\n" for i, s in registros)


# o exemplo do artigo de Henikoff & Henikoff (1994)
EXEMPLO_HH = [("s1", "GYVGS"), ("s2", "GFDGF"), ("s3", "GYDGF"), ("s4", "GYQGG")]
PESOS_HH = [1.3333333, 1.3333333, 1.0, 1.3333333]     # antes de normalizar


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()


# ===========================================================================
# leitura de MSA


class TesteLeituraMSA(TempDirCase):
    def test_fasta_alinhado(self):
        p = escrever(self.dir, "a.fasta", fasta(EXEMPLO_HH))
        msa = mh.ler_msa(p)
        self.assertEqual(msa.formato, "fasta")
        self.assertEqual(msa.n, 4)
        self.assertEqual(msa.largura, 5)
        self.assertEqual(msa.query, "GYVGS")
        self.assertEqual(msa.ids[0], "s1")

    def test_fasta_multilinha_e_minusculas(self):
        """FASTA quebrado em linhas e em minúsculas é o mesmo alinhamento."""
        p = escrever(self.dir, "b.fasta", ">x\nGYV\nGS\n>y\ngfdgf\n")
        msa = mh.ler_msa(p)
        self.assertEqual(msa.seqs, ["GYVGS", "GFDGF"])

    def test_a3m_remove_insercoes(self):
        """Minúscula em a3m é inserção: sai, e o alinhamento volta a fechar."""
        p = escrever(self.dir, "c.a3m",
                     ">query\nGYVGS\n>h1\nGYmmDGS\n>h2\nGF-GS\n")
        msa = mh.ler_msa(p)
        self.assertEqual(msa.formato, "a3m")
        self.assertEqual(msa.seqs, ["GYVGS", "GYDGS", "GF-GS"])
        self.assertTrue(any("inserções" in a for a in msa.avisos))

    def test_a2m_tratado_como_a3m(self):
        p = escrever(self.dir, "c.a2m", ">q\nGYVGS\n>h\nGYvVGS\n")
        self.assertEqual(mh.ler_msa(p).seqs, ["GYVGS", "GYVGS"])

    def test_stockholm_em_blocos(self):
        """jackhmmer quebra a sequência em blocos; juntar é obrigação nossa."""
        texto = (
            "# STOCKHOLM 1.0\n"
            "#=GF ID teste\n"
            "s1  GYV\n"
            "s2  GFD\n"
            "\n"
            "s1  GS\n"
            "s2  GF\n"
            "//\n"
        )
        p = escrever(self.dir, "d.sto", texto)
        msa = mh.ler_msa(p)
        self.assertEqual(msa.formato, "stockholm")
        self.assertEqual(msa.seqs, ["GYVGS", "GFDGF"])

    def test_stockholm_ponto_e_gap_equivalentes(self):
        p = escrever(self.dir, "e.sto", "# STOCKHOLM 1.0\ns1 GY.GS\ns2 GF-GF\n//\n")
        msa = mh.ler_msa(p)
        self.assertEqual(msa.seqs, ["GY-GS", "GF-GF"])

    def test_csv_do_boltz(self):
        p = escrever(self.dir, "f.csv", "key,sequence\nquery,GYVGS\nh1,GFDGF\n")
        msa = mh.ler_msa(p)
        self.assertEqual(msa.formato, "csv")
        self.assertEqual(msa.seqs, ["GYVGS", "GFDGF"])
        self.assertEqual(msa.ids, ["query", "h1"])

    def test_tsv_tambem(self):
        p = escrever(self.dir, "g.tsv", "id\tseq\nq\tGYVGS\nh\tGFDGF\n")
        self.assertEqual(mh.ler_msa(p).seqs, ["GYVGS", "GFDGF"])

    def test_quatro_formatos_dao_o_mesmo_perfil(self):
        """FASTA, a3m, Stockholm e CSV do mesmo alinhamento = mesmo resultado."""
        regs = [("q", "GYVGS"), ("h1", "GFDGF"), ("h2", "GYDGF"), ("h3", "GYQGG")]
        arquivos = [
            escrever(self.dir, "m.fasta", fasta(regs)),
            escrever(self.dir, "m.a3m", fasta(regs)),
            escrever(self.dir, "m.sto",
                     "# STOCKHOLM 1.0\n" + "".join(f"{i} {s}\n" for i, s in regs) + "//\n"),
            escrever(self.dir, "m.csv",
                     "key,sequence\n" + "".join(f"{i},{s}\n" for i, s in regs)),
        ]
        refs = [mh.analisar_msa(a)["colunas"] for a in arquivos]
        for outro in refs[1:]:
            for c1, c2 in zip(refs[0], outro):
                self.assertAlmostEqual(c1["jsd"], c2["jsd"], places=10)
                self.assertAlmostEqual(c1["entropia_bits"], c2["entropia_bits"],
                                       places=10)

    def test_extensao_desconhecida_usa_conteudo(self):
        p = escrever(self.dir, "sem_extensao", ">q\nGYVGS\n>h\nGFDGF\n")
        self.assertEqual(mh.ler_msa(p).n, 2)

    def test_caracteres_fora_do_alfabeto_viram_x(self):
        p = escrever(self.dir, "h.fasta", ">q\nGYVGS\n>h\nGBZGJ\n")
        msa = mh.ler_msa(p)
        self.assertEqual(msa.seqs[1], "GXXGX")

    def test_amostragem_mantem_a_query_e_e_deterministica(self):
        regs = [("q", "AAAAA")] + [(f"h{i}", "CCCCC") for i in range(50)]
        p = escrever(self.dir, "i.fasta", fasta(regs))
        a = mh.ler_msa(p, max_seqs=10, semente=7)
        b = mh.ler_msa(p, max_seqs=10, semente=7)
        self.assertEqual(a.n, 10)
        self.assertEqual(a.ids[0], "q")
        self.assertEqual(a.ids, b.ids)
        self.assertTrue(any("amostradas" in x for x in a.avisos))

    def test_ids_repetidos_ganham_sufixo(self):
        p = escrever(self.dir, "j.fasta", ">x\nGYVGS\n>x\nGFDGF\n")
        msa = mh.ler_msa(p)
        self.assertEqual(msa.ids, ["x", "x__1"])

    def test_sequencia_so_de_gap_e_descartada(self):
        p = escrever(self.dir, "k.fasta", ">q\nGYVGS\n>vazia\n-----\n>h\nGFDGF\n")
        msa = mh.ler_msa(p)
        self.assertEqual(msa.n, 2)
        self.assertEqual(msa.descartadas, ["vazia"])


class TesteErrosDeLeitura(TempDirCase):
    def test_arquivo_inexistente(self):
        with self.assertRaisesRegex(ValueError, "não encontrado"):
            mh.ler_msa(self.dir / "nao_existe.fasta")

    def test_arquivo_vazio(self):
        p = escrever(self.dir, "v.fasta", "")
        with self.assertRaisesRegex(ValueError, "vazio"):
            mh.ler_msa(p)

    def test_diretorio_em_vez_de_arquivo(self):
        with self.assertRaisesRegex(ValueError, "diretório"):
            mh.ler_msa(self.dir)

    def test_alinhamento_nao_retangular_diz_quem_esta_fora(self):
        p = escrever(self.dir, "r.fasta", ">q\nGYVGS\n>curta\nGYV\n>ok\nGFDGF\n")
        with self.assertRaises(ValueError) as ctx:
            mh.ler_msa(p)
        msg = str(ctx.exception)
        self.assertIn("não é retangular", msg)
        self.assertIn("curta", msg)
        self.assertIn("a3m", msg)             # sugere a causa provável

    def test_sem_cabecalho_fasta(self):
        p = escrever(self.dir, "s.fasta", "GYVGS\nGFDGF\n")
        with self.assertRaisesRegex(ValueError, "não é FASTA"):
            mh.ler_msa(p)

    def test_csv_sem_coluna_de_sequencia(self):
        p = escrever(self.dir, "t.csv", "a,b\n1,2\n")
        with self.assertRaisesRegex(ValueError, "coluna de sequência"):
            mh.ler_msa(p)

    def test_parquet_sem_engine_da_erro_tratado(self):
        """Sem pyarrow, o erro tem de ser ValueError com saída — não traceback.

        É o caso do .aligned.pqt do Chai num env que não tem pyarrow.
        """
        p = self.dir / "chai.aligned.pqt"
        p.write_bytes(b"PAR1 isto nao e um parquet valido")
        with self.assertRaises(ValueError) as ctx:
            mh.ler_msa(p)
        msg = str(ctx.exception)
        self.assertIn("parquet", msg)
        self.assertIn("FASTA alinhado", msg)       # diz o que fazer

    def test_sequencia_unica_e_erro_por_padrao(self):
        """ESMFold é single-sequence: o script tem de dizer isso, não calcular."""
        p = escrever(self.dir, "u.fasta", ">só_a_query\nGYVGS\n")
        with self.assertRaises(ValueError) as ctx:
            mh.analisar_msa(p)
        self.assertIn("UMA sequência", str(ctx.exception))
        self.assertIn("ESMFold", str(ctx.exception))

    def test_sequencia_unica_permitida_vira_aviso(self):
        p = escrever(self.dir, "u2.fasta", ">q\nGYVGS\n")
        r = mh.analisar_msa(p, permitir_unica=True)
        self.assertEqual(r["resumo"]["n_sequencias"], 1)
        self.assertIn("single-sequence", r["resumo"]["veredito"])
        self.assertEqual(r["pesos"], [1.0])

    def test_msa_raso_avisa(self):
        p = escrever(self.dir, "w.fasta", fasta(EXEMPLO_HH))
        r = mh.analisar_msa(p)
        self.assertTrue(any("raso" in a for a in r["resumo"]["avisos"]))


# ===========================================================================
# pesos de Henikoff


class TestePesosHenikoff(TempDirCase):
    def msa_exemplo(self):
        return mh.ler_msa(escrever(self.dir, "hh.fasta", fasta(EXEMPLO_HH)))

    def test_exemplo_publicado_soma(self):
        """Henikoff & Henikoff 1994: 1,333 / 1,333 / 1,000 / 1,333."""
        msa = self.msa_exemplo()
        pesos, info = mh.pesos_henikoff(msa, normalizacao="soma")
        esperado = [w * 4 / sum(PESOS_HH) for w in PESOS_HH]
        for obtido, esp in zip(pesos, esperado):
            self.assertAlmostEqual(obtido, esp, places=6)
        # normalizados pelo artigo: 0,267 0,267 0,200 0,267
        fracoes = [w / 4 for w in pesos]
        for obtido, esp in zip(fracoes, [0.26667, 0.26667, 0.2, 0.26667]):
            self.assertAlmostEqual(obtido, esp, places=4)

    def test_pesos_somam_n(self):
        msa = self.msa_exemplo()
        for norm in ("media", "soma"):
            pesos, _ = mh.pesos_henikoff(msa, normalizacao=norm)
            self.assertAlmostEqual(sum(pesos), msa.n, places=9)

    def test_media_e_soma_coincidem_sem_gap(self):
        """Sem gap toda sequência cobre todas as colunas: as duas normalizações
        têm de dar o mesmo (é o controle de que 'media' não distorce nada)."""
        msa = self.msa_exemplo()
        a, _ = mh.pesos_henikoff(msa, normalizacao="media")
        b, _ = mh.pesos_henikoff(msa, normalizacao="soma")
        for x, y in zip(a, b):
            self.assertAlmostEqual(x, y, places=9)

    def test_duplicata_divide_o_peso(self):
        """Duas cópias idênticas valem, juntas, o que uma valia."""
        regs = [("a", "ACDEF"), ("b", "ACDEF"), ("c", "WWWWW")]
        msa = mh.ler_msa(escrever(self.dir, "dup.fasta", fasta(regs)))
        pesos, _ = mh.pesos_henikoff(msa)
        self.assertAlmostEqual(pesos[0], pesos[1], places=9)
        self.assertGreater(pesos[2], pesos[0])
        self.assertAlmostEqual(pesos[0] + pesos[1], pesos[2], places=6)

    def test_redundancia_derruba_neff(self):
        """20 clones + 1 diferente: Neff tem de ficar MUITO abaixo de 21."""
        regs = [(f"clone{i}", "ACDEFGHIKL") for i in range(20)]
        regs.append(("distinta", "WWWWWWWWWW"))
        msa = mh.ler_msa(escrever(self.dir, "red.fasta", fasta(regs)))
        _, info = mh.pesos_henikoff(msa)
        self.assertLess(info["neff"], 6.0)
        self.assertGreater(info["redundancia"], 3.0)

    def test_sem_redundancia_neff_proximo_de_n(self):
        regs = [(f"s{i}", "".join(mh.ALFABETO[(i + k) % 20] for k in range(10)))
                for i in range(20)]
        msa = mh.ler_msa(escrever(self.dir, "div.fasta", fasta(regs)))
        _, info = mh.pesos_henikoff(msa)
        self.assertGreater(info["neff"], 19.0)

    def test_gap_nao_e_tipo_por_padrao(self):
        """Com gap como tipo, a sequência truncada ganha peso; sem, não."""
        regs = [("cheia1", "ACDEFGHIKL"), ("cheia2", "ACDEFGHIKL"),
                ("truncada", "ACDEF-----")]
        msa = mh.ler_msa(escrever(self.dir, "gap.fasta", fasta(regs)))
        sem, _ = mh.pesos_henikoff(msa, gap_como_tipo=False)
        com, _ = mh.pesos_henikoff(msa, gap_como_tipo=True)
        self.assertGreater(com[2], sem[2])

    def test_normalizacao_media_protege_parcial(self):
        """Um fragmento curto não deve ser punido por ser curto."""
        regs = [("cheia", "ACDEFGHIKL"), ("outra", "WYVTSRQPNM"),
                ("fragmento", "ACDE------")]
        msa = mh.ler_msa(escrever(self.dir, "par.fasta", fasta(regs)))
        media, _ = mh.pesos_henikoff(msa, normalizacao="media")
        soma, _ = mh.pesos_henikoff(msa, normalizacao="soma")
        self.assertGreater(media[2], soma[2])

    def test_so_colunas_query(self):
        regs = [("q", "AC--EF"), ("h1", "ACDDEF"), ("h2", "WCDDEW")]
        msa = mh.ler_msa(escrever(self.dir, "sq.fasta", fasta(regs)))
        _, info = mh.pesos_henikoff(msa, so_colunas_query=True)
        self.assertEqual(info["colunas_usadas"], 4)
        _, info2 = mh.pesos_henikoff(msa, so_colunas_query=False)
        self.assertEqual(info2["colunas_usadas"], 6)

    def test_coluna_toda_gap_nao_quebra(self):
        regs = [("q", "AC-EF"), ("h", "WC-EW")]
        msa = mh.ler_msa(escrever(self.dir, "tg.fasta", fasta(regs)))
        pesos, info = mh.pesos_henikoff(msa)
        self.assertAlmostEqual(sum(pesos), 2.0, places=9)
        self.assertFalse(info["uniforme"])

    def test_numpy_e_python_concordam(self):
        """O script usa numpy quando existe; os dois caminhos têm de bater."""
        regs = [("q", "ACDEFGHIKLMNPQRSTVWY")]
        for i in range(30):
            s = list("ACDEFGHIKLMNPQRSTVWY")
            for k in range(0, 20, (i % 5) + 2):
                s[k] = mh.ALFABETO[(i + k) % 20]
            if i % 4 == 0:
                s[0] = s[1] = mh.GAP
            regs.append((f"h{i}", "".join(s)))
        msa = mh.ler_msa(escrever(self.dir, "np.fasta", fasta(regs)))
        p_np, i_np = mh.pesos_henikoff(msa, usar_numpy=True)
        p_py, i_py = mh.pesos_henikoff(msa, usar_numpy=False)
        for a, b in zip(p_np, p_py):
            self.assertAlmostEqual(a, b, places=10)
        self.assertAlmostEqual(i_np["neff"], i_py["neff"], places=9)
        c_np = mh.perfil_colunas(msa, p_np, usar_numpy=True)
        c_py = mh.perfil_colunas(msa, p_py, usar_numpy=False)
        for x, y in zip(c_np, c_py):
            for campo in ("jsd", "entropia_bits", "ocupacao_pond", "neff_coluna",
                          "kl_bits", "conservacao_shannon", "frac_gap_pond"):
                self.assertAlmostEqual(x[campo], y[campo], places=10,
                                       msg=f"coluna {x['coluna']} campo {campo}")

    def test_normalizacao_invalida(self):
        msa = self.msa_exemplo()
        with self.assertRaisesRegex(ValueError, "media.*soma"):
            mh.pesos_henikoff(msa, normalizacao="qualquer")


# ===========================================================================
# escores por coluna


class TesteEscores(TempDirCase):
    def test_coluna_conservada_vs_variavel(self):
        regs = [(f"s{i}", "W" + mh.ALFABETO[i % 20]) for i in range(20)]
        msa = mh.ler_msa(escrever(self.dir, "cons.fasta", fasta(regs)))
        pesos, _ = mh.pesos_henikoff(msa)
        cols = mh.perfil_colunas(msa, pesos)
        conservada, variavel = cols[0], cols[1]
        self.assertAlmostEqual(conservada["entropia_bits"], 0.0, places=9)
        self.assertAlmostEqual(conservada["conservacao_shannon"], 1.0, places=9)
        self.assertAlmostEqual(variavel["entropia_bits"], math.log2(20), places=6)
        self.assertAlmostEqual(variavel["conservacao_shannon"], 0.0, places=6)
        self.assertGreater(conservada["jsd"], variavel["jsd"])
        self.assertAlmostEqual(variavel["jsd"], 0.0, delta=0.05)

    def test_jsd_limitada_entre_0_e_1(self):
        for aa in mh.ALFABETO:
            p = {aa: 1.0}
            d = mh.jsd(p)
            self.assertGreaterEqual(d, 0.0)
            self.assertLessEqual(d, 1.0)
        self.assertAlmostEqual(mh.jsd(dict(mh.FUNDO)), 0.0, places=9)

    def test_jsd_premia_residuo_raro(self):
        """W conservado informa mais que L conservado: W é raro no fundo."""
        self.assertGreater(mh.jsd({"W": 1.0}), mh.jsd({"L": 1.0}))

    def test_entropia_e_kl_coerentes(self):
        p = {"A": 0.5, "C": 0.5}
        self.assertAlmostEqual(mh.entropia_bits(p), 1.0, places=9)
        self.assertGreater(mh.kl_bits(p), 0.0)

    def test_ocupacao_e_neff_de_coluna(self):
        regs = [("q", "AA"), ("h1", "A-"), ("h2", "A-"), ("h3", "A-")]
        msa = mh.ler_msa(escrever(self.dir, "oc.fasta", fasta(regs)))
        pesos, _ = mh.pesos_henikoff(msa)
        cols = mh.perfil_colunas(msa, pesos, min_ocupacao=0.5, min_neff_coluna=3.0)
        self.assertAlmostEqual(cols[0]["ocupacao"], 1.0, places=9)
        self.assertAlmostEqual(cols[1]["ocupacao"], 0.25, places=9)
        self.assertEqual(cols[1]["profundidade"], 1)
        self.assertEqual(cols[1]["confiavel"], 0)       # rasa: não vale
        self.assertLess(cols[1]["neff_coluna"], 2.0)

    def test_coluna_rasa_nao_e_confiavel_mesmo_conservada(self):
        regs = [("q", "AA"), ("h1", "-A"), ("h2", "-A"), ("h3", "-A"), ("h4", "-A")]
        msa = mh.ler_msa(escrever(self.dir, "rasa.fasta", fasta(regs)))
        pesos, _ = mh.pesos_henikoff(msa)
        cols = mh.perfil_colunas(msa, pesos)
        self.assertAlmostEqual(cols[0]["conservacao_shannon"], 1.0, places=9)
        self.assertEqual(cols[0]["confiavel"], 0)
        self.assertEqual(cols[1]["confiavel"], 1)

    def test_janela_movel(self):
        v = [0.0, 0.0, 1.0, 0.0, 0.0]
        s = mh.janela_movel(v, meia_janela=1, peso_vizinhos=0.5)
        self.assertAlmostEqual(s[2], 0.5, places=9)      # metade de si + vizinhos 0
        self.assertGreater(s[1], 0.0)                    # vazou do vizinho
        self.assertEqual(mh.janela_movel(v, 0, 0.5), v)

    def test_pos_query_pula_gap_da_query(self):
        regs = [("q", "A-C"), ("h", "AWC")]
        msa = mh.ler_msa(escrever(self.dir, "pq.fasta", fasta(regs)))
        pesos, _ = mh.pesos_henikoff(msa)
        cols = mh.perfil_colunas(msa, pesos)
        self.assertEqual([c["pos_query"] for c in cols], [1, "", 2])

    def test_resumo_tem_veredito_e_efeito_dos_pesos(self):
        regs = [("q", "ACDEFGHIKL")] + [(f"c{i}", "ACDEFGHIKL") for i in range(30)]
        regs.append(("dif", "WWWWWWWWWW"))
        p = escrever(self.dir, "res.fasta", fasta(regs))
        r = mh.analisar_msa(p)
        s = r["resumo"]
        self.assertIn("veredito", s)
        self.assertIn("efeito_dos_pesos_jsd", s)
        self.assertIn("jsd_sem_peso", r["colunas"][0])
        # com 30 clones, a conservação sem peso é maior que a ponderada
        self.assertLess(s["jsd_media_colunas_query"], s["jsd_media_sem_peso"])

    def test_identidade(self):
        self.assertAlmostEqual(mh.identidade("ACDE", "ACDE"), 1.0)
        self.assertAlmostEqual(mh.identidade("ACDE", "ACDW"), 0.75)
        self.assertAlmostEqual(mh.identidade("AC--", "ACWW"), 1.0)
        self.assertAlmostEqual(mh.identidade("----", "ACWW"), 0.0)


# ===========================================================================
# estrutura: PDB, mmCIF, categorias


#   O ligante LIG tem dois átomos: C1 em (0,0,0) e C2 em (1,0,0). Todas as
#   distâncias abaixo são AO MAIS PRÓXIMO DOS DOIS — é o que o script mede, e
#   foi o erro da primeira versão deste arquivo medir só a partir do C1.
#   TRP3 fica a 4,99 Å do C2 (entra) e LEU4 a 5,10 Å (fica fora): o par que
#   testa o corte exato em 5 Å.
PDB_MINIMO = """\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 00.00           C
ATOM      2  CA  GLY A   2       3.000   0.000   0.000  1.00 00.00           C
ATOM      3  CA  TRP A   3       5.990   0.000   0.000  1.00 00.00           C
ATOM      4  CA  LEU A   4       6.100   0.000   0.000  1.00 00.00           C
ATOM      5  CA  MET A   5      20.000   0.000   0.000  1.00 00.00           C
ATOM      6  CA  CYS B   1       0.000   4.000   0.000  1.00 00.00           C
HETATM    7 ZN    ZN A 101       0.000  -4.000   0.000  1.00 00.00          ZN
HETATM    8  O   HOH W   1       2.000   0.000   0.000  1.00 00.00           O
HETATM    9  C1  GOL A 102       0.000   0.000   4.500  1.00 00.00           C
HETATM   10  FE  HEM A 103       0.000   0.000  -4.000  1.00 00.00          FE
HETATM   11  C1  NAG A 104      30.000   0.000   0.000  1.00 00.00           C
HETATM   12  C1  LIG L   1       0.000   0.000   0.000  1.00 00.00           C
HETATM   13  C2  LIG L   1       1.000   0.000   0.000  1.00 00.00           C
HETATM   14  H1  LIG L   1       0.500   0.500   0.500  1.00 00.00           H
END
"""


class TesteLeituraEstrutura(TempDirCase):
    def test_pdb_campos(self):
        p = escrever(self.dir, "m.pdb", PDB_MINIMO)
        atomos = bl.ler_estrutura(p)
        self.assertEqual(len(atomos), 14)
        a = atomos[0]
        self.assertEqual((a.nome, a.resname, a.cadeia, a.resnum), ("CA", "ALA", "A", 1))
        self.assertAlmostEqual(a.x, 0.0)
        zn = next(x for x in atomos if x.resname == "ZN")
        self.assertEqual(zn.elemento, "ZN")
        self.assertTrue(zn.hetero)

    def test_pdb_sem_coluna_de_elemento(self):
        """PDB velho/saída de ferramenta sem as colunas 77-78."""
        texto = (
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00\n"
            "HETATM    2 ZN    ZN A 101       1.000   0.000   0.000  1.00  0.00\n"
        )
        p = escrever(self.dir, "velho.pdb", texto)
        atomos = bl.ler_estrutura(p)
        self.assertEqual(atomos[0].elemento, "C")
        self.assertEqual(atomos[1].elemento, "ZN")

    def test_altloc_secundario_ignorado(self):
        texto = (
            "ATOM      1  CA AALA A   1       0.000   0.000   0.000  1.00 0.00           C\n"
            "ATOM      2  CA BALA A   1       9.000   0.000   0.000  1.00 0.00           C\n"
        )
        p = escrever(self.dir, "alt.pdb", texto)
        self.assertEqual(len(bl.ler_estrutura(p)), 1)

    def test_modelos_nmr(self):
        texto = (
            "MODEL        1\n"
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 0.00           C\n"
            "ENDMDL\n"
            "MODEL        2\n"
            "ATOM      1  CA  ALA A   1       9.000   0.000   0.000  1.00 0.00           C\n"
            "ENDMDL\n"
        )
        p = escrever(self.dir, "nmr.pdb", texto)
        self.assertAlmostEqual(bl.ler_estrutura(p, modelo=1)[0].x, 0.0)
        self.assertAlmostEqual(bl.ler_estrutura(p, modelo=2)[0].x, 9.0)
        with self.assertRaisesRegex(ValueError, "nenhum átomo"):
            bl.ler_estrutura(p, modelo=7)

    def test_cif_equivale_ao_pdb(self):
        """O mesmo conteúdo em mmCIF tem de dar os MESMOS contatos.

        É o teste que protege o parser de mmCIF: Boltz-2 e Chai-1/2 só entregam
        .cif, então um parser torto aqui invalidaria a análise inteira sem dar
        erro nenhum.
        """
        pdb = escrever(self.dir, "x.pdb", PDB_MINIMO)
        atomos_pdb = bl.ler_estrutura(pdb)
        linhas = [
            "data_teste", "#", "loop_",
            "_atom_site.group_PDB", "_atom_site.id", "_atom_site.type_symbol",
            "_atom_site.label_atom_id", "_atom_site.label_alt_id",
            "_atom_site.label_comp_id", "_atom_site.label_asym_id",
            "_atom_site.label_seq_id", "_atom_site.pdbx_PDB_ins_code",
            "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z",
            "_atom_site.occupancy", "_atom_site.B_iso_or_equiv",
            "_atom_site.auth_seq_id", "_atom_site.auth_asym_id",
            "_atom_site.pdbx_PDB_model_num",
        ]
        for a in atomos_pdb:
            linhas.append(
                f"{'HETATM' if a.hetero else 'ATOM'} {a.serial} {a.elemento} "
                f"{a.nome} . {a.resname} {a.cadeia} {a.resnum} ? "
                f"{a.x:.3f} {a.y:.3f} {a.z:.3f} 1.00 0.00 {a.resnum} {a.cadeia} 1"
            )
        cif = escrever(self.dir, "x.cif", "\n".join(linhas) + "\n")
        atomos_cif = bl.ler_estrutura(cif)
        self.assertEqual(len(atomos_cif), len(atomos_pdb))
        for a, b in zip(atomos_pdb, atomos_cif):
            self.assertEqual((a.nome, a.resname, a.cadeia, a.resnum, a.hetero),
                             (b.nome, b.resname, b.cadeia, b.resnum, b.hetero))
            self.assertAlmostEqual(a.x, b.x, places=3)
        lig_p, _, _ = bl.selecionar_ligante(atomos_pdb, resname="LIG")
        lig_c, _, _ = bl.selecionar_ligante(atomos_cif, resname="LIG")
        c_p, _ = bl.contatos(atomos_pdb, lig_p, raio=5.0)
        c_c, _ = bl.contatos(atomos_cif, lig_c, raio=5.0)
        self.assertEqual([x["residuo"] for x in c_p], [x["residuo"] for x in c_c])

    def test_cif_com_aspas_e_campos_ausentes(self):
        linhas = [
            "data_t", "loop_",
            "_atom_site.group_PDB", "_atom_site.id", "_atom_site.type_symbol",
            "_atom_site.label_atom_id", "_atom_site.label_comp_id",
            "_atom_site.auth_asym_id", "_atom_site.auth_seq_id",
            "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z",
            "ATOM 1 C \"C1'\" ALA A 1 0.000 0.000 0.000",
            "ATOM 2 C 'CA' GLY A 2 1.000 0.000 0.000",
        ]
        p = escrever(self.dir, "q.cif", "\n".join(linhas) + "\n")
        atomos = bl.ler_estrutura(p)
        self.assertEqual([a.nome for a in atomos], ["C1'", "CA"])

    def test_estrutura_inexistente_e_vazia(self):
        with self.assertRaisesRegex(ValueError, "não encontrada"):
            bl.ler_estrutura(self.dir / "nada.pdb")
        p = escrever(self.dir, "z.pdb", "")
        with self.assertRaisesRegex(ValueError, "vazia"):
            bl.ler_estrutura(p)


class TesteCategorias(unittest.TestCase):
    def test_categorias_principais(self):
        casos = {
            "ALA": "proteina", "MSE": "proteina", "HYP": "proteina",
            "ACE": "proteina",
            "DA": "acido_nucleico", "U": "acido_nucleico",
            "HOH": "agua", "WAT": "agua",
            "ZN": "ion", "MG": "ion", "CL": "ion", "SO4": "ion",
            "HEM": "cofator", "NAD": "cofator", "FAD": "cofator",
            "SAM": "cofator", "ATP": "cofator", "SF4": "cofator",
            "GOL": "aditivo", "EDO": "aditivo", "PEG": "aditivo",
            "NAG": "glicano", "MAN": "glicano",
            "LDA": "lipideo", "CLR": "lipideo",
            "LIG": "hetero_desconhecido", "UNK": "hetero_desconhecido",
            "063": "hetero_desconhecido",
        }
        for resname, esperado in casos.items():
            self.assertEqual(bl.categorizar(resname)[0], esperado, msg=resname)

    def test_ion_reconhecido_pelo_elemento(self):
        """Resíduo de 1 átomo com elemento metálico é íon, mesmo com nome exótico."""
        self.assertEqual(bl.categorizar("MN1", 1, ["MN"])[0], "ion")
        self.assertEqual(bl.categorizar("XYZ", 1, ["ZN"])[0], "ion")
        self.assertEqual(bl.categorizar("XYZ", 12, ["C"])[0], "hetero_desconhecido")

    def test_subtipo_de_residuo_modificado(self):
        self.assertEqual(bl.categorizar("ALA")[1], "residuo_padrao")
        self.assertEqual(bl.categorizar("MSE")[1], "residuo_modificado")


# ===========================================================================
# bolso


class TesteBolso(TempDirCase):
    def setUp(self):
        super().setUp()
        self.pdb = escrever(self.dir, "m.pdb", PDB_MINIMO)
        self.atomos = bl.ler_estrutura(self.pdb)

    def contatos_de(self, **kw):
        lig, _, _ = bl.selecionar_ligante(self.atomos, resname="LIG")
        return bl.contatos(self.atomos, lig, **kw)

    def test_corte_exato_em_5A(self):
        """TRP3 a 4,99 Å entra; LEU4 a 5,10 Å fica fora. É o ponto do teste."""
        lista, _ = self.contatos_de(raio=5.0)
        rotulos = {c["residuo"] for c in lista}
        self.assertIn("A/TRP3", rotulos)
        self.assertNotIn("A/LEU4", rotulos)
        self.assertNotIn("A/MET5", rotulos)
        self.assertNotIn("A/NAG104", rotulos)

    def test_lista_tudo_que_nao_e_proteina(self):
        """O pedido explícito: outras cadeias, cofatores, íons, águas."""
        lista, _ = self.contatos_de(raio=5.0)
        por_cat = {}
        for c in lista:
            por_cat.setdefault(c["categoria"], []).append(c["residuo"])
        self.assertIn("A/ZN101", por_cat["ion"])
        self.assertIn("W/HOH1", por_cat["agua"])
        self.assertIn("A/GOL102", por_cat["aditivo"])
        self.assertIn("A/HEM103", por_cat["cofator"])
        self.assertIn("B/CYS1", por_cat["proteina"])      # outra cadeia
        self.assertIn("A/ALA1", por_cat["proteina"])

    def test_distancia_minima_e_atomos(self):
        lista, pares = self.contatos_de(raio=5.0)
        d = {c["residuo"]: c for c in lista}
        self.assertAlmostEqual(d["A/ALA1"]["distancia_min"], 0.0, places=6)
        self.assertAlmostEqual(d["W/HOH1"]["distancia_min"], 1.0, places=6)
        self.assertAlmostEqual(d["A/TRP3"]["distancia_min"], 4.99, places=6)
        self.assertEqual(d["A/ZN101"]["atomo_mais_proximo"], "ZN")
        self.assertEqual(d["W/HOH1"]["atomo_ligante_mais_proximo"], "C2")
        # o CSV por átomo tem de ter uma linha por par dentro do raio
        self.assertTrue(all(p["distancia"] <= 5.0 for p in pares))
        self.assertEqual(pares, sorted(pares, key=lambda p: p["distancia"]))

    def test_proprio_ligante_fora_da_lista(self):
        lista, _ = self.contatos_de(raio=5.0)
        self.assertTrue(all(c["resname"] != "LIG" for c in lista))

    def test_hidrogenio_ignorado_por_padrao(self):
        lista, pares = self.contatos_de(raio=5.0)
        self.assertTrue(all(p["elemento_ligante"] != "H" for p in pares))
        lista_h, pares_h = self.contatos_de(raio=5.0, com_hidrogenios=True)
        self.assertTrue(any(p["elemento_ligante"] == "H" for p in pares_h))

    def test_sem_agua(self):
        lista, _ = self.contatos_de(raio=5.0, incluir_agua=False)
        self.assertTrue(all(c["categoria"] != "agua" for c in lista))

    def test_sem_agua_vale_tambem_no_csv_por_atomo(self):
        """Se a água saiu da lista de resíduos, não pode ficar no CSV de átomos."""
        lista, pares = self.contatos_de(raio=5.0, incluir_agua=False)
        self.assertTrue(all(c["resname"] != "HOH" for c in lista))
        self.assertTrue(all(p["resname"] != "HOH" for p in pares))
        lista2, pares2 = self.contatos_de(raio=5.0, incluir_agua=True)
        self.assertTrue(any(p["resname"] == "HOH" for p in pares2))

    def test_raio_muda_o_conjunto(self):
        p3, _ = self.contatos_de(raio=3.0)
        p5, _ = self.contatos_de(raio=5.0)
        p8, _ = self.contatos_de(raio=8.0)
        self.assertLess(len(p3), len(p5))
        self.assertLess(len(p5), len(p8))
        self.assertIn("A/NAG104", {c["residuo"] for c in self.contatos_de(raio=30.0)[0]})

    def test_raio_invalido(self):
        with self.assertRaisesRegex(ValueError, "raio"):
            self.contatos_de(raio=0.0)

    def test_resumo_marca_interface(self):
        lig, desc, _ = bl.selecionar_ligante(self.atomos, resname="LIG")
        lista, _ = bl.contatos(self.atomos, lig, raio=5.0)
        s = bl.resumir_bolso(lista, desc, 5.0, bl.inventario(self.atomos))
        self.assertEqual(s["cadeias_polimericas_em_contato"], ["A", "B"])
        self.assertTrue(s["bolso_de_interface"])
        self.assertEqual(s["n_residuos_proteicos"], 4)      # ALA1 GLY2 TRP3 + CYS1
        self.assertAlmostEqual(s["distancia_minima"], 0.0, places=6)

    def test_resumo_sem_interface(self):
        atomos = [a for a in self.atomos if a.cadeia != "B"]
        lig, desc, _ = bl.selecionar_ligante(atomos, resname="LIG")
        lista, _ = bl.contatos(atomos, lig, raio=5.0)
        s = bl.resumir_bolso(lista, desc, 5.0, bl.inventario(atomos))
        self.assertFalse(s["bolso_de_interface"])

    def test_pymol(self):
        lig, desc, _ = bl.selecionar_ligante(self.atomos, resname="LIG")
        lista, _ = bl.contatos(self.atomos, lig, raio=5.0)
        pml = bl.selecao_pymol(lista, desc, 5.0)
        self.assertIn("select bolso,", pml)
        self.assertIn("chain A and resi 1", pml)

    def test_grade_bate_com_forca_bruta(self):
        """A grade de células é otimização; o resultado tem de ser o da conta direta."""
        import itertools
        import random as _r
        rng = _r.Random(3)
        linhas = []
        for i in range(120):
            linhas.append(bl_linha(i + 1, "CA", "ALA", "A", i + 1,
                                   rng.uniform(-15, 15), rng.uniform(-15, 15),
                                   rng.uniform(-15, 15)))
        linhas.append(bl_linha(999, "C1", "LIG", "L", 1, 0.0, 0.0, 0.0, hetero=True))
        p = escrever(self.dir, "rand.pdb", "\n".join(linhas) + "\nEND\n")
        atomos = bl.ler_estrutura(p)
        lig = [a for a in atomos if a.resname == "LIG"]
        lista, _ = bl.contatos(atomos, lig, raio=5.0)
        esperado = set()
        for a, l in itertools.product(atomos, lig):
            if a.resname == "LIG":
                continue
            d = math.dist((a.x, a.y, a.z), (l.x, l.y, l.z))
            if d <= 5.0:
                esperado.add(a.rotulo_residuo)
        self.assertEqual({c["residuo"] for c in lista}, esperado)


def bl_linha(serial, nome, resname, cadeia, resnum, x, y, z, hetero=False):
    return cb._linha_pdb(serial, nome, resname, cadeia, resnum, x, y, z,
                         nome[0], hetero=hetero)


class TesteSelecaoDoLigante(TempDirCase):
    def setUp(self):
        super().setUp()
        self.atomos = bl.ler_estrutura(escrever(self.dir, "m.pdb", PDB_MINIMO))

    def test_por_resname(self):
        lig, desc, _ = bl.selecionar_ligante(self.atomos, resname="LIG")
        self.assertEqual(len(lig), 3)                   # 2 C + 1 H
        self.assertEqual(desc["resname"], "LIG")
        self.assertEqual(desc["como"], "seletor")

    def test_resname_minusculo_funciona(self):
        lig, _, _ = bl.selecionar_ligante(self.atomos, resname="lig")
        self.assertEqual(len(lig), 3)

    def test_auto_pega_o_maior_hetero(self):
        lig, desc, cands = bl.selecionar_ligante(self.atomos, auto=True)
        self.assertEqual(desc["resname"], "LIG")
        self.assertEqual(desc["como"], "auto")
        self.assertTrue(cands)

    def test_sem_seletor_e_sem_auto_e_erro(self):
        with self.assertRaisesRegex(ValueError, "diga qual é o ligante"):
            bl.selecionar_ligante(self.atomos)

    def test_resname_inexistente_lista_os_presentes(self):
        with self.assertRaises(ValueError) as ctx:
            bl.selecionar_ligante(self.atomos, resname="XXX")
        msg = str(ctx.exception)
        self.assertIn("nenhum resíduo casa", msg)
        self.assertIn("LIG", msg)                       # mostra o que existe

    def test_varias_copias_exige_escolha(self):
        texto = PDB_MINIMO.replace("END\n", "") + (
            "HETATM   20  C1  LIG M   1      40.000   0.000   0.000  1.00 0.00           C\n"
            "END\n"
        )
        atomos = bl.ler_estrutura(escrever(self.dir, "duas.pdb", texto))
        with self.assertRaises(ValueError) as ctx:
            bl.selecionar_ligante(atomos, resname="LIG")
        self.assertIn("2 cópias", str(ctx.exception))
        lig, desc, _ = bl.selecionar_ligante(atomos, resname="LIG", copia=2)
        self.assertEqual(desc["cadeia"], "M")
        lig2, desc2, _ = bl.selecionar_ligante(atomos, resname="LIG", cadeia="L")
        self.assertEqual(desc2["cadeia"], "L")

    def test_cadeia_de_proteina_como_ligante_e_erro_claro(self):
        linhas = [bl_linha(i + 1, "CA", "ALA", "A", i + 1, i * 4.0, 0, 0)
                  for i in range(40)]
        atomos = bl.ler_estrutura(escrever(self.dir, "prot.pdb",
                                           "\n".join(linhas) + "\nEND\n"))
        with self.assertRaisesRegex(ValueError, "ela é proteína, não ligante"):
            bl.selecionar_ligante(atomos, cadeia="A")

    def test_peptideo_curto_como_ligante(self):
        linhas = [bl_linha(i + 1, "CA", "ALA", "A", i + 1, i * 4.0, 0, 0)
                  for i in range(40)]
        linhas += [bl_linha(100 + i, "CA", "GLY", "P", i + 1, i * 4.0, 3.0, 0)
                   for i in range(5)]
        atomos = bl.ler_estrutura(escrever(self.dir, "pep.pdb",
                                           "\n".join(linhas) + "\nEND\n"))
        lig, desc, _ = bl.selecionar_ligante(atomos, cadeia="P")
        self.assertEqual(len(lig), 5)
        self.assertEqual(desc["subtipo"], "peptideo_curto")

    def test_auto_sem_hetero_nenhum(self):
        linhas = [bl_linha(i + 1, "CA", "ALA", "A", i + 1, i * 4.0, 0, 0)
                  for i in range(40)]
        atomos = bl.ler_estrutura(escrever(self.dir, "so_prot.pdb",
                                           "\n".join(linhas) + "\nEND\n"))
        with self.assertRaisesRegex(ValueError, "não achou ligante"):
            bl.selecionar_ligante(atomos, auto=True)

    def test_auto_ignora_agua_e_ion(self):
        texto = (
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 0.00           C\n"
            "HETATM    2  O   HOH W   1       3.000   0.000   0.000  1.00 0.00           O\n"
            "HETATM    3 ZN    ZN A 101       5.000   0.000   0.000  1.00 0.00          ZN\n"
            "HETATM    4  C1  GOL A 102       7.000   0.000   0.000  1.00 0.00           C\n"
            "HETATM    5  C2  GOL A 102       8.000   0.000   0.000  1.00 0.00           C\n"
        )
        atomos = bl.ler_estrutura(escrever(self.dir, "ai.pdb", texto))
        lig, desc, _ = bl.selecionar_ligante(atomos, auto=True)
        self.assertEqual(desc["resname"], "GOL")       # não água, não íon


class TesteSequenciaDeCadeia(TempDirCase):
    def test_sequencia_e_ordem(self):
        atomos = bl.ler_estrutura(escrever(self.dir, "m.pdb", PDB_MINIMO))
        seq, chaves = bl.sequencia_cadeia(atomos, "A")
        self.assertEqual(seq, "AGWLM")
        self.assertEqual([c[0] for c in chaves], [1, 2, 3, 4, 5])
        self.assertEqual(bl.sequencia_cadeia(atomos, "B")[0], "C")
        self.assertEqual(bl.sequencia_cadeia(atomos, "L")[0], "")

    def test_modificado_entra_como_o_aminoacido(self):
        texto = (
            "ATOM      1  CA  MSE A   1       0.000   0.000   0.000  1.00 0.00           C\n"
            "ATOM      2  CA  ALA A   2       3.800   0.000   0.000  1.00 0.00           C\n"
        )
        atomos = bl.ler_estrutura(escrever(self.dir, "mse.pdb", texto))
        self.assertEqual(bl.sequencia_cadeia(atomos, "A")[0], "MA")

    def test_fora_de_ordem_no_arquivo(self):
        texto = (
            "ATOM      1  CA  GLY A   2       3.800   0.000   0.000  1.00 0.00           C\n"
            "ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00 0.00           C\n"
        )
        atomos = bl.ler_estrutura(escrever(self.dir, "ord.pdb", texto))
        self.assertEqual(bl.sequencia_cadeia(atomos, "A")[0], "AG")

    def test_inventario_e_cadeias_polimericas(self):
        atomos = bl.ler_estrutura(escrever(self.dir, "m.pdb", PDB_MINIMO))
        inv = bl.inventario(atomos)
        self.assertEqual(inv["cadeias"]["A"]["tipo"], "proteina")
        self.assertEqual(sorted(bl.cadeias_polimericas(atomos, minimo_residuos=1)),
                         ["A", "B"])


# ===========================================================================
# alinhamento e mapeamento


class TesteAlinhamento(unittest.TestCase):
    def test_identicas(self):
        pares, ident = cb.alinhar("ACDEFGHIKL", "ACDEFGHIKL")
        self.assertEqual(len(pares), 10)
        self.assertAlmostEqual(ident, 1.0)
        self.assertEqual(pares[0], (0, 0))

    def test_estrutura_e_pedaco_do_meio(self):
        """Modelo truncado nas duas pontas: o caso normal de cofold/cristal."""
        query = "MMMMACDEFGHIKLMNPQWWWW"
        est = "ACDEFGHIKLMNPQ"
        pares, ident = cb.alinhar(query, est)
        self.assertAlmostEqual(ident, 1.0)
        self.assertEqual(len(pares), len(est))
        self.assertEqual(pares[0], (4, 0))       # começa depois dos 4 M
        self.assertEqual(pares[-1], (17, 13))

    def test_loop_faltando_na_estrutura(self):
        """Buraco de 6 resíduos no modelo não pode tirar o resto de registro."""
        query = "ACDEFGHIKLMNPQRSTVWY"
        est = "ACDEFG" + "TVWY"                  # faltam HIKLMNPQRS
        pares, ident = cb.alinhar(query, est)
        self.assertAlmostEqual(ident, 1.0)
        mapa = {ib: ia for ia, ib in pares}
        self.assertEqual(mapa[0], 0)
        self.assertEqual(mapa[6], 16)            # T volta na posição certa
        self.assertEqual(mapa[9], 19)

    def test_mutacoes_pontuais(self):
        query = "ACDEFGHIKLMNPQRSTVWY"
        est = "ACDWFGHIKLMNPQRSTVWY"
        pares, ident = cb.alinhar(query, est)
        self.assertEqual(len(pares), 20)
        self.assertAlmostEqual(ident, 0.95, places=6)

    def test_cadeia_errada_da_identidade_baixa(self):
        pares, ident = cb.alinhar("ACDEFGHIKLMNPQRSTVWY", "WWWWWWWWWWWWWWWWWWWW")
        self.assertLess(ident, 0.2)

    def test_vazias(self):
        self.assertEqual(cb.alinhar("", "ACDE"), ([], 0.0))
        self.assertEqual(cb.alinhar("ACDE", ""), ([], 0.0))

    def test_mapa_query_para_coluna(self):
        self.assertEqual(cb.mapa_query_para_coluna("A-C--D"), [0, 2, 5])


class TesteCruzamento(TempDirCase):
    """residuos_da_cadeia(): o ponto em que o MSA encontra a estrutura."""

    def monta(self, resnum_inicial=1, seq_est=None, seq_query=None, n_homologos=30):
        import random as _r
        rng = _r.Random(11)
        seq_est = seq_est or "ACDEFGHIKLMNPQRSTVWY"
        seq_query = seq_query or seq_est
        tres = {v: k for k, v in list(bl.AA3_PARA_1.items())[:20]}
        linhas = []
        for i, aa in enumerate(seq_est):
            linhas.append(bl_linha(i + 1, "CA", tres[aa], "A", resnum_inicial + i,
                                   i * 6.0, 0.0, 0.0))
        linhas.append(bl_linha(900, "C1", "LIG", "L", 1, 12.0, 4.0, 0.0, hetero=True))
        pdb = escrever(self.dir, f"c{resnum_inicial}.pdb",
                       "\n".join(linhas) + "\nEND\n")
        regs = [("query", seq_query)]
        for k in range(n_homologos):
            regs.append((f"h{k}", "".join(
                rng.choice(mh.ALFABETO) if (i not in (2, 3, 4) and rng.random() < 0.5)
                else aa for i, aa in enumerate(seq_query))))
        msa = escrever(self.dir, f"c{resnum_inicial}.fasta", fasta(regs))
        atomos = bl.ler_estrutura(pdb)
        lig, _, _ = bl.selecionar_ligante(atomos, resname="LIG")
        contatos, _ = bl.contatos(atomos, lig, raio=5.0)
        analise = mh.analisar_msa(msa)
        return analise, atomos, contatos

    def test_mapeamento_direto(self):
        analise, atomos, contatos = self.monta()
        linhas, diag = cb.residuos_da_cadeia(analise, atomos, "A", contatos)
        self.assertEqual(len(linhas), 20)
        self.assertAlmostEqual(diag["identidade_query_estrutura"], 1.0)
        self.assertEqual(diag["residuos_sem_coluna_no_msa"], 0)
        self.assertEqual([l["coluna_msa"] for l in linhas[:3]], [1, 2, 3])
        # o ligante está sobre o resíduo 3 (x=12): tem de estar no bolso
        no_bolso = [l["resnum"] for l in linhas if l["no_bolso"]]
        self.assertEqual(no_bolso, [3])
        self.assertAlmostEqual(
            [l for l in linhas if l["resnum"] == 3][0]["dist_min_ligante"], 4.0,
            places=3)

    def test_numeracao_deslocada(self):
        """Modelo numerado a partir de 101: o mapeamento não pode usar resnum."""
        analise, atomos, contatos = self.monta(resnum_inicial=101)
        linhas, diag = cb.residuos_da_cadeia(analise, atomos, "A", contatos)
        self.assertAlmostEqual(diag["identidade_query_estrutura"], 1.0)
        self.assertEqual(linhas[0]["resnum"], 101)
        self.assertEqual(linhas[0]["coluna_msa"], 1)
        self.assertEqual([l["resnum"] for l in linhas if l["no_bolso"]], [103])

    def test_msa_de_outra_cadeia_avisa(self):
        analise, atomos, contatos = self.monta(
            seq_est="ACDEFGHIKLMNPQRSTVWY", seq_query="WWWWWWWWWWWWWWWWWWWW")
        linhas, diag = cb.residuos_da_cadeia(analise, atomos, "A", contatos)
        self.assertLess(diag["identidade_query_estrutura"], 0.9)
        self.assertTrue(any("não é desta cadeia" in a for a in diag["avisos"]))

    def test_cadeia_sem_proteina(self):
        analise, atomos, contatos = self.monta()
        with self.assertRaisesRegex(ValueError, "não tem resíduo de aminoácido"):
            cb.residuos_da_cadeia(analise, atomos, "L", contatos)

    def test_estatistica_do_bolso(self):
        analise, atomos, contatos = self.monta()
        linhas, _ = cb.residuos_da_cadeia(analise, atomos, "A", contatos)
        est = cb.estatistica_do_bolso(linhas, escore="jsd", n_perm=500)
        self.assertEqual(est["n_bolso"], 1)
        self.assertGreater(est["n_resto"], 10)
        self.assertIn("amostra pequena", est["veredito"])


# ===========================================================================
# estatística


class TesteEstatistica(unittest.TestCase):
    def test_mann_whitney_separado(self):
        r = cb.mann_whitney([5, 6, 7, 8, 9, 10], [0, 1, 2, 3, 4])
        self.assertAlmostEqual(r["cliff_delta"], 1.0, places=9)
        self.assertLess(r["p"], 0.01)

    def test_mann_whitney_identico(self):
        r = cb.mann_whitney([1, 2, 3, 4], [1, 2, 3, 4])
        self.assertAlmostEqual(r["cliff_delta"], 0.0, places=9)
        self.assertGreater(r["p"], 0.5)

    def test_mann_whitney_valores_conhecidos(self):
        """U de um caso calculável à mão: x=[1,2], y=[3,4] -> U=0, δ=-1."""
        r = cb.mann_whitney([1, 2], [3, 4])
        self.assertAlmostEqual(r["U"], 0.0, places=9)
        self.assertAlmostEqual(r["cliff_delta"], -1.0, places=9)

    def test_mann_whitney_amostra_vazia(self):
        r = cb.mann_whitney([], [1, 2])
        self.assertTrue(math.isnan(r["p"]))

    def test_empates_nao_quebram(self):
        r = cb.mann_whitney([1, 1, 1, 1], [1, 1, 1, 1])
        self.assertAlmostEqual(r["cliff_delta"], 0.0, places=9)
        self.assertAlmostEqual(r["p"], 1.0, places=6)

    def test_permutacao_deterministica(self):
        x, y = [5, 6, 7, 8], [1, 2, 3, 4]
        a = cb.teste_permutacao(x, y, n_perm=500, semente=42)
        b = cb.teste_permutacao(x, y, n_perm=500, semente=42)
        self.assertEqual(a["p"], b["p"])
        self.assertLess(a["p"], 0.1)
        c = cb.teste_permutacao([1, 2, 3, 4], [1, 2, 3, 4], n_perm=500, semente=42)
        self.assertGreater(c["p"], 0.5)

    def test_permutacao_limites(self):
        """p nunca é 0 (correção +1) nem passa de 1."""
        r = cb.teste_permutacao(list(range(50, 70)), list(range(20)),
                                n_perm=200, semente=1)
        self.assertGreater(r["p"], 0.0)
        self.assertLessEqual(r["p"], 1.0)

    def test_spearman(self):
        x = [1, 2, 3, 4, 5]
        self.assertAlmostEqual(cb.spearman(x, x), 1.0, places=9)
        self.assertAlmostEqual(cb.spearman(x, x[::-1]), -1.0, places=9)
        self.assertTrue(math.isnan(cb.spearman([1, 2], [1, 2])))

    def test_mediana(self):
        self.assertEqual(cb._mediana([3, 1, 2]), 2)
        self.assertEqual(cb._mediana([4, 1, 2, 3]), 2.5)
        self.assertTrue(math.isnan(cb._mediana([])))


# ===========================================================================
# configuração


class TesteConfig(TempDirCase):
    def test_modelo_de_config_e_json_valido(self):
        p = self.dir / "cfg.json"
        p.write_text(json.dumps(cb.MODELO_CONFIG, indent=2, ensure_ascii=False))
        cfg = cb.ler_config(p)
        self.assertIn("msas", cfg)
        self.assertEqual(cfg["raio_angstrom"], 5.0)

    def test_placeholders_sao_detectados(self):
        """O modelo em branco NÃO pode passar pela validação."""
        erros = cb.validar_config(dict(cb.MODELO_CONFIG))
        self.assertTrue(any("estrutura" in e for e in erros))
        self.assertTrue(any("MSA não encontrado" in e or "não preenchido" in e
                            for e in erros))

    def test_config_parcial_herda_o_resto(self):
        p = self.dir / "parcial.json"
        p.write_text(json.dumps({"raio_angstrom": 4.5,
                                 "escores": {"janela": 1}}))
        cfg = cb.ler_config(p)
        self.assertEqual(cfg["raio_angstrom"], 4.5)
        self.assertEqual(cfg["escores"]["janela"], 1)
        self.assertEqual(cfg["escores"]["peso_janela"], 0.5)     # veio do modelo
        self.assertEqual(cfg["pesos"]["normalizacao"], "media")

    def test_json_invalido(self):
        p = self.dir / "ruim.json"
        p.write_text("{isto não é json}")
        with self.assertRaisesRegex(ValueError, "JSON inválido"):
            cb.ler_config(p)

    def test_config_inexistente(self):
        with self.assertRaisesRegex(ValueError, "não encontrada"):
            cb.ler_config(self.dir / "nada.json")

    def test_validacao_lista_todos_os_erros(self):
        cfg = json.loads(json.dumps(cb.MODELO_CONFIG))
        cfg["estrutura"] = str(escrever(self.dir, "e.pdb", PDB_MINIMO))
        cfg["raio_angstrom"] = -1
        cfg["escores"]["escore_principal"] = "inventado"
        cfg["pesos"]["normalizacao"] = "nada"
        cfg["msas"] = [{"ferramenta": "x", "cadeia": "", "arquivo": ""}]
        erros = cb.validar_config(cfg)
        self.assertTrue(any("raio_angstrom" in e for e in erros))
        self.assertTrue(any("escore_principal" in e for e in erros))
        self.assertTrue(any("normalizacao" in e for e in erros))
        self.assertTrue(any("'arquivo' não preenchido" in e for e in erros))
        self.assertTrue(any("'cadeia' não preenchida" in e for e in erros))

    def test_par_ferramenta_cadeia_repetido(self):
        f = escrever(self.dir, "m.fasta", fasta(EXEMPLO_HH))
        cfg = json.loads(json.dumps(cb.MODELO_CONFIG))
        cfg["estrutura"] = str(escrever(self.dir, "e.pdb", PDB_MINIMO))
        cfg["ligante"] = {"resname": "LIG"}
        cfg["msas"] = [{"ferramenta": "b", "cadeia": "A", "arquivo": str(f)},
                       {"ferramenta": "b", "cadeia": "A", "arquivo": str(f)}]
        self.assertTrue(any("repetido" in e for e in cb.validar_config(cfg)))

    def test_ligante_sem_seletor(self):
        f = escrever(self.dir, "m.fasta", fasta(EXEMPLO_HH))
        cfg = json.loads(json.dumps(cb.MODELO_CONFIG))
        cfg["estrutura"] = str(escrever(self.dir, "e.pdb", PDB_MINIMO))
        cfg["msas"] = [{"ferramenta": "b", "cadeia": "A", "arquivo": str(f)}]
        cfg["ligante"] = {"resname": None, "cadeia": None, "resnum": None,
                          "auto": False}
        self.assertTrue(any("ligante" in e for e in cb.validar_config(cfg)))

    def test_escore_sem_peso_exige_a_passagem_sem_peso(self):
        """Pedir 'jsd_sem_peso' com comparar_sem_peso desligado é incoerente."""
        f = escrever(self.dir, "m.fasta", fasta(EXEMPLO_HH))
        cfg = json.loads(json.dumps(cb.MODELO_CONFIG))
        cfg["estrutura"] = str(escrever(self.dir, "e.pdb", PDB_MINIMO))
        cfg["ligante"] = {"resname": "LIG"}
        cfg["msas"] = [{"ferramenta": "b", "cadeia": "A", "arquivo": str(f)}]
        cfg["escores"]["escore_principal"] = "jsd_sem_peso"
        cfg["escores"]["comparar_sem_peso"] = False
        self.assertTrue(any("comparar_sem_peso" in e for e in cb.validar_config(cfg)))
        cfg["escores"]["comparar_sem_peso"] = True
        self.assertFalse(cb.validar_config(cfg))

    def test_parse_do_argumento_msa(self):
        d = cb._parse_msa_arg("boltz2:A:/dados/msa.fasta")
        self.assertEqual((d["ferramenta"], d["cadeia"], d["arquivo"]),
                         ("boltz2", "A", "/dados/msa.fasta"))
        d2 = cb._parse_msa_arg("B:/dados/outro.a3m")
        self.assertEqual((d2["cadeia"], d2["arquivo"]), ("B", "/dados/outro.a3m"))
        self.assertEqual(d2["ferramenta"], "outro")
        with self.assertRaisesRegex(ValueError, "malformado"):
            cb._parse_msa_arg("so_um_campo")


# ===========================================================================
# pipeline completo


class TesteDemo(TempDirCase):
    """A demo é o autoteste do script: entrada conhecida, resposta conhecida."""

    def test_gerar_demo_escreve_tudo(self):
        d = cb.gerar_demo(self.dir / "entrada")
        for k in ("pdb", "msa_a", "msa_b", "config"):
            self.assertTrue(Path(d[k]).exists(), msg=k)
        self.assertEqual(len(d["seq_a"]), 30)
        self.assertEqual(d["bolso_esperado_a"], [11, 12, 13])

    def test_bolso_da_demo_e_o_esperado(self):
        d = cb.gerar_demo(self.dir / "entrada")
        r = bl.analisar_bolso(d["pdb"], resname="LIG", raio=5.0)
        # esperado: resíduos 11-13 da cadeia A, 1-3 da B, Zn, água e glicerol
        self.assertEqual(len(r["contatos"]), 9)
        self.assertEqual(
            sorted(c["resnum"] for c in r["contatos"] if c["cadeia"] == "A"
                   and c["categoria"] == "proteina"),
            [11, 12, 13])
        self.assertEqual(
            sorted(c["resnum"] for c in r["contatos"] if c["cadeia"] == "B"),
            [1, 2, 3])
        cats = r["resumo"]["por_categoria"]
        self.assertEqual(cats["ion"], 1)
        self.assertEqual(cats["agua"], 1)
        self.assertEqual(cats["aditivo"], 1)
        self.assertEqual(cats["proteina"], 6)
        self.assertTrue(r["resumo"]["bolso_de_interface"])
        self.assertEqual(r["resumo"]["cadeias_polimericas_em_contato"], ["A", "B"])

    def test_pipeline_da_demo_detecta_bolso_conservado(self):
        d = cb.gerar_demo(self.dir / "entrada")
        r = cb.executar(d["cfg"], self.dir / "saida")
        self.assertTrue((self.dir / "saida" / "relatorio.md").exists())
        self.assertTrue((self.dir / "saida" / "resumo.json").exists())
        cadeia_a = next(c for c in r["cadeias"] if c["cadeia"] == "A")
        est = cadeia_a["estatistica"]
        self.assertEqual(est["n_bolso"], 3)
        self.assertGreater(est["mediana_bolso"], est["mediana_resto"])
        self.assertLess(est["permutacao"]["p"], 0.05)
        self.assertIn("MAIS conservado", est["veredito"])
        # a cadeia B tem MSA aleatório: não pode dar sinal
        cadeia_b = next(c for c in r["cadeias"] if c["cadeia"] == "B")
        self.assertGreater(cadeia_b["estatistica"]["permutacao"]["p"], 0.05)

    def test_arquivos_de_saida_e_cabecalhos(self):
        d = cb.gerar_demo(self.dir / "entrada")
        saida = self.dir / "saida"
        cb.executar(d["cfg"], saida)
        esperados = [
            "relatorio.md", "resumo.json", "bolso.pml",
            "contatos_5A_residuos.csv", "contatos_5A_atomos.csv",
            "conservacao_colunas__demo__A.csv", "pesos_henikoff__demo__A.csv",
            "residuos__demo__A.csv", "residuos__demo__B.csv",
        ]
        for nome in esperados:
            self.assertTrue((saida / nome).exists(), msg=nome)
        import csv as _csv
        with (saida / "residuos__demo__A.csv").open() as fh:
            linhas = list(_csv.DictReader(fh))
        self.assertEqual(len(linhas), 30)
        self.assertEqual(sum(int(l["no_bolso"]) for l in linhas), 3)
        with (saida / "contatos_5A_residuos.csv").open() as fh:
            contatos = list(_csv.DictReader(fh))
        self.assertEqual(len(contatos), 9)
        self.assertIn("categoria", contatos[0])
        self.assertIn("mesma_cadeia_que_o_ligante", contatos[0])
        # o relatório tem de citar o íon, a água e o aditivo por nome
        rel = (saida / "relatorio.md").read_text()
        for trecho in ("ZN201", "HOH1", "GOL2", "bolso de interface",
                       "Henikoff"):
            self.assertIn(trecho, rel)

    def test_comparacao_entre_ferramentas(self):
        """Duas 'ferramentas' para a mesma cadeia -> tabela de concordância."""
        d = cb.gerar_demo(self.dir / "entrada")
        cfg = json.loads(json.dumps(d["cfg"]))
        cfg["msas"].append({"ferramenta": "outra", "cadeia": "A",
                            "arquivo": str(d["msa_a"]), "formato": "auto"})
        r = cb.executar(cfg, self.dir / "saida2")
        comp = r["comparacao_ferramentas"]
        self.assertTrue(comp)
        linha = next(c for c in comp if {c["ferramenta_a"], c["ferramenta_b"]}
                     == {"demo", "outra"})
        self.assertAlmostEqual(linha["spearman"], 1.0, places=6)  # mesmo arquivo
        self.assertTrue((self.dir / "saida2" / "comparacao_ferramentas.csv").exists())

    def test_sem_estrutura_so_msa(self):
        d = cb.gerar_demo(self.dir / "entrada")
        cfg = json.loads(json.dumps(d["cfg"]))
        cfg["estrutura"] = None
        r = cb.executar(cfg, self.dir / "saida3", sem_estrutura=True)
        self.assertNotIn("bolso", r)
        self.assertEqual(len(r["msas"]), 2)
        self.assertFalse(r["cadeias"])

    def test_msa_quebrado_nao_derruba_o_resto(self):
        """Um MSA ruim vira erro NAQUELA linha, e o resto do relatório sai."""
        d = cb.gerar_demo(self.dir / "entrada")
        ruim = escrever(self.dir, "ruim.fasta", ">q\nACDE\n>h\nAC\n")
        cfg = json.loads(json.dumps(d["cfg"]))
        cfg["msas"].append({"ferramenta": "quebrado", "cadeia": "A",
                            "arquivo": str(ruim), "formato": "auto"})
        r = cb.executar(cfg, self.dir / "saida4")
        com_erro = [m for m in r["msas"] if m.get("erro")]
        self.assertEqual(len(com_erro), 1)
        self.assertIn("retangular", com_erro[0]["erro"])
        self.assertEqual(len([m for m in r["msas"] if not m.get("erro")]), 2)

    def test_todos_os_msas_quebrados_e_erro(self):
        d = cb.gerar_demo(self.dir / "entrada")
        ruim = escrever(self.dir, "ruim.fasta", ">q\nACDE\n>h\nAC\n")
        cfg = json.loads(json.dumps(d["cfg"]))
        cfg["msas"] = [{"ferramenta": "x", "cadeia": "A", "arquivo": str(ruim)}]
        with self.assertRaisesRegex(ValueError, "nenhum MSA"):
            cb.executar(cfg, self.dir / "saida5")

    def test_config_incompleta_falha_antes_de_calcular(self):
        with self.assertRaisesRegex(ValueError, "configuração incompleta"):
            cb.executar(json.loads(json.dumps(cb.MODELO_CONFIG)), self.dir / "s6")

    def test_escores_alternativos_rodam(self):
        d = cb.gerar_demo(self.dir / "entrada")
        for escore in cb.ESCORES_VALIDOS:
            cfg = json.loads(json.dumps(d["cfg"]))
            cfg["escores"]["escore_principal"] = escore
            cfg["escores"]["n_permutacoes"] = 200
            r = cb.executar(cfg, self.dir / f"s_{escore}")
            self.assertTrue(r["cadeias"], msg=escore)

    def test_relatorio_sem_nada_no_raio(self):
        """Raio pequeno: o relatório tem de sair mesmo com bolso vazio."""
        d = cb.gerar_demo(self.dir / "entrada")
        cfg = json.loads(json.dumps(d["cfg"]))
        cfg["raio_angstrom"] = 0.5
        r = cb.executar(cfg, self.dir / "s7")
        self.assertEqual(r["bolso"]["n_itens_no_raio"], 0)
        self.assertIn("Nenhum cofator", (self.dir / "s7" / "relatorio.md").read_text())


class TesteCLI(TempDirCase):
    """As três CLIs, rodadas como o usuário roda: subprocess, código de saída."""

    def roda(self, script, *args):
        return subprocess.run(
            [PY, str(RAIZ / "scripts" / script), *map(str, args)],
            capture_output=True, text=True,
        )

    def test_cli_demo(self):
        r = self.roda("conservacao_bolso.py", "--demo", "--saida", self.dir / "d")
        self.assertEqual(r.returncode, 0, msg=r.stderr)
        self.assertIn("INTERFACE", r.stdout)
        self.assertIn("MAIS conservado", r.stdout)
        self.assertTrue((self.dir / "d" / "relatorio.md").exists())

    def test_cli_escrever_config(self):
        p = self.dir / "cfg.json"
        r = self.roda("conservacao_bolso.py", "--escrever-config", p)
        self.assertEqual(r.returncode, 0, msg=r.stderr)
        self.assertTrue(p.exists())
        cfg = json.loads(p.read_text())
        self.assertIn("msas", cfg)
        self.assertTrue(cb.validar_config(cfg))        # em branco: inválido

    def test_cli_com_config(self):
        d = cb.gerar_demo(self.dir / "entrada")
        r = self.roda("conservacao_bolso.py", "--config", d["config"],
                      "--saida", self.dir / "s")
        self.assertEqual(r.returncode, 0, msg=r.stderr)
        self.assertTrue((self.dir / "s" / "resumo.json").exists())

    def test_cli_argumentos_diretos(self):
        d = cb.gerar_demo(self.dir / "entrada")
        r = self.roda("conservacao_bolso.py",
                      "--estrutura", d["pdb"], "--ligante-resname", "LIG",
                      "--msa", f"boltz2:A:{d['msa_a']}",
                      "--msa", f"boltz2:B:{d['msa_b']}",
                      "--raio", "5.0", "--n-permutacoes", "500",
                      "--saida", self.dir / "s2")
        self.assertEqual(r.returncode, 0, msg=r.stderr)
        self.assertIn("boltz2__A", (self.dir / "s2" / "resumo.json").read_text())

    def test_cli_erro_sai_com_codigo_2(self):
        r = self.roda("conservacao_bolso.py", "--estrutura", "/nao/existe.pdb",
                      "--ligante-resname", "LIG", "--msa", "x:A:/nao/existe.fasta",
                      "--saida", self.dir / "s3")
        self.assertEqual(r.returncode, 2)
        self.assertIn("ERRO", r.stderr)

    def test_cli_ajuda(self):
        for script in ("msa_henikoff.py", "bolso_ligante.py", "conservacao_bolso.py"):
            r = self.roda(script, "--help")
            self.assertEqual(r.returncode, 0, msg=f"{script}: {r.stderr}")

    def test_cli_msa_henikoff(self):
        d = cb.gerar_demo(self.dir / "entrada")
        r = self.roda("msa_henikoff.py", "--msa", d["msa_a"],
                      "--rotulo", "boltz2", "--saida", self.dir / "m")
        self.assertEqual(r.returncode, 0, msg=r.stderr)
        self.assertIn("Neff", r.stdout)
        self.assertTrue((self.dir / "m" / "conservacao_colunas__boltz2.csv").exists())

    def test_cli_bolso_ligante(self):
        d = cb.gerar_demo(self.dir / "entrada")
        r = self.roda("bolso_ligante.py", "--estrutura", d["pdb"],
                      "--ligante-auto", "--raio", "5.0", "--saida", self.dir / "b")
        self.assertEqual(r.returncode, 0, msg=r.stderr)
        self.assertIn("ion", r.stdout)
        self.assertIn("INTERFACE", r.stdout)

    def test_cli_bolso_sem_ligante_definido(self):
        d = cb.gerar_demo(self.dir / "entrada")
        r = self.roda("bolso_ligante.py", "--estrutura", d["pdb"],
                      "--saida", self.dir / "b2")
        self.assertEqual(r.returncode, 2)
        self.assertIn("diga qual é o ligante", r.stderr)


# mmCIF no estilo do que o PDB (e o Boltz-2/Chai) realmente escreve, com as
# manhas que derrubam parser ingênuo:
#   - label_seq_id '.' nos não-polímeros (o número está em auth_seq_id)
#   - código de inserção (52 e 52A são DOIS resíduos)
#   - altLoc A/B na mesma side chain (fica a A)
#   - nome de átomo entre aspas com apóstrofo ("C1'") e de 4 caracteres
#   - auth_asym_id de 2 caracteres
#   - dois modelos (pdbx_PDB_model_num)
# Geometria com resposta conhecida, medida ao mais próximo entre C1 (0,0,0) e
# C2 (3,0,0) do ligante: TRP54 a EXATAMENTE 5,000 Å entra; LEU55 a 5,001 fica
# fora.
CIF_REALISTA = """\
data_TESTE
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM   1  N  N     . ALA AA 1 1  ? 0.000 1.000  4.000 1.00 20.0 52  ALA AA N     1
ATOM   2  C  CA    . ALA AA 1 1  ? 0.000 0.000  4.000 1.00 20.0 52  ALA AA CA    1
ATOM   3  C  CA    . ALA AA 1 2  A 3.000 0.000  4.900 1.00 20.0 52  ALA AA CA    1
ATOM   4  O  OG    A SER AA 1 3  ? 0.000 0.000  4.500 0.50 20.0 53  SER AA OG    1
ATOM   5  O  OG    B SER AA 1 3  ? 0.000 0.000 20.000 0.50 20.0 53  SER AA OG    1
ATOM   6  C  CA    . TRP AA 1 4  ? 0.000 0.000  5.000 1.00 20.0 54  TRP AA CA    1
ATOM   7  N  HE21  . TRP AA 1 4  ? 0.000 1.000  5.000 1.00 20.0 54  TRP AA HE21  1
ATOM   8  C  CA    . LEU AA 1 5  ? 0.000 0.000  5.001 1.00 20.0 55  LEU AA CA    1
ATOM   9  C  CA    . GLY B  2 1  ? 3.000 3.900  0.000 1.00 25.0 10  GLY B  CA    1
ATOM   10 C  "C1'" . DA  D  3 1  ? 9.000 0.000  0.000 1.00 30.0 1   DA  D  "C1'" 1
HETATM 11 ZN ZN    . ZN  C  4 .  ? 3.000 0.000 -2.500 1.00 15.0 501 ZN  AA ZN    1
HETATM 12 P  PA    . NAD E  5 .  ? -4.50 0.000  0.000 1.00 25.0 601 NAD AA PA    1
HETATM 13 C  C1    . NAG F  6 .  ? 6.000 0.000  3.000 1.00 25.0 1   NAG G  C1    1
HETATM 14 O  O     . HOH W  7 .  ? 1.500 1.500  0.000 1.00 30.0 701 HOH AA O     1
HETATM 15 C  C1    . LIG L  8 .  ? 0.000 0.000  0.000 1.00 18.0 401 LIG L  C1    1
HETATM 16 C  C2    . LIG L  8 .  ? 3.000 0.000  0.000 1.00 18.0 401 LIG L  C2    1
ATOM   17 C  CA    . ALA AA 1 1  ? 50.00 0.000  0.000 1.00 20.0 52  ALA AA CA    2
HETATM 18 C  C1    . LIG L  8 .  ? 50.00 0.000  0.000 1.00 18.0 401 LIG L  C1    2
#
"""


class TesteCifRealista(TempDirCase):
    """O parser de mmCIF contra um arquivo com as manhas dos arquivos de verdade."""

    def setUp(self):
        super().setUp()
        self.cif = escrever(self.dir, "realista.cif", CIF_REALISTA)
        self.atomos = bl.ler_estrutura(self.cif)

    def test_altloc_e_modelo(self):
        # 16 átomos no modelo 1, menos o altLoc B do SER53 = 15
        self.assertEqual(len(self.atomos), 15)
        ogs = [a for a in self.atomos if a.nome == "OG"]
        self.assertEqual(len(ogs), 1)
        self.assertAlmostEqual(ogs[0].z, 4.5, places=3)
        m2 = bl.ler_estrutura(self.cif, modelo=2)
        self.assertEqual(len(m2), 2)
        self.assertAlmostEqual(m2[0].x, 50.0, places=3)

    def test_auth_tem_prioridade(self):
        """auth_seq_id/auth_asym_id são a numeração que o usuário vê no PyMOL."""
        ala = next(a for a in self.atomos if a.resname == "ALA" and a.nome == "N")
        self.assertEqual((ala.cadeia, ala.resnum), ("AA", 52))
        zn = next(a for a in self.atomos if a.resname == "ZN")
        self.assertEqual((zn.cadeia, zn.resnum), ("AA", 501))   # não a label_asym 'C'
        nag = next(a for a in self.atomos if a.resname == "NAG")
        self.assertEqual(nag.cadeia, "G")

    def test_label_seq_id_ponto_nao_zera_o_hetero(self):
        lig = [a for a in self.atomos if a.resname == "LIG"]
        self.assertEqual({a.resnum for a in lig}, {401})

    def test_codigo_de_insercao_separa_residuos(self):
        """52 e 52A são dois resíduos: tratá-los como um some com um contato."""
        chaves = {a.chave_residuo for a in self.atomos if a.resname == "ALA"}
        self.assertEqual(len(chaves), 2)
        icodes = sorted((a.icode.strip() for a in self.atomos
                         if a.resname == "ALA"))
        self.assertIn("A", icodes)

    def test_nomes_de_atomo_com_aspas_e_quatro_caracteres(self):
        nomes = {a.nome for a in self.atomos}
        self.assertIn("C1'", nomes)
        self.assertIn("HE21", nomes)

    def test_bolso_do_cif_realista(self):
        lig, desc, _ = bl.selecionar_ligante(self.atomos, resname="LIG")
        lista, pares = bl.contatos(self.atomos, lig, raio=5.0)
        rotulos = {c["residuo"] for c in lista}
        esperado = {
            "AA/ALA52",      # 4,00 Å
            "AA/ALA52A",     # 4,90 Å — só existe se o icode for respeitado
            "AA/SER53",      # 4,50 Å (altLoc A)
            "AA/TRP54",      # 5,000 Å exatos: entra
            "B/GLY10",       # outra cadeia -> interface
            "AA/ZN501",      # íon
            "AA/NAD601",     # cofator
            "G/NAG1",        # glicano
            "AA/HOH701",     # água
        }
        self.assertEqual(rotulos, esperado)
        self.assertNotIn("AA/LEU55", rotulos)         # 5,001 Å: fora
        self.assertNotIn("D/DA1", rotulos)            # 9 Å: fora
        cat = {c["residuo"]: c["categoria"] for c in lista}
        self.assertEqual(cat["AA/ZN501"], "ion")
        self.assertEqual(cat["AA/NAD601"], "cofator")
        self.assertEqual(cat["G/NAG1"], "glicano")
        self.assertEqual(cat["AA/HOH701"], "agua")
        self.assertEqual(cat["B/GLY10"], "proteina")
        # o H do TRP não pode contar como átomo em contato
        self.assertTrue(all(p["elemento"] != "H" for p in pares))

    def test_resumo_do_cif_realista(self):
        lig, desc, _ = bl.selecionar_ligante(self.atomos, resname="LIG")
        lista, _ = bl.contatos(self.atomos, lig, raio=5.0)
        s = bl.resumir_bolso(lista, desc, 5.0, bl.inventario(self.atomos))
        self.assertEqual(s["cadeias_polimericas_em_contato"], ["AA", "B"])
        self.assertTrue(s["bolso_de_interface"])
        self.assertEqual(s["por_categoria"],
                         {"agua": 1, "cofator": 1, "glicano": 1, "ion": 1,
                          "proteina": 5})
        self.assertEqual(len(s["nao_proteicos_no_raio"]), 4)

    def test_sequencia_da_cadeia_com_insercao(self):
        seq, chaves = bl.sequencia_cadeia(self.atomos, "AA")
        self.assertEqual(seq, "AASWL")               # ALA52, ALA52A, SER53, TRP54, LEU55
        self.assertEqual(chaves[1][1].strip(), "A")  # o resíduo inserido

    def test_acido_nucleico_detectado(self):
        inv = bl.inventario(self.atomos)
        self.assertEqual(inv["cadeias"]["D"]["tipo"], "acido_nucleico")

    def test_pipeline_completo_no_cif_realista(self):
        """Fim a fim num arquivo realista, com MSA e código de inserção."""
        seq, _ = bl.sequencia_cadeia(self.atomos, "AA")
        import random as _r
        rng = _r.Random(5)
        regs = [("query", seq)]
        for k in range(25):
            regs.append((f"h{k}", "".join(
                rng.choice(mh.ALFABETO) if (i > 2 and rng.random() < 0.6) else aa
                for i, aa in enumerate(seq))))
        msa = escrever(self.dir, "msa_AA.fasta", fasta(regs))
        cfg = json.loads(json.dumps(cb.MODELO_CONFIG))
        cfg["estrutura"] = str(self.cif)
        cfg["ligante"] = {"resname": "LIG", "cadeia": None, "resnum": None,
                          "copia": None, "auto": False}
        cfg["msas"] = [{"ferramenta": "boltz2", "cadeia": "AA",
                        "arquivo": str(msa), "formato": "auto"}]
        cfg["escores"]["n_permutacoes"] = 300
        cfg["escores"]["so_colunas_confiaveis"] = False
        r = cb.executar(cfg, self.dir / "saida")
        self.assertTrue(r["bolso"]["bolso_de_interface"])
        c = r["cadeias"][0]
        self.assertAlmostEqual(c["diagnostico"]["identidade_query_estrutura"], 1.0)
        self.assertEqual(c["diagnostico"]["residuos_no_bolso"], 4)
        rel = (self.dir / "saida" / "relatorio.md").read_text()
        self.assertIn("NAD601", rel)
        self.assertIn("ZN501", rel)


class TesteNotebook(TempDirCase):
    """O notebook tem de RODAR, não só ter sintaxe válida.

    Executa as células de código em ordem, num só namespace, com a saída
    desviada para um diretório temporário e no modo de exemplo sintético — que
    é exatamente como o usuário vai rodá-lo na primeira vez.
    """

    NOTEBOOK = RAIZ / "notebooks" / "conservacao_msa_cofold.ipynb"

    def test_notebook_existe_e_e_json_valido(self):
        self.assertTrue(self.NOTEBOOK.exists(),
                        "rode: python scripts/build_notebook_conservacao.py")
        nb = json.loads(self.NOTEBOOK.read_text())
        self.assertEqual(nb["nbformat"], 4)
        self.assertTrue(any(c["cell_type"] == "code" for c in nb["cells"]))

    def test_celulas_de_codigo_compilam(self):
        nb = json.loads(self.NOTEBOOK.read_text())
        for i, c in enumerate(nb["cells"]):
            if c["cell_type"] == "code":
                compile("".join(c["source"]), f"<celula {i}>", "exec")

    def test_notebook_roda_inteiro_no_modo_exemplo(self):
        nb = json.loads(self.NOTEBOOK.read_text())
        destino = self.dir / "saida_nb"
        ns = {"__name__": "__main__", "display": lambda *a, **k: None}
        n = 0
        for i, c in enumerate(nb["cells"]):
            if c["cell_type"] != "code":
                continue
            src = "".join(c["source"])
            src = src.replace('SAIDA = Path.home() / "conservacao_cofold"',
                              f'SAIDA = Path(r"{destino}")')
            # a demo já é o padrão do notebook; menos permutações para ir rápido
            src = src.replace('"n_permutacoes": 10000,', '"n_permutacoes": 500,')
            try:
                exec(compile(src, f"<celula {i}>", "exec"), ns)
            except Exception as erro:            # pragma: no cover
                self.fail(f"célula {i} do notebook falhou: {type(erro).__name__}: {erro}")
            n += 1
        self.assertGreaterEqual(n, 10)
        self.assertTrue((destino / "relatorio.md").exists())
        self.assertTrue((destino / "contatos_5A_residuos.csv").exists())
        self.assertTrue(ns["USAR_EXEMPLO_SINTETICO"])
        # o notebook tem de ter chegado ao veredito certo da demo
        cadeia_a = next(r for r in ns["resultados_cadeia"] if r["cadeia"] == "A")
        self.assertIn("MAIS conservado", cadeia_a["estatistica"]["veredito"])

    def test_notebook_em_branco_nao_passa_pela_validacao(self):
        """Com USAR_EXEMPLO_SINTETICO=False e nada preenchido, tem de reclamar."""
        nb = json.loads(self.NOTEBOOK.read_text())
        ns = {"__name__": "__main__", "display": lambda *a, **k: None}
        saida = []
        for i, c in enumerate(nb["cells"][:8]):
            if c["cell_type"] != "code":
                continue
            src = "".join(c["source"])
            src = src.replace("USAR_EXEMPLO_SINTETICO = True",
                              "USAR_EXEMPLO_SINTETICO = False")
            src = src.replace('SAIDA = Path.home() / "conservacao_cofold"',
                              f'SAIDA = Path(r"{self.dir / "s"}")')
            import contextlib
            import io
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                exec(compile(src, f"<celula {i}>", "exec"), ns)
            saida.append(buf.getvalue())
        texto = "\n".join(saida)
        self.assertIn("CONFIGURAÇÃO INCOMPLETA", texto)
        self.assertIn("estrutura", texto)


if __name__ == "__main__":
    unittest.main(verbosity=2)
