# Conservação dos MSAs de cofold (pesos de Henikoff) + bolso a 5 Å

Avalia **quanto se pode confiar** no MSA que alimentou um modelo de cofold
(Boltz-2, Chai-1/2, ESMFold, OpenFold) e **se o bolso do ligante cai em região
conservada** — listando, junto, tudo que está a menos de 5 Å do ligante:
outras cadeias, cofatores, íons, glicanos, aditivos de cristalização e águas.

Só biblioteca padrão. `numpy` acelera; `pandas`/`matplotlib` só enfeitam o
notebook. Nada disso é exigido.

## Arquivos

| arquivo | o que faz |
|---|---|
| `notebooks/conservacao_msa_cofold.ipynb` | o notebook: campos de entrada genéricos, tabelas, gráficos, exportação |
| `scripts/msa_henikoff.py` | leitura do MSA, pesos de Henikoff, escores por coluna |
| `scripts/bolso_ligante.py` | PDB/mmCIF, seleção do ligante, tudo a ≤ N Å com categoria |
| `scripts/conservacao_bolso.py` | junta os dois: mapeia MSA↔estrutura, testa bolso × resto, escreve o relatório |
| `scripts/build_notebook_conservacao.py` | gera o notebook (edite aqui, não o `.ipynb`) |
| `tests/test_conservacao_msa.py` | 147 testes, incluindo rodar o notebook inteiro |

## Uso

```bash
# 0) conferir que tudo funciona, sem ter nenhum dado ainda
python scripts/conservacao_bolso.py --demo --saida /tmp/demo_conservacao

# 1) escrever a configuração em branco e preencher
python scripts/conservacao_bolso.py --escrever-config config_conservacao.json
python scripts/conservacao_bolso.py --config config_conservacao.json --saida resultados/

# ou direto na linha de comando
python scripts/conservacao_bolso.py \
    --estrutura modelo_boltz2.cif --ligante-resname LIG --raio 5.0 \
    --msa boltz2:A:msa_alvo.fasta --msa boltz2:B:msa_e3.fasta \
    --saida resultados/

# as duas partes também rodam sozinhas
python scripts/msa_henikoff.py --msa msa.a3m --rotulo boltz2 --saida conservacao/
python scripts/bolso_ligante.py --estrutura modelo.cif --ligante-auto --raio 5.0

# testes
python tests/test_conservacao_msa.py
```

No notebook, `USAR_EXEMPLO_SINTETICO = True` (padrão) gera um caso com resposta
conhecida e roda o notebook inteiro — serve para validar a instalação antes de
os resultados chegarem.

## De onde vem o MSA de cada ferramenta

| Ferramenta | O que ela usa | Formato |
|---|---|---|
| Boltz-2 | o `.a3m`/`.csv` do YAML, ou o que o `--use_msa_server` baixou | a3m, csv |
| Chai-1 / Chai-2 | `*.aligned.pqt`, ou o a3m do ColabFold | parquet, a3m |
| OpenFold / AlphaFold | o diretório `msas/`: `*.sto` (jackhmmer) e `*.a3m` (hhblits) | sto, a3m |
| ESMFold | **nenhum** — é single-sequence por construção | — |

Os quatro formatos são lidos direto, além de **FASTA alinhado**, que é o
formato recomendado para exportar. O ESMFold não tem conservação a medir: o
script **avisa** em vez de produzir número (`--permitir-sequencia-unica` força,
mas o veredito diz que não se aplica).

Parquet exige `pyarrow` no env; sem ele, o erro diz para exportar FASTA.

## As decisões que não são óbvias

1. **Por que pesar.** Esses MSAs vêm de bancos redundantes. Dez mil ortólogos
   de primata não são dez mil observações independentes. Henikoff & Henikoff
   (1994): numa coluna com `r` tipos, um resíduo do tipo `a` que aparece `s`
   vezes contribui `1/(r·s)` para o peso da sequência. `N/Neff` é o fator de
   redundância, e a tabela "efeito dos pesos" mostra a conservação com e sem
   peso — a diferença é a inflação que a redundância causava.

2. **Gap não é o 21º tipo** (padrão). Senão um bloco de sequências truncadas
   vira "tipo raro" e ganha peso alto por não ter sido sequenciado.
   (`--gap-como-tipo` inverte.)

3. **Normalização `media`** (padrão): a soma das contribuições é dividida pelo
   número de colunas que a sequência cobre. MSA de servidor é cheio de hit
   parcial, e na soma crua um fragmento de 40 resíduos em 400 colunas fica com
   peso 10x menor só por ser curto. (`--normalizacao soma` dá o clássico.)

4. **Escore principal: JSD** contra o fundo de Robinson & Robinson, suavizada
   em janela (Capra & Singh 2007). Ela sabe que uma coluna toda de triptofano
   informa mais que uma coluna toda de leucina. Entropia, conservação
   normalizada e KL saem também, por coluna, no CSV.

5. **Dois testes estatísticos.** Mann-Whitney (aproximação normal) é grosseiro
   com bolso de 8 resíduos; o teste de permutação é exato por construção e
   determinístico (semente). O Cliff δ aparece junto porque tamanho de efeito
   e significância não são a mesma pergunta.

6. **Distância entre átomos pesados**, altLoc principal, primeiro MODEL. É a
   convenção de contato, e os modelos de cofold não têm hidrogênio.

## Como ler o resultado

- **`Neff/resíduo < 1`** → MSA raso. É a mesma faixa em que o cofold erra mais:
  conservação e estrutura são, as duas, pouco confiáveis pelo mesmo motivo.
  Não trate uma como validação da outra.
- **`coluna_confiavel = 0`** → coluna rasa. "100% conservada" com três
  sequências não-gap é acidente.
- **`identidade query × estrutura < 90%`** → aquele MSA provavelmente não é
  dessa cadeia. Todo o mapeamento abaixo disso é ficção. É o erro mais fácil
  de cometer quando se tem um MSA por entidade.
- **Bolso de interface** (≥ 2 cadeias poliméricas a 5 Å) → conservação baixa é
  **esperada**: a interface induzida por um PROTAC é neo-interface, não sítio
  ancestral. Num bolso de uma cadeia só, conservação baixa é suspeita — o
  modelo provavelmente encaixou o ligante num bolso de superfície qualquer.

Isto **não** mede a qualidade do modelo (pLDDT/ipTM/PAE: ver
`select_ternary_pose.py`) nem decide se o ternário existe. Conservação do
bolso é evidência independente, para somar.

## Saídas

| arquivo | conteúdo |
|---|---|
| `relatorio.md` | o relatório legível, com os vereditos |
| `resumo.json` | tudo em forma de máquina |
| `contatos_5A_residuos.csv` | a lista do bolso: um item por resíduo/cofator/íon/água, com categoria e distância |
| `contatos_5A_atomos.csv` | par de átomos a par de átomos |
| `bolso.pml` | seleção para conferir no PyMOL |
| `conservacao_colunas__*.csv` | por coluna: entropia, JSD, KL, ocupação, Neff |
| `pesos_henikoff__*.csv` | peso de cada sequência, com id e identidade com a query |
| `residuos__*.csv` | por resíduo da cadeia: conservação + se está no bolso |
| `comparacao_ferramentas.csv` | concordância entre ferramentas (ρ de Spearman) |
