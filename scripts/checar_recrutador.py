#!/usr/bin/env python
"""
O recrutador é quimicamente compatível com a E3 que o projeto declara?

Roda no env `mdtools` (RDKit). Instantâneo.

    python checar_recrutador.py --e3 VHL --sdf recruiter_VHL_ligand_087.sdf
    python checar_recrutador.py --e3 VHL --smiles "O=C1CC..."
    python checar_recrutador.py --e3 VHL --protac-smi protac.smi

Por que esta verificação não existia, e o que custou
---------------------------------------------------
O WP1 escolheu o recrutador por SCORE DE DOCKING no sítio da E3, mais
enterramento. Nenhuma das duas medidas sabe que E3 a molécula de fato recruta:
uma biblioteca de blocos de PROTAC (Chemspace) contém ligantes de CRBN
(análogos de talidomida) E de VHL, e um análogo de talidomida pode pontuar bem
encaixado no bolso da VHL sem ter razão nenhuma para se ligar nela.

Foi o que aconteceu. O `protac.smi` do candidato SC0013__WH023 — o PROTAC que
consumiu sete dias de PRosettaC no "track da VHL" — carrega glutarimida +
ftalimida fluorada, ou seja, uma talidomida: o recrutador canônico da CRBN.
Nenhuma das três assinaturas de ligante de VHL está presente. A corrida inteira
dockeou a proteína VHL contra um ligante de CRBN.

E isso explica o resultado sem precisar de mais nada: um "recrutador" que não
se liga à E3 não produz interface cooperativa, e 103 clusters com mediana de
41 Å é o que a geometria devolve quando a única coisa que prende as duas
proteínas é a restrição de distância do linker.

Calibração das evidências
-------------------------
A identificação POSITIVA é a forte: glutarimida ligada a ftalimida é
talidomida, e talidomida liga CRBN. Isso basta.

A ausência das assinaturas de VHL é evidência forte mas não absoluta: a
esmagadora maioria dos ligantes de VHL deriva do VH032 e tem
(2S,4R)-4-hidroxiprolina, mas existem exceções na literatura. Por isso o script
reporta os dois lados em vez de decidir por ausência.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

# SMARTS por E3. A chave é a E3 que a molécula RECRUTA.
ASSINATURAS = {
    "CRBN": [
        ("glutarimida",   "O=C1CCC(C(=O)N1)N",  True),
        ("ftalimida",     "O=C1c2ccccc2C(=O)N1", False),
        ("isoindolinona", "O=C1Cc2ccccc2N1",     False),
    ],
    "VHL": [
        ("4-hidroxiprolina",     "OC1CC(NC1)C(=O)N",            True),
        ("terc-leucina amida",   "CC(C)(C)C(NC(=O)C)C(=O)N",    False),
        ("tiazol 4-metil",       "Cc1scnc1",                    False),
        ("benzilamida terminal", "c1ccccc1CNC(=O)",             False),
    ],
}
# o terceiro campo diz se a assinatura é SUFICIENTE sozinha para identificar


def motivos(mol):
    achados = {}
    for e3, lista in ASSINATURAS.items():
        achados[e3] = [(nome, suf) for nome, sma, suf in lista
                       if (q := Chem.MolFromSmarts(sma)) is not None
                       and mol.HasSubstructMatch(q)]
    return achados


def ler_biblioteca(caminho: Path, coluna: str | None):
    """[(id, mol)] de um .sdf, .smi/.txt, .csv ou DIRETÓRIO de .sdf."""
    itens, falhas = [], 0
    # O prep/ligands guarda um .sdf por ligante (ligand_087.sdf), não um
    # arquivo único. Varrer o diretório é o que responde "a biblioteca tem
    # ligante de VHL?" — pedir ao operador um `for` no shell daria uma
    # contagem por arquivo em vez da contagem agregada, que é a pergunta.
    if caminho.is_dir():
        for f in sorted(caminho.glob("*.sdf")):
            m = next(iter(Chem.SDMolSupplier(str(f), removeHs=True)), None)
            if m is None:
                falhas += 1
                continue
            itens.append((f.stem, m))
        return itens, falhas
    suf = caminho.suffix.lower()
    if suf in (".sdf", ".sd"):
        for i, m in enumerate(Chem.SDMolSupplier(str(caminho), removeHs=True)):
            if m is None:
                falhas += 1
                continue
            nome = m.GetProp("_Name") if m.HasProp("_Name") else f"mol{i}"
            itens.append((nome or f"mol{i}", m))
        return itens, falhas
    if suf == ".csv":
        import csv as _csv
        with open(caminho, newline="") as fh:
            leitor = _csv.DictReader(fh)
            cols = leitor.fieldnames or []
            col = coluna or next(
                (c for c in cols
                 if c and c.lower() in ("smiles", "canonical_smiles",
                                        "protac_smiles", "smi")), None)
            if col is None:
                raise SystemExit(
                    f"não achei coluna de SMILES em {caminho.name}.\n"
                    f"  colunas: {cols}\n  use --coluna <nome>")
            idc = next((c for c in cols if c and "id" in c.lower()), None)
            for i, linha in enumerate(leitor):
                m = Chem.MolFromSmiles((linha.get(col) or "").strip())
                if m is None:
                    falhas += 1
                    continue
                itens.append((linha.get(idc) or f"linha{i}", m))
        return itens, falhas
    # .smi / .txt: SMILES no primeiro campo, nome no segundo se houver
    for i, linha in enumerate(caminho.read_text(errors="ignore").splitlines()):
        campos = linha.split()
        if not campos:
            continue
        m = Chem.MolFromSmiles(campos[0])
        if m is None:
            falhas += 1
            continue
        itens.append((campos[1] if len(campos) > 1 else f"linha{i}", m))
    return itens, falhas


def ranking(scores: Path, col_id: str, col_score: str):
    """{id: (posição, score)} — posição 1 é o melhor (score mais negativo)."""
    import csv as _csv
    linhas = []
    with open(scores) as fh:
        leitor = _csv.DictReader(fh)
        if col_id not in (leitor.fieldnames or []) or \
                col_score not in (leitor.fieldnames or []):
            raise SystemExit(
                f"{scores.name} não tem '{col_id}'/'{col_score}'.\n"
                f"  colunas: {leitor.fieldnames}\n"
                f"  use --coluna-id e --coluna-score")
        for l in leitor:
            try:
                linhas.append((l[col_id], float(l[col_score])))
            except (TypeError, ValueError):
                continue
    linhas.sort(key=lambda t: t[1])      # mais negativo = melhor
    return {i: (pos, s) for pos, (i, s) in enumerate(linhas, 1)}, len(linhas)


def varrer(caminho: Path, e3: str, coluna: str | None,
           scores: Path | None = None, col_id: str = "ligand_id",
           col_score: str = "best_score"):
    """A pergunta decisiva: a biblioteca TEM ligante da E3 que queremos?

    Refazer a triagem com filtro de quimiotipo só funciona se houver o que
    filtrar. Se a biblioteca não contém nenhum ligante de VHL, nenhuma
    quantidade de triagem nova encontra um — e a única saída é usar um
    recrutador conhecido da literatura. Uma busca de segundos separa as duas
    situações, e sem ela a escolha entre os caminhos seria palpite.
    """
    if not caminho.exists():
        raise SystemExit(f"não achei {caminho}")
    itens, falhas = ler_biblioteca(caminho, coluna)
    print(f"biblioteca: {caminho.name}")
    print(f"  {len(itens)} moléculas lidas"
          + (f", {falhas} não parsearam" if falhas else ""))
    if not itens:
        raise SystemExit("  nenhuma molécula legível")

    conta = {a: [] for a in ASSINATURAS}
    for nome, m in itens:
        for alvo, lista in motivos(m).items():
            if any(s for _, s in lista):
                conta[alvo].append(nome)
    print()
    for alvo in sorted(conta):
        n = len(conta[alvo])
        pct = 100 * n / len(itens)
        print(f"  quimiotipo de {alvo:<5}: {n:>6} de {len(itens)} ({pct:.1f}%)")
        if n:
            print(f"      ex.: {', '.join(str(x) for x in conta[alvo][:4])}")
    # --- o quimiotipo CONTRA o ranking da triagem -------------------------
    # A pergunta de método: existindo ligante do quimiotipo certo, o score o
    # teria encontrado? Se ele existe e ficou no fim do ranking, o problema
    # não é a biblioteca — é a função de score não discriminar o bolso, e isso
    # é um achado, não um contratempo.
    if scores and scores.exists():
        pos, total = ranking(scores.expanduser(), col_id, col_score)
        print(f"\n  RANKING da triagem ({scores.name}, {total} ligantes):")
        for alvo in sorted(conta):
            com_pos = sorted(((pos[n][0], pos[n][1], n)
                              for n in conta[alvo] if n in pos))
            if not com_pos:
                print(f"    {alvo:<5}: nenhum dos quimiotipos está no CSV")
                continue
            # Com o NOME e não só a posição: "o melhor de quimiotipo VHL está
            # em #4" é interessante, mas o passo seguinte precisa saber QUAL
            # ligante é esse para olhar o exit vector dele.
            print(f"    {alvo:<5}:")
            for pp, ss, nn in com_pos[:6]:
                print(f"        #{pp:<4} {nn:<14} {ss:>7.1f}")

    print()
    if conta.get(e3):
        print(f"  A biblioteca TEM {len(conta[e3])} candidato(s) com quimiotipo"
              f" de {e3}.")
        print(f"  Refazer a triagem restringindo a eles é viável.")
        return 0
    print(f"  *** A biblioteca NÃO TEM nenhuma molécula com quimiotipo"
          f" de {e3}.")
    print(f"  *** Refazer a triagem não resolve: não há o que a triagem possa")
    print(f"  *** encontrar. O caminho é usar um recrutador conhecido de {e3}")
    print(f"  *** (para VHL, o scaffold VH032 / (2S,4R)-4-hidroxiprolina, que")
    print(f"  *** tem estrutura cristalográfica e dezenas de PROTACs")
    print(f"  *** publicados), em vez de procurar um.")
    return 2


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--e3", required=True, help="a E3 que o projeto declara")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--sdf", type=Path)
    g.add_argument("--smiles")
    g.add_argument("--protac-smi", type=Path,
                   help="o protac.smi do job: confere a molécula inteira")
    g.add_argument("--biblioteca", type=Path,
                   help="varre uma biblioteca e conta quantas moléculas têm "
                        "quimiotipo de cada E3. Aceita .sdf, .smi, .csv ou um "
                        "DIRETÓRIO de .sdf (é assim que o prep/ligands guarda)")
    ap.add_argument("--coluna", default=None,
                    help="nome da coluna de SMILES, para --biblioteca em CSV")
    ap.add_argument("--scores", type=Path, default=None,
                    help="CSV da triagem (ligand_id, best_score). Junta o "
                         "quimiotipo com o RANKING: responde se o score teria "
                         "achado o ligante certo, caso ele exista")
    ap.add_argument("--coluna-id", default="ligand_id")
    ap.add_argument("--coluna-score", default="best_score")
    args = ap.parse_args()

    if args.biblioteca:
        return varrer(args.biblioteca.expanduser(), args.e3.upper(),
                      args.coluna, args.scores, args.coluna_id,
                      args.coluna_score)

    if args.sdf:
        mol = next(iter(Chem.SDMolSupplier(str(args.sdf.expanduser()),
                                           removeHs=True)), None)
        fonte = args.sdf.name
    elif args.protac_smi:
        txt = args.protac_smi.expanduser().read_text().split()
        mol = Chem.MolFromSmiles(txt[0]) if txt else None
        fonte = args.protac_smi.name
    else:
        mol = Chem.MolFromSmiles(args.smiles)
        fonte = "--smiles"
    if mol is None:
        raise SystemExit(f"não consegui ler a molécula de {fonte}")

    e3 = args.e3.upper()
    print(f"molécula: {fonte} — {mol.GetNumAtoms()} átomos pesados")
    print(f"E3 declarada pelo projeto: {e3}")
    ach = motivos(mol)
    print()
    for alvo, lista in ach.items():
        if lista:
            print(f"  assinaturas de {alvo}: "
                  + ", ".join(n + (" [suficiente]" if s else "")
                              for n, s in lista))
        else:
            print(f"  assinaturas de {alvo}: NENHUMA")

    identificadas = {a for a, l in ach.items() if any(s for _, s in l)}
    print()
    if e3 in identificadas and len(identificadas) == 1:
        print(f"  OK: a molécula tem assinatura suficiente de {e3} e de mais"
              f" nenhuma E3.")
        return 0
    outras = identificadas - {e3}
    if outras:
        print(f"  *** INCOMPATÍVEL: assinatura suficiente de"
              f" {', '.join(sorted(outras))},")
        print(f"  *** e o projeto declara {e3}.")
        print(f"  ***")
        print(f"  *** Score de docking não sabe qual E3 a molécula recruta. Uma")
        print(f"  *** biblioteca de blocos de PROTAC tem ligantes das duas, e um")
        print(f"  *** análogo de talidomida encaixa no bolso da VHL sem ter razão")
        print(f"  *** para se ligar nela. Seguir daqui modela um complexo que a")
        print(f"  *** química não sustenta.")
        return 2
    print(f"  *** SEM ASSINATURA de {e3}, e sem assinatura de outra E3.")
    print(f"  *** Isto não prova que a molécula não se liga à {e3} — há")
    print(f"  *** ligantes fora dos quimiotipos canônicos — mas exige")
    print(f"  *** justificativa explícita antes de gastar dias de máquina.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
