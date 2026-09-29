#!/usr/bin/env python
"""
O preparo e a triagem da E3 existem, e a fase 4 conseguiria rodar com eles?

Roda no env `mdtools`. Segundos, e não escreve nada.

    python check_wp1.py --screening ~/.../screening --prep ~/.../prep \\
                        --e3 VHL CRBN --anchors-sdf <chemspace_anchors.sdf>

Por que ele importa as funções do wp1_select_recruiter
-----------------------------------------------------
A pergunta não é "existe uma pasta VHL". É "o `wp1_select_recruiter.py`
encontraria o que precisa". São coisas diferentes: aquele script escolhe o CSV
por `*round2*.csv` com fallback para `*.csv`, pega o PRIMEIRO subdiretório que
começa com o nome da E3 em ordem alfabética, e resolve cada pose por uma lista
de seis padrões de glob.

Reimplementar essa busca aqui daria um check que passa e uma fase 4 que falha —
que é exatamente o erro que este projeto já cometeu uma vez, com um teste que
comparava um arquivo com ele mesmo e dizia OK para o que a ferramenta recusava.
Então aqui as funções são IMPORTADAS, não reescritas: se o consumidor mudar, o
check muda com ele.

O que ele verifica, na ordem em que a fase 4 precisa
---------------------------------------------------
    1. CSV de triagem, e as colunas de id e score que o leitor reconhece
    2. pasta de preparo, e o `*_receptor.pdb` dentro dela
    3. as CADEIAS que o receptor realmente tem (o config pedir 'A' num arquivo
       que só tem 'C' já custou 104 configs errados neste projeto)
    4. o ligante de referência — sem ele não há revalidação do redocking
    5. as POSES dos primeiros ligantes do ranking, resolvidas pelo mesmo glob
       do consumidor. É aqui que a fase 4 descarta em silêncio.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from wp1_select_recruiter import achar_pose, ler_ranking  # noqa: E402


def cadeias_e_atomos(pdb: Path):
    cads, n = {}, 0
    for l in pdb.read_text(errors="ignore").splitlines():
        if l.startswith("ATOM"):
            n += 1
            cads[l[21]] = cads.get(l[21], 0) + 1
    return cads, n


def arvore(d: Path, prof: int = 2, limite: int = 25):
    """O que existe de verdade, quando o esperado não está lá."""
    if not d.exists():
        return [f"    {d} não existe"]
    linhas, n = [], 0
    for p in sorted(d.rglob("*")):
        rel = p.relative_to(d)
        if len(rel.parts) > prof:
            continue
        linhas.append(f"    {rel}{'/' if p.is_dir() else ''}")
        n += 1
        if n >= limite:
            linhas.append(f"    ... (cortado em {limite})")
            break
    return linhas or [f"    {d} está vazio"]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--screening", type=Path, required=True)
    ap.add_argument("--prep", type=Path, required=True)
    ap.add_argument("--e3", nargs="*", default=["VHL", "CRBN"])
    ap.add_argument("--anchors-sdf", type=Path)
    ap.add_argument("--top-n", type=int, default=10,
                    help="quantos ligantes do topo do ranking testar a pose")
    args = ap.parse_args()

    screening = args.screening.expanduser()
    prep = args.prep.expanduser()

    print("=" * 66)
    print("Check do WP1 — o que a fase 4 encontraria")
    print("=" * 66)
    print(f"  screening: {screening}")
    print(f"  prep:      {prep}")
    if not screening.exists():
        print(f"\n  [BLOQUEIO] a pasta de triagem não existe.")
        print("  O que existe no nível acima:")
        for l in arvore(screening.parent, 1):
            print(l)
    if not prep.exists():
        print(f"\n  [BLOQUEIO] a pasta de preparo não existe.")
        for l in arvore(prep.parent, 1):
            print(l)

    if args.anchors_sdf:
        a = args.anchors_sdf.expanduser()
        if a.exists():
            n = a.read_text(errors="ignore").count("$$$$")
            print(f"  âncoras Chemspace: {a.name} ({n} moléculas)")
        else:
            print(f"  [aviso] âncoras Chemspace ausentes ({a}) — a fase 4 roda,"
                  f" mas sem identificar o recrutador no catálogo")

    veredito = {}
    for e3 in args.e3:
        print(f"\n{'=' * 66}\n{e3}\n{'=' * 66}")
        faltas = []

        # --- 1. triagem ----------------------------------------------------
        dir_e3 = screening / e3
        csvs = sorted(dir_e3.glob("*round2*.csv")) or sorted(dir_e3.glob("*.csv"))
        if not csvs:
            print(f"  [BLOQUEIO] nenhum CSV de triagem em {dir_e3}")
            faltas.append("triagem (docking da biblioteca)")
            print("  o que existe por perto:")
            for l in arvore(dir_e3 if dir_e3.exists() else screening, 2):
                print(l)
            df = None
        else:
            escolhido = csvs[0]
            rodada2 = "round2" in escolhido.name
            try:
                df, col_id, col_sc = ler_ranking(escolhido)
                print(f"  [ok] triagem: {escolhido.name} — {len(df)} ligantes, "
                      f"id='{col_id}', score='{col_sc}'")
                if not rodada2:
                    print("       (é o CSV de fallback, não um *round2*: a fase 4"
                          " usa este mesmo)")
                if len(csvs) > 1:
                    print(f"       {len(csvs)} CSVs na pasta; a fase 4 usa o "
                          f"primeiro em ordem alfabética")
            except Exception as e:
                print(f"  [BLOQUEIO] o CSV existe mas o leitor falhou: "
                      f"{str(e)[:90]}")
                faltas.append("CSV de triagem legível")
                df = None

        # --- 2 e 3. preparo, receptor, cadeias -----------------------------
        subdirs = sorted([d for d in prep.iterdir()
                          if d.is_dir() and d.name.startswith(e3)]) \
            if prep.exists() else []
        rec_pdb = None
        if not subdirs:
            print(f"  [BLOQUEIO] nenhuma pasta {e3}* em {prep}")
            faltas.append(f"preparo da {e3} (receptor limpo)")
            if prep.exists():
                print("  as pastas que existem:")
                for d in sorted(p.name for p in prep.iterdir() if p.is_dir()):
                    print(f"    {d}/")
        else:
            usada = subdirs[0]
            print(f"  [ok] preparo: {usada.name}"
                  + (f"  ({len(subdirs)} pastas {e3}*; a fase 4 usa esta)"
                     if len(subdirs) > 1 else ""))
            recs = sorted(usada.glob("*_receptor.pdb"))
            if not recs:
                print(f"  [BLOQUEIO] nenhum *_receptor.pdb em {usada}")
                faltas.append("receptor preparado (*_receptor.pdb)")
                print("  o que existe lá:")
                for l in arvore(usada, 1):
                    print(l)
            else:
                rec_pdb = recs[0]
                cads, n_at = cadeias_e_atomos(rec_pdb)
                print(f"  [ok] receptor: {rec_pdb.name} ({n_at} átomos ATOM)")
                print(f"       cadeias presentes: "
                      + ", ".join(f"'{c}' ({v} átomos)"
                                  for c, v in sorted(cads.items())))
                if len(cads) > 1:
                    print("       [ATENÇÃO] mais de uma cadeia: o config do "
                          "PRosettaC precisa dizer QUAL, e o run_prosettac.sh")
                    print("       só corrige automaticamente quando há uma só.")

            # --- 4. ligante de referência ---------------------------------
            refs = sorted(usada.glob("*ref_ligand*"))
            if refs:
                tipos = {p.suffix for p in refs}
                print(f"  [ok] ligante de referência: "
                      + ", ".join(p.name for p in refs[:3]))
                if ".sdf" not in tipos:
                    print("       (só .pdb — a fase 4 converte para .sdf com o "
                          "obabel antes de revalidar)")
            else:
                print(f"  [aviso] nenhum *ref_ligand* em {usada}: a revalidação "
                      f"do redocking é pulada")
            poses_redock = sorted(usada.glob("redock*.pdbqt"))
            print(f"  {'[ok]' if poses_redock else '[aviso]'} poses de redocking:"
                  f" {len(poses_redock)}"
                  + ("" if poses_redock else " — sem elas o portão de RMSD ≤ 2 Å"
                                             " não é revalidado"))

        # --- 5. as poses dos ligantes do topo ------------------------------
        if df is not None and len(df):
            achadas, exemplos = 0, []
            for _, linha in df.head(args.top_n).iterrows():
                lig = str(linha[col_id])
                p = achar_pose(dir_e3, lig) or (achar_pose(prep, lig)
                                                if prep.exists() else None)
                if p is not None:
                    achadas += 1
                    if len(exemplos) < 2:
                        exemplos.append(f"{lig} -> {p.name}")
                elif len(exemplos) < 2:
                    exemplos.append(f"{lig} -> NÃO ACHEI")
            marca = "[ok]" if achadas == args.top_n else (
                "[BLOQUEIO]" if achadas == 0 else "[aviso]")
            print(f"  {marca} poses resolvidas: {achadas} de {args.top_n} "
                  f"do topo do ranking")
            for e in exemplos:
                print(f"       {e}")
            if achadas == 0:
                faltas.append("poses do docking (a fase 4 descartaria tudo)")
                print("       A fase 4 procura nestes padrões, nesta ordem:")
                print("         **/round2/<id>/*.pdbqt, **/round2/<id>*.pdbqt,")
                print("         **/<id>/run*.pdbqt, **/<id>*.pdbqt,")
                print("         **/<id>/*.sdf, **/<id>*.sdf")

        veredito[e3] = faltas

    print(f"\n{'=' * 66}")
    print("VEREDITO")
    print("=" * 66)
    pronta = [e for e, f in veredito.items() if not f]
    for e3, faltas in veredito.items():
        if not faltas:
            print(f"  {e3}: a fase 4 RODA com o que está no disco")
        else:
            print(f"  {e3}: falta — " + "; ".join(faltas))
    print()
    if pronta:
        # NÃO sugerir "use a que está pronta": qual E3 usar é decisão de
        # projeto, não de conveniência de arquivo. A CRBN estava pronta e não
        # produziu geometria ternária — trocar para a VHL foi justamente a
        # decisão tomada depois disso.
        print(f"  Prontas no disco: {', '.join(pronta)}")
        faltando = [e for e in veredito if veredito[e]]
        if faltando:
            print(f"  NÃO prontas: {', '.join(faltando)} — se a E3 do seu track")
            print(f"  está aqui, o que falta é a etapa manual do WP1 para ela,")
            print(f"  listada acima item por item.")
        print(f"\n  Para lançar o track de uma E3 pronta:")
        print(f"    WP1_E3_LIST=\"<E3>\" em config/pipeline_vhl.conf, e --from 4")
    else:
        print("  Nenhuma E3 está pronta. O que falta é a etapa manual do WP1:")
        print("    baixar o PDB (VHL: 6GFZ/6GFY), limpar em ChimeraX, extrair o")
        print("    ligante co-cristalizado, preparar o receptor em PDBQT e docar")
        print("    a biblioteca Chemspace de âncoras contra o sítio.")
        print("  Diga-me o que apareceu acima e eu escrevo essa etapa como fase")
        print("  do driver — ela é a mesma para VHL e CRBN, só mudam o PDB e o")
        print("  ligante de referência.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
