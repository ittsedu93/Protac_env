#!/usr/bin/env python
"""
Deixa os .sdf dos heads legíveis para o PRosettaC — e PROVA que ficaram.

Roda no env `mdtools` (RDKit). Segundos.

    python prosettac_fix_heads.py --dir <run_dir>             # diagnostica
    python prosettac_fix_heads.py --dir <run_dir> --aplicar    # conserta

O problema que isto resolve
---------------------------
O PRosettaC protona o head (gera `<head>_H.sdf`) e depois remapeia a âncora
casando o head ORIGINAL contra a versão protonada:

    OldSdf = Chem.SDMolSupplier(old_sdf, sanitize=False)[0]
    try:    Chem.SanitizeMol(OldSdf)
    except: pass                       # <- a falha é engolida
    NewMatch = NewSdf.GetSubstructMatch(OldSdf)
    if len(NewMatch) == 0: return -1

Saindo o match vazio, a função devolve -1 e o -1 só aparece cinquenta linhas
adiante, como índice negativo:

    anchor_b = HeadB.GetConformer().GetAtomPosition(Anchors[1])
    OverflowError: can't convert negative value to unsigned int

O head fica nesse estado quando o .sdf traz a contabilidade de hidrogênio
presa — ligações aromáticas (tipo 4) que o RDKit não consegue kekulizar, ou o
campo de valência do bloco de átomos fixando `noImplicit` com zero H. O
sintoma no SMILES são os colchetes:

    Oc1[c]c([C][C]Nc2nnn[nH]2)[c]c(Oc2[c][c][c][c][c]2)[c]1   <- head
    Oc1cc(CCNc2nnn[nH]2)cc(Oc2ccccc2)c1                       <- _H.sdf

Um átomo de consulta que declara 0 H não casa com o mesmo átomo carregando 1 H
na versão protonada. Cada átomo assim mata o casamento inteiro.

Por que a primeira versão deste script mentiu
---------------------------------------------
Ela comparava o arquivo com uma RELEITURA DELE MESMO. As duas leituras herdam
os mesmos campos tortos, então o match casava sempre — e o script dizia "OK"
para um head que o PRosettaC recusava na linha seguinte. O alvo do casamento é
a versão PROTONADA, e é contra ela que este script testa: contra o `_H.sdf`
real quando ele existe, e contra `Chem.AddHs(head)` — que é o que o PRosettaC
vai gerar — sempre.

Nada aqui é heurístico: cada candidato a conserto é ESCRITO, RELIDO do jeito
que o PRosettaC lê, casado contra a versão protonada e conferido átomo por
átomo contra as coordenadas originais. Só sobrevive o que passa nos três.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

DESVIO_MAX_A = 1e-3   # o conserto mexe em ligação e H, nunca em coordenada


# --------------------------------------------------------------------------
# 1. o que o PRosettaC vê
# --------------------------------------------------------------------------
def como_o_prosettac_le(caminho: Path):
    """A mesma leitura do translate_anchors, com o mesmo except silencioso."""
    mol = Chem.SDMolSupplier(str(caminho), sanitize=False)[0]
    if mol is None:
        return None, "o RDKit não leu o arquivo nem sem sanitizar"
    erro = None
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        erro = str(e).split("\n")[0][:110]
    return mol, erro


def atomos_presos(mol: Chem.Mol) -> list[str]:
    """Átomos que declaram 0 H com `noImplicit` — os que matam o casamento."""
    presos = []
    for a in mol.GetAtoms():
        if a.GetAtomicNum() == 1:
            continue
        if a.GetNoImplicit() and a.GetNumExplicitHs() == 0:
            presos.append(f"{a.GetSymbol()}{a.GetIdx()}")
    return presos


def protonado_de_referencia(mol: Chem.Mol):
    """A versão protonada, que é o ALVO do casamento dentro do PRosettaC."""
    try:
        m = Chem.Mol(mol)
        Chem.SanitizeMol(m)
        return Chem.AddHs(m, addCoords=bool(m.GetNumConformers()))
    except Exception:
        return None


def casa(caminho: Path, alvo: Chem.Mol | None = None):
    """(n_atomos_casados, erro_de_sanitizacao, smiles_como_lido, alvo_usado).

    `alvo=None` faz o alvo ser a protonação do próprio arquivo — exatamente o
    `_H.sdf` que o PRosettaC gera a partir dele.
    """
    velho, erro = como_o_prosettac_le(caminho)
    if velho is None:
        return 0, erro, None, None
    smi = Chem.MolToSmiles(velho)
    novo = alvo if alvo is not None else protonado_de_referencia(velho)
    if novo is None:
        return 0, erro or "não consegui protonar para comparar", smi, None
    return len(novo.GetSubstructMatch(velho)), erro, smi, novo


# --------------------------------------------------------------------------
# 2. candidatos a conserto, do menos invasivo ao mais
# --------------------------------------------------------------------------
def _sem_contabilidade_de_H(bruto: Chem.Mol, so_carbono: bool):
    """Solta a contagem de H, devolvendo-a ao modelo de valência do RDKit.

    É o conserto do caso real: o arquivo fixa 0 H em átomos que têm H, e o
    RDKit obedece. Zerando `noImplicit` ele volta a contar sozinho. Em
    `so_carbono` os N e O ficam como estão — num tetrazol o [nH] é a única
    coisa que permite kekulizar, e apagá-lo troca um problema por outro.
    """
    mol = Chem.RWMol(bruto)
    for a in mol.GetAtoms():
        if so_carbono and a.GetAtomicNum() != 6:
            continue
        a.SetNoImplicit(False)
        a.SetNumExplicitHs(0)
        a.SetNumRadicalElectrons(0)
    m = mol.GetMol()
    m.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(m)
    return m


def _com_NH_no_anel(bruto: Chem.Mol):
    """Atribui o NH que falta a um N aromático de anel, testando cada um.

    O caso é o tetrazol da warhead: ligações marcadas aromáticas, nenhum N com
    H, e o RDKit não tem como adivinhar qual dos quatro N carrega o próton.
    Testa-se cada candidato e fica o primeiro que sanitiza — determinístico, e
    sem inventar: ou existe um tautômero válido, ou o arquivo está errado de
    outro jeito.

    Não dá para exigir GetIsAromatic() aqui: sem sanitizar o RDKit ainda não
    percebeu aromaticidade nos ÁTOMOS. O que é observável é símbolo, grau, H
    explícitos e pertencer a anel.
    """
    base = Chem.Mol(bruto)
    base.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(base)
    candidatos = [a.GetIdx() for a in base.GetAtoms()
                  if a.GetSymbol() == "N" and a.GetDegree() == 2
                  and a.GetNumExplicitHs() == 0 and a.IsInRing()]
    for idx in candidatos:
        mol = Chem.RWMol(base)
        at = mol.GetAtomWithIdx(idx)
        at.SetNumExplicitHs(1)
        at.SetNoImplicit(True)
        m = mol.GetMol()
        try:
            m.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(m)
        except Exception:
            continue
        yield m, (f"NH atribuído ao N {idx} "
                  f"(de {len(candidatos)} candidatos de anel)")


def candidatos(caminho: Path):
    """Gera (mol, nota). Quem escolhe é a verificação, não esta função."""
    bruto = Chem.SDMolSupplier(str(caminho), sanitize=False)[0]
    if bruto is None:
        return

    # (a) o arquivo já sanitiza por conta própria? reescrevê-lo a partir de uma
    #     leitura sanitizada normaliza aromaticidade e sai kekulizado.
    limpo = Chem.SDMolSupplier(str(caminho), sanitize=True)[0]
    if limpo is not None:
        yield limpo, "reescrito a partir de uma leitura sanitizada"

    # (b) e (c) soltar a contagem de H — o conserto do caso real
    for so_c in (False, True):
        try:
            yield (_sem_contabilidade_de_H(bruto, so_c),
                   "contagem de H devolvida ao RDKit"
                   + (" (só nos carbonos)" if so_c else ""))
        except Exception:
            pass

    # (d) o tetrazol sem NH: atribuir o próton e, depois, soltar os carbonos
    for m, nota in _com_NH_no_anel(bruto):
        yield m, nota
        try:
            yield _sem_contabilidade_de_H(m, True), nota + " + carbonos soltos"
        except Exception:
            pass


# --------------------------------------------------------------------------
# 3. verificação: escreve, relê como o PRosettaC, casa, confere coordenadas
# --------------------------------------------------------------------------
def coordenadas(mol: Chem.Mol):
    if not mol.GetNumConformers():
        return {}
    c = mol.GetConformer()
    return {i: tuple(round(v, 4) for v in (c.GetAtomPosition(i).x,
                                           c.GetAtomPosition(i).y,
                                           c.GetAtomPosition(i).z))
            for i, a in enumerate(mol.GetAtoms()) if a.GetAtomicNum() != 1}


def desvio_de_coordenadas(antes: dict, depois: dict) -> float:
    """Maior deslocamento de átomo pesado, em Å. Zero é o esperado."""
    if not antes or not depois:
        return 0.0
    comuns = set(antes) & set(depois)
    if len(comuns) != len(antes):
        return float("inf")   # mudou o número de átomos pesados: não serve
    return max(max(abs(a - b) for a, b in zip(antes[i], depois[i]))
               for i in comuns)


def experimenta(caminho: Path, mol: Chem.Mol, coords_originais: dict,
                temp: Path):
    """Escreve o candidato num arquivo temporário e testa tudo nele."""
    with Chem.SDWriter(str(temp)) as w:
        w.write(mol)
    n, erro, smi, _ = casa(temp)
    relido, _ = como_o_prosettac_le(temp)
    desvio = desvio_de_coordenadas(coords_originais,
                                   coordenadas(relido) if relido else {})
    presos = atomos_presos(relido) if relido else ["ilegível"]
    ok = (n > 0 and erro is None and not presos and desvio <= DESVIO_MAX_A)
    return ok, dict(n=n, erro=erro, smi=smi, desvio=desvio, presos=presos)


# --------------------------------------------------------------------------
def heads_do_config(cfg: Path, d: Path) -> list[Path]:
    for l in cfg.read_text().splitlines():
        if l.startswith("Heads: "):
            return [d / h for h in l.split(": ", 1)[1].split()]
    raise SystemExit(f"não achei a linha Heads: em {cfg}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--config", default="prosetta_config.txt")
    ap.add_argument("--aplicar", action="store_true")
    args = ap.parse_args()

    d = args.dir.expanduser().resolve()
    heads = heads_do_config(d / args.config, d)

    problemas = 0
    for h in heads:
        print(f"\n{h.name}")
        if not h.exists():
            print("  [ERRO] não existe")
            problemas += 1
            continue

        n, erro, smi, _ = casa(h)
        bruto, _ = como_o_prosettac_le(h)
        presos = atomos_presos(bruto) if bruto else []
        if smi:
            print(f"  como o PRosettaC lê:  {smi[:78]}")
        if erro:
            print(f"  sanitização FALHOU:   {erro}")
        if presos:
            print(f"  átomos com 0 H fixado: {len(presos)} "
                  f"({', '.join(presos[:8])}{'...' if len(presos) > 8 else ''})")

        # o _H.sdf real, quando existe, é o alvo verdadeiro — mais severo
        h_sdf = h.with_name(h.stem + "_H.sdf")
        if h_sdf.exists():
            alvo = Chem.SDMolSupplier(str(h_sdf), sanitize=True)[0]
            if alvo is not None:
                n_real, _, _, _ = casa(h, alvo)
                print(f"  contra {h_sdf.name} (o que o PRosettaC gerou): "
                      f"{n_real} átomos")
                n = min(n, n_real)

        if n > 0 and erro is None and not presos:
            print(f"  OK — o GetSubstructMatch casa {n} átomos contra a "
                  f"versão protonada")
            continue

        print("  [PROBLEMA] o casamento contra a versão protonada sai VAZIO.")
        print("             translate_anchors devolve -1, e o -1 vira o")
        print("             OverflowError cinquenta linhas adiante.")
        problemas += 1
        if not args.aplicar:
            continue

        coords = coordenadas(bruto) if bruto else {}
        temp = h.with_suffix(".sdf.tentativa")
        escolhido = None
        for mol, nota in candidatos(h):
            ok, r = experimenta(h, mol, coords, temp)
            estado = "serve" if ok else "não serve"
            motivo = ""
            if not ok:
                if r["n"] == 0:
                    motivo = "casamento vazio"
                elif r["erro"]:
                    motivo = f"não sanitiza: {r['erro'][:40]}"
                elif r["presos"]:
                    motivo = f"{len(r['presos'])} átomos ainda com 0 H fixado"
                else:
                    motivo = f"coordenadas mudaram {r['desvio']:.3f} Å"
                motivo = f" — {motivo}"
            print(f"    tentativa: {nota}: {estado}{motivo}")
            if ok:
                escolhido = (mol, nota, r)
                break
        temp.unlink(missing_ok=True)

        if escolhido is None:
            print("  [ERRO] nenhum conserto passou na verificação. O head não")
            print("         descreve uma molécula válida — reveja a origem do")
            print("         arquivo (o docking que o gerou).")
            continue

        mol, nota, r = escolhido
        bak = h.with_suffix(".sdf.bak")
        if not bak.exists():
            shutil.copy2(h, bak)
        with Chem.SDWriter(str(h)) as w:
            w.write(mol)
        print(f"  CORRIGIDO ({nota})")
        print(f"    agora o PRosettaC lê:  {r['smi'][:78]}")
        print(f"    casa {r['n']} átomos | coordenadas intactas "
              f"(desvio {r['desvio']:.4f} Å) | original em {bak.name}")

        # O _H.sdf velho nasceu do arquivo torto. Deixá-lo aqui é pedir para o
        # PRosettaC reaproveitá-lo e casar o head novo contra o lixo antigo.
        if h_sdf.exists():
            h_sdf.replace(h_sdf.with_suffix(".sdf.velho"))
            print(f"    {h_sdf.name} era do arquivo antigo: movido para "
                  f"{h_sdf.name}.velho, o PRosettaC gera outro")
        problemas -= 1

    print()
    if problemas and not args.aplicar:
        print(f"  {problemas} head(s) com problema. Para corrigir:")
        print(f"    python {Path(__file__).name} --dir {d} --aplicar")
        raise SystemExit(1)
    if problemas:
        raise SystemExit(f"  {problemas} head(s) seguem com problema")
    print("  todos os heads casam contra a versão protonada — "
          "translate_anchors não vai devolver -1")


if __name__ == "__main__":
    main()
