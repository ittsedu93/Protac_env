#!/usr/bin/env python
"""
Fase 10a — prepara o sistema de MD do complexo TERNÁRIO (nível iii).

Roda no env `mdtools` (RDKit). Segundos.

    python md_prepare_ternario.py --candidato SC0013__WH023 \\
                                  --outdir ~/pipeline_vhl/md_ternario
    python md_prepare_ternario.py --modelo <pose.pdb> --protac-smi <protac.smi> \\
                                  --id SC0013__WH023 --outdir ...

Nível (iii) é o que o WP3 do projeto chama de "the complete ternary complex":
E3 ligase + recrutador-linker-warhead + PCSK9. É a pergunta final do trabalho —
a interface ternária modelada pelo PRosettaC sobrevive 200 ns de solvente, ou
ela só existe porque o docking a construiu?

O que muda do nível (ii) para o (iii)
-------------------------------------
No nível (ii) as coordenadas do PROTAC eram CONSTRUÍDAS (embebidas com o
recrutador fixo na pose do docking), porque não havia estrutura do complexo.
Aqui elas já existem: o modelo do PRosettaC traz E3, alvo e PROTAC juntos, numa
geometria que o Rosetta refinou. O trabalho deixa de ser construir e passa a ser
SEPARAR — proteínas de um lado, ligante de outro — sem mexer em coordenada
nenhuma.

E por isso o risco também muda. O do nível (ii) era geometria inventada; o daqui
é suposição sobre o formato do arquivo. Este script não supõe: ele INSPECIONA o
modelo, imprime a composição, e para com o conteúdo real na tela quando o que
encontra não é o que o resto do pipeline precisa.

O que ele confere antes de escrever qualquer coisa
--------------------------------------------------
    duas cadeias de proteína, com contagem de resíduos plausível
    um ligante HETATM, contíguo, com o número de átomos pesados do PROTAC
    nenhum átomo do ligante coincidindo com átomo de proteína
    o SMILES casando com as coordenadas, átomo por átomo
    a carga formal do PROTAC, que é o que vai no `acpype -n`

A carga está aqui porque é o erro mais caro desta etapa: errada, ela não dá
mensagem nenhuma e contamina a trajetória inteira.

Por que o SMILES é obrigatório
------------------------------
O PDB do PRosettaC não guarda ordem de ligação nem carga. Entregar esse PDB
cru ao acpype faz o antechamber perceber as ligações só pela geometria — e uma
amida refinada pelo Rosetta pode sair como ligação simples, um anel aromático
como ciclo-hexeno. O resultado é uma parametrização plausível de uma molécula
que não é o PROTAC. Aqui as ordens de ligação vêm do SMILES
(AssignBondOrdersFromTemplate), e o PDB entregue ao acpype sai do RDKit COM
registros CONECT e COM hidrogênios — o mesmo caminho do nível (ii).

O CONTRATO com o md_run.sh
--------------------------
O md_run.sh lê deste md_sistema.json, pelo nome, as chaves:

    candidate_id   carga_formal   resname   lig_pdb   receptor_pdb

e espera um `protac.sdf` no diretório (o fix_receptor_for_md.py e o
split_chain_gaps.py o usam para medir distâncias ao ligante). Faltar qualquer
uma delas é a falha que este projeto já teve duas vezes: arquivo presente,
config correto, fiação faltando. Por isso elas são escritas aqui e CONFERIDAS
no fim, contra a lista de chaves que o md_run.sh de fato lê.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
# Batizar o resíduo do ligante é a MESMA operação do nível (ii), e este projeto
# já aprendeu o que acontece com duas implementações da mesma coisa: elas
# divergem. O nome do resíduo viaja do PDB para o .gro, para o make_ndx e para
# a análise — se os dois níveis batizarem diferente, a análise do (iii) procura
# um resname que não existe.
from md_prepare import LIG_RESNAME, nomear_residuo  # noqa: E402

RDLogger.DisableLog("rdApp.*")

# Abaixo disto não existe química: é o mesmo espaço ocupado duas vezes.
SOBREPOSICAO_A = 1.2
# Entre os dois, é choque — ruim, mas a minimização resolve.
CHOQUE_A = 2.0

# As chaves que o scripts/md_run.sh lê com `jq_ <chave>`. Mudar o md_run.sh sem
# mudar esta lista é pedir o mesmo erro de fiação de novo.
CHAVES_DO_MD_RUN = ("candidate_id", "carga_formal", "resname",
                    "lig_pdb", "receptor_pdb")

SOLVENTE_CRISTAL = {"HOH", "WAT", "NA", "CL", "SO4", "GOL", "EDO", "ZN", "MG"}


def ler_pdb(p: Path):
    """(ATOM por cadeia, HETATM por resname) — sem tocar em coordenadas."""
    prot, het = {}, {}
    for l in p.read_text(errors="ignore").splitlines():
        if l.startswith("ATOM"):
            prot.setdefault(l[21], []).append(l)
        elif l.startswith("HETATM"):
            rn = l[17:20].strip()
            if rn in SOLVENTE_CRISTAL:
                continue
            het.setdefault(rn, []).append(l)
    return prot, het


def residuos(linhas):
    return len({l[22:27] for l in linhas})


def coords(linhas):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in linhas])


def eh_hidrogenio(linha: str) -> bool:
    """Elemento da coluna 77-78; se vier vazia, cai no nome do átomo."""
    el = linha[76:78].strip().upper()
    if el:
        return el == "H"
    nome = linha[12:16].strip()
    return bool(re.match(r"^\d*H", nome))


def extensao(xyz: np.ndarray) -> float:
    """Maior aresta da caixa envolvente, em Å."""
    return float((xyz.max(axis=0) - xyz.min(axis=0)).max())


def custo_da_caixa(xyz_tudo: np.ndarray, padding_nm: float = 1.3):
    """Quantos átomos o sistema solvatado terá, e quanta RAM isso pede.

    A máquina é compartilhada. Dizer o tamanho ANTES de solvatar é mais
    educado que descobrir o custo quando o GROMACS já tomou a placa.
    """
    aresta_nm = extensao(xyz_tudo) / 10.0 + 2 * padding_nm
    volume = aresta_nm ** 3
    n_soluto = len(xyz_tudo)
    # 33,4 moléculas de água por nm3 a 310 K, 3 átomos cada; o soluto desloca
    # cerca de 1,2 nm3 por 100 átomos pesados
    volume_livre = max(volume - n_soluto * 0.012, 0.0)
    n_agua = int(volume_livre * 33.4)
    n_total = n_soluto + 3 * n_agua
    # ORDEM DE GRANDEZA, não medida: ~0,3 GB de base mais ~4 kB por átomo
    # entre coordenadas, forças, listas de vizinhos e a malha do PME. Serve
    # para distinguir "cabe" de "não cabe", não para reservar memória.
    ram_gb = 0.3 + n_total * 4.0e-6
    return {"aresta_nm": round(aresta_nm, 2), "n_atomos_soluto": n_soluto,
            "n_aguas_estimado": n_agua, "n_atomos_total_estimado": n_total,
            "ram_estimada_gb": round(ram_gb, 1)}


def ligante_com_ordens(lig_linhas, smiles: str):
    """Mol com as COORDENADAS do modelo e as ORDENS DE LIGAÇÃO do SMILES."""
    bloco = "\n".join(lig_linhas) + "\nEND\n"
    bruto = Chem.MolFromPDBBlock(bloco, removeHs=False, sanitize=False,
                                 proximityBonding=True)
    if bruto is None:
        raise SystemExit("  o RDKit não leu o bloco HETATM do modelo")
    bruto = Chem.RemoveHs(bruto, sanitize=False)
    modelo = Chem.MolFromSmiles(smiles)
    if modelo is None:
        raise SystemExit(f"  SMILES inválido: {smiles}")
    try:
        mol = AllChem.AssignBondOrdersFromTemplate(modelo, bruto)
    except Exception as e:
        raise SystemExit(
            f"  [NÃO CASOU] o SMILES do PROTAC não casa com as coordenadas do\n"
            f"  modelo ({e}).\n"
            f"  Sem esse casamento as ordens de ligação teriam de sair da\n"
            f"  geometria, e o antechamber erraria amidas e aromáticos — uma\n"
            f"  parametrização bonita de uma molécula que não é o PROTAC.\n"
            f"  Confira se o protac.smi é o deste candidato.")
    Chem.SanitizeMol(mol)
    # Hs explícitos com coordenada: é o que o acpype/antechamber precisa, e é o
    # mesmo caminho do nível (ii)
    mol = Chem.AddHs(mol, addCoords=True)
    return mol


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--modelo", type=Path,
                    help="PDB do ternário (ex.: Results/cluster1/*.pdb)")
    ap.add_argument("--candidato",
                    help="acha o melhor cluster deste candidato sozinho")
    ap.add_argument("--id", dest="cand_id",
                    help="candidate_id, quando se passa --modelo direto")
    ap.add_argument("--prosettac-dir", type=Path,
                    help="raiz dos jobs do PRosettaC (default: wp3/prosettac "
                         "sob $PIPELINE_OUT)")
    ap.add_argument("--protac-smi", type=Path,
                    help="protac.smi do candidato (obrigatório: dele saem as "
                         "ordens de ligação e a carga formal)")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--aceitar-divergencia", action="store_true",
                    help="seguir mesmo com o ligante diferente do SMILES")
    ap.add_argument("--ram-limite-gb", type=float, default=8.0,
                    help="recusa se o sistema solvatado passar disto (a "
                         "workstation é compartilhada; default 24)")
    args = ap.parse_args()

    # --- achar o modelo ---------------------------------------------------
    modelo = args.modelo
    cand_id = args.cand_id or args.candidato
    if modelo is None:
        if not args.candidato:
            raise SystemExit("passe --modelo ou --candidato")
        raiz = args.prosettac_dir
        if raiz is None:
            raise SystemExit("passe --prosettac-dir (ou --modelo direto)")
        base = raiz.expanduser() / args.candidato
        res = base / "Results"
        if not res.is_dir():
            raise SystemExit(f"não achei {res} — o PRosettaC terminou?")
        # o cluster com mais membros é o representativo; empate, o cluster1
        clusters = [(len(list(d.rglob("*.pdb"))), d) for d in res.iterdir()
                    if d.is_dir()]
        clusters.sort(key=lambda t: (-t[0], t[1].name))
        if not clusters or not clusters[0][0]:
            raise SystemExit(f"nenhum .pdb em {res}")
        n_mem, d = clusters[0]
        modelo = sorted(d.rglob("*.pdb"))[0]
        print(f"cluster escolhido: {d.name} ({n_mem} membros)")
        if not args.protac_smi:
            smi = base / "protac.smi"
            if smi.exists():
                args.protac_smi = smi
    modelo = modelo.expanduser().resolve()
    if not modelo.exists():
        raise SystemExit(f"não achei {modelo}")
    if not cand_id:
        cand_id = modelo.parent.parent.parent.name
    print(f"modelo: {modelo}")
    print(f"candidato: {cand_id}")

    if not args.protac_smi or not args.protac_smi.expanduser().exists():
        raise SystemExit(
            "\n*** falta o --protac-smi.\n"
            "*** Dele saem DUAS coisas que o PDB do PRosettaC não guarda: as\n"
            "*** ordens de ligação do PROTAC e a carga formal que vai no\n"
            "*** `acpype -n`. Sem elas o antechamber adivinha pela geometria,\n"
            "*** e carga errada não dá mensagem — contamina a trajetória\n"
            "*** inteira. O arquivo é o wp3/prosettac/<candidato>/protac.smi.")
    smiles = args.protac_smi.expanduser().read_text().split()[0]

    # --- inspecionar, antes de qualquer decisão ---------------------------
    prot, het = ler_pdb(modelo)
    print("\ncomposição do arquivo:")
    for c, linhas in sorted(prot.items()):
        print(f"  cadeia {c}: {residuos(linhas):>4} resíduos, "
              f"{len(linhas):>5} átomos")
    for rn, linhas in sorted(het.items()):
        print(f"  HETATM {rn}: {len(linhas)} átomos")
    if not prot:
        raise SystemExit("\nnenhum registro ATOM: isto não é um complexo")
    if len(prot) < 2:
        raise SystemExit(
            f"\nsó uma cadeia de proteína ({list(prot)}). O nível (iii) exige "
            f"E3 E alvo —\no que está aqui é um complexo binário.")
    if not het:
        raise SystemExit(
            "\nnenhum HETATM: o PROTAC não está no arquivo, ou veio como ATOM."
            "\nO conteúdo acima é o que o arquivo tem; sem o ligante não há o"
            " que parametrizar.")

    # --- o ligante --------------------------------------------------------
    if len(het) > 1:
        print(f"\n  [ATENÇÃO] {len(het)} resíduos HETATM distintos: "
              f"{', '.join(het)}")
        print("  O PRosettaC pode escrever o PROTAC como DUAS cabeças + linker."
              "\n  Elas serão unidas num só ligante — confira a contagem de "
              "átomos abaixo.")
    lig_linhas = [l for linhas in het.values() for l in linhas]
    n_lig = len(lig_linhas)

    # nenhum átomo do ligante pode coincidir com proteína: isso denunciaria
    # um arquivo com modelos sobrepostos, que o GROMACS aceitaria e explodiria
    xyz_lig = coords(lig_linhas)
    xyz_prot = np.vstack([coords(l) for l in prot.values()])
    d2 = ((xyz_lig[:, None, :] - xyz_prot[None, :, :]) ** 2).sum(-1)
    dmin = float(np.sqrt(d2.min()))
    print(f"\n  menor distância ligante–proteína: {dmin:.2f} Å")
    # Os dois limiares são física, não gosto. Um contato não-ligado entre
    # átomos pesados mede 2,8 Å numa estrutura refinada; 2,0 Å é um choque que
    # a minimização resolve; abaixo de SOBREPOSICAO_A não existe química — é
    # um arquivo com dois modelos somados, e minimizar isso não conserta.
    # O limiar antigo era 0,5 Å, que só pegava coincidência literal: um modelo
    # com átomos a 0,6 Å passava, e estourava o GROMACS dez minutos depois.
    if dmin < SOBREPOSICAO_A:
        raise SystemExit(
            f"  átomos sobrepostos ({dmin:.2f} Å < {SOBREPOSICAO_A} Å). Isto"
            f" não é choque,\n  é o mesmo espaço ocupado duas vezes: o arquivo"
            f" tem mais de um modelo\n  somado. Conserta-se na origem, não na"
            f" minimização.")
    if dmin < CHOQUE_A:
        print(f"      [ATENÇÃO] choque de {dmin:.2f} Å (contato normal: 2,8 Å)."
              f"\n      A minimização deve resolver, mas confira o Fmax do em2:"
              f"\n      se ele ficar alto, a geometria de partida é a causa.")

    # --- conferir contra o SMILES ----------------------------------------
    m_ref = Chem.MolFromSmiles(smiles)
    if m_ref is None:
        raise SystemExit(f"  SMILES inválido em {args.protac_smi}")
    esperado = m_ref.GetNumAtoms()
    carga = Chem.GetFormalCharge(m_ref)
    pesados = len([l for l in lig_linhas if not eh_hidrogenio(l)])
    print(f"  PROTAC do .smi: {esperado} átomos pesados, carga formal {carga:+d}")
    print(f"  ligante no modelo: {pesados} átomos pesados (de {n_lig} com H)")
    if pesados != esperado:
        msg = (f"\n  [DIVERGÊNCIA] o modelo tem {pesados} átomos pesados e o "
               f"SMILES pede {esperado}.")
        if not args.aceitar_divergencia:
            raise SystemExit(
                msg + "\n  Parametrizar um ligante que não é o PROTAC produz "
                "uma trajetória\n  bonita de uma molécula errada. Confira o "
                "modelo, ou passe\n  --aceitar-divergencia se souber por que "
                "diferem.")
        print(msg + "  (aceito por --aceitar-divergencia)")

    # ordens de ligação do SMILES sobre as coordenadas do modelo
    print("\n  casando o SMILES com as coordenadas do modelo")
    mol = ligante_com_ordens(lig_linhas, smiles)
    print(f"      casou: {mol.GetNumAtoms()} átomos com H, "
          f"{mol.GetNumBonds()} ligações")

    # --- escrever ---------------------------------------------------------
    # resolve(): o md_run.sh faz `cd "$MD_DIR"` antes de usar estes caminhos,
    # e um caminho relativo no JSON apontaria para dentro do próprio md_dir.
    out = args.outdir.expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)

    rec = out / "receptor_ternario.pdb"
    with open(rec, "w") as fh:
        for c in sorted(prot):
            fh.writelines(l + "\n" for l in prot[c])
            fh.write("TER\n")
        fh.write("END\n")

    # O ligante para o acpype sai do RDKit: um resname só, CONECT presente,
    # hidrogênios explícitos. É o mesmo arquivo que o nível (ii) entregava, e
    # batizado pela MESMA função.
    nomear_residuo(mol, LIG_RESNAME)
    lig = out / "protac_ternario.pdb"
    Chem.MolToPDBFile(mol, str(lig))

    # protac.sdf: o fix_receptor_for_md.py e o split_chain_gaps.py o leem para
    # medir distância ao ligante. O nome é o que eles esperam.
    sdf = out / "protac.sdf"
    with Chem.SDWriter(str(sdf)) as w:
        w.write(mol)

    complexo = out / "complexo_ternario.pdb"
    with open(complexo, "w") as fh:
        fh.write(rec.read_text().replace("END\n", ""))
        fh.write("\n".join(l for l in lig.read_text().splitlines()
                           if l.startswith(("HETATM", "CONECT"))) + "\nEND\n")

    # --- o custo, antes de tomar a máquina --------------------------------
    xyz_tudo = np.vstack([xyz_prot, xyz_lig])
    custo = custo_da_caixa(xyz_tudo)
    print(f"\n  tamanho do sistema solvatado (caixa cúbica, 1,3 nm de folga):")
    print(f"      maior extensão do complexo: {extensao(xyz_tudo)/10:.1f} nm")
    print(f"      aresta da caixa: {custo['aresta_nm']} nm")
    print(f"      ~{custo['n_aguas_estimado']} águas, "
          f"~{custo['n_atomos_total_estimado']} átomos no total")
    print(f"      RAM estimada: ~{custo['ram_estimada_gb']} GB  "
          f"(ordem de grandeza: 0,3 GB + 4 kB/átomo)")
    print(f"      A RAM é o menor dos custos. O que a MD realmente OCUPA numa"
          f"\n      máquina compartilhada é a GPU, por dias seguidos: o"
          f" md_run.sh\n      usa `-nb gpu -pme gpu -bonded gpu -update gpu`"
          f" e a placa fica\n      dedicada. Combine o horário com quem mais"
          f" usa a workstation.")
    if custo["ram_estimada_gb"] > args.ram_limite_gb:
        raise SystemExit(
            f"\n*** O sistema pede ~{custo['ram_estimada_gb']} GB e o limite"
            f" é {args.ram_limite_gb} GB.\n"
            f"*** A workstation é compartilhada: tomar a memória dela inteira"
            f" derruba\n*** o trabalho de outras pessoas, e esse custo não é"
            f" nosso para gastar.\n*** \n"
            f"*** Uma aresta de {custo['aresta_nm']} nm para um complexo de"
            f" {extensao(xyz_tudo)/10:.1f} nm costuma\n"
            f"*** significar que o modelo tem cadeia estendida ou mais de um\n"
            f"*** complexo no arquivo — confira o receptor antes de aumentar o\n"
            f"*** limite. Um ternário E3+PCSK9 real mede ~10 nm, dá ~120 mil\n"
            f"*** átomos e pede ~0,8 GB — bem abaixo do limite.\n"
            f"*** \n"
            f"*** Para seguir de propósito: --ram-limite-gb"
            f" {custo['ram_estimada_gb'] + 1:.0f}")

    info = {
        # --- o contrato com o md_run.sh ---
        "candidate_id": cand_id,
        "carga_formal": carga,
        "resname": LIG_RESNAME,
        "lig_pdb": str(lig),
        "receptor_pdb": str(rec),
        # --- o resto é registro para a tese ---
        "nivel": "iii",
        "nivel_descricao": "ternário completo — E3 + recrutador-linker-"
                           "warhead + PCSK9",
        "modelo_origem": str(modelo),
        "protac_smiles": smiles,
        "cadeias_proteina": {c: residuos(l) for c, l in sorted(prot.items())},
        "hetatm_originais": {k: len(v) for k, v in het.items()},
        "n_atomos_ligante": mol.GetNumAtoms(),
        "n_atomos_pesados_esperado": esperado,
        "n_atomos_pesados_no_modelo": pesados,
        "menor_distancia_ligante_proteina_A": round(dmin, 2),
        "lig_sdf": str(sdf),
        "complexo_pdb": str(complexo),
        "custo_estimado": custo,
    }
    faltando = [k for k in CHAVES_DO_MD_RUN if info.get(k) in (None, "")]
    if faltando:
        raise SystemExit(f"  [BUG] faltam chaves que o md_run.sh lê: {faltando}")
    (out / "md_sistema.json").write_text(json.dumps(info, indent=2,
                                                    ensure_ascii=False))

    print(f"\n  escrito em {out}:")
    for p in (rec, lig, sdf, complexo, out / "md_sistema.json"):
        print(f"    {p.name}")
    print(f"\n  contrato com o md_run.sh conferido: "
          f"{', '.join(CHAVES_DO_MD_RUN)}")
    print(f"  carga formal do PROTAC: {carga:+d}  (vai em `acpype -n`)")


if __name__ == "__main__":
    main()
