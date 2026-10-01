#!/usr/bin/env python
"""
Fase 10c — analisa a MD do complexo TERNÁRIO e decide se ele é estável.

Roda no env `mdtools` (MDAnalysis). Minutos.

    python md_analyze_ternario.py --md-dir ~/pipeline_vhl/md_ternario

A pergunta é outra, e por isso o script é outro
-----------------------------------------------
O nível (ii) perguntava se a warhead gruda na PRÓPRIA E3 na ausência da PCSK9 —
uma interação espúria, que reprova o candidato. O nível (iii) pergunta o que o
trabalho existe para responder:

    a interface ternária que o PRosettaC modelou sobrevive ao solvente,
    ou ela só existia porque o docking a construiu?

E, respondida essa, a pergunta biológica que vem depois:

    sobra alguma lisina da PCSK9 exposta e virada para a E3?

Sem lisina acessível não há transferência de ubiquitina, e um complexo
ternário perfeitamente estável sem lisina exposta não degrada nada. Esta é a
medida que distingue "o modelo é estável" de "o modelo é produtivo".

O que ele mede (as análises que o WP3 pede)
-------------------------------------------
    RMSD            do conjunto, de cada proteína, e o da INTERFACE
    RMSF            por resíduo, em cada cadeia
    contatos        E3 <-> PCSK9, e PROTAC <-> cada uma das duas
    ligações de H   atravessando a interface
    torções         do linker: ele explora conformações ou travou?
    SASA das Lys    da PCSK9, e a distância de cada uma à E3
    agrupamento     as conformações visitadas são uma família ou várias?

O RMSD da interface é o número central, e ele é diferente do RMSD global.
Superpõe-se a E3 e mede-se o deslocamento da PCSK9: é isso que diz se o alvo
ESCORREGOU em relação à ligase. Um RMSD global de 4 Å com interface de 1,5 Å é
movimento de domínio distal; interface de 8 Å com global de 4 Å é o complexo
se desfazendo. Os dois juntos informam, cada um sozinho engana.

Honestidade sobre a lisina
--------------------------
A zona de transferência de ubiquitina fica no E2, carregado pelo complexo
CRL (Cul2-Rbx1 para a VHL), e esse complexo NÃO está no nosso modelo. Então
aqui não se prediz transferência: mede-se acessibilidade da lisina e a
geometria dela em relação à E3 que temos. É triagem, não predição, e a tese
precisa dizer isso com essas palavras.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
# Carregar a trajetória, achar a topologia legível, descobrir o nome do resíduo
# do ligante e achar os términos artificiais são os MESMOS problemas do nível
# (ii), já resolvidos lá — inclusive o .tpr versão 138 que o MDAnalysis não lê.
# Reimplementar aqui seria a terceira vez que este projeto duplica uma busca e
# vê as duas cópias divergirem.
from md_analyze import carregar, descobrir_resname  # noqa: E402

# --- os critérios, explícitos, para citar na tese -------------------------
CRITERIOS = {
    # A interface é o portão. 3,0 Å de deslocamento da PCSK9 em relação à E3
    # é o limite do que ainda é a mesma pose; acima disso o modelo mudou de
    # hipótese no meio da simulação.
    "rmsd_interface_max_A": 3.0,
    # Contatos mantidos em relação ao frame inicial, que é a pose do PRosettaC.
    "retencao_interface_min": 0.5,
    "fracao_frames_interface_min": 0.7,
    # O PROTAC tem de permanecer ancorado nas DUAS pontas: perder uma é perder
    # o motivo de ele existir.
    "retencao_protac_por_cadeia_min": 0.4,
    # Lisina produtiva: exposta ao solvente e do lado da E3.
    "sasa_lisina_min_A2": 30.0,
    "dist_lisina_e3_max_A": 20.0,
    "n_lisinas_produtivas_min": 1,
}

# Raios de Bondi, em Å. O probe de 1,4 Å é o raio da água.
RAIOS = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80,
         "P": 1.80, "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98}
PROBE = 1.4
RAIO_PADRAO = 1.70

CORTE_CONTATO_A = 4.5
# Interface = CA a esta distância da outra proteína, medida no frame 0.
CORTE_INTERFACE_A = 10.0
# Ligação de H: H...aceptor <= 2,5 Å e ângulo D-H...A >= 120 graus. É o
# critério geométrico padrão, e está escrito aqui para poder ser citado.
HB_DIST_A = 2.5
HB_ANGULO_GRAUS = 120.0


def elemento(nome: str) -> str:
    """Elemento a partir do nome do átomo: um .gro não guarda elementos."""
    n = nome.strip().upper().lstrip("0123456789")
    for e in ("CL", "BR"):
        if n.startswith(e):
            return e
    return n[:1] if n else "C"


def raio(nome: str) -> float:
    return RAIOS.get(elemento(nome), RAIO_PADRAO)


# ===========================================================================
# quem é E3 e quem é PCSK9
# ===========================================================================
def rotular_cadeias(u, receptor_pdb: Path):
    """Cadeia de cada resíduo de proteína da trajetória.

    O .gro não guarda cadeia, e a trajetória veio de um .gro. Mas o PDB que o
    pdb2gmx consumiu guarda — e o pdb2gmx preserva a ORDEM dos resíduos. Então
    a cadeia sai por posição, do arquivo que originou a topologia.

    Isto é verificado, não suposto: as contagens têm de bater e os nomes de
    resíduo têm de concordar. Não batendo, o script para e imprime os dois
    números, porque atribuir a cadeia errada inverteria E3 e alvo em TODAS as
    medidas abaixo sem dar nenhum sinal.
    """
    linhas = [l for l in receptor_pdb.read_text(errors="ignore").splitlines()
              if l.startswith("ATOM")]
    ordem, vistos = [], set()
    for l in linhas:
        chave = (l[21], l[22:27])
        if chave not in vistos:
            vistos.add(chave)
            ordem.append((l[21], l[17:20].strip()))

    residuos = u.select_atoms("protein").residues
    if len(ordem) != len(residuos):
        raise SystemExit(
            f"\n*** O receptor tem {len(ordem)} resíduos e a trajetória tem "
            f"{len(residuos)}.\n"
            f"*** Sem essa correspondência não se sabe quais resíduos são da "
            f"E3 e\n*** quais são da PCSK9 — e trocar as duas inverteria todas "
            f"as medidas\n*** sem dar sinal nenhum.\n"
            f"*** O receptor usado foi {receptor_pdb.name}. Se o md_run.sh "
            f"cortou\n*** lacunas (receptor_capped.pdb) ou completou resíduos\n"
            f"*** (receptor_fixed.pdb), passe ESSE arquivo em --receptor.")

    # HIS vira HISE/HISD/HSE no campo de força; GLU/ASP protonados idem
    def compativel(a, b):
        return a == b or a[:2] == b[:2]

    divergem = sum(1 for (_, rn), r in zip(ordem, residuos)
                   if not compativel(rn, r.resname))
    if divergem > 0.05 * len(ordem):
        exemplos = [f"{rn}/{r.resname}" for (_, rn), r in zip(ordem, residuos)
                    if not compativel(rn, r.resname)][:6]
        raise SystemExit(
            f"\n*** {divergem} de {len(ordem)} resíduos não concordam entre o "
            f"receptor e a\n*** trajetória (ex.: {', '.join(exemplos)}). A "
            f"ordem não é a mesma, então\n*** a atribuição de cadeia seria "
            f"inventada.")

    rotulo = {}
    for (cad, _), r in zip(ordem, residuos):
        rotulo[int(r.resid)] = cad
    return rotulo, divergem


def selecao_de_cadeia(u, rotulo, cadeia):
    resids = sorted(r for r, c in rotulo.items() if c == cadeia)
    if not resids:
        return u.select_atoms("protein and name CA and resid 999999")
    return u.select_atoms("protein and resid " + " ".join(map(str, resids)))


# ===========================================================================
# SASA — Shrake-Rupley, em numpy
# ===========================================================================
def pontos_esfera(n=256):
    """Pontos quase uniformes na esfera unitária (espiral de Fibonacci)."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(phi)], axis=1)


_ESFERA = pontos_esfera()


def sasa(xyz, raios, alvos, vizinhanca=12.0):
    """SASA por átomo (Å²) para os índices em `alvos`.

    Shrake-Rupley: cada átomo ganha uma esfera de pontos no raio r+probe, e
    conta-se a fração de pontos que nenhum vizinho cobre. Implementado aqui em
    vez de chamado do `gmx sasa` porque o gmx exige grupos do make_ndx, e a
    numeração desses grupos já custou a este projeto uma etapa inteira.

    Conferido contra os dois casos que têm solução fechada, porque um número
    que vai para a tese não se aceita por parecer razoável:

        átomo isolado      120,76 Å² contra 4*pi*(r+probe)^2 = 120,76  (0,00%)
        par ligado a 1,54  151,43 Å² contra as duas calotas = 150,76   (0,44%)

    O erro de 0,44% é a discretização dos 256 pontos da esfera.
    """
    raios = np.asarray(raios, dtype=float) + PROBE
    out = np.zeros(len(alvos))
    for k, i in enumerate(alvos):
        d = np.linalg.norm(xyz - xyz[i], axis=1)
        viz = np.where((d > 0) & (d < raios[i] + raios.max()) &
                       (d < vizinhanca + raios[i]))[0]
        pts = xyz[i] + _ESFERA * raios[i]
        if len(viz):
            dd = np.linalg.norm(pts[:, None, :] - xyz[viz][None, :, :], axis=2)
            livre = (dd >= raios[viz][None, :]).all(axis=1)
        else:
            livre = np.ones(len(pts), dtype=bool)
        out[k] = 4 * np.pi * raios[i] ** 2 * livre.mean()
    return out


# ===========================================================================
# ligações de hidrogênio, pelo critério geométrico
# ===========================================================================
def doadores_e_aceptores(u, grupo):
    """(pares doador-H, aceptores) do grupo, pelas distâncias do frame 0.

    Um .gro não traz ligações, então 'H ligado a N/O' sai da geometria: H a
    menos de 1,2 Å de um N ou O está ligado a ele. Não é heurística frouxa —
    é a distância de uma ligação N-H (1,01 Å) ou O-H (0,96 Å), e o segundo
    vizinho mais próximo está a mais de 1,8 Å.
    """
    pesados = grupo.select_atoms("not name H*")
    nos = np.array([i for i, a in enumerate(pesados)
                    if elemento(a.name) in ("N", "O")], dtype=int)
    hs = grupo.select_atoms("name H*")
    pares = []
    if len(hs) and len(nos):
        xyz_no = pesados.positions[nos]
        d = np.linalg.norm(hs.positions[:, None, :] - xyz_no[None, :, :], axis=2)
        for ih, linha in enumerate(d):
            j = int(np.argmin(linha))
            if linha[j] <= 1.2:
                pares.append((int(pesados[nos[j]].index), int(hs[ih].index)))
    aceptores = [int(pesados[i].index) for i in nos]
    return pares, aceptores


def contar_pares(xyz_a, xyz_b, corte=CORTE_CONTATO_A):
    """Pares de átomos a menos de `corte`, por lista de células.

    A matriz cheia custa n_a x n_b por frame. Num ternário real isso é
    2500 x 7000 = 17 milhões de distâncias, mil vezes — minutos virando horas
    para contar pares que são milhares. O `capped_distance` do MDAnalysis usa
    células e devolve SÓ os pares dentro do corte; o resultado é idêntico.
    """
    from MDAnalysis.lib.distances import capped_distance
    if not len(xyz_a) or not len(xyz_b):
        return 0
    pares = capped_distance(xyz_a, xyz_b, max_cutoff=corte,
                            return_distances=False)
    return int(len(pares))


def contar_hbonds(u, pares_a, acept_a, pares_b, acept_b):
    """Ligações de H atravessando a interface, nos dois sentidos."""
    from MDAnalysis.lib.distances import capped_distance
    pos = u.atoms.positions
    total = 0
    for pares, acept in ((pares_a, acept_b), (pares_b, acept_a)):
        if not pares or not acept:
            continue
        xyz_d = pos[[d for d, _ in pares]]
        xyz_h = pos[[h for _, h in pares]]
        xyz_a = pos[acept]
        vizinhos = capped_distance(xyz_h, xyz_a, max_cutoff=HB_DIST_A,
                                   return_distances=False)
        if not len(vizinhos):
            continue
        ii, jj = vizinhos[:, 0], vizinhos[:, 1]
        v1 = xyz_d[ii] - xyz_h[ii]
        v2 = xyz_a[jj] - xyz_h[ii]
        cos = (v1 * v2).sum(1) / (np.linalg.norm(v1, axis=1)
                                  * np.linalg.norm(v2, axis=1) + 1e-9)
        ang = np.degrees(np.arccos(np.clip(cos, -1, 1)))
        total += int((ang >= HB_ANGULO_GRAUS).sum())
    return total


# ===========================================================================
# torções do linker
# ===========================================================================
def ligacoes_do_ligante(lig, sdf: Path | None):
    """Pares ligados do ligante — do SDF quando ele CONFERE, senão da geometria.

    O SDF escrito pelo md_prepare_ternario tem ordens de ligação e anéis
    exatos, vindos do SMILES. Mas a ordem dos átomos na trajetória passou pelo
    antechamber, e supor que ela se manteve seria exatamente o tipo de suposição
    que custou caro a este projeto. Então compara-se a sequência de elementos:
    batendo, usam-se as ligações do SDF; não batendo, caem-se para a geometria
    e isso é dito em voz alta.
    """
    pesados = lig.select_atoms("not name H*")
    els = [elemento(a.name) for a in pesados]
    if sdf and sdf.exists():
        try:
            from rdkit import Chem
            m = Chem.RemoveHs(Chem.SDMolSupplier(str(sdf), removeHs=True)[0])
            els_sdf = [a.GetSymbol().upper() for a in m.GetAtoms()]
            if els_sdf == els:
                print("      ligações do ligante: do protac.sdf (ordens e "
                      "anéis exatos)")
                return ([(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
                         for b in m.GetBonds()],
                        {(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
                         for b in m.GetBonds() if b.IsInRing()},
                        pesados)
            print(f"      [ATENÇÃO] a ordem dos átomos do ligante na "
                  f"trajetória não é a\n      do protac.sdf "
                  f"({len(els)} vs {len(els_sdf)} pesados) — o antechamber "
                  f"reordenou.\n      As ligações saem da geometria, e anel "
                  f"vs. não-anel é detectado por busca.")
        except Exception as e:
            print(f"      [ATENÇÃO] não li o protac.sdf ({e}); ligações pela "
                  f"geometria")

    xyz = pesados.positions
    d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=2)
    pares = [(int(i), int(j)) for i in range(len(xyz))
             for j in range(i + 1, len(xyz)) if d[i, j] <= 1.95]
    # anel: tirando a ligação, os dois átomos continuam conectados
    adj = {i: set() for i in range(len(xyz))}
    for i, j in pares:
        adj[i].add(j); adj[j].add(i)
    aneis = set()
    for i, j in pares:
        vistos, pilha = {i}, [i]
        while pilha:
            a = pilha.pop()
            for b in adj[a]:
                if (a, b) in ((i, j), (j, i)):
                    continue
                if b not in vistos:
                    vistos.add(b); pilha.append(b)
        if j in vistos:
            aneis.add((i, j))
    return pares, aneis, pesados


def torcoes_rotaveis(lig, sdf: Path | None):
    """Quádruplas de índices globais para cada ligação simples rotável."""
    pares, aneis, pesados = ligacoes_do_ligante(lig, sdf)
    adj = {i: set() for i in range(len(pesados))}
    for i, j in pares:
        adj[i].add(j); adj[j].add(i)
    quadras = []
    for i, j in pares:
        if (i, j) in aneis or (j, i) in aneis:
            continue
        vi = [a for a in adj[i] if a != j]
        vj = [b for b in adj[j] if b != i]
        if not vi or not vj:          # ponta de cadeia: não há torção
            continue
        quadras.append((int(pesados[vi[0]].index), int(pesados[i].index),
                        int(pesados[j].index), int(pesados[vj[0]].index)))
    return quadras


def diedro(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / (np.linalg.norm(b1) + 1e-9)
    v = b0 - (b0 @ b1) * b1
    w = b2 - (b2 @ b1) * b1
    return np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w))


def desvio_circular(angulos):
    """Desvio padrão circular, em graus: ângulo é periódico e a média comum
    diz que -179 e +179 estão a 358 graus de distância."""
    a = np.radians(np.asarray(angulos))
    R = np.hypot(np.mean(np.cos(a)), np.mean(np.sin(a)))
    R = min(max(R, 1e-9), 1.0)
    return float(np.degrees(np.sqrt(-2 * np.log(R))))


# ===========================================================================
# agrupamento das conformações (gromos)
# ===========================================================================
def agrupar(coords, corte=2.0):
    """Gromos: o frame com mais vizinhos vira cluster, sai, repete.

    `coords` é (n_frames, n_atomos, 3) já superposto. Devolve as populações.
    """
    n = len(coords)
    pl = coords.reshape(n, -1)
    d = np.sqrt(((pl[:, None, :] - pl[None, :, :]) ** 2).sum(-1) / coords.shape[1])
    restantes = set(range(n))
    pops = []
    while restantes:
        idx = sorted(restantes)
        sub = d[np.ix_(idx, idx)]
        viz = [(int(((sub[k] <= corte)).sum()), idx[k]) for k in range(len(idx))]
        viz.sort(reverse=True)
        tam, centro = viz[0]
        # O centro entra no proprio cluster SEMPRE. Sem isto, uma matriz com
        # NaN (que acontece quando a selecao da interface sai vazia) deixa
        # `membros` vazio, `restantes` nunca encolhe, e o laco nao termina —
        # foi exatamente assim que o primeiro teste travou, sem mensagem.
        membros = {j for j in idx if d[centro, j] <= corte} | {centro}
        pops.append(len(membros))
        restantes -= membros
    return pops


def superpor(movel, ref):
    """Kabsch: devolve `movel` rodado/transladado sobre `ref`."""
    cm, cr = movel.mean(0), ref.mean(0)
    A, B = movel - cm, ref - cr
    V, _, Wt = np.linalg.svd(A.T @ B)
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1, 1, d]) @ Wt
    return A @ R + cr, R, cm, cr


# ===========================================================================
# a análise de uma réplica
# ===========================================================================
def analisar_replica(u, rotulo, cad_e3, cad_alvo, resname, sdf, passo_sasa):
    from MDAnalysis.analysis.distances import distance_array

    e3 = selecao_de_cadeia(u, rotulo, cad_e3)
    alvo = selecao_de_cadeia(u, rotulo, cad_alvo)
    lig = u.select_atoms(f"resname {resname}")
    e3_ca = e3.select_atoms("name CA")
    alvo_ca = alvo.select_atoms("name CA")
    if not len(e3_ca) or not len(alvo_ca) or not len(lig):
        raise SystemExit("seleção vazia: E3, alvo ou ligante não encontrados")

    # --- referências no frame 0, que é a pose do PRosettaC ----------------
    u.trajectory[0]
    ref_e3 = e3_ca.positions.copy()
    ref_alvo = alvo_ca.positions.copy()
    ref_lig = lig.positions.copy()
    # resíduos de interface: CA a menos de CORTE_INTERFACE_A da outra proteína.
    # A matriz cheia aqui é aceitável porque é UMA vez, só com CA.
    d0 = distance_array(e3_ca.positions, alvo_ca.positions)
    ie3 = np.where(d0.min(axis=1) <= CORTE_INTERFACE_A)[0]
    ialvo = np.where(d0.min(axis=0) <= CORTE_INTERFACE_A)[0]
    contatos0 = contar_pares(e3.positions, alvo.positions)

    # SEM interface não há nível (iii). Isto tem de parar aqui, e com o número
    # na tela: sem a guarda, `ie3` vazio faz o RMSD da interface sair NaN, o
    # agrupamento receber uma matriz de NaN, e o veredito comparar NaN com o
    # critério — que dá False e REPROVA, dizendo "a interface não sobreviveu"
    # de um modelo cujas proteínas nunca se tocaram. Diagnóstico errado com
    # cara de resultado é pior que erro.
    if len(ie3) < 3 or len(ialvo) < 3:
        raise SystemExit(
            f"\n*** As duas proteínas NÃO se tocam no primeiro frame.\n"
            f"*** CA da E3 a menos de {CORTE_INTERFACE_A} Å do alvo: "
            f"{len(ie3)}; do alvo à E3: {len(ialvo)}.\n"
            f"*** Menor distância CA-CA: {float(d0.min()):.1f} Å.\n"
            f"*** \n"
            f"*** Um complexo ternário tem interface. Sem ela, ou o modelo do\n"
            f"*** PRosettaC não é ternário, ou a PBC não foi corrigida e as\n"
            f"*** duas proteínas foram envolvidas em lados opostos da caixa —\n"
            f"*** rode antes: bash scripts/md_fix_pbc.sh <dir_md>")

    pares_e3, acept_e3 = doadores_e_aceptores(u, e3)
    pares_al, acept_al = doadores_e_aceptores(u, alvo)

    quadras = torcoes_rotaveis(lig, sdf)
    print(f"      {len(quadras)} torções rotáveis no PROTAC")

    # lisinas do ALVO: é a PCSK9 que precisa ser ubiquitinada
    lys_nz = alvo.select_atoms("resname LYS LYN and name NZ")
    lys_sc = {int(a.resid): alvo.select_atoms(
                  f"resid {int(a.resid)} and name CD CE NZ") for a in lys_nz}
    print(f"      {len(lys_nz)} lisinas na cadeia do alvo ({cad_alvo})")

    serie = {k: [] for k in ("rmsd_e3", "rmsd_alvo", "rmsd_interface",
                             "rmsd_protac", "contatos", "hbonds",
                             "contatos_lig_e3", "contatos_lig_alvo")}
    # Guardar os CA aqui custa 610 x n_frames x 3 floats e poupa DUAS
    # passagens extras pela trajetória — num .xtc de 200 ns cada passagem é
    # leitura de disco, não aritmética.
    traj_ca = {cad_e3: [], cad_alvo: []}
    torsoes = [[] for _ in quadras]
    sasa_lys = {r: [] for r in lys_sc}
    dist_lys = {r: [] for r in lys_sc}
    conf_interface = []

    todos = u.atoms
    raios_todos = np.array([raio(a.name) for a in todos])
    idx_sasa_alvo = None

    for nf, _ in enumerate(u.trajectory):
        pos_e3 = e3_ca.positions
        pos_alvo = alvo_ca.positions
        traj_ca[cad_e3].append(pos_e3.copy())
        traj_ca[cad_alvo].append(pos_alvo.copy())

        # RMSD de cada proteína sobre si mesma: movimento interno
        sup, _, _, _ = superpor(pos_e3, ref_e3)
        serie["rmsd_e3"].append(float(np.sqrt(((sup - ref_e3) ** 2).sum(1).mean())))
        sup, _, _, _ = superpor(pos_alvo, ref_alvo)
        serie["rmsd_alvo"].append(
            float(np.sqrt(((sup - ref_alvo) ** 2).sum(1).mean())))

        # O NÚMERO CENTRAL: superpõe na E3 e mede o quanto o ALVO andou.
        # É este que diz se a PCSK9 escorregou em relação à ligase.
        _, R, cm, cr = superpor(pos_e3[ie3], ref_e3[ie3])
        alvo_no_quadro = (pos_alvo[ialvo] - cm) @ R + cr
        serie["rmsd_interface"].append(float(np.sqrt(
            ((alvo_no_quadro - ref_alvo[ialvo]) ** 2).sum(1).mean())))
        conf_interface.append(alvo_no_quadro.copy())

        sup, _, _, _ = superpor(lig.positions, ref_lig)
        serie["rmsd_protac"].append(
            float(np.sqrt(((sup - ref_lig) ** 2).sum(1).mean())))

        serie["contatos"].append(contar_pares(e3.positions, alvo.positions))
        serie["hbonds"].append(contar_hbonds(u, pares_e3, acept_e3,
                                             pares_al, acept_al))
        serie["contatos_lig_e3"].append(
            contar_pares(lig.positions, e3.positions))
        serie["contatos_lig_alvo"].append(
            contar_pares(lig.positions, alvo.positions))

        for k, (a, b, c, dd) in enumerate(quadras):
            p = todos.positions
            torsoes[k].append(diedro(p[a], p[b], p[c], p[dd]))

        # SASA é o passo caro: um frame a cada `passo_sasa`
        if nf % passo_sasa == 0 and lys_sc:
            xyz = todos.positions
            if idx_sasa_alvo is None:
                idx_sasa_alvo = {r: [int(a.index) for a in g]
                                 for r, g in lys_sc.items()}
            centro_e3 = e3.positions.mean(0)
            for r, idxs in idx_sasa_alvo.items():
                sasa_lys[r].append(float(sasa(xyz, raios_todos, idxs).sum()))
                nz = [i for i in idxs
                      if elemento(todos[i].name) == "N"]
                ponto = xyz[nz[0]] if nz else xyz[idxs[-1]]
                dist_lys[r].append(float(np.linalg.norm(ponto - centro_e3)))

    # --- RMSF por resíduo, em cada cadeia ---------------------------------
    # Superpõe-se cada frame na MÉDIA, não no frame 0: o RMSF mede flutuação
    # em torno da posição média, e usar o frame 0 como referência mistura
    # flutuação com a deriva do conjunto.
    rmsf = {}
    for nome, sel in ((cad_e3, e3_ca), (cad_alvo, alvo_ca)):
        traj = np.array(traj_ca[nome])
        media = traj.mean(0)
        sup = np.array([superpor(f, media)[0] for f in traj])
        rmsf[nome] = {int(a.resid): float(v) for a, v in zip(
            sel, np.sqrt(((sup - sup.mean(0)) ** 2).sum(2).mean(0)))}

    n = len(serie["contatos"])
    c = np.array(serie["contatos"], dtype=float)
    c0 = max(contatos0, 1)
    cl_e3 = np.array(serie["contatos_lig_e3"], dtype=float)
    cl_al = np.array(serie["contatos_lig_alvo"], dtype=float)

    # lisinas produtivas: expostas E do lado da E3
    lys_resumo = []
    for r in sorted(sasa_lys):
        if not sasa_lys[r]:
            continue
        s, dd = float(np.mean(sasa_lys[r])), float(np.mean(dist_lys[r]))
        lys_resumo.append({
            "resid": r, "sasa_media_A2": round(s, 1),
            "dist_media_a_E3_A": round(dd, 1),
            "exposta": s >= CRITERIOS["sasa_lisina_min_A2"],
            "do_lado_da_E3": dd <= CRITERIOS["dist_lisina_e3_max_A"],
            "produtiva": (s >= CRITERIOS["sasa_lisina_min_A2"]
                          and dd <= CRITERIOS["dist_lisina_e3_max_A"]),
        })
    n_prod = sum(1 for l in lys_resumo if l["produtiva"])

    tor_desvios = [desvio_circular(t) for t in torsoes] if quadras else []
    pops = agrupar(np.array(conf_interface)) if len(conf_interface) > 1 else [n]

    return {
        "n_frames": n,
        "rmsd_e3_media": float(np.mean(serie["rmsd_e3"])),
        "rmsd_alvo_media": float(np.mean(serie["rmsd_alvo"])),
        "rmsd_interface_media": float(np.mean(serie["rmsd_interface"])),
        "rmsd_interface_final": float(serie["rmsd_interface"][-1]),
        "rmsd_protac_media": float(np.mean(serie["rmsd_protac"])),
        "n_residuos_interface_e3": int(len(ie3)),
        "n_residuos_interface_alvo": int(len(ialvo)),
        "contatos_inicial": contatos0,
        "contatos_media": float(c.mean()),
        "retencao_interface": float(c.mean() / c0),
        "fracao_frames_interface": float(
            (c >= CRITERIOS["retencao_interface_min"] * c0).mean()),
        "hbonds_media": float(np.mean(serie["hbonds"])),
        "hbonds_inicial": int(serie["hbonds"][0]),
        "retencao_protac_e3": float(cl_e3.mean() / max(cl_e3[0], 1)),
        "retencao_protac_alvo": float(cl_al.mean() / max(cl_al[0], 1)),
        "n_torcoes": len(quadras),
        "torcoes_travadas": int(sum(1 for d in tor_desvios if d < 30.0)),
        "torcao_desvio_mediano_graus": float(np.median(tor_desvios))
                                        if tor_desvios else 0.0,
        "lisinas": lys_resumo,
        "n_lisinas_produtivas": n_prod,
        "clusters_interface": pops[:10],
        "fracao_maior_cluster": float(max(pops) / n) if n else 0.0,
        "rmsf_por_residuo": rmsf,
        "_series": {k: list(map(float, v)) for k, v in serie.items()},
    }


# ===========================================================================
def veredito(m):
    """Regras explícitas, uma linha por critério, para citar na tese."""
    testes = [
        ("RMSD da interface", m["rmsd_interface_media"],
         f"<= {CRITERIOS['rmsd_interface_max_A']} Å",
         m["rmsd_interface_media"] <= CRITERIOS["rmsd_interface_max_A"]),
        ("retenção da interface", m["retencao_interface"],
         f">= {CRITERIOS['retencao_interface_min']}",
         m["retencao_interface"] >= CRITERIOS["retencao_interface_min"]),
        ("frames com interface", m["fracao_frames_interface"],
         f">= {CRITERIOS['fracao_frames_interface_min']}",
         m["fracao_frames_interface"] >= CRITERIOS["fracao_frames_interface_min"]),
        ("PROTAC ancorado na E3", m["retencao_protac_e3"],
         f">= {CRITERIOS['retencao_protac_por_cadeia_min']}",
         m["retencao_protac_e3"] >= CRITERIOS["retencao_protac_por_cadeia_min"]),
        ("PROTAC ancorado no alvo", m["retencao_protac_alvo"],
         f">= {CRITERIOS['retencao_protac_por_cadeia_min']}",
         m["retencao_protac_alvo"] >= CRITERIOS["retencao_protac_por_cadeia_min"]),
        ("lisinas produtivas", m["n_lisinas_produtivas"],
         f">= {CRITERIOS['n_lisinas_produtivas_min']}",
         m["n_lisinas_produtivas"] >= CRITERIOS["n_lisinas_produtivas_min"]),
    ]
    return testes, all(t[3] for t in testes)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--md-dir", type=Path, required=True)
    ap.add_argument("--receptor", type=Path,
                    help="PDB com as cadeias, que originou a topologia "
                         "(default: o receptor_capped/fixed/ternario do md-dir)")
    ap.add_argument("--resname", default=None)
    ap.add_argument("--e3-chain", help="cadeia da E3 (default: a menor)")
    ap.add_argument("--passo-sasa", type=int, default=10,
                    help="SASA a cada N frames (é o passo caro; default 10)")
    args = ap.parse_args()

    md = args.md_dir.expanduser().resolve()
    info_p = md / "md_sistema.json"
    if not info_p.exists():
        raise SystemExit(f"não achei {info_p} — rode md_prepare_ternario.py")
    info = json.loads(info_p.read_text())
    if info.get("nivel") != "iii":
        print(f"  [ATENÇÃO] este md_sistema.json diz nível "
              f"'{info.get('nivel', '?')}', não 'iii'. Esta análise pressupõe "
              f"DUAS proteínas;\n  para o nível (ii) o script é o md_analyze.py.")

    # O receptor certo é o que o pdb2gmx consumiu, e o md_run.sh o conserta em
    # dois passos. Procurar na ordem inversa da cirurgia é procurar o último.
    rec = args.receptor
    if rec is None:
        for nome in ("receptor_capped.pdb", "receptor_fixed.pdb",
                     "receptor_ternario.pdb"):
            if (md / nome).exists():
                rec = md / nome
                break
    if rec is None or not rec.exists():
        raise SystemExit(f"não achei o PDB do receptor em {md}")
    print(f"receptor com as cadeias: {rec.name}")

    reps = sorted(d.name for d in md.glob("rep*") if d.is_dir())
    if not reps:
        raise SystemExit(f"nenhuma réplica rep* em {md}")

    linhas, detalhe = [], {}
    for rep in reps:
        print(f"\n[{rep}]")
        u = carregar(md, rep)
        if u is None:
            print("      sem prod.xtc — pulando")
            continue
        resname = descobrir_resname(u, args.resname or info.get("resname", "PTC"))
        rotulo, divergem = rotular_cadeias(u, rec)
        cadeias = sorted(set(rotulo.values()))
        if len(cadeias) < 2:
            raise SystemExit(
                f"o receptor tem {len(cadeias)} cadeia(s) ({cadeias}). O nível "
                f"(iii) exige duas.")
        if len(cadeias) > 2:
            print(f"      [ATENÇÃO] {len(cadeias)} cadeias ({cadeias}); as "
                  f"duas maiores serão E3 e alvo")
        tamanhos = {c: sum(1 for v in rotulo.values() if v == c)
                    for c in cadeias}
        duas = sorted(tamanhos, key=lambda c: -tamanhos[c])[:2]
        # A E3 é a MENOR das duas: VHL tem ~160 resíduos, a PCSK9 ~450. Dizer
        # isso em voz alta é mais seguro que confiar no ID da cadeia, que vem
        # do PDB de origem e já mudou de nome duas vezes neste projeto.
        cad_e3 = args.e3_chain or min(duas, key=lambda c: tamanhos[c])
        cad_alvo = [c for c in duas if c != cad_e3][0]
        print(f"      E3 = cadeia {cad_e3} ({tamanhos[cad_e3]} resíduos) | "
              f"alvo = cadeia {cad_alvo} ({tamanhos[cad_alvo]} resíduos)")
        if tamanhos[cad_e3] > tamanhos[cad_alvo]:
            print(f"      [CONFIRA] a E3 ficou MAIOR que o alvo. Se isso está "
                  f"trocado,\n      passe --e3-chain {cad_alvo}.")

        m = analisar_replica(u, rotulo, cad_e3, cad_alvo, resname,
                             md / "protac.sdf", args.passo_sasa)
        m["replica"] = rep
        m["cadeia_e3"], m["cadeia_alvo"] = cad_e3, cad_alvo
        detalhe[rep] = m
        linhas.append({k: v for k, v in m.items()
                       if not k.startswith("_")
                       and k not in ("lisinas", "rmsf_por_residuo",
                                     "clusters_interface")})

    if not linhas:
        raise SystemExit("nenhuma réplica analisada")

    df = pd.DataFrame(linhas).set_index("replica")
    csv = md / "ternario_metricas.csv"
    df.to_csv(csv)
    (md / "ternario_detalhe.json").write_text(
        json.dumps(detalhe, indent=2, ensure_ascii=False))

    print("\n" + "=" * 70)
    print("MD NÍVEL (iii) — complexo ternário")
    print("=" * 70)
    colunas = ["rmsd_interface_media", "retencao_interface", "hbonds_media",
               "retencao_protac_e3", "retencao_protac_alvo",
               "n_lisinas_produtivas", "fracao_maior_cluster"]
    print(df[colunas].round(2).to_string())

    media = {k: float(df[k].mean()) for k in df.columns
             if pd.api.types.is_numeric_dtype(df[k])}
    testes, passou = veredito(media)
    print("\nveredito (média das réplicas):")
    for nome, valor, regra, ok in testes:
        print(f"  [{'OK ' if ok else 'NÃO'}] {nome:<26} {valor:>8.2f}  {regra}")

    print("\nlisinas da PCSK9, por réplica:")
    for rep, m in detalhe.items():
        prod = [l for l in m["lisinas"] if l["produtiva"]]
        print(f"  {rep}: {len(prod)} produtiva(s) de {len(m['lisinas'])}"
              + (f" — resíduos {[l['resid'] for l in prod]}" if prod else ""))
    print("\n  Produtiva = SASA >= "
          f"{CRITERIOS['sasa_lisina_min_A2']} Å² e a <= "
          f"{CRITERIOS['dist_lisina_e3_max_A']} Å do centro da E3.")
    print("  ISTO NÃO É PREDIÇÃO DE TRANSFERÊNCIA. A zona de transferência de")
    print("  ubiquitina fica no E2, carregado pelo complexo CRL (Cul2-Rbx1 para")
    print("  a VHL), que NÃO está neste modelo. É triagem de acessibilidade, e")
    print("  a tese precisa dizer isso com estas palavras.")

    print(f"\n{'APROVADO' if passou else 'REPROVADO'} nos critérios do nível (iii)")
    if not passou:
        print("  Reprovar aqui é resultado: a interface que o PRosettaC modelou")
        print("  não sobreviveu ao solvente, e isso é uma conclusão sobre o")
        print("  candidato — não uma falha do pipeline.")
    print(f"\n  {csv}")
    print(f"  {md / 'ternario_detalhe.json'}")
    return 0 if passou else 6


if __name__ == "__main__":
    raise SystemExit(main())
