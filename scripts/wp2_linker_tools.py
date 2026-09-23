#!/usr/bin/env python
"""
WP2 revisado — seleção de linkers, colocação no exit vector e contabilidade dos
pontos de conjugação, em acordo com o que o WP3 espera receber.

Três defeitos do WP2 original que este módulo corrige
------------------------------------------------------

1. `has_terminal_group()` aceitava linker MONOfuncional.
   Ela devolve True com UM match de `[NH2]`/`[OH]`/etc. Um linker de PROTAC é
   bifuncional por definição: uma ponta vai no recrutador, a outra no warhead.
   Um linker com uma única alça passa no filtro e só falha lá na frente, na
   montagem do WP3. Pior, `HasSubstructMatch('[NH2]')` casa amina NÃO terminal
   (`CC(N)CC` casa), então "grupo terminal" não era terminal.

2. `exit_vector_score()` media confôrmero solto no espaço.
   `EmbedMultipleConfs` gera coordenadas num referencial arbitrário. A função
   calculava clash contra o receptor e ângulo com a direção de saída sem nunca
   ter POSICIONADO o linker no exit vector — as distâncias eram entre nuvens de
   pontos sem relação espacial. Aqui o confôrmero é ancorado no ponto de saída
   e alinhado à direção antes de qualquer medida, e a rotação livre em torno do
   eixo é amostrada.

3. A contabilidade do cap destruía a montagem do WP3, em silêncio.
   `cap_free_terminus()` troca TODO `[#0]` por metil. Depois disso
   `assemble_full_protac()` procura `[#0]`, não acha, e o RDKit devolve a
   molécula INALTERADA dentro de uma tupla de um elemento — sem exceção, e
   `SanitizeMol` passa. O resultado é um lote inteiro de "PROTACs" que são só o
   recrutador-linker sem warhead nenhum, e nada avisa.
   Aqui os dois pontos ficam rotulados (`[1*]` recrutador, `[2*]` warhead) e a
   montagem é verificada.
"""

from __future__ import annotations

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

# Rótulos dos pontos de conjugação do linker
RECRUITER_ISOTOPE = 1     # [1*] — ponta que vai no recrutador (WP2)
WARHEAD_ISOTOPE = 2       # [2*] — ponta que vai no warhead    (WP3)


# ---------------------------------------------------------------------------
# 1. Pontos de conjugação
# ---------------------------------------------------------------------------
# (nome, SMARTS, índice do átomo âncora dentro do match)
ATTACHMENT_PATTERNS = [
    ("amina_primaria",   "[NX3;H2;!$(N[!#6])][CX4]",      0),
    ("hidroxila",        "[OX2;H1][CX4]",                 0),
    ("acido_carboxilico", "[CX3](=[OX1])[OX2;H1]",        0),
    ("azida",            "[NX1-0]=[NX2+1]=[NX1-1]",       2),
    ("alcino_terminal",  "[CX2;H1]#[CX2]",                0),
    ("halogeneto_alquila", "[CX4][F,Cl,Br,I]",            1),
    ("dummy",            "[#0]",                          0),
]


def _heavy_degree(atom) -> int:
    return sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() > 1)


def find_attachment_points(mol, require_terminal: bool = True):
    """Pontos onde o linker pode ser conjugado.

    require_terminal: exige que o átomo âncora esteja na periferia (um único
    vizinho pesado), que é o que distingue uma alça de ponta de um grupo
    funcional no meio da cadeia. O ácido carboxílico é tratado à parte: a
    âncora é o carbono do carboxilo, que tem dois vizinhos pesados por
    construção.
    """
    points = {}
    for name, smarts, anchor_pos in ATTACHMENT_PATTERNS:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue
        for match in mol.GetSubstructMatches(patt):
            idx = match[anchor_pos]
            atom = mol.GetAtomWithIdx(idx)
            if require_terminal and name != "acido_carboxilico":
                if _heavy_degree(atom) > 1:
                    continue
            points.setdefault(idx, name)
    return [{"atom_idx": i, "type": t} for i, t in sorted(points.items())]


def topological_span(mol, i: int, j: int) -> int:
    return int(Chem.GetDistanceMatrix(mol)[i][j])


def is_bifunctional(mol, min_span: int = 3):
    """True se houver >= 2 pontos de conjugação genuinamente separados.

    min_span evita contar duas alças do mesmo grupo (ex.: os dois oxigênios de
    um carboxilo) como se fossem as duas pontas do linker.
    """
    pts = find_attachment_points(mol)
    for a in range(len(pts)):
        for b in range(a + 1, len(pts)):
            if topological_span(mol, pts[a]["atom_idx"], pts[b]["atom_idx"]) >= min_span:
                return True, (pts[a], pts[b])
    return False, None


def filter_linkers(mols, atom_count_range=(5, 40), rotb_max=15,
                   min_span=3, verbose=True):
    """Substitui `filter_linkers` do WP2: agora exige bifuncionalidade."""
    kept, rejected = [], {"atomos": 0, "rotb": 0, "monofuncional": 0}
    for m in mols:
        if m is None:
            continue
        n = m.GetNumHeavyAtoms()
        if not (atom_count_range[0] <= n <= atom_count_range[1]):
            rejected["atomos"] += 1
            continue
        if rdMolDescriptors.CalcNumRotatableBonds(m) > rotb_max:
            rejected["rotb"] += 1
            continue
        ok, pair = is_bifunctional(m, min_span)
        if not ok:
            rejected["monofuncional"] += 1
            continue
        m.SetProp("attach_1_idx", str(pair[0]["atom_idx"]))
        m.SetProp("attach_1_type", pair[0]["type"])
        m.SetProp("attach_2_idx", str(pair[1]["atom_idx"]))
        m.SetProp("attach_2_type", pair[1]["type"])
        kept.append(m)
    if verbose:
        print(f"{len(kept)}/{len(mols)} linkers aprovados")
        print(f"  rejeitados: {rejected['atomos']} por contagem de átomos, "
              f"{rejected['rotb']} por flexibilidade, "
              f"{rejected['monofuncional']} por não serem bifuncionais")
    return kept


# ---------------------------------------------------------------------------
# 2. Colocação do confôrmero no exit vector
# ---------------------------------------------------------------------------
def _rotation_between(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Matriz que leva o vetor a sobre o vetor b (Rodrigues)."""
    a = a / (np.linalg.norm(a) + 1e-12)
    b = b / (np.linalg.norm(b) + 1e-12)
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    if np.linalg.norm(v) < 1e-9:
        return np.eye(3) if c > 0 else -np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * (1.0 / (1.0 + c))


def _axis_rotation(axis: np.ndarray, angle: float) -> np.ndarray:
    axis = axis / (np.linalg.norm(axis) + 1e-12)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)


def place_conformer_at_exit(mol, conf_id: int, attach_idx: int,
                            exit_point, exit_direction, spin: float = 0.0):
    """Coordenadas do confôrmero ancoradas no exit vector do recrutador.

    O átomo de conjugação vai para `exit_point`; a direção (âncora -> resto do
    linker) é alinhada com `exit_direction`; `spin` gira em torno desse eixo,
    que é o grau de liberdade que sobra depois do alinhamento.
    """
    pos = np.array(mol.GetConformer(conf_id).GetPositions(), dtype=float)
    anchor = pos[attach_idx].copy()

    neigh = [n.GetIdx() for n in mol.GetAtomWithIdx(attach_idx).GetNeighbors()
             if n.GetAtomicNum() > 1]
    ref_target = pos[neigh[0]] if neigh else pos.mean(axis=0)
    body = ref_target - anchor
    if np.linalg.norm(body) < 1e-9:
        body = pos.mean(axis=0) - anchor

    R = _rotation_between(body, np.asarray(exit_direction, dtype=float))
    out = (pos - anchor) @ R.T
    if abs(spin) > 1e-9:
        out = out @ _axis_rotation(np.asarray(exit_direction, dtype=float), spin).T
    return out + np.asarray(exit_point, dtype=float)


def score_conformer_at_exit(mol, conf_id: int, attach_idx: int,
                            exit_point, exit_direction, receptor_coords,
                            clash_dist: float = 2.5, n_spins: int = 12):
    """Avalia um confôrmero JÁ POSICIONADO no exit vector.

    Amostra `n_spins` rotações em torno do eixo de saída e fica com a de menor
    penalidade de clash. Devolve None se o linker não couber em nenhuma.
    """
    rec = np.asarray(receptor_coords, dtype=float) if receptor_coords is not None else None
    best = None
    for k in range(max(1, n_spins)):
        spin = 2 * np.pi * k / max(1, n_spins)
        coords = place_conformer_at_exit(mol, conf_id, attach_idx,
                                         exit_point, exit_direction, spin)
        anchor = coords[attach_idx]
        extension = float(np.max(np.linalg.norm(coords - anchor, axis=1)))

        centroid = coords.mean(axis=0)
        v = centroid - np.asarray(exit_point, dtype=float)
        d = np.asarray(exit_direction, dtype=float)
        cos_angle = float(np.dot(v, d) /
                          (np.linalg.norm(v) * np.linalg.norm(d) + 1e-9))

        if rec is not None and len(rec):
            dist = np.linalg.norm(coords[:, None, :] - rec[None, :, :], axis=2)
            min_dist = float(dist.min())
            n_clashes = int((dist < clash_dist).sum())
        else:
            min_dist, n_clashes = float("inf"), 0

        cand = {"conf_id": conf_id, "spin_rad": round(spin, 3),
                "extension_A": round(extension, 2),
                "cos_angle_to_exit": round(cos_angle, 3),
                "min_dist_to_receptor_A": (round(min_dist, 2)
                                           if np.isfinite(min_dist) else None),
                "n_clashes": n_clashes, "clash": n_clashes > 0}
        if best is None or (cand["n_clashes"], -cand["cos_angle_to_exit"]) < \
                           (best["n_clashes"], -best["cos_angle_to_exit"]):
            best = cand
    return best


def evaluate_linker(mol, attach_idx: int, exit_point, exit_direction,
                    receptor_coords, n_confs: int = 50, n_spins: int = 12,
                    seed: int = 0xC0FFEE):
    """Gera confôrmeros e devolve a melhor colocação de cada um, ordenada."""
    molh = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    cids = AllChem.EmbedMultipleConfs(molh, numConfs=n_confs, params=params)
    if not cids:
        return molh, []
    AllChem.MMFFOptimizeMoleculeConfs(molh, maxIters=1000, mmffVariant="MMFF94s")

    results = []
    for cid in cids:
        s = score_conformer_at_exit(molh, cid, attach_idx, exit_point,
                                    exit_direction, receptor_coords,
                                    n_spins=n_spins)
        if s:
            results.append(s)
    results.sort(key=lambda r: (r["n_clashes"], -r["cos_angle_to_exit"]))
    return molh, results


# ---------------------------------------------------------------------------
# 3. Contabilidade dos pontos de conjugação
# ---------------------------------------------------------------------------
def label_attachment_points(mol, recruiter_idx: int, warhead_idx: int):
    """Troca os dois pontos por dummies rotulados: [1*] recrutador, [2*] warhead.

    A partir daqui cada ponta é rastreável, e a montagem do WP3 pede
    explicitamente o `[2*]` em vez de um `[#0]` qualquer.
    """
    rw = Chem.RWMol(mol)
    for idx, iso in ((recruiter_idx, RECRUITER_ISOTOPE),
                     (warhead_idx, WARHEAD_ISOTOPE)):
        atom = rw.GetAtomWithIdx(idx)
        atom.SetAtomicNum(0)
        atom.SetIsotope(iso)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(0)
        atom.SetFormalCharge(0)
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return out


def find_dummy(mol, isotope: int):
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetIsotope() == isotope:
            return atom.GetIdx()
    return None


def cap_for_md(mol, isotope: int = WARHEAD_ISOTOPE, cap_smiles: str = "C"):
    """Capeia UMA ponta com um placeholder, para a MD de nível (i).

    Devolve a molécula capeada SEM tocar na original: o objeto que segue para o
    WP3 tem de manter o `[2*]` intacto. Foi a confusão entre esses dois papéis
    que produzia PROTACs sem warhead.
    """
    idx = find_dummy(mol, isotope)
    if idx is None:
        raise ValueError(
            f"nenhum dummy [{isotope}*] em {Chem.MolToSmiles(mol)} — esta "
            f"molécula já foi capeada? Capeie sempre uma CÓPIA.")
    capped = AllChem.ReplaceSubstructs(
        Chem.Mol(mol), Chem.MolFromSmarts(f"[{isotope}#0]"),
        Chem.MolFromSmiles(cap_smiles), replacementConnectionPoint=0)[0]
    Chem.SanitizeMol(capped)
    return capped


def assemble_protac(recruiter_linker, warhead_mol,
                    isotope: int = WARHEAD_ISOTOPE):
    """Conjuga o warhead no `[2*]` do recrutador-linker, com verificação.

    Substitui `assemble_full_protac()` do notebook. A diferença que importa:
    se o ponto de conjugação não existir, aqui ISSO É UM ERRO. No original o
    RDKit devolvia a molécula intacta e o pipeline seguia com um "PROTAC" sem
    warhead nenhum.

    O warhead precisa ter o átomo de conjugação no índice 0 — que é como o
    generate_pcsk9_warheads.py grava os SDFs.
    """
    idx = find_dummy(recruiter_linker, isotope)
    if idx is None:
        raise ValueError(
            f"recrutador-linker sem ponto [{isotope}*]: "
            f"{Chem.MolToSmiles(recruiter_linker)}. Se ele já foi capeado para "
            f"a MD, use a cópia não capeada.")

    n_before = recruiter_linker.GetNumHeavyAtoms()
    combo = AllChem.ReplaceSubstructs(
        recruiter_linker, Chem.MolFromSmarts(f"[{isotope}#0]"),
        warhead_mol, replacementConnectionPoint=0)
    full = combo[0]
    Chem.SanitizeMol(full)

    # GetNumHeavyAtoms exclui dummies (número atômico 0), então o [2*] que sai
    # não entra na conta: o total é simplesmente a soma dos dois lados.
    expected = n_before + warhead_mol.GetNumHeavyAtoms()
    if full.GetNumHeavyAtoms() != expected:
        raise ValueError(
            f"montagem inconsistente: {full.GetNumHeavyAtoms()} átomos pesados, "
            f"esperados {expected}")
    if find_dummy(full, isotope) is not None:
        raise ValueError("o ponto de conjugação sobreviveu à montagem")
    return full


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("=== 1. filtro bifuncional ===")
    tests = {
        "PEG3 diamina (bifuncional)": "NCCOCCOCCN",
        "PEG3 mono-amina (só uma alça)": "NCCOCCOCC",
        "amino-acido (2 alcas)": "NCCCCCC(=O)O",
        "amina interna (nao terminal)": "CCC(N)CC",
        "alquil diol": "OCCCCCCO",
    }
    mols = []
    for name, smi in tests.items():
        m = Chem.MolFromSmiles(smi)
        m.SetProp("_Name", name)
        ok, pair = is_bifunctional(m)
        print(f"  {name:32s} bifuncional={str(ok):5s} "
              f"{'pontos: ' + str([p['type'] for p in pair]) if ok else ''}")
        mols.append(m)
    print()
    kept = filter_linkers(mols)

    print("\n=== 2. montagem verificada ===")
    linker = Chem.MolFromSmiles("NCCOCCOCCN")
    pts = find_attachment_points(linker)
    labeled = label_attachment_points(linker, pts[0]["atom_idx"], pts[1]["atom_idx"])
    print("  linker rotulado:", Chem.MolToSmiles(labeled))

    capped = cap_for_md(labeled)
    print("  capeado p/ MD :", Chem.MolToSmiles(capped), "(cópia)")
    print("  original      :", Chem.MolToSmiles(labeled), "(intacto)")

    warhead = Chem.MolFromSmiles("NCCc1cc(Oc2ccccc2)cc(C(F)(F)F)c1")
    full = assemble_protac(labeled, warhead)
    print("  PROTAC montado:", Chem.MolToSmiles(full))

    print("\n=== 3. o modo de falha silencioso, agora ruidoso ===")
    try:
        assemble_protac(capped, warhead)
        print("  ERRO: deveria ter levantado exceção")
    except ValueError as e:
        print(f"  capeado + montar -> ValueError: {str(e)[:70]}...")

    print("\n=== 4. colocação no exit vector ===")
    molh, res = evaluate_linker(
        labeled, attach_idx=pts[0]["atom_idx"],
        exit_point=[0.0, 0.0, 0.0], exit_direction=[0.0, 1.0, 0.0],
        receptor_coords=np.array([[0.0, -3.0, 0.0], [2.0, -3.0, 0.0]]),
        n_confs=10, n_spins=8)
    print(f"  {len(res)} confôrmeros avaliados; melhor:")
    for r in res[:3]:
        print(f"    ext {r['extension_A']:5.1f} Å | cos {r['cos_angle_to_exit']:6.3f} "
              f"| clashes {r['n_clashes']} | dmin {r['min_dist_to_receptor_A']}")

    print("\nautoteste OK")
