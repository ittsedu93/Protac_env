#!/usr/bin/env python
"""
RMSD in-place, com simetria — para validação de redocking.

NÃO use rdMolAlign.GetBestRMS para validar redocking
-----------------------------------------------------
GetBestRMS SUPERPÕE as moléculas antes de medir. Num redocking isso invalida o
teste: a pose e o ligante cristalográfico já estão no mesmo referencial, e o que
se quer saber é o quanto a pose se afastou *no sítio*. Uma pose ancorada no
bolsão errado, a 25 Å do lugar certo, é reportada como RMSD 0,00 (autoteste no
fim deste arquivo):

    rmsd_inplace : 25.00   <- correto
    GetBestRMS   :  0.00   <- "validado" no bolsão errado
    CalcRMS      : 25.00   <- correto

Pior: GetBestRMS MODIFICA a molécula probe, deixando as coordenadas dela
superpostas à referência. Quem chama GetBestRMS e depois reaproveita o objeto da
pose (salvar em disco, medir contatos, alimentar o PRosettaC) passa a trabalhar
com a pose deslocada, não com a que o docking produziu. Se precisar chamar
GetBestRMS, passe uma cópia: Chem.Mol(pose).

rdMolAlign.CalcRMS mede in-place e com simetria, e resolve o caso geral. Este
módulo existe pelo que a CalcRMS não faz: restringir o RMSD a um SUBCONJUNTO de
átomos. Isso importa no 6U26, cujo ligante `063` tem um braço de ~15 Å solto no
solvente: exigir que o docking reproduza esse braço reprovaria um protocolo que
acerta o núcleo ancorado. A validação honesta mede o núcleo enterrado e reporta
o valor do ligante inteiro à parte.
"""

from __future__ import annotations

import numpy as np
from rdkit import Chem


def _heavy_atom_mol(mol):
    return Chem.RemoveHs(Chem.Mol(mol))


def symmetry_mappings(probe, ref, max_matches: int = 5000):
    """Mapeamentos probe->ref equivalentes por simetria do grafo molecular."""
    matches = ref.GetSubstructMatches(probe, uniquify=False,
                                      useChirality=False, maxMatches=max_matches)
    if not matches:
        # grafos diferentes (percepção de ligação divergente): tenta pelo ref
        matches = probe.GetSubstructMatches(ref, uniquify=False,
                                            useChirality=False,
                                            maxMatches=max_matches)
        if not matches:
            raise ValueError(
                "probe e referência não casam como subestrutura — a percepção "
                "de ligações divergiu. Gere os dois lados pela mesma rota "
                "(ex.: ambos via obabel a partir do mesmo SDF)."
            )
        matches = [tuple(np.argsort(m)) for m in matches]
    return matches


def rmsd_inplace(probe, ref, atom_indices: list[int] | None = None,
                 heavy_only: bool = True) -> float:
    """RMSD entre probe e ref SEM superposição, considerando simetria.

    atom_indices: se dado, restringe o cálculo a esses índices da *probe*
        (use para medir só o núcleo enterrado, ignorando um braço flexível
        exposto ao solvente que o docking não tem como reproduzir).
    """
    if heavy_only:
        probe, ref = _heavy_atom_mol(probe), _heavy_atom_mol(ref)
    if probe.GetNumAtoms() != ref.GetNumAtoms():
        raise ValueError(
            f"contagem de átomos difere: probe {probe.GetNumAtoms()}, "
            f"ref {ref.GetNumAtoms()}")

    pp = probe.GetConformer().GetPositions()
    rp = ref.GetConformer().GetPositions()

    sel = list(range(probe.GetNumAtoms())) if atom_indices is None else list(atom_indices)
    if not sel:
        raise ValueError("seleção de átomos vazia")

    best = np.inf
    for mapping in symmetry_mappings(probe, ref):
        m = np.asarray(mapping)
        diff = pp[sel] - rp[m[sel]]
        best = min(best, float(np.sqrt((diff ** 2).sum(axis=1).mean())))
    return best


def buried_atom_indices(mol, receptor_xyz: np.ndarray,
                        threshold: int = 20, radius: float = 6.0) -> list[int]:
    """Índices (átomos pesados) com pelo menos `threshold` átomos de receptor a
    `radius` Å — o núcleo ancorado, por oposição ao braço no solvente."""
    m = _heavy_atom_mol(mol)
    pos = m.GetConformer().GetPositions()
    d = np.linalg.norm(pos[:, None, :] - receptor_xyz[None, :, :], axis=2)
    return [int(i) for i in np.where((d <= radius).sum(axis=1) >= threshold)[0]]


if __name__ == "__main__":
    # auto-teste: demonstra a diferença contra GetBestRMS
    from rdkit.Chem import AllChem, rdMolAlign

    m = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1CCN"))
    AllChem.EmbedMolecule(m, randomSeed=1)
    a, b = Chem.Mol(m), Chem.Mol(m)
    cb = b.GetConformer()
    for i in range(b.GetNumAtoms()):
        p = cb.GetAtomPosition(i)
        cb.SetAtomPosition(i, (p.x + 25.0, p.y, p.z))

    print("pose deslocada 25 Å do cristal:")
    # medir PRIMEIRO, porque GetBestRMS altera as coordenadas do probe
    print(f"  rmsd_inplace : {rmsd_inplace(a, b):.2f}  <- correto")
    print(f"  GetBestRMS   : {rdMolAlign.GetBestRMS(Chem.Mol(a), b):.2f}"
          f"  <- superpõe: reporta 0 para o bolsão errado")
    print(f"  CalcRMS      : {rdMolAlign.CalcRMS(Chem.Mol(a), b):.2f}"
          f"  <- idem nesta versão")

    assert rmsd_inplace(a, b) > 24.0, "deveria detectar o deslocamento"
    assert rmsd_inplace(a, a) < 1e-6, "molécula contra si mesma deve dar 0"

    # a mutação silenciosa do probe por GetBestRMS
    before = a.GetConformer().GetPositions().mean(axis=0).copy()
    rdMolAlign.GetBestRMS(a, b)
    after = a.GetConformer().GetPositions().mean(axis=0)
    shift = float(np.linalg.norm(after - before))
    print(f"\nGetBestRMS moveu o centroide do probe em {shift:.1f} Å "
          f"(mutação silenciosa)")
    assert shift > 24.0

    print("\nautoteste OK")
