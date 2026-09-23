#!/usr/bin/env python
"""
Motores de docking: AutoDock Vina (API Python) e Uni-Dock (CLI, GPU).

Por que os dois
---------------
Nesta workstation os envs conda pertencem a outro usuário, então instalar o
`vina` no `pf_vs` falha com EnvironmentNotWritableError. O **Uni-Dock já está
instalado** no `pf_vs` e usa a mesma função de scoring do Vina
(`--scoring vina`), em GPU — então não há nada a instalar.

A diferença não é só conveniência. O Vina docka um ligante por chamada; o
Uni-Dock recebe um **lote** de ligantes numa chamada só (`--gpu_batch`). Para
45 warheads × 30 runs isso é a diferença entre 1350 processos sequenciais e 30
chamadas em lote. O `staged_screening` explora isso invertendo os laços: o laço
externo passa a ser o seed, não o ligante.

Uso: as duas funções `dock_batch_*` recebem e devolvem a mesma estrutura, então
o chamador só escolhe o motor.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

VINA = "vina"
UNIDOCK = "unidock"


# ---------------------------------------------------------------------------
# utilidades
# ---------------------------------------------------------------------------
def find_unidock(explicit: str | None = None) -> str:
    """Localiza o executável do Uni-Dock."""
    candidatos = [explicit] if explicit else []
    candidatos += [
        shutil.which("unidock"),
        "/home/soberano/miniconda3/envs/pf_vs/bin/unidock",
        str(Path.home() / "envs/pf_vs/bin/unidock"),
    ]
    for c in candidatos:
        if c and Path(c).exists():
            return c
    raise SystemExit(
        "não achei o executável `unidock`. Ative o env pf_vs, ou passe o "
        "caminho em --unidock-exe. Para conferir:\n"
        "  ls /home/soberano/miniconda3/envs/pf_vs/bin/ | grep -i unidock")


def parse_pdbqt_score(pdbqt: Path) -> float | None:
    """Melhor score do arquivo de pose (primeiro REMARK VINA RESULT)."""
    try:
        for line in Path(pdbqt).read_text().splitlines():
            if line.startswith("REMARK VINA RESULT"):
                return float(line.split()[3])
    except (OSError, ValueError, IndexError):
        return None
    return None


def check_engine(engine: str, unidock_exe: str | None = None) -> str:
    """Confirma que o motor está utilizável antes de gastar horas."""
    if engine == VINA:
        try:
            import vina  # noqa: F401
        except ImportError as exc:
            raise SystemExit(
                f"o módulo `vina` não está disponível neste interpretador "
                f"({exc}). Use --engine unidock, que não precisa de instalação "
                f"nesta máquina.")
        return "vina (API Python)"

    exe = find_unidock(unidock_exe)
    try:
        r = subprocess.run([exe, "--version"], capture_output=True, text=True,
                           timeout=60)
        versao = (r.stdout or r.stderr).strip().splitlines()
        return f"unidock: {versao[0] if versao else exe}"
    except (subprocess.SubprocessError, OSError) as exc:
        raise SystemExit(f"`{exe}` não executou: {exc}")


# ---------------------------------------------------------------------------
# Vina (um ligante por chamada)
# ---------------------------------------------------------------------------
def dock_batch_vina(receptor_pdbqt, ligands: dict, center, box_size, seed,
                    exhaustiveness, out_dir: Path, n_poses: int = 1, **_):
    """ligands: {id: caminho_pdbqt}. Devolve {id: (score, caminho_da_pose)}."""
    from vina import Vina

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    resultados = {}

    v = Vina(sf_name="vina", seed=seed, cpu=0, verbosity=0)
    v.set_receptor(str(receptor_pdbqt))
    v.compute_vina_maps(center=list(center), box_size=list(box_size))

    for lid, lpath in ligands.items():
        pose = out_dir / f"{lid}_out.pdbqt"
        if pose.exists():                                   # retomada
            resultados[lid] = (parse_pdbqt_score(pose), pose)
            continue
        v.set_ligand_from_file(str(lpath))
        v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
        v.write_poses(str(pose), n_poses=n_poses, overwrite=True)
        resultados[lid] = (v.energies(n_poses=n_poses)[0][0], pose)
    return resultados


# ---------------------------------------------------------------------------
# Uni-Dock (lote na GPU)
# ---------------------------------------------------------------------------
def dock_batch_unidock(receptor_pdbqt, ligands: dict, center, box_size, seed,
                       exhaustiveness, out_dir: Path, n_poses: int = 1,
                       unidock_exe: str | None = None, batch_size: int = 40,
                       search_mode: str | None = None, verbose: bool = False):
    """Mesma assinatura de dock_batch_vina, mas docka em lote na GPU.

    `batch_size` fatia o lote: 45 ligantes de uma vez podem estourar a memória
    da GPU, e o erro nesse caso não é óbvio.
    """
    exe = find_unidock(unidock_exe)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pendentes = {lid: p for lid, p in ligands.items()
                 if not (out_dir / f"{Path(p).stem}_out.pdbqt").exists()}

    itens = list(pendentes.items())
    for i in range(0, len(itens), batch_size):
        lote = itens[i:i + batch_size]
        cmd = [
            exe,
            "--receptor", str(receptor_pdbqt),
            "--center_x", str(center[0]),
            "--center_y", str(center[1]),
            "--center_z", str(center[2]),
            "--size_x", str(box_size[0]),
            "--size_y", str(box_size[1]),
            "--size_z", str(box_size[2]),
            "--scoring", "vina",
            "--num_modes", str(n_poses),
            "--seed", str(seed),
            "--dir", str(out_dir),
            "--gpu_batch", *[str(p) for _, p in lote],
        ]
        if search_mode:
            cmd += ["--search_mode", search_mode]
        else:
            cmd += ["--exhaustiveness", str(exhaustiveness)]

        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            raise SystemExit(
                f"unidock falhou (código {r.returncode}).\n"
                f"Comando:\n  {' '.join(cmd[:20])} ... "
                f"[{len(lote)} ligantes]\n"
                f"stderr:\n{r.stderr[-1500:]}\n\n"
                f"Se as flags não baterem com a sua build, rode "
                f"`{exe} --help` e me diga a saída — a CLI do Uni-Dock variou "
                f"entre versões.")
        if verbose and r.stdout:
            print(r.stdout[-400:])

    resultados = {}
    for lid, lpath in ligands.items():
        pose = out_dir / f"{Path(lpath).stem}_out.pdbqt"
        resultados[lid] = (parse_pdbqt_score(pose) if pose.exists() else None,
                           pose)
    return resultados


# ---------------------------------------------------------------------------
def dock_batch(engine: str, **kwargs):
    return (dock_batch_vina if engine == VINA else dock_batch_unidock)(**kwargs)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Verifica os motores disponíveis.")
    ap.add_argument("--unidock-exe")
    args = ap.parse_args()

    for eng in (VINA, UNIDOCK):
        try:
            print(f"  {eng:8s} OK  — {check_engine(eng, args.unidock_exe)}")
        except SystemExit as e:
            print(f"  {eng:8s} indisponível: {str(e).splitlines()[0]}")
