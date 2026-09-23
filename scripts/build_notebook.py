#!/usr/bin/env python
"""
Monta o notebook integrado protac_wp1_wp3_pipeline_v2.ipynb.

Junta, na ordem certa de execução:
  - as células do notebook original que continuam válidas
  - as correções do WP1 (portão de redocking) e do WP2 (linkers, montagem)
  - a etapa 3.0 nova (ancoragem dos warheads na PCSK9 via co-cristal 6U26)

Gerado por script para ser reproduzível: para mudar o notebook, edite aqui e
rode de novo.

    python scripts/build_notebook.py
"""

import json
from pathlib import Path

CELLS = []


def md(src: str):
    CELLS.append({"cell_type": "markdown", "metadata": {},
                  "source": src.strip("\n").splitlines(keepends=True)})


def code(src: str):
    CELLS.append({"cell_type": "code", "execution_count": None, "metadata": {},
                  "outputs": [], "source": src.strip("\n").splitlines(keepends=True)})


# ===========================================================================
md(r"""
# PROTAC PCSK9 — WP1–WP3 integrado (v2)

Pipeline completo: recrutador VHL/CRBN → linker → warhead PCSK9 → complexo
ternário → MD.

**Este notebook monta o pipeline; as etapas pesadas rodam fora do kernel.**
Elas aparecem como comandos prontos para colar no terminal, em `nohup`/SLURM.
Os envs são diferentes e isso importa: o kernel deste notebook é `mdtools`, mas
o **AutoDock Vina só existe no `pf_vs`**.

## Ordem de execução

Siga `docs/RUNBOOK.md`. Resumo:

| Fase | O quê | Onde |
|---|---|---|
| 0 | **[PORTÃO]** revalidar o redocking do WP1 | terminal `mdtools` |
| 1 | gerar warheads + PDBQT | terminal `mdtools` → `pf_vs` |
| 2 | sítio na PCSK9, **[PORTÃO]** validar, triagem | terminal `mdtools` → `pf_vs` |
| 3 | WP2: linkers, sub-complexos, **[PORTÃO]** validar | notebook |
| 4 | WP3: montar PROTACs, PRosettaC/AF3, MD | notebook + workstation |

## O que mudou em relação à v1

Quatro defeitos corrigidos. Os dois primeiros invalidavam resultados **sem dar
erro**:

1. **`cap_free_terminus` + `assemble_full_protac`** — o cap consumia o `[#0]`;
   a montagem não achava o ponto e o RDKit devolvia a molécula **inalterada**,
   sem exceção. O lote virava recrutador-linker **sem warhead**.
2. **`redocking_rmsd`** — `GetBestRMS` **superpõe antes de medir**: pose a 25 Å
   do sítio correto passava como RMSD 0,00. Também **altera as coordenadas do
   probe**.
3. **`filter_linkers`** — aceitava linker **monofuncional**.
4. **`exit_vector_score`** — media clash **sem nunca posicionar** o linker no
   exit vector.

## Antes de começar

```bash
# o env pf_vs pode não ter pip próprio; instale-o primeiro e use -m pip
conda install -n pf_vs -y -c conda-forge pip
conda run -n pf_vs python -m pip install vina==1.2.7
conda run -n pf_vs python -c "import vina; print('vina ok', vina.__version__)"
```

**Não rode isso numa célula do notebook.** `conda activate` não muda o
interpretador de um kernel que já está no ar, e num cell o IPython transforma
`conda ...` em `%conda`, que devolve `CondaError: Run 'conda init' before
'conda activate'`. Terminal do VS Code, sempre.
""")

# --------------------------------------------------------------------------
md("## Setup — rode esta célula primeiro, sempre")

code(r'''
import json
import subprocess
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import (AllChem, Descriptors, rdMolAlign, rdMolDescriptors)

# ---- ferramentas desta máquina (ajuste se mudar de host) ----
CHIMERAX_EXE = "/usr/bin/chimerax"
OBABEL_EXE = "/home/soberano/miniconda3/envs/obabel_env/bin/obabel"
GMX_EXE = "/usr/local/gromacs/bin/gmx"
ACPYPE_EXE = "/home/soberano/miniconda3/envs/mdtools/bin/acpype"
PF_VS_PYTHON = "/home/soberano/miniconda3/envs/pf_vs/bin/python"
BOLTZ_ENV = "gustavo_boltz-2"
PROSETTAC_DIR = Path("/mnt/hd2tb/Documentos/PRosettaC")

# ---- este repositório (detectado, não fixo: funciona em qualquer clone) ----
def _find_repo() -> Path:
    marcador = Path("scripts") / "rmsd_inplace.py"
    candidatos = [Path.cwd(), *Path.cwd().parents,
                  Path.home() / "Protac_env", Path.home() / "protac_env"]
    for c in candidatos:
        if (c / marcador).exists():
            return c
    raise FileNotFoundError(
        "não achei o repositório (procurei por scripts/rmsd_inplace.py). "
        "Defina REPO_DIR à mão na célula de Setup.")


REPO_DIR = _find_repo()
SCRIPTS_DIR = REPO_DIR / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

from rmsd_inplace import rmsd_inplace, buried_atom_indices
from wp2_linker_tools import (
    find_attachment_points, is_bifunctional, filter_linkers,
    label_attachment_points, cap_for_md, assemble_protac,
    evaluate_linker, place_conformer_at_exit, find_dummy,
    RECRUITER_ISOTOPE, WARHEAD_ISOTOPE,
)

# ---- pasta de trabalho ----
WORK_DIR = Path.home() / "PRosettaC_runs" / "vhl_crbn_pcsk9_protac"
WORK_DIR.mkdir(parents=True, exist_ok=True)

print("REPO_DIR :", REPO_DIR)
print("WORK_DIR :", WORK_DIR)
for exe, nome in ((CHIMERAX_EXE, "chimerax"), (OBABEL_EXE, "obabel"),
                  (GMX_EXE, "gmx"), (ACPYPE_EXE, "acpype")):
    if not Path(exe).exists():
        print(f"  [aviso] {nome} não está em {exe} — ajuste o caminho acima")
''')

# ===========================================================================
md(r"""
---
# Work Package 1 — recrutador de E3 ligase (VHL/CRBN) por SBVS

Triagem virtual baseada em estrutura da coleção Enamine dedicada a PROTAC
contra VHL e CRBN, com validação do protocolo por redocking do ligante
co-cristalizado, triagem em dois estágios (30 → 100 runs independentes) e
identificação do exit vector para conjugação do linker.
""")

md("## 1.1 Configuração")

code(r'''
E3_TARGETS = {
    "VHL": {
        "pdb_ids": ["6GFZ", "6GFY"],
        "chain_receptor": "A",
        "ref_ligand_resname": None,   # preencher após inspecionar o HETATM
    },
    "CRBN": {
        "pdb_ids": ["4CI1", "4CI3", "4TZ4"],
        "chain_receptor": "A",
        "ref_ligand_resname": None,   # talidomida/pomalidomida/lenalidomida
    },
}

# Catálogo Enamine para PROTAC — baixe manualmente e aponte aqui
ENAMINE_PROTAC_BB_DIR = WORK_DIR / "enamine_protac_building_blocks"

FILTERS = dict(mw_max=500.0, hbd_max=5, hba_max=10, rotb_max=10,
               charge_range=(-1, 1))

REDOCK_RMSD_CUTOFF = 2.0
ROUND1_N_RUNS = 30
ROUND2_N_RUNS = 100
EXHAUSTIVENESS = 32
BOX_PAD_A = 4.0
''')

md("## 1.2 Estruturas 3D de VHL e CRBN (PDB)")

code(r'''
def fetch_pdb(pdb_id: str, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{pdb_id.upper()}.pdb"
    if not out_path.exists():
        urllib.request.urlretrieve(
            f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb", out_path)
    return out_path


pdb_paths = {}
for e3, info in E3_TARGETS.items():
    pdb_paths[e3] = {pid: fetch_pdb(pid, WORK_DIR / "pdb" / e3)
                     for pid in info["pdb_ids"]}
    print(e3, list(pdb_paths[e3]))
''')

md(r"""
## 1.3 Preparo do receptor (ChimeraX) e PDBQT

O ligante co-cristalizado é exportado à parte — é ele que valida o protocolo
no passo 1.6.
""")

code(r'''
def write_chimerax_prep_script(pdb_path: Path, ligand_resname: str, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    receptor_pdb = out_dir / f"{pdb_path.stem}_receptor.pdb"
    ligand_pdb = out_dir / f"{pdb_path.stem}_ref_ligand.pdb"
    script = out_dir / f"{pdb_path.stem}_prep.cxc"
    script.write_text(
        f"open {pdb_path}\n"
        f"select :{ligand_resname}\n"
        f"save {ligand_pdb} selectedOnly true\n"
        f"select ~:{ligand_resname} & ~solvent & ~ions\n"
        f"addh\n"
        f"save {receptor_pdb} selectedOnly true\n"
        f"exit\n"
    )
    return script, receptor_pdb, ligand_pdb


def prepare_receptor_pdbqt(receptor_pdb: Path, out_pdbqt: Path):
    """obabel: simples e testado. O preparo de RECEPTOR no Meeko varia por
    versão instalada, então não dependemos dele aqui."""
    subprocess.run([OBABEL_EXE, str(receptor_pdb), "-O", str(out_pdbqt),
                    "-xr", "-p", "7.4"], check=True)
    return out_pdbqt


def ref_ligand_to_sdf(ligand_pdb: Path) -> Path:
    """Referência do redocking. Converta a pose pela MESMA rota, senão a
    percepção de ligações diverge e o RMSD não consegue casar os grafos."""
    sdf = ligand_pdb.with_suffix(".sdf")
    subprocess.run([OBABEL_EXE, str(ligand_pdb), "-O", str(sdf), "-h"], check=True)
    return sdf


# Exemplo (repita para cada pdb_id):
# script, rec_pdb, lig_pdb = write_chimerax_prep_script(
#     pdb_paths["VHL"]["6GFZ"], E3_TARGETS["VHL"]["ref_ligand_resname"],
#     WORK_DIR / "prep" / "VHL_6GFZ")
# subprocess.run([CHIMERAX_EXE, "--nogui", str(script)], check=True)
# prepare_receptor_pdbqt(rec_pdb, rec_pdb.with_suffix(".pdbqt"))
# ref_ligand_to_sdf(lig_pdb)
''')

md("## 1.4 Biblioteca Enamine — carregar e filtrar")

code(r'''
def load_library(sdf_dir: Path):
    mols = []
    for sdf in sorted(Path(sdf_dir).glob("*.sdf")):
        for mol in Chem.SDMolSupplier(str(sdf), removeHs=False):
            if mol is not None:
                mols.append(mol)
    return mols


def passes_filters(mol, mw_max, hbd_max, hba_max, rotb_max, charge_range):
    return (Descriptors.MolWt(mol) <= mw_max
            and rdMolDescriptors.CalcNumHBD(mol) <= hbd_max
            and rdMolDescriptors.CalcNumHBA(mol) <= hba_max
            and rdMolDescriptors.CalcNumRotatableBonds(mol) <= rotb_max
            and charge_range[0] <= Chem.GetFormalCharge(mol) <= charge_range[1])


def filter_library(mols, filters):
    kept = [m for m in mols if passes_filters(m, **filters)]
    print(f"{len(kept)}/{len(mols)} moléculas passaram nos filtros")
    return kept


# library = load_library(ENAMINE_PROTAC_BB_DIR)
# library_filtered = filter_library(library, FILTERS)
''')

md(r"""
## 1.5 Preparo de ligantes e docking

`run_vina_once` faz `from vina import Vina`, que **só existe no env `pf_vs`**.
Para rodar a triagem, empacote como script e chame
`conda run -n pf_vs python ...`, ou abra este notebook com o kernel `pf_vs`.
""")

code(r'''
def embed_3d(mol):
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def ligand_to_pdbqt(mol, out_pdbqt: Path):
    tmp_sdf = out_pdbqt.with_suffix(".sdf")
    with Chem.SDWriter(str(tmp_sdf)) as w:
        w.write(mol)
    subprocess.run(
        [PF_VS_PYTHON, "-c",
         "from meeko import MoleculePreparation, PDBQTWriterLegacy; "
         "from rdkit import Chem; "
         f"mol = next(Chem.SDMolSupplier('{tmp_sdf}', removeHs=False)); "
         "prep = MoleculePreparation(); setups = prep.prepare(mol); "
         "s = PDBQTWriterLegacy.write_string(setups[0])[0]; "
         f"open('{out_pdbqt}', 'w').write(s)"],
        check=True)
    return out_pdbqt


def run_vina_once(receptor_pdbqt, ligand_pdbqt, center, box_size, seed,
                  exhaustiveness, n_poses, out_pose_pdbqt):
    from vina import Vina           # env pf_vs
    v = Vina(sf_name="vina", seed=seed, cpu=0)
    v.set_receptor(str(receptor_pdbqt))
    v.set_ligand_from_file(str(ligand_pdbqt))
    v.compute_vina_maps(center=list(center), box_size=list(box_size))
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    v.write_poses(str(out_pose_pdbqt), n_poses=n_poses, overwrite=True)
    return v.energies(n_poses=n_poses)[:, 0].tolist()


def independent_docking_runs(receptor_pdbqt, ligand_pdbqt, center, box_size,
                             n_runs, exhaustiveness, out_dir, seed_start=1):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for i in range(n_runs):
        seed = seed_start + i
        pose = out_dir / f"run{seed:03d}.pdbqt"
        scores = run_vina_once(receptor_pdbqt, ligand_pdbqt, center, box_size,
                               seed, exhaustiveness, 1, pose)
        rows.append({"run": i + 1, "seed": seed, "best_score": scores[0],
                     "pose": pose})
    return pd.DataFrame(rows)
''')

md(r"""
## 1.6 Validação do protocolo — **CORRIGIDO**

A versão original usava `rdMolAlign.GetBestRMS`, que **superpõe as moléculas
antes de medir**. Uma pose ancorada a 25 Å do sítio correto era reportada como
RMSD 0,00 — o portão aprovava qualquer coisa. Demonstração:

```bash
python scripts/rmsd_inplace.py
```
```
rmsd_inplace : 25.00   <- correto
GetBestRMS   :  0.00   <- "VALIDADO" no bolsão errado
CalcRMS      : 25.00   <- correto
```

Segundo efeito: `GetBestRMS` **altera as coordenadas do probe**. Quem media e
depois reaproveitava o objeto da pose seguia com ela deslocada.
""")

code(r'''
def redocking_rmsd(pose_path, ref_ligand_path, receptor_pdb=None,
                   core_threshold: int = 20) -> dict:
    """Substitui a versão que usava GetBestRMS. Mede sem superposição, com
    simetria. Com `receptor_pdb`, devolve também o RMSD restrito ao núcleo
    enterrado — o número que decide quando a referência tem partes flexíveis
    expostas ao solvente que o docking não tem como reproduzir."""
    def _load(p):
        p = Path(p)
        if p.suffix.lower() == ".pdbqt":
            sdf = p.with_suffix(".conv.sdf")
            subprocess.run([OBABEL_EXE, str(p), "-O", str(sdf)],
                           check=True, capture_output=True)
            p = sdf
        if p.suffix.lower() in (".sdf", ".mol"):
            return Chem.MolFromMolFile(str(p), removeHs=True)
        if p.suffix.lower() == ".pdb":
            return Chem.MolFromPDBFile(str(p), removeHs=True)
        raise ValueError(f"formato não suportado: {p.suffix}")

    probe, ref = _load(pose_path), _load(ref_ligand_path)
    if probe is None or ref is None:
        raise ValueError(f"falha ao ler {pose_path} ou {ref_ligand_path}")

    out = {"rmsd_full_A": rmsd_inplace(probe, ref)}
    if receptor_pdb:
        rec = np.array([(float(l[30:38]), float(l[38:46]), float(l[46:54]))
                        for l in Path(receptor_pdb).read_text().splitlines()
                        if l.startswith("ATOM") and (l[76:78].strip() or "C") != "H"])
        core = buried_atom_indices(ref, rec, threshold=core_threshold)
        if len(core) >= 3:
            out["rmsd_core_A"] = rmsd_inplace(probe, ref, atom_indices=core)
            out["n_core_atoms"] = len(core)
    return out


def validate_docking_protocol(receptor_pdbqt, ref_ligand_pdbqt, ref_ligand_sdf,
                              center, box_size, out_dir, receptor_pdb=None):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pose = out_dir / "redock_pose.pdbqt"
    scores = run_vina_once(receptor_pdbqt, ref_ligand_pdbqt, center, box_size,
                           seed=1, exhaustiveness=EXHAUSTIVENESS, n_poses=1,
                           out_pose_pdbqt=pose)
    r = redocking_rmsd(pose, ref_ligand_sdf, receptor_pdb)
    decisivo = r.get("rmsd_core_A", r["rmsd_full_A"])
    ok = decisivo <= REDOCK_RMSD_CUTOFF
    print(f"score {scores[0]:.2f} kcal/mol | RMSD {decisivo:.2f} Å "
          f"(corte {REDOCK_RMSD_CUTOFF}) -> {'VALIDADO' if ok else 'NAO VALIDADO'}")
    return ok, r
''')

md(r"""
### 1.6b **[PORTÃO]** Revalidar as poses que já passaram

Se o WP1 já rodou com a função antiga, as conclusões precisam ser reconferidas.
Este passo **não redoca** — só remede o que está no disco.

Preencha `WP1_ALVOS` com os caminhos reais e rode os comandos no terminal.
""")

code(r'''
WP1_ALVOS = {
    # "VHL_6GFZ": dict(
    #     receptor   = WORK_DIR / "prep/VHL_6GFZ/6GFZ_receptor.pdb",
    #     ref_ligand = WORK_DIR / "prep/VHL_6GFZ/6GFZ_ref_ligand.sdf",
    #     poses      = str(WORK_DIR / "prep/VHL_6GFZ/redock*.pdbqt"),
    # ),
    # "CRBN_4CI1": dict(receptor=..., ref_ligand=..., poses=...),
}

if not WP1_ALVOS:
    print("WP1_ALVOS vazio — preencha com os caminhos das poses já geradas.")
for nome, cfg in WP1_ALVOS.items():
    print(f"python {SCRIPTS_DIR}/revalidate_redocking.py \\")
    print(f"    --receptor-pdb {cfg['receptor']} \\")
    print(f"    --ref-ligand   {cfg['ref_ligand']} \\")
    print(f"    --poses        '{cfg['poses']}' \\")
    print(f"    --label        {nome} \\")
    print(f"    --out          {WORK_DIR}/wp1_revalidation_{nome}.csv\n")
''')

code(r'''
def consolidar_revalidacao(work_dir: Path = None):
    """Rode depois dos comandos acima."""
    work_dir = Path(work_dir or WORK_DIR)
    csvs = sorted(work_dir.glob("wp1_revalidation_*.csv"))
    if not csvs:
        print("nenhum CSV de revalidação ainda")
        return None
    df = pd.concat([pd.read_csv(c) for c in csvs], ignore_index=True)
    col = "rmsd_core_A" if df["rmsd_core_A"].notna().any() else "rmsd_full_A"
    resumo = df.groupby("label").agg(poses=("pose", "count"),
                                     passaram=("passa_corte", "sum"),
                                     melhor_rmsd=(col, "min"))
    print(resumo)
    reprovados = resumo[resumo["passaram"] == 0]
    if len(reprovados):
        print(f"\nATENÇÃO: sem nenhuma pose válida em {list(reprovados.index)}")
        print("A triagem do WP1 nesses alvos não tem protocolo validado.")
    return df


# consolidar_revalidacao()
''')

md("## 1.7 Triagem em dois estágios (30 → 100 runs)")

code(r'''
def score_summary(df: pd.DataFrame) -> dict:
    return {"best_score": df["best_score"].min(),
            "mean_score": df["best_score"].mean(),
            "std_score": df["best_score"].std(),   # consistência entre runs
            "n_runs": len(df)}


def staged_screening(ligands: dict, receptor_pdbqt, center, box_size, out_dir,
                     round1_cutoff_score):
    out_dir = Path(out_dir)
    rows1 = []
    for lig_id, lig in ligands.items():
        df = independent_docking_runs(receptor_pdbqt, lig, center, box_size,
                                      ROUND1_N_RUNS, EXHAUSTIVENESS,
                                      out_dir / "round1" / lig_id)
        s = score_summary(df); s["ligand_id"] = lig_id
        rows1.append(s)
    round1 = pd.DataFrame(rows1).sort_values("best_score")

    survivors = round1[round1["best_score"] <= round1_cutoff_score]["ligand_id"].tolist()
    print(f"Round 1: {len(survivors)}/{len(ligands)} abaixo de {round1_cutoff_score}")

    rows2 = []
    for lig_id in survivors:
        df = independent_docking_runs(receptor_pdbqt, ligands[lig_id], center,
                                      box_size, ROUND2_N_RUNS, EXHAUSTIVENESS,
                                      out_dir / "round2" / lig_id)
        s = score_summary(df); s["ligand_id"] = lig_id
        rows2.append(s)
    return round1, pd.DataFrame(rows2).sort_values("best_score")
''')

md(r"""
## 1.8 Inspeção estrutural e exit vector

O exit vector do recrutador é anotado **manualmente** aqui: é a posição
acessível ao solvente por onde o linker vai sair. Guarde `exit_point` (as
coordenadas do átomo) e `exit_direction` (o vetor que aponta para fora) — o
WP2 usa os dois.
""")

code(r'''
def write_pymol_inspection_script(receptor_pdb: Path, pose_pdbqts: list,
                                  out_pml: Path):
    lines = [f"load {receptor_pdb}, receptor"]
    for i, pose in enumerate(pose_pdbqts):
        pose_pdb = Path(pose).with_suffix(".pdb")
        subprocess.run([OBABEL_EXE, str(pose), "-O", str(pose_pdb)], check=True)
        lines.append(f"load {pose_pdb}, hit_{i:02d}")
    lines += ["hide everything", "show cartoon, receptor", "show sticks, hit_*",
              "util.cbag hit_*", "set stick_radius, 0.15"]
    Path(out_pml).write_text("\n".join(lines) + "\n")
    print("Abra com: pymol", out_pml)


# Anote aqui, após a inspeção (um por recrutador priorizado):
RECRUITER_EXIT_VECTORS = {
    # "VHL_hit01": dict(exit_point=[x, y, z], exit_direction=[dx, dy, dz],
    #                   atom_note="carbono do anel apontando para o solvente"),
}
''')

md("**Milestone WP1:** recrutadores VHL/CRBN priorizados, com protocolo "
   "validado pelo critério correto e exit vector anotado.")

# ===========================================================================
md(r"""
---
# Work Package 2 — linker e integração recrutador–linker

Seleção de linkers da coleção Enamine por contagem de átomos, flexibilidade e
**bifuncionalidade**, avaliação geométrica a partir do exit vector do WP1,
docking restrito e modelagem complementar em Boltz-2.
""")

md("## 2.1 Configuração")

code(r'''
ENAMINE_LINKER_DIR = WORK_DIR / "enamine_protac_linkers"

LINKER_FILTERS = dict(atom_count_range=(5, 40), rotb_max=15, min_span=3)
N_CONFORMERS = 50
N_SPINS = 12              # rotações amostradas em torno do eixo de saída
CLASH_DISTANCE_A = 2.5

WP2_DIR = WORK_DIR / "wp2_subcomplexes"
WP2_DIR.mkdir(parents=True, exist_ok=True)
''')

md(r"""
## 2.2 Filtrar linkers — **CORRIGIDO**

A versão original usava `has_terminal_group()`, que devolve True com **um**
match de `[NH2]`/`[OH]`. Um linker de PROTAC é bifuncional por definição: uma
ponta no recrutador, a outra no warhead. Além disso `[NH2]` casa amina **não
terminal** (`CC(N)CC` casa), então "grupo terminal" não era terminal.

`filter_linkers` de `wp2_linker_tools` exige dois pontos de conjugação com
separação topológica mínima, e grava em cada molécula as props
`attach_1_idx` / `attach_2_idx`.
""")

code(r'''
# linkers = load_library(ENAMINE_LINKER_DIR)
# linkers_ok = filter_linkers(linkers, **LINKER_FILTERS)

# Demonstração do filtro:
for smi, rotulo in [("NCCOCCOCCN", "PEG diamina"),
                    ("NCCOCCOCC",  "PEG mono-amina"),
                    ("CCC(N)CC",   "amina interna")]:
    ok, par = is_bifunctional(Chem.MolFromSmiles(smi))
    print(f"{rotulo:18s} {smi:14s} bifuncional={ok}")
''')

md(r"""
## 2.3 Geometria a partir do exit vector — **CORRIGIDO**

A versão original calculava clash contra o receptor sobre as coordenadas do
`EmbedMultipleConfs`, que estão num referencial arbitrário — o linker **nunca
era posicionado** no exit vector, então as distâncias eram entre nuvens de
pontos sem relação espacial.

Agora cada confôrmero é ancorado no ponto de saída, alinhado à direção, e a
rotação livre em torno do eixo é amostrada (`N_SPINS`).
""")

code(r'''
def receptor_coords_from_pdb(receptor_pdb: Path) -> np.ndarray:
    return np.array([(float(l[30:38]), float(l[38:46]), float(l[46:54]))
                     for l in Path(receptor_pdb).read_text().splitlines()
                     if l.startswith("ATOM") and (l[76:78].strip() or "C") != "H"])


def avaliar_linkers(linkers, exit_point, exit_direction, receptor_pdb,
                    n_confs=N_CONFORMERS, n_spins=N_SPINS):
    rec = receptor_coords_from_pdb(receptor_pdb)
    linhas = []
    for m in linkers:
        attach = int(m.GetProp("attach_1_idx"))
        _, res = evaluate_linker(m, attach, exit_point, exit_direction, rec,
                                 n_confs=n_confs, n_spins=n_spins)
        if not res:
            continue
        melhor = res[0]
        sem_clash = [r for r in res if r["n_clashes"] == 0]
        linhas.append({"smiles": Chem.MolToSmiles(m),
                       "n_heavy": m.GetNumHeavyAtoms(),
                       "extension_A": melhor["extension_A"],
                       "cos_to_exit": melhor["cos_angle_to_exit"],
                       "n_clashes_melhor": melhor["n_clashes"],
                       "min_dist_A": melhor["min_dist_to_receptor_A"],
                       "frac_sem_clash": round(len(sem_clash) / len(res), 3)})
    df = pd.DataFrame(linhas)
    if len(df):
        df = df.sort_values(["frac_sem_clash", "cos_to_exit"], ascending=False)
    return df


# ev = RECRUITER_EXIT_VECTORS["VHL_hit01"]
# df_linkers = avaliar_linkers(linkers_ok, ev["exit_point"], ev["exit_direction"],
#                              WORK_DIR / "prep/VHL_6GFZ/6GFZ_receptor.pdb")
''')

md(r"""
**Priorize por `frac_sem_clash`, não pelo melhor confôrmero isolado.** Um
linker que encaixa em 1 de 50 confôrmeros está geometricamente forçado, por
melhor que pontue esse confôrmero.
""")

md("## 2.4 Docking restrito à pose validada do WP1")

code(r'''
def dock_linker_constrained(receptor_pdbqt, linker_pdbqt, wp1_exit_point,
                            box_size, out_dir, n_runs=ROUND1_N_RUNS):
    """Box pequeno centrado no exit vector já validado — não faz nova busca de
    sítio. Aqui o critério é adequação espacial, não afinidade pela E3."""
    return independent_docking_runs(receptor_pdbqt, linker_pdbqt,
                                    center=wp1_exit_point, box_size=box_size,
                                    n_runs=n_runs, exhaustiveness=EXHAUSTIVENESS,
                                    out_dir=out_dir)
''')

md(r"""
## 2.5 Montar o sub-complexo — **CORRIGIDO**

Aqui estava o defeito mais grave do pipeline. `cap_free_terminus()` trocava
**todo** `[#0]` por metil; depois disso `assemble_full_protac()` procurava
`[#0]`, não achava, e o RDKit devolvia a molécula **inalterada** numa tupla de
um elemento — sem exceção, com `SanitizeMol` passando. O lote inteiro virava
recrutador-linker **sem warhead nenhum**, em silêncio.

Agora as pontas são rotuladas (`[1*]` recrutador, `[2*]` warhead) e são
gravadas **duas** versões: a do WP3, com `[2*]` intacto, e uma **cópia**
capeada só para a MD de nível (i).
""")

code(r'''
def preparar_subcomplexo(linker, recruiter_mol, sc_id: str, out_dir: Path = None):
    """recruiter_mol precisa ter o átomo de conjugação no índice 0 (mesma
    convenção dos warheads)."""
    out_dir = Path(out_dir or WP2_DIR) / sc_id
    out_dir.mkdir(parents=True, exist_ok=True)

    a1 = int(linker.GetProp("attach_1_idx"))
    a2 = int(linker.GetProp("attach_2_idx"))
    rotulado = label_attachment_points(linker, a1, a2)

    sub = AllChem.ReplaceSubstructs(
        rotulado, Chem.MolFromSmarts(f"[{RECRUITER_ISOTOPE}#0]"),
        recruiter_mol, replacementConnectionPoint=0)[0]
    Chem.SanitizeMol(sub)

    if find_dummy(sub, WARHEAD_ISOTOPE) is None:
        raise ValueError(f"{sc_id}: o [2*] sumiu ao conjugar o recrutador")

    with Chem.SDWriter(str(out_dir / f"{sc_id}.sdf")) as w:
        w.write(sub)                       # -> WP3, [2*] intacto
    with Chem.SDWriter(str(out_dir / f"{sc_id}_capped_md.sdf")) as w:
        w.write(cap_for_md(sub))           # -> MD nível (i), cópia capeada

    return {"subcomplex_id": sc_id, "smiles": Chem.MolToSmiles(sub),
            "path_wp3": out_dir / f"{sc_id}.sdf",
            "path_md": out_dir / f"{sc_id}_capped_md.sdf"}
''')

md(r"""
### 2.6 **[PORTÃO]** Validar os sub-complexos

Rode antes de entregar qualquer coisa ao WP3. É o que impede o lote de PROTACs
sem warhead de existir.
""")

code(r'''
def validar_subcomplexos(wp2_dir: Path = None):
    wp2_dir = Path(wp2_dir or WP2_DIR)
    arquivos = [p for p in wp2_dir.rglob("*.sdf") if "_capped_md" not in p.name]
    if not arquivos:
        print(f"nenhum sub-complexo em {wp2_dir}")
        return None
    linhas = []
    for p in arquivos:
        m = Chem.MolFromMolFile(str(p), removeHs=False, sanitize=False)
        if m is None:
            linhas.append({"arquivo": p.name, "ok": False, "motivo": "ilegível"})
            continue
        tem_wh = find_dummy(m, WARHEAD_ISOTOPE) is not None
        tem_rec = find_dummy(m, RECRUITER_ISOTOPE) is not None
        linhas.append({"arquivo": p.name, "ok": tem_wh and not tem_rec,
                       "motivo": ("pronto" if tem_wh and not tem_rec else
                                  "[2*] ausente — foi capeado?" if not tem_wh else
                                  "[1*] sobrou — recrutador não conjugado")})
    df = pd.DataFrame(linhas)
    print(f"{df['ok'].sum()}/{len(df)} sub-complexos prontos para o WP3")
    if not df["ok"].all():
        print(df[~df["ok"]].to_string(index=False))
    return df


# validar_subcomplexos()
''')

md(r"""
## 2.7 Boltz-2 — modelagem complementar

`filter_boltz_models` da v1 caía em `df.columns[-1]` quando não achava
`confidence_score`, filtrando por uma coluna qualquer sem avisar. Aqui a
coluna é obrigatória.
""")

code(r'''
def write_boltz_input_yaml(protein_seq: str, ligand_smiles: str, out_yaml: Path):
    Path(out_yaml).parent.mkdir(parents=True, exist_ok=True)
    Path(out_yaml).write_text(
        "version: 1\nsequences:\n"
        "  - protein:\n      id: A\n"
        f'      sequence: "{protein_seq}"\n'
        "  - ligand:\n      id: B\n"
        f'      smiles: "{ligand_smiles}"\n')
    return out_yaml


def boltz_predict_cmd(input_yaml: Path, out_dir: Path) -> str:
    return (f"conda run -n {BOLTZ_ENV} boltz predict {input_yaml} "
            f"--out_dir {out_dir} --use_msa_server --output_format pdb")


def filter_boltz_models(confidence_json_paths, confidence_cutoff=0.6,
                        col="confidence_score"):
    linhas = [{"model": Path(p).parent.name, **json.loads(Path(p).read_text())}
              for p in confidence_json_paths]
    df = pd.DataFrame(linhas)
    if col not in df.columns:
        raise KeyError(f"coluna '{col}' ausente. Disponíveis: {list(df.columns)}. "
                       f"Escolha explicitamente, sem fallback silencioso.")
    sobrev = df[df[col] >= confidence_cutoff]
    print(f"{len(sobrev)}/{len(df)} modelos acima de {confidence_cutoff}")
    return df, sobrev


def rmsd_to_reference(model_pdb: Path, reference_pdb: Path, selection="protein"):
    """Comparação modelo x referência em referenciais diferentes: aqui a
    superposição é o que se quer (ao contrário do redocking)."""
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms
    u, ref = mda.Universe(str(model_pdb)), mda.Universe(str(reference_pdb))
    return rms.rmsd(u.select_atoms(selection).positions,
                    ref.select_atoms(selection).positions, superposition=True)
''')

md("**Milestone WP2:** sub-complexos recrutador–linker validados, com `[2*]` "
   "livre, prontos para conjugação com o warhead.")

# ===========================================================================
md(r"""
---
# Work Package 3 — warhead PCSK9, complexo ternário e MD

Warheads da série feniletilamina (Figura 3) conjugados aos sub-complexos do
WP2; modelagem do ternário por PRosettaC e AlphaFold 3; MD em três níveis de
complexidade.
""")

md("## 3.1 Configuração — warheads PCSK9 e co-cristal 6U26")

code(r'''
# --- saída de generate_pcsk9_warheads.py ---
WARHEADS_DIR = Path.home() / "PCSK9_warheads"
WARHEADS_CSV = WARHEADS_DIR / "pcsk9_warheads.csv"
WARHEADS_SDF_DIR = WARHEADS_DIR / "sdf"        # átomo 0 = N de conjugação
WARHEADS_PDBQT_DIR = WARHEADS_DIR / "pdbqt"

# --- PCSK9: co-cristal 6U26 (PCSK9 + composto 16, HET 063, 1,53 Å) ---
#   cadeia A = pró-domínio (61-152); cadeia B = catalítico + CHRD (153-682)
#   O ligante 063 ancora na cadeia B, a 23,4 Å da interface EGF(A)/LDLR:
#   NÃO é o sítio clássico de inibição PCSK9-LDLR. Decisão consciente.
PCSK9_PDB_FILE = Path.home() / "structures" / "6U26_1.pdb"
PCSK9_CHAIN = "B"
PCSK9_REF_LIGAND = "063"
PCSK9_DOCKING_DIR = WORK_DIR / "pcsk9_docking"
PCSK9_SITE_JSON = PCSK9_DOCKING_DIR / "pcsk9_site.json"

# --- MD (inalterado da metodologia) ---
MD_LEVELS = ["E3_recruiter_linker", "E3_recruiter_linker_warhead", "full_ternary"]
N_REPLICATES = 3
PRODUCTION_NS = 200
EQUIL_NPT_NS = 5
DT_PS = 0.002
COORD_SAVE_PS = 200
TEMP_K = 310
PRESSURE_ATM = 1.0
BOX_DIST_NM = 1.3
RVDW_RCOULOMB_NM = 0.9

WP3_DIR = WORK_DIR / "wp3"
WP3_DIR.mkdir(parents=True, exist_ok=True)

print("cristal:", PCSK9_PDB_FILE,
      "(ok)" if PCSK9_PDB_FILE.exists() else "*** ponha o 6U26_1.pdb aqui ***")
''')

md(r"""
### 3.1a Gerar e priorizar a série de warheads

Rode uma vez no terminal (env `mdtools`, ~20 min):

```bash
python scripts/generate_pcsk9_warheads.py --outdir ~/PCSK9_warheads --n-confs 50
```

180 moléculas únicas: R1 (8) × R2 (6) × R3 (4), deduplicadas por SMILES
canônico. O gerador **não descarta nada** — marca em `filter_flags`. A
priorização fica abaixo, explícita, para poder ser justificada.
""")

code(r'''
if not WARHEADS_CSV.exists():
    print(f"ainda não existe: {WARHEADS_CSV}")
    print(f"rode: python {SCRIPTS_DIR}/generate_pcsk9_warheads.py "
          f"--outdir {WARHEADS_DIR} --n-confs 50")
    wh = serie_A = serie_B = None
else:
    wh = pd.read_csv(WARHEADS_CSV)
    print(f"{len(wh)} warheads enumerados")
    prioritized = wh[wh["filter_flags"] == "ok"].copy()

    # Série A — R3 = H: N primário, ponto de conjugação limpo. Principal.
    serie_A = prioritized[prioritized["r3"] == "H"]
    # Série B — R3 farmacofórico mantido; o linker ocupa o H restante e o N
    # vira terciário. Só depois de confirmar que o N-H não é essencial.
    serie_B = prioritized[prioritized["r3"].isin(["Ms", "Tet", "CH2Tet"])]

    print(f"Série A (R3=H):      {len(serie_A)}")
    print(f"Série B (R3 mantido): {len(serie_B)}")
    display(serie_A[["warhead_id", "r1", "r2", "mw", "clogp",
                     "net_charge_ph74"]].head(10))
''')

md(r"""
## 3.0 **(rode antes do 3.2)** Ancorar os warheads na PCSK9

**Por que esta etapa não existia na v1:** a metodologia partia dos ligantes
GMF, que já tinham pose conhecida. Os warheads da Figura 3 são moléculas novas,
e o PRosettaC exige `Heads` **já posicionados** dentro da estrutura.

O co-cristal 6U26 resolve os três pontos:

| precisa de | 6U26 fornece |
|---|---|
| onde dockar | bolsão do ligante `063`, medido da estrutura |
| como validar | redocking do `063`, corte RMSD ≤ 2,0 Å |
| para onde o linker sai | braço PEG/guanidina, 17 átomos sem vizinho proteico |
""")

md("### 3.0a Receptor e sítio (terminal, env `mdtools`, minutos)")

code(r'''
print("# 1) inspeção — confira as cadeias antes de qualquer coisa:")
print(f"python {SCRIPTS_DIR}/prep_pcsk9_receptor.py \\")
print(f"    --pdb-file {PCSK9_PDB_FILE} \\")
print(f"    --outdir {PCSK9_DOCKING_DIR} --inspect-only")

print("\n# 2) sítio + receptor + ligante de referência:")
print(f"python {SCRIPTS_DIR}/prep_pcsk9_receptor.py \\")
print(f"    --pdb-file {PCSK9_PDB_FILE} \\")
print(f"    --outdir {PCSK9_DOCKING_DIR} \\")
print(f"    --target-chain {PCSK9_CHAIN} --keep-chains A B \\")
print(f"    --site-mode ligand --ref-ligand-resname {PCSK9_REF_LIGAND} \\")
print(f"    --burial-threshold 20")
''')

md(r"""
`--keep-chains A B` mantém o pró-domínio: ele fica a 9,1 Å do ligante e compõe
a borda do bolsão.

`--burial-threshold 20` não é cosmético. O `063` tem 76 átomos pesados, ~17
deles num braço de 15 Å no solvente. Com o ligante inteiro o box sai com ~28 Å
de aresta e o docking vira busca cega; restrito ao núcleo fica em
~20 × 18 × 18 Å, tamanho certo para 230–420 Da.

**Inspecione o box no ChimeraX antes de gastar horas.**
""")

code(r'''
site = None
if PCSK9_SITE_JSON.exists():
    site = json.loads(PCSK9_SITE_JSON.read_text())
    print("centro         :", site["center"])
    print("box (Å)        :", site["box_size"])
    print("vetor de saída :", site["exit_vector"],
          f"(braço de {site['exit_arm_length_A']} Å)")
    print("contatos       :", len(site["contact_residues"]), "resíduos")
    print("origem         :", site["provenance"])
    print(f"\nchimerax {site['receptor_pdb']} {site['ref_ligand_sdf']}")
else:
    print(f"ainda não existe: {PCSK9_SITE_JSON}")
    print("rode os comandos da célula 3.0a no terminal primeiro.")
''')

md(r"""
### 3.0b **[PORTÃO]** Validar o protocolo (terminal, env `pf_vs`)

`--validate-only` faz só o redocking do `063` e para. Se o RMSD de núcleo não
fechar em ≤ 2,0 Å, não adianta dockar 45 warheads — o protocolo não reproduz
nem a pose que o cristal já entregou. O script recusa seguir.

Reprovando, mexa nesta ordem: (a) protonação do receptor, (b) tamanho do box,
(c) `--exhaustiveness`.
""")

code(r'''
print("conda activate pf_vs")
print(f"python {SCRIPTS_DIR}/dock_warheads_pcsk9.py \\")
print(f"    --site {PCSK9_SITE_JSON} \\")
print(f"    --warheads-pdbqt {WARHEADS_PDBQT_DIR} \\")
print(f"    --outdir {PCSK9_DOCKING_DIR / 'docking'} --validate-only")
''')

md(r"""
### 3.0c Triagem em dois estágios (terminal, env `pf_vs`, horas)

30 runs por warhead → corte por score → 100 runs nos sobreviventes.
**Retomável**: poses já gravadas são puladas.
""")

code(r'''
print(f"nohup python {SCRIPTS_DIR}/dock_warheads_pcsk9.py \\")
print(f"    --site {PCSK9_SITE_JSON} \\")
print(f"    --warheads-pdbqt {WARHEADS_PDBQT_DIR} \\")
print(f"    --warheads-csv {WARHEADS_CSV} --series A \\")
print(f"    --outdir {PCSK9_DOCKING_DIR / 'docking'} \\")
print(f"    > {PCSK9_DOCKING_DIR}/docking.log 2>&1 &")
''')

code(r'''
def load_docking_results(docking_dir: Path = None):
    docking_dir = Path(docking_dir or PCSK9_DOCKING_DIR / "docking")
    r2 = pd.read_csv(docking_dir / "round2_scores.csv")
    heads = pd.DataFrame(json.loads(
        (docking_dir / "heads_manifest.json").read_text()))
    return (r2.merge(heads, left_on="ligand_id", right_on="warhead_id")
              .sort_values("best_score"))


# results = load_docking_results()
# display(results[["warhead_id", "best_score", "std_score", "head_sdf"]].head(15))
''')

md(r"""
### 3.0d Conferir as poses antes de gastar dias de PRosettaC

Duas checagens que custam minutos.
""")

code(r'''
def check_pose_consistency(results: pd.DataFrame, sd_max: float = 1.0):
    """Desvio alto entre runs independentes = pose instável. Score bom com
    dispersão alta é candidato a falso positivo, não a hit."""
    instaveis = results[results["std_score"] > sd_max]
    print(f"{len(instaveis)}/{len(results)} com desvio > {sd_max} kcal/mol")
    return instaveis


def check_exit_vector_accessible(head_sdf: Path, site: dict, receptor_pdb: Path,
                                 min_clearance: float = 3.0):
    """O N de conjugação aponta para o solvente, na direção do braço do 063?
    Se ficar enterrado, não há por onde o linker sair e o candidato morre no
    WP2 — melhor descobrir agora."""
    mol = Chem.MolFromMolFile(str(head_sdf), removeHs=True)
    if mol is None:
        return None
    match = mol.GetSubstructMatches(Chem.MolFromSmarts("c[CH2][CH2][NX3]"))
    if not match:
        return None
    n_pos = np.array(mol.GetConformer().GetAtomPosition(match[0][-1]))

    rec = receptor_coords_from_pdb(receptor_pdb)
    clearance = float(np.linalg.norm(rec - n_pos, axis=1).min())
    v = n_pos - np.array(site["center"])
    cos = float(np.dot(v / (np.linalg.norm(v) + 1e-9),
                       np.array(site["exit_vector"])))
    return {"n_clearance_A": round(clearance, 2),
            "cos_to_exit_vector": round(cos, 3),
            "solvent_exposed": clearance >= min_clearance,
            "aligned_with_063_arm": cos > 0.3}
''')

md(r"""
## 3.2 Montar o PROTAC completo — **CORRIGIDO**

`assemble_protac()` procura o `[2*]` rotulado, confere a contagem de átomos e
**levanta exceção** se o ponto não existir. A versão antiga recebia de volta a
molécula inalterada e seguia em silêncio.

Os warheads já saem do gerador com o **átomo 0 sendo o N de conjugação**, que é
o que `replacementConnectionPoint=0` espera.
""")

code(r'''
def load_warhead_for_assembly(warhead_id: str, sdf_dir: Path = None):
    sdf_dir = Path(sdf_dir or WARHEADS_SDF_DIR)
    mol = Chem.MolFromMolFile(str(sdf_dir / f"{warhead_id}.sdf"))
    if mol is None:
        raise FileNotFoundError(f"não li o SDF de {warhead_id}")
    assert mol.GetAtomWithIdx(0).GetSymbol() == "N", \
        f"{warhead_id}: átomo 0 deveria ser o N de conjugação"
    return mol


def build_protac_matrix(subcomplex_mols: dict, warhead_ids: list,
                        out_dir: Path = None):
    out_dir = Path(out_dir or WP3_DIR / "protacs")
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for sc_id, sc_mol in subcomplex_mols.items():
        for wid in warhead_ids:
            try:
                full = assemble_protac(sc_mol, load_warhead_for_assembly(wid))
            except ValueError as exc:
                print(f"  [falha] {sc_id} x {wid}: {exc}")
                continue
            cand_id = f"{sc_id}__{wid}"
            smi = Chem.MolToSmiles(full)
            cand_dir = out_dir / cand_id
            cand_dir.mkdir(parents=True, exist_ok=True)
            (cand_dir / "protac.smi").write_text(f"{smi}\t{cand_id}\n")
            rows.append({"candidate_id": cand_id, "subcomplex_id": sc_id,
                         "warhead_id": wid, "protac_smiles": smi,
                         "protac_mw": round(Descriptors.MolWt(full), 1),
                         "protac_rotb": rdMolDescriptors.CalcNumRotatableBonds(full),
                         "protac_smi_path": cand_dir / "protac.smi"})
    df = pd.DataFrame(rows)
    print(f"{len(df)} PROTACs montados em {out_dir}")
    return df


# Carregue as versões com [2*] livre — NUNCA os *_capped_md.sdf:
# subcomplex_mols = {p.stem: Chem.MolFromMolFile(str(p))
#                    for p in sorted(WP2_DIR.rglob("*.sdf"))
#                    if "_capped_md" not in p.name}
# protac_df = build_protac_matrix(subcomplex_mols, serie_A["warhead_id"].tolist())
''')

md(r"""
### 3.2a Jobs do PRosettaC

`Anchor atoms` = índice, em cada `Heads`, do átomo que liga ao linker. Nos SDFs
deste pipeline ele é o **primeiro**, mas sua build pode contar de 0 ou de 1.
**Rode um candidato sozinho e confira** antes de liberar o lote.
""")

code(r'''
WARHEAD_ANCHOR_SERIAL = 1      # troque para 0 se a sua build for 0-based


def generate_prosettac_config(struct_a, chain_a, struct_b, chain_b, head_a,
                              head_b, anchor_a, anchor_b, protac_smi,
                              out_dir: Path, full_run=False):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg = out_dir / "prosetta_config.txt"
    cfg.write_text(
        f"Structures: {struct_a} {struct_b}\n"
        f"Chains: {chain_a} {chain_b}\n"
        f"Heads: {head_a} {head_b}\n"
        f"Anchor atoms: {anchor_a} {anchor_b}\n"
        f"Protac: {protac_smi}\n"
        f"Full: {full_run}\n"
        f"ClusterName: SLURM\n")     # sem isto o default é PBS, que não existe aqui
    return cfg


def emit_prosettac_jobs(protac_df, struct_e3, chain_e3, struct_pcsk9,
                        chain_pcsk9, head_e3_sdf, e3_anchor,
                        head_warhead_dir, out_root: Path = None,
                        full_run=False):
    out_root = Path(out_root or WP3_DIR / "prosettac")
    out_root.mkdir(parents=True, exist_ok=True)
    cmds = []
    for r in protac_df.itertuples():
        cfg = generate_prosettac_config(
            struct_e3, chain_e3, struct_pcsk9, chain_pcsk9, head_e3_sdf,
            Path(head_warhead_dir) / f"{r.warhead_id}_in_pcsk9.sdf",
            e3_anchor, WARHEAD_ANCHOR_SERIAL, r.protac_smi_path,
            out_root / r.candidate_id, full_run)
        cmds.append(f"cd {cfg.parent} && nohup {PROSETTAC_DIR}/run_prosettac.sh "
                    f"{cfg.parent} {cfg.name} > run_prosettac.log 2>&1 &")
    script = out_root / "launch_all_prosettac.sh"
    script.write_text("#!/usr/bin/env bash\nset -u\n" + "\n".join(cmds) + "\n")
    script.chmod(0o755)
    print(f"{len(cmds)} jobs -> {script}")
    print("RODE UM PRIMEIRO e confira o config antes de liberar o lote.")
    return script
''')

md("### 3.2b Carga líquida de cada PROTAC (para `acpype -n` e `genion`)")

code(r'''
def protac_net_charge(protac_smiles: str, warhead_id: str, wh_df=None) -> int:
    """Carga em pH 7,4 = carga do warhead (do manifesto) + carga formal do
    resto. O N de conjugação deixa de ser básico quando vira amida/sulfonamida
    no acoplamento — confirme caso a caso."""
    wh_df = wh_df if wh_df is not None else wh
    q_wh = int(wh_df.loc[wh_df["warhead_id"] == warhead_id,
                         "net_charge_ph74"].iloc[0])
    return Chem.GetFormalCharge(Chem.MolFromSmiles(protac_smiles)) + q_wh
''')

md(r"""
## 3.3 AlphaFold 3 — submissão manual

AF3 não roda local nesta máquina. As células abaixo só montam o JSON no schema
do AlphaFold Server; confira o schema atual antes de submeter.
""")

code(r'''
def build_alphafold3_job(job_name, protein_sequences: dict, ligand_smiles: str,
                         out_json: Path):
    Path(out_json).parent.mkdir(parents=True, exist_ok=True)
    job = {"name": job_name, "modelSeeds": [1],
           "sequences": ([{"proteinChain": {"sequence": s, "count": 1}}
                          for s in protein_sequences.values()]
                         + [{"ligand": {"smiles": ligand_smiles, "count": 1}}])}
    Path(out_json).write_text(json.dumps([job], indent=2))
    return out_json


def emit_af3_jobs(protac_df, e3_sequences: dict, pcsk9_sequence: str,
                  out_dir: Path = None, top_n: int = 10):
    out_dir = Path(out_dir or WP3_DIR / "af3")
    seqs = dict(e3_sequences); seqs["PCSK9"] = pcsk9_sequence
    paths = [build_alphafold3_job(r.candidate_id, seqs, r.protac_smiles,
                                  out_dir / f"{r.candidate_id}.json")
             for r in protac_df.head(top_n).itertuples()]
    print(f"{len(paths)} JSONs em {out_dir} — submeta em alphafoldserver.com")
    return paths
''')

md("## 3.4 Filtrar por métricas nativas de confiança")

code(r'''
def filter_prosettac_models(score_sc_path: Path, score_cutoff=0.0):
    rows, header = [], None
    for line in Path(score_sc_path).read_text().splitlines():
        if not line.startswith("SCORE:"):
            continue
        parts = line.split()
        if header is None:
            header = parts
            continue
        rows.append(parts)
    df = pd.DataFrame(rows, columns=header)
    num = [c for c in df.columns if c not in ("SCORE:", "description")]
    df[num] = df[num].apply(pd.to_numeric, errors="coerce")
    col = "total_score" if "total_score" in df.columns else num[0]
    return df[df[col] <= score_cutoff].sort_values(col)


def filter_af3_models(confidence_json_paths, iptm_cutoff=0.6, pae_cutoff=10.0):
    df = pd.DataFrame([{"model": Path(p).parent.name,
                        **{k: json.loads(Path(p).read_text()).get(k)
                           for k in ("iptm", "ptm", "pae_mean")}}
                       for p in confidence_json_paths])
    return df[(df["iptm"] >= iptm_cutoff) & (df["pae_mean"] <= pae_cutoff)]
''')

md("## 3.5 Checagens estruturais do complexo ternário")

code(r'''
def interface_contacts(universe_path, sel_a, sel_b, cutoff=4.5) -> int:
    import MDAnalysis as mda
    from MDAnalysis.analysis.distances import distance_array
    u = mda.Universe(str(universe_path))
    d = distance_array(u.select_atoms(sel_a).positions,
                       u.select_atoms(sel_b).positions)
    return int((d <= cutoff).sum())


def structural_checks(model_pdb: Path, e3_chain: str, pcsk9_chain: str,
                      warhead_resname: str, reference_e3_recruiter_pdb: Path):
    return {
        "warhead_pcsk9_contacts": interface_contacts(
            model_pdb, f"resname {warhead_resname}", f"chainID {pcsk9_chain}"),
        "e3_recruiter_rmsd_to_reference": rmsd_to_reference(
            model_pdb, reference_e3_recruiter_pdb, f"chainID {e3_chain}"),
        # contato direto E3-PCSK9 sem o PROTAC mediando = geometria improdutiva
        "spurious_e3_pcsk9_direct_contacts": interface_contacts(
            model_pdb, f"chainID {e3_chain}", f"chainID {pcsk9_chain}"),
    }
''')

md(r"""
## 3.6 MD (GROMACS) — três níveis, 3 réplicas cada

- nível (i) `E3_recruiter_linker` → use os `*_capped_md.sdf` do WP2
- nível (ii) `E3_recruiter_linker_warhead` → **rode antes do ternário**; é onde
  aparece interação espúria warhead–E3 sem a PCSK9
- nível (iii) `full_ternary`

AMBER ff14SB + GAFF2/AM1-BCC (acpype), TIP3P, box cúbico 13 Å, cutoff 9 Å, PME.
""")

code(r'''
def parametrize_ligand_gaff2(ligand_pdb: Path, resname: str, net_charge: int,
                             out_dir: Path):
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    cmd = [ACPYPE_EXE, "-i", str(ligand_pdb), "-b", resname,
           "-n", str(net_charge), "-a", "gaff2", "-c", "bcc"]
    print("Terminal (env mdtools):", " ".join(cmd), " (cwd =", out_dir, ")")
    print("A carga vem de protac_net_charge() — errada aqui, contamina a MD "
          "inteira em silêncio.")
    return cmd


def write_mdp_files(out_dir: Path):
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    base = (f"cutoff-scheme = Verlet\ncoulombtype = PME\n"
            f"rvdw = {RVDW_RCOULOMB_NM}\nrcoulomb = {RVDW_RCOULOMB_NM}\n")
    tc = (f"tcoupl = Nose-Hoover\ntc-grps = Protein_Ligand Water_and_ions\n"
          f"tau-t = 1.0 1.0\nref-t = {TEMP_K} {TEMP_K}\n")
    pc = (f"pcoupl = Parrinello-Rahman\npcoupltype = isotropic\ntau-p = 2.0\n"
          f"ref-p = {PRESSURE_ATM}\ncompressibility = 4.5e-5\n")

    em = "integrator = steep\nemtol = 1000.0\nemstep = 0.01\nnsteps = 5000\n" + base
    files = {
        "em_steep.mdp": em,
        "em_cg.mdp": em.replace("integrator = steep", "integrator = cg"),
        "nvt.mdp": (f"integrator = md\ndt = {DT_PS}\nnsteps = 50000\n" + base + tc
                    + f"pcoupl = no\ngen-vel = yes\ngen-temp = {TEMP_K}\n"
                      "define = -DPOSRES\n"),
        "npt.mdp": (f"integrator = md\ndt = {DT_PS}\n"
                    f"nsteps = {int(EQUIL_NPT_NS * 1000 / DT_PS)}\n"
                    + base + tc + pc + "gen-vel = no\ndefine = -DPOSRES\n"),
        "production.mdp": (f"integrator = md\ndt = {DT_PS}\n"
                           f"nsteps = {int(PRODUCTION_NS * 1000 / DT_PS)}\n"
                           + base + tc + pc
                           + f"gen-vel = no\n"
                             f"nstxout-compressed = {int(COORD_SAVE_PS / DT_PS)}\n"
                             f"define =\n"),
    }
    for name, content in files.items():
        (out_dir / name).write_text(content)
    print("mdp escritos em", out_dir)
    return {n: out_dir / n for n in files}


def gromacs_prep_commands(system_pdb: Path, out_dir: Path = None):
    cmds = [
        f"{GMX_EXE} pdb2gmx -f {system_pdb} -ff amber14sb -water tip3p "
        f"-o processed.gro -p topol.top",
        f"{GMX_EXE} editconf -f processed.gro -o boxed.gro -c -d {BOX_DIST_NM} -bt cubic",
        f"{GMX_EXE} solvate -cp boxed.gro -cs spc216.gro -p topol.top -o solvated.gro",
        f"{GMX_EXE} grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr",
        f"{GMX_EXE} genion -s ions.tpr -o neutral.gro -p topol.top "
        f"-pname NA -nname CL -neutral",
    ]
    for c in cmds:
        print(c)
    return cmds


def replicate_seeds(level: str, n_replicates=N_REPLICATES, base_seed=1000):
    return {f"{level}_rep{i+1}": base_seed + i * 97 for i in range(n_replicates)}
''')

md("## 3.7 Análise de trajetória (MDAnalysis)")

code(r'''
def load_trajectory(tpr_or_gro: Path, xtc: Path):
    import MDAnalysis as mda
    return mda.Universe(str(tpr_or_gro), str(xtc))


def traj_rmsd(u, selection, ref_frame=0):
    from MDAnalysis.analysis import rms
    R = rms.RMSD(u, select=selection, ref_frame=ref_frame); R.run()
    return pd.DataFrame(R.results.rmsd, columns=["frame", "time_ps", "rmsd_A"])


def traj_rmsf(u, selection):
    from MDAnalysis.analysis import rms
    ag = u.select_atoms(selection)
    R = rms.RMSF(ag).run()
    return pd.DataFrame({"atom": [a.name for a in ag],
                         "resid": [a.resid for a in ag],
                         "rmsf_A": R.results.rmsf})


def traj_hbonds(u, sel_donor, sel_acceptor):
    from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
    hb = HydrogenBondAnalysis(u, donors_sel=sel_donor, hydrogens_sel=None,
                              acceptors_sel=sel_acceptor)
    hb.run()
    return hb.count_by_time()


def traj_contact_frequency(u, sel_a, sel_b, cutoff=4.5):
    from MDAnalysis.analysis.distances import distance_array
    linhas = []
    for ts in u.trajectory:
        d = distance_array(u.select_atoms(sel_a).positions,
                           u.select_atoms(sel_b).positions)
        linhas.append({"time_ps": ts.time, "n_contacts": int((d <= cutoff).sum())})
    return pd.DataFrame(linhas)


def linker_torsions(u, dihedral_atom_groups):
    from MDAnalysis.analysis.dihedrals import Dihedral
    ags = [u.atoms[list(idxs)] for idxs in dihedral_atom_groups]
    dih = Dihedral(ags).run()
    return pd.DataFrame(dih.results.angles,
                        columns=[f"dih_{i}" for i in range(len(dihedral_atom_groups))])


def cluster_trajectory(u, selection, n_clusters=5):
    from sklearn.cluster import KMeans
    coords = np.array([u.select_atoms(selection).positions.flatten()
                       for _ in u.trajectory])
    return KMeans(n_clusters=n_clusters, random_state=0,
                  n_init=10).fit_predict(coords)


def sasa_lysines_gmx(tpr: Path, xtc: Path, out_dir: Path):
    """Exposição das lisinas da PCSK9 ao longo do ternário. Informação
    estrutural complementar — não é evidência de ubiquitinação por si só."""
    cmd = [GMX_EXE, "sasa", "-s", str(tpr), "-f", str(xtc),
           "-o", str(Path(out_dir) / "lys_sasa.xvg"),
           "-output", str(Path(out_dir) / "lys_sasa_residue.xvg")]
    print("Terminal:", " ".join(cmd), "(selecione o grupo das lisinas)")
    return cmd
''')

md(r"""
## 3.8 Integração final e priorização

Junta SBVS (WP1), Boltz-2 (WP2), MD binária/pré-ternária e PRosettaC/AF3 do
ternário completo. Os PROTACs com maior consistência estrutural, estabilidade
de interação e viabilidade sintética seguem para síntese.
""")

code(r'''
def integrate_prioritization(wp3_prosettac: pd.DataFrame, wp3_af3: pd.DataFrame,
                             md_stability: pd.DataFrame,
                             wp1_scores: pd.DataFrame = None,
                             wp2_boltz: pd.DataFrame = None) -> pd.DataFrame:
    merged = (wp3_prosettac.merge(wp3_af3, on="candidate_id", how="inner")
                           .merge(md_stability, on="candidate_id", how="inner"))
    for extra, chave in ((wp1_scores, "candidate_id"), (wp2_boltz, "candidate_id")):
        if extra is not None and chave in extra.columns:
            merged = merged.merge(extra, on=chave, how="left")
    ordenar = [c for c in ("structural_consistency", "interaction_stability")
               if c in merged.columns]
    return merged.sort_values(ordenar, ascending=False) if ordenar else merged
''')

md(r"""
**Milestone WP3:** conjunto priorizado e validado por MD de candidatos PROTAC
(recrutador VHL/CRBN – linker – warhead PCSK9) com estabilidade do complexo
ternário demonstrada.

---

## Checklist final

- [ ] Fase 0: `revalidate_redocking.py` rodado; VHL/CRBN passam pelo RMSD correto
- [ ] Fase 1: 180 warheads gerados; PDBQT criados
- [ ] Fase 2: sítio do 6U26 definido; redocking do `063` ≤ 2,0 Å; triagem feita
- [ ] Fase 3: linkers bifuncionais; `validar_subcomplexos()` sem falhas
- [ ] Fase 4: 1 job PRosettaC conferido antes do lote; `acpype -n` com a carga certa
- [ ] MD nível (ii) rodado antes do ternário completo
""")

# ===========================================================================
nb = {
    "cells": CELLS,
    "metadata": {
        "kernelspec": {"display_name": "mdtools", "language": "python",
                       "name": "python3"},
        "language_info": {"name": "python", "version": "3.11"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}

out = Path(__file__).resolve().parent.parent / "notebooks" / \
    "protac_wp1_wp3_pipeline_v2.ipynb"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(nb, indent=1, ensure_ascii=False))
n_code = sum(1 for c in CELLS if c["cell_type"] == "code")
print(f"{out}\n{len(CELLS)} células ({n_code} de código, "
      f"{len(CELLS) - n_code} markdown)")
