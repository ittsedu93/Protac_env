# ---------------------------------------------------------------------------
# WP1 (revalidação do portão de docking) e WP2 (seleção de linkers revisada)
# para o protac_wp1_wp3_pipeline.ipynb.
#
# Formato `# %%` — o VS Code renderiza cada bloco como célula.
#
# Estas células substituem:
#   - célula 1.6  `redocking_rmsd` / `validate_docking_protocol`
#   - célula 2.2  `has_terminal_group` / `filter_linkers`
#   - célula 2.3  `exit_vector_score`
#   - célula 2.5  `cap_free_terminus`
#   - célula 3.2  `assemble_full_protac`
# ---------------------------------------------------------------------------

# %%
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

SCRIPTS_DIR = Path.home() / "Protac_env" / "scripts"   # ajuste se necessário
sys.path.insert(0, str(SCRIPTS_DIR))

from rmsd_inplace import rmsd_inplace, buried_atom_indices       # noqa: E402
from wp2_linker_tools import (                                    # noqa: E402
    find_attachment_points, is_bifunctional, filter_linkers,
    label_attachment_points, cap_for_md, assemble_protac,
    evaluate_linker, place_conformer_at_exit, find_dummy,
    RECRUITER_ISOTOPE, WARHEAD_ISOTOPE,
)


# %% [markdown]
# # WP1 — revalidação do portão de redocking
#
# `redocking_rmsd()` (célula 1.6) usa `rdMolAlign.GetBestRMS`, que **superpõe
# as moléculas antes de medir**. Uma pose ancorada a 25 Å do sítio correto é
# reportada como RMSD 0,00. Como é esse portão que autoriza toda a triagem do
# WP1, as conclusões de VHL e CRBN precisam ser reconferidas.
#
# Segundo efeito: `GetBestRMS` **altera as coordenadas do probe**. Quem chamou e
# depois reaproveitou o objeto da pose trabalhou com ela deslocada.

# %%
def redocking_rmsd_correto(pose_path, ref_ligand_path, receptor_pdb=None,
                           core_threshold: int = 20):
    """Substitui `redocking_rmsd()` da célula 1.6.

    Mede sem superposição e com simetria. Se `receptor_pdb` for dado, devolve
    também o RMSD restrito ao núcleo enterrado — o número que decide quando o
    ligante de referência tem partes flexíveis expostas ao solvente que o
    docking não tem como reproduzir.
    """
    def _load(p):
        p = Path(p)
        if p.suffix.lower() in (".sdf", ".mol"):
            return Chem.MolFromMolFile(str(p), removeHs=True)
        if p.suffix.lower() == ".pdb":
            return Chem.MolFromPDBFile(str(p), removeHs=True)
        raise ValueError(f"converta {p.suffix} para SDF antes (obabel)")

    probe, ref = _load(pose_path), _load(ref_ligand_path)
    if probe is None or ref is None:
        raise ValueError(f"falha ao ler {pose_path} ou {ref_ligand_path}")

    out = {"rmsd_full_A": rmsd_inplace(probe, ref)}
    if receptor_pdb:
        rec = np.array([(float(l[30:38]), float(l[38:46]), float(l[46:54]))
                        for l in Path(receptor_pdb).read_text().splitlines()
                        if l.startswith("ATOM")
                        and (l[76:78].strip() or "C") != "H"])
        core = buried_atom_indices(ref, rec, threshold=core_threshold)
        if len(core) >= 3:
            out["rmsd_core_A"] = rmsd_inplace(probe, ref, atom_indices=core)
            out["n_core_atoms"] = len(core)
    return out


# %% [markdown]
# ### Reconferir VHL e CRBN
#
# Para cada alvo do WP1, aponte o receptor, o ligante co-cristalizado e as poses
# de redocking que já estão no disco. Nada é dockado de novo — só remedido.

# %%
WP1_ALVOS = {
    # "VHL_6GFZ": dict(
    #     receptor=WORK_DIR / "prep/VHL_6GFZ/6GFZ_receptor.pdb",
    #     ref_ligand=WORK_DIR / "prep/VHL_6GFZ/6GFZ_ref_ligand.sdf",
    #     poses=str(WORK_DIR / "prep/VHL_6GFZ/redock_pose.pdbqt"),
    # ),
    # "CRBN_4CI1": dict(...),
}

for nome, cfg in WP1_ALVOS.items():
    print(f"python {SCRIPTS_DIR}/revalidate_redocking.py \\")
    print(f"    --receptor-pdb {cfg['receptor']} \\")
    print(f"    --ref-ligand   {cfg['ref_ligand']} \\")
    print(f"    --poses        '{cfg['poses']}' \\")
    print(f"    --label        {nome} \\")
    print(f"    --out          {WORK_DIR}/wp1_revalidation_{nome}.csv\n")

# %%
# Consolidar os CSVs depois de rodar:
def consolidar_revalidacao(work_dir: Path = None):
    work_dir = work_dir or WORK_DIR
    csvs = sorted(Path(work_dir).glob("wp1_revalidation_*.csv"))
    if not csvs:
        print("nenhum CSV de revalidação encontrado")
        return None
    df = pd.concat([pd.read_csv(c) for c in csvs], ignore_index=True)
    resumo = df.groupby("label").agg(
        poses=("pose", "count"),
        passaram=("passa_corte", "sum"),
        melhor_rmsd=("rmsd_core_A", "min"),
    )
    print(resumo)
    reprovados = resumo[resumo["passaram"] == 0]
    if len(reprovados):
        print(f"\nATENÇÃO: {len(reprovados)} alvo(s) sem nenhuma pose válida:")
        print(f"  {list(reprovados.index)}")
        print("  A triagem do WP1 nesses alvos não tem protocolo validado.")
    return df


# %% [markdown]
# # WP2 — seleção de linkers revisada
#
# Três defeitos corrigidos, na ordem de gravidade:
#
# **1. O cap destruía a montagem do WP3, em silêncio.** `cap_free_terminus()`
# troca todo `[#0]` por metil. Depois disso `assemble_full_protac()` procura
# `[#0]`, não acha, e o RDKit devolve a molécula **inalterada** numa tupla de um
# elemento — sem exceção, e `SanitizeMol` passa. O lote inteiro vira
# recrutador-linker sem warhead nenhum, e nada avisa. Agora as pontas são
# rotuladas (`[1*]` recrutador, `[2*]` warhead) e a montagem é verificada.
#
# **2. `filter_linkers()` aceitava linker monofuncional.** `has_terminal_group()`
# devolve True com UM match. Um linker de PROTAC é bifuncional por definição.
# Além disso `[NH2]` casa amina não terminal (`CC(N)CC` casa), então "grupo
# terminal" não era terminal.
#
# **3. `exit_vector_score()` media confôrmero solto no espaço.**
# `EmbedMultipleConfs` gera coordenadas num referencial arbitrário; a função
# calculava clash contra o receptor sem nunca ter posicionado o linker no exit
# vector. As distâncias eram entre nuvens de pontos sem relação espacial.

# %% [markdown]
# ## 2.1 Configuração

# %%
ENAMINE_LINKER_DIR = WORK_DIR / "enamine_protac_linkers"

LINKER_FILTERS = dict(atom_count_range=(5, 40), rotb_max=15, min_span=3)
N_CONFORMERS = 50
N_SPINS = 12                # rotações amostradas em torno do eixo de saída
CLASH_DISTANCE_A = 2.5      # < 2,5 Å entre pesados = clash

WP2_DIR = WORK_DIR / "wp2_subcomplexes"
WP2_DIR.mkdir(parents=True, exist_ok=True)


# %% [markdown]
# ## 2.2 Filtrar linkers (agora exigindo bifuncionalidade)

# %%
# linkers = load_library(ENAMINE_LINKER_DIR)          # do WP1, célula 1.4
# linkers_ok = filter_linkers(linkers, **LINKER_FILTERS)
#
# Cada molécula aprovada leva as props attach_1_idx/attach_2_idx com os dois
# pontos de conjugação já identificados.


# %% [markdown]
# ## 2.3 Avaliação geométrica a partir do exit vector do recrutador
#
# `exit_point` e `exit_direction` vêm do WP1 (átomo do exit vector do recrutador
# na pose validada, e a direção que aponta para o solvente).

# %%
def receptor_coords_from_pdb(receptor_pdb: Path) -> np.ndarray:
    return np.array([(float(l[30:38]), float(l[38:46]), float(l[46:54]))
                     for l in Path(receptor_pdb).read_text().splitlines()
                     if l.startswith("ATOM") and (l[76:78].strip() or "C") != "H"])


def avaliar_linkers(linkers, exit_point, exit_direction, receptor_pdb,
                    n_confs=N_CONFORMERS, n_spins=N_SPINS):
    """Para cada linker: posiciona no exit vector e devolve a melhor colocação."""
    rec = receptor_coords_from_pdb(receptor_pdb)
    linhas = []
    for m in linkers:
        attach = int(m.GetProp("attach_1_idx"))
        molh, res = evaluate_linker(m, attach, exit_point, exit_direction,
                                    rec, n_confs=n_confs, n_spins=n_spins)
        if not res:
            print(f"  [sem confôrmero] {Chem.MolToSmiles(m)}")
            continue
        melhor = res[0]
        sem_clash = [r for r in res if r["n_clashes"] == 0]
        linhas.append({
            "smiles": Chem.MolToSmiles(m),
            "n_heavy": m.GetNumHeavyAtoms(),
            "extension_A": melhor["extension_A"],
            "cos_to_exit": melhor["cos_angle_to_exit"],
            "n_clashes_melhor": melhor["n_clashes"],
            "min_dist_A": melhor["min_dist_to_receptor_A"],
            "confs_sem_clash": len(sem_clash),
            "frac_sem_clash": round(len(sem_clash) / len(res), 3),
        })
    df = pd.DataFrame(linhas)
    if len(df):
        df = df.sort_values(["n_clashes_melhor", "cos_to_exit"],
                            ascending=[True, False])
    return df


# %% [markdown]
# `frac_sem_clash` é o critério mais informativo: um linker que só encaixa em
# 1 de 50 confôrmeros está geometricamente forçado, mesmo que esse único
# confôrmero pontue bem. Prefira os que acomodam o exit vector em boa parte do
# seu espaço conformacional.

# %% [markdown]
# ## 2.4 Rotular as duas pontas e montar o sub-complexo
#
# A partir daqui cada ponta é rastreável. **Nunca capeie o objeto que segue para
# o WP3** — capeie uma cópia.

# %%
def preparar_subcomplexo(linker, recruiter_mol, out_dir: Path, sc_id: str):
    """Rotula as pontas, conjuga o recrutador em [1*] e guarda as duas versões:
    a que vai para o WP3 (com [2*] livre) e a capeada para a MD de nível (i)."""
    out_dir = Path(out_dir) / sc_id
    out_dir.mkdir(parents=True, exist_ok=True)

    a1, a2 = int(linker.GetProp("attach_1_idx")), int(linker.GetProp("attach_2_idx"))
    rotulado = label_attachment_points(linker, a1, a2)

    # recrutador entra no [1*]; recruiter_mol precisa ter o átomo de conjugação
    # no índice 0, mesma convenção dos warheads
    sub = AllChem.ReplaceSubstructs(
        rotulado, Chem.MolFromSmarts(f"[{RECRUITER_ISOTOPE}#0]"),
        recruiter_mol, replacementConnectionPoint=0)[0]
    Chem.SanitizeMol(sub)

    if find_dummy(sub, WARHEAD_ISOTOPE) is None:
        raise ValueError(f"{sc_id}: o [2*] sumiu na conjugação do recrutador")

    # versão para o WP3 — [2*] intacto
    with Chem.SDWriter(str(out_dir / f"{sc_id}.sdf")) as w:
        w.write(sub)
    # versão para a MD de nível (i) — cópia capeada
    with Chem.SDWriter(str(out_dir / f"{sc_id}_capped_md.sdf")) as w:
        w.write(cap_for_md(sub))

    return {"subcomplex_id": sc_id, "smiles": Chem.MolToSmiles(sub),
            "path_wp3": out_dir / f"{sc_id}.sdf",
            "path_md": out_dir / f"{sc_id}_capped_md.sdf"}


# %% [markdown]
# ## 2.5 Checagem antes de entregar ao WP3
#
# Rode isto sobre a pasta de sub-complexos. É o que impede o lote de PROTACs sem
# warhead de existir.

# %%
def validar_subcomplexos(wp2_dir: Path = None):
    wp2_dir = Path(wp2_dir or WP2_DIR)
    arquivos = [p for p in wp2_dir.rglob("*.sdf") if "_capped_md" not in p.name]
    if not arquivos:
        print(f"nenhum SDF de sub-complexo em {wp2_dir}")
        return None

    linhas = []
    for p in arquivos:
        m = Chem.MolFromMolFile(str(p), removeHs=False, sanitize=False)
        if m is None:
            linhas.append({"arquivo": p.name, "ok": False,
                           "motivo": "ilegível"})
            continue
        tem_wh = find_dummy(m, WARHEAD_ISOTOPE) is not None
        tem_rec = find_dummy(m, RECRUITER_ISOTOPE) is not None
        linhas.append({
            "arquivo": p.name, "ok": tem_wh and not tem_rec,
            "motivo": ("pronto" if tem_wh and not tem_rec else
                       "[2*] ausente — foi capeado?" if not tem_wh else
                       "[1*] sobrou — recrutador não foi conjugado"),
        })
    df = pd.DataFrame(linhas)
    print(f"{df['ok'].sum()}/{len(df)} sub-complexos prontos para o WP3")
    if not df["ok"].all():
        print(df[~df["ok"]].to_string(index=False))
    return df


# %% [markdown]
# ## 2.6 Boltz-2 — só sobre os que passaram
#
# `filter_boltz_models()` (célula 2.6) cai em `df.columns[-1]` quando não acha
# `confidence_score`, ou seja, pode filtrar por uma coluna qualquer sem avisar.
# Confira o nome real da métrica na saída da sua versão do Boltz antes de usar.

# %%
def filter_boltz_models_estrito(confidence_json_paths, confidence_cutoff=0.6,
                                col="confidence_score"):
    import json
    linhas = [{"model": Path(p).parent.name, **json.loads(Path(p).read_text())}
              for p in confidence_json_paths]
    df = pd.DataFrame(linhas)
    if col not in df.columns:
        raise KeyError(
            f"coluna '{col}' ausente. Disponíveis: {list(df.columns)}. "
            f"Escolha explicitamente em vez de aceitar um fallback silencioso.")
    sobreviventes = df[df[col] >= confidence_cutoff]
    print(f"{len(sobreviventes)}/{len(df)} modelos acima de {confidence_cutoff}")
    return df, sobreviventes
