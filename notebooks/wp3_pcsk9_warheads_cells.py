# ---------------------------------------------------------------------------
# Células de integração dos warheads PCSK9 (feniletilamina) na Seção 3 do
# protac_wp1_wp3_pipeline.ipynb.
#
# Arquivo em formato `# %%` — o VS Code renderiza cada bloco como célula, então
# dá para rodar aqui mesmo ou copiar bloco a bloco para o notebook.
#
# Pré-requisito: rodar antes, uma vez,
#     python scripts/generate_pcsk9_warheads.py --outdir ~/PCSK9_warheads
# ---------------------------------------------------------------------------

# %% [markdown]
# ## 3.1 Configuração — warheads PCSK9 (substitui os ligantes GMF)

# %%
from pathlib import Path
import json
import subprocess

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# --- saída do generate_pcsk9_warheads.py ---
WARHEADS_DIR = Path.home() / "PCSK9_warheads"
WARHEADS_CSV = WARHEADS_DIR / "pcsk9_warheads.csv"
WARHEADS_SDF_DIR = WARHEADS_DIR / "sdf"          # atomo 0 = N de conjugação
WARHEADS_PDBQT_DIR = WARHEADS_DIR / "pdbqt"      # criado por prepare_pdbqt.sh

# --- alvo: PCSK9 ---
# 2P4E = PCSK9 madura (pró-domínio + catalítico); 3BPS = PCSK9:LDLR-EGF(A).
# Escolha conforme o sítio que os warheads devem ocupar e confira as cadeias.
PCSK9_PDB_IDS = ["2P4E", "3BPS"]
PCSK9_CHAIN = "A"                                 # TODO: conferir no PDB baixado

# --- sub-complexos recrutador-linker validados no WP2 ---
WP2_SUBCOMPLEXES_DIR = WORK_DIR / "wp2_subcomplexes"

# --- parâmetros MD (inalterados da metodologia) ---
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


# %% [markdown]
# ### 3.1a Carregar o manifesto e priorizar a série de warheads
#
# O gerador NÃO descarta nada — ele marca. A priorização fica aqui, explícita,
# para poder ser justificada na tese.

# %%
wh = pd.read_csv(WARHEADS_CSV)
print(f"{len(wh)} warheads enumerados")
print(wh.groupby(["r3"]).size().rename("n"))

# Critérios de priorização (ajuste e documente):
#  1. dentro da faixa típica de warhead (filter_flags == 'ok')
#  2. série de N primário (R3 = H) => ponto de conjugação limpo para o linker
#  3. diversidade: pelo menos um representante de cada R1 e cada R2
prioritized = wh[(wh["filter_flags"] == "ok")].copy()

# Série A — conjugação no N livre (recomendada como principal)
serie_A = prioritized[prioritized["r3"] == "H"]

# Série B — R3 farmacofórico mantido (metilsulfonil / tetrazol); o linker entra
# no H restante e o N vira terciário. Só siga com estes se a inspeção do pose em
# PCSK9 mostrar que o N-H não faz ligação de hidrogênio essencial.
serie_B = prioritized[prioritized["r3"].isin(["Ms", "Tet", "CH2Tet"])]

print(f"\nSérie A (R3=H, N primário): {len(serie_A)}")
print(f"Série B (R3 mantido):       {len(serie_B)}")
display(serie_A[["warhead_id", "r1", "r2", "r3", "mw", "clogp",
                 "net_charge_ph74"]].head(10))


# %% [markdown]
# ## 3.0 (NOVO — rode antes do 3.2) Ancorar os warheads na PCSK9
#
# **Por que esta etapa não existia:** a metodologia original partia dos ligantes
# GMF, que já vinham com pose conhecida na PCSK9. Os warheads da Figura 3 são
# moléculas novas — o PRosettaC precisa de `HeadB.sdf` **já posicionado** dentro
# da estrutura da PCSK9, então é preciso produzir essa pose primeiro.
#
# Reaproveite as funções do WP1 (`independent_docking_runs`, `score_summary`,
# `validate_docking_protocol`) apontando o receptor para a PCSK9.

# %%
# 1) baixar e preparar a PCSK9 (mesmo caminho do WP1)
pcsk9_paths = {pid: fetch_pdb(pid, WORK_DIR / "pdb" / "PCSK9") for pid in PCSK9_PDB_IDS}
print(pcsk9_paths)

# 2) definir o sítio: centro/box do bolsão que os warheads devem ocupar.
#    Sem ligante co-cristalizado não há redocking para validar — declare o sítio
#    a partir da literatura/análise de cavidade e registre essa escolha.
PCSK9_SITE_CENTER = (0.0, 0.0, 0.0)   # TODO: preencher (ChimeraX/fpocket)
PCSK9_BOX_SIZE = (24.0, 24.0, 24.0)

# 3) docking em dois estágios, igual ao WP1
# warhead_pdbqts = {r.warhead_id: WARHEADS_PDBQT_DIR / f"{r.warhead_id}.pdbqt"
#                   for r in serie_A.itertuples()}
# r1_pcsk9, r2_pcsk9 = staged_screening(
#     warhead_pdbqts, receptor_pdbqt=pcsk9_receptor_pdbqt,
#     center=PCSK9_SITE_CENTER, box_size=PCSK9_BOX_SIZE,
#     out_dir=WP3_DIR / "warhead_docking", round1_cutoff_score=-7.0)
#
# 4) converter a melhor pose de cada warhead priorizado em SDF — este é o
#    HeadB.sdf do PRosettaC (pose na PCSK9, não o confôrmero isolado):
# subprocess.run([OBABEL_EXE, str(best_pose_pdbqt), "-O", str(head_b_sdf)], check=True)


# %% [markdown]
# ## 3.2 Montar o PROTAC completo (recrutador–linker + warhead PCSK9)
#
# `assemble_full_protac()` do notebook funciona sem alteração com estes
# warheads: o gerador renumera cada molécula para que **o átomo 0 seja o
# nitrogênio de conjugação**, que é exatamente o que
# `ReplaceSubstructs(..., replacementConnectionPoint=0)` espera.

# %%
def load_warhead_for_assembly(warhead_id: str, sdf_dir: Path = WARHEADS_SDF_DIR):
    """Lê o warhead 3D com o átomo de conjugação já no índice 0."""
    mol = Chem.MolFromMolFile(str(Path(sdf_dir) / f"{warhead_id}.sdf"))
    if mol is None:
        raise FileNotFoundError(f"não consegui ler o SDF de {warhead_id}")
    assert mol.GetAtomWithIdx(0).GetSymbol() == "N", \
        f"{warhead_id}: átomo 0 deveria ser o N de conjugação"
    return mol


def build_protac_matrix(subcomplex_mols: dict, warhead_ids: list, out_dir: Path):
    """Produto cartesiano sub-complexo x warhead -> SMILES do PROTAC completo.

    subcomplex_mols: {subcomplex_id: mol do recrutador-linker com cap dummy [*]}
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for sc_id, sc_mol in subcomplex_mols.items():
        for wid in warhead_ids:
            wh_mol = load_warhead_for_assembly(wid)
            try:
                full = assemble_full_protac(sc_mol, wh_mol)
            except Exception as exc:                      # valência, cap ausente...
                print(f"  [falha] {sc_id} x {wid}: {exc}")
                continue
            cand_id = f"{sc_id}__{wid}"
            smi = Chem.MolToSmiles(full)
            cand_dir = out_dir / cand_id
            cand_dir.mkdir(parents=True, exist_ok=True)
            (cand_dir / "protac.smi").write_text(f"{smi}\t{cand_id}\n")
            rows.append({
                "candidate_id": cand_id,
                "subcomplex_id": sc_id,
                "warhead_id": wid,
                "protac_smiles": smi,
                "protac_mw": Descriptors.MolWt(full),
                "protac_rotb": rdMolDescriptors.CalcNumRotatableBonds(full),
                "protac_smi_path": cand_dir / "protac.smi",
            })
    df = pd.DataFrame(rows)
    print(f"{len(df)} PROTACs montados em {out_dir}")
    return df


# Exemplo de uso (preencha subcomplex_mols com as saídas reais do WP2):
# subcomplex_mols = {p.stem: Chem.MolFromMolFile(str(p))
#                    for p in sorted(WP2_SUBCOMPLEXES_DIR.glob("*.sdf"))}
# protac_df = build_protac_matrix(subcomplex_mols,
#                                 serie_A["warhead_id"].tolist(),
#                                 WP3_DIR / "protacs")


# %% [markdown]
# ### 3.2a Gerar um `prosetta_config.txt` por candidato
#
# `Anchor atoms` = índice, em cada arquivo `Heads`, do átomo que liga ao linker.
# Nos SDFs gerados esse átomo está na **primeira posição**. O PRosettaC da sua
# instalação pode contar a partir de 0 ou de 1 — confira num caso isolado antes
# de disparar o lote inteiro.

# %%
WARHEAD_ANCHOR_SERIAL = 1     # 1-based; troque para 0 se a sua build usar 0-based

def emit_prosettac_jobs(protac_df: pd.DataFrame,
                        struct_e3: Path, chain_e3: str,
                        struct_pcsk9: Path, chain_pcsk9: str,
                        head_e3_sdf: Path, e3_anchor: int,
                        head_warhead_sdf_dir: Path,
                        out_root: Path, full_run: bool = False):
    cmds = []
    for r in protac_df.itertuples():
        cfg = generate_prosettac_config(
            struct_a=struct_e3, chain_a=chain_e3,
            struct_b=struct_pcsk9, chain_b=chain_pcsk9,
            head_a=head_e3_sdf,
            head_b=Path(head_warhead_sdf_dir) / f"{r.warhead_id}_in_pcsk9.sdf",
            anchor_a=e3_anchor, anchor_b=WARHEAD_ANCHOR_SERIAL,
            protac_smi=r.protac_smi_path,
            out_dir=out_root / r.candidate_id,
            full_run=full_run,
        )
        cmds.append(
            f"cd {cfg.parent} && nohup {PROSETTAC_DIR}/run_prosettac.sh "
            f"{cfg.parent} {cfg.name} > run_prosettac.log 2>&1 &"
        )
    script = out_root / "launch_all_prosettac.sh"
    script.write_text("#!/usr/bin/env bash\nset -u\n" + "\n".join(cmds) + "\n")
    script.chmod(0o755)
    print(f"{len(cmds)} jobs. Dispare no terminal da workstation:\n  bash {script}")
    print("ATENÇÃO: isso sobe todos de uma vez. Rode 1 candidato primeiro "
          "para conferir o config, depois libere o lote (ou use SLURM).")
    return script


# %% [markdown]
# ### 3.2b Carga líquida de cada PROTAC (para `acpype -n` e `genion`)

# %%
def protac_net_charge(protac_smiles: str, warhead_id: str, wh_df: pd.DataFrame) -> int:
    """Carga do PROTAC em pH 7,4 = carga do warhead (já estimada no manifesto)
    + carga do recrutador-linker. O N de conjugação deixa de ser básico quando
    vira amida/sulfonamida no acoplamento — confirme caso a caso."""
    q_wh = int(wh_df.loc[wh_df["warhead_id"] == warhead_id, "net_charge_ph74"].iloc[0])
    q_formal = Chem.GetFormalCharge(Chem.MolFromSmiles(protac_smiles))
    return q_formal + q_wh


# %% [markdown]
# ## 3.3 AlphaFold 3 — um job por candidato priorizado

# %%
def emit_af3_jobs(protac_df: pd.DataFrame, e3_sequences: dict,
                  pcsk9_sequence: str, out_dir: Path, top_n: int = 10):
    out_dir.mkdir(parents=True, exist_ok=True)
    seqs = dict(e3_sequences)
    seqs["PCSK9"] = pcsk9_sequence
    paths = []
    for r in protac_df.head(top_n).itertuples():
        paths.append(build_alphafold3_job(
            job_name=r.candidate_id,
            protein_sequences=seqs,
            ligand_smiles=r.protac_smiles,
            out_json=out_dir / f"{r.candidate_id}.json",
        ))
    print(f"{len(paths)} JSONs em {out_dir} — submeta em https://alphafoldserver.com")
    return paths
