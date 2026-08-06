#!/usr/bin/env python3

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
import streamlit.components.v1 as components
import yaml

try:
    import py3Dmol
except ImportError:
    py3Dmol = None

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Draw
except ImportError:
    Chem = None
    Descriptors = None
    Draw = None


# ============================================================
# PATHS
# ============================================================

ROOT = Path(__file__).resolve().parent

CONFIG = ROOT / "config" / "pipeline.yaml"
RUNNER = ROOT / "run_pipeline.py"
RUNS = ROOT / "runs"

PHASE5 = (
    ROOT
    / "data/processed/phase5/pilot/MCULE-2304697456/results"
)

PHASE6 = (
    ROOT
    / "data/processed/phase6/comparative_validation"
)

SHORTLIST = PHASE6 / "phase6_candidate_shortlist.csv"

CANDIDATES = [
    "MCULE-2304697456",
    "MCULE-1302802546",
    "MCULE-5294193665",
    "MCULE-4176465450",
]

REFERENCE = "MCULE-2304697456"


# ============================================================
# PAGE
# ============================================================

st.set_page_config(
    page_title="DockForgeAI",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
<style>
.block-container {
    padding-top: 1.3rem;
    padding-bottom: 3rem;
}

[data-testid="stMetric"] {
    border: 1px solid rgba(120,120,120,0.20);
    border-radius: 12px;
    padding: 14px;
    background: rgba(120,120,120,0.06);
}

.status-card {
    border: 1px solid rgba(40,150,90,0.35);
    background: rgba(40,150,90,0.10);
    border-radius: 10px;
    padding: 12px 16px;
}
</style>
""",
    unsafe_allow_html=True,
)


# ============================================================
# HELPERS
# ============================================================

@st.cache_data
def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()

    try:
        return pd.read_csv(path)
    except Exception:
        return pd.DataFrame()


@st.cache_data
def read_text(path: Path) -> str:
    if not path.exists():
        return ""

    return path.read_text(
        encoding="utf-8",
        errors="replace",
    )


@st.cache_data
def load_xvg(path: Path) -> pd.DataFrame:
    rows = []

    if not path.exists():
        return pd.DataFrame(columns=["time", "value"])

    for line in read_text(path).splitlines():
        line = line.strip()

        if not line or line.startswith(("#", "@")):
            continue

        fields = line.split()

        if len(fields) < 2:
            continue

        try:
            rows.append(
                [float(fields[0]), float(fields[1])]
            )
        except ValueError:
            continue

    return pd.DataFrame(rows, columns=["time", "value"])


def results_dir(candidate: str) -> Path:
    if candidate == REFERENCE:
        return PHASE5

    return PHASE6 / candidate / "results"


def summary_path(candidate: str) -> Path:
    if candidate == REFERENCE:
        return results_dir(candidate) / "phase5_summary.md"

    return results_dir(candidate) / "phase6_candidate_summary.md"


def result_file(candidate: str, metric: str) -> Path:
    reference_files = {
        "ligand_rmsd": "ligand_rmsd_5ns.xvg",
        "chain_a": "chain_A_backbone_rmsd_5ns.xvg",
        "chain_b": "chain_B_backbone_rmsd_5ns.xvg",
        "hbonds": "protein_ligand_hbonds_5ns.xvg",
        "contacts": "protein_ligand_contacts_5ns.xvg",
        "distance": "protein_ligand_min_distance_5ns.xvg",
    }

    pilot_files = {
        "ligand_rmsd": "ligand_rmsd_1ns.xvg",
        "chain_a": "chain_A_backbone_rmsd_1ns.xvg",
        "chain_b": "chain_B_backbone_rmsd_1ns.xvg",
        "hbonds": "protein_ligand_hbonds_1ns.xvg",
        "contacts": "protein_ligand_contacts_1ns.xvg",
        "distance": "protein_ligand_min_distance_1ns.xvg",
    }

    filenames = (
        reference_files
        if candidate == REFERENCE
        else pilot_files
    )

    return results_dir(candidate) / filenames[metric]


def structure_path(
    candidate: str,
    structure_type: str,
) -> Path:
    folder = results_dir(candidate)

    if candidate == REFERENCE:
        filename = (
            "final_5ns_complex.pdb"
            if structure_type == "Final"
            else "representative_4p5ns_complex.pdb"
        )
    else:
        filename = (
            "final_1ns_complex.pdb"
            if structure_type == "Final"
            else "representative_0p75ns_complex.pdb"
        )

    return folder / filename


def to_ns(
    dataframe: pd.DataFrame,
    metric: str,
) -> pd.DataFrame:
    data = dataframe.copy()

    if data.empty:
        return data

    if metric in {"ligand_rmsd", "chain_a", "chain_b"}:
        data["time_ns"] = data["time"]
    else:
        data["time_ns"] = data["time"] / 1000.0

    return data


def last_half(dataframe: pd.DataFrame) -> pd.DataFrame:
    if dataframe.empty:
        return dataframe

    cutoff = dataframe["time"].max() / 2

    return dataframe[dataframe["time"] >= cutoff]


def candidate_metrics(candidate: str) -> dict[str, float]:
    ligand = load_xvg(
        result_file(candidate, "ligand_rmsd")
    )
    hbonds = load_xvg(
        result_file(candidate, "hbonds")
    )
    contacts = load_xvg(
        result_file(candidate, "contacts")
    )
    distance = load_xvg(
        result_file(candidate, "distance")
    )

    ligand_last = last_half(ligand)

    return {
        "ligand_rmsd": (
            float(ligand_last["value"].mean())
            if not ligand_last.empty
            else np.nan
        ),
        "ligand_rmsd_max": (
            float(ligand["value"].max())
            if not ligand.empty
            else np.nan
        ),
        "mean_hbonds": (
            float(hbonds["value"].mean())
            if not hbonds.empty
            else np.nan
        ),
        "hbond_occupancy": (
            float((hbonds["value"] >= 1).mean() * 100)
            if not hbonds.empty
            else np.nan
        ),
        "mean_contacts": (
            float(contacts["value"].mean())
            if not contacts.empty
            else np.nan
        ),
        "contact_occupancy": (
            float((contacts["value"] >= 1).mean() * 100)
            if not contacts.empty
            else np.nan
        ),
        "minimum_distance": (
            float(distance["value"].mean())
            if not distance.empty
            else np.nan
        ),
    }


@st.cache_data
def comparison_table() -> pd.DataFrame:
    shortlist = read_csv(SHORTLIST)
    rows = []

    for candidate in CANDIDATES:
        row = {
            "compound_id": candidate,
            **candidate_metrics(candidate),
        }

        if (
            not shortlist.empty
            and "compound_id" in shortlist.columns
        ):
            match = shortlist[
                shortlist["compound_id"].astype(str)
                == candidate
            ]

            if not match.empty:
                for column in [
                    "consensus_rank",
                    "predicted_pki",
                    "predicted_ki_nm",
                    "vina_score_kcal_mol",
                    "weighted_consensus_score",
                    "pains_flag",
                    "brenk_flag",
                    "md_status",
                ]:
                    if column in match.columns:
                        row[column] = match.iloc[0][column]

        rows.append(row)

    return pd.DataFrame(rows)


def shell(command: list[str]) -> tuple[int, str]:
    try:
        result = subprocess.run(
            command,
            cwd=ROOT,
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
    except Exception as error:
        return 1, str(error)

    return (
        result.returncode,
        (result.stdout or "") + (result.stderr or ""),
    )


def render_3d(path: Path) -> None:
    if py3Dmol is None:
        st.warning("py3Dmol is not installed.")
        return

    if not path.exists():
        st.warning(f"Structure not found: {path}")
        return

    viewer = py3Dmol.view(
        width=950,
        height=620,
    )

    viewer.addModel(read_text(path), "pdb")

    viewer.setStyle(
        {"protein": True},
        {
            "cartoon": {
                "color": "spectrum",
                "opacity": 0.85,
            }
        },
    )

    for residue in ["LIG", "MOL"]:
        viewer.setStyle(
            {"resn": residue},
            {
                "stick": {"radius": 0.25},
                "sphere": {"scale": 0.28},
            },
        )

    viewer.zoomTo()
    viewer.setBackgroundColor("white")

    components.html(
        viewer._make_html(),
        height=650,
        scrolling=False,
    )


def stage_names() -> list[str]:
    if not CONFIG.exists():
        return []

    config = yaml.safe_load(read_text(CONFIG)) or {}

    return [
        str(stage["name"])
        for stage in config.get("stages", [])
        if isinstance(stage, dict) and stage.get("name")
    ]


def execute_pipeline(
    dry_run: bool,
    resume: bool,
    from_stage: str | None,
) -> None:
    command = [
        sys.executable,
        str(RUNNER),
        "--config",
        str(CONFIG),
    ]

    if dry_run:
        command.append("--dry-run")

    if resume:
        command.append("--resume")

    if from_stage:
        command.extend(["--from-stage", from_stage])

    st.code(" ".join(command), language="bash")

    output = st.empty()
    status = st.empty()
    collected = []

    process = subprocess.Popen(
        command,
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )

    assert process.stdout is not None

    status.info("Pipeline running...")

    for line in process.stdout:
        collected.append(line)
        output.code(
            "".join(collected[-200:]),
            language="text",
        )

    code = process.wait()

    if code == 0:
        status.success("Pipeline completed successfully.")
    else:
        status.error(f"Pipeline failed with code {code}.")


# ============================================================
# SIDEBAR
# ============================================================

st.sidebar.title("🧬 DockForgeAI")

st.sidebar.markdown(
    """
**Target:** HIV-1 protease  
**Mode:** Prospective lead prioritization  
**Current release:** Portfolio v1.0
"""
)

active_candidate = st.sidebar.selectbox(
    "Active candidate",
    CANDIDATES,
)

st.sidebar.divider()

st.sidebar.markdown(
    """
### Evidence layers

- Machine-learning prediction
- Applicability domain
- Drug-likeness assessment
- PAINS/Brenk alerts
- AutoDock Vina docking
- Interaction profiling
- Consensus ranking
- Molecular dynamics
"""
)

st.sidebar.warning(
    "These simulations support computational prioritization. "
    "They do not prove experimental inhibition."
)


# ============================================================
# HEADER
# ============================================================

st.title("🧬 DockForgeAI")

st.caption(
    "AI-assisted computational drug discovery, molecular docking "
    "and molecular-dynamics validation platform"
)


(
    dashboard_tab,
    workflow_tab,
    candidate_tab,
    md_tab,
    viewer_tab,
    pipeline_tab,
    logs_tab,
    repository_tab,
) = st.tabs(
    [
        "Dashboard",
        "Workflow",
        "Candidate Explorer",
        "MD Comparison",
        "3D Viewer",
        "Run Pipeline",
        "Runs & Logs",
        "Repository",
    ]
)


# ============================================================
# DASHBOARD
# ============================================================

with dashboard_tab:
    shortlist = read_csv(SHORTLIST)
    comparison = comparison_table()

    st.header("Project dashboard")

    columns = st.columns(5)

    columns[0].metric(
        "Shortlisted candidates",
        len(shortlist),
    )

    completed = 0

    if (
        not shortlist.empty
        and "md_status" in shortlist.columns
    ):
        completed = int(
            shortlist["md_status"]
            .astype(str)
            .str.contains("complete", case=False)
            .sum()
        )

    columns[1].metric("MD completed", completed)
    columns[2].metric("Reference simulation", "5 ns")
    columns[3].metric("Comparative pilots", "3 × 1 ns")

    git_code, git_text = shell(
        ["git", "status", "--porcelain"]
    )

    columns[4].metric(
        "Repository",
        (
            "Clean"
            if git_code == 0 and not git_text.strip()
            else "Changes"
        ),
    )

    left, right = st.columns(2)

    with left:
        st.subheader("Consensus ranking")

        if (
            not shortlist.empty
            and "weighted_consensus_score"
            in shortlist.columns
        ):
            plot = shortlist.sort_values(
                "weighted_consensus_score"
            ).copy()

            plot["label"] = (
                plot["compound_id"]
                .astype(str)
                .str.replace("MCULE-", "", regex=False)
            )

            figure = px.bar(
                plot,
                x="weighted_consensus_score",
                y="label",
                orientation="h",
                text="weighted_consensus_score",
                labels={
                    "weighted_consensus_score": (
                        "Consensus score"
                    ),
                    "label": "Candidate",
                },
            )

            figure.update_traces(
                texttemplate="%{text:.3f}",
                textposition="outside",
            )

            figure.update_layout(height=410)

            st.plotly_chart(
                figure,
                width="stretch",
            )

    with right:
        st.subheader("Predicted activity vs docking")

        required = {
            "compound_id",
            "predicted_pki",
            "vina_score_kcal_mol",
        }

        if (
            not shortlist.empty
            and required.issubset(shortlist.columns)
        ):
            plot = shortlist.copy()

            plot["label"] = (
                plot["compound_id"]
                .astype(str)
                .str.replace("MCULE-", "", regex=False)
            )

            figure = px.scatter(
                plot,
                x="vina_score_kcal_mol",
                y="predicted_pki",
                size="weighted_consensus_score",
                text="label",
                hover_name="compound_id",
                labels={
                    "vina_score_kcal_mol": (
                        "Vina score (kcal/mol)"
                    ),
                    "predicted_pki": "Predicted pKi",
                },
            )

            figure.update_traces(
                textposition="top center"
            )

            figure.update_layout(height=410)

            st.plotly_chart(
                figure,
                width="stretch",
            )

    st.subheader("Comparative MD summary")

    st.dataframe(
        comparison,
        width="stretch",
        hide_index=True,
    )


# ============================================================
# WORKFLOW
# ============================================================

with workflow_tab:
    st.header("DockForgeAI workflow")

    graph = """
digraph pipeline {
    rankdir=LR;
    graph [bgcolor="transparent"];
    node [
        shape=box,
        style="rounded,filled",
        fillcolor="#EAF2FF",
        color="#4F73A5",
        fontname="Helvetica"
    ];
    edge [color="#6B7280", penwidth=1.5];

    A [label="Compound Library"];
    B [label="Structure Standardization"];
    C [label="ML pKi / Ki Prediction"];
    D [label="Applicability Domain"];
    E [label="PAINS / Brenk / ADMET"];
    F [label="AutoDock Vina"];
    G [label="Interaction Profiling"];
    H [label="Consensus Ranking"];
    I [label="Top Candidates"];
    J [label="Molecular Dynamics"];
    K [label="RMSD / H-Bonds / Contacts"];
    L [label="Lead Prioritization"];

    A -> B -> C -> D -> E -> F -> G -> H -> I -> J -> K -> L;
}
"""

    st.graphviz_chart(
        graph,
        width="stretch",
    )

    st.markdown(
        """
### Pipeline purpose

DockForgeAI accepts a selected protein target and chemical
library, predicts activity, filters unreliable or problematic
compounds, generates docking poses, integrates independent evidence
streams and validates shortlisted candidates with molecular dynamics.
"""
    )


# ============================================================
# CANDIDATE EXPLORER
# ============================================================

with candidate_tab:
    selected = st.selectbox(
        "Select candidate",
        CANDIDATES,
        index=CANDIDATES.index(active_candidate),
    )

    metrics = candidate_metrics(selected)
    shortlist = read_csv(SHORTLIST)

    row = pd.DataFrame()

    if (
        not shortlist.empty
        and "compound_id" in shortlist.columns
    ):
        row = shortlist[
            shortlist["compound_id"].astype(str)
            == selected
        ]

    st.markdown(
        '<div class="status-card"><b>Completed MD evidence package</b></div>',
        unsafe_allow_html=True,
    )

    st.write("")

    cards = st.columns(6)

    rank = "Reference"

    if (
        not row.empty
        and "consensus_rank" in row.columns
    ):
        rank = int(row.iloc[0]["consensus_rank"])

    cards[0].metric("Consensus rank", rank)

    predicted_pki = "—"
    vina = "—"

    if not row.empty:
        if "predicted_pki" in row.columns:
            predicted_pki = f"{row.iloc[0]['predicted_pki']:.3f}"

        if "vina_score_kcal_mol" in row.columns:
            vina = f"{row.iloc[0]['vina_score_kcal_mol']:.3f}"

    cards[1].metric("Predicted pKi", predicted_pki)
    cards[2].metric("Vina score", vina)

    cards[3].metric(
        "Ligand RMSD",
        f"{metrics['ligand_rmsd']:.3f} nm",
    )

    cards[4].metric(
        "Mean H-bonds",
        f"{metrics['mean_hbonds']:.2f}",
    )

    cards[5].metric(
        "H-bond occupancy",
        f"{metrics['hbond_occupancy']:.1f}%",
    )

    smiles_column = None

    for possible in [
        "canonical_smiles",
        "smiles",
        "standardized_smiles",
    ]:
        if possible in shortlist.columns:
            smiles_column = possible
            break

    left, right = st.columns([0.8, 1.4])

    with left:
        st.subheader("Ligand structure")

        if (
            not row.empty
            and smiles_column
            and Chem is not None
            and Draw is not None
        ):
            smiles = str(row.iloc[0][smiles_column])
            molecule = Chem.MolFromSmiles(smiles)

            if molecule is not None:
                image = Draw.MolToImage(
                    molecule,
                    size=(600, 420),
                )

                st.image(
                    image,
                    width="stretch",
                )

                descriptors = pd.DataFrame(
                    [
                        {
                            "Property": "Molecular weight",
                            "Value": round(
                                Descriptors.MolWt(molecule),
                                2,
                            ),
                        },
                        {
                            "Property": "LogP",
                            "Value": round(
                                Descriptors.MolLogP(molecule),
                                2,
                            ),
                        },
                        {
                            "Property": "H-bond donors",
                            "Value": (
                                Descriptors.NumHDonors(molecule)
                            ),
                        },
                        {
                            "Property": "H-bond acceptors",
                            "Value": (
                                Descriptors.NumHAcceptors(
                                    molecule
                                )
                            ),
                        },
                    ]
                )

                st.dataframe(
                    descriptors,
                    width="stretch",
                    hide_index=True,
                )
        else:
            st.info(
                "RDKit or a compatible SMILES column is unavailable."
            )

    with right:
        st.subheader("Ligand RMSD")

        data = to_ns(
            load_xvg(
                result_file(selected, "ligand_rmsd")
            ),
            "ligand_rmsd",
        )

        figure = px.line(
            data,
            x="time_ns",
            y="value",
            labels={
                "time_ns": "Time (ns)",
                "value": "Ligand RMSD (nm)",
            },
        )

        figure.update_layout(height=430)

        st.plotly_chart(
            figure,
            width="stretch",
        )

    chart_left, chart_right = st.columns(2)

    with chart_left:
        st.subheader("Hydrogen bonds")

        data = to_ns(
            load_xvg(result_file(selected, "hbonds")),
            "hbonds",
        )

        figure = px.line(
            data,
            x="time_ns",
            y="value",
            labels={
                "time_ns": "Time (ns)",
                "value": "Hydrogen bonds",
            },
        )

        figure.update_layout(height=360)

        st.plotly_chart(
            figure,
            width="stretch",
        )

    with chart_right:
        st.subheader("Protein-ligand contacts")

        data = to_ns(
            load_xvg(result_file(selected, "contacts")),
            "contacts",
        )

        figure = px.line(
            data,
            x="time_ns",
            y="value",
            labels={
                "time_ns": "Time (ns)",
                "value": "Contacts within 0.4 nm",
            },
        )

        figure.update_layout(height=360)

        st.plotly_chart(
            figure,
            width="stretch",
        )

    st.subheader("Scientific summary")
    st.markdown(read_text(summary_path(selected)))


# ============================================================
# MD COMPARISON
# ============================================================

with md_tab:
    st.header("Comparative molecular-dynamics analysis")

    selected_candidates = st.multiselect(
        "Candidates",
        CANDIDATES,
        default=CANDIDATES,
    )

    selection = st.selectbox(
        "Metric",
        [
            "Ligand RMSD",
            "Chain A RMSD",
            "Chain B RMSD",
            "Hydrogen bonds",
            "Contacts",
            "Minimum distance",
        ],
    )

    mapping = {
        "Ligand RMSD": "ligand_rmsd",
        "Chain A RMSD": "chain_a",
        "Chain B RMSD": "chain_b",
        "Hydrogen bonds": "hbonds",
        "Contacts": "contacts",
        "Minimum distance": "distance",
    }

    metric = mapping[selection]
    figure = go.Figure()

    for candidate in selected_candidates:
        data = to_ns(
            load_xvg(result_file(candidate, metric)),
            metric,
        )

        if data.empty:
            continue

        figure.add_trace(
            go.Scatter(
                x=data["time_ns"],
                y=data["value"],
                mode="lines",
                name=candidate,
            )
        )

    figure.update_layout(
        height=540,
        xaxis_title="Time (ns)",
        yaxis_title=selection,
        hovermode="x unified",
    )

    st.plotly_chart(
        figure,
        width="stretch",
    )

    comparison = comparison_table()

    plot_left, plot_right = st.columns(2)

    with plot_left:
        figure = px.bar(
            comparison.sort_values("ligand_rmsd"),
            x="compound_id",
            y="ligand_rmsd",
            labels={
                "compound_id": "Candidate",
                "ligand_rmsd": (
                    "Stable-period ligand RMSD (nm)"
                ),
            },
        )

        figure.update_layout(height=420)

        st.plotly_chart(
            figure,
            width="stretch",
        )

    with plot_right:
        figure = px.bar(
            comparison.sort_values(
                "hbond_occupancy"
            ),
            x="compound_id",
            y="hbond_occupancy",
            labels={
                "compound_id": "Candidate",
                "hbond_occupancy": (
                    "H-bond occupancy (%)"
                ),
            },
        )

        figure.update_layout(height=420)

        st.plotly_chart(
            figure,
            width="stretch",
        )


# ============================================================
# 3D VIEWER
# ============================================================

with viewer_tab:
    st.header("Interactive protein-ligand structure")

    candidate = st.selectbox(
        "Candidate",
        CANDIDATES,
        key="viewer_candidate",
    )

    structure_type = st.radio(
        "Structure",
        ["Representative", "Final"],
        horizontal=True,
    )

    path = structure_path(
        candidate,
        structure_type,
    )

    st.caption(str(path))
    render_3d(path)

    st.info(
        "Drag to rotate, scroll to zoom and right-drag to translate."
    )


# ============================================================
# PIPELINE
# ============================================================

with pipeline_tab:
    st.header("Pipeline execution")

    dry_run = st.checkbox(
        "Dry run",
        value=True,
    )

    resume = st.checkbox(
        "Resume completed stages",
        value=False,
    )

    stages = stage_names()

    selected_stage = st.selectbox(
        "Start stage",
        ["Beginning"] + stages,
    )

    from_stage = (
        None
        if selected_stage == "Beginning"
        else selected_stage
    )

    if st.button(
        "▶ Run DockForgeAI",
        type="primary",
        width="stretch",
    ):
        execute_pipeline(
            dry_run,
            resume,
            from_stage,
        )

    with st.expander("View pipeline configuration"):
        st.code(
            read_text(CONFIG),
            language="yaml",
        )


# ============================================================
# LOGS
# ============================================================

with logs_tab:
    st.header("Runs and logs")

    if not RUNS.exists():
        st.info("No run records found.")
    else:
        run_dirs = sorted(
            [
                path
                for path in RUNS.iterdir()
                if path.is_dir()
            ],
            reverse=True,
        )

        if run_dirs:
            selected_run = st.selectbox(
                "Run",
                [path.name for path in run_dirs],
            )

            run_dir = RUNS / selected_run
            manifest = run_dir / "run_manifest.json"

            if manifest.exists():
                st.json(json.loads(read_text(manifest)))

            log_dir = run_dir / "logs"

            if log_dir.exists():
                logs = sorted(log_dir.glob("*.log"))

                if logs:
                    selected_log = st.selectbox(
                        "Log",
                        [path.name for path in logs],
                    )

                    st.code(
                        read_text(log_dir / selected_log),
                        language="text",
                    )


# ============================================================
# REPOSITORY
# ============================================================

with repository_tab:
    st.header("Repository monitoring")

    code, output = shell(["git", "status"])

    if code == 0:
        st.code(output, language="text")
    else:
        st.error(output)

    st.subheader("Recent commits")

    code, output = shell(
        ["git", "log", "--oneline", "-15"]
    )

    if code == 0:
        st.code(output, language="text")
    else:
        st.error(output)


st.divider()

st.caption(
    "DockForgeAI v1.0 — computational prioritization only. "
    "Experimental validation is required before biological conclusions."
)
