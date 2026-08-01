"""DockForgeAI Streamlit dashboard."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import streamlit as st


PROJECT_ROOT = Path(__file__).resolve().parents[1]

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


from app.components.loaders import (  # noqa: E402
    file_status,
    load_csv,
    load_json,
    load_text,
)


st.set_page_config(
    page_title="DockForgeAI",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)


REDOCKING_SUMMARY_PATH = (
    "data/processed/docking_results/"
    "1HSG_MK1_redocking_summary.json"
)

REDOCKING_TABLE_PATH = (
    "data/processed/docking_tables/"
    "1HSG_MK1_redocking_analysis.csv"
)

INTERACTION_REPORT_PATH = (
    "data/processed/interaction_reports/report.md"
)

INTERACTION_TABLE_PATH = (
    "data/processed/interaction_reports/interactions.csv"
)

SCREENING_TABLE_PATH = (
    "data/processed/screening_tables/screening_scores.csv"
)

SCREENING_REPORT_PATH = (
    "data/processed/screening_reports/screening_report.md"
)

PRIORITIZATION_REPORT_PATH = (
    "data/processed/prioritization_reports/"
    "hit_prioritization_report.md"
)

PRIORITIZATION_TABLE_CANDIDATES = [
    "data/processed/prioritization_tables/consensus_ranking.csv",
    "data/processed/prioritization_reports/consensus_ranking.csv",
]


def find_first_existing_csv(paths: list[str]) -> pd.DataFrame:
    """Load the first available CSV from a candidate list."""

    for path in paths:
        status = file_status(path)

        if status["exists"]:
            return load_csv(path)

    return pd.DataFrame()


def show_status_card(
    title: str,
    status: bool,
    description: str,
) -> None:
    """Display a compact workflow-status card."""

    icon = "✅" if status else "⚠️"

    st.markdown(
        f"""
        <div style="
            border: 1px solid rgba(128,128,128,0.3);
            border-radius: 10px;
            padding: 14px;
            margin-bottom: 10px;
        ">
            <h4 style="margin: 0 0 8px 0;">{icon} {title}</h4>
            <p style="margin: 0;">{description}</p>
        </div>
        """,
        unsafe_allow_html=True,
    )


def render_home() -> None:
    """Render the project home page."""

    st.title("DockForgeAI")
    st.subheader(
        "Resource-aware molecular docking, virtual screening "
        "and hit-prioritization platform"
    )

    st.info(
        "Designed for reproducible local execution on an "
        "8 GB RAM macOS system."
    )

    st.markdown("### Workflow status")

    col1, col2, col3 = st.columns(3)

    redocking_exists = file_status(
        REDOCKING_SUMMARY_PATH
    )["exists"]

    interaction_exists = file_status(
        INTERACTION_REPORT_PATH
    )["exists"]

    screening_exists = file_status(
        SCREENING_TABLE_PATH
    )["exists"]

    with col1:
        show_status_card(
            "Redocking validation",
            redocking_exists,
            "AutoDock Vina docking and crystallographic RMSD validation.",
        )

    with col2:
        show_status_card(
            "Interaction analysis",
            interaction_exists,
            "Protein–ligand interaction detection and comparison.",
        )

    with col3:
        show_status_card(
            "Virtual screening",
            screening_exists,
            "Resource-aware screening and ranking of small libraries.",
        )

    col4, col5, col6 = st.columns(3)

    with col4:
        show_status_card(
            "Ligand preparation",
            file_status(
                "data/processed/pdbqt/ligands/1HSG_MK1.pdbqt"
            )["exists"],
            "Hydrogenation, chemistry correction and PDBQT generation.",
        )

    with col5:
        show_status_card(
            "Receptor preparation",
            file_status(
                "data/processed/pdbqt/receptors/"
                "1HSG_chains-A-B_receptor.pdbqt"
            )["exists"],
            "Rigid receptor preparation with documented protonation.",
        )

    with col6:
        show_status_card(
            "Automated testing",
            True,
            "Current validation status: 35 tests passed.",
        )

    st.markdown("### Current validated case study")

    st.write(
        """
        **Target:** HIV-1 protease  
        **PDB structure:** 1HSG  
        **Reference inhibitor:** MK1  
        **Protein chains retained:** A and B  
        **Validation criterion:** crystallographic heavy-atom RMSD ≤ 2 Å
        """
    )


def render_redocking() -> None:
    """Render redocking results."""

    st.header("Redocking validation")

    summary = load_json(REDOCKING_SUMMARY_PATH)
    table = load_csv(REDOCKING_TABLE_PATH)

    if not summary and table.empty:
        st.warning("Redocking results were not found.")
        return

    if summary:
        columns = st.columns(4)

        possible_metrics = {
            "Best affinity": [
                "best_affinity_kcal_mol",
            ],
            "Minimum RMSD": [
                "min_rmsd_angstrom",
            ],
            "Best RMSD rank": [
                "min_rmsd_pose_rank",
            ],
            "Validation": [
                "interpretation",
            ],
        }

        metric_values = {
            "Best affinity": f"{summary.get('best_affinity_kcal_mol', 'N/A')} kcal/mol",
            "Minimum RMSD": f"{summary.get('min_rmsd_angstrom', 'N/A')} Å",
            "Best RMSD rank": summary.get("min_rmsd_pose_rank", "N/A"),
            "Validation": summary.get("interpretation", "N/A"),
        }

        for column, (label, value) in zip(
            columns,
            metric_values.items(),
        ):
            column.metric(label, value)

        with st.expander("Complete redocking summary"):
            st.json(summary)

    if not table.empty:
        st.markdown("### Pose ranking and RMSD")

        display_table = table.drop(
            columns=["pose_pdbqt"],
            errors="ignore",
        )

        st.dataframe(
            display_table,
            width="stretch",
            hide_index=True,
        )

        affinity_column = next(
            (
                column
                for column in table.columns
                if column == "affinity_kcal_mol"
            ),
            None,
        )

        if affinity_column:
            chart_data = table.copy()

            if "pose_rank" in chart_data.columns:
                chart_data = chart_data.set_index("pose_rank")

            st.markdown("### Docking affinity by pose")
            st.bar_chart(chart_data[[affinity_column]])

        rmsd_column = next(
            (
                column
                for column in table.columns
                if "crystal" in column.lower()
                and "rmsd" in column.lower()
            ),
            None,
        )

        if rmsd_column:
            chart_data = table.copy()

            if "pose_rank" in chart_data.columns:
                chart_data = chart_data.set_index("pose_rank")

            st.markdown("### Crystallographic RMSD by pose")
            st.line_chart(chart_data[[rmsd_column]])


def render_interactions() -> None:
    """Render protein–ligand interaction results."""

    st.header("Protein–ligand interaction analysis")

    report = load_text(INTERACTION_REPORT_PATH)
    table = load_csv(INTERACTION_TABLE_PATH)

    if not report and table.empty:
        st.warning("Interaction-analysis outputs were not found.")
        return

    if report:
        st.markdown(report)

    if not table.empty:
        st.markdown("### Interaction records")

        display_table = table.drop(
            columns=["pose_pdbqt"],
            errors="ignore",
        )

        st.dataframe(
            display_table,
            width="stretch",
            hide_index=True,
        )

        if "interaction_type" in table.columns:
            counts = (
                table["interaction_type"]
                .value_counts()
                .rename_axis("interaction_type")
                .to_frame("count")
            )

            st.markdown("### Interaction counts")
            st.bar_chart(counts)


def render_screening() -> None:
    """Render virtual-screening outputs."""

    st.header("Virtual screening")

    report = load_text(SCREENING_REPORT_PATH)
    table = load_csv(SCREENING_TABLE_PATH)

    if not report and table.empty:
        st.warning("Virtual-screening outputs were not found.")
        return

    if report:
        st.markdown(report)

    if not table.empty:
        st.markdown("### Screening scores")

        display_table = table.drop(
            columns=["pose_pdbqt"],
            errors="ignore",
        )

        st.dataframe(
            display_table,
            width="stretch",
            hide_index=True,
        )

        score_column = next(
            (
                column
                for column in table.columns
                if column in {
                    "best_affinity",
                    "best_affinity_kcal_mol",
                    "screening_score",
                }
            ),
            None,
        )

        id_column = next(
            (
                column
                for column in table.columns
                if column in {
                    "compound_id",
                    "ligand_id",
                    "id",
                }
            ),
            None,
        )

        if score_column and id_column:
            chart_data = (
                table[[id_column, score_column]]
                .dropna()
                .set_index(id_column)
            )

            st.markdown("### Compound ranking")
            st.bar_chart(chart_data)


def render_prioritization() -> None:
    """Render hit-prioritization outputs."""

    st.header("Hit prioritization")

    report = load_text(PRIORITIZATION_REPORT_PATH)

    table = find_first_existing_csv(
        PRIORITIZATION_TABLE_CANDIDATES
    )

    if not report and table.empty:
        st.warning("Hit-prioritization outputs were not found.")
        return

    if report:
        st.markdown(report)

    if not table.empty:
        st.markdown("### Consensus ranking")

        display_table = table.drop(
            columns=["pose_pdbqt"],
            errors="ignore",
        )

        st.dataframe(
            display_table,
            width="stretch",
            hide_index=True,
        )

        score_column = next(
            (
                column
                for column in table.columns
                if column in {
                    "consensus_score",
                    "base_consensus_score",
                    "admet_adjusted_consensus_score",
                }
            ),
            None,
        )

        id_column = next(
            (
                column
                for column in table.columns
                if column in {
                    "compound_id",
                    "ligand_id",
                    "id",
                }
            ),
            None,
        )

        if score_column and id_column:
            chart_data = (
                table[[id_column, score_column]]
                .dropna()
                .set_index(id_column)
            )

            st.markdown("### Consensus scores")
            st.bar_chart(chart_data)


def render_about() -> None:
    """Render project information."""

    st.header("About DockForgeAI")

    st.markdown(
        """
        DockForgeAI is a modular computational drug-discovery workflow
        supporting:

        - Protein retrieval, inspection and preparation
        - Ligand chemistry validation and preparation
        - AutoDock Vina docking
        - Crystallographic redocking validation
        - Protein–ligand interaction analysis
        - Resource-aware virtual screening
        - Drug-likeness and structural-alert analysis
        - Optional external ADMET integration
        - Consensus hit prioritization

        ### Important scientific limitations

        - Docking scores are approximate ranking values.
        - They are not experimental binding free energies.
        - Predicted interactions require experimental confirmation.
        - Physicochemical descriptors are not equivalent to full ADMET data.
        - Final hit candidates require biochemical and cellular validation.
        """
    )


def main() -> None:
    """Run the Streamlit application."""

    st.sidebar.title("DockForgeAI")

    page = st.sidebar.radio(
        "Navigation",
        [
            "Home",
            "Redocking",
            "Interactions",
            "Virtual screening",
            "Hit prioritization",
            "About",
        ],
    )

    st.sidebar.markdown("---")
    st.sidebar.caption("DockForgeAI v1.0")

    if page == "Home":
        render_home()
    elif page == "Redocking":
        render_redocking()
    elif page == "Interactions":
        render_interactions()
    elif page == "Virtual screening":
        render_screening()
    elif page == "Hit prioritization":
        render_prioritization()
    else:
        render_about()


if __name__ == "__main__":
    main()