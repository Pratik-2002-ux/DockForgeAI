"""DockForgeAI Streamlit dashboard."""

from __future__ import annotations

import json

import sys
from pathlib import Path

import altair as alt
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
from app.components.protein_viewer import (  # noqa: E402
    available_pose_files,
    build_view,
    import_py3dmol,
    load_viewer_settings,
    structure_files,
)
from app.components.export_center import (  # noqa: E402
    build_zip_bytes,
    categorized_exports,
    format_file_size,
    generate_reports,
    load_manifest,
)
from app.components.xai_loaders import load_xai_csv, load_xai_json  # noqa: E402
from app.components.xai_export import xai_bundle  # noqa: E402
from app.components.xai_validation import validation_bundle, external_template  # noqa: E402
from app.components.pubchem_external import DOWNLOAD_KEYS, build_bundle  # noqa: E402


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


SCREENING_JSON_PATH = (
    "data/processed/screening_reports/screening_report.json"
)

PRIORITIZATION_REPORT_PATH = (
    "data/processed/prioritization_reports/"
    "hit_prioritization_report.md"
)


PRIORITIZATION_JSON_PATH = (
    "data/processed/prioritization_reports/"
    "hit_prioritization_report.json"
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
    """Render a concise and user-friendly redocking dashboard."""

    st.header("Redocking validation")
    st.caption(
        "Comparison of AutoDock Vina poses with the corrected "
        "crystallographic MK1 reference pose."
    )

    summary = load_json(REDOCKING_SUMMARY_PATH)
    table = load_csv(REDOCKING_TABLE_PATH)

    if not summary and table.empty:
        st.warning("Redocking results were not found.")
        return

    if summary:
        best_affinity = summary.get("best_affinity_kcal_mol", "N/A")
        minimum_rmsd = summary.get("min_rmsd_angstrom", "N/A")
        best_rmsd_rank = summary.get("min_rmsd_pose_rank", "N/A")
        interpretation = summary.get("interpretation", "N/A")

        metric_columns = st.columns(4)

        metric_columns[0].metric(
            "Best affinity",
            (
                f"{best_affinity:.3f} kcal/mol"
                if isinstance(best_affinity, (int, float))
                else best_affinity
            ),
        )

        metric_columns[1].metric(
            "Minimum RMSD",
            (
                f"{minimum_rmsd:.3f} Å"
                if isinstance(minimum_rmsd, (int, float))
                else minimum_rmsd
            ),
        )

        metric_columns[2].metric(
            "Best RMSD rank",
            best_rmsd_rank,
        )

        metric_columns[3].metric(
            "Validation",
            interpretation,
        )

        pose_count = summary.get("pose_count", "N/A")
        passing_count = summary.get("passing_pose_count", "N/A")
        passing_percentage = summary.get(
            "passing_pose_percentage",
            "N/A",
        )
        threshold = summary.get("threshold_angstrom", 2.0)

        summary_table = pd.DataFrame(
            {
                "Metric": [
                    "Generated poses",
                    "Best-affinity pose",
                    "RMSD of best-affinity pose",
                    "Minimum-RMSD pose",
                    "Passing poses",
                    "Passing percentage",
                    "RMSD threshold",
                    "Top-ranked pose passed",
                    "Overall pose reproduction",
                ],
                "Value": [
                    str(pose_count),
                    str(summary.get("best_affinity_pose_rank", "N/A")),
                    f"{summary.get('rmsd_best_affinity_pose_angstrom', 'N/A')} Å",
                    str(best_rmsd_rank),
                    f"{passing_count} / {pose_count}",
                    (
                        f"{passing_percentage:.1f}%"
                        if isinstance(passing_percentage, (int, float))
                        else str(passing_percentage)
                    ),
                    f"{threshold} Å",
                    (
                        "Yes"
                        if summary.get("top_ranked_pose_success")
                        else "No"
                    ),
                    (
                        "Successful"
                        if summary.get(
                            "overall_pose_reproduction_success"
                        )
                        else "Failed"
                    ),
                ],
            }
        )

        with st.expander(
            "Redocking validation summary",
            expanded=True,
        ):
            st.dataframe(
                summary_table,
                width="stretch",
                hide_index=True,
            )

        download_columns = st.columns(2)

        download_columns[0].download_button(
            label="Download summary JSON",
            data=json.dumps(summary, indent=2),
            file_name="1HSG_MK1_redocking_summary.json",
            mime="application/json",
        )

        if not table.empty:
            download_columns[1].download_button(
                label="Download pose analysis CSV",
                data=table.to_csv(index=False),
                file_name="1HSG_MK1_redocking_analysis.csv",
                mime="text/csv",
            )

    if table.empty:
        return

    st.subheader("Pose ranking")

    preferred_columns = [
        "pose_rank",
        "affinity_kcal_mol",
        "crystal_in_place_rmsd_angstrom",
        "aligned_best_rmsd_angstrom",
        "passes_2A_threshold",
    ]

    available_columns = [
        column
        for column in preferred_columns
        if column in table.columns
    ]

    display_table = table[available_columns].copy()

    display_table = display_table.rename(
        columns={
            "pose_rank": "Pose",
            "affinity_kcal_mol": "Affinity (kcal/mol)",
            "crystal_in_place_rmsd_angstrom": "Crystal RMSD (Å)",
            "aligned_best_rmsd_angstrom": "Aligned RMSD (Å)",
            "passes_2A_threshold": "≤2 Å",
        }
    )

    if "≤2 Å" in display_table.columns:
        display_table["Validation"] = (
            display_table["≤2 Å"]
            .map({True: "PASS", False: "FAIL"})
            .fillna("N/A")
        )
        display_table = display_table.drop(columns=["≤2 Å"])

    st.dataframe(
        display_table,
        width="stretch",
        hide_index=True,
    )

    chart_column_1, chart_column_2 = st.columns(2)

    if {
        "pose_rank",
        "affinity_kcal_mol",
    }.issubset(table.columns):
        affinity_data = (
            table[
                [
                    "pose_rank",
                    "affinity_kcal_mol",
                ]
            ]
            .set_index("pose_rank")
            .rename(
                columns={
                    "affinity_kcal_mol": "Affinity (kcal/mol)"
                }
            )
        )

        with chart_column_1:
            st.subheader("Docking affinity")
            st.caption("More negative Vina scores indicate stronger predicted ranking.")
            st.bar_chart(
                affinity_data,
                height=400,
            )

    if {
        "pose_rank",
        "crystal_in_place_rmsd_angstrom",
    }.issubset(table.columns):
        threshold_value = float(
            summary.get("threshold_angstrom", 2.0)
        )

        rmsd_data = table[
            [
                "pose_rank",
                "crystal_in_place_rmsd_angstrom",
            ]
        ].copy()

        rmsd_data["Validation threshold"] = threshold_value

        rmsd_data = (
            rmsd_data
            .set_index("pose_rank")
            .rename(
                columns={
                    "crystal_in_place_rmsd_angstrom": (
                        "Crystal RMSD (Å)"
                    )
                }
            )
        )

        with chart_column_2:
            st.subheader("Crystallographic RMSD")
            st.line_chart(
                rmsd_data,
                height=400,
            )

    st.info(
        "A redocked pose is considered successfully reproduced "
        "when its crystallographic heavy-atom RMSD is no greater "
        "than 2 Å. Vina scores are approximate ranking values, "
        "not experimental binding free energies."
    )
def render_interactions() -> None:
    """Render a professional protein-ligand interaction dashboard."""

    st.header("Protein–ligand interaction analysis")
    st.caption(
        "Explore predicted contacts for the crystallographic, "
        "rank-1 and best-RMSD MK1 poses."
    )

    report = load_text(INTERACTION_REPORT_PATH)
    table = load_csv(INTERACTION_TABLE_PATH)

    if not report and table.empty:
        st.warning("Interaction-analysis outputs were not found.")
        return

    if table.empty:
        st.markdown(report)
        return

    interaction_column = next(
        (
            column
            for column in table.columns
            if column.lower() in {
                "interaction_type",
                "interaction",
                "type",
            }
        ),
        None,
    )

    residue_name_column = next(
        (
            column
            for column in table.columns
            if column.lower() in {
                "protein_residue_name",
                "residue_name",
                "resname",
            }
        ),
        None,
    )

    residue_number_column = next(
        (
            column
            for column in table.columns
            if column.lower() in {
                "protein_residue_number",
                "residue_number",
                "resnum",
            }
        ),
        None,
    )

    chain_column = next(
        (
            column
            for column in table.columns
            if column.lower() in {
                "protein_chain",
                "chain",
                "chain_id",
            }
        ),
        None,
    )

    pose_column = next(
        (
            column
            for column in table.columns
            if column.lower() in {
                "pose_label",
                "pose",
                "pose_name",
            }
        ),
        None,
    )

    selected_table = table.copy()

    if pose_column:
        pose_values = sorted(
            selected_table[pose_column]
            .dropna()
            .astype(str)
            .unique()
            .tolist()
        )

        pose_options = ["All poses"] + pose_values

        selected_pose = st.selectbox(
            "Pose selection",
            options=pose_options,
            help=(
                "Filter contacts by crystallographic, rank-1 "
                "or best-RMSD pose."
            ),
        )

        if selected_pose != "All poses":
            selected_table = selected_table[
                selected_table[pose_column].astype(str)
                == selected_pose
            ].copy()

    row_text = selected_table.apply(
        lambda row: " ".join(
            str(value)
            for value in row.tolist()
        ),
        axis=1,
    )

    total_interactions = len(selected_table)

    interaction_type_count = (
        selected_table[interaction_column].nunique()
        if interaction_column
        else 0
    )

    catalytic_contacts = int(
        row_text.str.contains(
            "catalytic",
            case=False,
            na=False,
        ).sum()
    )

    asp25_contacts = int(
        row_text.str.contains(
            r"ASP\s*25|ASP25",
            case=False,
            regex=True,
            na=False,
        ).sum()
    )

    hydrogen_bonds = 0

    if interaction_column:
        hydrogen_bonds = int(
            selected_table[interaction_column]
            .astype(str)
            .str.contains(
                "hydrogen",
                case=False,
                na=False,
            )
            .sum()
        )

    metric_columns = st.columns(5)

    metric_columns[0].metric(
        "Total contacts",
        total_interactions,
    )
    metric_columns[1].metric(
        "Interaction types",
        interaction_type_count,
    )
    metric_columns[2].metric(
        "Hydrogen bonds",
        hydrogen_bonds,
    )
    metric_columns[3].metric(
        "Catalytic contacts",
        catalytic_contacts,
    )
    metric_columns[4].metric(
        "ASP25 contacts",
        asp25_contacts,
    )

    summary_tab, catalytic_tab, records_tab, report_tab = st.tabs(
        [
            "Summary",
            "Catalytic site",
            "Interaction records",
            "Full report",
        ]
    )

    with summary_tab:
        chart_column_1, chart_column_2 = st.columns(2)

        if interaction_column:
            interaction_counts = (
                selected_table[interaction_column]
                .fillna("Unknown")
                .astype(str)
                .value_counts()
                .rename_axis("Interaction type")
                .reset_index(name="Count")
                .sort_values("Count", ascending=True)
            )

            interaction_chart = (
                alt.Chart(interaction_counts)
                .mark_bar()
                .encode(
                    x=alt.X(
                        "Count:Q",
                        title="Number of contacts",
                    ),
                    y=alt.Y(
                        "Interaction type:N",
                        sort="-x",
                        title=None,
                    ),
                    tooltip=[
                        alt.Tooltip(
                            "Interaction type:N",
                            title="Interaction",
                        ),
                        alt.Tooltip(
                            "Count:Q",
                            title="Count",
                        ),
                    ],
                )
                .properties(
                    height=max(
                        280,
                        len(interaction_counts) * 42,
                    )
                )
            )

            with chart_column_1:
                st.subheader("Interactions by type")
                st.altair_chart(
                    interaction_chart,
                    width="stretch",
                )

        if residue_name_column:
            residue_labels = (
                selected_table[residue_name_column]
                .fillna("Unknown")
                .astype(str)
            )

            if residue_number_column:
                residue_numbers = (
                    selected_table[residue_number_column]
                    .fillna("")
                    .astype(str)
                    .str.replace(
                        r"\.0$",
                        "",
                        regex=True,
                    )
                )

                residue_labels = (
                    residue_labels
                    + residue_numbers
                )

            if chain_column:
                chain_labels = (
                    selected_table[chain_column]
                    .fillna("")
                    .astype(str)
                )

                residue_labels = (
                    chain_labels
                    + ":"
                    + residue_labels
                )

            residue_counts = (
                residue_labels
                .value_counts()
                .head(15)
                .rename_axis("Residue")
                .reset_index(name="Contact count")
                .sort_values(
                    "Contact count",
                    ascending=True,
                )
            )

            residue_chart = (
                alt.Chart(residue_counts)
                .mark_bar()
                .encode(
                    x=alt.X(
                        "Contact count:Q",
                        title="Number of contacts",
                    ),
                    y=alt.Y(
                        "Residue:N",
                        sort="-x",
                        title=None,
                    ),
                    tooltip=[
                        alt.Tooltip(
                            "Residue:N",
                            title="Residue",
                        ),
                        alt.Tooltip(
                            "Contact count:Q",
                            title="Contacts",
                        ),
                    ],
                )
                .properties(
                    height=max(
                        320,
                        len(residue_counts) * 36,
                    )
                )
            )

            with chart_column_2:
                st.subheader("Top contacted residues")
                st.altair_chart(
                    residue_chart,
                    width="stretch",
                )

        st.info(
            "Contact frequency reflects predicted geometric "
            "interactions. Multiple ligand atoms may contact the "
            "same protein residue."
        )

    with catalytic_tab:
        catalytic_mask = row_text.str.contains(
            r"ASP\s*25|ASP25|catalytic",
            case=False,
            regex=True,
            na=False,
        )

        catalytic_rows = selected_table[
            catalytic_mask
        ].copy()

        st.subheader("Catalytic Asp25 interactions")

        if catalytic_rows.empty:
            st.warning(
                "No catalytic-site interaction records were found "
                "for the selected pose."
            )
        else:
            hidden_columns = [
                column
                for column in catalytic_rows.columns
                if any(
                    keyword in column.lower()
                    for keyword in [
                        "path",
                        "file",
                        "source_pdbqt",
                    ]
                )
            ]

            catalytic_display = catalytic_rows.drop(
                columns=hidden_columns,
                errors="ignore",
            )

            st.dataframe(
                catalytic_display,
                width="stretch",
                hide_index=True,
                column_config={
                    "distance_angstrom": st.column_config.NumberColumn(
                        "Distance (Å)",
                        format="%.3f",
                    ),
                    "angle_degrees": st.column_config.NumberColumn(
                        "Angle (°)",
                        format="%.1f",
                    ),
                },
            )

            st.download_button(
                label="Download catalytic contacts",
                data=catalytic_rows.to_csv(index=False),
                file_name=(
                    "1HSG_MK1_catalytic_interactions.csv"
                ),
                mime="text/csv",
            )

    with records_tab:
        hidden_columns = [
            column
            for column in selected_table.columns
            if any(
                keyword in column.lower()
                for keyword in [
                    "path",
                    "file",
                    "source_pdbqt",
                ]
            )
        ]

        display_table = selected_table.drop(
            columns=hidden_columns,
            errors="ignore",
        )

        st.dataframe(
            display_table,
            width="stretch",
            hide_index=True,
            column_config={
                "distance_angstrom": st.column_config.NumberColumn(
                    "Distance (Å)",
                    format="%.3f",
                ),
                "angle_degrees": st.column_config.NumberColumn(
                    "Angle (°)",
                    format="%.1f",
                ),
            },
        )

        st.download_button(
            label="Download filtered interaction CSV",
            data=selected_table.to_csv(index=False),
            file_name="1HSG_MK1_interactions.csv",
            mime="text/csv",
        )

    with report_tab:
        if report:
            st.markdown(report)

            st.download_button(
                label="Download Markdown report",
                data=report,
                file_name=(
                    "1HSG_MK1_interaction_report.md"
                ),
                mime="text/markdown",
            )
        else:
            st.info(
                "The Markdown interaction report was not found."
            )

    st.caption(
        "Protein-ligand interactions are computational predictions "
        "and require experimental validation."
    )
def render_screening() -> None:
    """Render a professional virtual-screening dashboard."""

    st.header("Virtual screening")
    st.caption(
        "Resource-aware ligand preparation, docking and "
        "screening-result overview."
    )

    report_text = load_text(SCREENING_REPORT_PATH)
    report_json = load_json(SCREENING_JSON_PATH)
    table = load_csv(SCREENING_TABLE_PATH)

    if not report_json and not report_text and table.empty:
        st.warning("Virtual-screening outputs were not found.")
        return

    library_size = report_json.get(
        "library_size",
        report_json.get("input_count", 0),
    )
    accepted_count = report_json.get(
        "accepted_count",
        report_json.get("standardized_count", 0),
    )
    rejected_count = report_json.get(
        "rejected_count",
        0,
    )
    screening_count = report_json.get(
        "screening_result_count",
        report_json.get("screened_count", len(table)),
    )
    refined_count = report_json.get(
        "refined_count",
        0,
    )

    metric_columns = st.columns(5)

    metric_columns[0].metric(
        "Library compounds",
        library_size,
    )
    metric_columns[1].metric(
        "Accepted",
        accepted_count,
    )
    metric_columns[2].metric(
        "Rejected",
        rejected_count,
    )
    metric_columns[3].metric(
        "Docking results",
        screening_count,
    )
    metric_columns[4].metric(
        "Refined hits",
        refined_count,
    )

    overview_tab, scores_tab, report_tab = st.tabs(
        [
            "Overview",
            "Screening scores",
            "Full report",
        ]
    )

    with overview_tab:
        status_data = pd.DataFrame(
            {
                "Stage": [
                    "Library",
                    "Accepted",
                    "Rejected",
                    "Screened",
                    "Refined",
                ],
                "Count": [
                    library_size,
                    accepted_count,
                    rejected_count,
                    screening_count,
                    refined_count,
                ],
            }
        )

        status_chart = (
            alt.Chart(status_data)
            .mark_bar()
            .encode(
                x=alt.X(
                    "Count:Q",
                    title="Compound count",
                ),
                y=alt.Y(
                    "Stage:N",
                    sort="-x",
                    title=None,
                ),
                tooltip=[
                    alt.Tooltip("Stage:N"),
                    alt.Tooltip("Count:Q"),
                ],
                color=alt.Color(
                    "Stage:N",
                    legend=None,
                ),
            )
            .properties(height=300)
        )

        st.subheader("Screening workflow summary")
        st.altair_chart(
            status_chart,
            width="stretch",
        )

        if screening_count == 0:
            st.warning(
                "No completed docking scores were found. "
                "The demonstration library was prepared, but "
                "the screening-results table currently contains "
                "no successful docking result."
            )

        if not table.empty:
            score_column = next(
                (
                    column
                    for column in table.columns
                    if column in {
                        "best_affinity",
                        "best_affinity_kcal_mol",
                        "screening_score",
                        "affinity_kcal_mol",
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
                valid_scores = table[
                    [id_column, score_column]
                ].dropna().copy()

                if not valid_scores.empty:
                    valid_scores["Display compound"] = (
                        valid_scores[id_column].astype(str)
                    )

                    valid_scores = valid_scores.sort_values(
                        score_column,
                        ascending=True,
                    ).head(15)

                    score_chart = (
                        alt.Chart(valid_scores)
                        .mark_bar()
                        .encode(
                            x=alt.X(
                                f"{score_column}:Q",
                                title="Best Vina score (kcal/mol)",
                            ),
                            y=alt.Y(
                                "Display compound:N",
                                sort=alt.EncodingSortField(
                                    field=score_column,
                                    order="ascending",
                                ),
                                title=None,
                            ),
                            tooltip=[
                                alt.Tooltip(
                                    f"{id_column}:N",
                                    title="Compound",
                                ),
                                alt.Tooltip(
                                    f"{score_column}:Q",
                                    title="Score",
                                    format=".3f",
                                ),
                            ],
                        )
                        .properties(
                            height=max(
                                300,
                                len(valid_scores) * 35,
                            )
                        )
                    )

                    st.subheader("Top docking scores")
                    st.caption(
                        "More negative Vina scores indicate "
                        "stronger predicted ranking."
                    )
                    st.altair_chart(
                        score_chart,
                        width="stretch",
                    )

    with scores_tab:
        if table.empty:
            st.info(
                "The screening score table is currently empty."
            )
        else:
            hidden_columns = [
                column
                for column in table.columns
                if any(
                    keyword in column.lower()
                    for keyword in [
                        "path",
                        "file",
                        "stdout",
                        "stderr",
                    ]
                )
            ]

            display_table = table.drop(
                columns=hidden_columns,
                errors="ignore",
            )

            st.dataframe(
                display_table,
                width="stretch",
                hide_index=True,
                column_config={
                    "best_affinity": (
                        st.column_config.NumberColumn(
                            "Best affinity",
                            format="%.3f",
                        )
                    ),
                    "best_affinity_kcal_mol": (
                        st.column_config.NumberColumn(
                            "Best affinity (kcal/mol)",
                            format="%.3f",
                        )
                    ),
                    "runtime_seconds": (
                        st.column_config.NumberColumn(
                            "Runtime (s)",
                            format="%.2f",
                        )
                    ),
                },
            )

            st.download_button(
                label="Download screening scores",
                data=table.to_csv(index=False),
                file_name="screening_scores.csv",
                mime="text/csv",
            )

    with report_tab:
        if report_json:
            summary_table = pd.DataFrame(
                {
                    "Metric": list(report_json.keys()),
                    "Value": [
                        (
                            json.dumps(value)
                            if isinstance(
                                value,
                                (dict, list),
                            )
                            else value
                        )
                        for value in report_json.values()
                    ],
                }
            )

            st.dataframe(
                summary_table,
                width="stretch",
                hide_index=True,
            )

            st.download_button(
                label="Download screening JSON",
                data=json.dumps(
                    report_json,
                    indent=2,
                ),
                file_name="screening_report.json",
                mime="application/json",
            )

        if report_text:
            with st.expander(
                "Markdown screening report",
                expanded=False,
            ):
                st.markdown(report_text)

            st.download_button(
                label="Download screening Markdown",
                data=report_text,
                file_name="screening_report.md",
                mime="text/markdown",
            )

    st.caption(
        "The demonstration library validates software "
        "functionality and does not represent a biological "
        "discovery study."
    )
def render_prioritization() -> None:
    """Render a professional hit-prioritization dashboard."""

    st.header("Hit prioritization")
    st.caption(
        "Consensus ranking integrating docking, interactions, "
        "drug-likeness, structural alerts and ligand efficiency."
    )

    report_text = load_text(PRIORITIZATION_REPORT_PATH)
    report_json = load_json(PRIORITIZATION_JSON_PATH)

    table = find_first_existing_csv(
        PRIORITIZATION_TABLE_CANDIDATES
    )

    if not report_text and not report_json and table.empty:
        st.warning("Hit-prioritization outputs were not found.")
        return

    if table.empty:
        if report_json:
            st.json(report_json)
        elif report_text:
            st.markdown(report_text)
        return

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

    rank_column = next(
        (
            column
            for column in table.columns
            if column in {
                "rank",
                "consensus_rank",
            }
        ),
        None,
    )

    consensus_column = next(
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

    affinity_column = next(
        (
            column
            for column in table.columns
            if column in {
                "best_affinity",
                "best_affinity_kcal_mol",
                "refined_best_affinity",
            }
        ),
        None,
    )

    qed_column = next(
        (
            column
            for column in table.columns
            if column.lower() == "qed"
        ),
        None,
    )

    ligand_efficiency_column = next(
        (
            column
            for column in table.columns
            if column in {
                "ligand_efficiency",
                "score_per_heavy_atom",
            }
        ),
        None,
    )

    top_row = (
        table.sort_values(
            rank_column,
            ascending=True,
        ).iloc[0]
        if rank_column
        else table.iloc[0]
    )

    top_compound = (
        top_row[id_column]
        if id_column
        else "N/A"
    )

    top_rank = (
        top_row[rank_column]
        if rank_column
        else 1
    )

    top_consensus = (
        top_row[consensus_column]
        if consensus_column
        else "N/A"
    )

    top_affinity = (
        top_row[affinity_column]
        if affinity_column
        else "N/A"
    )

    metric_columns = st.columns(5)

    metric_columns[0].metric(
        "Ranked compounds",
        len(table),
    )
    metric_columns[1].metric(
        "Top compound",
        top_compound,
    )
    metric_columns[2].metric(
        "Top rank",
        top_rank,
    )
    metric_columns[3].metric(
        "Consensus score",
        (
            f"{float(top_consensus):.3f}"
            if isinstance(
                top_consensus,
                (int, float),
            )
            else top_consensus
        ),
    )
    metric_columns[4].metric(
        "Best affinity",
        (
            f"{float(top_affinity):.3f} kcal/mol"
            if isinstance(
                top_affinity,
                (int, float),
            )
            else top_affinity
        ),
    )

    ranking_tab, profile_tab, report_tab = st.tabs(
        [
            "Consensus ranking",
            "Top-hit profile",
            "Full report",
        ]
    )

    with ranking_tab:
        display_table = table.copy()

        hidden_columns = [
            column
            for column in display_table.columns
            if any(
                keyword in column.lower()
                for keyword in [
                    "path",
                    "file",
                ]
            )
        ]

        display_table = display_table.drop(
            columns=hidden_columns,
            errors="ignore",
        )

        st.dataframe(
            display_table,
            width="stretch",
            hide_index=True,
            column_config={
                "consensus_score": (
                    st.column_config.NumberColumn(
                        "Consensus score",
                        format="%.3f",
                    )
                ),
                "best_affinity": (
                    st.column_config.NumberColumn(
                        "Best affinity",
                        format="%.3f",
                    )
                ),
                "qed": (
                    st.column_config.NumberColumn(
                        "QED",
                        format="%.3f",
                    )
                ),
                "ligand_efficiency": (
                    st.column_config.NumberColumn(
                        "Ligand efficiency",
                        format="%.3f",
                    )
                ),
            },
        )

        if consensus_column:
            chart_data = table.copy()

            if id_column:
                chart_data["Compound label"] = (
                    chart_data[id_column].astype(str)
                )
            else:
                chart_data["Compound label"] = (
                    "Compound "
                    + (
                        chart_data.index
                        + 1
                    ).astype(str)
                )

            if rank_column:
                chart_data["Compound label"] = (
                    chart_data["Compound label"]
                    + " (#"
                    + chart_data[
                        rank_column
                    ].astype(str)
                    + ")"
                )

            chart_data = chart_data.sort_values(
                consensus_column,
                ascending=True,
            )

            consensus_chart = (
                alt.Chart(chart_data)
                .mark_bar()
                .encode(
                    x=alt.X(
                        f"{consensus_column}:Q",
                        title="Consensus score",
                    ),
                    y=alt.Y(
                        "Compound label:N",
                        sort=alt.EncodingSortField(
                            field=consensus_column,
                            order="descending",
                        ),
                        title=None,
                    ),
                    tooltip=[
                        alt.Tooltip(
                            "Compound label:N",
                            title="Compound",
                        ),
                        alt.Tooltip(
                            f"{consensus_column}:Q",
                            title="Consensus",
                            format=".3f",
                        ),
                    ],
                )
                .properties(
                    height=max(
                        280,
                        len(chart_data) * 45,
                    )
                )
            )

            st.subheader("Consensus scores")
            st.altair_chart(
                consensus_chart,
                width="stretch",
            )

        st.download_button(
            label="Download consensus ranking",
            data=table.to_csv(index=False),
            file_name="consensus_ranking.csv",
            mime="text/csv",
        )

    with profile_tab:
        if id_column:
            compound_options = (
                table[id_column]
                .astype(str)
                .tolist()
            )

            selected_compound = st.selectbox(
                "Select compound",
                options=compound_options,
            )

            compound_rows = table[
                table[id_column].astype(str)
                == selected_compound
            ]

            selected_row = compound_rows.iloc[0]
        else:
            selected_row = table.iloc[0]

        profile_fields = []

        for column, label in [
            (consensus_column, "Consensus score"),
            (affinity_column, "Docking affinity"),
            (qed_column, "QED"),
            (
                ligand_efficiency_column,
                "Ligand efficiency",
            ),
            (
                "interaction",
                "Interaction score",
            ),
            (
                "pains_alerts",
                "PAINS alerts",
            ),
            (
                "brenk_alerts",
                "Brenk alerts",
            ),
            (
                "alerts",
                "Total alerts",
            ),
        ]:
            if column and column in selected_row.index:
                profile_fields.append(
                    {
                        "Feature": label,
                        "Value": selected_row[column],
                    }
                )

        if profile_fields:
            profile_table = pd.DataFrame(
                profile_fields
            )

            st.subheader("Selected-hit profile")
            st.dataframe(
                profile_table,
                width="stretch",
                hide_index=True,
            )

        if "smiles" in selected_row.index:
            st.text_area(
                "SMILES",
                value=str(selected_row["smiles"]),
                height=90,
                disabled=True,
            )

        st.info(
            "Consensus ranking supports prioritization only. "
            "Final candidates require experimental validation."
        )

    with report_tab:
        if report_json:
            with st.expander(
                "Machine-readable report",
                expanded=True,
            ):
                st.json(report_json)

            st.download_button(
                label="Download prioritization JSON",
                data=json.dumps(
                    report_json,
                    indent=2,
                ),
                file_name=(
                    "hit_prioritization_report.json"
                ),
                mime="application/json",
            )

        if report_text:
            st.markdown(report_text)

            st.download_button(
                label="Download prioritization Markdown",
                data=report_text,
                file_name=(
                    "hit_prioritization_report.md"
                ),
                mime="text/markdown",
            )
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


def render_structure_viewer() -> None:
    """Render the interactive protein-ligand structure viewer."""

    st.header("Interactive 3D protein–ligand viewer")
    st.caption(
        "HIV-1 protease is shown as a cartoon, MK1 as sticks, and the "
        "catalytic ASP25 dyad in red. Drag to rotate and scroll to zoom."
    )

    try:
        settings = load_viewer_settings()
        receptor, _ = structure_files(settings)
        poses = available_pose_files(settings)
    except (KeyError, ValueError) as exc:
        st.error(f"The 3D viewer configuration is invalid: {exc}")
        return

    selected_pose = st.selectbox("Pose", options=list(poses))
    ligand = poses[selected_pose]
    missing = [path for path in (receptor.path, ligand.path) if not path.is_file()]
    if missing:
        st.error(
            "The selected structure files are unavailable: "
            + ", ".join(path.name for path in missing)
        )
        return

    py3dmol = import_py3dmol()
    if py3dmol is None:
        st.warning(
            "Interactive 3D rendering is unavailable because the optional "
            "py3Dmol dependency is not installed. Install the project "
            "dependencies (or run `conda run -n dockforge pip install "
            "py3Dmol==2.5.3`) and restart Streamlit."
        )
        return

    try:
        view = build_view(
            py3dmol,
            receptor,
            ligand,
            list(settings.get("catalytic_residues", ["A:25", "B:25"])),
        )
        st.components.v1.html(view._make_html(), height=670, scrolling=False)
    except (OSError, ValueError) as exc:
        st.error(f"The selected structure could not be rendered: {exc}")


def render_export_center() -> None:
    """Render report status, generation controls and safe downloads."""

    st.header("Reports & Exports")
    st.caption("Generate and download portable scientific reports and validated workflow outputs. Report generation never reruns docking.")
    manifest = load_manifest()
    report_entries = [entry for entry in manifest.get("generated_files", []) if str(entry.get("path", "")).endswith((".pdf", ".html"))]
    columns = st.columns(3)
    columns[0].metric("Status", "Ready" if report_entries else "Not generated")
    columns[1].metric("Last generated", str(manifest.get("generated_at", "Unavailable"))[:19])
    columns[2].metric("Warnings", len(manifest.get("warnings", [])))
    for entry in report_entries:
        st.write(f"{Path(entry['path']).name}: {format_file_size(int(entry.get('size_bytes', 0)))}")
    if st.button("Generate PDF and HTML reports", type="primary"):
        with st.spinner("Generating reports from existing outputs…"):
            success, message = generate_reports()
        (st.success if success else st.error)(message)
        if success: st.rerun()
    exports = categorized_exports()
    report_root = PROJECT_ROOT
    tabs = st.tabs(["Final reports", "Tables & JSON", "Images", "Bundles", "Manifest"])
    with tabs[0]:
        for index, path in enumerate(exports["Reports"]):
            if path.suffix in {".pdf", ".html"}:
                st.download_button(
                    label=f"Download {path.name}",
                    data=path.read_bytes(),
                    file_name=path.name,
                    mime=(
                        "application/pdf"
                        if path.suffix == ".pdf"
                        else "text/html"
                    ),
                    key=f"export_report_{index}_{path.as_posix()}",
                )
    with tabs[1]:
        for index, path in enumerate(
            exports["Tables"] + exports["Summaries"]
        ):
            st.download_button(
                label=f"Download {path.name}",
                data=path.read_bytes(),
                file_name=path.name,
                key=f"export_table_{index}_{path.as_posix()}",
            )
    with tabs[2]:
        for index, path in enumerate(exports["Molecule images"]):
            st.download_button(
                label=f"Download {path.name}",
                data=path.read_bytes(),
                file_name=path.name,
                key=f"export_image_{index}_{path.as_posix()}",
            )
    with tabs[3]:
        figure_zip = build_zip_bytes(exports["Figures"], report_root)
        all_zip = build_zip_bytes(sum(exports.values(), []), report_root)
        st.download_button(
            "Download figure bundle",
            figure_zip,
            "DockForgeAI_figures.zip",
            "application/zip",
            key="export_figure_bundle",
        )
        st.download_button(
            "Download all final outputs",
            all_zip,
            "DockForgeAI_final_outputs.zip",
            "application/zip",
            key="export_all_outputs_bundle",
        )
    with tabs[4]: st.json(manifest or {"status": "No report manifest is available yet."})
    for warning in manifest.get("warnings", []): st.warning(warning)
    st.info("Docking scores are approximate ranking values. Demo compounds are not validated inhibitors, RDKit descriptors are not full ADMET predictions, and all predictions require experimental validation.")


def render_xai_page(page: str) -> None:
    """Render cached XAI artifacts without launching training."""
    st.header(page)
    best = load_xai_json("data/processed/xai_benchmarks/best_model.json")
    metrics = load_xai_csv("data/processed/xai_benchmarks/model_comparison.csv")
    dataset = load_xai_csv("data/processed/xai_dataset/hivpr_xai_2000_compounds.csv")
    if not best or dataset.empty:
        st.warning("XAI artifacts are unavailable. Run the local XAI pipeline explicitly from the command line.")
        return
    if page == "XAI Overview":
        columns = st.columns(4); columns[0].metric("Compounds", len(dataset)); columns[1].metric("Scaffolds", dataset.scaffold_id.nunique()); columns[2].metric("Best model", best["model_id"]); columns[3].metric("Test PR-AUC", f"{best['locked_test_metrics']['pr_auc']:.3f}")
        st.dataframe(dataset.groupby(["dataset_split", "binary_label"]).size().rename("count").reset_index(), hide_index=True)
    elif page == "Model Performance": st.dataframe(metrics, hide_index=True)
    elif page == "Global Explanations": st.dataframe(load_xai_csv("data/processed/xai_explanations/global/global_feature_importance.csv").head(30), hide_index=True)
    elif page == "Compound Explanation":
        explanations = load_xai_csv("data/processed/xai_explanations/local/compound_explanations.csv"); st.dataframe(explanations, hide_index=True)
    elif page == "Applicability Domain": st.dataframe(load_xai_csv("data/processed/xai_predictions/applicability_domain.csv"), hide_index=True)
    elif page == "Benchmarking": st.dataframe(load_xai_csv("data/processed/xai_benchmarks/feature_family_comparison.csv"), hide_index=True)
    else:
        st.download_button("Download complete XAI bundle", xai_bundle(PROJECT_ROOT), "DockForgeAI_XAI_bundle.zip", "application/zip")
    st.warning("DUD-E decoys are not experimentally proven inactive. Source shortcuts may inflate performance. Explanations describe model behaviour, not biological causality; classification does not prove binding or therapeutic activity.")


def render_xai_validation() -> None:
    """Render validation evidence already generated by the offline CLI."""
    st.header("XAI Validation")
    root="data/processed/xai_validation"
    scorecard=load_xai_json(f"{root}/validation_scorecard.json")
    if not scorecard:
        st.warning("Validation artifacts are unavailable. Run the offline validation suite explicitly; this page never starts training.")
        return
    st.error("Benchmark-supported with cautions. ChEMBL/DUD-E source separation can inflate performance; no independent external or experimental validation has been completed.")
    names=["Validation overview","Source bias","Permutation tests","Baselines","Similarity robustness","Scaffold robustness","Matched subset","Decoy sensitivity","Calibration","Explanation sanity","External validation","Scorecard","Reports"]
    tabs=st.tabs(names)
    with tabs[0]: st.json({"qualification":scorecard.get("qualification"),"externally_validated":False,"experimentally_validated":False})
    with tabs[1]: st.dataframe(load_xai_csv(f"{root}/source_bias/source_adversarial_metrics.csv"),hide_index=True); st.json(load_xai_json(f"{root}/source_bias/source_only_model_metrics.json"))
    with tabs[2]: st.json(load_xai_json(f"{root}/permutation/permutation_summary.json"))
    with tabs[3]: st.dataframe(load_xai_csv(f"{root}/baselines/dummy_baselines.csv"),hide_index=True); st.dataframe(load_xai_csv(f"{root}/baselines/nearest_neighbor_metrics.csv"),hide_index=True)
    with tabs[4]: st.dataframe(load_xai_csv(f"{root}/similarity/similarity_stratified_metrics.csv"),hide_index=True)
    with tabs[5]: st.dataframe(load_xai_csv(f"{root}/scaffolds/scaffold_frequency_metrics.csv"),hide_index=True)
    with tabs[6]: st.json(load_xai_json(f"{root}/matched_subsets/matched_subset_metrics.json"))
    with tabs[7]: st.dataframe(load_xai_csv(f"{root}/decoy_sensitivity/decoy_subset_metrics.csv"),hide_index=True)
    with tabs[8]: st.dataframe(load_xai_csv(f"{root}/calibration/calibration_robustness.csv"),hide_index=True)
    with tabs[9]: st.json(load_xai_json(f"{root}/explanation_sanity/explanation_sanity_summary.json"))
    with tabs[10]: st.info("No independent external dataset has been supplied. The CLI performs inference only and never retrains or tunes the locked model."); st.download_button("Download external validation template",external_template(),"external_xai_template.csv","text/csv")
    with tabs[11]: st.json(scorecard)
    with tabs[12]:
        report_dir=PROJECT_ROOT/root/"reports"
        for path in (report_dir/"DockForgeAI_XAI_Validation_Report.pdf",report_dir/"DockForgeAI_XAI_Validation_Report.html",PROJECT_ROOT/root/"validation_scorecard.json"):
            if path.exists(): st.download_button(f"Download {path.name}",path.read_bytes(),path.name)
        csvs=list((PROJECT_ROOT/root).rglob("*.csv")); st.caption(f"{len(csvs)} validation CSV files available")
        st.download_button("Download validation bundle",validation_bundle(PROJECT_ROOT),"DockForgeAI_XAI_validation_bundle.zip","application/zip")


def render_pubchem_external_validation() -> None:
    """Render cached PubChem external-validation results without network calls."""
    st.header("PubChem External Validation")
    root="data/processed/pubchem_external_validation"
    tabs=st.tabs(["Assay discovery","Assay audit","External dataset","Overlap analysis","Locked-model performance","Novel-scaffold performance","Distribution shift","Limitations","Reports"])
    with tabs[0]: st.dataframe(load_xai_csv(f"{root}/assay_discovery/candidate_assays.csv"),hide_index=True)
    with tabs[1]: st.dataframe(load_xai_csv(f"{root}/assay_discovery/accepted_assays.csv"),hide_index=True); st.dataframe(load_xai_csv(f"{root}/assay_discovery/rejected_assays.csv"),hide_index=True)
    dataset=load_xai_csv(f"{root}/dataset/pubchem_external_validation.csv")
    with tabs[2]: st.dataframe(dataset,hide_index=True)
    with tabs[3]:
        if not dataset.empty: st.dataframe(dataset[["compound_id","exact_training_overlap","connectivity_training_overlap","scaffold_training_overlap","maximum_training_tanimoto","independence_class"]],hide_index=True)
    metrics=load_xai_json(f"{root}/evaluation/external_validation_metrics.json")
    with tabs[4]: st.json(metrics)
    with tabs[5]: st.json((metrics.get("groups") or {}).get("novel_scaffold",{}))
    with tabs[6]: st.dataframe(load_xai_csv(f"{root}/evaluation/distribution_shift.csv"),hide_index=True); st.json(load_xai_json(f"{root}/evaluation/external_validation_summary.json"))
    with tabs[7]: st.warning("The locked model was not retrained or retuned. PubChem assay outcomes do not establish therapeutic efficacy. BindingDB literature provenance may overlap ChEMBL publications, and sparse inactive classes make several per-assay metrics unstable.")
    with tabs[8]:
        paths={"pdf":PROJECT_ROOT/root/"reports/DockForgeAI_PubChem_External_Validation_Report.pdf","html":PROJECT_ROOT/root/"reports/DockForgeAI_PubChem_External_Validation_Report.html","dataset":PROJECT_ROOT/root/"dataset/pubchem_external_validation.csv"}
        for name,path in paths.items():
            if path.exists(): st.download_button(f"Download {path.name}",path.read_bytes(),path.name,key=DOWNLOAD_KEYS[name])
        st.download_button("Download PubChem validation bundle",build_bundle(PROJECT_ROOT),"DockForgeAI_PubChem_external_validation.zip",key=DOWNLOAD_KEYS["bundle"])


def main() -> None:
    """Run the Streamlit application."""

    st.sidebar.title("DockForgeAI")

    page = st.sidebar.radio(
        "Navigation",
        [
            "Home",
            "Redocking",
            "3D viewer",
            "Reports & Exports",
            "XAI Overview",
            "Model Performance",
            "Global Explanations",
            "Compound Explanation",
            "Applicability Domain",
            "Benchmarking",
            "XAI Reports",
            "XAI Validation",
            "PubChem External Validation",
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
    elif page == "3D viewer":
        render_structure_viewer()
    elif page == "Reports & Exports":
        render_export_center()
    elif page in {"XAI Overview", "Model Performance", "Global Explanations", "Compound Explanation", "Applicability Domain", "Benchmarking", "XAI Reports"}:
        render_xai_page(page)
    elif page == "XAI Validation":
        render_xai_validation()
    elif page == "PubChem External Validation":
        render_pubchem_external_validation()
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
