"""DockForgeAI command-line entry point."""

from __future__ import annotations

import argparse
import logging
import sys

from src.proteins.cleaner import (
    ProteinCleaningError,
    clean_protein,
)
from src.proteins.downloader import (
    ProteinDownloadError,
    retrieve_protein,
)
from src.proteins.inspector import (
    ProteinInspectionError,
    inspect_protein_by_id,
)


def configure_logging() -> None:
    """Configure application logging."""

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def parse_comma_separated(value: str | None) -> list[str]:
    """Convert a comma-separated CLI value into a clean list."""

    if not value:
        return []

    return [
        item.strip()
        for item in value.split(",")
        if item.strip()
    ]


def build_parser() -> argparse.ArgumentParser:
    """Create the command-line argument parser."""

    parser = argparse.ArgumentParser(
        prog="DockForgeAI",
        description=(
            "Download, inspect and clean protein structures "
            "for molecular docking."
        ),
    )

    parser.add_argument(
        "--pdb-id",
        required=True,
        help="Four-character RCSB PDB identifier, such as 1HSG.",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace previously generated files.",
    )

    parser.add_argument(
        "--inspect",
        action="store_true",
        help="Inspect chains, residues, ligands, waters and metals.",
    )

    parser.add_argument(
        "--clean",
        action="store_true",
        help="Create an intermediate cleaned receptor PDB.",
    )

    parser.add_argument(
        "--chains",
        default="",
        help=(
            "Comma-separated chains to retain, for example A,B. "
            "Required when --clean is used."
        ),
    )

    parser.add_argument(
        "--extract-ligands",
        default="",
        help=(
            "Comma-separated residue names to extract, "
            "for example MK1."
        ),
    )

    parser.add_argument(
        "--retain-hetero",
        default="",
        help=(
            "Comma-separated hetero-residue names to retain "
            "in the receptor."
        ),
    )

    return parser


def print_inspection_summary(
    report: dict,
    report_path: str,
) -> None:
    """Display a concise inspection summary."""

    metadata = report.get("experimental_metadata", {})

    print("\nProtein inspection summary")
    print("-" * 55)
    print(f"PDB ID:             {report['structure_id']}")
    print(f"Models:             {report['number_of_models']}")
    print(f"Chains:             {report['number_of_chains']}")
    print(f"Protein residues:   {report['protein_residues']}")
    print(f"Total atoms:        {report['total_atoms']}")
    print(f"Water molecules:    {report['water_molecules']}")
    print(f"Ligand instances:   {report['ligand_instances']}")
    print(f"Metal instances:    {report['metal_instances']}")

    methods = metadata.get("experimental_methods") or []
    resolution = metadata.get("resolution_angstrom") or []

    print(
        "Experimental method: "
        + (
            ", ".join(methods)
            if methods
            else "Not available"
        )
    )

    print(
        "Resolution:          "
        + (
            str(resolution)
            if resolution
            else "Not available"
        )
    )

    print("\nChain summary:")

    for chain in report["chains"]:
        print(
            f"  Chain {chain['chain_id']}: "
            f"{chain['protein_residues']} residues, "
            f"{chain['atoms']} atoms, "
            f"{chain['waters']} waters"
        )

    print("\nDetected ligands/cofactors:")

    ligand_summary = report.get("ligand_summary", {})

    if ligand_summary:
        for name, count in ligand_summary.items():
            print(f"  {name}: {count}")
    else:
        print("  None")

    print("\nDetected metals:")

    metal_summary = report.get("metal_summary", {})

    if metal_summary:
        for name, count in metal_summary.items():
            print(f"  {name}: {count}")
    else:
        print("  None")

    print(f"\nFull report: {report_path}")


def print_cleaning_summary(report: dict) -> None:
    """Display protein-cleaning results."""

    print("\nProtein cleaning summary")
    print("-" * 55)
    print(f"PDB ID:               {report['pdb_id']}")
    print(
        "Selected chains:      "
        + ", ".join(report["selected_chains"])
    )
    print(
        f"Protein residues kept: "
        f"{report['protein_residues_retained']}"
    )
    print(
        f"Water molecules removed: "
        f"{report['waters_removed']}"
    )
    print(
        f"Extracted ligands:     "
        f"{len(report['extracted_ligands'])}"
    )
    print(
        f"Retained hetero-residues: "
        f"{len(report['retained_hetero_residues'])}"
    )
    print(
        f"Removed hetero-residues: "
        f"{len(report['removed_hetero_residues'])}"
    )

    print(
        f"\nClean receptor: "
        f"{report['clean_receptor_file']}"
    )

    if report["extracted_ligands"]:
        print("\nExtracted ligand files:")

        for ligand in report["extracted_ligands"]:
            print(
                f"  {ligand['name']} "
                f"chain {ligand['chain']} "
                f"residue {ligand['residue_number']}: "
                f"{ligand['output_file']}"
            )

    print(f"\nPreparation report: {report['report_file']}")


def main() -> int:
    """Run DockForgeAI."""

    configure_logging()
    parser = build_parser()
    arguments = parser.parse_args()

    try:
        files = retrieve_protein(
            pdb_id=arguments.pdb_id,
            overwrite=False,
        )

        print("\nProtein retrieval completed successfully.")
        print(f"Structure: {files['structure']}")
        print(f"Metadata:  {files['metadata']}")

        if arguments.inspect:
            inspection_report, inspection_path = (
                inspect_protein_by_id(arguments.pdb_id)
            )

            print_inspection_summary(
                report=inspection_report,
                report_path=str(inspection_path),
            )

        if arguments.clean:
            chains = parse_comma_separated(
                arguments.chains
            )

            if not chains:
                raise ValueError(
                    "--chains is required when --clean is used. "
                    "Example: --chains A,B"
                )

            cleaning_report = clean_protein(
                pdb_id=arguments.pdb_id,
                selected_chains=chains,
                extract_ligands=parse_comma_separated(
                    arguments.extract_ligands
                ),
                retain_hetero=parse_comma_separated(
                    arguments.retain_hetero
                ),
                overwrite=arguments.overwrite,
            )

            print_cleaning_summary(cleaning_report)

    except (
        ValueError,
        FileExistsError,
        FileNotFoundError,
        KeyError,
        ProteinDownloadError,
        ProteinInspectionError,
        ProteinCleaningError,
    ) as error:
        logging.error("%s", error)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())