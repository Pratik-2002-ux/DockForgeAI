"""DockForgeAI command-line entry point."""

from __future__ import annotations

import argparse
import logging
import sys

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


def build_parser() -> argparse.ArgumentParser:
    """Create the command-line argument parser."""

    parser = argparse.ArgumentParser(
        prog="DockForgeAI",
        description=(
            "Download and inspect protein structures for "
            "molecular docking."
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
        help="Replace previously downloaded files.",
    )

    parser.add_argument(
        "--inspect",
        action="store_true",
        help="Inspect chains, residues, ligands, waters and metals.",
    )

    return parser


def print_inspection_summary(
    report: dict,
    report_path: str,
) -> None:
    """Display a concise protein inspection summary."""

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
        f"Experimental method: "
        f"{', '.join(methods) if methods else 'Not available'}"
    )

    print(
        f"Resolution:          "
        f"{resolution if resolution else 'Not available'}"
    )

    print("\nChain summary:")

    for chain in report["chains"]:
        print(
            f"  Chain {chain['chain_id']}: "
            f"{chain['protein_residues']} residues, "
            f"{chain['atoms']} atoms, "
            f"{chain['waters']} waters"
        )

    ligand_summary = report.get("ligand_summary", {})

    print("\nDetected ligands/cofactors:")

    if ligand_summary:
        for name, count in ligand_summary.items():
            print(f"  {name}: {count}")
    else:
        print("  None")

    metal_summary = report.get("metal_summary", {})

    print("\nDetected metals:")

    if metal_summary:
        for name, count in metal_summary.items():
            print(f"  {name}: {count}")
    else:
        print("  None")

    print(f"\nFull report: {report_path}")


def main() -> int:
    """Run DockForgeAI."""

    configure_logging()
    parser = build_parser()
    arguments = parser.parse_args()

    try:
        files = retrieve_protein(
            pdb_id=arguments.pdb_id,
            overwrite=arguments.overwrite,
        )

        print("\nProtein retrieval completed successfully.")
        print(f"Structure: {files['structure']}")
        print(f"Metadata:  {files['metadata']}")

        if arguments.inspect:
            report, report_path = inspect_protein_by_id(
                arguments.pdb_id
            )

            print_inspection_summary(
                report=report,
                report_path=str(report_path),
            )

    except (
        ValueError,
        FileNotFoundError,
        KeyError,
        ProteinDownloadError,
        ProteinInspectionError,
    ) as error:
        logging.error("%s", error)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())