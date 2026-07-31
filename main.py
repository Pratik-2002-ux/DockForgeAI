"""DockForgeAI command-line entry point."""

from __future__ import annotations

import argparse
import logging
import sys

from src.proteins.downloader import (
    ProteinDownloadError,
    retrieve_protein,
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
            "Download and process protein structures for "
            "molecular docking."
        ),
    )

    parser.add_argument(
        "--pdb-id",
        required=True,
        help="Four-character RCSB PDB identifier, for example 1HSG.",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace previously downloaded files.",
    )

    return parser


def main() -> int:
    """Run the DockForgeAI command-line application."""

    configure_logging()
    parser = build_parser()
    arguments = parser.parse_args()

    try:
        files = retrieve_protein(
            pdb_id=arguments.pdb_id,
            overwrite=arguments.overwrite,
        )
    except (
        ValueError,
        FileNotFoundError,
        KeyError,
        ProteinDownloadError,
    ) as error:
        logging.error("%s", error)
        return 1

    print("\nProtein retrieval completed successfully.")
    print(f"Structure: {files['structure']}")
    print(f"Metadata:  {files['metadata']}")

    return 0


if __name__ == "__main__":
    sys.exit(main())