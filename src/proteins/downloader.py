"""Protein structure download utilities for DockForgeAI.

This module:
1. Validates PDB identifiers.
2. Downloads protein structures from RCSB PDB.
3. Retrieves entry metadata using the RCSB Data API.
4. Saves files using a reproducible directory structure.
"""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Any

import requests
import yaml


LOGGER = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SETTINGS_PATH = PROJECT_ROOT / "config" / "settings.yaml"


class ProteinDownloadError(RuntimeError):
    """Raised when a protein structure or metadata cannot be downloaded."""


def load_settings(settings_path: Path = SETTINGS_PATH) -> dict[str, Any]:
    """Load project settings from the YAML configuration file."""

    if not settings_path.exists():
        raise FileNotFoundError(
            f"Settings file was not found: {settings_path}"
        )

    with settings_path.open("r", encoding="utf-8") as file:
        settings = yaml.safe_load(file)

    if not isinstance(settings, dict):
        raise ValueError("The settings file must contain a YAML dictionary.")

    return settings


def validate_pdb_id(pdb_id: str) -> str:
    """Validate and standardize a PDB identifier.

    A conventional PDB identifier contains exactly four alphanumeric
    characters.

    Parameters
    ----------
    pdb_id:
        User-provided PDB identifier.

    Returns
    -------
    str
        Uppercase validated PDB identifier.

    Raises
    ------
    ValueError
        If the PDB identifier is invalid.
    """

    normalized_id = pdb_id.strip().upper()

    if not re.fullmatch(r"[A-Z0-9]{4}", normalized_id):
        raise ValueError(
            "Invalid PDB ID. Enter exactly four letters or numbers, "
            "for example: 1HSG, 4G1T or 6LU7."
        )

    return normalized_id


def create_output_directories(
    settings: dict[str, Any],
) -> tuple[Path, Path]:
    """Create and return the protein and metadata output directories."""

    raw_directory = (
        PROJECT_ROOT / settings["directories"]["raw_proteins"]
    )
    metadata_directory = (
        PROJECT_ROOT / settings["directories"]["protein_metadata"]
    )

    raw_directory.mkdir(parents=True, exist_ok=True)
    metadata_directory.mkdir(parents=True, exist_ok=True)

    return raw_directory, metadata_directory


def download_structure(
    pdb_id: str,
    overwrite: bool = False,
) -> Path:
    """Download a PDB structure file from RCSB.

    Parameters
    ----------
    pdb_id:
        Four-character PDB identifier.
    overwrite:
        Download the file again when it already exists.

    Returns
    -------
    Path
        Location of the downloaded PDB file.
    """

    settings = load_settings()
    pdb_id = validate_pdb_id(pdb_id)

    raw_directory, _ = create_output_directories(settings)
    output_path = raw_directory / f"{pdb_id}.pdb"

    if output_path.exists() and output_path.stat().st_size > 0:
        if not overwrite:
            LOGGER.info("Structure already exists: %s", output_path)
            return output_path

    url_template = settings["rcsb"]["structure_url"]
    url = url_template.format(pdb_id=pdb_id)

    timeout = settings["network"]["timeout_seconds"]

    LOGGER.info("Downloading structure %s from %s", pdb_id, url)

    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
    except requests.Timeout as error:
        raise ProteinDownloadError(
            f"Download timed out for PDB ID {pdb_id}."
        ) from error
    except requests.RequestException as error:
        raise ProteinDownloadError(
            f"Could not download PDB structure {pdb_id}: {error}"
        ) from error

    content = response.text

    if not content.strip():
        raise ProteinDownloadError(
            f"RCSB returned an empty structure file for {pdb_id}."
        )

    if "<html" in content.lower():
        raise ProteinDownloadError(
            f"RCSB returned an HTML error page instead of a "
            f"structure file for {pdb_id}."
        )

    output_path.write_text(content, encoding="utf-8")

    LOGGER.info("Structure saved to %s", output_path)

    return output_path


def download_metadata(
    pdb_id: str,
    overwrite: bool = False,
) -> Path:
    """Download RCSB entry metadata and save it as JSON."""

    settings = load_settings()
    pdb_id = validate_pdb_id(pdb_id)

    _, metadata_directory = create_output_directories(settings)
    output_path = metadata_directory / f"{pdb_id}_metadata.json"

    if output_path.exists() and output_path.stat().st_size > 0:
        if not overwrite:
            LOGGER.info("Metadata already exists: %s", output_path)
            return output_path

    url_template = settings["rcsb"]["metadata_url"]
    url = url_template.format(pdb_id=pdb_id)

    timeout = settings["network"]["timeout_seconds"]

    LOGGER.info("Downloading metadata for %s", pdb_id)

    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
    except requests.Timeout as error:
        raise ProteinDownloadError(
            f"Metadata request timed out for PDB ID {pdb_id}."
        ) from error
    except requests.RequestException as error:
        raise ProteinDownloadError(
            f"Could not retrieve metadata for {pdb_id}: {error}"
        ) from error

    try:
        metadata = response.json()
    except requests.JSONDecodeError as error:
        raise ProteinDownloadError(
            f"RCSB returned invalid metadata for {pdb_id}."
        ) from error

    with output_path.open("w", encoding="utf-8") as file:
        json.dump(metadata, file, indent=2)

    LOGGER.info("Metadata saved to %s", output_path)

    return output_path


def retrieve_protein(
    pdb_id: str,
    overwrite: bool = False,
) -> dict[str, Path]:
    """Download both the protein structure and metadata."""

    structure_path = download_structure(
        pdb_id=pdb_id,
        overwrite=overwrite,
    )

    metadata_path = download_metadata(
        pdb_id=pdb_id,
        overwrite=overwrite,
    )

    return {
        "structure": structure_path,
        "metadata": metadata_path,
    }