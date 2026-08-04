"""Calculate transparent distance-based protein-ligand interaction features."""
from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import pandas as pd

from src.docking.phase4_common import atomic_json, provenance, sha256, write_checksums, write_reports

AROMATIC = {"PHE", "TYR", "TRP", "HIS"}; POSITIVE = {"LYS", "ARG", "HIS"}; NEGATIVE = {"ASP", "GLU"}; METALS = {"ZN", "MG", "MN", "FE", "CA", "CU", "CO", "NI"}


def _atoms(path: Path, ligand: bool = False) -> list[dict[str, Any]]:
    atoms: list[dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM", "HETATM")): continue
        try: x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError: continue
        if ligand:
            fields = line.split(); atom_type = fields[-1] if fields else ""; charge = float(fields[-2]) if len(fields) > 2 else 0.0
            atoms.append({"atom": line[12:16].strip(), "element": atom_type[:1].upper(), "atom_type": atom_type, "charge": charge, "x": x, "y": y, "z": z})
        else:
            atoms.append({"atom": line[12:16].strip(), "element": (line[76:78].strip() or line[12:14].strip()).upper(), "residue": line[17:20].strip(), "chain": line[21].strip(), "residue_number": line[22:26].strip(), "x": x, "y": y, "z": z})
    return atoms


def _distance(left: dict[str, Any], right: dict[str, Any]) -> float:
    return math.sqrt(sum((left[key] - right[key]) ** 2 for key in ("x", "y", "z")))


def analyze_poses(receptor: Path, docking_dir: Path, output_dir: Path, root: Path, overwrite: bool = False) -> dict[str, Any]:
    """Analyze every computed best pose using documented geometric rules."""

    if output_dir.exists() and any(output_dir.iterdir()) and not overwrite: raise FileExistsError("Analysis outputs exist; use --overwrite")
    output_dir.mkdir(parents=True, exist_ok=True)
    protein = _atoms(receptor); contacts: list[dict[str, Any]] = []; interactions: list[dict[str, Any]] = []; summaries: list[dict[str, Any]] = []
    pose_paths = sorted(docking_dir.glob("*/best_pose.pdbqt"))
    for pose in pose_paths:
        compound_id = pose.parent.name; ligand = _atoms(pose, ligand=True); counts = {key: 0 for key in ["hydrogen_bond", "hydrophobic_contact", "salt_bridge", "pi_stacking", "pi_cation", "metal_interaction"]}; contacting_residues: set[str] = set()
        for latom in ligand:
            for patom in protein:
                distance = _distance(latom, patom)
                if distance > 5.5: continue
                residue_id = f"{patom['chain']}:{patom['residue']}{patom['residue_number']}"
                if distance <= 4.0:
                    contacts.append({"compound_id": compound_id, "ligand_atom": latom["atom"], "protein_atom": patom["atom"], "protein_residue": residue_id, "distance_angstrom": round(distance, 3)}); contacting_residues.add(residue_id)
                kinds: list[str] = []
                if distance <= 3.5 and latom["element"] in {"N", "O", "S"} and patom["element"] in {"N", "O", "S"}: kinds.append("hydrogen_bond")
                if distance <= 4.0 and latom["element"] == "C" and patom["element"] == "C": kinds.append("hydrophobic_contact")
                if distance <= 4.0 and ((latom["charge"] > 0.3 and patom["residue"] in NEGATIVE) or (latom["charge"] < -0.3 and patom["residue"] in POSITIVE)): kinds.append("salt_bridge")
                if distance <= 5.5 and latom["atom_type"] == "A" and patom["residue"] in AROMATIC: kinds.append("pi_stacking")
                if distance <= 5.5 and latom["atom_type"] == "A" and patom["residue"] in POSITIVE: kinds.append("pi_cation")
                if distance <= 3.0 and patom["element"] in METALS and latom["element"] in {"N", "O", "S"}: kinds.append("metal_interaction")
                for kind in kinds:
                    counts[kind] += 1; interactions.append({"compound_id": compound_id, "interaction_type": kind, "ligand_atom": latom["atom"], "protein_atom": patom["atom"], "protein_residue": residue_id, "distance_angstrom": round(distance, 3), "method": "distance_and_atom-type_rule"})
        summaries.append({"compound_id": compound_id, "ligand_atom_count": len(ligand), "contacting_residue_count": len(contacting_residues), **{f"{key}_count": value for key, value in counts.items()}, "interaction_fingerprint": ";".join(f"{key}:{value}" for key, value in counts.items())})
    interaction_frame = pd.DataFrame(interactions, columns=["compound_id", "interaction_type", "ligand_atom", "protein_atom", "protein_residue", "distance_angstrom", "method"]); contact_frame = pd.DataFrame(contacts, columns=["compound_id", "ligand_atom", "protein_atom", "protein_residue", "distance_angstrom"]); summary_frame = pd.DataFrame(summaries)
    interaction_frame.to_csv(output_dir / "interaction_table.csv", index=False); contact_frame.to_csv(output_dir / "protein_ligand_contacts.csv", index=False); summary_frame.to_csv(output_dir / "pose_summary.csv", index=False)
    manifest = {"schema_version": 1, "module": "4.4_pose_analysis", "method": "deterministic distance and PDBQT atom-type rules", "limitations": ["Hydrogen-bond geometry does not include donor-H-acceptor angles.", "Pi interactions use atom-level aromatic proximity rather than ring-plane geometry.", "Counts are computational annotations, not experimental interactions."], "pose_count": len(pose_paths), "interaction_count": len(interaction_frame), "contact_count": len(contact_frame), "thresholds_angstrom": {"hydrogen_bond": 3.5, "hydrophobic": 4.0, "salt_bridge": 4.0, "pi": 5.5, "metal": 3.0}, **provenance(root)}
    write_reports("DockForgeAI Pose Analysis", {"Poses": len(pose_paths), "Contacts": len(contact_frame), "Interaction annotations": len(interaction_frame), "Method": manifest["method"], "Limitations": " | ".join(manifest["limitations"])}, output_dir / "pose_analysis_report.html", output_dir / "pose_analysis_report.pdf")
    outputs = [output_dir / name for name in ["interaction_table.csv", "pose_summary.csv", "protein_ligand_contacts.csv", "pose_analysis_report.html", "pose_analysis_report.pdf"]]; manifest["output_checksums"] = {path.name: sha256(path) for path in outputs}; atomic_json(manifest, output_dir / "pose_analysis_manifest.json"); write_checksums([*outputs, output_dir / "pose_analysis_manifest.json"], output_dir / "pose_analysis_checksums.json"); return manifest
