# MK1 redocking report

## Project and software
- Project: DockForgeAI 1.0.0
- Vina executable: /Library/anaconda3/envs/dockforge/bin/vina
- Vina version: AutoDock Vina v1.2.7
- Runtime seconds: 34.203

## Inputs
- Receptor PDBQT: /Users/pratikramchandrachaudhari/DockForgeAI/data/processed/pdbqt/receptors/1HSG_chains-A-B_receptor.pdbqt
- Ligand PDBQT: data/processed/pdbqt/ligands/1HSG_MK1.pdbqt
- Docking box: /Users/pratikramchandrachaudhari/DockForgeAI/data/processed/docking_boxes/1HSG_MK1_box.json
- Docking config: /Users/pratikramchandrachaudhari/DockForgeAI/config/docking/1HSG_MK1_redocking.txt
- RMSD threshold: 2.00 Å

## Settings
- Exhaustiveness: n/a
- Num modes: n/a
- Energy range: n/a
- CPU: n/a
- Seed: n/a

## Score and RMSD table

| pose_rank | affinity_kcal_mol | vina_rmsd_lower_bound_angstrom | vina_rmsd_upper_bound_angstrom | crystal_in_place_rmsd_angstrom | passes_2A_threshold |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 | -11.575 | 0.0 | 0.0 | 0.0 | True |
| 2 | -10.601 | 2.042 | 11.172 | 0.0 | True |
| 3 | -10.497 | 1.997 | 10.508 | 0.0 | True |
| 4 | -10.37 | 1.527 | 4.315 | 0.0 | True |
| 5 | -10.21 | 2.684 | 4.656 | 0.0 | True |
| 6 | -10.159 | 2.06 | 10.686 | 0.0 | True |
| 7 | -10.089 | 1.871 | 11.045 | 0.0 | True |
| 8 | -10.054 | 2.324 | 4.274 | 0.0 | True |
| 9 | -9.91 | 2.943 | 5.82 | 0.0 | True |
| 10 | -9.832 | 2.442 | 11.019 | 0.0 | True |

## Validation outcome
- Overall pose reproduction success: True
- Top-ranked pose success: True
- Best-energy pose rank: 1
- Best-energy pose RMSD: 0.0 Å
- Best-RMSD pose rank: 1
- Best-RMSD pose affinity: -11.575 kcal/mol
- Interpretation: Validated strongly

## Scientific limitations
- Vina scores are approximate docking scores, not experimental free energies.
- The receptor protonation assumption was documented in /Users/pratikramchandrachaudhari/DockForgeAI/data/processed/preparation_reports/1HSG_receptor_preparation.json.
- RMSD is computed in-place against the corrected crystallographic reference without prior minimization or alignment.

## Reproducibility
- Vina command: /Library/anaconda3/envs/dockforge/bin/vina --config /Users/pratikramchandrachaudhari/DockForgeAI/config/docking/1HSG_MK1_redocking.txt --out /Users/pratikramchandrachaudhari/DockForgeAI/data/processed/docking_runs/1HSG_MK1/1HSG_MK1_redocked.pdbqt.tmp
- Analysis CSV: /Users/pratikramchandrachaudhari/DockForgeAI/data/processed/docking_tables/1HSG_MK1_redocking_analysis.csv
- Summary JSON: /Users/pratikramchandrachaudhari/DockForgeAI/data/processed/docking_results/1HSG_MK1_redocking_summary.json
