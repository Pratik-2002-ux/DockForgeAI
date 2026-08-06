# Phase 6 Comparative MD Candidate Summary

## Candidate
MCULE-1302802546

## Prioritization evidence
- Consensus rank: 2
- Predicted pKi: 8.946
- Predicted Ki: 1.132 nM
- Vina score: -9.498 kcal/mol
- Applicability-domain status: in-domain
- PAINS flag: false
- Brenk flag: true
- Weighted consensus score: 0.871319

## Simulation
- Energy minimization
- 100 ps restrained NVT equilibration
- 100 ps restrained NPT equilibration
- 1 ns unrestrained molecular dynamics
- 101 saved production frames

## Thermodynamic stability
- NVT mean temperature: 299.761 K
- NPT mean temperature: 300.032 K
- NPT mean density: 1015.16 kg/m3
- No numerical instability or actual LINCS warning was detected.

## Structural stability
- Protein-backbone final RMSD: approximately 0.102 nm
- Protein-backbone maximum RMSD: approximately 0.116 nm
- Protein-backbone last-half RMSD: approximately 0.094 ± 0.008 nm
- Chain A last-half RMSD: approximately 0.087 ± 0.008 nm
- Chain B last-half RMSD: approximately 0.087 ± 0.009 nm
- Ligand final RMSD: approximately 0.136 nm
- Ligand maximum RMSD: approximately 0.178 nm
- Ligand last-half RMSD: approximately 0.130 ± 0.014 nm

## Hydrogen bonding
- Mean protein-ligand hydrogen bonds: approximately 1.366 per frame
- Maximum simultaneous hydrogen bonds: 3
- At least one hydrogen bond was present in approximately 98.0% of frames
- Last-half mean hydrogen bonds: approximately 1.627
- At least one hydrogen bond was present in approximately 98.0% of final-half frames

## Contact persistence
- Mean minimum protein-ligand distance: 0.195 nm
- Maximum minimum protein-ligand distance: 0.227 nm
- Mean contacts within 0.4 nm: 392.9 atom pairs per frame
- Minimum contacts within 0.4 nm: 300
- Maximum contacts within 0.4 nm: 499
- Frames with at least one contact: 100.0%

## Interpretation
MCULE-1302802546 remained close to its starting docked pose during the 1 ns pilot and showed highly persistent protein-ligand hydrogen bonding. The result supports continued comparative computational evaluation. The Brenk flag remains a structural-review caution and does not by itself establish toxicity or exclude the compound. These computational results do not establish experimental inhibition, binding affinity, safety, or biological activity.
