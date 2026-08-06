# Phase 6 Comparative MD Candidate Summary

## Candidate
MCULE-4176465450

## Prioritization evidence
- Consensus rank: 4
- Predicted pKi: 9.562
- Predicted Ki: 0.274 nM
- Vina score: -7.644 kcal/mol
- Applicability-domain status: in-domain
- PAINS flag: false
- Brenk flag: false
- Weighted consensus score: 0.828162

## Simulation
- Energy minimization
- 100 ps restrained NVT equilibration
- 100 ps restrained NPT equilibration
- 1 ns unrestrained molecular dynamics
- 101 saved production frames

## Thermodynamic stability
- NVT mean temperature: 299.941 K
- NPT mean temperature: 300.065 K
- NPT mean density: 1015.19 kg/m3
- Production mean temperature: 300.024 K
- Production mean pressure: 0.888 bar
- Production mean density: 1014.89 kg/m3
- No numerical instability or actual LINCS warning was detected.

## Structural stability
- Chain A final backbone RMSD: approximately 0.111 nm
- Chain A maximum backbone RMSD: approximately 0.138 nm
- Chain A last-half RMSD: approximately 0.108 ± 0.009 nm
- Chain B final backbone RMSD: approximately 0.090 nm
- Chain B maximum backbone RMSD: approximately 0.123 nm
- Chain B last-half RMSD: approximately 0.096 ± 0.008 nm
- Ligand final RMSD: approximately 0.109 nm
- Ligand maximum RMSD: approximately 0.137 nm
- Ligand last-half RMSD: approximately 0.117 ± 0.012 nm
- The whole-dimer RMSD was excluded because periodic-boundary imaging produced a confirmed artifact.

## Contact persistence
- Mean minimum protein-ligand distance: approximately 0.202 nm
- Maximum minimum protein-ligand distance: approximately 0.227 nm
- Last-half mean minimum distance: approximately 0.203 nm
- Mean contacts within 0.4 nm: approximately 263 atom pairs per frame
- Minimum contacts within 0.4 nm: 79
- Maximum contacts within 0.4 nm: 408
- Protein-ligand contacts were present in 100% of frames
- Last-half mean contacts: approximately 252 atom pairs per frame

## Hydrogen bonding
- Mean protein-ligand hydrogen bonds: approximately 1.020 per frame
- Maximum simultaneous hydrogen bonds: 3
- At least one hydrogen bond was present in approximately 98.0% of frames
- Last-half mean hydrogen bonds: approximately 1.000
- At least one hydrogen bond was present in approximately 98.0% of final-half frames

## Interpretation
MCULE-4176465450 showed the lowest ligand RMSD among the three 1 ns comparative candidates and retained close protein contacts throughout the simulation. Hydrogen bonding was highly persistent, although the mean number of hydrogen bonds was lower than for MCULE-5294193665. The compound remained structurally stable during this short technical pilot. These computational results do not establish experimental inhibition, binding affinity, safety, or biological activity.
