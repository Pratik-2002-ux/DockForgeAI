# Phase 6 Comparative MD Candidate Summary

## Candidate
MCULE-5294193665

## Prioritization evidence
- Consensus rank: 3
- Predicted pKi: 8.589
- Predicted Ki: 2.575 nM
- Vina score: -8.972 kcal/mol
- Applicability-domain status: in-domain
- PAINS flag: false
- Brenk flag: true
- Weighted consensus score: 0.833716

## Simulation
- Energy minimization
- 100 ps restrained NVT equilibration
- 100 ps restrained NPT equilibration
- 1 ns unrestrained molecular dynamics
- 101 saved production frames

## Thermodynamic stability
- NVT mean temperature: 299.639 K
- NPT mean temperature: 299.953 K
- NPT mean density: 1015.7 kg/m3
- Production mean temperature: 300.037 K
- Production mean pressure: 2.047 bar
- Production mean density: 1014.95 kg/m3
- No numerical instability or actual LINCS warning was detected.

## Structural stability
- Chain A final backbone RMSD: approximately 0.094 nm
- Chain A maximum backbone RMSD: approximately 0.127 nm
- Chain A last-half RMSD: approximately 0.106 ± 0.007 nm
- Chain B final backbone RMSD: approximately 0.100 nm
- Chain B maximum backbone RMSD: approximately 0.147 nm
- Chain B last-half RMSD: approximately 0.103 ± 0.009 nm
- Ligand final RMSD: approximately 0.206 nm
- Ligand maximum RMSD: approximately 0.206 nm
- Ligand last-half RMSD: approximately 0.160 ± 0.018 nm
- Mean chain COM distance: approximately 2.334 nm
- Chain COM distance remained between approximately 2.210 and 2.396 nm.
- The whole-dimer RMSD was excluded because periodic-boundary imaging produced a confirmed artifact.

## Contact persistence
- Mean minimum protein-ligand distance: approximately 0.195 nm
- Maximum minimum protein-ligand distance: approximately 0.234 nm
- Last-half mean minimum distance: approximately 0.198 nm
- Mean contacts within 0.4 nm: approximately 293 atom pairs per frame
- Minimum contacts within 0.4 nm: 66
- Maximum contacts within 0.4 nm: 555
- Protein-ligand contacts were present in 100% of frames
- Last-half mean contacts: approximately 170 atom pairs per frame

## Hydrogen bonding
- Mean protein-ligand hydrogen bonds: approximately 1.921 per frame
- Maximum simultaneous hydrogen bonds: 5
- At least one hydrogen bond was present in approximately 97.0% of frames
- Last-half mean hydrogen bonds: approximately 1.804
- At least one hydrogen bond was present in approximately 94.1% of final-half frames

## Interpretation
MCULE-5294193665 remained associated with the HIV-1 protease binding region during the 1 ns pilot. Both protease chains remained individually stable, the ligand retained close protein contacts in every frame, and hydrogen bonding was highly persistent. The ligand showed slightly greater positional deviation than MCULE-1302802546. The Brenk flag remains a structural-review caution. These computational results do not establish experimental inhibition, binding affinity, safety, or biological activity.
