# Phase 5 Pilot MD Summary

## Candidate
MCULE-2304697456

## System
HIV-1 protease dimer with one ligand, TIP3P water, and approximately 0.15 M NaCl.

## Simulation
- Energy minimization
- 100 ps restrained NVT equilibration
- 100 ps restrained NPT equilibration
- 5 ns unrestrained molecular dynamics
- 37,638 atoms in the solvated system
- 501 saved frames across 5 ns

## Thermodynamic Stability
- Mean temperature remained near 300 K.
- Mean density remained near 1015 kg/m3.
- No numerical instability or actual LINCS warning was detected.

## Structural Stability
- Chain A backbone RMSD remained below approximately 0.175 nm.
- Chain B backbone RMSD remained below approximately 0.130 nm.
- Ligand RMSD remained below approximately 0.186 nm.
- The protease dimer COM distance remained near 2.22 nm.
- Protein-ligand hydrogen bonds averaged approximately 0.369 per frame.
- At least one hydrogen bond was present in approximately 35.5% of frames.
- At least one hydrogen bond was present in approximately 37.6% of frames during the final 1 ns.
- The mean minimum protein-ligand atom distance was approximately 0.192 nm, with a maximum of approximately 0.246 nm.
- Protein-ligand contacts within 0.4 nm were present in 100% of trajectory frames.
- The mean number of atom-pair contacts within 0.4 nm was approximately 292 per frame.

## Interpretation
MCULE-2304697456 remained close to its initial docked pose during the 5 ns technical pilot. The protein chains, dimer organization, and ligand pose showed low structural deviation. The result supports continued computational evaluation but does not establish experimental inhibition, binding affinity, toxicity, or biological activity.
