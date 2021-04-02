# LamAnalysis
Code for Simulation and Analysis of a Coarse-Grained di-Block Copolymer Lamellar Phase Nanocomposite

## LAMMPS scripts
This folder contains the minimal scripts needed to reproduce what was done for the article.

The scripts are separated into two folders: "pure", for the pristine polymer systems, and "pnc", as in polymer nanocomposites for the systems containing nanodimers.
Only one of the particle mass fractions is represented in the script names, as the scripts are all the same regardless of the number of particles.

The subfolders include the following scripts sorted by the order they were run:

1. Initial: Thermalisation of the systems.
2. Initial2: Pressurisation of the systems.
3. Equilibration: Equilibration of the systems. Sometimes this run was not enough for systems to reach equilibrium and extra simulation time was assigned for it.
4. Production: Generation of trajectories (raw data) for subsequent analysis.

In the case of systems containing nanodimers, the first two steps were run for neutral particles only, while for the other particles equilibration started from the coordinates of neutral particles after the "Initial2" step.

For pristine polymers an expansion-contraction cycle was run between steps 2 and 3 to determine the phase behaviour and the Pressure for equilibration and production runs of all systems.
