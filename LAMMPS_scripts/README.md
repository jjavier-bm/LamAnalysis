# LamAnalysis
Code for Simulation and Analysis of a Coarse-Grained di-Block Copolymer Lamellar Phase Nanocomposite

## LAMMPS scripts
This folder contains the minimal scripts needed to reproduce what was done for the article.

The scripts are separated into two folders: `pure`, for the pristine polymer systems, and `pnc`, as in polymer nanocomposites for the systems containing nanodimers.
Only one of the particle mass fractions is represented in the `pnc` script names, as the scripts are all the same regardless of the number of particles.

The subfolders include the following scripts (listed by the order they were run):

- Pure polymer:
  1. **Initial**: Thermalisation at an estimated pressure.
  2. **Initial2**: Expansion of systems (lowering pressure).
  3. **Fast and slow Cycles**: Two independent runs to produce the Pressure vs Density curve and determine the pressure for the desired density value.
  4. **Equilibration**: Equilibration of the systems at desired pressure and density. Sometimes this run was not enough for systems to reach equilibrium and extra simulation time was assigned for it.
  5. **Production**: Generation of trajectories (raw data) for subsequent analysis.

- Polymer nanocomposite (`wpcx`, with x>0, meaning mass fraction = x wt%):
  1. **Initial**: Thermalisation and Pressurisation of the systems (desired pressure previously determined for pure polymer).
  2. **Equilibration**: Equilibration of the systems at desired pressure and density. Sometimes this run was not enough for systems to reach equilibrium and extra simulation time was assigned for it.
  3. **Production**: Generation of trajectories (raw data) for subsequent analysis.

  Note: In the case of systems containing nanodimers, the first two steps were run for neutral particles only, while for the other particles equilibration started from the coordinates of neutral particles after step 2.
