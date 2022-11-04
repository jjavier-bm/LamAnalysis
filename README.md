# LamAnalysis
**Code for Simulation and Analysis of a Coarse-Grained di-Block Copolymer Lamellar Phase Nanocomposite**

[![DOI](http://img.shields.io/badge/DOI-10.3390/polym13091524-3A145B.svg)](https://doi.org/10.3390/polym13091524)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4660283.svg)](https://doi.org/10.5281/zenodo.4660283)

## About
This repository includes the code used for the research published as Burgos-Mármol, J.J.; Patti, A. Molecular Dynamics of Janus Nanodimers Dispersed in
Lamellar Phases of a Block Copolymer. *Polymers* **2021**, 13, 1524. doi: [10.3390/polym13091524](https://doi.org/10.3390/polym13091524). The content of this repository can be listed as follows:

- A ``FORTRAN`` file that generates LAMMPS scripts and LAMMPS data files.
- A ``FORTRAN`` file for structural analysis.
- A ``FORTRAN`` file for dynamical analysis.
- A ``Python`` script to remove duplicated configurations from linear trajectories and create a single long log trajectory from 100 shorter log trajectories.
- A series of LAMMPS scripts

## Cite this work
If you make use of any of the codes in this repository for your own research, make sure you include the following two references in the corresponding publications:

- <ins>The original article</ins>: Burgos-Mármol, J.J.; Patti, A. "Molecular Dynamics of Janus Nanodimers Dispersed in Lamellar Phases of a Block Copolymer". *Polymers* **2021**, 13, 1524. doi: [10.3390/polym13091524](https://doi.org/10.3390/polym13091524).
- <ins>This repository</ins>: Burgos-Mármol, J.J. LamAnalysis: Code for Simulation and Analysis of a Coarse-Grained di-Block Copolymer Lamellar Phase Nanocomposite, *Zenodo* **2021**. doi: [10.5281/zenodo.4660283](https://doi.org/10.5281/zenodo.4660283).

You can download both references from this repository in BibTeX format (`LamAnalysis.bib`).

## License
BSD 3-Clause License

Copyright (c) 2021, J. Javier Burgos-Mármol (Github username: jjavier-bm)
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
