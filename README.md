cpmd2eXYZ
==============================
Python tool to convert a [CPMD](https://github.com/CPMD-code) simulation to an extended XYZ trajectory file

---------
Ab initio molecular dynamics (AIMD) from CPMD
---------
An AIMD simulation performed with the CPMD code generates:
- an output file
- a file ENERGIES which includes temperature (in K) and energies (in Ha) along the trajectory
- a trajectory file TRAJEC.xyz in the XYZ format if the XYZ keyword is present in the input file
- an FTRAJECTORY file including the forces (in Ha/Bohr) if the FORCES keyword is present in the input file

> The present tool works if and only if the following line exists in the &CPMD section of the input file:
>
> TRAJECTORY XYZ FORCES SAMPLE

Extended XYZ format
----------
In comparison to the simpler XYZ format, its [extended version](https://github.com/libAtoms/extxyz) stores useful informations on the second line of each configuration (a.k.a. the comment line), such as the cell parameters, the energy or the periodic boundary conditions.

> Some machine-learning packages, e.g. [Allegro](https://github.com/mir-group/allegro), support and recommend the use of the extended XYZ format.

> The present tool will extract the required informations from the CPMD files listed above and generate a trajectory in the extended XYZ format. Coordinates will be stored in Å and forces in ev/Å.

Using the tool
----------
Simply provide the path to the CPMD output file, i.e.,

./cpmd2eXYZ.py path/to/your/cpmd_output

This will create the EXTTRAJ.xyz file.
