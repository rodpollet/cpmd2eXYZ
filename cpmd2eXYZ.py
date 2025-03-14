#!/usr/bin/env python3
"""
Convert CPMD output files to extended XYZ trajectory format
"""

import os
import sys
import argparse
from contextlib import suppress

import numpy as np
from scipy import constants

# Constants for unit conversions
BOHR_TO_ANGSTROM = constants.value("atomic unit of length") * 1e10
HARTREE_TO_EV = constants.value("atomic unit of energy") / constants.elementary_charge
FORCE_CONVERSION = HARTREE_TO_EV / BOHR_TO_ANGSTROM

def read_cpmd_output(path_out):
    """Read atom types, PBC and lattice vectors from CPMD output file"""
    atom_types = []
    pbc = True
    lattice_vectors = []
    nstep = 1

    with open(path_out, 'r') as f:
        for line in f:
            if "*** ATOMS ***" in line:
                # Read atom types block
                while True:
                    next_line = next(f)
                    if "******" in next_line:
                        break
                    if "BOHR" not in next_line:
                        atom_types.append(next_line.split()[1])

            elif "ISOLATED SYSTEM CALCULATION" in line:
                # no PBC
                pbc = False

            elif "LATTICE VECTOR" in line:
                # Read lattice vectors (in bohr)
                lattice_vectors.append(line.split()[3:6])

            elif "TRAJECTORIES ARE SAVED ON FILE EVERY" in line:
                # Read number of steps for snapshots
                nstep = line.split()[6]

            elif "INITIALIZATION TIME" in line:
                break  # Stop reading after initialization

    return atom_types, pbc, np.array(lattice_vectors, dtype=np.float64), int(nstep)

def read_trajectory_data(path_trj, aref, num_atoms):
    """Read coordinate and force trajectory data"""
    coordinates = []
    if aref:
        forces = []
    
    with open(path_trj, 'r') as f:
        current_block = []
        for line in f:
            if "NEW DATA" not in line:
                current_block.append(line)
                if len(current_block) == num_atoms:
                    # Process complete block
                    if aref:
                        arr = np.loadtxt(current_block, dtype=np.float64,
                                         usecols=(1,2,3,7,8,9))
                        coords, frcs = np.split(arr, 2, axis=1)
                        forces.append(frcs)
                    else:
                        coords = np.loadtxt(current_block, dtype=np.float64,
                                            usecols=(1,2,3))
                    coordinates.append(coords)
                    current_block = []

    if aref:
        return np.array(coordinates), np.array(forces)
    else:
        return np.array(coordinates)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("CPMD_output", help="output file from a CPMD run")
    parser.add_argument("-j", "--jump", type=int, default=1,
                        help="do not keep every configuration")
    args = parser.parse_args()
    print(f"Output file is {args.CPMD_output}")
    print(f"Will keep every {args.jump} configurations")

    njump = args.jump

    path_out = args.CPMD_output
    base_dir = os.path.dirname(path_out)
    path_energy = os.path.join(base_dir, "ENERGIES")
    path_force = os.path.join(base_dir, "FTRAJECTORY")
    path_noforce = os.path.join(base_dir, "TRAJECTORY")

    if os.path.isfile(path_force):
        print("Forces and coordinates found!")
        areforces = True
        path_traj = path_force
        propert = 'species:S:1:pos:R:3:forces:R:3'
    elif os.path.isfile(path_noforce):
        print("Only coordinates found!")
        areforces = False
        path_traj = path_noforce
        propert = 'species:S:1:pos:R:3'
    else:
        sys.exit("Neither forces nor coordinates found!")

    # Read basic simulation information
    atom_types, pbc, lattice_vectors, nstep = read_cpmd_output(path_out)
    
    # Process atom type information
    atom_names, atom_type_indices, atom_counts = np.unique(
        atom_types, return_inverse=True, return_counts=True
    )
    
    # Set periodic boundary conditions
    if pbc:
        periodic = 'pbc="T T T"'
    else:
        periodic = 'pbc="F F F"'

    # Convert lattice vectors to angstrom and prepare for output
    lattice_vectors = lattice_vectors * BOHR_TO_ANGSTROM
    lattice_str = ' '.join(map(str, lattice_vectors.flatten()))
    num_atoms = sum(atom_counts)

    # Read trajectory data
    if areforces:
        coordinates, forces = read_trajectory_data(path_traj, areforces,
                                                   num_atoms)
    else:
        coordinates = read_trajectory_data(path_traj, areforces, num_atoms)

    num_frames = len(coordinates)  # size of the 1st dimension = number of frames

    # Convert units
    coordinates *= BOHR_TO_ANGSTROM
    if areforces:
        forces *= FORCE_CONVERSION

    # Read energies (every nstep to match the trajectory)
    energy_data = np.loadtxt(path_energy, dtype=np.float64, usecols=(0, 3))
    energies = energy_data[::nstep, 1] * HARTREE_TO_EV

    # Verify consistency
    if len(energies) != num_frames:
        raise ValueError("Mismatch between number of energy records and trajectory frames")

    # Write extended XYZ file
    with open("EXTTRAJ.xyz", 'w', encoding='utf-8') as outfile:
        for frame_idx in range(0,num_frames,njump):
            # Write header
            outfile.write(f"{num_atoms}\n")
            outfile.write(
                f'Lattice="{lattice_str}" Properties={propert} '
                f'energy={energies[frame_idx]:.16f} {periodic}\n'
            )

            # Write atom data
            for atom_idx in range(num_atoms):
                species = atom_types[atom_idx]
                pos = coordinates[frame_idx, atom_idx]
                pos_str   = ' '.join(f"{x:>25.16f}" for x in pos)
                if areforces:
                    force = forces[frame_idx, atom_idx]
                    force_str = ' '.join(f"{f:>25.16f}" for f in force)
                    outfile.write(f"{species:2} {pos_str} {force_str}\n")
                else:
                    outfile.write(f"{species:2} {pos_str}\n")

if __name__ == "__main__":
    with suppress(KeyboardInterrupt):
        main()
