#!/usr/bin/env python3
"""
Convert CPMD output files to extended XYZ trajectory format
"""

import os
import sys
from contextlib import suppress

import numpy as np
from scipy import constants

# Constants for unit conversions
BOHR_TO_ANGSTROM = constants.value("atomic unit of length") * 1e10
HARTREE_TO_EV = constants.value("atomic unit of energy") / constants.elementary_charge
FORCE_CONVERSION = HARTREE_TO_EV / BOHR_TO_ANGSTROM

def read_cpmd_output(path_out):
    """Read atom types and lattice vectors from CPMD output file"""
    atom_types = []
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

            elif "LATTICE VECTOR" in line:
                # Read lattice vectors (in bohr)
                lattice_vectors.append(line.split()[3:6])

            elif "TRAJECTORIES ARE SAVED ON FILE EVERY" in line:
                # Read number of steps for snapshots
                nstep = line.split()[6]

            elif "INITIALIZATION TIME" in line:
                break  # Stop reading after initialization

    return atom_types, np.array(lattice_vectors, dtype=np.float64), int(nstep)

def read_trajectory_data(path_frc, num_atoms):
    """Read coordinate and force trajectory data"""
    coordinates = []
    forces = []
    
    with open(path_frc, 'r') as f:
        current_block = []
        for line in f:
            if "NEW DATA" not in line:
                current_block.append(line)
                if len(current_block) == num_atoms:
                    # Process complete block
                    arr = np.loadtxt(current_block, dtype=np.float64, usecols=(1,2,3,7,8,9))
                    coords, frcs = np.split(arr, 2, axis=1)
                    coordinates.append(coords)
                    forces.append(frcs)
                    current_block = []

    return np.array(coordinates), np.array(forces)

def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: python cpmd2eXYZ.py <CPMD_output_file>")

    path_out = sys.argv[1]
    base_dir = os.path.dirname(path_out)
    path_energy = os.path.join(base_dir, "ENERGIES")
    path_force = os.path.join(base_dir, "FTRAJECTORY")

    # Read basic simulation information
    atom_types, lattice_vectors, nstep = read_cpmd_output(path_out)
    
    # Process atom type information
    atom_names, atom_type_indices, atom_counts = np.unique(
        atom_types, return_inverse=True, return_counts=True
    )
    
    # Convert lattice vectors to angstrom and prepare for output
    lattice_vectors = lattice_vectors * BOHR_TO_ANGSTROM
    lattice_str = ' '.join(map(str, lattice_vectors.flatten()))
    num_atoms = sum(atom_counts)

    # Read trajectory data
    coordinates, forces = read_trajectory_data(path_force, num_atoms)
    num_frames = len(coordinates)

    # Convert units
    coordinates *= BOHR_TO_ANGSTROM
    forces *= FORCE_CONVERSION

    # Read energies (every nstep to match the trajectory)
    energy_data = np.loadtxt(path_energy, dtype=np.float64, usecols=(0, 3))
    energies = energy_data[::nstep, 1] * HARTREE_TO_EV

    # Verify consistency
    if len(energies) != num_frames:
        raise ValueError("Mismatch between number of energy records and trajectory frames")

    # Write extended XYZ file
    with open("EXTTRAJ.xyz", 'w', encoding='utf-8') as outfile:
        for frame_idx in range(num_frames):
            # Write header
            outfile.write(f"{num_atoms}\n")
            outfile.write(
                f'Lattice="{lattice_str}" Properties=species:S:1:pos:R:3:forces:R:3 '
                f'energy={energies[frame_idx]:.16f} pbc="T T T"\n'
            )

            # Write atom data
            for atom_idx in range(num_atoms):
                species = atom_types[atom_idx]
                pos = coordinates[frame_idx, atom_idx]
                force = forces[frame_idx, atom_idx]
                
                pos_str   = ' '.join(f"{x:>25.16f}" for x in pos)
                force_str = ' '.join(f"{f:>25.16f}" for f in force)

                outfile.write(f"{species:2} {pos_str} {force_str}\n")

if __name__ == "__main__":
    with suppress(KeyboardInterrupt):
        main()
