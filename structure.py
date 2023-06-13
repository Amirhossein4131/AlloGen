# %%
import random
import os

from copy import deepcopy

from ase import Atoms
from ase.build import bulk
from ase.visualize import view
from ase.io import read, write


# %%
def fcc_supercell(supercell_size, atom_type, lattice_constant):
    # Create a bulk FCC structure with the given lattice constant and atom type
    bulk_fcc = bulk(atom_type, crystalstructure='fcc', a=lattice_constant, orthorhombic=True)

    # Create a supercell of the bulk FCC structure with the given size
    supercell = bulk_fcc.repeat(supercell_size)

    return supercell

# %%
def create_directory(directory):
    try:
        os.mkdir(directory)
        print(f"Directory '{directory}' created successfully.")
    except FileExistsError:
        print(f"Directory '{directory}' already exists.")

# %%
def binary_configs(matrix, substitutional_type, percentage, num_configs, output_prefix, lmpdir, dbdir):
    """ Substitutes atoms in a matrix randomly and creates as many configurations as requested. """
    
    # Create a list to store the resulting supercells
    configs = []

    # Create the lammps director
    create_directory(lmpdir)
    create_directory(dbdir)

    for i in range(num_configs):

        # a list to stotre points individually
        configs_lmp = []

        # Create a deep copy of the input supercell
        new_supercell = deepcopy(matrix)

        # Determine the number of atoms to be replaced based on the percentage
        num_to_replace = int(len(new_supercell) * percentage)

        # Select random atoms to be replaced
        atoms_to_replace = random.sample(list(new_supercell), num_to_replace)

        # Replace selected atoms with the substitutional type
        for atom in atoms_to_replace:
            atom.symbol = substitutional_type

        # Add the resulting supercell to the list
        configs.append(new_supercell)
        configs_lmp.append(new_supercell)   
        
        # Write the supercell as an XYZ/lmp file with the specified output prefix and index number
        write(f'{lmpdir}/{output_prefix}_{i}.lmp', configs_lmp, format="lammps-data")
    
    write(f'{dbdir}/{output_prefix}.xyz', configs, format="xyz")


    return configs_lmp

# %%
#Ni_matrix = fcc_supercell(2, "Ni", 3.52)

# %%
#NiAl10 = binary_configs(Ni_matrix, "Al", 0.10, 10, "NiAl10", lmpdir="lammps-data")



