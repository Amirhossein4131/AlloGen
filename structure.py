# %%
import random
import os
from pathlib import Path

from copy import deepcopy

from ase import Atoms
from ase.build import bulk
from ase.visualize import view
from ase.io import read, write


# %%
def reset_db_creation_pipeline():
    os.system("rm -r LMP/db LMP/lammps-data")

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
        pass

# %%
def binary_configs(matrix, substitutional_type, percentage, num_configs, output_prefix, lmpdir, dbdir):
    """ Substitutes atoms in a matrix randomly and creates as many configurations as requested. """
    # Create a list to store the resulting supercells
    configs = []

    # Create the lammps director
    create_directory('LMP/'+lmpdir)
    create_directory('LMP/'+dbdir)

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
        write(f'LMP/{lmpdir}/{output_prefix}_{i}.lmp', configs_lmp, format="lammps-data")
    
    write(f'LMP/{dbdir}/{output_prefix}.xyz', configs, format="xyz")
    return configs_lmp

# %%
def binary_data_creator (elm1, elm2, celldim, lc, conf_num):
    # define matrix
    matrix = fcc_supercell(celldim, elm2, lc)

    # Create configuration strings
    config_range = []
    elm1_perc = 10
    for i in range (8):
        config_range.append(elm1+elm2+f"{elm1_perc}")
        elm1_perc += 10
    
    # create lammps data
    elm1_comp = 0.10
    for i in range (len(config_range)):
        binary_configs(matrix, elm1, elm1_comp, conf_num, "%s"%(config_range[i]) , lmpdir="lammps-data", dbdir="db")    
        elm1_comp += 0.10

# %%
def ternary_configs(matrix, type1, substitutional_type2, substitutional_type3, percentage2, percentage3, num_configs, output_prefix, lmpdir, dbdir):
    """Substitutes atoms in a matrix randomly and creates as many ternary configurations as requested."""
    # Create a list to store the resulting supercells
    configs = []

    # Create the LAMMPS and database directories if they don't exist
    Path('LMP/'+lmpdir).mkdir(parents=True, exist_ok=True)
    Path('LMP/'+dbdir).mkdir(parents=True, exist_ok=True)

    for i in range(num_configs):

        configs_lmp = []
        # Create a deep copy of the input supercell
        new_supercell = deepcopy(matrix)

        # Determine the number of atoms to be replaced based on the percentages
        num_to_replace2 = int(len(new_supercell) * percentage2)
        num_to_replace3 = int(len(new_supercell) * percentage3)

        # Select random atoms to be replaced for substitutional type 2
        atoms_to_replace2 = random.sample(list(new_supercell), num_to_replace2)

        # Replace selected atoms with substitutional type 2
        for atom in atoms_to_replace2:
            atom.symbol = substitutional_type2

        # Select random atoms to be replaced for substitutional type 3, only if the atom is of type 1
        atoms_to_replace3 = [atom for atom in new_supercell if atom.symbol == type1]
        atoms_to_replace3 = random.sample(atoms_to_replace3, num_to_replace3)

        # Replace selected atoms of type 1 with substitutional type 3
        for atom in atoms_to_replace3:
            atom.symbol = substitutional_type3

        # Add the resulting supercell to the list
        configs.append(new_supercell)
        configs_lmp.append(new_supercell)

        # Write the supercell as an XYZ/LAMMPS file with the specified output prefix and index number
        write(f'LMP/{lmpdir}/{output_prefix}_{i}.lmp', configs_lmp, format="lammps-data")

    write(f'LMP/{dbdir}/{output_prefix}.xyz', configs, format="xyz")
    return configs_lmp




def ternary_data_creator(elm1, elm2, elm3, celldim, lc, conf_num):
    """ Creatres lmp input datafiles for ternary alloys for different compositions"""
    # define matrix
    matrix = fcc_supercell(celldim, elm1, lc)

    # compositions
    e1 = [0.8, 0.1, 0.1, 0.2, 0.2, 0.2, 0.4, 0.4, 0.6]
    e2 = [0.1, 0.1, 0.8, 0.6, 0.4, 0.2, 0.4, 0.2, 0.2]
    e3 = [0.1, 0.8, 0.1, 0.2, 0.4, 0.6, 0.2, 0.4, 0.2]

    # Create configuration strings
    config_range = []
    for i in range (len(e1)):
        config_range.append(elm1+f"{int(e1[i]*100)}"+elm2+f"{int(e2[i]*100)}"+elm3+f"{int(e3[i]*100)}")

    # write the lmp data files
    for j in range (len(config_range)):
        ternary_configs(matrix, elm1, elm2, elm3, e2[j], e3[j], conf_num, f"{config_range[j]}", lmpdir="lammps-data", dbdir="db")
