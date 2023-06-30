# %%
import random
import os
from pathlib import Path

from copy import deepcopy

from ase import Atoms
from ase.build import bulk
from ase.visualize import view
from ase.io import read, write


# %%/
def reset_db_creation_pipeline():
    os.system("rm -r LMP/Binary LMP/Ternary LMP/relaxed-structures LMP/finalDB LMP/cif-files")

def replace_line(file_name, line_num, text):
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()

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
def binary_configs(matrix, substitutional_type, elm2, percentage,
                    num_configs, output_prefix, lmpdir, dbdir, m1, m2):
    """ Substitutes atoms in a matrix randomly and creates as many configurations as requested. """
    # Create a list to store the resulting supercells
    configs = []

    # Create the lammps director
    create_directory('LMP/Binary')
    Path('LMP/Binary/'+lmpdir).mkdir(parents=True, exist_ok=True)
    Path('LMP/Binary/'+dbdir).mkdir(parents=True, exist_ok=True)

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
        write(f'LMP/Binary/{lmpdir}/{output_prefix}_{i}.lmp', configs_lmp, format="lammps-data")
        replace_line(f'LMP/Binary/{lmpdir}/{output_prefix}_{i}.lmp', 8, f'\nMasses\n\n1 {m1}\n2 {m2}\n')
        replace_line(f'LMP/Binary/{lmpdir}/{output_prefix}_{i}.lmp', 7, '0.0 0.0 0.0 xy xz yz')
        replace_line(f'LMP/Binary/{lmpdir}/{output_prefix}_{i}.lmp', 0, f'{substitutional_type} {elm2}')

            
    write(f'LMP/Binary/{dbdir}/{output_prefix}.xyz', configs, format="xyz")
    return configs_lmp

# %%
def binary_data_creator (elm1, elm2, celldim, lc, conf_num, m1, m2, mat):
    # define matrix
    matrix = fcc_supercell(celldim, elm2, lc)

    # Create configuration strings
    config_range = []
    # elm1_perc = 5
    # for i in range (18):
    #     config_range.append(elm1+elm2+f"{elm1_perc}")
    #     elm1_perc += 5
    # print(config_range)

    for i in range (50):
        c1 = round(random.uniform(0.02, 0.98), 2)
        config_range.append(c1)
       

    # create lammps data

    for i in range (len(config_range)):
        c1 = int(config_range[i]*100)
        c2 = int(100 - c1)
        binary_configs(matrix, elm1, elm2, config_range[i], conf_num,
                        elm1+f"{c1}"+elm2+f"{c2}", m1=m1, m2=m2, lmpdir=f"{mat}/lammps-data", dbdir=f"{mat}/db")    
        

# %%
def ternary_configs(matrix, type1, substitutional_type2, substitutional_type3,
                     percentage2, percentage3, num_configs, output_prefix, lmpdir, dbdir, m1, m2, m3):
    """Substitutes atoms in a matrix randomly and creates as many ternary configurations as requested."""
    # Create a list to store the resulting supercells
    configs = []

    # Create the LAMMPS and database directories if they don't exist
    Path('LMP/Ternary').mkdir(parents=True, exist_ok=True)
    Path('LMP/Ternary/'+lmpdir).mkdir(parents=True, exist_ok=True)
    Path('LMP/Ternary/'+dbdir).mkdir(parents=True, exist_ok=True)

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
        write(f'LMP/Ternary/{lmpdir}/{output_prefix}_{i}.lmp', configs_lmp, format="lammps-data")
        replace_line(f'LMP/Ternary/{lmpdir}/{output_prefix}_{i}.lmp', 8, f'\nMasses\n\n1 {m1}\n2 {m2}\n3 {m3}\n')
        replace_line(f'LMP/Ternary/{lmpdir}/{output_prefix}_{i}.lmp', 7, '0.0 0.0 0.0 xy xz yz')
        replace_line(f'LMP/Ternary/{lmpdir}/{output_prefix}_{i}.lmp', 0, f'{type1} {substitutional_type2} {substitutional_type3}')

    write(f'LMP/Ternary/{dbdir}/{output_prefix}.xyz', configs, format="xyz")
    return configs_lmp




def ternary_data_creator(elm1, elm2, elm3, celldim, lc, conf_num, m1, m2, m3):
    """ Creatres lmp input datafiles for ternary alloys for different compositions"""
    # define matrix
    matrix = fcc_supercell(celldim, elm1, lc)

    # compositions
    #e1 = [0.8, 0.1, 0.1, 0.2, 0.2, 0.2, 0.4, 0.4, 0.6, 0.82, 0.85, 0.90, 0.12, 0.08, 0.13, 0.11, 0.13, 0.09, 0.23, 0.25, 0.18, 0.22, 0.23, 0.18, 0.17, 0.21, 0.24, 0.44, 0.38, 0.45, 0.42, 0.36, 0.45, 0.62, 0.65, 0.58]
    #e2 = [0.1, 0.1, 0.8, 0.6, 0.4, 0.2, 0.4, 0.2, 0.2, 0.08, 0.05, 0.07, 0.07, 0.09, 0.12, 0.84, 0.83, 0.76, 0.66, 0.58, 0.63, 0.44, 0.38, 0.46, 0.21, 0.17, 0.23, 0.45, 0.44, 0.35, 0.21, 0.24, 0.18, 0.18, 0.16, 0.24]
    #e3 = [0.1, 0.8, 0.1, 0.2, 0.4, 0.6, 0.2, 0.4, 0.2, 0.10, 0.10, 0.03, 0.81, 0.83, 0.75, 0.05, 0.04, 0.15, 0.11, 0.17, 0.19, 0.36, 0.39, 0.36, 0.62, 0.62, 0.53, 0.11, 0.18, 0.20, 0.37, 0.40, 0.37, 0.20, 0.09, 0.18]
    
    e1 = []
    e2 = []
    e3 = []
    
    for i in range (350):
        c1 = random.uniform(0.02, 0.98)
        c2 = random.uniform(0.02, 1-c1)
        c3 = 1-c1-c2
    
        e1.append(c1)
        e2.append(c2)
        e3.append(c3)

    # Create configuration strings
    config_range = []
    for i in range (len(e1)):
        config_range.append(elm1+f"{int(e1[i]*100)}"+elm2+f"{int(e2[i]*100)}"+elm3+f"{int(e3[i]*100)}")

    print(config_range)


    # write the lmp data files
    for j in range (len(e1)):
        ternary_configs(matrix, elm1, elm2, elm3, e2[j], e1[j], conf_num,
                         f"{config_range[j]}", m1=m1, m2=m2, m3=m3, lmpdir="lammps-data", dbdir="db")
