# %%
import os
import tempfile
import numpy as np
from ase.io import read, write
from pymatgen.io.cif import CifWriter
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import lammpsdata
from structure import *
import json
import subprocess
import re

lammps = "lmp"


# %%
def modify_file(input_file, output_file, search_lines, modification_lines):
    # Open the input file for reading
    with open(input_file, 'r') as file:
        # Read the content of the input file
        content = file.readlines()

    # Modify the lines that match the search criteria
    modified_content = []
    for line in content:
        if any(search_line in line for search_line in search_lines):
            # Modify the line based on the modification lines
            modified_line = line
            for modification_line in modification_lines:
                modified_line = modified_line.replace(modification_line[0], modification_line[1])
            modified_content.append(modified_line)
        else:
            modified_content.append(line)

    # Create a temporary file to write the modified content
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
        # Write the modified content to the temporary file
        temp_file.writelines(modified_content)

    # Rename the temporary file to the output file
    os.rename(temp_file.name, output_file)

    

# %%
def lmp_energy_calculator(pot, pot_name, elm1, elm2, alloytype, elm3=""):
    """minimises the structures and calculates the energy"""
    # dict to write info to and needed directories

    Path('LMP/relaxed-structures').mkdir(parents=True, exist_ok=True)
    Path('LMP/finalDB').mkdir(parents=True, exist_ok=True)

    if os.path.isfile('LMP/finalDB/energy.json'):
        # Load the existing data from the JSON file
        with open('LMP/finalDB/energy.json', "r") as file:
            prop = json.load(file)
            #print (prop)
        os.system('rm LMP/finalDB/energy.json')
    else:
        prop = {}

    # list of available structures
    folder_path = f'LMP/{alloytype}/lammps-data'
    files = os.listdir(folder_path)
    file_names = []
    for file_name in files:
        file_names.append(file_name)
    
    # modify read_data and run the minimisation
    input_file = 'LMP/lammps-inputs/in.minimize'
    output_file = f'LMP/lammps-inputs/'

    # Specify the lines to search for in the input file
    search_lines = ['pair_style', 'pair_coeff', 'read_data', 'wd']

    # minimise and save prop to file
    for name in file_names:
        modification_lines = [
            ('pair_style', f'pair_style {pot}'),
            ('pair_coeff', f'pair_coeff * * LMP/potentials/{pot_name} {elm1} {elm2} {elm3}'),
            ("read_data", f"read_data {folder_path}/{name}"), 
            ('wd', f'write_data LMP/relaxed-structures/{name}')]
        
        modify_file(input_file, output_file+f'in.{name}_min', search_lines, modification_lines)
        
        # run the simulation and get the energy
        os.system(f"{lammps} -in LMP/lammps-inputs/in.{name}_min > /dev/null")
        command = """grep -A 1 "Energy initial" log.lammps | tail -n 1 | awk '{print $3}'"""
        energy = subprocess.check_output(command, shell=True, text=True)
        prop.update({f"{name.split('.')[0]}":[float(energy.split('\n')[0])]})
        print(f"{name} DONE")
            
    os.system("rm LMP/lammps-inputs/*_min")
    os.system("rm log.lammps")

    with open("LMP/finalDB/energy.json", "a") as file:
        # Write the dictionary to the file in JSON format
        json.dump(prop, file)


def lmp_elastic_calculator(pot, pot_name, alloytype):
    """minimises the structures and calculates the energy"""
    # dict to write info to and needed directories

    Path('LMP/finalDB').mkdir(parents=True, exist_ok=True)

    if os.path.isfile('LMP/finalDB/elastic.json'):
        # Load the existing data from the JSON file
        with open('LMP/finalDB/elastic.json', "r") as file:
            prop = json.load(file)
            #print (prop)
        os.system('rm LMP/finalDB/elastic.json')
    else:
        prop = {}

    # list of available structures
    folder_path = f'LMP/{alloytype}/lammps-data'
    files = os.listdir(folder_path)
    file_names = []
    for file_name in files:
        file_names.append(file_name)
    
    # modify read_data and run in.elastic
    input_file_init = 'LMP/lammps-inputs/elastic/init.mod'
    input_file_pot = 'LMP/lammps-inputs/elastic/potential.mod'
    output_file = f'LMP/lammps-inputs/'

    # Specify the lines to search for in the input file
    search_lines_init = ['read_data']
    search_lines_pot = ['pair_style', 'pair_coeff']

    # minimise and save prop to file
    for name in file_names:
        
        with open(f'{folder_path}/{name}') as file:
            lines = file.readlines()
            elms = lines[0]
        
        modification_init = [
            ("read_data", f"read_data {folder_path}/{name}")]
        
        modification_pot = [
            ('pair_style', f'pair_style {pot}'),
            ('pair_coeff', f'pair_coeff * * LMP/potentials/{pot_name} {elms}')]
        
        # modify init
        modify_file(input_file_init, output_file+f'init.mod', search_lines_init, modification_init)
        # modify potential
        modify_file(input_file_pot, output_file+f'potential.mod', search_lines_pot, modification_pot)

        # run the simulation and get the energy
        os.system(f"{lammps} -in LMP/lammps-inputs/elastic/in.elastic > /dev/null")
        os.system(f"rm -r LMP/lammps-inputs/init.mod LMP/lammps-inputs/potential.mod")
        command = "grep -oP '(?<=^cdvae )\[.*?\]' log.lammps"
        elastic_vector = subprocess.check_output(command, shell=True, text=True).strip().replace('\n', '')
        elastic_vector = eval(elastic_vector)
        prop.update({f"{name.split('.')[0]}":elastic_vector})
        print(f"{name} DONE")
            
    os.system("rm log.lammps restart.equil")

    with open("LMP/finalDB/elastic.json", "a") as file:
        # Write the dictionary to the file in JSON format
        json.dump(prop, file)
    
    
    

# %%
def extract_element_types(string):
    element_types = re.findall(r'[A-Z][a-z]*', string)
    return element_types

def lammps_data_to_cif(type1, type2, type3=None):
    """Converts relaxed lammps datafiles to cif files"""
    # Path to folders
    Path('LMP/cif-files').mkdir(parents=True, exist_ok=True)
    folder_path = f'LMP/relaxed-structures/'
    cif_path = f'LMP/cif-files/'
    # list of available structures
    files = os.listdir(folder_path)
    file_names = []
    for file_name in files:
        file_names.append(file_name)
    
    for name in file_names:
        # file prefix name
        d = name.split(".")[0]
        
        # Read LAMMPS data file
        atoms = lammpsdata.read_lammps_data(folder_path+name, style='atomic')

        # Manually map LAMMPS atom types to element symbols
        element_symbols = extract_element_types(f'{name}')

        atom_types = atoms.get_array('type')
        atom_symbols = [element_symbols[atom_type-1] for atom_type in atom_types]
        atoms.set_chemical_symbols(atom_symbols)

        # Convert ASE Atoms to PyMatGen Structure
        structure = AseAtomsAdaptor().get_structure(atoms)

        # Write CIF file using CifWriter
        cif_writer = CifWriter(structure)
        cif_writer.write_file(cif_path+name.split('.')[0]+'.cif')


#%%

def convert_cif_to_lammps(directory):
    # Get a list of all CIF files in the directory
    cif_files = [f for f in os.listdir(directory) if f.endswith('.cif')]

    # Convert each CIF file to LAMMPS data format
    for cif_file in cif_files:
        cif_path = os.path.join(directory, cif_file)

        # Load CIF file
        atoms = read(cif_path)

        # Save LAMMPS data file
        lammps_file = os.path.splitext(cif_path)[0] + '.data'
        write(lammps_file, atoms, format='lammps-data')

        # Update the LAMMPS data file with masses
        update_lammps_data_file(lammps_file, atoms)

        print(f"Converted {cif_file} to {lammps_file}")






def update_lammps_data_file(lammps_file, atoms):
    # Retrieve the masses and atom types from the CIF file
    masses = atoms.get_masses()
    atom_types = atoms.get_atomic_numbers()

    # Read the content of the LAMMPS data file
    with open(lammps_file, 'r') as f:
        lines = f.readlines()

    # Find the line number where the atom section starts
    atom_section_line = None
    for i, line in enumerate(lines):
        if line.strip() == "Atoms":
            atom_section_line = i
            break

    if atom_section_line is not None:
        # Insert mass information before the atom section
        mass_info = ["0 0 0 xy xz yz \n Masses \n\n"]
        unique_atom_types, indices = np.unique(atom_types, return_index=True)
        type = 1
        symbols = []
        for atom_type, index in zip(unique_atom_types, indices):
            atom_index = index + 1
            mass = masses[index]
            symbol = atoms.get_chemical_symbols()[index]
            mass_info.append(f"{type} {mass}  #{symbol}\n")
            type += 1
            symbols.append(symbol+" ")
        mass_info.append("\n")
        lines.insert(atom_section_line, "".join(mass_info))
        lines[0] = "".join(symbols)
        #print(atom_section_line)

    # Overwrite the LAMMPS data file with the updated content
    with open(lammps_file, 'w') as f:
        f.writelines(lines)
    
    return symbols 


# Specify the directory containing CIF files
directory = './LMP/generated_structures/eval_gen_cif/'

# Convert CIF files to LAMMPS data format
convert_cif_to_lammps(directory)