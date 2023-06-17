# %%
import os
import tempfile
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


