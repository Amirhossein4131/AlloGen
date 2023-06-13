# %%
import os
import tempfile
from ase.io import read, write
from pymatgen.io.cif import CifWriter
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import lammpsdata


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
def create_tmp_min(material, pot, pot_name, elm1, elm2, mass1, mass2, mass3):
    # Define input and output file paths
    input_file = 'lammps-inputs/in.minimize'
    output_file = f'lammps-inputs/in.{material}_min_temp'

    # Specify the lines to search for in the input file
    search_lines = ['pair_style', 'pair_coeff', 'm1', 'm2', 'm3']

    # Specify the lines to modify in the input file
    modification_lines = [
        ('pair_style', f'pair_style {pot}'),
        ('pair_coeff', f'pair_coeff * * ../potentials/{pot_name} {elm1} {elm2}'),
        ('m1', f'{mass1}'),
        ('m2', f'{mass2}'),
        ('m3', f'{mass3}')
    ]

    # Call the function to modify the input file
    modify_file(input_file, output_file, search_lines, modification_lines)


# %%
def lmp_calculator():
    os.system(f"{lammps} -in lammps-inputs/in.NiFe_min_temp")

lmp_calculator()

# %%

def lammps_data_to_cif(data_file, cif_file, type1, type2):
    # Read LAMMPS data file
    atoms = lammpsdata.read_lammps_data(data_file, style='atomic')

    # Manually map LAMMPS atom types to element symbols
    element_symbols = [f'{type1}', f'{type2}']
    atom_types = atoms.get_array('type')
    atom_symbols = [element_symbols[atom_type-1] for atom_type in atom_types]
    atoms.set_chemical_symbols(atom_symbols)

    # Convert ASE Atoms to PyMatGen Structure
    structure = AseAtomsAdaptor().get_structure(atoms)

    # Write CIF file using CifWriter
    cif_writer = CifWriter(structure)
    cif_writer.write_file(cif_file)

lammps_data_to_cif("./lmp-data-files/NiFe10_0.lmp", "NiFe_2.cif")


# %%
create_tmp_min("NiFe", "eam/alloy", "NiFeCr.eam.alloy", "Ni", "Fe", "mass 1 2", "mass 2 3", " ")


