import json
import os

def read_json(filename):
    with open(filename, 'r') as file:
        data = json.load(file)
    return data

def append_value_to_cif(cif_file, json_energy, json_elastic, mat, numerator, material_id, pretty_formula):
    # Read the JSON file
    energy = read_json(json_energy)[mat][0]
    ealstic = read_json(json_elastic)[mat][3]

    # Read the CIF file
    with open(cif_file, 'r') as file:
        cif_lines = file.readlines()

    # Append the value to the first line
    if cif_lines:
        cif_lines[0] = f'''{numerator},{material_id},{energy}, {ealstic},{pretty_formula},"{cif_lines[0]} '''
        cif_lines[-1] = f'''{cif_lines[-1]}" '''

    # Save the modified CIF file
    with open(cif_file, 'w') as file:
        file.writelines(cif_lines)

def prepare_cif():
    # list of available structures
    num = 1
    folder_path = f'LMP/cif-files'
    files = os.listdir(folder_path)
    file_names = []
    for file_name in files:
        file_names.append(file_name)

    for name in file_names:
        json_energy_path = './LMP/finalDB/energy.json'
        json_elastic_path = './LMP/finalDB/elastic.json'
        cif_file_path = f'./LMP/cif-files/{name}'
        mat = name.split(".")[0]
        pf = name.split("_")[0]
        append_value_to_cif(cif_file_path, json_energy_path, json_elastic_path, mat, num, mat, pf)
        num += 1



def create_db(db_type):
    prepare_cif()
    directory_path = './LMP/cif-files'
    output_file = f'./LMP/finalDB/{db_type}.csv'
    # Create the output text file and write the header line
    with open(output_file, 'w') as file:
        header = ",material_id,formation_energy_per_atom,ealstic_vector,pretty_formula,cif\n"
        file.write(header)

        # Find all .cif files in the directory
        cif_files = [file for file in os.listdir(directory_path) if file.endswith('.cif')]

        # Append the contents of each .cif file to the output text file
        for cif_file in cif_files:
            with open(os.path.join(directory_path, cif_file), 'r') as cif:
                file.write(f"{cif.read()}\n")
    
    print(f"Data appended to {output_file} successfully!")