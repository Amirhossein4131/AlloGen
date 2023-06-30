import seaborn as sns
import matplotlib.pyplot as plt
from create_db import *
import os
import re

def plot_distribution(val):
    # Create an empty list to store the values
    data_list = []
    file_names = []

    num = 1
    folder_path = f'LMP/cif-files'
    files = os.listdir(folder_path)
        
    for file_name in files:
        file_names.append(file_name)
        
    for name in file_names:
        json_energy_path = './LMP/finalDB/energy.json'
        json_elastic_path = './LMP/finalDB/elastic.json'
        mat = name.split(".")[0]
        matt = name.split("_")[0]
        #if matt == "Ni10Fe80Cr10":
        #    data_list.append(read_json(json_elastic_path)[mat][-1])
        data_list.append(read_json(json_elastic_path)[mat][2])
    

    # Plot the distribution using Seaborn
    sns.histplot(data_list, kde=True, color="orange")

    # Add labels and title to the plot
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.title(f'Distribution of {val}')

    # Show the plot
    plt.show()


#plot_distribution('C44 - All')

def list_cif_files(directory):
    cif_files = []
    for filename in os.listdir(directory):
        if filename.endswith('.cif'):
            name_without_extension = os.path.splitext(filename)[0]
            cif_files.append(name_without_extension)
    return cif_files




dire = "/home/amirhossein/Desktop/repos/AlloGen/LMP/cif-files/"
#dire = "./LMP/cif-files/"

cifs = list_cif_files(dire)

def extract_chemical_formula(cif_file_path):
    with open(cif_file_path, 'r') as file:
        for line in file:
            if line.startswith("_chemical_formula_sum"):
                match = re.search(r"'(.+)'", line)
                if match:
                    formula = match.group(1)
                    elements = re.findall(r'([A-Z][a-z]*)(\d+)?', formula)
                    element_dict = {element: int(count) if count else 1 for element, count in elements}
                    return element_dict
    return None


import pandas as pd
import plotly.express as px

structures =[]
color = []
elm1 = [] 
elm2 = []
elm3 = []

for cif in cifs:
    color.append(1)
    #print(cif)
    cif_file_path = dire + f"{cif}" + ".cif"

    element_dict = extract_chemical_formula(cif_file_path)
    #print(type(element_dict))
    if isinstance(element_dict, dict):
        structures.append(cif)
        for elm in element_dict.keys():
            if elm == 'Ni':
                elm1.append(element_dict["Ni"])
            elif elm == 'Fe':
                elm2.append(element_dict["Fe"])
            elif elm == 'Cr':
                elm3.append(element_dict["Cr"])

        if len(element_dict.keys()) == 2:
            if "Ni" not in element_dict.keys():
                elm1.append(0)
            elif "Fe" not in element_dict.keys():
                elm2.append(0)
            elif "Cr" not in element_dict.keys():
                elm3.append(0)
    else:
        pass

data = {
    "structure": structures, 
    "Ni": elm1,
    "Fe": elm2,
    "Cr": elm3,
}

df = pd.DataFrame(data, columns=['structure', 'Ni', 'Fe', 'Cr'])

print(len(elm1))
# print(df)

fig = px.scatter_ternary(df, a="Ni", b="Fe", c="Cr")
fig.show()
