import seaborn as sns
import matplotlib.pyplot as plt
from create_db import *

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


plot_distribution('C44 - All')



