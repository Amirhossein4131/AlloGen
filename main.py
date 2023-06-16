# %%
import os
from structure import *
from calculator import *

# %%
reset_db_creation_pipeline()
#binary_data_creator("Ni", "Fe", celldim=2, lc=3.52, conf_num=2, m1='20', m2='50')
#ternary_data_creator("Ni", "Fe", "Cr", celldim=3, lc=3.52, conf_num=2, m1='58.69', m2='55.84', m3='51.99')

# %%
#create_tmp_min("LMP/eam/alloy", "NiFeCr.eam.alloy", "Ni", "Fe", "mass 1 2", "mass 2 3", " ")
#lmp_calculator()

# %%
#lmp_energy_calculator("eam/alloy", "NiFeCr.eam.alloy", "Ni", "Fe", alloytype='Ternary', elm3='Cr')
#lmp_energy_calculator("eam/alloy", "NiFeCr.eam.alloy", "Ni", "Fe", alloytype='Binary')

# %%
#lammps_data_to_cif('Ni', 'Fe', 'Cr')
#lammps_data_to_cif('Ni', 'Fe')