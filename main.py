# %%
import os
from structure import *
from calculator import *
from create_db import *

# %%
reset_db_creation_pipeline()
# binary_data_creator(elm1="Ni", elm2="Fe", celldim=2, lc=3.52, conf_num=15, m1='58.69', m2='55.84', mat="NiFe")
# binary_data_creator(elm1="Ni", elm2="Cr", celldim=2, lc=3.52, conf_num=15, m1='58.69',  m2='51.99', mat="NiCr")
# binary_data_creator(elm1="Fe", elm2="Cr", celldim=2, lc=3.52, conf_num=15, m1='55.84',  m2='51.99', mat="FeCr")
# ternary_data_creator("Ni", "Fe", "Cr", celldim=3, lc=3.52, conf_num=15, m1='58.69', m2='55.84', m3='51.99')


# # %%
# lmp_energy_calculator("eam/alloy", "NiFeCr.eam.alloy", "Ni", "Fe", alloytype='Ternary', elm3='Cr')
# lmp_energy_calculator("eam/alloy", "NiFeCr.eam.alloy", elm1="Ni", elm2="Fe", alloytype='Binary/NiFe')
# lmp_energy_calculator("eam/alloy", "NiFeCr.eam.alloy", elm1="Ni", elm2="Cr", alloytype='Binary/NiCr')
# lmp_energy_calculator("eam/alloy", "NiFeCr.eam.alloy", elm1="Fe", elm2="Cr", alloytype='Binary/FeCr')

# # %%
# lammps_data_to_cif('Ni', 'Fe', 'Cr')

# create_db("train")
