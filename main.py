# %%
import os

from structure import *
from calculator import *
from create_db import *

# %%
reset_db_creation_pipeline()
binary_data_creator(elm1="Ni", elm2="Fe", celldim=3, lc=3.52, conf_num=1, m1='58.69', m2='55.84', mat="NiFe")
binary_data_creator(elm1="Ni", elm2="Cr", celldim=3, lc=3.52, conf_num=20, m1='58.69',  m2='51.99', mat="NiCr")
binary_data_creator(elm1="Fe", elm2="Cr", celldim=3, lc=3.52, conf_num=20, m1='55.84',  m2='51.99', mat="FeCr")
ternary_data_creator("Ni", "Fe", "Cr", celldim=3, lc=3.52, conf_num=20, m1='58.69', m2='55.84', m3='51.99')


# # # # %%
lmp_energy_calculator("eam/alloy", "NiFeCr.eam.alloy", alloytype='Ternary')
lmp_energy_calculator("eam/alloy", "NiFeCr.eam.alloy", alloytype='Binary/NiFe')
lmp_energy_calculator("eam/alloy", "NiFeCr.eam.alloy", alloytype='Binary/NiCr')
lmp_energy_calculator("eam/alloy", "NiFeCr.eam.alloy", alloytype='Binary/FeCr')

lmp_elastic_calculator("eam/alloy", "NiFeCr.eam.alloy", alloytype='Ternary')
lmp_elastic_calculator("eam/alloy", "NiFeCr.eam.alloy", alloytype='Binary/NiFe')
lmp_elastic_calculator("eam/alloy", "NiFeCr.eam.alloy", alloytype='Binary/NiCr')
lmp_elastic_calculator("eam/alloy", "NiFeCr.eam.alloy", alloytype='Binary/FeCr')

# # # %%
lammps_data_to_cif('LMP/')

create_db("train")

retrain_property_calculator("LMP/opt", "opt/")
