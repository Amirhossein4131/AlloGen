# %%
import os
from structure import *
from calculator import *

# %%
#binary_data_creator("Ni", "Fe", celldim=2, lc=3.52, conf_num=2)
#ternary_data_creator("Ni", "Fe", "Cr", celldim=3, lc=3.52, conf_num=2)

# %%
#reset_db_creation_pipeline()

# %%
#create_tmp_min("eam/alloy", "NiFeCr.eam.alloy", "Ni", "Fe", "mass 1 2", "mass 2 3", " ")
#lmp_calculator()

# %%
lmp_energy_calculator("eam/alloy", "NiFeCr.eam.alloy", "Ni", "Fe", " ", "mass 1 2", "mass 2 3", " ")