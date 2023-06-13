# %%
from structure import *

# %%
Ni_matrix = fcc_supercell(2, "Ni", 3.52)

# %%
NiFe_db = ["NiFe5", "NiFe10", "NiFe15", "NiFe20", "NiFe25", "NiFe30", "NiFe35", "NiFe40", "NiFe45", "NiFe50"
            , "NiFe55", "NiFe60", "NiFe65", "NiFe70", "NiFe75", "NiFe80", "NiFe85", "NiFe90", "NiFe95", "Fe"]

Fe_comp  = 0.05

for i in range (2):

    binary_configs(Ni_matrix, "Fe", Fe_comp, 10, "%s"%(NiFe_db[i]) , lmpdir="lammps-data", dbdir="db")    
    Fe_comp += 0.05

# %%
NiCr_db = ["Ni", "NiCr5", "NiCr10", "NiCr15", "NiCr20", "NiCr25", "NiCr30", "NiCr35", "NiCr40", "NiCr45", "NiCr50"
            , "NiCr55", "NiCr60", "NiCr65", "NiCr70", "NiCr75", "NiCr80", "NiCr85", "NiCr90", "NiCr95", "Cr"]

Cr_comp  = 0

for i in range (len(NiCr_db)):

    binary_configs(Ni_matrix, "Cr", Cr_comp, 10, "%s"%(NiCr_db[i]) , lmpdir="lammps-data", dbdir="db")    
    Cr_comp += 0.05

# %%
Fe_matrix = fcc_supercell(2, "Fe", 3.52)

FeCr_db = ["FeCr5", "FeCr10", "FeCr15", "FeCr20", "FeCr25", "FeCr30", "FeCr35", "FeCr40", "FeCr45", "FeCr50"
            , "FeCr55", "FeCr60", "FeCr65", "FeCr70", "FeCr75", "FeCr80", "FeCr85", "FeCr90", "FeCr95"]

Cr_comp  = 0.05

for i in range (len(FeCr_db)):

    binary_configs(Fe_matrix, "Cr", Cr_comp, 10, "%s"%(FeCr_db[i]) , lmpdir="lammps-data", dbdir="db")    
    Cr_comp += 0.05


