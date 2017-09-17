import sys
snap_name = sys.argv[1]
aexp      = sys.argv[2]
redshift  = sys.argv[3]
age_yr    = sys.argv[4]
rho_crit  = sys.argv[5]
delta_vir = sys.argv[6]


import MultiDark as md
box = md.MultiDarkSimulation(L_box=1000.0)

import time
t0=time.time()

box.convert_to_emerge_input_catalog_to_h5_format_light(snap_name, aexp, redshift, age_yr, rho_crit, delta_vir)
#box.convert_to_emerge_input_catalog_to_h5_format_ultralight(snap_name, aexp, redshift, age_yr, rho_crit, delta_vir)

print time.time()-t0

