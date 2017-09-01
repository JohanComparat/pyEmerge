import sys
snap_name = sys.argv[1]

import MultiDark as md
box = md.MultiDarkSimulation(L_box=1000.0)

import time
t0=time.time()

box.convert_to_emerge_input_catalog_to_h5_format(snap_name)

print time.time()-t0

