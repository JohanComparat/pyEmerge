import sys

path_2_snap = sys.argv[1]
path_2_output = sys.argv[2]

import MultiDark as md
box = md.MultiDarkSimulation(L_box=1000.0)

import time
t0=time.time()

box.transform_rockstar_hlist_catalog_into_emerge_input_catalog(path_2_snap, path_2_output, mmin=1e13, option='complete')

print time.time()-t0

