import sys
import os
import EmergeMultiDark as md
box = md.MultiDarkSimulation(L_box=1000.0)
import numpy as n
import time
t0=time.time()
import glob 

path_2_snaps = n.array(glob.glob(os.path.join(os.environ['MD10'], 'cluster_data', 'hlist*.data')))
path_2_snaps.sort()
path_2_h5_files = n.array([os.path.join(os.environ['MD10'], 'cluster_h5', os.path.basename(path_2_snap[:-5]) + '.hdf5') for path_2_snap in path_2_snaps ])
aexps = n.array([ float(os.path.basename(path_2_snap)[6:-5]) for path_2_snap in path_2_snaps ])

for path_snap, path_h5, aexp in zip(path_2_snaps, path_2_h5_files, aexps):
	print(path_snap, path_h5, aexp)
	box.convert_to_h5_format_cluster(path_snap, path_h5, aexp, 1./aexp-1.)
	print(time.time()-t0)

