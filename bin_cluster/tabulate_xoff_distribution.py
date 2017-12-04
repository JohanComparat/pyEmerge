import numpy as n
import glob
import h5py
import os
import time
import sys
from scipy.stats import scoreatpercentile

h5_file = os.path.join(os.environ['MD10'], "cluster_h5", "hlist_1.00000.hdf5")

out_file = os.path.join(os.environ['MD10'], 'scaling_relations', 'hlist_1.00000_xoff.txt' )

f1 = h5py.File(h5_file,  "r")
redshift = f1.attrs['redshift']
xoff = f1['/halo_properties/Xoff'].value
f1.close()

n.savetxt(out_file, scoreatpercentile(xoff, [0, 25, 50, 75, 100]) )