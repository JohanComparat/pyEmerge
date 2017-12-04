"""
Routine to add Moster et al. 2013 stellar masses

python3 lc_add_Ms_Mo13.py 115 MD10
"""
import sys
ii = int(sys.argv[1])
env = sys.argv[2] # 'MD10'
status = sys.argv[3]

import h5py    # HDF5 support
import os
import glob
import numpy as n

h5_dir = os.path.join(os.environ[env], 'cluster_h5/' )
input_list = n.array(glob.glob(os.path.join(h5_dir, "hlist_?.?????.hdf5")))
input_list.sort()

xoff_000, xoff_025, xoff_050, xoff_075, xoff_100 = n.loadtxt( os.path.join(os.environ['MD10'], 'scaling_relations', 'hlist_1.00000_xoff.txt' ) )

file_1 = input_list[ii]
print(file_1)
f1 = h5py.File(file_1,  "r+")

xoff = f1['/halo_properties/Xoff'].value

cool_class = n.ones_like(xoff)
cool_class[cool_class>xoff_025] = cool_class[cool_class>xoff_025]+1
cool_class[cool_class>xoff_050] = cool_class[cool_class>xoff_050]+1
cool_class[cool_class>xoff_075] = cool_class[cool_class>xoff_075]+1

if status=='create':
  halo_data = f1.create_group('cluster_data')

  ds = halo_data.create_dataset('cool_class', data = cool_class.astype('int') )
  ds.attrs['units'] = 'coolness class'
  ds.attrs['long_name'] = '1: very relaxed cool core, 2: relaxed cool core, 3: disturbed non cool core, 4: very disturbed non cool core' 

if status=='update':
  f1['/cluster_data/cool_class'][:] = cool_class.astype('int')

f1.close()


