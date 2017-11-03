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

h5_dir = os.path.join(os.environ[env], 'h5' )
input_list = n.array(glob.glob(os.path.join(h5_dir, "hlist_?.?????_emerge.hdf5")))
input_list.sort()

file_1 = input_list[ii]
print(file_1)
f1 = h5py.File(file_1,  "r+")
z = f1.attrs['redshift']

from scipy.stats import norm
# original
# meanSM= lambda Mh, z : n.log10(Mh * 2. * 2. * ( 0.0351 - 0.0247 * z/(1.+z)) / ((Mh/ (10**(11.59 + 1.195 * z/(1.+z))) )**(- 1.376 + 0.826 * z/(1.+z)) + ( Mh /(10**(11.59 + 1.195 * z/(1.+z))) )**(0.608 + 0.329 *z/(1.+z)) ) )

meanSM= lambda Mh, z : n.log10(Mh * 2. * ( 0.0351 - 0.0247 * z/(1.+z)) / ((Mh/ (10**(11.79 + 1.5 * z/(1.+z))) )**(- 0.9 + 0.5  * z/(1.+z)) + ( Mh /(10**(11.79 + 1.5 * z/(1.+z))) )**(0.67 + 0.2 * z/(1.+z)) ) )

mean_SM = meanSM(f1['/halo_properties/mvir'].value/0.6777, z)

fun = lambda mmm : norm.rvs( loc = mmm, scale = 0.15 )

Mgal_mvir_Mo13 = fun(mean_SM) # n.array(pool.starmap( fun, mean_SM ))

#print( "res  mgal", Mgal_mvir_Mo13)
#print( "diff mgal - mvir", n.mean(mean_SM-Mgal_mvir_Mo13) )
#print( "mean, std magl - mh",n.mean(mean_SM-Mgal_mvir_Mo13), n.std(mean_SM-Mgal_mvir_Mo13))

if status=='create':
  halo_data = f1.create_group('moster_2013_data')

  ds = halo_data.create_dataset('stellar_mass', data = 10**Mgal_mvir_Mo13 )
  ds.attrs['units'] = '\log_{10}(M/M_\odot)'
  ds.attrs['long_name'] = 'log10 of the stellar mass' 

if status=='update':
  f1['/moster_2013_data/stellar_mass'][:] = 10**Mgal_mvir_Mo13

f1.close()


