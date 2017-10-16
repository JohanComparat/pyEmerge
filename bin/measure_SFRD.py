import numpy as n
import glob
import h5py
import os
import time
import sys

h5_files = n.array(glob.glob(os.path.join(os.environ['MD10'], "h5", "hlist_?.?????_emerge.hdf5")))
h5_files.sort()

bins = n.arange(6,13,0.1)
xb = (bins[1:] + bins[:-1]) / 2.

def measureSMF(h5_file, volume=1000.**3., update=True):
  f1 = h5py.File(h5_file,  "r+")
  sfr = f1['/emerge_data/star_formation_rate'].value
  print( h5_file, len(sfr) )
  if len(sfr)>0:
	ok = (sfr>0)&(sfr<1000)
    sfrd = n.sum(sfr)/volume
  
    if update:
      print('updates')
      f1['/star_formation_rate_density/sfrd'][:] = sfrd

    else:
      print('creates')
      stellar_mass_function_data = f1.create_group('star_formation_rate_density')

      ds = stellar_mass_function_data.create_dataset('sfrd', data = sfrd )
      ds.attrs['units'] = r'$M_\odot h^{-3} Mpc^{3}$'
      ds.attrs['long_name'] = 'SFRD' 

  f1.close()

for h5_file in h5_files:
  try:
    measureSMF(h5_file)
  except( ValueError, KeyError ):
    pass