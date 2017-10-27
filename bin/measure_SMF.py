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
  mass = f1['/emerge_data/stellar_mass'].value
  sel = (mass>0) & (mass!=n.inf) & (n.isnan(mass)==False)
  print( h5_file, len(mass), len(mass[sel]), len(mass[sel])>0 )

  if len(mass[sel])>0:
    counts, bb = n.histogram(n.log10(mass[sel]), bins=bins)
    dN_dVdlogM = counts*0.6777**3./(bins[1:]-bins[:-1])/volume/n.log(10)
  
    if update:
      print('updates')
      #print(f1['/stellar_mass_function/stellar_mass_low'].value)
      f1['/stellar_mass_function/stellar_mass_low'][:] = bins[:-1]
      f1['/stellar_mass_function/stellar_mass_up'][:] = bins[1:] 
      f1['/stellar_mass_function/counts'][:] = counts 
      f1['/stellar_mass_function/dN_dVdlogM'][:] = dN_dVdlogM  

    else:
      print('creates')
      stellar_mass_function_data = f1.create_group('stellar_mass_function')

      ds = stellar_mass_function_data.create_dataset('stellar_mass_low', data = bins[:-1] )
      ds.attrs['units'] = r'$M_\odot$'
      ds.attrs['long_name'] = r'$M_\odot$' 

      ds = stellar_mass_function_data.create_dataset('stellar_mass_up', data = bins[1:] )
      ds.attrs['units'] = r'$M_\odot$'
      ds.attrs['long_name'] = r'$M_\odot$' 

      ds = stellar_mass_function_data.create_dataset('dN_dVdlogM', data = dN_dVdlogM )
      ds.attrs['units'] = r'$ Mpc^{-3} dex^{-1}$'
      ds.attrs['long_name'] = r'$dN / (dV/, dlogM) $' 

      ds = stellar_mass_function_data.create_dataset('counts', data = counts )
      ds.attrs['units'] = r'count'
      ds.attrs['long_name'] = r'galaxy counts' 

  f1.close()

for h5_file in h5_files:
  #try:
  measureSMF(h5_file, update=True)
  #except( ValueError ):
  #pass
