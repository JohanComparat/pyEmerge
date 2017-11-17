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
hh = 0.6777

def measureSMF(h5_file, volume=1000.**3./hh**3., update=True):
  f1 = h5py.File(h5_file,  "r+")
  mass = f1['/moster_2013_data/stellar_mass'].value 
  print(f1['/agn_properties/agn_activity'].value)
  sel = (mass>0) & (mass!=n.inf) & (n.isnan(mass)==False)&(f1['/agn_properties/agn_activity'].value==1)
  print( h5_file, len(mass), len(mass[sel]), len(mass[sel])>0 )

  if len(mass[sel])>0:
    counts, bb = n.histogram(n.log10(mass[sel]), bins=bins)
    dN_dVdlogM = counts/(bins[1:]-bins[:-1])/volume/n.log(10)
  
    if update:
      print('updates')
      f1['/stellar_mass_function_moster_2013/counts_AGN'][:] = counts 
      f1['/stellar_mass_function_moster_2013/dN_dVdlogM_AGN'][:] = dN_dVdlogM  

    else:
      print('creates')
     
      ds = f1['stellar_mass_function_moster_2013'].create_dataset('dN_dVdlogM_AGN', data = dN_dVdlogM )
      ds.attrs['units'] = r'$ Mpc^{-3} dex^{-1}$'
      ds.attrs['long_name'] = r'$dN / (dV/, dlogM) $' 

      ds = f1['stellar_mass_function_moster_2013'].create_dataset('counts_AGN', data = counts )
      ds.attrs['units'] = r'count'
      ds.attrs['long_name'] = r'AGN counts' 

  f1.close()
  
#measureSMF(h5_files[50], update=True)
#measureSMF(h5_files[65], update=True)
for h5_file in h5_files[::-1]:
  ##try:
  measureSMF(h5_file, update=False)
  ##except( ValueError ):
  ##pass
