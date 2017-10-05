import numpy as n
import glob
import h5py
import os

from scipy.interpolate import interp1d

import XrayLuminosity
xr = XrayLuminosity.XrayLuminosity()

h5_files = n.array(glob.glob(os.path.join(os.environ['MD10'], "h5", "hlist_?.?????_emerge.hdf5")))
h5_files.sort()

# compare the stellar mass function measured to the Ilbert function
# take the AGN HGMF model

def tabulate_duty_cycle(h5_file, update=True):
  print(h5_file)
  f1 = h5py.File(h5_file,  "r+")
  redshift = f1.attrs['redshift']
  # Loads mass function
  logMs_low    = f1['stellar_mass_function/stellar_mass_low'].value
  logMs_up     = f1['stellar_mass_function/stellar_mass_up'].value
  counts       = f1['stellar_mass_function/counts'].value
  dN_dVdlogM = f1['stellar_mass_function/dN_dVdlogM'].value 
  # interpolates the model of the AGN host galaxy mass function
  AGN_HGMF = interp1d(n.arange(5.5,13.5,0.01), n.array([xr.Phi_stellar_mass(logMs_i, redshift) for logMs_i in n.arange(5.5,13.5,0.01)]))
  # interpolates the duty cycle in the interesting region
  maxMS = n.max(logMs_up[(counts>1)])
  minMS = n.min(logMs_low[(counts>1)])
  x_SMF = (logMs_low+ logMs_up)/2.
  sel=(x_SMF>minMS)&(x_SMF<maxMS)
  duty_cycle = AGN_HGMF(x_SMF[sel]) / dN_dVdlogM[sel]
  # writes the results
  if update:
    print('updates')
    #print(f1['/stellar_mass_function/stellar_mass_low'].value)
    f1['/agn_model/stellar_mass'][:] = AGN_HGMF.x 
    f1['/agn_model/HGMF'][:] = AGN_HGMF.y
    f1['/agn_model/stellar_mass_duty_cycle'][:] = x_SMF[sel]
    f1['/agn_model/duty_cycle'][:] = duty_cycle
  else:
    print('creates')
    model_agn_data = f1.create_group('agn_model')
    ds = model_agn_data.create_dataset('stellar_mass', data = AGN_HGMF.x )
    ds.attrs['units'] = r'$h^{-1} M_\odot$'
    ds.attrs['long_name'] = r'$h^{-1} M_\odot$' 
    ds = model_agn_data.create_dataset('HGMF', data = AGN_HGMF.y  )
    ds.attrs['units'] = r'$ Mpc^{-3} dex^{-1}$'
    ds.attrs['long_name'] = r'AGN host galaxy stellar mass function' 
    ds = model_agn_data.create_dataset('stellar_mass_duty_cycle', data = x_SMF[sel] )
    ds.attrs['units'] = r'$h^{-1} M_\odot$'
    ds.attrs['long_name'] = r'$h^{-1} M_\odot$' 
    ds = model_agn_data.create_dataset('duty_cycle', data = duty_cycle  )
    ds.attrs['units'] = r'fraction'
    ds.attrs['long_name'] = r'duty cycle' 

  f1.close()

for el in h5_files[20:]:
  try:
    tabulate_duty_cycle(el, update=False)
  except( ValueError, KeyError ):
    pass
  

