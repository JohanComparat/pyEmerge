import numpy as n
import glob
import h5py
import os
import time
import sys

import XrayLuminosity
xr = XrayLuminosity.XrayLuminosity()

h5_files = n.array(glob.glob(os.path.join(os.environ['MD10'], "h5", "hlist_?.?????_emerge.hdf5")))
h5_files.sort()

#bins = n.arange(6,13,0.1)
#xb = (bins[1:] + bins[:-1]) / 2.
#hh=0.6777


def measure_HGMF(h5_file, update=True):
	f1 = h5py.File(h5_file,  "r+")
	redshift = f1.attrs['redshift']

	logMs_low  = f1['/stellar_mass_function/stellar_mass_low'].value
	logMs_up   = f1['/stellar_mass_function/stellar_mass_up'].value 
	counts     = f1['/stellar_mass_function/counts'].value
	dN_dVdlogM = f1['/stellar_mass_function/dN_dVdlogM'].value 


	# interpolates the duty cycle in the interesting region
	maxMS = n.max(logMs_up[(counts>1)])
	minMS = n.min(logMs_low[(counts>1)])
	x_SMF = (logMs_low+ logMs_up)/2.

	sel=(x_SMF>minMS)&(x_SMF<maxMS)

	AGN_HGMF = n.array([xr.Phi_stellar_mass(logMs_i, redshift) for logMs_i in x_SMF])
	duty_cycle = n.zeros_like(dN_dVdlogM)
	duty_cycle[sel] = AGN_HGMF[sel] / dN_dVdlogM[sel]
	print("HGMF", AGN_HGMF[sel])
	print("SMF", dN_dVdlogM[sel])
	print("DC",duty_cycle[sel],n.min(duty_cycle[sel]), n.max(duty_cycle[sel]))
	
	if update:
		print('updates')
		f1['/stellar_mass_function/AGN_HGMF'][:] = AGN_HGMF  
		f1['/stellar_mass_function/duty_cycle'][:] = duty_cycle  

	else:
		print('creates')

		ds = f1['/stellar_mass_function'].create_dataset('AGN_HGMF', data = AGN_HGMF )
		ds.attrs['units'] = r'$ Mpc^{-3} dex^{-1}$'
		ds.attrs['long_name'] = r'$dN / (dV/, dlogM) $' 
		ds = f1['/stellar_mass_function'].create_dataset('duty_cycle', data = duty_cycle )
		ds.attrs['units'] = r'fraction'
		ds.attrs['long_name'] = r'Duty cycle, fraction of active galaxies at a given stellar mass' 

	f1.close()

for h5_file in h5_files[::-1][:3]:
	#try:
	measure_HGMF(h5_file, update=True)
	#except( ValueError ):
	#pass
