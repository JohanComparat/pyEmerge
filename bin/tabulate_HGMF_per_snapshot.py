import numpy as n
import glob
import h5py
import os
import time
import sys

from scipy.integrate import quad

f_z = lambda z : n.piecewise(z, [z <= 1.1, z > 1.1], [ lambda z : (1.+z)**(5.82), lambda z : (1. + 1.1)**(5.82) * ((1.+z)/(1.+1.1))**(2.36)])

f_Mstar =lambda logM, z : (10**(logM - 10.99) )**(0.24) * n.e**( - 10**(logM - 10.99) )
	
def f_lambda_sar( logM, z, log_lambda_SAR ):
	lambda_SAR_var = 10**( log_lambda_SAR - 33.8 + 0.48 * (logM - 11.) )		
	g1z = 1.01 - 0.58 * (z - 1.1)
	return 1. / ( lambda_SAR_var**(g1z) + lambda_SAR_var**(3.72) )

psi_log = lambda log_lambda_SAR, logM, z : 10**(- 6.86) * f_lambda_sar( logM, z, log_lambda_SAR ) * f_Mstar(logM, z) * f_z(z)	

psi = lambda lambda_SAR, mass, z : psi_log(n.log10(lambda_SAR), n.log10(mass), z)

stellar_mass = 10**(10.8)
redshift=2.

integrand = lambda lambda_SAR, mass, z : psi( lambda_SAR, mass, z) / lambda_SAR / n.log(10)
fel = n.array([quad(integrand, 10**32, 10**33, args=(stellar_mass, redshift))[0], quad(integrand, 10**33, 10**34, args=(stellar_mass, redshift))[0], quad(integrand, 10**34, 10**35, args=(stellar_mass, redshift))[0], quad(integrand, 10**35, 10**36, args=(stellar_mass, redshift))[0]])
print(n.log10(fel), n.log10(n.sum(fel)))



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

	logMs_low  = f1['/stellar_mass_function_moster_2013/stellar_mass_low'].value
	logMs_up   = f1['/stellar_mass_function_moster_2013/stellar_mass_up'].value 
	counts     = f1['/stellar_mass_function_moster_2013/counts'].value
	dN_dVdlogM = f1['/stellar_mass_function_moster_2013/dN_dVdlogM'].value 


	# interpolates the duty cycle in the interesting region
	maxMS = n.max(logMs_up[(counts>1)])
	minMS = n.min(logMs_low[(counts>1)])
	x_SMF = (logMs_low+ logMs_up)/2.

	sel=(x_SMF>minMS)&(x_SMF<maxMS)

	AGN_HGMF = n.array([xr.Phi_stellar_mass(logMs_i, redshift) for logMs_i in x_SMF])
	duty_cycle = n.zeros_like(dN_dVdlogM)
	duty_cycle[sel] = AGN_HGMF[sel] / dN_dVdlogM[sel]
	duty_cycle[duty_cycle>1] = n.ones_like(duty_cycle[duty_cycle>1])
	print("HGMF", AGN_HGMF[sel])
	print("SMF", dN_dVdlogM[sel])
	print("DC",duty_cycle[sel],n.min(duty_cycle[sel]), n.max(duty_cycle[sel]))
	
	if update:
		print('updates')
		f1['/stellar_mass_function_moster_2013/AGN_HGMF'][:] = AGN_HGMF  
		f1['/stellar_mass_function_moster_2013/duty_cycle'][:] = duty_cycle  

	else:
		print('creates')

		ds = f1['/stellar_mass_function_moster_2013'].create_dataset('AGN_HGMF', data = AGN_HGMF )
		ds.attrs['units'] = r'$ Mpc^{-3} dex^{-1}$'
		ds.attrs['long_name'] = r'$dN / (dV/, dlogM) $' 
		ds = f1['/stellar_mass_function_moster_2013'].create_dataset('duty_cycle', data = duty_cycle )
		ds.attrs['units'] = r'fraction'
		ds.attrs['long_name'] = r'Duty cycle, fraction of active galaxies at a given stellar mass' 

	f1.close()

for h5_file in h5_files[::-1]:
	#try:
	measure_HGMF(h5_file, update=False)
	#except( ValueError ):
	#pass
