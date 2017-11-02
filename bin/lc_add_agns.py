"""
Routine to add sky coordinates

L1 L2 L3   u11 u12 u13   u21 u22 u23   u31 u32 u33   (periodicity)
'5.9161', '0.4140', '0.4082', '5', '3', '1', '1', '1', '0', '0', '1', '0', '(1)'
'2.4495', '0.7071', '0.5774', '2', '1', '1', '1', '1', '0', '0', '1', '0', '(1)'

python3 lc_add_sky_coordinates.py 1.3

"""


import sys
env = 'MD10'# sys.argv[1] # 'MD10'
L_box = 1000./0.6777 # float(sys.argv[2]) / 0.6777
z_max = 1.3 # float(sys.argv[3])
positions_group_name = 'remaped_position_L3'
x_obs, y_obs, z_obs = 0., 0.7071/2.*L_box, 0.5774/2.*L_box
strech_factor_los = 2.4495
strech_factor_y = 0.7071
strech_factor_z = 0.5774

import h5py    # HDF5 support
import os
import glob
import numpy as n

f_z = lambda z : n.piecewise(z, [z <= 1.1, z > 1.1], [ lambda z : (1.+z)**(5.82), lambda z : (1. + 1.1)**(5.82) * ((1.+z)/(1.+1.1))**(2.36)])

f_Mstar =lambda logM, z : (10**(logM - 10.99) )**(0.24) * n.e**( - 10**(logM - 10.99) )
	
def f_lambda_sar( logM, z, log_lambda_SAR ):
	log_lambda_SAR_var = 10**( log_lambda_SAR - 33.8 + 0.48 * (logM - 11.) )		
	g1z = 1.01 - 0.58 * (z - 1.1)
	return 1. / ( log_lambda_SAR_var**(g1z) + log_lambda_SAR_var**(3.72) )
		


#import XrayLuminosity
#xr = XrayLuminosity.XrayLuminosity()

from scipy.interpolate import interp1d
from scipy.misc import derivative

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)


print("reads the light cone")
file_lc = os.path.join(os.environ[env], 'h5_lc', 'lc_'+positions_group_name+'_obs_0-0_0-7_0-5.hdf5')
f = h5py.File(file_lc,  "r")

# interpolation grid for the Bongiorno et al. 2016 model
dl=0.01
log_lambda_SAR_values = n.arange(32-dl,36+2*dl,dl)

# set the columns
agn_log_lambda_SAR = n.zeros_like(f['/sky_position/redshift_R'].value)
agn_random_number = n.random.random(len(f['/sky_position/redshift_R'].value))
 
# objects in the light cone with reasonnable parameters
ok = (f['/sky_position/selection'].value) & (f['/emerge_data/stellar_mass'].value>0) & (f['/emerge_data/stellar_mass'].value<1e14)
# area subtended by the light cone
area_fraction = (n.max(f['/sky_position/RA'].value[ok])-n.min(f['/sky_position/RA'].value[ok]))*(n.max(f['/sky_position/DEC'].value[ok])-n.min(f['/sky_position/DEC'].value[ok])) * n.pi / 129600.
# derivative of the comoving volume
dvdz = lambda z : derivative(cosmoMD.comoving_volume, z) * area_fraction 

dl=0.01
log_lambda_SAR_values = n.arange(32-dl,36+2*dl,dl)
zs = f['/sky_position/redshift_R'].value[ok]
fz = f_z( zs )
logM = n.log10(f['/emerge_data/stellar_mass'].value[ok])
fMstar = f_Mstar( logM, zs )

#flsar = n.array([f_lambda_sar( logM, zs, log_lambda_SAR_val ) for log_lambda_SAR_val in log_lambda_SAR_values])

psi_log = n.array([10**(- 6.86) * fMstar * fz  * f_lambda_sar( logM, zs, log_lambda_SAR_val ) for log_lambda_SAR_val in log_lambda_SAR_values])

norm = n.sum(psi_log, axis=1)
probas = psi_log	 / norm
	
import time
t0 = time.time()
# computes the probability function per object on the grid of lambda_sar
for ii in range(len(f['/sky_position/RA'].value[ok])):
	probas_un = psi_log(log_lambda_SAR_values, n.log10(f['/emerge_data/stellar_mass'].value[ok][ii]), f['/sky_position/redshift_R'].value[ok][ii]) * dvdz( f['/sky_position/redshift_R'].value[ok][ii]).value * dl 
	norm = n.sum(probas)
	probas = probas_un / norm
	prb_rev = probas[::-1]
	itp_prb = interp1d( n.hstack((n.array([ n.sum(prb_rev[:jj]) for jj in range(len(prb_rev)) ]), 1.)), n.hstack((log_lambda_SAR_values[::-1], log_lambda_SAR_values[0]-dl )) )
	agn_log_lambda_SAR[ok][ii] = itp_prb(n.array([agn_random_number[ok][ii]]))[0]
	print( (time.time()-t0)/(ii+1) )



halo_data = f.create_group('agn_properties')

ds = halo_data.create_dataset('log_lambda_sar', data = agn_log_lambda_SAR )
ds.attrs['units'] = '\log_{10}(\lambda_{SAR}/[])'
ds.attrs['long_name'] = 'log10 of the specific accretion rate' 

ds = halo_data.create_dataset('random_number_lambda_sar', data = agn_random_number )

#random_number_activity = n.random.random(len(f['/sky_position/redshift_R'].value))
#ds = halo_data.create_dataset('random_number_activity', data = random_number_activity )

f.close()
