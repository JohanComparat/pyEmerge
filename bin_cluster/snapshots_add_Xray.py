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

from scipy.stats import lognorm
from scipy.stats import norm

import astropy.cosmoMDlogy as co
import astropy.units as u
import astropy.constants as cc
from astropy.cosmoMDlogy import FlatLambdaCDM
cosmoMDMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115)

from scipy.interpolate import interp1d

"""

definitions
-----------
	- Planck flat LCDM cosmoMDlogy
	- :math:`m = ln(M_{500} / (10^{15} M_\odot))`
	- :math:`m_{gas} = ln(M_{gas, 500} / (10^{15} M_\odot))` is the gas mass within r500
	- :math:`m_{lens} = ln(M_{lens, 500} / (10^{15} M_\odot))` is the spherical mass estimate from lensing corresponding to an
idealized shear profile without statistical noise
	- :math:`l = ln(L_{500} / (E(z) 10^{44} erg/s))` where L is defined as the cluster rest-frame luminosity in the 0.1 - 2.4 keV band.
	- :math:`t = ln(kT_{500} / keV)` where kT is the emission weighted temperature measured in annulus from 0.15 to 1.0 r500
	- :math:`E(z) = H(z)/H_0`
	- :math:`\epsilon = ln(E(z))`

Parameters
----------

	* normalization for parameter X, :math:`N_X`
	* slope for E(z) for parameter X, :math:`slope_E_X`
	* slope for M500 for parameter X, :math:`slope_{M500}_X`
  

Workflow
--------
	- Select relaxed clusters from the DM point of view. to be defined how ... with T/U ? Look at the publications from Sembolini, Yepes, Knebe ...
	- using M15, add gas density profile and temperature using scaling relations

"""

N_Mgas = 31.98
N_kT 	= 2.18
N_L 	= 103.7
N_Lce 	= 102.66

slope_E_Mgas 	= -0.11
slope_E_kT 	= 0.61
slope_E_L 		= 1.20
slope_E_Lce 	= 1.82

slope_M500_Mgas= 1.04
slope_M500_kT 	= 0.66
slope_M500_L 	= 1.26
slope_M500_Lce = 1.06

scatter_Mgas = 0.086
scatter_kT = 0.18
scatter_L = 0.24
scatter_Lce = 0.17

E035 = cosmoMD.efunc(0.35)

# converts logM500 to clusters observables
m500_to_qty = lambda logM500, z, slope_efunc, slope_m500, normalization : n.e**normalization * (cosmoMD.efunc(z)/E035)**(slope_efunc) * (10**(logM500-n.log10(6)-14))**(slope_m500)

logM500_to_logMgas 	= lambda logM500, z : m500_to_qty( logM500, z, slope_E_Mgas, slope_M500_Mgas, N_Mgas)
logM500_to_kT 		= lambda logM500, z : m500_to_qty( logM500, z, slope_E_kT, slope_M500_kT, N_kT)
logM500_to_L 		= lambda logM500, z : m500_to_qty( logM500, z, slope_E_L, slope_M500_L, N_L)
logM500_to_Lce		= lambda logM500, z : m500_to_qty( logM500, z, slope_E_Lce, slope_M500_Lce, N_Lce)


file_1 = input_list[ii]
print(file_1)
f1 = h5py.File(file_1,  "r+")

z = f1.attrs['redshift']
log_m500c = n.log10(f1['/halo_properties/M500c'].value)

rds = (n.random.rand(len(log_m500c))-0.5)*2.                   
norm.rvs(loc=0,scale=0.086,size=5)

Mean_Mgas	= n.log10(logM500_to_logMgas	(log_m500c, z))
scatter_Mgas = norm.cdf(rds*3*scatter_Mgas, scale=scatter_Mgas)
VAL_Mgas    = Mean_Mgas + scatter_Mgas

Mean_kT 		= logM500_to_kT		(log_m500c, z)
scatter_kT = norm.cdf(rds, scale=scatter_kT)
VAL_kT    = Mean_kT + scatter_kT

Mean_L 		= logM500_to_L 		(log_m500c, z)
scatter_L = norm.cdf(rds, scale=scatter_L)
VAL_L    = Mean_L + scatter_L

Mean_Lce		= logM500_to_Lce	(log_m500c, z)
scatter_Lce = norm.cdf(rds, scale=scatter_Lce)
VAL_Lce    = Mean_Lce + scatter_Lce



if status=='create':
  halo_data = f1.create_group('cluster_data')

  ds = halo_data.create_dataset('cool_class', data = cool_class.astype('int') )
  ds.attrs['units'] = 'coolness class'
  ds.attrs['long_name'] = '1: very relaxed cool core, 2: relaxed cool core, 3: disturbed non cool core, 4: very disturbed non cool core' 

if status=='update':
  f1['/cluster_data/cool_class'][:] = cool_class.astype('int')

f1.close()


