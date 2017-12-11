"""

+ + + + + + + HEADER + + + + + + + + +
file_name lc_remaped_position_L3_.hdf5
HDF5_Version 1.8.18
h5py_version 2.7.1


+ + + + + + + DATA   + + + + + + + + + +
========================================
agn_properties <HDF5 group "/agn_properties" (2 members)>
- - - - - - - - - - - - - - - - - - - - 
agn_activity <HDF5 dataset "agn_activity": shape (23191107,), type "<f8">
log_lambda_sar <HDF5 dataset "log_lambda_sar": shape (23191107,), type "<f8">
========================================
emerge_data <HDF5 group "/emerge_data" (3 members)>
- - - - - - - - - - - - - - - - - - - - 
dMdt <HDF5 dataset "dMdt": shape (23191107,), type "<f8">
mvir_dot <HDF5 dataset "mvir_dot": shape (23191107,), type "<f8">
rvir_dot <HDF5 dataset "rvir_dot": shape (23191107,), type "<f8">
========================================
halo_position <HDF5 group "/halo_position" (7 members)>
- - - - - - - - - - - - - - - - - - - - 
vx <HDF5 dataset "vx": shape (23191107,), type "<f8">
vy <HDF5 dataset "vy": shape (23191107,), type "<f8">
vz <HDF5 dataset "vz": shape (23191107,), type "<f8">
x <HDF5 dataset "x": shape (23191107,), type "<f8">
y <HDF5 dataset "y": shape (23191107,), type "<f8">
z <HDF5 dataset "z": shape (23191107,), type "<f8">
z_snap <HDF5 dataset "z_snap": shape (23191107,), type "<f8">
========================================
halo_properties <HDF5 group "/halo_properties" (7 members)>
- - - - - - - - - - - - - - - - - - - - 
Mpeak <HDF5 dataset "Mpeak": shape (23191107,), type "<f8">
Vmax <HDF5 dataset "Vmax": shape (23191107,), type "<f8">
id <HDF5 dataset "id": shape (23191107,), type "<f8">
mvir <HDF5 dataset "mvir": shape (23191107,), type "<f8">
pid <HDF5 dataset "pid": shape (23191107,), type "<f8">
rs <HDF5 dataset "rs": shape (23191107,), type "<f8">
rvir <HDF5 dataset "rvir": shape (23191107,), type "<f8">
========================================
moster_2013_data <HDF5 group "/moster_2013_data" (1 members)>
- - - - - - - - - - - - - - - - - - - - 
stellar_mass <HDF5 dataset "stellar_mass": shape (23191107,), type "<f8">
========================================
sky_position <HDF5 group "/sky_position" (5 members)>
- - - - - - - - - - - - - - - - - - - - 
DEC <HDF5 dataset "DEC": shape (23191107,), type "<f8">
RA <HDF5 dataset "RA": shape (23191107,), type "<f8">
redshift_R <HDF5 dataset "redshift_R": shape (23191107,), type "<f8">
redshift_S <HDF5 dataset "redshift_S": shape (23191107,), type "<f8">
selection <HDF5 dataset "selection": shape (23191107,), type "|b1">

TO DELETE PREVIOUS DATA : 

path_to_lc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_remaped_position_L3.hdf5'

f = h5py.File(path_to_lc, 'a')
del f['/agn_properties/rxay_flux_05_20']
del f['/agn_properties/logNH']
f.close()

"""
"""
Convert to observed fluxes

intrinsic extinction. Thin / thick obscuration

Follows Buchner et al. 2016

"""


import h5py    # HDF5 support
import os
import glob
import numpy as n
from scipy.interpolate import interp1d
import sys
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)

from scipy.special import erf
ricci_ct_f = lambda z: 0.22 + 0.18 * z**0.4
#fraction_ricci = lambda lsar, z : ricci_ct_f(z)+(0.8-ricci_ct_f(z))*(0.5+0.5*erf((-lsar+34.7)/0.4))
#print(fraction_ricci(n.array([32, 32.7, 33, 34, 35, 36, 37]), 0.))
fraction_ricci = lambda lsar, z : ricci_ct_f(z)+(0.8-ricci_ct_f(z))*(0.5+0.5*erf((-lsar+35.2)/0.4))
print(fraction_ricci(n.array([32, 32.7, 33, 34, 35, 36, 37]), 0.))

#status = 'create'
status = 'update'

model_NH = 'ricci_2017'
# model_NH = 'buchner_2017'

path_to_lc = sys.argv[1]
#path_to_lc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_remaped_position_L3.hdf5'

f = h5py.File(path_to_lc, 'r+')

is_gal = (f['/sky_position/selection'].value)&(f['/sky_position/redshift_R'].value<3.)
is_agn = (is_gal)&(f['/agn_properties/agn_activity'].value==1)

n_gal = len(f['/sky_position/redshift_S'].value[is_gal])

n_agn = len(f['/sky_position/redshift_S'].value[is_agn])

z = f['/sky_position/redshift_S'].value[is_agn]
logm = n.log10(f['/moster_2013_data/stellar_mass'].value[is_agn])
lsar = f['/agn_properties/log_lambda_sar'].value[is_agn]
print('lsar',lsar, n.min(lsar), n.max(lsar))
lx = logm + lsar 

logNH = n.random.uniform(20, 22, n_agn)
obs_type = n.zeros(n_agn)

if model_NH == 'ricci_2017':
	# obscuration, Ricci + 2017
	randomNH = n.random.rand(n_agn)
	# thin obscuration 20 - 22
	frac_thin = fraction_ricci(lsar, z)
	print('frac thin min', n.min(frac_thin))
	thinest = (randomNH >= frac_thin)
	# 22% of thick, 24-26
	thick = (randomNH < 0.22)
	# obscured 22-24
	obscured = (thinest==False)&(thick==False)
	# assigns logNH values randomly :
	print(n_agn, len(thick.nonzero()[0]), len(obscured.nonzero()[0]), len(thinest.nonzero()[0]))
	logNH[thick] = n.random.uniform(24, 26, len(logNH[thick]))
	obs_type[thick] = n.ones_like(logNH[thick])*2
	logNH[obscured] = n.random.uniform(22, 24, len(logNH[obscured]))
	obs_type[obscured] =  n.ones_like(logNH[obscured])
	

if model_NH == 'buchner_2017':
	# randomObscuration = n.random.rand(n_gal)		
	# obscured = (xr.obscured_fraction_optical_Merloni2015(lsar + stellar_mass) < randomObscuration )

	# obscuration Buchner et al. 2015 + 2017
	# add the log NH of the logNH_host
	# 35 % have a thick obscuration 24 - 26
	# 65 % have a thin obscuration that depends on stellar mass and Xray luminosity
	logNH = n.random.uniform(20, 22, n_agn)
	obs_type = n.zeros(n_agn)

	# 35% of thick, 24-26
	randomNH = n.random.rand(n_agn)
	thick_obscuration = (randomNH < 0.35)
	thin_obscuration = (randomNH >= 0.35)
	#med_obscuration = (randomNH >= 0.6)

	logNH[thick_obscuration] = n.random.uniform(24, 26, len(logNH[thick_obscuration]))
	obs_type[thick_obscuration] = n.ones_like(logNH[thick_obscuration])*2

	# the thin : about 40 % are thin whatever happens: 22-24
	logNH_host_mean = 21.7 + (logm- 9.5)*0.38
	logNH_host = n.random.normal(logNH_host_mean, 0.5)
	L_transition  = 21.62
	logNH[(thin_obscuration)&(logNH_host>= L_transition)] = n.random.uniform(22, 24, len(logNH[(thin_obscuration)&(logNH_host>=L_transition)]))
	obs_type[(thin_obscuration)&(logNH_host>= L_transition)] =  n.ones_like(logNH[(thin_obscuration)&(logNH_host>=L_transition)])
	# a few more are thin depending on their Xray luminosity: 22-24
	#rest = (thin_obscuration)&(logNH_host<22)
	#randomNH2 = n.random.rand(n_agn)
	#rest_obscured = (rest)&(randomNH2 < obscured_fraction_interpolated(lsar + logm))
	#logNH[(rest_obscured)] = random.uniform(22, 24, len(logNH[(rest_obscured)]))
	#obs_type[(rest_obscured)] =  n.ones_like(logNH[(rest_obscured)])

	#logNH_host[logNH_host<=20] = n.random.uniform(20, 22, len(logNH_host[logNH_host<=20]))
	#logNH_host[logNH_host>=26] = n.random.uniform(24, 26, len(logNH_host[logNH_host>=26]))
	#logNH = logNH_host

obscuration_z_grid, obscuration_nh_grid, obscuration_fraction_obs_erosita = n.loadtxt( os.path.join( os.environ['GIT_NBODY_NPT'], "data", "AGN", "fraction_observed_by_erosita_due_2_obscuration.txt"), unpack=True)
nh_vals = 10**n.arange(-2,4,0.05)
z_vals = 10**n.arange(-3,0.68,0.025)

obscuration_interpolation_grid = n.array([ 
  interp1d(
    n.hstack((obscuration_nh_grid[ (obscuration_z_grid==zz) ], 26.)), 
    n.hstack((obscuration_fraction_obs_erosita[( obscuration_z_grid==zz) ], obscuration_fraction_obs_erosita[( obscuration_z_grid==zz) ][-1]))
		) 
  for zz in z_vals])

index = n.searchsorted(z_vals, z)

percent_observed = n.array([ obscuration_interpolation_grid[ind]( nh) for ind, nh in zip(index, logNH)])

# lx is 2-10 keV
# converts in the erosita band 0.5-2 keV
lx_absorbed_05_20 = n.log10(10**lx * percent_observed)

d_L = cosmoMD.luminosity_distance(z)
dl_cm = (d_L.to(u.cm)).value
adjusting_factor = 0.
fx_05_20 = 10**(lx_absorbed_05_20-adjusting_factor) / (4 * n.pi * dl_cm**2.)

fx_05_20_out = n.ones_like(f['/sky_position/redshift_S'].value)*-9999.
logNH_out = n.ones_like(f['/sky_position/redshift_S'].value)*-9999.
fx_05_20_out[is_agn] = fx_05_20
logNH_out[is_agn] = logNH

if status == 'create':
  f['/agn_properties'].create_dataset('rxay_flux_05_20', data = fx_05_20_out )
  f['/agn_properties'].create_dataset('logNH', data = logNH_out )

if status == 'update':
  f['/agn_properties/rxay_flux_05_20'][:] = fx_05_20_out 
  f['/agn_properties/logNH'][:] = logNH_out 

#ds = f['/agn_properties'].create_dataset('rxay_flux_05_20', data = fx_05_20 )
#ds.attrs['units'] = 'erg/cm2/s'
#ds.attrs['long_name'] = 'X ray flux in the 0.5-2 keV band' 

#ds = f['/agn_properties'].create_dataset('logNH', data = logNH )
#ds.attrs['units'] = 'logNH'
#ds.attrs['long_name'] = 'logNH' 

f.close()
