import h5py    # HDF5 support
import os
import glob
import numpy as n
from scipy.interpolate import interp1d
import sys
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)

#status = 'create'
status = 'update'

path_to_lc = sys.argv[1]
#path_to_lc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_cluster_remaped_position_L3.hdf5'

f = h5py.File(path_to_lc, 'r+')

is_gal = (f['/sky_position/selection'].value)&(f['/sky_position/redshift_R'].value<3.)

z = f['/sky_position/redshift_S'].value[is_gal]
lx = f['/cluster_data/log_LceX_05_24'].value[is_gal] 

percent_observed = 1.
lx_absorbed_05_24 = n.log10(10**lx * percent_observed)

d_L = cosmoMD.luminosity_distance(z)
dl_cm = (d_L.to(u.cm)).value
adjusting_factor = 0.35 # accounts for absorption for now !
fx_05_24 = 10**(lx_absorbed_05_24-adjusting_factor) / (4 * n.pi * dl_cm**2.)

fx_05_24_out = n.ones_like(f['/sky_position/redshift_S'].value)*-9999.
fx_05_24_out[is_gal] = fx_05_24

if status == 'create':
  f['/cluster_data'].create_dataset('rxay_flux_05_24', data = fx_05_24_out )
  
if status == 'update':
  f['/cluster_data/rxay_flux_05_24'][:] = fx_05_24_out 

f.close()
