from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115)#, Ob0=0.048206)

import glob
import os
import time
import numpy as n
import sys

# specific functions
from scipy.stats import norm
from scipy.integrate import quad
from scipy.interpolate import interp1d
import random
from scipy.stats import norm

import h5py    # HDF5 support

status = 'create'
#status = 'update'

lc_id = sys.argv[1]
path_to_lc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_'+lc_id+'.hdf5'


if lc_id=='L3':
  z_max = 1.3
  area = 2*2*6.7529257176359 * 8.269819492449505

if lc_id=='L3_z1':
  z_max = 3.
  area = 2*2*3.3764628588325674*4.134909746242654

if lc_id=='L6':
  z_max = 10.
  area = 2*2*1.9766516114702513*2.0047373031569915
  
if lc_id=='L15':
  z_max = 3.
  area = 2*2*14.323944878104827*20.257311381848154
  
f = h5py.File(path_to_lc, 'r+')

is_gal = (f['/sky_position/selection'].value)&(f['/sky_position/redshift_R'].value<3.)
N_halos = len(f['/sky_position/redshift_S'].value)
n_gal = len(f['/sky_position/redshift_S'].value[is_gal])

z = f['/sky_position/redshift_S'].value

EBE_zmin, EBE_zmax, EBE_SGC_long, EBE_SGC_nom, EBE_NGC = n.loadtxt(os.path.join(os.environ['GIT_EMERGE'], "data/NZ/eBOSS-ELG-Raichoor-17-deg2.txt"), unpack=True)

elg_selection =  (n.ones(N_halos)==0)

# ELG parameters
# ELG select on Mvir
all_mvir = n.log10(f1['halo_properties/mvir'].value)
mh_mean, mh_scatter = 12.2, 0.25

for zmin,zmax,N_p_deg2 in zip(EBE_zmin, EBE_zmax, EBE_SGC_long):
	z_sel = (is_gal)&(z>=zmin)&(z<zmax)
	N_elg  = int(area * N_p_deg2 ) # 2
	if N_elg > 100 :
		mh_bins = n.arange(mh_mean -2*mh_scatter, mh_mean +2*mh_scatter+0.05, 0.05)
		mh_bins_pos = 0.5*(mh_bins[1:]+mh_bins[:-1])
		proba = lambda x : norm.pdf(x, loc=mh_mean,scale=mh_scatter)
		proba_norm = proba(mh_bins_pos).sum()
		N_2_select_per_bin = (N_elg*proba(mh_bins_pos)/proba_norm).astype('int')
		for id_bin in range(len(mh_bins)-1):
			id_in_bin =(z_sel) & (all_mvir > mh_bins[id_bin]) &( all_mvir < mh_bins[id_bin+1]) 
			N_avail = len(id_in_bin.nonzero()[0])
			rds = n.random.rand(len(all_MS))
			bin_selection = (id_in_bin)&(rds < N_2_select_per_bin[id_bin]*1./N_avail)
			elg_selection = (bin_selection)|(elg_selection)



if status == 'create' :
  halo_data = f1.create_group('cosmo_4most')
  ds = halo_data.create_dataset('is_ELG_eBOSS', data = elg_selection )
  
if status=='update' :
  f1['/cosmo_4most/is_ELG_eBOSS'][:]  = elg_selection 
  
f1.close()
