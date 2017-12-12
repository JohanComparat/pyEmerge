"""
identifies all cluster members within 10 r500c in the galaxy light cone
"""
import glob
#import astropy.io.fits as fits
import os
import time
import numpy as n
import sys
import scipy.spatial.ckdtree as t
import h5py    # HDF5 support

from scipy.stats import norm
fun = lambda mmm, sss : norm.rvs( loc = mmm, scale = sss )

from scipy.interpolate import interp1d
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115)#, Ob0=0.048206)

lc_id = sys.argv[1]

status = 'create'
#status = 'update'

path_to_glc =  '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_'+lc_id+'.hdf5'
f = h5py.File(path_to_glc, 'r+')

is_gal = (f['/sky_position/selection'].value)&(f['/cluster_galaxies/cluster_id'].value > 0 ) & (f['/cluster_galaxies/cluster_flux_05_24'].value> 0) &(f['/cluster_galaxies/d_cluster_center'].value<2.5)

red_sequence_flag = n.ones_like(f['/cluster_galaxies/d_cluster_center']).astype('int')

r_bins = n.arange(0,2.6,0.5)

t0 = time.time()

cl_ids = n.array(list(set(f['/cluster_galaxies/cluster_id'].value[is_gal])))
print(len(cl_ids),'clusters')
for cl_id in  cl_ids:
	cluster = (f['/cluster_galaxies/cluster_id'].value==cl_id)#&(is_gal)
	NALL =  len(cluster.nonzero()[0])
	print(NALL, time.time())
	if NALL>3:
		z_mean = n.mean(f['/sky_position/redshift_R'].value[cluster])
		N_total = len(cluster.nonzero()[0])
		x = 0.5*(r_bins[:-1] + r_bins[1:])
		frac_old = ( x**(-0.25) - x/100. - 0.47 )*(1.+ z_mean)**2./3.2
		#print(cl_id, z_mean, N_total, frac_old)
		for ii, (rmin, rmax) in enumerate(zip(r_bins[:-1], r_bins[1:])):
			sel=(f['/cluster_galaxies/d_cluster_center'].value>rmin) &(f['/cluster_galaxies/d_cluster_center'].value<rmax) &(cluster)
			NGAL = len(sel.nonzero()[0])
			#print(NGAL)
			if NGAL>0:
				rds = n.random.rand(NGAL)
				red_sequence_flag[sel][rds < frac_old[ii]] = n.ones_like(red_sequence_flag[sel][rds < frac_old[ii]]).astype('int')
 
#sys.exit()

if status == 'create' :
  f['/cluster_galaxies'].create_dataset('red_sequence_flag', data = red_sequence_flag )
  
if status=='update' :
  f['/cluster_galaxies/red_sequence_flag'][:]  = red_sequence_flag 

f.close()