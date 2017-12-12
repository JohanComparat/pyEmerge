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

#status = 'create'
status = 'update'

path_to_glc =  '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_'+lc_id+'.hdf5'
fg = h5py.File(path_to_glc, 'r+')

path_to_clc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_cluster_Dietrich_2017/lc_cluster_remaped_position_'+lc_id+'.hdf5'
fc = h5py.File(path_to_clc, 'r')
  
x, y, z = fc['/halo_position/x'].value, fc['/halo_position/y'].value, fc['/halo_position/z'].value
radial_distance_2_member = 2.5 * fc['/halo_properties/rvir'].value /1000.
position_BCG = n.transpose([x, y, z])

t0=time.time()
treeData = t.cKDTree(n.transpose([ fg['/halo_position/x'].value, fg['/halo_position/y'].value, fg['/halo_position/z'].value ]),10.0)
out = n.array([ treeData.query_ball_point(position_BCG[id_bcg], radial_distance_2_member[id_bcg]) for id_bcg in range(len(position_BCG))])
print('galaxy tree created', time.time()-t0)

radial_separation = n.zeros_like(fg['/halo_position/x'].value)
clus_id = n.zeros_like(fg['/halo_position/x'].value).astype('int')
clus_flux_05_24 = n.zeros_like(fg['/halo_position/x'].value)
bcg_flag = n.zeros_like(fg['/halo_position/x'].value).astype('int')
#bcg_ids = []

t0=time.time()
for jj in n.arange(len(position_BCG)):
  clus_id[out[jj]] = (fc['/halo_properties/id'].value).astype('int')[jj]*n.ones_like(out[jj]).astype('int')
  clus_flux_05_24[out[jj]] = fc['/cluster_data/rxay_flux_05_24'].value[jj]*n.ones_like(out[jj])
  radial_distance = n.sum((position_BCG[jj]-treeData.data[out[jj]])**2., axis=1)**0.5 / radial_distance_2_member[jj] * 2.5
  radial_separation[out[jj]] =  radial_distance
  #bcg_ids.append( out[jj][n.argmin(radial_distance)] )
  #print(n.min(radial_separation[out[jj]]), bcg_ids[-1])

print('clusters assigned', time.time()-t0)

#bcg_ids = n.array(bcg_ids)
#bcg_flag[bcg_ids] = n.ones_like(bcg_flag[bcg_ids])

# relation fitted on magneticum light cones. Works for Vmax > 100 (80 is ok, lower it is wrong)
Mag_abs_i = lambda x : - 24.6 - 5.9*( x*(1+x*x) - 0.3)
Mag_abs_r = lambda x : - 24.1 - 5.9*( x*(1+x*x) - 0.3)
sigma_Mi = lambda x : x**(-3.5) * 24. * 2.**0.5 / 5.

#in_clus=(clus_id>0)

#log_vmax = n.log10(fg['/halo_properties/Vmax'].value[in_clus])
#mean_Mr = Mag_abs_r(log_vmax-2.5)
#scatter = sigma_Mi(log_vmax)

#dist_mods = cosmoMD.distmod(fg['/sky_position/redshift_R'].value[in_clus])
#mag_abs = n.array([fun(mm,ss) for mm, ss in zip(mean_Mr, scatter)])
#mag_r = n.zeros_like(fg['/halo_position/x'].value)
#mag_r[in_clus] = mag_abs + dist_mods.value

#fg['/cluster_galaxies'].create_dataset('mag_r', data = mag_r )
#fg['/cluster_galaxies'].create_dataset('bcg_flag', data = bcg_flag )
fg['/cluster_galaxies'].create_dataset('d_cluster_center', data = radial_separation)

if status == 'create' :
  halo_data = fg.create_group('cluster_galaxies')
  ds = halo_data.create_dataset('cluster_id', data = clus_id )
  ds = halo_data.create_dataset('cluster_flux_05_24', data = clus_flux_05_24 )
  #ds = halo_data.create_dataset('mag_r', data = mag_r )
  ds = halo_data.create_dataset('bcg_flag', data = bcg_flag )
    
if status=='update' :
  fg['/cluster_galaxies/cluster_id'][:]  = clus_id 
  fg['/cluster_galaxies/cluster_flux_05_24'][:]  = clus_flux_05_24 
  #fg['/cluster_galaxies/mag_r'][:]  = mag_r 
  fg['/cluster_galaxies/bcg_flag'][:] = bcg_flag 
  fg['/cluster_galaxies/d_cluster_center'][:] = radial_separation

fc.close()
fg.close()