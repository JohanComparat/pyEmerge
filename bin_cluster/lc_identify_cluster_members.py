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

lc_id = sys.argv[1]

#status = 'create'
status = 'update'

path_to_glc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_remaped_position_'+lc_id+'.hdf5'
fg = h5py.File(path_to_glc, 'r+')

path_to_clc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_cluster_remaped_position_'+lc_id+'.hdf5'
fc = h5py.File(path_to_clc, 'r')
  
x, y, z = fc['/halo_position/x'].value, fc['/halo_position/y'].value, fc['/halo_position/z'].value
radial_distance_2_member = 10. * fc['/halo_properties/rvir'].value /1000.
position_BCG = n.transpose([x, y, z])

t0=time.time()
treeData = t.cKDTree(n.transpose([ fg['/halo_position/x'].value, fg['/halo_position/y'].value, fg['/halo_position/z'].value ]),10.0)
out = n.array([ treeData.query_ball_point(position_BCG[id_bcg], radial_distance_2_member[id_bcg]) for id_bcg in range(len(position_BCG))])
print('galaxy tree created', time.time()-t0)

clus_id = n.zeros_like(fg['/halo_position/x'].value).astype('int')
clus_flux_05_24 = n.zeros_like(fg['/halo_position/x'].value)

t0=time.time()
for jj in n.arange(len(position_BCG)):
  clus_id[out[jj]] = (fc['/halo_properties/id'].value).astype('int')[jj]*n.ones_like(out[jj]).astype('int')
  clus_flux_05_24[out[jj]] = fc['/cluster_data/rxay_flux_05_24'].value[jj]*n.ones_like(out[jj])

print('clusters assigned', time.time()-t0)

if status == 'create' :
  halo_data = fg.create_group('cluster_galaxies')
  ds = halo_data.create_dataset('cluster_id', data = clus_id )
  ds = halo_data.create_dataset('cluster_flux_05_24', data = clus_flux_05_24 )
  
if status=='update' :
  fg['/cluster_galaxies/cluster_id'][:]  = clus_id 
  fg['/cluster_galaxies/cluster_flux_05_24'][:]  = clus_flux_05_24 

fc.close()
fg.close()