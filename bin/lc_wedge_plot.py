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


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p

plotDir = os.path.join(os.environ['HOME'], 'wwwDir', "eRoMok", "h5")

def get_positions(path_to_lc, area, z_max=3., ra_max=1.):
  f = h5py.File(path_to_lc, 'r')
  is_gal = (abs(f['/sky_position/RA'].value)<ra_max) & (f['/sky_position/selection'].value) & (f['/sky_position/redshift_R'].value<z_max)
  is_agn = (is_gal) & (f['/agn_properties/agn_activity'].value==1) & (f['/agn_properties/rxay_flux_05_20'].value>0)
  #n_gal = len(f['/sky_position/redshift_S'].value[is_gal])
  #n_agn = len(f['/sky_position/redshift_S'].value[is_agn])
  #zR_g = f['/sky_position/redshift_R'].value[is_gal]
  #zS_g = f['/sky_position/redshift_S'].value[is_gal]
  #dec_g = f['/sky_position/DEC'].value[is_gal]*n.pi/180.
  zR_a = f['/sky_position/redshift_R'].value[is_agn]
  zS_a = f['/sky_position/redshift_S'].value[is_agn]
  dec_a = f['/sky_position/DEC'].value[is_agn]*n.pi/180.
  f.close()
  # return zR_g, zS_g, dec_g, zR_a, zS_a, dec_a 
  return zR_a, zS_a, dec_a 

p.figure(2, (10,4))
p.axes([0,0,1,1])
path_to_lc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_remaped_position_L15.hdf5'
area = 14.323944878104827*2. * 2*20.257311381848154
#zR_g, zS_g, dec_g, zR_a, zS_a, dec_a = get_positions(path_to_lc, area, z_max=3.)
zR_a, zS_a, dec_a = get_positions(path_to_lc, area, z_max=3.)
p.plot(zR_a, zR_a*n.tan(dec_a), 'r,', alpha=0.2, rasterized = True, label = 'L15 z<0.54 1160deg2'  )

path_to_lc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_remaped_position_L3.hdf5'
area = 6.7529257176359*2. * 2* 8.269819492449505
#zR_g, zS_g, dec_g, zR_a, zS_a, dec_a = get_positions(path_to_lc, area, z_max=1.1)
zR_a, zS_a, dec_a = get_positions(path_to_lc, area, z_max=1.1)
p.plot(zR_a, zR_a*n.tan(dec_a), 'b,', alpha=0.2, rasterized=True, label = 'L3 z<1.08, 223deg2')

path_to_lc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_remaped_position_L6.hdf5'
area = 1.9766516114702513*2. * 2*2.0047373031569915
#zR_g, zS_g, dec_g, zR_a, zS_a, dec_a  = get_positions(path_to_lc, area, z_max=3.)
zR_a, zS_a, dec_a  = get_positions(path_to_lc, area, z_max=3.)
p.plot(zR_a, zR_a*n.tan(dec_a), 'k,', alpha=0.2, rasterized=True, label = 'L6 z<3, 15deg2')

p.xlabel('redshift')
p.ylabel('DEC')
p.legend(frameon=False, loc=0)
p.xscale('log')
p.xlim((0,3))
p.ylim((-0.2, 0.2))
#p.title('Mocks')
#p.grid()
p.axis('off')
p.savefig(os.path.join(plotDir, "wedges_AGN.jpg"))
p.clf()


