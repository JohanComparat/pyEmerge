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
import astropy.io.fits as fits

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)


def write_fits_lcC(path_to_lc, out_filename, z_min, z_max, dec_max, ra_max):
  f = h5py.File(path_to_lc, 'r')

  is_gal = (f['/sky_position/selection'].value)&(f['/cluster_galaxies/cluster_id'].value > 0 ) & (f['/cluster_galaxies/cluster_flux_05_24'].value> 0) & (f['/sky_position/redshift_R'].value<z_max)

  hdu_cols  = fits.ColDefs([
  #fits.Column(name='is_agn',format='K',array= (f['/agn_properties/agn_activity'].value[is_gal]).astype('int'), unit='=1 is AGN' ) 
  fits.Column(name='log_xray_flux_05_24',format='D',    array= n.log10(f['/cluster_galaxies/cluster_flux_05_24'].value[is_gal]), unit='log10(Xray flux 0.5-2.4 keV [erg/cm2/s])' ) 
  ,fits.Column(name='Mpeak',format='D',    array= f['/halo_properties/Mpeak'].value[is_gal], unit='Msun'  ) 
  ,fits.Column(name='Vmax',format='D',    array= f['/halo_properties/Vmax'].value[is_gal], unit='km/s'  ) 
  ,fits.Column(name='mvir',format='D',    array= f['/halo_properties/mvir'].value[is_gal], unit='Msun'  ) 
  #,fits.Column(name='id',format='K',    array= (f['/halo_properties/id'].value[is_gal]).astype('int'), unit='MDPL ID integer'  ) 
  #,fits.Column(name='pid',format='K',    array= (f['/halo_properties/pid'].value[is_gal]).astype('int'), unit='MDPL PID integer'   ) 
  ,fits.Column(name='log_stellar_mass',format='D', array= n.log10(f['/moster_2013_data/stellar_mass'].value[is_gal]) , unit='log10(stellar_mass/[Msun])'  ) 
  ,fits.Column(name='RA',format='D',    array= f['/sky_position/RA'].value[is_gal] , unit='RA/[deg]'  ) 
  ,fits.Column(name='DEC',format='D',    array= f['/sky_position/DEC'].value[is_gal], unit='DEC/[deg]'   ) 
  ,fits.Column(name='redshift_R',format='D',    array= f['/sky_position/redshift_R'].value[is_gal], unit='real space redshift'  ) 
  ,fits.Column(name='redshift_S',format='D',    array= f['/sky_position/redshift_S'].value[is_gal], unit='redshift space redshift'  ) 
  ])

  f.close()

  tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
  #define the header
  prihdr = fits.Header()
  prihdr['author'] = 'JC'
  prihdr['DEC_max'] = dec_max
  prihdr['DEC_max'] = - dec_max 
  prihdr['RA_max'] = ra_max
  prihdr['RA_max'] = - ra_max
  prihdr['z_min'] = z_min
  prihdr['z_max'] = z_max

  prihdu = fits.PrimaryHDU(header=prihdr)
  #writes the file
  thdulist = fits.HDUList([prihdu, tb_hdu])

  print( out_filename )
  if os.path.isfile(out_filename):
    os.system("rm "+out_filename)
  thdulist.writeto(out_filename)



path_to_lc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_remaped_position_L3.hdf5'
out_filename = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_cluster_remaped_position_L3.fits'
z_min = 0.
z_max = 1.08
dec_max = 8.269819492449505
ra_max = 6.7529257176359
write_fits_lcC(path_to_lc, out_filename, z_min, z_max, dec_max, ra_max)

path_to_lc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_remaped_position_L3_z1.hdf5'
out_filename = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_cluster_remaped_position_L3_z1.fits'
z_min = 1.08
z_max = 3.0
dec_max = 4.134909746242654
ra_max = 3.3764628588325674
write_fits_lcC(path_to_lc, out_filename, z_min, z_max, dec_max, ra_max)

path_to_lc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_remaped_position_L6.hdf5'
out_filename = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_cluster_remaped_position_L6.fits'
z_min = 0.
z_max = 3.0
dec_max = 2.0047373031569915
ra_max = 1.9766516114702513
write_fits_lcC(path_to_lc, out_filename, z_min, z_max, dec_max, ra_max)

path_to_lc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_remaped_position_L15.hdf5'
out_filename = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_cluster_remaped_position_L15.fits'
z_min = 0.
z_max = 0.54
dec_max = 20.257311381848154
ra_max = 14.323944878104827
write_fits_lcC(path_to_lc, out_filename, z_min, z_max, dec_max, ra_max)

# z< 0.5423857379098544 |ra [deg]|< 14.323944878104827 |dec [deg]|< 20.257311381848154

# L3 characteristics :
# z< 1.0889947373832305 |ra [deg]|< 6.7529257176359 |dec [deg]|< 8.269819492449505
# N points: 8037075 
# 
# L3_z1 characteristics
# z< 3.8309961826584344 |ra [deg]|< 3.3764628588325674 |dec [deg]|< 4.134909746242654
# N points: 8511571
# 
# L6 characteristics
# z< 6.697087333514605 |ra [deg]|< 1.9766516114702513 |dec [deg]|< 2.0047373031569915
# N points: 3287299
