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


"""
import magnitude_library
photo = magnitude_library.Photo()


import h5py    # HDF5 support
import os
import glob
import numpy as n
from scipy.interpolate import interp1d
import astropy.io.fits as fits
import sys
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115)#, Ob0=0.048206)

arch_dir = '/data36s/comparat/archetypes/v2/'


zmins = n.array([0.001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ])
zmaxs = n.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2 ])

gal_type = 'LRG'
Nspec_max = 2000.
sn_min = 2.0

dir_4most = '/data27s/4most/comparat/mocks/'
cat_filename = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_cluster_remaped_position_L3.fits'
out_filename = os.path.join(dir_4most, '4MOST_CLUSTER_L3_r_lt_225.fits')
dec_max = 8.269819492449505
ra_max = 6.7529257176359


hdus =  fits.open(cat_filename)


template_names = n.ones_like(hdus[1].data['redshift_R'] ).astype('|S64')
  
for zmin, zmax in zip(zmins, zmaxs):
  #zmin = zmins[0]
  #zmax = zmaxs[0]
  name = "sdss_"+gal_type+"_zmin_"+str(zmin)+"_zmax_"+str(zmax)+"_Nlt_"+str(Nspec_max)
  path_2_arch = os.path.join(arch_dir, gal_type, "archetypes_"+name+"_snMin"+str(sn_min)+".txt") 
  DATA = n.loadtxt( path_2_arch )
  tpl_name = os.path.basename(path_2_arch)[16:-4]+".fits"
  out_template_name = os.path.join(dir_4most, 'templates', tpl_name)

  wl_i = DATA[0]*(1+0.5*(zmin+zmax))
  fl_i = DATA[1]*1e-17

  wl_j = n.hstack((1000., wl_i[0]-1, wl_i, wl_i[-1]+1, 30000))
  fl_j = n.hstack((0., 0., fl_i, 0., 0.))

  distMod = 0. #cosmoMD.distmod( 0.5*(zmin+zmax)).value 
  observed_magnitudes, arr = photo.computeMagnitudes(interp1d(wl_j, fl_j), distMod) 
  rescale_by = 10**((14. + 48.6)/-2.5) / 10**((observed_magnitudes[2] + 48.6)/-2.5) 
  normed_magnitudes, arr = photo.computeMagnitudes(interp1d(wl_j, fl_j*rescale_by), distMod) 
  print(normed_magnitudes)
  # writes template

  hdu_cols  = fits.ColDefs([fits.Column(name='LAMBDA', format='D', array=wl_j, unit='Angstrom'), fits.Column(name='FLUX', format='D', array=fl_j*rescale_by, unit='erg/cm2/s/A')])

  tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
  #define the header
  prihdr = fits.Header()
  prihdr['author'] = 'JC'
  prihdr['ABMAG'] = n.round(normed_magnitudes[2],1)
  prihdu = fits.PrimaryHDU(header=prihdr)
  #writes the file
  thdulist = fits.HDUList([prihdu, tb_hdu])

  print( out_template_name )
  if os.path.isfile(out_template_name):
    os.system("rm "+out_template_name)
    
  thdulist.writeto(out_template_name)

  gal = (hdus[1].data['redshift_R'] > n.round(zmin,1) )& ( hdus[1].data['redshift_R'] <= n.round(zmax,1) )
  template_names[gal] = tpl_name


ruleset = n.array(['Clusters' for ii in n.arange(len(hdus[1].data['redshift_R']))])

bright = (hdus[1].data['mag_r']<22.5)

hdu_cols  = fits.ColDefs([
fits.Column(name='OBJECT_ID', format='15A', array=n.arange(len(hdus[1].data['redshift_R'])).astype('str')[bright] )
,fits.Column(name='IDNUM', format='1J', array=n.arange(len(hdus[1].data['redshift_R']))[bright] )
,fits.Column(name='RA', format='1D', array=hdus[1].data['RA'][bright] )
,fits.Column(name='DEC', format='1D', array=hdus[1].data['DEC'][bright] )
,fits.Column(name='PRIORITY', format='1I', array=n.ones_like(hdus[1].data['redshift_R']).astype('int')[bright]*100 )
#,fits.Column(name='ORIG_TEXP_D', format='1E', array=n.ones_like(hdus[1].data['redshift_R']).astype('int')[bright]*10.  )
#,fits.Column(name='ORIG_TEXP_G', format='1E', array=n.ones_like(hdus[1].data['redshift_R']).astype('int')[bright]*10.  )
#,fits.Column(name='ORIG_TEXP_B', format='1E', array=n.ones_like(hdus[1].data['redshift_R']).astype('int')[bright]*10.  )
,fits.Column(name='RESOLUTION', format='1B', array=n.ones_like(hdus[1].data['redshift_R']).astype('int')[bright] )
,fits.Column(name='R_MAG', format='1E', array=hdus[1].data['mag_r'][bright] )
,fits.Column(name='TEMPLATE', format='64A', array=template_names [bright]) 
,fits.Column(name='REDSHIFT', format='1E', array=hdus[1].data['redshift_S'][bright] )
,fits.Column(name='RULESET', format='9A', array=ruleset[bright] )
,fits.Column(name='SUB_CAT_ID', format='1J', array=hdus[1].data['cluster_id'][bright] )
,fits.Column(name='cluster_log_xray_flux_05_20'  , format='1D', array=hdus[1].data['cluster_log_xray_flux_05_20' ][bright] )
,fits.Column(name='Vmax'                         , format='1D', array=hdus[1].data['Vmax'                        ][bright] )
,fits.Column(name='mvir'                         , format='1D', array=hdus[1].data['mvir'                        ][bright] )
,fits.Column(name='log_stellar_mass'             , format='1D', array=hdus[1].data['log_stellar_mass'            ][bright] )
,fits.Column(name='d_cluster_center'             , format='1D', array=hdus[1].data['d_cluster_center'            ][bright] )
,fits.Column(name='is_agn'                       , format='1D', array=hdus[1].data['is_agn'                      ][bright] )
,fits.Column(name='agn_logNH'                    , format='1D', array=hdus[1].data['agn_logNH'                   ][bright] )
,fits.Column(name='agn_log_lambda_sar'           , format='1D', array=hdus[1].data['agn_log_lambda_sar'          ][bright] )
,fits.Column(name='agn_log_xray_flux_05_20'      , format='1D', array=hdus[1].data['agn_log_xray_flux_05_20'     ][bright] )
,fits.Column(name='agn_log_xray_luminosity_2_10' , format='1D', array=hdus[1].data['agn_log_xray_luminosity_2_10'][bright] )
])

tb_hdu = fits.BinTableHDU.from_columns( hdu_cols  )
#define the header
prihdr = fits.Header()
prihdr['author'] = 'JC'

prihdu = fits.PrimaryHDU(header=prihdr)
#writes the file
thdulist = fits.HDUList([prihdu, tb_hdu])

print( out_filename )
if os.path.isfile(out_filename):
  os.system("rm "+out_filename)
thdulist.writeto(out_filename)

sys.exit()

def read_template_columns(filename):
	h = fits.open(filename)
	return os.path.basename(filename), h[1].header['ABMAG'], h[1].data['LAMBDA'], h[1].data['FLUX']


# shift to observed frame at mean redshift
# norm r mag 14 
# zero-padded left and right
# integrate r mag and write template in the 4MOST format  
  

def write_fits_lcC(path_to_lc, out_filename, z_min, z_max, dec_max, ra_max):
  f = h5py.File(path_to_lc, 'r')

  is_gal = (f['/sky_position/selection'].value)&(f['/cluster_galaxies/cluster_id'].value > 0 ) & (f['/cluster_galaxies/cluster_flux_05_24'].value> 0) & (f['/sky_position/redshift_R'].value<z_max)

  hdu_cols  = fits.ColDefs([
  fits.Column(name='cluster_log_xray_flux_05_20',format='D',    array= n.log10(f['/cluster_galaxies/cluster_flux_05_24'].value[is_gal]), unit='cluster log10(Xray flux 0.5-2 keV [erg/cm2/s])' ) 
  ,fits.Column(name='cluster_id',format='K',    array= f['/cluster_galaxies/cluster_id'].value[is_gal], unit='cluster unique identifier' ) 
  #,fits.Column(name='id',format='K',    array= f['/halo_properties/id'].value[is_gal], unit='halo unique id'  ) 
  ,fits.Column(name='is_subhalo',format='L',    array= (f['/halo_properties/pid'].value[is_gal]==-1), unit='True if galaxy is in a subhalo'  ) 
  ,fits.Column(name='mvir_dot',format='D',    array= f['/emerge_data/mvir_dot'].value[is_gal], unit='Msun/yr'  ) 
  ,fits.Column(name='dMdt',format='D',    array= f['/emerge_data/dMdt'].value[is_gal], unit='Msun/yr'  ) 
  ,fits.Column(name='Mpeak',format='D',    array= f['/halo_properties/Mpeak'].value[is_gal], unit='Msun'  ) 
  ,fits.Column(name='Vmax',format='D',    array= f['/halo_properties/Vmax'].value[is_gal], unit='km/s'  ) 
  ,fits.Column(name='mvir',format='D',    array= f['/halo_properties/mvir'].value[is_gal], unit='Msun'  ) 
  ,fits.Column(name='log_stellar_mass',format='D', array= n.log10(f['/moster_2013_data/stellar_mass'].value[is_gal]) , unit='log10(stellar_mass/[Msun])'  ) 
  ,fits.Column(name='RA',format='D',    array= f['/sky_position/RA'].value[is_gal] , unit='RA/[deg]'  ) 
  ,fits.Column(name='DEC',format='D',    array= f['/sky_position/DEC'].value[is_gal], unit='DEC/[deg]'   ) 
  ,fits.Column(name='redshift_R',format='D',    array= f['/sky_position/redshift_R'].value[is_gal], unit='real space redshift'  ) 
  ,fits.Column(name='redshift_S',format='D',    array= f['/sky_position/redshift_S'].value[is_gal], unit='redshift space redshift'  ) 
  ,fits.Column(name='mag_r',format='D',    array= f['/cluster_galaxies/mag_r'].value[is_gal], unit='sdss r band magnitude for cluster members' ) 
  ,fits.Column(name='is_bcg',format='L',    array= (f['/cluster_galaxies/bcg_flag'].value[is_gal]==1), unit='True if galaxy is the BCG'  ) 
  ,fits.Column(name='d_cluster_center',format='D',    array= f['/cluster_galaxies/d_cluster_center'].value[is_gal], unit='distance to cluster most massive halo center [rvir]'  ) 
  ,fits.Column(name='is_agn',format='L',    array= (f['/agn_properties/agn_activity'].value[is_gal]==1), unit='True if agn' ) 
  ,fits.Column(name='agn_logNH',format='D',    array= f['/agn_properties/logNH'].value[is_gal], unit='log10(NH/[cm2])' ) 
  ,fits.Column(name='agn_log_lambda_sar',format='D',    array= f['/agn_properties/log_lambda_sar'].value[is_gal], unit='log10(LSAR/[erg/s/Msun])' ) 
  ,fits.Column(name='agn_log_xray_flux_05_20',format='D',    array= n.log10(f['/agn_properties/rxay_flux_05_20'].value[is_gal]), unit='log10(Xray flux 0.5-2 keV [erg/cm2/s])' ) 
  ,fits.Column(name='agn_log_xray_luminosity_2_10',format='D',    array= f['/agn_properties/log_lambda_sar'].value[is_gal] + n.log10(f['/moster_2013_data/stellar_mass'].value[is_gal]), unit='log10(LX [2-10keV] / erg/s])' ) 
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

path_to_lc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_agn_Bongiorno_Ricci_Zevol/lc_remaped_position_L3_z1.hdf5'
out_filename = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_cluster_remaped_position_L3_z1.fits'
z_min = 1.08
z_max = 3.0
dec_max = 4.134909746242654
ra_max = 3.3764628588325674
write_fits_lcC(path_to_lc, out_filename, z_min, z_max, dec_max, ra_max)

path_to_lc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_agn_Bongiorno_Ricci_Zevol/lc_remaped_position_L3.hdf5'
out_filename = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_cluster_remaped_position_L3.fits'
z_min = 0.
z_max = 1.08
dec_max = 8.269819492449505
ra_max = 6.7529257176359
write_fits_lcC(path_to_lc, out_filename, z_min, z_max, dec_max, ra_max)

path_to_lc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_agn_Bongiorno_Ricci_Zevol/lc_remaped_position_L6.hdf5'
out_filename = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_cluster_remaped_position_L6.fits'
z_min = 0.
z_max = 3.0
dec_max = 2.0047373031569915
ra_max = 1.9766516114702513
write_fits_lcC(path_to_lc, out_filename, z_min, z_max, dec_max, ra_max)

path_to_lc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_agn_Bongiorno_Ricci_Zevol/lc_remaped_position_L15.hdf5'
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

#log_xray_flux_05_24 > n.log10(1.5e-14)

#is_gal = (f['/sky_position/selection'].value)&(f['/cluster_galaxies/cluster_id'].value > 0 ) & (f['/cluster_galaxies/cluster_flux_05_24'].value> 0) &  (f['/cluster_galaxies/cluster_flux_05_24'].value > 1.5e-14 ) & (f['/cluster_galaxies/mag_r'].value<23)

#cluIDS = list(set(f['/cluster_galaxies/cluster_id'].value[is_gal]))

#hdu_cols  = fits.ColDefs([
#fits.Column(name='cluster_log_xray_flux_05_20',format='D',    array= n.log10(f['/cluster_galaxies/cluster_flux_05_24'].value[is_gal]), unit='cluster log10(Xray flux 0.5-2 keV [erg/cm2/s])' ) 
#,fits.Column(name='cluster_id',format='K',    array= n.log10(f['/cluster_galaxies/cluster_id'].value[is_gal]), unit='cluster unique identifier' ) 
#,fits.Column(name='Mpeak',format='D',    array= f['/halo_properties/Mpeak'].value[is_gal], unit='Msun'  ) 
#,fits.Column(name='Vmax',format='D',    array= f['/halo_properties/Vmax'].value[is_gal], unit='km/s'  ) 
#,fits.Column(name='mvir',format='D',    array= f['/halo_properties/mvir'].value[is_gal], unit='Msun'  ) 
#,fits.Column(name='log_stellar_mass',format='D', array= n.log10(f['/moster_2013_data/stellar_mass'].value[is_gal]) , unit='log10(stellar_mass/[Msun])'  ) 
#,fits.Column(name='RA',format='D',    array= f['/sky_position/RA'].value[is_gal] , unit='RA/[deg]'  ) 
#,fits.Column(name='DEC',format='D',    array= f['/sky_position/DEC'].value[is_gal], unit='DEC/[deg]'   ) 
#,fits.Column(name='redshift_R',format='D',    array= f['/sky_position/redshift_R'].value[is_gal], unit='real space redshift'  ) 
#,fits.Column(name='redshift_S',format='D',    array= f['/sky_position/redshift_S'].value[is_gal], unit='redshift space redshift'  ) 
#,fits.Column(name='mag_r',format='D',    array= f['/cluster_galaxies/mag_r'].value[is_gal], unit='sdss r band magnitude for cluster members' ) 
#,fits.Column(name='is_agn',format='L',    array= (f['/agn_properties/agn_activity'].value[is_gal]==1), unit='True if agn' ) 
#,fits.Column(name='agn_logNH',format='D',    array= f['/agn_properties/logNH'].value[is_gal], unit='log10(NH/[cm2])' ) 
#,fits.Column(name='agn_log_lambda_sar',format='D',    array= f['/agn_properties/log_lambda_sar'].value[is_gal], unit='log10(LSAR/[erg/s/Msun])' ) 
#,fits.Column(name='agn_log_xray_flux_05_20',format='D',    array= n.log10(f['/agn_properties/rxay_flux_05_20'].value[is_gal]), unit='log10(Xray flux 0.5-2 keV [erg/cm2/s])' ) 
#,fits.Column(name='agn_log_xray_luminosity_2_10',format='D',    array= f['/agn_properties/log_lambda_sar'].value[is_gal] + n.log10(f['/moster_2013_data/stellar_mass'].value[is_gal]), unit='log10(LX [2-10keV] / erg/s])' ) 
#])