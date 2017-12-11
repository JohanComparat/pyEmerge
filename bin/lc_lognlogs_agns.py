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

plotDir = os.path.join(os.environ['HOME'], 'wwwDir', "eRoMok", "logNlogS")


from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)

def get_lognlogs(path_to_lc, area, z_max=3.):
  f = h5py.File(path_to_lc, 'r+')
  is_gal = (f['/sky_position/selection'].value)&(f['/sky_position/redshift_R'].value<z_max)
  is_agn = (f['/sky_position/selection'].value)&(f['/agn_properties/agn_activity'].value==1)&(f['/agn_properties/rxay_flux_05_20'].value>0)
  n_gal = len(f['/sky_position/redshift_S'].value[is_gal])
  n_agn = len(f['/sky_position/redshift_S'].value[is_agn])
  z = f['/sky_position/redshift_S'].value[is_agn]
  #logm = n.log10(f['/moster_2013_data/stellar_mass'].value[is_agn])
  #lsar = f['/agn_properties/log_lambda_sar'].value[is_agn]
  #lx = logm + lsar
  log_f_05_20 = n.log10(f['/agn_properties/rxay_flux_05_20'].value[is_agn]) #- 0.6
  f.close()
  out = n.histogram(log_f_05_20, bins = n.arange(-18, -8., 0.2))
  # cumulative number density per square degrees
  x_out = 0.5*(out[1][1:] + out[1][:-1])
  N_out = n.array([n.sum(out[0][ii:]) for ii in range(len(out[0])) ])
  c_out = n.array([n.sum(out[0][ii:]) for ii in range(len(out[0])) ]) / area
  c_out_up = (1 + N_out**(-0.5)) * c_out
  c_out_low = (1 - N_out**(-0.5)) * c_out
  c_err = (n.log10(c_out_up) - n.log10(c_out_low))/2.
  return x_out, c_out, c_err

p.figure(1, (6,6))

path_to_lc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_L3.hdf5'
area = 6.7529257176359*2. * 2* 8.269819492449505
x_out, c_out, c_err = get_lognlogs(path_to_lc, area, z_max=1.1)
#p.plot(x_out, n.log10(c_out), lw=2, rasterized = True, label = 'z<1.08' )
p.errorbar(x_out, n.log10(c_out), yerr = c_err, rasterized = True, label = 'L3 z<1.08, 223deg2'  )
x_out_a, c_out_a, c_err_a = x_out, c_out, c_err 


#path_to_lc=='/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_remaped_position_L3_z1.hdf5'
#area = 3.3764628588325674*2. * 2* 4.134909746242654
#x_out, c_out, c_err = get_lognlogs(path_to_lc, area, z_max=3.)
#p.errorbar(x_out, n.log10(c_out), yerr = c_err, rasterized = True, label = 'L3 1.08<z<3.' )
#p.plot(x_out, n.log10(c_out+c_out_a), ls='dashed', label='total')

path_to_lc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_L6.hdf5'
area = 1.9766516114702513*2. * 2*2.0047373031569915
x_out, c_out, c_err = get_lognlogs(path_to_lc, area, z_max=3.)
p.errorbar(x_out, n.log10(c_out), yerr = c_err, rasterized = True, label = 'L6 z<3., 15deg2'  )

#p.plot(x_out-0.1, n.log10(c_out), 'k', lw=2, rasterized = True, label = 'L3 lc-0.1' )
#p.plot(x_out, n.log10(c_out*(1-frac_err_13deg2)), 'k--', lw=1, rasterized = True, label = 'v0.6, 13.3deg2 scatter' )
#p.plot(x_out, n.log10(c_out*(1+frac_err_13deg2)), 'k--', lw=1, rasterized = True)
#p.plot(x_out, n.log10(c_out*(1-frac_err_3deg2)), 'r--', lw=1, rasterized = True, label = 'v0.6, 3.5deg2 scatter' )
#p.plot(x_out, n.log10(c_out*(1+frac_err_3deg2)), 'r--', lw=1, rasterized = True)
#p.plot(x_out_0, n.log10(c_out_0), 'm--', rasterized = True, label = 'Planck mock v0.0' )

path_to_lc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_L15.hdf5'
area = 14.323944878104827*2. * 2*20.257311381848154
x_out, c_out, c_err = get_lognlogs(path_to_lc, area, z_max=3.)
p.errorbar(x_out, n.log10(c_out), yerr = c_err, rasterized = True, label = 'L15 z<0.54 1160deg2'  )

path_2_logNlogS_data = os.path.join(os.environ["DARKSIM_DIR"], 'observations', 'logNlogS', 'logNlogS_Georgakakis_08_AGN.data')
x_data, y_data, yerr = n.loadtxt(path_2_logNlogS_data, unpack=True)
p.fill_between(x_data, y1 = n.log10(y_data-yerr), y2=n.log10(y_data+yerr), color='b' , rasterized = True, alpha=0.5, label = 'Georgakakis 08' )
#p.plot(x_data, n.log10(y_data))
path_2_logNlogS_data = os.path.join(os.environ["DARKSIM_DIR"], 'observations', 'logNlogS', 'logNlogS_Merloni_12_AGN.data')
x_data, y_data  = n.loadtxt(path_2_logNlogS_data, unpack=True)
p.plot(x_data, n.log10(y_data), label = 'Merloni 12' )

p.axhline(7, ls='dashed')
p.xlabel('log(F[0.5-2 keV])')
p.ylabel('log(>F) [/deg2]')
p.legend(frameon=False, loc=0)
#p.yscale('log')
p.xlim((-17, -12))
p.ylim((-2, 4.))
#p.title('Mocks')
p.grid()
p.savefig(os.path.join(plotDir, "logN_logS_AGN.jpg"))
p.clf()


