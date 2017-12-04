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
  is_gal = (f['/sky_position/selection'].value)&(f['/sky_position/redshift_R'].value<z_max)&(f['/agn_properties/rxay_flux_05_20'].value>0)
  n_gal = len(f['/sky_position/redshift_S'].value[is_gal])
  z = f['/sky_position/redshift_S'].value[is_gal]
  log_f_05_20 = n.log10(f['/cluster_data/rxay_flux_05_24'].value[is_gal]) #- 0.6
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

path_to_lc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_cluster_remaped_position_L3.hdf5'
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

path_to_lc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_cluster_remaped_position_L6.hdf5'
area = 1.9766516114702513*2. * 2*2.0047373031569915
x_out, c_out, c_err = get_lognlogs(path_to_lc, area, z_max=3.)
p.errorbar(x_out, n.log10(c_out), yerr = c_err, rasterized = True, label = 'L6 z<3., 15deg2'  )

#p.plot(x_out-0.1, n.log10(c_out), 'k', lw=2, rasterized = True, label = 'L3 lc-0.1' )
#p.plot(x_out, n.log10(c_out*(1-frac_err_13deg2)), 'k--', lw=1, rasterized = True, label = 'v0.6, 13.3deg2 scatter' )
#p.plot(x_out, n.log10(c_out*(1+frac_err_13deg2)), 'k--', lw=1, rasterized = True)
#p.plot(x_out, n.log10(c_out*(1-frac_err_3deg2)), 'r--', lw=1, rasterized = True, label = 'v0.6, 3.5deg2 scatter' )
#p.plot(x_out, n.log10(c_out*(1+frac_err_3deg2)), 'r--', lw=1, rasterized = True)
#p.plot(x_out_0, n.log10(c_out_0), 'm--', rasterized = True, label = 'Planck mock v0.0' )

path_to_lc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_cluster_remaped_position_L15.hdf5'
area = 14.323944878104827*2. * 2*20.257311381848154
x_out, c_out, c_err = get_lognlogs(path_to_lc, area, z_max=3.)
p.errorbar(x_out, n.log10(c_out), yerr = c_err, rasterized = True, label = 'L15 z<0.54 1160deg2'  )

path_2_logNlogS_data = os.path.join(os.environ["DARKSIM_DIR"], 'observations', 'logNlogS', 'logNlogS_Finoguenov_cosmos_2007_clusters.data')
x_data, y_data, y_data_min, y_data_max = n.loadtxt(path_2_logNlogS_data, unpack=True)
p.fill_between(n.log10(x_data), y1 = n.log10(y_data_min), y2=n.log10(y_data_max), color='b' , rasterized = True, alpha=0.5, label = 'F07 COSMOS' )
#p.plot(x_data, n.log10(y_data))
path_2_logNlogS_data = os.path.join(os.environ["DARKSIM_DIR"], 'observations', 'logNlogS', 'logNlogS_Finoguenov_ecdfs_2015_clusters.data')
x_data, y_data, y_data_min, y_data_max  = n.loadtxt(path_2_logNlogS_data, unpack=True)
p.fill_between(n.log10(x_data), y1 = n.log10(y_data_min), y2=n.log10(y_data_max), color='r' , rasterized = True, alpha=0.5, label = 'F15 CDFS' )

p.axhline(7, ls='dashed')
p.xlabel('log(F[0.5-2 keV])')
p.ylabel('log(>F) [/deg2]')
p.legend(frameon=False, loc=0)
#p.yscale('log')
p.xlim((-17, -12))
p.ylim((-2, 4.))
#p.title('Mocks')
p.grid()
p.savefig(os.path.join(plotDir, "logN_logS_CLUSTERS.jpg"))
p.clf()


