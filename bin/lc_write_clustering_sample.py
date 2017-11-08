"""
PYTHON 2

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

f = h5py.File('/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_remaped_position_L3_.hdf5', 'r')

is_gal = (f['/sky_position/selection'].value)
is_agn = (f['/sky_position/selection'].value)&(f['/agn_properties/agn_activity'].value==1)

n_gal = len(f['/sky_position/redshift_S'].value[is_gal])

n_agn = len(f['/sky_position/redshift_S'].value[is_agn])

z = f['/sky_position/redshift_S'].value[is_agn]
logm = n.log10(f['/moster_2013_data/stellar_mass'].value[is_agn])
lsar = f['/agn_properties/log_lambda_sar'].value[is_agn]
lx = logm + lsar

log_f_05_20 = n.log10(f['/agn_properties/rxay_flux_05_20'].value) 

area = 6.7529257176359*2. * 2* 8.269819492449505
topdir = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_remaped_position_L3/'

raR, decR = n.loadtxt(topdir + 'random-ra-dec.txt', unpack=True)

def write_samp(zmax,lxmin, out_name='lc_remaped_position_L3_z_lt_03_lx_gt_438.ascii'):
	zmin=0.08
	sel = (is_agn)&(f['/sky_position/redshift_S'].value>zmin)&(f['/sky_position/redshift_S'].value<zmax)&(n.log10(f['/moster_2013_data/stellar_mass'].value)+f['/agn_properties/log_lambda_sar'].value>lxmin)
	
	n.savetxt(out_name, n.transpose([f['/sky_position/RA'].value[sel], f['/sky_position/DEC'].value[sel], f['/sky_position/redshift_S'].value[sel], n.ones_like(f['/sky_position/redshift_S'].value[sel])]) )
	print(zmax, lxmin, len(f['/sky_position/RA'].value[sel]))


	N_data = len(f['/sky_position/RA'].value[sel]) 
	N_rds = len(raR) 
	print("D,R=",N_data, N_rds)
	dz=0.05
	zs=np.arange(zmin, zmax + dz, dz)
	nn,bb = np.histogram(f['/sky_position/redshift_S'].value[sel], bins=zs)#, weights=1./w_col.array)
	nz=interp1d((zs[1:]+zs[:-1])/2.,nn)
	rdsz=[]
	for i in range(1,len(zs)-1,1):
		inter=np.random.uniform(low=zs[i]-dz/2., high=zs[i]+dz/2., size=int( 1000* nz( zs[i] )))
		rdsz.append(inter)

	rds=np.hstack((rdsz))
	np.random.shuffle(rds)
	RR=rds[:N_rds]#-dz/2.
	print("RR=",len(rds), len(RR))

	n.savetxt(out_name[:-4]+'random', n.transpose([raR, decR, RR, np.ones_like(RR) ]))

	p.figure(1, (6,6))
	p.plot(f['/sky_position/redshift_S'].value[sel], n.log10(f['/halo_properties/mvir'].value[sel]), 'k,', rasterized = True )
	p.axvline(0.08, ls='dashed')
	p.ylabel('mvir')
	p.xlabel('redshift')
	p.legend(frameon=False, loc=0)
	#p.yscale('log')
	p.xlim((0,1.2))
	#p.ylim((40, 46))
	p.title('200deg2 mock')
	p.grid()
	p.savefig(os.path.join(plotDir, "HOD_z_"+str(zmax)+"_lx_"+str(lxmin)+".jpg"))
	p.clf()
	return sel

#p.figure(1, (6,6))
#p.plot(f['/sky_position/redshift_S'].value[sel], n.log10(f['/halo_properties/mvir'].value[sel]), 'k,', rasterized = True )
#p.axvline(0.08, ls='dashed')
#p.ylabel('mvir')
#p.xlabel('redshift')
#p.legend(frameon=False, loc=0)
##p.yscale('log')
#p.xlim((0,1.2))
#p.ylim((40, 46))
#p.title('200deg2 mock')
#p.grid()
#p.savefig(os.path.join(plotDir, "HOD_z_"+str(zmax)+"_lx_"+str(lxmin)+".jpg"))
#p.clf()

sel = write_samp(0.3, 44.0, out_name=topdir+'lc_L3_z_lt_03_lx_gt_440.ascii')
sel = write_samp(0.3, 43.5, out_name=topdir+'lc_L3_z_lt_03_lx_gt_435.ascii')
sel = write_samp(0.3, 43., out_name=topdir+'lc_L3_z_lt_03_lx_gt_430.ascii')
sel = write_samp(0.3, 42.5, out_name=topdir+'lc_L3_z_lt_03_lx_gt_425.ascii')
sel = write_samp(0.3, 42., out_name=topdir+'lc_L3_z_lt_03_lx_gt_420.ascii')
sel = write_samp(0.3, 41.5, out_name=topdir+'lc_L3_z_lt_03_lx_gt_415.ascii')

sel = write_samp(0.4, 44., out_name=topdir+'lc_L3_z_lt_04_lx_gt_440.ascii')
sel = write_samp(0.4, 43., out_name=topdir+'lc_L3_z_lt_04_lx_gt_430.ascii')
sel = write_samp(0.4, 42.5, out_name=topdir+'lc_L3_z_lt_04_lx_gt_425.ascii')
sel = write_samp(0.4, 42., out_name=topdir+'lc_L3_z_lt_04_lx_gt_420.ascii')
sel = write_samp(0.4, 41.5, out_name=topdir+'lc_L3_z_lt_04_lx_gt_415.ascii')

# create a mangle mask and make randoms in it	
#  6.7529257176359*2. * 2* 8.269819492449505
import pymangle
'/data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_remaped_position_L3/mask.ply'

p.figure(1, (6,6))
p.plot(z, lx, 'k,', rasterized = True )
p.plot(z[log_f_05_20>-12.7], lx[log_f_05_20>-12.7], 'r+', rasterized = True )
p.axvline(0.08, ls='dashed')
p.ylabel('log(LX)')
p.xlabel('redshift')
p.legend(frameon=False, loc=0)
#p.yscale('log')
p.xlim((0,1.2))
p.ylim((40, 46))
p.title('200deg2 mock')
p.grid()
p.savefig(os.path.join(plotDir, "z_lx_AGN.jpg"))
p.clf()

