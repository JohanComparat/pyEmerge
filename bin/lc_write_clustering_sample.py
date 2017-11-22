"""
Writes clustering samples: ra, dec, z for a set of LX cuts
Based on the MDPL lightcones.
"""
import h5py    # HDF5 support
import os
import glob
import numpy as n
from scipy.interpolate import interp1d

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p



from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115)

Lname='L3'

def write_samples(Lname):
	path_2_lc = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_remaped_position_'+Lname+'.hdf5'
	topdir = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_remaped_position_'+Lname+'/'
	plotDir = os.path.join(os.environ['HOME'], 'wwwDir', "eRoMok", "h5", "clustering_AGN", Lname)
	if os.path.isdir(plotDir)==False:
		os.mkdir(plotDir)

	f = h5py.File(path_2_lc, 'r')

	is_gal = (f['/sky_position/selection'].value)
	is_agn = (f['/sky_position/selection'].value)&(f['/agn_properties/agn_activity'].value==1)

	n_gal = len(f['/sky_position/redshift_S'].value[is_gal])

	n_agn = len(f['/sky_position/redshift_S'].value[is_agn])

	z = f['/sky_position/redshift_S'].value[is_agn]
	logm = n.log10(f['/moster_2013_data/stellar_mass'].value[is_agn])
	lsar = f['/agn_properties/log_lambda_sar'].value[is_agn]
	lx = logm + lsar

	log_f_05_20 = n.log10(f['/agn_properties/rxay_flux_05_20'].value) 

	raR, decR = n.loadtxt(topdir + 'random-ra-dec.txt', unpack=True)

	def write_samp(zmax,lxmin, out_name = 'lc_remaped_position_'+Lname+'_z_lt_03_lx_gt_438.ascii'):
		zmin=0.001
		sel = (is_agn)&(f['/sky_position/redshift_S'].value>zmin)&(f['/sky_position/redshift_S'].value<zmax)&(n.log10(f['/moster_2013_data/stellar_mass'].value)+f['/agn_properties/log_lambda_sar'].value>lxmin)
		
		n.savetxt(out_name, n.transpose([f['/sky_position/RA'].value[sel], f['/sky_position/DEC'].value[sel], f['/sky_position/redshift_S'].value[sel], n.ones_like(f['/sky_position/redshift_S'].value[sel])]) )
		print(zmax, lxmin, len(f['/sky_position/RA'].value[sel]))


		N_data = len(f['/sky_position/RA'].value[sel]) 
		N_rds = 20*N_data # len(raR) 
		print("D,R=",N_data, N_rds)
		dz=0.05
		zs=n.arange(zmin, zmax + dz, dz)
		nn,bb = n.histogram(f['/sky_position/redshift_S'].value[sel], bins=zs)#, weights=1./w_col.array)
		nz=interp1d((zs[1:]+zs[:-1])/2.,nn)
		rdsz=[]
		for i in range(1,len(zs)-1,1):
			inter=n.random.uniform(low=zs[i]-dz/2., high=zs[i]+dz/2., size=int( 1000* nz( zs[i] )))
			rdsz.append(inter)

		rds=n.hstack((rdsz))
		n.random.shuffle(rds)
		RR=rds[:N_rds]#-dz/2.
		print("RR=",len(rds), len(RR))
		n.random.shuffle(raR)
		n.random.shuffle(decR)
		
		n.savetxt(out_name[:-5]+'random', n.transpose([raR[:N_rds], decR[:N_rds], RR, n.ones_like(RR) ]))

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

	sel = write_samp(0.3, 44.0, out_name=topdir+'lc_'+Lname+'_z_lt_03_lx_gt_440.ascii')
	sel = write_samp(0.3, 43.5, out_name=topdir+'lc_'+Lname+'_z_lt_03_lx_gt_435.ascii')
	sel = write_samp(0.3, 43., out_name=topdir+'lc_'+Lname+'_z_lt_03_lx_gt_430.ascii')
	sel = write_samp(0.3, 42.5, out_name=topdir+'lc_'+Lname+'_z_lt_03_lx_gt_425.ascii')
	sel = write_samp(0.3, 42., out_name=topdir+'lc_'+Lname+'_z_lt_03_lx_gt_420.ascii')
	sel = write_samp(0.3, 41.5, out_name=topdir+'lc_'+Lname+'_z_lt_03_lx_gt_415.ascii')

	sel = write_samp(0.4, 44., out_name=topdir+'lc_'+Lname+'_z_lt_04_lx_gt_440.ascii')
	sel = write_samp(0.4, 43.5, out_name=topdir+'lc_'+Lname+'_z_lt_04_lx_gt_435.ascii')
	sel = write_samp(0.4, 43., out_name=topdir+'lc_'+Lname+'_z_lt_04_lx_gt_430.ascii')
	sel = write_samp(0.4, 42.5, out_name=topdir+'lc_'+Lname+'_z_lt_04_lx_gt_425.ascii')
	sel = write_samp(0.4, 42., out_name=topdir+'lc_'+Lname+'_z_lt_04_lx_gt_420.ascii')
	sel = write_samp(0.4, 41.5, out_name=topdir+'lc_'+Lname+'_z_lt_04_lx_gt_415.ascii')


	#p.figure(1, (6,6))
	#p.plot(z, lx, 'k,', rasterized = True )
	#p.plot(z[log_f_05_20>-12.7], lx[log_f_05_20>-12.7], 'r+', rasterized = True )
	##p.axvline(0.08, ls='dashed')
	#p.ylabel('log(LX)')
	#p.xlabel('redshift')
	#p.legend(frameon=False, loc=0)
	##p.yscale('log')
	#p.xlim((0,1.2))
	#p.ylim((40, 46))
	#p.title('200deg2 mock')
	#p.grid()
	#p.savefig(os.path.join(plotDir, Lname+"_z_lx_AGN.jpg"))
	#p.clf()

#write_samples("L3")
write_samples("L6")
write_samples("L15")
