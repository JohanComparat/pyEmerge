"""
identifies COSMO tracers
"""
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115)#, Ob0=0.048206)

import glob
import os
import time
import numpy as n
import sys

# specific functions
from scipy.stats import norm
from scipy.integrate import quad
from scipy.interpolate import interp1d
import random
from scipy.stats import norm

import h5py    # HDF5 support



def get_dn_dv(file_name = os.path.join(os.environ['GIT_EMERGE'], "data/NZ/4FS_scenario2.17Jul2017.CoLoRe.BGhiz.nz")):
	print( file_name)
	zmean, dN_dz = n.loadtxt(file_name, unpack=True)
	print(zmean[1:]-zmean[:-1])
	dz_all = zmean[1:]-zmean[:-1]
	dz = dz_all[0] #0.025 # 
	zmin = zmean-dz/2.
	zmax = zmean+dz/2.
	print( n.sum(dN_dz * dz))
	dV_per_deg2 = (cosmoMD.comoving_volume(zmax) - cosmoMD.comoving_volume(zmin))*n.pi/129600.
	dN_dV_data = dN_dz * dz / dV_per_deg2
	dN_dV = interp1d(n.hstack((-1., zmin[0], zmean, zmax[-1], 10.)), n.hstack((0., 0., dN_dV_data, 0., 0.)) )
	return dN_dV

dN_dV_lrg1=get_dn_dv(file_name = os.path.join(os.environ['GIT_EMERGE'], "data/NZ/4FS_scenario2.17Jul2017.CoLoRe.BGlowz.nz"))
dN_dV_lrg2=get_dn_dv(file_name = os.path.join(os.environ['GIT_EMERGE'], "data/NZ/4FS_scenario2.17Jul2017.CoLoRe.BGhiz.nz"))
dN_dV_elg = get_dn_dv(os.path.join(os.environ['GIT_EMERGE'], "data/NZ/4FS_scenario2.17Jul2017.CoLoRe.ELG.nz"))
dN_dV_qso = get_dn_dv(os.path.join(os.environ['GIT_EMERGE'], "data/NZ/4FS_scenario2.17Jul2017.CoLoRe.QSO.nz"))

ii=0
status = 'create'
#status = 'update'

h5_dir = os.path.join(os.environ[env], 'h5' )
input_list = n.array(glob.glob(os.path.join(h5_dir, "hlist_?.?????_emerge.hdf5")))
input_list.sort()

file_1 = input_list[ii]
print(file_1)
f1 = h5py.File(file_1,  "r+")

z_snap = f1.attrs['redshift']
volume = 1000.**3.

all_MS = n.log10(f1['/moster_2013_data/stellar_mass'].value)
rds = n.random.rand(len(all_MS))
# target bits
N_lrg1 = int(volume * dN_dV_lrg1(z_snap)) # 0
N_lrg2 = int(volume * dN_dV_lrg2(z_snap)) # 1
N_elg  = int(volume * dN_dV_elg(z_snap) ) # 2
N_qso  = int(volume * dN_dV_qso(z_snap) ) # 3
#N_Lya  = int(volume * dN_dV_qso(z_snap) ) # 4



if N_lrg1 > 100 :
	# LRG: select on Ms
	
	all_MS_sort_id = n.argsort(all_MS)
	min_mass = all_MS[all_MS_sort_id[-N_lrg1*2-1]]

	mass_selection = (all_MS>min_mass)
	
	lrgs = id_superset[1][ (id_superset[0]==num) & (rds>0.5) ]
	
	for num in list(set(id_superset[0])):
		file_out = os.path.join(os.environ["MD10"], 'work_agn', 'out_'+el['snap_name']+'_SAM_Nb_'+str(num)+'_4MOST_S8_BG1.fits')
		
		print( file_out, lrgs, len(lrgs))
		hdu_cols  = fits.ColDefs([
		fits.Column(name='line_number',format='K', array= lrgs )])
		tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
		prihdr = fits.Header()
		prihdr['HIERARCH nameSnapshot'] = el['snap_name']
		prihdr['batchN'] = num
		prihdr['tracer'] = 'BG1'
		prihdr['author'] = 'JC'
		prihdu = fits.PrimaryHDU(header=prihdr)
		#writes the file
		thdulist = fits.HDUList([prihdu, tb_hdu])
		if os.path.isfile(file_out):
			os.system("rm "+file_out)
		thdulist.writeto(file_out)

if N_lrg2 > 100 :
	# LRG: select on Ms
	MS = n.array([fits.open(file)[1].data['stellar_mass_Mo13_mvir'] for file in fileList_snap_MS])

	all_MS = n.hstack((MS))
	all_MS_sort_id = n.argsort(all_MS)
	min_mass = all_MS[all_MS_sort_id[-N_lrg2*4-1]]

	id_superset = n.where(MS>min_mass)
	rds = n.random.rand(len(id_superset[0]))

	for num in list(set(id_superset[0])):
		file_out = os.path.join(os.environ["MD10"], 'work_agn', 'out_'+el['snap_name']+'_SAM_Nb_'+str(num)+'_4MOST_S8_BG2.fits')
		lrgs = id_superset[1][ (id_superset[0]==num) & (rds>0.75) ]
		print( file_out, lrgs, len(lrgs))
		hdu_cols  = fits.ColDefs([
		fits.Column(name='line_number',format='K', array= lrgs )])
		tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
		prihdr = fits.Header()
		prihdr['HIERARCH nameSnapshot'] = el['snap_name']
		prihdr['batchN'] = num
		prihdr['tracer'] = 'BG2'
		prihdr['author'] = 'JC'
		prihdu = fits.PrimaryHDU(header=prihdr)
		#writes the file
		thdulist = fits.HDUList([prihdu, tb_hdu])
		if os.path.isfile(file_out):
			os.system("rm "+file_out)
		thdulist.writeto(file_out)


if N_elg > 100 :
	MH = n.array([fits.open(file)[1].data['mvir'] for file in fileList_snap])
	# ELG select on Mvir
	#p_elg=[10**(12.2),0.25]
	mh_mean, mh_scatter = 12.2, 0.25
	mh_bins = n.arange(mh_mean -2*mh_scatter, mh_mean +2*mh_scatter+0.1, 0.1)
	mh_bins_pos = 0.5*(mh_bins[1:]+mh_bins[:-1])
	proba = lambda x : norm.pdf(x, loc=mh_mean,scale=mh_scatter)
	proba_norm = proba(mh_bins_pos).sum()
	N_2_select_per_bin = (N_elg*proba(mh_bins_pos)/proba_norm).astype('int')

	id_0 = []
	id_1 = []
	for id_bin in range(len(mh_bins)-1):
		id_superset = n.where( (MH > mh_bins[id_bin]) &( MH < mh_bins[id_bin+1]) )
		N_avail = id_superset[0].shape[0]
		rds = n.random.rand(len(id_superset[0]))
		bin_selection = (rds < N_2_select_per_bin[id_bin]*1./N_avail)
		id_0.append(id_superset[0][bin_selection])
		id_1.append(id_superset[1][bin_selection])

	id_0 = n.hstack((id_0))
	id_1 = n.hstack((id_1))

	for num in list(set(id_0)):
		file_out = os.path.join(os.environ["MD10"], 'work_agn', 'out_'+el['snap_name']+'_SAM_Nb_'+str(num)+'_4MOST_S8_ELG.fits')
		elgs = id_1[ (id_0==num) ]
		print( file_out, elgs, len(elgs))
		hdu_cols  = fits.ColDefs([
		fits.Column(name='line_number',format='K', array= elgs )])
		tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
		prihdr = fits.Header()
		prihdr['HIERARCH nameSnapshot'] = el['snap_name']
		prihdr['batchN'] = num
		prihdr['tracer'] = 'ELG'
		prihdr['author'] = 'JC'
		prihdu = fits.PrimaryHDU(header=prihdr)
		#writes the file
		thdulist = fits.HDUList([prihdu, tb_hdu])
		if os.path.isfile(file_out):
			os.system("rm "+file_out)
		thdulist.writeto(file_out)

if N_qso > 100 :
	MH = n.array([fits.open(file)[1].data['mvir'] for file in fileList_snap])
	# QSO select on Mvir
	#p_qso=[10**(12.7),0.25]
	mh_mean, mh_scatter = 12.7, 0.25
	mh_bins = n.arange(mh_mean -2*mh_scatter, mh_mean +2*mh_scatter+0.1, 0.1)
	mh_bins_pos = 0.5*(mh_bins[1:]+mh_bins[:-1])
	proba_qso = lambda x : norm.pdf(x, loc=mh_mean,scale=mh_scatter)
	proba_norm_qso = proba_qso(mh_bins_pos).sum()
	N_2_select_per_bin = (N_qso*proba_qso(mh_bins_pos)/proba_norm_qso).astype('int')

	id_0 = []
	id_1 = []
	for id_bin in range(len(mh_bins)-1):
		id_superset = n.where( (MH > mh_bins[id_bin]) &( MH < mh_bins[id_bin+1]) )
		N_avail = id_superset[0].shape[0]
		rds = n.random.rand(len(id_superset[0]))
		bin_selection = (rds < N_2_select_per_bin[id_bin]*1./N_avail)
		id_0.append(id_superset[0][bin_selection])
		id_1.append(id_superset[1][bin_selection])

	id_0 = n.hstack((id_0))
	id_1 = n.hstack((id_1))

	for num in list(set(id_0)):
		file_out = os.path.join(os.environ["MD10"], 'work_agn', 'out_'+el['snap_name']+'_SAM_Nb_'+str(num)+'_4MOST_S8_QSO.fits')
		qsos = id_1[ (id_0==num) ]
		print( file_out, qsos, len(qsos))
		hdu_cols  = fits.ColDefs([
		fits.Column(name='line_number',format='K', array= qsos )])
		tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
		prihdr = fits.Header()
		prihdr['HIERARCH nameSnapshot'] = el['snap_name']
		prihdr['batchN'] = num
		prihdr['tracer'] = 'QSO'
		prihdr['author'] = 'JC'
		prihdu = fits.PrimaryHDU(header=prihdr)
		#writes the file
		thdulist = fits.HDUList([prihdu, tb_hdu])
		if os.path.isfile(file_out):
			os.system("rm "+file_out)
		thdulist.writeto(file_out)


if status == 'create' :
  halo_data = fg.create_group('cosmo_4most')
  ds = halo_data.create_dataset('is_BG_lz', data = clus_flux_05_24 )
  ds = halo_data.create_dataset('is_BG_hz', data = clus_flux_05_24 )
  ds = halo_data.create_dataset('is_ELG', data = clus_id )
  ds = halo_data.create_dataset('is_QSO', data = mag_r )
  ds = halo_data.create_dataset('is_LyaQSO', data = mag_r )
  
if status=='update' :
  fg['/cluster_galaxies/cluster_id'][:]  = clus_id 
  fg['/cluster_galaxies/cluster_flux_05_24'][:]  = clus_flux_05_24 
  #fg['/cluster_galaxies/mag_r'][:]  = mag_r 
  fg['/cluster_galaxies/bcg_flag'][:] = bcg_flag 
  fg['/cluster_galaxies/d_cluster_center'][:] = radial_separation

fc.close()
fg.close()