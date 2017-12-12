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
	#print(zmean[1:]-zmean[:-1])
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
EBE_zmin, EBE_zmax, EBE_SGC_long, EBE_SGC_nom, EBE_NGC = n.loadtxt(os.path.join(os.environ['GIT_EMERGE'], "data/NZ/eBOSS-ELG-Raichoor-17.txt"), unpack=True)
dN_dV_elg2 = interp1d( n.hstack((-1., EBE_zmin[0], (EBE_zmin+ EBE_zmax)*0.5, EBE_zmax[-1])), n.hstack((0., 0., 0.7**3.*EBE_SGC_long*0.0001, 0.)) )

ii = int(sys.argv[1]) # -10
env = sys.argv[2] # 'MD10'
status = sys.argv[3] # 'create'
#status = 'update'

h5_dir = os.path.join(os.environ[env], 'h5' )        
input_list = n.array(glob.glob(os.path.join(h5_dir, "hlist_?.?????_emerge.hdf5")))
input_list.sort()

file_1 = input_list[ii]
print(file_1)
f1 = h5py.File(file_1,  "r+")

z_snap = f1.attrs['redshift']
volume = 1000.**3./0.6777**3.

all_MS = n.log10(f1['/moster_2013_data/stellar_mass'].value)
all_mvir = n.log10(f1['halo_properties/mvir'].value)
N_halos = len(all_MS)

# target bits
N_lrg1 = int(volume * dN_dV_lrg1(z_snap)) # 0
N_lrg2 = int(volume * dN_dV_lrg2(z_snap)) # 1
N_elg  = int(volume * dN_dV_elg(z_snap) ) # 2
N_qso  = int(volume * dN_dV_qso(z_snap) ) # 3
N_elg2  = int(volume * dN_dV_elg2(z_snap) ) # 2
#N_Lya  = int(volume * dN_dV_qso(z_snap) ) # 4

lrg1_selection = (n.ones(N_halos)==0)
lrg2_selection = (n.ones(N_halos)==0)
elg_selection =  (n.ones(N_halos)==0)
elg2_selection =  (n.ones(N_halos)==0)
qso_selection =  (n.ones(N_halos)==0)
lya_selection =  (n.ones(N_halos)==0)

print(N_lrg1, N_lrg2, N_elg, N_elg2, N_qso)
if N_lrg1 > 10 or N_lrg2 > 10 :
	all_MS_sort_id = n.argsort(all_MS)
	
if N_lrg1 > 10 :	
	# LRG1: select on Ms
	min_mass = all_MS[all_MS_sort_id[-N_lrg1*3-1]]
	mass_selection = (all_MS>min_mass)
	rds = n.random.rand(len(all_MS))
	lrg1_selection = (mass_selection) & (rds < N_lrg1*1./len(mass_selection.nonzero()[0]))
	#len(lrg1_selection.nonzero()[0]), N_lrg1

if N_lrg2 > 10 :
	# LRG2: select on Ms, less strict. No overlap with LRG1.
	min_mass_2 = all_MS[all_MS_sort_id[-N_lrg2*3-1]]
	mass_selection_2 = (all_MS>min_mass_2) & (lrg1_selection==False)
	rds = n.random.rand(len(all_MS))
	lrg2_selection = (mass_selection_2) & (rds < N_lrg2*1./len(mass_selection_2.nonzero()[0]))
	#print(len(mass_selection_2.nonzero()[0]), len(lrg2_selection.nonzero()[0]), N_lrg2, min_mass_2)

if N_elg > 10 :
	# ELG parameters
	# ELG select on Mvir
	mh_mean, mh_scatter = 12.2, 0.25
	mh_bins = n.arange(mh_mean -2*mh_scatter, mh_mean +2*mh_scatter+0.05, 0.05)
	mh_bins_pos = 0.5*(mh_bins[1:]+mh_bins[:-1])
	proba = lambda x : norm.pdf(x, loc=mh_mean,scale=mh_scatter)
	proba_norm = proba(mh_bins_pos).sum()
	N_2_select_per_bin = (N_elg*proba(mh_bins_pos)/proba_norm).astype('int')
	for id_bin in range(len(mh_bins)-1):
		id_in_bin =(all_mvir > mh_bins[id_bin]) &( all_mvir < mh_bins[id_bin+1]) 
		N_avail = len(id_in_bin.nonzero()[0])
		rds = n.random.rand(len(all_MS))
		bin_selection = (id_in_bin)&(rds < N_2_select_per_bin[id_bin]*1./N_avail)
		elg_selection = (bin_selection)|(elg_selection)

if N_elg2 > 10 :
	# ELG parameters
	# ELG select on Mvir
	mh_mean, mh_scatter = 12.2, 0.25
	mh_bins = n.arange(mh_mean -2*mh_scatter, mh_mean +2*mh_scatter+0.05, 0.05)
	mh_bins_pos = 0.5*(mh_bins[1:]+mh_bins[:-1])
	proba = lambda x : norm.pdf(x, loc=mh_mean,scale=mh_scatter)
	proba_norm = proba(mh_bins_pos).sum()
	N_2_select_per_bin = (N_elg2*proba(mh_bins_pos)/proba_norm).astype('int')
	for id_bin in range(len(mh_bins)-1):
		id_in_bin =(all_mvir > mh_bins[id_bin]) &( all_mvir < mh_bins[id_bin+1]) 
		N_avail = len(id_in_bin.nonzero()[0])
		rds = n.random.rand(len(all_MS))
		bin_selection = (id_in_bin)&(rds < N_2_select_per_bin[id_bin]*1./N_avail)
		elg2_selection = (bin_selection)|(elg2_selection)

if N_qso > 10 :
	# QSO parameters
	# QSO select on Mvir
	mh_mean, mh_scatter = 12.75, 0.25
	mh_bins = n.arange(mh_mean -2*mh_scatter, mh_mean +2*mh_scatter+0.05, 0.05)
	mh_bins_pos = 0.5*(mh_bins[1:]+mh_bins[:-1])
	proba = lambda x : norm.pdf(x, loc=mh_mean,scale=mh_scatter)
	proba_norm = proba(mh_bins_pos).sum()
	N_2_select_per_bin = (N_qso*proba(mh_bins_pos)/proba_norm).astype('int')
	for id_bin in range(len(mh_bins)-1):
		id_in_bin =(all_mvir > mh_bins[id_bin]) &( all_mvir < mh_bins[id_bin+1]) 
		N_avail = len(id_in_bin.nonzero()[0])
		rds = n.random.rand(len(all_MS))
		bin_selection = (id_in_bin)&(rds < N_2_select_per_bin[id_bin]*1./N_avail)
		qso_selection = (bin_selection)|(qso_selection)

ds = f1['/cosmo_4most'].create_dataset('is_ELG_eBOSS', data = elg2_selection )

if status == 'create' :
  halo_data = f1.create_group('cosmo_4most')
  ds = halo_data.create_dataset('is_BG_lz', data = lrg1_selection )
  ds = halo_data.create_dataset('is_BG_hz', data = lrg2_selection )
  ds = halo_data.create_dataset('is_ELG', data = elg_selection )
  ds = halo_data.create_dataset('is_ELG_eBOSS', data = elg2_selection )
  ds = halo_data.create_dataset('is_QSO', data = qso_selection )
  ds = halo_data.create_dataset('is_Lya', data = lya_selection )
  
if status=='update' :
  f1['/cosmo_4most/is_BG_lz'][:]  = lrg1_selection 
  f1['/cosmo_4most/is_BG_hz'][:]  = lrg2_selection 
  f1['/cosmo_4most/is_ELG'][:]  = elg_selection 
  f1['/cosmo_4most/is_ELG_eBOSS'][:]  = elg2_selection 
  f1['/cosmo_4most/is_QSO'][:]  = qso_selection 
  f1['/cosmo_4most/is_Lya'][:]  = lya_selection 
  
f1.close()
