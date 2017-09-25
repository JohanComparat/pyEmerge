# procedure for existing halos, after the first snapshot was initiated :
import sys
# only input parameter is the numbering of the snapshot. These have to be 
#processed in sequence, cannot be done in parallel ... First 1, 2, 3,  ...
#ii = int(sys.argv[1])
#print('snapshot' ii)
import time
t0 = time.time()
from multiprocessing import Pool
#p=Pool(12)

import h5py    
import os
import glob
import numpy as n
import EmergeStellarMass as sm
model = sm.StellarMass()
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115)#, Ob0=0.048206)
import astropy.constants as constants

# generic functions
# =================

f_loss = lambda t : 0.05*n.log( 1 + t / (1.4*10**6))
t_dyn = lambda rvir, mvir : (rvir**3./(9.797465327217671e-24*mvir))**0.5

def tau_quenching( m_star, tdyn, tau_0=4.282, tau_s=0.363):
	out = n.zeros_like(m_star)
	case_1 = (m_star < 1e10 )
	out[case_1] = tdyn[case_1] * tau_0
	case_2 = (case_1==False)
	out[case_2] = tdyn[case_2] * tau_0 * (m_star[case_2] * 10.**(-10.))**(tau_s)
	return out

class EmergeIterate():
	"""
	Loads iterates one step with the Emerge model.
	:param ii: index of the snapshot of interest
	:param env: environment variable of the box. In this dir must be a sub dir 
'h5' with the 'hlist_?.?????_emerge.hdf5' data files in it
	:param L_box: length of the box in Mpc/h 
	
	Running the iteration
	---------------------
	
	ipython3 
	
	import EmergeIterate
	iterate = EmergeIterate.EmergeIterate(12, 'MD10')
	iterate.open_snapshots()
	iterate.map_halos_between_snapshots()
	iterate.init_new_quantities()
	
	if len((iterate.mask_f1_new_halos).nonzero()[0]) > 0 :
		iterate.compute_qtys_new_halos()
	if len((iterate.mask_f0_evolving_11_halos).nonzero()[0]) > 0 :
		iterate.compute_qtys_evolving_halos() 
	if len(self.mask_f1_in_a_merging.nonzero()[0]) > 0 :
		iterate.compute_qtys_merging_halos()
	self.compute_qtys_merging_halos()
	self.write_results()
	
	"""
	def __init__(self, ii, env, L_box=1000.0 ):
		self.ii = ii
		self.env = env
		self.L_box = L_box # box length

	def open_snapshots(self):
		"""
		Opens the files into the class as f0 and f1  
		"""
		h5_dir = os.path.join(os.environ[self.env], 'h5' )
		input_list = n.array(glob.glob(os.path.join(h5_dir, "hlist_?.?????_emerge.hdf5")))
		input_list.sort()
		file_0 = input_list[self.ii-1]
		file_1 = input_list[self.ii]
		
		self.f0 = h5py.File(file_0,  "r")
		self.f0_scale = os.path.basename(file_0).split('_')[1]
		self.positions_f0 = n.arange(len(self.f0['/halo_properties/id'].value))
		
		self.f1 = h5py.File(file_1,  "r+")
		self.f1_scale = os.path.basename(file_1).split('_')[1]
		self.positions_f1 = n.arange(len(self.f1['/halo_properties/id'].value))

	def map_halos_between_snapshots(self):
		"""
		id mapping for halos present in the previous snapshot
		
		Creates 6 arrays to do the mapping in the different cases 
		 * mask_f1_new_halos
		 * mask_f0_evolving_11_halos
		 * mask_f1_evolving_11_halos
		 * f1_id_with_multiple_progenitors
		 * mask_f1_in_a_merging
		 * mask_f0_in_a_merging
		
		"""
		#f0_desc_id_unique_list_all_descendents = n.unique(self.f0['/halo_properties/desc_id'].value)
		f1_id_unique_list_descendents_detected_at_next_scale = n.intersect1d(n.unique(self.f0['/halo_properties/desc_id'].value), self.f1['/halo_properties/id'].value)
		mask_f0_to_propagate = n.in1d(self.f0['/halo_properties/desc_id'].value, f1_id_unique_list_descendents_detected_at_next_scale)
		# mask_f0_lost = (mask_f0_to_propagate == False )
		# evolving halos are given after applying this boolean mask to a f1 quantity :
		mask_f1_evolved_from_previous = n.in1d( self.f1['/halo_properties/id'].value,  f1_id_unique_list_descendents_detected_at_next_scale )
		# new halos are given after applying this boolean mask to a f1 quantity 
		# new halos in f1, not present in f0
		self.mask_f1_new_halos = (mask_f1_evolved_from_previous==False)
		print('new halos', len(self.mask_f1_new_halos.nonzero()[0]))
		# halos descending :
		# mask_f0_to_propagate
		# mask_f1_evolved_from_previous
		s = pd.Series(self.f0['/halo_properties/desc_id'].value[mask_f0_to_propagate])
		self.f1_id_with_multiple_progenitors = s[s.duplicated()].get_values()
		# also = f0_desc_id merging into 1 halo in f1
		# merging systems [many halos in f0 into a single f1 halo]
		self.mask_f1_in_a_merging = n.in1d( self.f1['/halo_properties/id'].value,  self.f1_id_with_multiple_progenitors )
		self.mask_f0_in_a_merging = n.in1d( self.f0['/halo_properties/desc_id'].value,  self.f1_id_with_multiple_progenitors )

		# halos mapped 1::1 between snapshots	
		self.mask_f0_evolving_11_halos = ( mask_f0_to_propagate ) & ( self.mask_f0_in_a_merging == False )
		self.mask_f1_evolving_11_halos = ( mask_f1_evolved_from_previous ) & ( self.mask_f1_in_a_merging == False )
		print('11 mapping', len(self.mask_f0_evolving_11_halos.nonzero()[0]), len(self.mask_f1_evolving_11_halos.nonzero()[0]))
		print('merging systems', len(self.f1_id_with_multiple_progenitors))

		#for dede in self.f1_id_with_multiple_progenitors :
			#sel = self.f0['/halo_properties/desc_id'].value == dede
			#print( dede )
			#print('desc id', self.f0['/halo_properties/desc_id'].value[sel])
			#print('id', self.f0['/halo_properties/id'].value[sel])
			#print('pid', self.f0['/halo_properties/pid'].value[sel])
			#print('mvir', self.f0['/halo_properties/mvir'].value[sel])
			#print('futur merger mmpid', self.f0['/halo_properties/Future_merger_MMP_ID'].value[sel])
			#print('time to future merger', self.f0['/halo_properties/Time_to_future_merger'].value[sel])
			#print('Ms', self.f0['/emerge_data/stellar_mass'].value[sel])
			#print('SFR',self.f0['/emerge_data/star_formation_rate'].value[sel])
			#print('mCIM',self.f0['/emerge_data/m_icm'].value[sel])
			#print('=================================')

	def init_new_quantities(self):
		"""
		Quantities computed for every halos are initialized to 0
		 * mvir_dot            
		 * rvir_dot            
		 * dMdt                 
		 * dmdt_star           
		 * dmdt_star_accretion 
		 * stellar_mass      
		 * star_formation_rate   
		 * m_icm                 
		 * t_dynamical [in years] 
		
		Along with the iteration, these quantities will be updated accordingly
		"""
		self.mvir_dot   =                n.zeros_like(self.f1['/halo_properties/mvir'].value)
		self.rvir_dot   =                n.zeros_like(self.f1['/halo_properties/mvir'].value)
		self.dMdt   =                    n.zeros_like(self.f1['/halo_properties/mvir'].value)
		self.dmdt_star   =               n.zeros_like(self.f1['/halo_properties/mvir'].value)
		self.dmdt_star_accretion   =     n.zeros_like(self.f1['/halo_properties/mvir'].value)
		self.stellar_mass   =            n.zeros_like(self.f1['/halo_properties/mvir'].value)
		self.star_formation_rate   =     n.zeros_like(self.f1['/halo_properties/mvir'].value)
		self.m_icm =                     n.zeros_like(self.f1['/halo_properties/mvir'].value)
		self.t_dynamical = t_dyn( self.f1['/halo_properties/rvir'].value, self.f1['/halo_properties/mvir'].value )
		

	def compute_qtys_new_halos(self):
		"""
		Creates a new galaxy along with the new halo.
		Integrates since the start of the Universe.
		
		Updates the initiated quantities with the values of interest.
		"""
		# evaluate equation (4)
		self.mvir_dot[self.mask_f1_new_halos] = self.f1['/halo_properties/mvir'].value[self.mask_f1_new_halos] / self.f1.attrs['age_yr']
		self.rvir_dot[self.mask_f1_new_halos] = self.f1['/halo_properties/rvir'].value[self.mask_f1_new_halos] / self.f1.attrs['age_yr']
		# no pseudo evolution correction
		self.dMdt[self.mask_f1_new_halos] = self.mvir_dot[self.mask_f1_new_halos] 
		# evaluate equation (1)
		self.dmdt_star[self.mask_f1_new_halos] = model.f_b * self.dMdt[self.mask_f1_new_halos] * model.epsilon(self.f1['/halo_properties/mvir'].value[self.mask_f1_new_halos], self.f1.attrs['redshift'] * n.ones_like(self.f1['/halo_properties/mvir'].value[self.mask_f1_new_halos]))
		# evaluate accretion: 0 in this first step
		# self.dmdt_star_accretion[self.mask_f1_new_halos] = n.zeros_like(self.dmdt_star[self.mask_f1_new_halos])
		# evaluate equation (11)
		f_lost = f_loss(self.f1.attrs['age_yr']) # equation (12)
		# evaluate stellar mass 
		self.star_formation_rate[self.mask_f1_new_halos] = self.dmdt_star[self.mask_f1_new_halos] * (1. - f_lost) + self.dmdt_star_accretion[self.mask_f1_new_halos] 
		self.stellar_mass[self.mask_f1_new_halos] = self.star_formation_rate[self.mask_f1_new_halos] * self.f1.attrs['age_yr'] 
		# intra-cluster mass is currently 0
		# self.m_icm[self.mask_f1_new_halos] = n.zeros_like(self.stellar_mass[self.mask_f1_new_halos])

	def compute_qtys_evolving_halos(self): 
		"""
		update the quantities for evolving halos, present in f0 and f1.
		
		masks :
		 * mask_f1_evolving_11_halos 
		 * mask_f0_evolving_11_halos
		 
		subcases :
		 * quenching : (mvir < Mpeak) & (Mpeak_scale < f1_scale)
		   * case 1. ( age >= t_mpeak ) & ( age_yr < t_mpeak  + t_quench)  
		   * case 2. (age_yr >= t_mpeak  + t_quench) 
		 * stripping, case 1 : (dMdt < 0), then all mass goes to ICM, m=0, mdot=0 
		 * stripping, case 2 : after reaching its peak mass, if M < 0.122 * Mpeak, then all mass goes to ICM, m=0, mdot=0
		"""
		# computing dMdt for the halo
		self.mvir_dot[self.mask_f1_evolving_11_halos] = (self.f1['/halo_properties/mvir'].value[self.mask_f1_evolving_11_halos]-self.f0['/halo_properties/mvir'].value[self.mask_f0_evolving_11_halos]) / (self.f1.attrs['age_yr'] - self.f0.attrs['age_yr'])
		self.rvir_dot[self.mask_f1_evolving_11_halos] = (self.f1['/halo_properties/rvir'].value[self.mask_f1_evolving_11_halos]-self.f0['/halo_properties/rvir'].value[self.mask_f0_evolving_11_halos]) / (self.f1.attrs['age_yr'] - self.f0.attrs['age_yr'])
		c = self.f1['/halo_properties/rvir'].value[self.mask_f1_evolving_11_halos] / self.f1['/halo_properties/rs'].value[self.mask_f1_evolving_11_halos]
		rho_nfw = self.f1['/halo_properties/mvir'].value[self.mask_f1_evolving_11_halos] / (self.f1['/halo_properties/rs'].value[self.mask_f1_evolving_11_halos]**3. * 4. * n.pi * c * (1+c)**2. * (n.log(1.+c)-c/(1.+c)))
		pseudo_evolution_correction = 4.*n.pi*self.f1['/halo_properties/rvir'].value[self.mask_f1_evolving_11_halos] *self.f1['/halo_properties/rvir'].value[self.mask_f1_evolving_11_halos] * self.rvir_dot[self.mask_f1_evolving_11_halos] * rho_nfw
		self.dMdt[self.mask_f1_evolving_11_halos] = self.mvir_dot[self.mask_f1_evolving_11_halos] - pseudo_evolution_correction
		
		# initialize the ICM mass to the previous value
		self.m_icm[self.mask_f1_evolving_11_halos] = self.f0['/emerge_data/m_icm'].value[self.mask_f0_evolving_11_halos]
		
		# Direct estimates of stellar mass and SFR
		self.dmdt_star[self.mask_f1_evolving_11_halos] = model.f_b * self.dMdt[self.mask_f1_evolving_11_halos] * model.epsilon(self.f1['/halo_properties/mvir'].value[self.mask_f1_evolving_11_halos], self.f1.attrs['redshift'] * n.ones_like(self.f1['/halo_properties/mvir'].value[self.mask_f1_evolving_11_halos]))
		# evaluate accretion: 0 in this first step
		# dmdt_star_accretion = n.zeros_like(dmdt_star[self.mask_f1_evolving_11_halos])
		# evaluate equation (11)
		f_lost = f_loss(self.f1.attrs['age_yr']-self.f0.attrs['age_yr'])
		# evaluate stellar mass 
		self.star_formation_rate[self.mask_f1_evolving_11_halos] = self.dmdt_star[self.mask_f1_evolving_11_halos] * (1. - f_lost) + self.dmdt_star_accretion[self.mask_f1_evolving_11_halos] 
		self.stellar_mass[self.mask_f1_evolving_11_halos] = self.star_formation_rate[self.mask_f1_evolving_11_halos] * (self.f1.attrs['age_yr']-self.f0.attrs['age_yr']) + self.f0['/emerge_data/stellar_mass'].value[self.mask_f0_evolving_11_halos]
		# Variations due to stripping, merging and quenching
	
		# quenching
		quenching = (self.f1['/halo_properties/mvir'].value[self.mask_f1_evolving_11_halos] < self.f1['/halo_properties/Mpeak'].value[self.mask_f1_evolving_11_halos]) & (self.f1['/halo_properties/Mpeak_scale'].value[self.mask_f1_evolving_11_halos] < float(self.f1_scale))
		t_quench = tau_quenching( self.f0['/emerge_data/stellar_mass'].value[self.mask_f0_evolving_11_halos], self.t_dynamical[self.mask_f1_evolving_11_halos] )
		t_mpeak = cosmoMD.age( 1. / self.f1['/halo_properties/Mpeak_scale'].value[self.mask_f1_evolving_11_halos] - 1. ).to(u.yr).value
		
		# case 1. mdot = mdot at tpeak
		quench_1 = (quenching) & (self.f1.attrs['age_yr'] >= t_mpeak ) & ( self.f1.attrs['age_yr'] < t_mpeak + t_quench) 
		if len(quench_1.nonzero()[0])>0:
			print("quenching1")
			self.star_formation_rate[self.mask_f1_evolving_11_halos][quench_1] = n.ones_like(self.star_formation_rate[self.mask_f1_evolving_11_halos][quench_1])*self.f0['/emerge_data/stellar_mass'].value[self.mask_f0_evolving_11_halos][quench_1]
			self.stellar_mass[self.mask_f1_evolving_11_halos][quench_1] = self.star_formation_rate[self.mask_f1_evolving_11_halos][quench_1] * (self.f1.attrs['age_yr']-self.f0.attrs['age_yr']) + self.f0['/emerge_data/stellar_mass'].value[self.mask_f0_evolving_11_halos][quench_1]

		# case 2. m dot =0
		quench_2 = (quenching) &(self.f1.attrs['age_yr'] >= t_mpeak + t_quench )
		if len(quench_2.nonzero()[0])>0:
			print("quenching2")
			self.star_formation_rate[self.mask_f1_evolving_11_halos][quench_2] = n.zeros_like(self.star_formation_rate[self.mask_f1_evolving_11_halos][quench_2])
			self.stellar_mass[self.mask_f1_evolving_11_halos][quench_2] = self.f0['/emerge_data/stellar_mass'].value[self.mask_f0_evolving_11_halos][quench_2]
		
		# stripping, case 1
		# negative growth value self.dMdt[self.mask_f1_evolving_11_halos] => 0
		stripping_1 = (self.dMdt[self.mask_f1_evolving_11_halos] < 0) 
		# stripping, case 2
		# after reaching its peak mass,
		# if M < 0.122 * Mpeak, all mass goes to ICM, m=0, mdot=0
		stripping_2 = (self.f1['/halo_properties/mvir'].value[self.mask_f1_evolving_11_halos] < 0.122*self.f1['/halo_properties/Mpeak'].value[self.mask_f1_evolving_11_halos]) & (self.f1['/halo_properties/Mpeak_scale'].value[self.mask_f1_evolving_11_halos] < float(self.f1_scale))
		# both cases together
		stripping = (stripping_1) | (stripping_1)
		if len(stripping.nonzero()[0])>0:
			print("stripping")
			self.m_icm[self.mask_f1_evolving_11_halos][stripping] += self.f0['/emerge_data/stellar_mass'].value[self.mask_f0_evolving_11_halos][stripping]
			self.stellar_mass[self.mask_f1_evolving_11_halos][stripping] = n.zeros_like(self.stellar_mass[self.mask_f1_evolving_11_halos][stripping])
			self.star_formation_rate[self.mask_f1_evolving_11_halos][stripping] = n.zeros_like(self.star_formation_rate[self.mask_f1_evolving_11_halos][stripping])
		

	def get_position_merger_players(self, merger_id):
		"""
		Given the identifier of the merger
		:param merger_id: id of the parent halo of the merger at the later time. One integer.
		
		Outputs the position on the f0 and f1 arrays of the hosts and of the merging systems
		
		returns :
		 position_f1_host [int], position_f0_host [int], position_f0_merging [list]
		"""
		# about the host at t1
		#print(merger_id)
		mask_f1_host = (self.f1['/halo_properties/id'].value == merger_id)
		#print(mask_f1_host)
		position_f1_host = self.positions_f1[mask_f1_host]
		#print(position_f1_host)
		# about the host and merging subhalos at t0
		mask_f0_all = (self.f0['/halo_properties/desc_id'].value == merger_id)
		#print(mask_f0_all)
		id_f0_all = self.f0['/halo_properties/id'].value[mask_f0_all]
		#print(id_f0_all)
		# the host at t1 is flagged at t0 as the most massive progenitor
		#print(n.unique(self.f0['/halo_properties/Future_merger_MMP_ID'].value[mask_f0_all]))
		#print(n.in1d(id_f0_all, n.unique(self.f0['/halo_properties/Future_merger_MMP_ID'].value[mask_f0_all])))
		#print(id_f0_all[n.in1d(id_f0_all, n.unique(self.f0['/halo_properties/Future_merger_MMP_ID'].value[mask_f0_all]))])
		f0_host_id = id_f0_all[n.in1d(id_f0_all, n.unique(self.f0['/halo_properties/Future_merger_MMP_ID'].value[mask_f0_all]))][0] 
		mask_f0_host = (mask_f0_all) & (self.f0['/halo_properties/id'].value == f0_host_id)
		mask_f0_merging = (mask_f0_all) & (self.f0['/halo_properties/id'].value != f0_host_id)
		position_f0_host    = self.positions_f0[mask_f0_host]
		position_f0_merging = self.positions_f0[mask_f0_merging]
		return position_f1_host, position_f0_host, position_f0_merging

	def merging_single_system(self, merger_id):
		"""
		:param merger_id: id of the parent halo of the merger at the later time. One integer.		
		
		Merging goes as follows. Assume escape fraction: f_esc = 0.388, then 
		 * m_star_satellite x f_esc     goes to m_host_ICM
		 * m_star_satellite x (1-f_esc) goes to m_star_host
		
		returns :
		 parameters of the emerge model of the galaxies undergoing merger at this point.
		 [ mvir_dot, rvir_dot, dMdt, dmdt_star, dmdt_star_accretion, stellar_mass, star_formation_rate, m_icm ]
		
		"""
		position_f1_host, position_f0_host, position_f0_merging = self.get_position_merger_players(merger_id)
		mvir_dot = (self.f1['/halo_properties/mvir'].value[position_f1_host]-self.f0['/halo_properties/mvir'].value[position_f0_host]) / (self.f1.attrs['age_yr'] - self.f0.attrs['age_yr'])
		rvir_dot = (self.f1['/halo_properties/rvir'].value[position_f1_host]-self.f0['/halo_properties/rvir'].value[position_f0_host]) / (self.f1.attrs['age_yr'] - self.f0.attrs['age_yr'])
		c = self.f1['/halo_properties/rvir'].value[position_f1_host] / self.f1['/halo_properties/rs'].value[position_f1_host]
		rho_nfw = self.f1['/halo_properties/mvir'].value[position_f1_host] / (self.f1['/halo_properties/rs'].value[position_f1_host]**3. * 4. * n.pi * c * (1+c)**2. * (n.log(1.+c)-c/(1.+c)))
		pseudo_evolution_correction = 4.*n.pi*self.f1['/halo_properties/rvir'].value[position_f1_host] *self.f1['/halo_properties/rvir'].value[position_f1_host] * rvir_dot * rho_nfw
		dMdt = mvir_dot - pseudo_evolution_correction
		# initialize the ICM mass to the previous value
		m_icm = self.f0['/emerge_data/m_icm'].value[position_f0_host]
		# Direct estimates of stellar mass and SFR
		dmdt_star = model.f_b * dMdt * model.epsilon(self.f1['/halo_properties/mvir'].value[position_f1_host], self.f1.attrs['redshift'] * n.ones_like(self.f1['/halo_properties/mvir'].value[position_f1_host]))
		# evaluate accretion: 0 in this first step
		dmdt_star_accretion = n.zeros_like(dmdt_star)
		# evaluate equation (11)
		f_lost = f_loss(self.f1.attrs['age_yr']-self.f0.attrs['age_yr'])
		# evaluate stellar mass 
		star_formation_rate = dmdt_star * (1. - f_lost) + dmdt_star_accretion 
		stellar_mass = star_formation_rate * (self.f1.attrs['age_yr']-self.f0.attrs['age_yr']) + self.f0['/emerge_data/stellar_mass'].value[position_f0_host]
		# merging
		# m_star_sat x f_esc => m_host_ICM
		# m_star_sat x (1-f_esc) => m_star_host
		# f_esc = 0.388
		#Time_to_future_merger: Time (in Gyr) until the given halo merges into a larger halo.  (-1 if no future merger happens)
		#Future_merger_MMP_ID: most-massive progenitor of the halo into which the given halo merges. (-1 if the main progenitor of the future merger halo does not exist at the given scale factor.)

		stellar_mass += (1.-0.388)*n.sum(self.f0['/emerge_data/stellar_mass'].value[position_f0_merging])
		m_icm += 0.388*n.sum(self.f0['/emerge_data/stellar_mass'].value[position_f0_merging])

		return mvir_dot, rvir_dot, dMdt, dmdt_star, dmdt_star_accretion, stellar_mass, star_formation_rate, m_icm
	
	def merging_set_of_system(self, merger_ids):
		"""
		Loops over self.merging_single_system over a list of ids and returns a merged output array
		"""
		return n.hstack(( n.array([self.merging_single_system(merger_id) for merger_id in merger_ids]) ))
	
	def compute_qtys_merging_halos(self): 
		"""
		computes all quantities for merging halos
		"""
		pool = Pool(processes=12)
		self.out3 = pool.map(self.merging_set_of_system, self.f1['/halo_properties/id'].value[ self.mask_f1_in_a_merging ])
		#self.out3 = p.starmap(self.merging_set_of_system, self.f1['/halo_properties/id'].value[ self.mask_f1_in_a_merging ])
		
	def write_results(self):
		"""
		After computing all quantities, you need to write the results in the h5 file.
		"""
		emerge_data = self.f1.create_group('emerge_data')

		#emerge_data.attrs['f_lost'] = f_lost

		ds = emerge_data.create_dataset('mvir_dot', data = self.mvir_dot )
		ds.attrs['units'] = r'$h^{-1} M_\odot / yr$'
		ds.attrs['long_name'] = r'$d M_{vir} / dt$' 

		ds = emerge_data.create_dataset('rvir_dot', data = self.rvir_dot )
		ds.attrs['units'] = r'$h^{-1} kpc / yr$'
		ds.attrs['long_name'] = r'$d r_{vir} / dt$' 

		ds = emerge_data.create_dataset('dMdt', data = self.dMdt )
		ds.attrs['units'] = r'$h^{-1} M_\odot / yr$'
		ds.attrs['long_name'] = r'$\langle d M / dt \rangle$ (4)' 

		ds = emerge_data.create_dataset('dmdt_star', data = self.dmdt_star )
		ds.attrs['units'] = r'$h^{-1} M_\odot / yr$'
		ds.attrs['long_name'] = r'$ d m_* / dt $ (1)' 

		ds = emerge_data.create_dataset('dmdt_star_accretion', data = 
		self.dmdt_star_accretion )
		ds.attrs['units'] = r'$h^{-1} M_\odot / yr$'
		ds.attrs['long_name'] = r'$ d m_{acc} / dt $ ' 

		ds = emerge_data.create_dataset('star_formation_rate', data = 
		self.star_formation_rate )
		ds.attrs['units'] = r'$h^{-1} M_\odot / yr$'
		ds.attrs['long_name'] = r'$ d m / dt $ ' 

		ds = emerge_data.create_dataset('stellar_mass', data = self.stellar_mass )
		ds.attrs['units'] = r'$h^{-1} M_\odot $'
		ds.attrs['long_name'] = r'$ m_* $ (11)' 

		ds = emerge_data.create_dataset('m_icm', data = self.m_icm )
		ds.attrs['units'] = r'$h^{-1} M_\odot $'
		ds.attrs['long_name'] = r'$ m_{ICM}$ ' 

		self.f0.close()
		self.f1.close()

		print("Results written")



if __name__ == '__main__':
	iterate = EmergeIterate.EmergeIterate(22, 'MD10')
	iterate.open_snapshots()
	iterate.map_halos_between_snapshots()
	iterate.init_new_quantities()
	
	if len((iterate.mask_f1_new_halos).nonzero()[0]) > 0 :
		iterate.compute_qtys_new_halos()
	if len((iterate.mask_f0_evolving_11_halos).nonzero()[0]) > 0 :
		iterate.compute_qtys_evolving_halos() 
	if len(iterate.mask_f1_in_a_merging.nonzero()[0]) > 0 :
		iterate.compute_qtys_merging_halos()
		
	# iterate.write_results()
