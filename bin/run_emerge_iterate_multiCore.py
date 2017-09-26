import numpy as n 
import sys
import time

import EmergeIterate

from multiprocessing import Pool
def doIT(n_proc):
	t0=time.time()
	pool = Pool(n_proc)

	iterate = EmergeIterate.EmergeIterate(22, 'MD10')
	iterate.open_snapshots()
	iterate.map_halos_between_snapshots()
	iterate.init_new_quantities()

	# computes the new quantitiess
	if len((iterate.mask_f1_new_halos).nonzero()[0]) > 0 :
		DATA = n.transpose([
			iterate.f1['/halo_properties/mvir'].value[iterate.mask_f1_new_halos], 
			iterate.f1['/halo_properties/rvir'].value[iterate.mask_f1_new_halos], 
			iterate.f1.attrs['redshift']*n.ones_like(iterate.f1['/halo_properties/mvir'].value[iterate.mask_f1_new_halos]), 
			iterate.f1.attrs['age_yr']*n.ones_like(iterate.f1['/halo_properties/mvir'].value[iterate.mask_f1_new_halos]) 
			])

		mvir_dot, rvir_dot, dMdt, dmdt_star, star_formation_rate, stellar_mass = n.array(pool.starmap(EmergeIterate.compute_qtys_new_halos_pk, DATA)).T
		# updates the initiated array with the results
		iterate.mvir_dot[iterate.mask_f1_new_halos]   			=  mvir_dot
		iterate.rvir_dot[iterate.mask_f1_new_halos]   			=  rvir_dot
		iterate.dMdt[iterate.mask_f1_new_halos]   				=  dMdt
		iterate.dmdt_star[iterate.mask_f1_new_halos]  			=  dmdt_star
		iterate.star_formation_rate[iterate.mask_f1_new_halos]	=  star_formation_rate
		iterate.stellar_mass[iterate.mask_f1_new_halos]   		=  stellar_mass

	if len((iterate.mask_f0_evolving_11_halos).nonzero()[0]) > 0 :
		uns = n.ones_like(iterate.f1['/halo_properties/mvir'].value[iterate.mask_f1_evolving_11_halos])
		DATA = n.transpose([
		iterate.f0['/halo_properties/mvir'].value[iterate.mask_f0_evolving_11_halos]
		, iterate.f1['/halo_properties/mvir'].value[iterate.mask_f1_evolving_11_halos]
		, iterate.f0.attrs['age_yr'] * uns
		, iterate.f1.attrs['age_yr'] * uns
		, iterate.f0['/halo_properties/rvir'].value[iterate.mask_f0_evolving_11_halos]
		, iterate.f1['/halo_properties/rvir'].value[iterate.mask_f1_evolving_11_halos]
		, iterate.f1.attrs['redshift'] * uns
		, iterate.t_dynamical[iterate.mask_f1_evolving_11_halos]
		, iterate.f1['/halo_properties/rs'].value[iterate.mask_f1_evolving_11_halos]
		, iterate.f1['/halo_properties/Mpeak'].value[iterate.mask_f1_evolving_11_halos]
		, iterate.f1['/halo_properties/Mpeak_scale'].value[iterate.mask_f1_evolving_11_halos]
		, float(iterate.f1_scale) * uns
		, iterate.f0['/emerge_data/m_icm'].value[iterate.mask_f0_evolving_11_halos]
		, iterate.f0['/emerge_data/stellar_mass'].value[iterate.mask_f0_evolving_11_halos]
		, iterate.f0['/emerge_data/star_formation_rate'].value[iterate.mask_f0_evolving_11_halos]
		])
		
		mvir_dot, rvir_dot, dMdt, dmdt_star, star_formation_rate, stellar_mass, m_icm = n.array( pool.starmap( EmergeIterate.compute_qtys_evolving_halos_pk, DATA)).T
		
		iterate.mvir_dot[iterate.mask_f1_evolving_11_halos]             = mvir_dot
		iterate.rvir_dot[iterate.mask_f1_evolving_11_halos]             = rvir_dot
		iterate.dMdt[iterate.mask_f1_evolving_11_halos]                 = dMdt
		iterate.dmdt_star[iterate.mask_f1_evolving_11_halos]            = dmdt_star
		iterate.star_formation_rate[iterate.mask_f1_evolving_11_halos]  = star_formation_rate
		iterate.stellar_mass[iterate.mask_f1_evolving_11_halos]         = stellar_mass
		iterate.m_icm[iterate.mask_f1_evolving_11_halos]                = m_icm
		

	#if len(iterate.mask_f1_in_a_merging.nonzero()[0]) > 0 :
	#      iterate.compute_qtys_merging_halos()
	print('dt/_proc=',(time.time()-t0)/n_proc)
	
doIT(2)
doIT(4)
doIT(6)
doIT(8)
doIT(10)
doIT(12)

