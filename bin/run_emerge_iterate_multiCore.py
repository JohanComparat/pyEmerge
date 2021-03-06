import time
t0 = time.time()
import sys

# environment variable to the simulation directory
sim_dir = sys.argv[1]
# snapshot id (within a sorted list of ids)
snap_id = int(sys.argv[2])
# python3 run_emerge_iterate_multiCore.py MD10 22
hh = 0.6777

import numpy as n 
import pandas as pd

import EmergeIterate

from multiprocessing import Pool
t0=time.time()
n_proc=12
pool = Pool(n_proc)

# start the emerge method 
iterate = EmergeIterate.EmergeIterate(snap_id, sim_dir)
iterate.open_snapshots()
iterate.init_new_quantities()

# options 
run_new_halos = True
run_evolving_halos = True
run_merging_halos = True

#################################################3
#################################################3
#################################################3
#                    new halos
#################################################3
#################################################3
#################################################3
t1=time.time()
print('elapsed time', t1-t0, 'seconds')
print('---------------------------------')
print(iterate.f1.attrs['aexp'])
if run_new_halos:
	print('N f1 halos', len(iterate.f1['/halo_properties/id'].value))
	print('N f0 halos', len(iterate.f0['/halo_properties/id'].value))
	f1_new_halos = (n.in1d(iterate.f1['/halo_properties/id'].value, iterate.f0['/halo_properties/desc_id'].value)==False)
	# new halos [f1_new_halos] to be fed to EmergeIterate.compute_qtys_new_halos_pk
	print('f1 new halos', len(iterate.f1['/halo_properties/desc_id'].value[f1_new_halos]))
	print('computing galaxies in new halos')
	if len((f1_new_halos).nonzero()[0]) > 0 :
		DATA = n.transpose([
			iterate.f1['/halo_properties/mvir'].value[f1_new_halos]*hh**(-1), 
			iterate.f1['/halo_properties/rvir'].value[f1_new_halos]*hh**(-1), 
			iterate.f1.attrs['redshift']*n.ones_like(iterate.f1['/halo_properties/mvir'].value[f1_new_halos]), 
			iterate.f1.attrs['age_yr']*n.ones_like(iterate.f1['/halo_properties/mvir'].value[f1_new_halos]) 
			])

		mvir_dot, rvir_dot, dMdt, dmdt_star, star_formation_rate, stellar_mass = n.array(pool.starmap(EmergeIterate.compute_qtys_new_halos_pk, DATA)).T
		# updates the initiated array with the results
		#iterate.mvir_dot[f1_new_halos]   			=  mvir_dot
		#iterate.rvir_dot[f1_new_halos]   			=  rvir_dot
		#iterate.dMdt[f1_new_halos]   				=  dMdt
		#iterate.dmdt_star[f1_new_halos]  			=  dmdt_star
		#iterate.star_formation_rate[f1_new_halos]	=  star_formation_rate
		#iterate.stellar_mass[f1_new_halos]   		=  stellar_mass

		iterate.f1['/emerge_data/mvir_dot'][f1_new_halos] = mvir_dot 
		iterate.f1['/emerge_data/rvir_dot'][f1_new_halos] = rvir_dot 
		iterate.f1['/emerge_data/dMdt'][f1_new_halos] = dMdt 
		iterate.f1['/emerge_data/dmdt_star'][f1_new_halos] = dmdt_star 
		#iterate.f1['/emerge_data/dmdt_star_accretion'][f1_new_halos] = dmdt_star_accretion 
		iterate.f1['/emerge_data/star_formation_rate'][f1_new_halos] = star_formation_rate 
		iterate.f1['/emerge_data/stellar_mass'][f1_new_halos] = stellar_mass 
		#iterate.f1['/emerge_data/m_icm'][f1_new_halos] = m_icm 
		print("Results updated")


t2=time.time()
print('elapsed time', t2-t1, 'seconds to create new galaxies')
print('---------------------------------')
#################################################3
#################################################3
#################################################3
#              1:1 evolving halos
#################################################3
#################################################3
#################################################3
if run_evolving_halos:
	print('evolving halos 1:1')
	f1_evolved_halos = (f1_new_halos==False)

	f0_propagated_halos = n.in1d(iterate.f0['/halo_properties/desc_id'].value, iterate.f1['/halo_properties/id'].value)
	f0_lost_halos = (f0_propagated_halos==False)

	# parmi f0_propagated_halos
	# cherchez les desc_id qui apparaissent plusieurs fois
	s = pd.Series(iterate.f0['/halo_properties/desc_id'].value[f0_propagated_halos])
	duplicated_desc_ids = s[s.duplicated()].get_values()
	f0_with_desc_id_duplicates = n.in1d( iterate.f0['/halo_properties/desc_id'].value,  duplicated_desc_ids )
	f0_in_a_merging = f0_with_desc_id_duplicates & f0_propagated_halos
	f0_not_merging = (f0_with_desc_id_duplicates == False) & f0_propagated_halos

	print('f0 lost halos', len(iterate.f0['/halo_properties/desc_id'].value[f0_lost_halos]))
	print('f0 propagated halos', len(iterate.f0['/halo_properties/desc_id'].value[f0_propagated_halos]))
	print('f0 merging halos', len(iterate.f0['/halo_properties/id'].value[f0_in_a_merging]))
	print('f0 not merging halos', len(iterate.f0['/halo_properties/id'].value[f0_not_merging]))

	print('f1 halos with progenitor', len(iterate.f1['/halo_properties/id'].value[f1_evolved_halos]))


	# mapper les 12964 sur les 12981 avec selection tel que
	# iterate.f1['/halo_properties/id'].value[f1_evolved_halos & selection] = iterate.f0['/halo_properties/desc_id'].value[f0_not_merging]
	# comme les indices sont ordonnes, il suffit de retirer ceux qui apparaissent plusieurs fois
	# selection = mask_f1_in_a_merging == False
	mask_f1_in_a_merging = n.in1d( iterate.f1['/halo_properties/id'].value,  duplicated_desc_ids )
	f1_evolved_halos_no_merger = (f1_evolved_halos) & (mask_f1_in_a_merging==False)
	f1_evolved_halos_with_merger = (f1_evolved_halos) & (mask_f1_in_a_merging)

	# evolving halos to be fed to EmergeIterate.compute_qtys_evolving_halos_pk
	# 1:1 with f0_not_merging <=> f1_evolved_halos_no_merger
	print("f0 and f1 connected halos without merger",len(f0_not_merging.nonzero()[0]), len(f1_evolved_halos_no_merger.nonzero()[0]))


	if len((f0_not_merging).nonzero()[0]) > 0 :
		uns = n.ones_like(iterate.f1['/halo_properties/mvir'].value[f1_evolved_halos_no_merger])
		DATA = n.transpose([
		iterate.f0['/halo_properties/mvir'].value[f0_not_merging]*hh**(-1)
		, iterate.f1['/halo_properties/mvir'].value[f1_evolved_halos_no_merger]*hh**(-1)
		, iterate.f0.attrs['age_yr'] * uns
		, iterate.f1.attrs['age_yr'] * uns
		, iterate.f0['/halo_properties/rvir'].value[f0_not_merging]*hh**(-1)
		, iterate.f1['/halo_properties/rvir'].value[f1_evolved_halos_no_merger]*hh**(-1)
		, iterate.f1.attrs['redshift'] * uns
		, iterate.t_dynamical[f1_evolved_halos_no_merger]
		, iterate.f1['/halo_properties/rs'].value[f1_evolved_halos_no_merger]*hh**(-1)
		, iterate.f1['/halo_properties/Mpeak'].value[f1_evolved_halos_no_merger]*hh**(-1)
		, iterate.f1['/halo_properties/Mpeak_scale'].value[f1_evolved_halos_no_merger]
		, float(iterate.f1_scale) * uns
		, iterate.f0['/emerge_data/m_icm'].value[f0_not_merging]*hh**(-1)
		, iterate.f0['/emerge_data/stellar_mass'].value[f0_not_merging]
		, iterate.f0['/emerge_data/star_formation_rate'].value[f0_not_merging]
		])
		
		mvir_dot, rvir_dot, dMdt, dmdt_star, star_formation_rate, stellar_mass, m_icm = n.array( pool.starmap( EmergeIterate.compute_qtys_evolving_halos_pk, DATA)).T
		
		#iterate.mvir_dot[f1_evolved_halos_no_merger]             = mvir_dot
		#iterate.rvir_dot[f1_evolved_halos_no_merger]             = rvir_dot
		#iterate.dMdt[f1_evolved_halos_no_merger]                 = dMdt
		#iterate.dmdt_star[f1_evolved_halos_no_merger]            = dmdt_star
		#iterate.star_formation_rate[f1_evolved_halos_no_merger]  = star_formation_rate
		#iterate.stellar_mass[f1_evolved_halos_no_merger]         = stellar_mass
		#iterate.m_icm[f1_evolved_halos_no_merger]                = m_icm
		
		iterate.f1['/emerge_data/mvir_dot'][f1_evolved_halos_no_merger] = mvir_dot 
		iterate.f1['/emerge_data/rvir_dot'][f1_evolved_halos_no_merger] = rvir_dot 
		iterate.f1['/emerge_data/dMdt'][f1_evolved_halos_no_merger] = dMdt 
		iterate.f1['/emerge_data/dmdt_star'][f1_evolved_halos_no_merger] = dmdt_star 
		#iterate.f1['/emerge_data/dmdt_star_accretion'][f1_evolved_halos_no_merger] = dmdt_star_accretion 
		iterate.f1['/emerge_data/star_formation_rate'][f1_evolved_halos_no_merger] = star_formation_rate 
		iterate.f1['/emerge_data/stellar_mass'][f1_evolved_halos_no_merger] = stellar_mass 
		iterate.f1['/emerge_data/m_icm'][f1_evolved_halos_no_merger] = m_icm 
		print("Results updated")


t3=time.time()
print('elapsed time', t3-t2,'seconds to evolve 1:1 galaxies')
print('---------------------------------')


#################################################3
#################################################3
#################################################3
#              merging halos
#################################################3
#################################################3
#################################################3
# merging halos to be fed to EmergeIterate.get_position_merger_players_pk
# X=>1 with f0_in_a_merging => f1_evolved_halos_with_merger
#print("merging halos in f0 and f1",len(f0_in_a_merging.nonzero()[0]), len(f1_evolved_halos_with_merger.nonzero()[0])))
# 34 => 17
if run_merging_halos:
	merger_ids = iterate.f1['/halo_properties/id'].value[f1_evolved_halos_with_merger]

	if len(merger_ids) > 0 :
		positions_f0 = n.arange(len(iterate.f0['/halo_properties/desc_id'].value[f0_in_a_merging]))
		sum_stellar_mass_guests=[] #n.zeros_like(merger_ids)
		positions_hosts=[]

		for merger_id in merger_ids:
			f0_corr_mid = (iterate.f0['/halo_properties/desc_id'].value[f0_in_a_merging] == merger_id)
			id_host_f0 = n.argmax(iterate.f0['/halo_properties/mvir'].value[f0_in_a_merging][f0_corr_mid])
			positions_hosts.append(positions_f0[f0_corr_mid][id_host_f0])
			id_guests_f0 = n.delete(n.arange(len(iterate.f0['/halo_properties/mvir'].value[f0_in_a_merging][f0_corr_mid])), id_host_f0)
			sum_stellar_mass_guests.append( n.sum(iterate.f0['/emerge_data/stellar_mass'].value[f0_in_a_merging][f0_corr_mid][id_guests_f0]))

		hosts_f0 = n.array(positions_hosts)
		sum_stellar_mass_guests = n.array(sum_stellar_mass_guests)

		uns = n.ones_like(iterate.f1['/halo_properties/mvir'].value[f1_evolved_halos_with_merger])
		DATA = n.transpose([
		iterate.f0['/halo_properties/mvir'].value[f0_in_a_merging][hosts_f0]*hh**(-1)
		, iterate.f1['/halo_properties/mvir'].value[f1_evolved_halos_with_merger]*hh**(-1)
		, iterate.f0.attrs['age_yr'] * uns
		, iterate.f1.attrs['age_yr'] * uns
		, iterate.f0['/halo_properties/rvir'].value[f0_in_a_merging][hosts_f0]*hh**(-1)
		, iterate.f1['/halo_properties/rvir'].value[f1_evolved_halos_with_merger]*hh**(-1)
		, iterate.f1.attrs['redshift'] * uns
		, iterate.t_dynamical[f1_evolved_halos_with_merger]
		, iterate.f1['/halo_properties/rs'].value[f1_evolved_halos_with_merger]*hh**(-1)
		, iterate.f1['/halo_properties/Mpeak'].value[f1_evolved_halos_with_merger]*hh**(-1)
		, iterate.f1['/halo_properties/Mpeak_scale'].value[f1_evolved_halos_with_merger]
		, float(iterate.f1_scale) * uns
		, iterate.f0['/emerge_data/m_icm'].value[f0_in_a_merging][hosts_f0]*hh**(-1)
		, iterate.f0['/emerge_data/stellar_mass'].value[f0_in_a_merging][hosts_f0]
		, iterate.f0['/emerge_data/star_formation_rate'].value[f0_in_a_merging][hosts_f0]
		, sum_stellar_mass_guests
		])
		
		mvir_dot, rvir_dot, dMdt, dmdt_star, star_formation_rate, stellar_mass, m_icm = n.array( pool.starmap( EmergeIterate.merge_system, DATA)).T
		
		iterate.mvir_dot[f1_evolved_halos_with_merger]             = mvir_dot
		iterate.rvir_dot[f1_evolved_halos_with_merger]             = rvir_dot
		iterate.dMdt[f1_evolved_halos_with_merger]                 = dMdt
		iterate.dmdt_star[f1_evolved_halos_with_merger]            = dmdt_star
		iterate.star_formation_rate[f1_evolved_halos_with_merger]  = star_formation_rate
		iterate.stellar_mass[f1_evolved_halos_with_merger]         = stellar_mass
		iterate.m_icm[f1_evolved_halos_with_merger]                = m_icm
		
		iterate.f1['/emerge_data/mvir_dot'][f1_evolved_halos_with_merger] = mvir_dot 
		iterate.f1['/emerge_data/rvir_dot'][f1_evolved_halos_with_merger] = rvir_dot 
		iterate.f1['/emerge_data/dMdt'][f1_evolved_halos_with_merger] = dMdt 
		iterate.f1['/emerge_data/dmdt_star'][f1_evolved_halos_with_merger] = dmdt_star 
		#iterate.f1['/emerge_data/dmdt_star_accretion'][f1_evolved_halos_with_merger] = dmdt_star_accretion 
		iterate.f1['/emerge_data/star_formation_rate'][f1_evolved_halos_with_merger] = star_formation_rate 
		iterate.f1['/emerge_data/stellar_mass'][f1_evolved_halos_with_merger] = stellar_mass 
		iterate.f1['/emerge_data/m_icm'][f1_evolved_halos_with_merger] = m_icm 
		print("Results updated")


# setting to 0 the negative SFRs and negative masses
sfr_neg = (iterate.f1['/emerge_data/star_formation_rate'].value < 0)
mass_neg = (iterate.f1['/emerge_data/stellar_mass'].value < 0)
print("neg sfr", len(sfr_neg.nonzero()[0]), 'neg masses', len(mass_neg.nonzero()[0]) )
iterate.f1['/emerge_data/star_formation_rate'][sfr_neg] = n.zeros_like(iterate.f1['/emerge_data/star_formation_rate'][sfr_neg]) 
iterate.f1['/emerge_data/stellar_mass'][mass_neg] = n.zeros_like(iterate.f1['/emerge_data/stellar_mass'][mass_neg]) 

t4=time.time()
print('elapsed time', t4-t3,'seconds to evolve 1:1 halos')
print('---------------------------------')
#t5=time.time()
#print('elapsed time', t5-t4, 'seconds getting merging masses')
#print('---------------------------------')
#t6=time.time()
#print('elapsed time merger', t6-t5, 'seconds to do the mergers')
#print('---------------------------------')
print('total time ',time.time()-t0,'seconds, total time/ n proc=',(time.time()-t0)/n_proc)

iterate.f0.close()
iterate.f1.close()

