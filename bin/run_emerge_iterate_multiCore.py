import time
import EmergeIterate
from multiprocessing import Pool
iterate = EmergeIterate.EmergeIterate(22, 'MD10')
iterate.open_snapshots()
iterate.map_halos_between_snapshots()
iterate.init_new_quantities()

# computes the new quantitiess
pool = Pool(processes=12)
DATA = n.transpose([iterate.f1['/halo_properties/mvir'].value[iterate.mask_f1_new_halos], rvir=iterate.f1['/halo_properties/rvir'].value[iterate.mask_f1_new_halos], iterate.f1.attrs['redshift']*n.ones_like(iterate.f1['/halo_properties/mvir'].value[iterate.mask_f1_new_halos]), iterate.f1.attrs['age_yr']*n.ones_like(iterate.f1['/halo_properties/mvir'].value[iterate.mask_f1_new_halos]) ])
out = p.starmap(iterate.compute_qtys_new_halos_pk, DATA)
mvir_dot, rvir_dot, dMdt, dmdt_star, star_formation_rate, stellar_mass = out

#, f_b=model.f_b, epsilon = model.epsilon(mvir, redshift * n.ones_like(mvir)), f_lost = f_loss(iterate.f1.attrs['age_yr']))
# updates the initiated array with the results
iterate.mvir_dot[iterate.mask_f1_new_halos]   			=  mvir_dot
iterate.rvir_dot[iterate.mask_f1_new_halos]   			=  rvir_dot
iterate.dMdt[iterate.mask_f1_new_halos]   				=  dMdt
iterate.dmdt_star[iterate.mask_f1_new_halos]  			=  dmdt_star
iterate.star_formation_rate[iterate.mask_f1_new_halos]	=  star_formation_rate
iterate.stellar_mass[iterate.mask_f1_new_halos]   		=  stellar_mass


#iterate.compute_qtys_new_halos()
	
if len((iterate.mask_f0_evolving_11_halos).nonzero()[0]) > 0 :
	iterate.compute_qtys_evolving_halos() 
if len(iterate.mask_f1_in_a_merging.nonzero()[0]) > 0 :
	iterate.compute_qtys_merging_halos()
	
# iterate.write_results()


#ipython3 
import time
import EmergeIterate
from multiprocessing import Pool
pool = Pool(processes=12)
f=lambda x:x*x
iterate = EmergeIterate.EmergeIterate(22, 'MD10')
iterate.open_snapshots()
iterate.map_halos_between_snapshots()
iterate.init_new_quantities()

merger_ids = iterate.f1['/halo_properties/id'].value[ iterate.mask_f1_in_a_merging]

prepare_emerge_run

merger_ids

one python job per 1e6 merger_ids

python run_emerge_per_batch index_min index_max

h5/tmp/hlist_xxxx_N_0

t0 = time.time()
out1 =iterate.merging_set_of_system(merger_ids)
t1=time.time()-t0
print( pool.starmap(f, merger_ids.T ) )
#out2 = pool.map(iterate.merging_single_system, merger_ids )
t2=time.time()-t1
print(out1, out2)
print(t2-t1, t1-t0)


"""
merger_id = merger_ids[0] # 23938863

position_f1_host, position_f0_host, position_f0_merging = iterate.get_position_merger_players(merger_id)

iterate.merging_single_system(merger_id)


if len((iterate.mask_f1_new_halos).nonzero()[0]) > 0 :
      iterate.compute_qtys_new_halos()


if len((iterate.mask_f0_evolving_11_halos).nonzero()[0]) > 0 :
      iterate.compute_qtys_evolving_halos() 


if len(iterate.mask_f1_in_a_merging.nonzero()[0]) > 0 :
      iterate.compute_qtys_merging_halos()
#iterate.write_results()

pool = Pool(processes=12)
out = pool.map(iterate.merging_single_system,  ])
"""
