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

t0 = time.time()
out1 =iterate.merging_set_of_system(merger_ids)
t1=time.time()-t0
print( pool.map(f, merger_ids ) )
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
