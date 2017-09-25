#ipython3 

import EmergeIterate
iterate = EmergeIterate.EmergeIterate(22, 'MD10')
iterate.open_snapshots()
iterate.map_halos_between_snapshots()
iterate.init_new_quantities()


merger_id = 23938863
position_f1_host, position_f0_host, position_f0_merging = iterate.get_position_merger_players(merger_id)
iterate.merging_single_system(merger_id)



if len((iterate.mask_f1_new_halos).nonzero()[0]) > 0 :
      iterate.compute_qtys_new_halos()


if len((iterate.mask_f0_evolving_11_halos).nonzero()[0]) > 0 :
      iterate.compute_qtys_evolving_halos() 


if len(iterate.mask_f1_in_a_merging.nonzero()[0]) > 0 :
      iterate.compute_qtys_merging_halos()

#self.write_results()
