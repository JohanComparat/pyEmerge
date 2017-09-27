import glob
import os
import numpy as n

h5_dir = os.path.join(os.environ['MD10'], 'h5' )
input_list = n.array(glob.glob(os.path.join(h5_dir, "hlist_?.?????_emerge.hdf5")))
input_list.sort()

input_list[30:]
hlist_0.28290_emerge.hdf5

MD10 63  -- 113
#!/bin/bash
#python3 run_emerge_iterate_multiCore.py MD10 60
#python3 run_emerge_iterate_multiCore.py MD10 61
#python3 run_emerge_iterate_multiCore.py MD10 62
python3 run_emerge_iterate_multiCore.py MD10 63
python3 run_emerge_iterate_multiCore.py MD10 64
python3 run_emerge_iterate_multiCore.py MD10 65
python3 run_emerge_iterate_multiCore.py MD10 66
python3 run_emerge_iterate_multiCore.py MD10 67
python3 run_emerge_iterate_multiCore.py MD10 68
python3 run_emerge_iterate_multiCore.py MD10 69
python3 run_emerge_iterate_multiCore.py MD10 70
python3 run_emerge_iterate_multiCore.py MD10 71
python3 run_emerge_iterate_multiCore.py MD10 72
python3 run_emerge_iterate_multiCore.py MD10 73
python3 run_emerge_iterate_multiCore.py MD10 74
python3 run_emerge_iterate_multiCore.py MD10 75
python3 run_emerge_iterate_multiCore.py MD10 76
python3 run_emerge_iterate_multiCore.py MD10 77
python3 run_emerge_iterate_multiCore.py MD10 78
python3 run_emerge_iterate_multiCore.py MD10 79
python3 run_emerge_iterate_multiCore.py MD10 80
python3 run_emerge_iterate_multiCore.py MD10 81
python3 run_emerge_iterate_multiCore.py MD10 82
python3 run_emerge_iterate_multiCore.py MD10 83
python3 run_emerge_iterate_multiCore.py MD10 84
python3 run_emerge_iterate_multiCore.py MD10 85
python3 run_emerge_iterate_multiCore.py MD10 86
python3 run_emerge_iterate_multiCore.py MD10 87
python3 run_emerge_iterate_multiCore.py MD10 88
python3 run_emerge_iterate_multiCore.py MD10 89
python3 run_emerge_iterate_multiCore.py MD10 90
python3 run_emerge_iterate_multiCore.py MD10 91
python3 run_emerge_iterate_multiCore.py MD10 92
python3 run_emerge_iterate_multiCore.py MD10 93
python3 run_emerge_iterate_multiCore.py MD10 94
python3 run_emerge_iterate_multiCore.py MD10 95
python3 run_emerge_iterate_multiCore.py MD10 96
python3 run_emerge_iterate_multiCore.py MD10 97
python3 run_emerge_iterate_multiCore.py MD10 98
python3 run_emerge_iterate_multiCore.py MD10 99
python3 run_emerge_iterate_multiCore.py MD10 100
python3 run_emerge_iterate_multiCore.py MD10 101
python3 run_emerge_iterate_multiCore.py MD10 102
python3 run_emerge_iterate_multiCore.py MD10 103
python3 run_emerge_iterate_multiCore.py MD10 104
python3 run_emerge_iterate_multiCore.py MD10 105
python3 run_emerge_iterate_multiCore.py MD10 106
python3 run_emerge_iterate_multiCore.py MD10 107
python3 run_emerge_iterate_multiCore.py MD10 108
python3 run_emerge_iterate_multiCore.py MD10 109
python3 run_emerge_iterate_multiCore.py MD10 110
python3 run_emerge_iterate_multiCore.py MD10 111
python3 run_emerge_iterate_multiCore.py MD10 112
python3 run_emerge_iterate_multiCore.py MD10 113

h5_dir = os.path.join(os.environ['MD04'], 'h5' )
input_list = n.array(glob.glob(os.path.join(h5_dir, "hlist_?.?????_emerge.hdf5")))
input_list.sort()

MD04 29 -- 103

python3 run_emerge_iterate_multiCore.py MD04 0
python3 run_emerge_iterate_multiCore.py MD04 1
python3 run_emerge_iterate_multiCore.py MD04 2
python3 run_emerge_iterate_multiCore.py MD04 3
python3 run_emerge_iterate_multiCore.py MD04 4
python3 run_emerge_iterate_multiCore.py MD04 5
python3 run_emerge_iterate_multiCore.py MD04 6
python3 run_emerge_iterate_multiCore.py MD04 7
python3 run_emerge_iterate_multiCore.py MD04 8
python3 run_emerge_iterate_multiCore.py MD04 9


hlist_0.23800_emerge.hdf5


import numpy as n 
import sys
import time

import EmergeIterate

t0=time.time()

iterate = EmergeIterate.EmergeIterate(22, 'MD10')
iterate.open_snapshots()
iterate.init_new_quantities()
#iterate.map_halos_between_snapshots()

import pandas as pd

print('N f1 halos', len(iterate.f1['/halo_properties/id'].value))
print('N f0 halos', len(iterate.f0['/halo_properties/id'].value))
f1_new_halos = (n.in1d(iterate.f1['/halo_properties/id'].value, iterate.f0['/halo_properties/desc_id'].value)==False)
# new halos [f1_new_halos] to be fed to EmergeIterate.compute_qtys_new_halos_pk
 
f1_evolved_halos = (new_halos==False)

f0_propagated_halos = n.in1d(iterate.f0['/halo_properties/desc_id'].value, iterate.f1['/halo_properties/id'].value)
f0_lost_halos = (propagated_halos==False)

print('lost halos', len(iterate.f0['/halo_properties/desc_id'].value[f0_lost_halos]))
print('propagated halos', len(iterate.f0['/halo_properties/desc_id'].value[f0_propagated_halos]))

print('new halos', len(iterate.f1['/halo_properties/desc_id'].value[f1_new_halos]))
print('evolved halos', len(iterate.f1['/halo_properties/id'].value[f1_evolved_halos]))


# parmi f0_propagated_halos
# cherchez les desc_id qui apparaissent plusieurs fois
s = pd.Series(iterate.f0['/halo_properties/desc_id'].value[f0_propagated_halos])
duplicated_desc_ids = s[s.duplicated()].get_values()
f0_with_desc_id_duplicates = n.in1d( iterate.f0['/halo_properties/desc_id'].value,  duplicated_desc_ids )
f0_in_a_merging = f0_with_desc_id_duplicates & f0_propagated_halos
f0_not_merging = (f0_with_desc_id_duplicates == False) & f0_propagated_halos

print('merging halos', len(iterate.f0['/halo_properties/id'].value[f0_in_a_merging]))
print('not merging halos', len(iterate.f0['/halo_properties/id'].value[f0_not_merging]))

# mapper les 12964 sur les 12981 avec selection tel que
# iterate.f1['/halo_properties/id'].value[f1_evolved_halos & selection] = iterate.f0['/halo_properties/desc_id'].value[f0_not_merging]
# comme les indices sont ordonnes, il suffit de retirer ceux qui apparaissent plusieurs fois
# selection = mask_f1_in_a_merging == False
mask_f1_in_a_merging = n.in1d( iterate.f1['/halo_properties/id'].value,  duplicated_desc_ids )
f1_evolved_halos_no_merger = (f1_evolved_halos) & (mask_f1_in_a_merging==False)
f1_evolved_halos_with_merger = (f1_evolved_halos) & (mask_f1_in_a_merging)

# evolving halos to be fed to EmergeIterate.compute_qtys_evolving_halos_pk
# 1:1 with f0_not_merging <=> f1_evolved_halos_no_merger
print("not merging halos in f0 and f1",len(f0_not_merging.nonzero()[0]), len(f1_evolved_halos_no_merger.nonzero()[0]))

# merging halos to be fed to EmergeIterate.get_position_merger_players_pk
# X=>1 with f0_in_a_merging => f1_evolved_halos_with_merger
print("merging halos in f0 and f1",len(f0_in_a_merging.nonzero()[0]), len(f1_evolved_halos_with_merger.nonzero()[0])))
# 34 => 17
iterate.f0['/halo_properties/desc_id'].value[f0_in_a_merging]
iterate.f0['/halo_properties/desc_id'].value[f0_in_a_merging].reshape(17,2)
iterate.f0['/halo_properties/id'].value[f0_in_a_merging]
iterate.f0['/halo_properties/Future_merger_MMP_ID'].value[f0_in_a_merging]

iterate.f0['/halo_properties/mvir'].value[f0_in_a_merging]

merger_ids = iterate.f1['/halo_properties/id'].value[f1_evolved_halos_with_merger]

for id_host_f1, merger_id in enumerate(merger_ids):
	f0_corr_mid = (iterate.f0['/halo_properties/desc_id'].value[f0_in_a_merging] == merger_id)
	id_host_f0 = n.argmax(iterate.f0['/halo_properties/mvir'].value[f0_in_a_merging][f0_corr_mid])
	id_guests_f0 = n.delete(n.arange(len(iterate.f0['/halo_properties/mvir'].value[f0_in_a_merging][f0_corr_mid])), id_host_f0)
	print(merger_id, id_host_f0, id_guests_f0)
    # apply those indexes as follow on f0 and f1 
	#f1 [f1_evolved_halos_with_merger][merger_id]
	#f0 [f0_in_a_merging][f0_corr_mid][id_host_f0]
	#f0 [f0_in_a_merging][f0_corr_mid][id_guests_f0]
	sum_stellar_mass_guests = n.sum(iterate.f0['/emerge_data/stellar_mass'].value[f0_in_a_merging][f0_corr_mid][id_guests_f0])

EmergeIterate.merge_system(mvir_f0, mvir_f1, age_f0, age_f1, rvir_f0, rvir_f1, redshift, t_dynamical, rs_f1, mpeak_f1, mpeak_scale_f1,  f1_scale, m_icm_f0, stellar_mass_f0, star_formation_rate_f0, sum_stellar_mass_guests)


f0_is_evolving = (iterate.f0['/halo_properties/Future_merger_MMP_ID'].value[f0_in_a_merging] == iterate.f0['/halo_properties/id'].value[f0_in_a_merging])
f0_is_merging = (iterate.f0['/halo_properties/Future_merger_MMP_ID'].value[f0_in_a_merging] != iterate.f0['/halo_properties/id'].value[f0_in_a_merging])

print('is evolving', len(f0_is_evolving.nonzero()[0]))
print('is merging', len(f0_is_merging.nonzero()[0]))

f1_evolved_halos_with_merger

#test = iterate.f1['/halo_properties/id'].value[(f1_evolved_halos )&(mask_f1_in_a_merging==False)] == iterate.f0['/halo_properties/desc_id'].value[f0_not_merging]
#print(len(iterate.f1['/halo_properties/id'].value[(f1_evolved_halos )&(mask_f1_in_a_merging==False)][test]))

iterate.f0['/halo_properties/desc_id'].value[f0_not_merging]
iterate.f1['/halo_properties/id'].value[f1_evolved_halos]

iterate.f0['/halo_properties/desc_id'].value[f0_in_a_merging]

iterate.f1_id_with_multiple_progenitors = s[s.duplicated()].get_values()

iterate.mask_f1_in_a_merging = n.in1d( iterate.f1['/halo_properties/id'].value,  iterate.f1_id_with_multiple_progenitors )

f0_in_a_merging = n.in1d( iterate.f0['/halo_properties/desc_id'].value,  duplicated_desc_ids ) & f0_propagated_halos




f0_is_evolving = (iterate.f0['/halo_properties/Future_merger_MMP_ID'].value == iterate.f0['/halo_properties/id'].value)&(propagated_halos)
f0_is_merging = (iterate.f0['/halo_properties/Future_merger_MMP_ID'].value != iterate.f0['/halo_properties/id'].value)&(propagated_halos)

print('is evolving', len(f0_is_evolving.nonzero()[0]))
print('is merging', len(f0_is_merging.nonzero()[0]))


f1_id_unique_list_descendents_detected_at_next_scale = n.intersect1d(n.unique(iterate.f0['/halo_properties/desc_id'].value), iterate.f1['/halo_properties/id'].value)
mask_f0_to_propagate = n.in1d(iterate.f0['/halo_properties/desc_id'].value, f1_id_unique_list_descendents_detected_at_next_scale)
mask_f0_lost = (mask_f0_to_propagate == False )

# evolving halos are given after applying this boolean mask to a f1 quantity :
mask_f1_evolved_from_previous = n.in1d( iterate.f1['/halo_properties/id'].value,  f1_id_unique_list_descendents_detected_at_next_scale )

# new halos are given after applying this boolean mask to a f1 quantity 
# new halos in f1, not present in f0
iterate.mask_f1_new_halos = (mask_f1_evolved_from_previous==False)
print('new halos', len(iterate.mask_f1_new_halos.nonzero()[0]))
# f1 ids not in f0 desc id

# f1 ids present in desc id
# in the following list, ids only appear once 
all_ids_f1 = iterate.f1['/halo_properties/id'].value[mask_f1_evolved_from_previous]
#total number = N_f1

# all f0 halos that become halos in f1
#mask_f0_to_propagate
#total number : N_f0 > N_f1
 
# a single appearance in desc id => evolving halo 
# multiple appearance in desc id => merging halos

# halos descending :
# mask_f0_to_propagate
# mask_f1_evolved_from_previous
s = pd.Series(iterate.f0['/halo_properties/desc_id'].value[mask_f0_to_propagate])

iterate.f1_id_with_multiple_progenitors = s[s.duplicated()].get_values()

# also = f0_desc_id merging into 1 halo in f1
# merging systems [many halos in f0 into a single f1 halo]
iterate.mask_f1_in_a_merging = n.in1d( iterate.f1['/halo_properties/id'].value,  iterate.f1_id_with_multiple_progenitors )
iterate.mask_f0_in_a_merging = n.in1d( iterate.f0['/halo_properties/desc_id'].value,  iterate.f1_id_with_multiple_progenitors )

# halos mapped 1::1 between snapshots	
iterate.mask_f0_evolving_11_halos = ( mask_f0_to_propagate ) & ( iterate.mask_f0_in_a_merging == False )
iterate.mask_f1_evolving_11_halos = ( mask_f1_evolved_from_previous ) & ( iterate.mask_f1_in_a_merging == False )
print('11 mapping', len(iterate.mask_f0_evolving_11_halos.nonzero()[0]), len(iterate.mask_f1_evolving_11_halos.nonzero()[0]))
print('merging systems', len(iterate.f1_id_with_multiple_progenitors))

s = pd.Series(iterate.f0['/halo_properties/desc_id'].value[mask_f0_to_propagate])
iterate.f1_id_with_multiple_progenitors = s[s.duplicated()].get_values()
# also = f0_desc_id merging into 1 halo in f1
# merging systems [many halos in f0 into a single f1 halo]
iterate.mask_f1_in_a_merging = n.in1d( iterate.f1['/halo_properties/id'].value,  iterate.f1_id_with_multiple_progenitors )

iterate.mask_f0_in_a_merging = n.in1d( iterate.f0['/halo_properties/desc_id'].value,  iterate.f1_id_with_multiple_progenitors )

# now we need to figure out at f0 what are hosts and what are merging entities

positions_f0 = iterate.positions_f0[iterate.mask_f0_in_a_merging]
ids_f0 = iterate.f0['/halo_properties/desc_id'].value[iterate.mask_f0_in_a_merging]
desc_id_f0 = iterate.f0['/halo_properties/id'].value[iterate.mask_f0_in_a_merging]
Future_merger_MMP_ID_f0 = iterate.f0['/halo_properties/Future_merger_MMP_ID'].value[iterate.mask_f0_in_a_merging]


iterate.f1['/halo_properties/id'].value[ iterate.mask_f1_in_a_merging ]

id_f0 = iterate.f0['/halo_properties/id'].value[ iterate.mask_f0_in_a_merging ]
desc_id_f0 = iterate.f0['/halo_properties/desc_id'].value[ iterate.mask_f0_in_a_merging ]
Future_merger_MMP_ID_f0 = iterate.f0['/halo_properties/Future_merger_MMP_ID'].value[ iterate.mask_f0_in_a_merging ]

is_evolving = (iterate.f0['/halo_properties/Future_merger_MMP_ID'].value == iterate.f0['/halo_properties/id'].value)
is_merging = (iterate.f0['/halo_properties/Future_merger_MMP_ID'].value != iterate.f0['/halo_properties/id'].value)

len(is_evolving.nonzero()[0])
len(is_merging.nonzero()[0])

mask2_f0_in_a_merging_hosts = (id_f0 == Future_merger_MMP_ID_f0)

mask2_f0_in_a_merging_guest = (id_f0 != Future_merger_MMP_ID_f0)

desc_id_f0[mask2_f0_in_a_merging_hosts]

# hosts of a merger
ids_f1 = iterate.f1['/halo_properties/id'].value[iterate.mask_f1_in_a_merging]
positions_f1 = iterate.positions_f1[iterate.mask_f1_in_a_merging]

