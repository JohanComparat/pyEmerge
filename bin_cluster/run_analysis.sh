# Running sequence cluster in the simulations. Galaxies in the clusters are in the other light cone.

cd $GIT_EMERGE/bin_cluster

# Revised implementaion of the Moster et al. 2013 + revised parameters to match Planck cosmology.
# mass + agn functions are ok until redshift 3, then it is wrong.
sh run_MSMo13_md10.sh
# uses snapshots_add_Ms_Mo13.py

# splits the cluster population in 4 categories according on how disturbed they are.
python3.4 tabulate_xoff_distribution.py
sh run_add_coolcore.sh
# uses snapshots_add_coolness.py

# adds remaped coordinates
sh run_remap_MD10_f15.sh
# uses remap_lc.py

# adds Xray using basic scaling relations Matz 2016 or Dietrich 2017
sh run_xray.sh

# create light cone shells
sh lc_create_shells_run.sh
# uses lc_create_shells.py

# merges the shells into a single light cone file
python3.4 lc_merge_shells.py L3
python3.4 lc_merge_shells.py L3_z1
python3.4 lc_merge_shells.py L6
python3.4 lc_merge_shells.py L15

# Adds ra, dec, z 
python3.4 lc_add_sky_coordinates.py remaped_position_L6
python3.4 lc_add_sky_coordinates.py remaped_position_L3
python3.4 lc_add_sky_coordinates.py remaped_position_L3_z1
python3.4 lc_add_sky_coordinates.py remaped_position_L15
 
# L3 characteristics :
# z< 1.0889947373832305 |ra [deg]|< 6.7529257176359 |dec [deg]|< 8.269819492449505
# N points: 8037075 
# 
# L3_z1 characteristics
# z< 3.8309961826584344 |ra [deg]|< 3.3764628588325674 |dec [deg]|< 4.134909746242654
# N points: 8511571
# 
# L6 characteristics
# z< 6.697087333514605 |ra [deg]|< 1.9766516114702513 |dec [deg]|< 2.0047373031569915
# N points: 3287299
#
# L15 characteristics
# z< 0.5423857379098544 |ra [deg]|< 14.323944878104827 |dec [deg]|< 20.257311381848154
# N points: 8237518

# converts the Bongiorno luminosity into the eRosita band assuming a NH distribution and following Buchner et al. 2017.
# 
python3.4 lc_add_clusters.py /data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_cluster_remaped_position_L6.hdf5
python3.4 lc_add_clusters.py /data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_cluster_remaped_position_L3.hdf5
python3.4 lc_add_clusters.py /data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_cluster_remaped_position_L3_z1.hdf5
python3.4 lc_add_clusters.py /data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_cluster_remaped_position_L15.hdf5

python3.4 lc_lognlogs_clusters.py
# results are shown here : 
# http://www.mpe.mpg.de/~comparat/eRoMok/logNlogS/


python3.4 lc_identify_cluster_members.py L6
python3.4 lc_identify_cluster_members.py L3
python3.4 lc_identify_cluster_members.py L15
python3.4 lc_identify_cluster_members.py L3_z1

# identifies statistically the fraction of galaxies with higher age as a function of its distance to the cluster radius
python3.4 lc_add_red_sequence.py L3
python3.4 lc_add_red_sequence.py L3_z1
python3.4 lc_add_red_sequence.py L6
python3.4 lc_add_red_sequence.py L15

# create fits files for each light cone. Modify the script according to which LC you want to generate
python3.4 lc_convert_2_fits.py

# create a 4MOST mock by linking to templates :
python3.4 fits_add_templates.py

###################3
###################3
###################3

###################3
###################3
###################3

###################3
###################3
###################3

###################3
###################3###################3
# now creates clustering catalogs equivalent to SPIDERS
python3 lc_create_random_ra_dec.py
python3 lc_write_clustering_sample.py
python3 lc_write_CUTE_param_files.py

# and compute the clustering 
cd /data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_remaped_position_L3/
sh  run_clustering.sh
cd /data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_remaped_position_L6/
sh  run_clustering.sh
cd /data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_remaped_position_L15/
sh  run_clustering.sh



