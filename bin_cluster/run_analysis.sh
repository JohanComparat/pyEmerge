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





#########################################
# LIGHT CONES
#########################################
# create the shells of the light cone

sh lc_create_shells_run.sh
sh lc_create_shells_run_z1.sh
sh lc_create_shells_run_L6.sh
sh lc_create_shells_run_L15.sh
# remains a problem for L6, some logSAR are not written
# based on  lc_create_shells.py

# merges the shells into a single light cone file
python3 lc_merge_shells.py L3
python3 lc_merge_shells.py L3_z1
python3 lc_merge_shells.py L6
python3 lc_merge_shells.py L15

# Adds ra, dec, z 
python3 lc_add_sky_coordinates.py remaped_position_L6
python3 lc_add_sky_coordinates.py remaped_position_L3
python3 lc_add_sky_coordinates.py remaped_position_L3_z1
python3 lc_add_sky_coordinates.py remaped_position_L15
 
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
python3 lc_add_agns.py /data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_remaped_position_L6.hdf5
python3 lc_add_agns.py /data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_remaped_position_L3.hdf5
python3 lc_add_agns.py /data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_remaped_position_L3_z1.hdf5
python3 lc_add_agns.py /data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_remaped_position_L15.hdf5

python3 lc_lognlogs_agns.py
# results are shown here : 
# http://www.mpe.mpg.de/~comparat/eRoMok/logNlogS/
# logs logn is still too bright by 0.6 dex. Probably an issue when applying the NH attenuation.

# create fits files for each light cone. Modify the script according to which LC you want to generate
python3 lc_convert_2_fits.py


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

# Adds cluster related columns 
python3 lc_add_clusters.py  
# Adds galaxy related columns 
python3 lc_add_galaxies.py  




