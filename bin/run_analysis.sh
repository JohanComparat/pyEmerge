# Running sequence to obtain EMERGE galaxies in the MultiDark simulation

module load git
module load anaconda
export PYTHONPATH="${PYTHONPATH}:/u/joco/software/pyEmerge/python/"

#########################################
# snapshot pipeline
#########################################

# download all hlist files generated by rockstar
ls /ptmp/joco/MD/MD_1.0Gpc/hlists/hlist_?.?????.list

# generate smaller files without all the unnecessary columns
# saves only halos with mvir > M of 100 particles
# the columns saved are limited to the stricti minimum ;
#  ('id'                     , 'i4' )
# ,('desc_scale'             , 'f4' )
# ,('desc_id'                , 'i4' )
# ,('pid'                    , 'i4' )
# ,('mvir'                   , 'f4' )
# ,('rvir'                   , 'f4' )
# ,('rs'                     , 'f4' )
# ,('Mpeak'                  , 'f4' )
# ,('Mpeak_scale'            , 'f4' )
# ,('Acc_Rate_1Tdyn'         , 'f4' )
# ,('Time_to_future_merger'  , 'f4' )
# ,('Future_merger_MMP_ID'   , 'i4' )
python generate_and_submit_MD10-write-smallFile.py  
# it generates a job per hlist file and submits it to the queue at RZG
# it uses the script :
MD10-write-smallFile.py 
# it writes files here :
ls /ptmp/joco/MD/MD_1.0Gpc/emerge/hlist_?.?????.data

# generate hdf5 files to access easily the data

f['/halo_data/id'].value

# f.attrs['aexp']      = float(aexp)#[snap_id]
# f.attrs['redshift']  = float(redshift)#[snap_id]
# f.attrs['age_yr']    = float(age_yr)#[snap_id] 
# f.attrs['rho_crit']  = float(rho_crit)#[snap_id] 
# f.attrs['delta_vir'] = float(delta_vir)#[snap_id]
# 
# halo_data = f.create_group('halo_data')
# halo_data.attrs['N_halos'] =  N_halo

python generate_and_submit_convert-2-h5.py  
# it generates a job per hlist file and submits it to the queue at RZG
# it uses the script :
convert-2-h5.py  
# it writes files here :
ls -lh /ptmp/joco/MD/MD_1.0Gpc/h5/hlist_?.?????_emerge.hdf5

# initialize the emerge columns in the files 
sh h5_init_emerge_data_with_0s_run_md04.sh  
sh h5_init_emerge_data_with_0s_run_md10.sh
# using the script
h5_init_emerge_data_with_0s.py  


# add the emerge information: star formation rate, stellar mass and icm mass.
# all following snapshots on ds52 this way :
sh run_emerge_iterate_multiCore_run_md04.sh  
sh run_emerge_iterate_multiCore_run_md10.sh
# using this script
run_emerge_iterate_multiCore.py

# remap the coordinates to have a cuboid
sh run_remap_MD04.s
sh run_remap_MD10.s
# executes the remapping for the MD04 box
# 12 cores is the maximum on ds52
# snapshots 72, 76, 98, 99, 111 fail at this step
# BLACKLIST REMAP FAILs :
# hlist_0.27060_emerge.hdf5
# hlist_0.43090_emerge.hdf5
# hlist_0.71730_emerge.hdf5
# hlist_0.93570_emerge.hdf5

# Revised implementaion of the Moster et al. 2013 + revised parameters to match Planck cosmology.
# mass + agn functions are ok until redshift 3, then it is wrong.
sh run_MSMo13_md10.sh
# uses snapshots_add_Ms_Mo13.py

# DATA MODEL FOR SNAPSHOTS
python3 print_data_structure.py 22 MD10
python3 print_data_structure.py 22 MD04

# stellar mass function
python3 measure_SMF.py
python3 measure_SMF_Mo13.py
# model for the host galaxy mass function for AGN
#tabulates the duty cycles and AGN host galaxy function for the Mo13 implementation
python3 tabulate_HGMF_per_snapshot.py

sh run_LSAR_md10.sh
# uses snapshots_add_LSAR_Bo16.py
# BLACKLIST LSAR FAILS
# hlist_0.28920_emerge.hdf5    

970000 2956779 332.835337638855
Traceback (most recent call last):
  File "snapshots_add_LSAR_Bo16.py", line 69, in <module>
    log_lSAR[ii0:ii1] = log_lambda_SAR_values[n.array([n.min(n.where(cmat.T[jj]==True)) for jj in n.arange(len(cmat.T)) ])]
  File "snCapshots_add_LSAR_Bo16.py", line 69, in <listcomp>
    log_lSAR[ii0:ii1] = log_lambda_SAR_values[n.array([n.min(n.where(cmat.T[jj]==True)) for jj in n.arange(len(cmat.T)) ])]
  File "/home/comparat/.local/lib/python3.4/site-packages/numpy/core/fromnumeric.py", line 2372, in amin
    out=out, **kwargs)
  File "/home/comparat/.local/lib/python3.4/site-packages/numpy/core/_methods.py", line 29, in _amin
    return umr_minimum(a, axis, None, out, keepdims)
ValueError: zero-size array to reduction operation minimum which has no identity

sh run_AGN_activity_md10.sh
# uses snapshots_add_AGN_activity_Bo16.py

# star formation rate density. WAS FOR EMERGE. NOT USED ANYMORE.
#python measure_SFRD.py
#python plot_SFRD.py

python3 plot_SMF.py
# plots are here: http://www.mpe.mpg.de/~comparat/eRoMok/h5/stellar_mass_function/

#########################################
# LIGHT CONES
#########################################
# create the shells of the light cone

sh lc_create_shells_run.sh
sh lc_create_shells_run_z1.sh
sh lc_create_shells_run_L6.sh
# remains a problem for L6, some logSAR are not written
# based on  lc_create_shells.py

# merges the shells into a single light cone file
python3 lc_merge_shells.py L3
python3 lc_merge_shells.py L3_z1
python3 lc_merge_shells.py L6

# Adds ra, dec, z 
python3 lc_add_sky_coordinates.py remaped_position_L6
python3 lc_add_sky_coordinates.py remaped_position_L3
python3 lc_add_sky_coordinates.py remaped_position_L3_z1
 
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

# converts the Bongiorno luminosity into the eRosita band assuming a NH distribution and following Buchner et al. 2017.
# 
python3 lc_add_agns.py /data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_remaped_position_L6.hdf5
python3 lc_add_agns.py /data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_remaped_position_L3.hdf5
python3 lc_add_agns.py /data17s/darksim/MD/MD_1.0Gpc/h5_lc/lc_remaped_position_L3_z1.hdf5

python3 lc_lognlogs_agns.py
# results are shown here : 
# http://www.mpe.mpg.de/~comparat/eRoMok/logNlogS/
# logs logn is still too bright by 0.6 dex. Probably an issue when applying the NH attenuation.

# now creates clustering catalogs equivalent to SPIDERS
python3 create_random_ra_dec.py
python3 lc_write_clustering_sample.py
# and compute the clustering 
cd /data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_remaped_position_L3/
sh  run_clustering.sh

# Adds cluster related columns 
python3 lc_add_clusters.py  
# Adds galaxy related columns 
python3 lc_add_galaxies.py  




