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

# Revised implementaion of the Moster et al. 2013 + revised parameters to match Planck cosmology.
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
python3 plot_SMF.py

sh run_LSAR_md10.sh
# uses snapshots_add_LSAR_Bo16.py
#ONGOING
#ONGOING
#ONGOING
# Some value error for snapshot 37, ... occuring 

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

# star formation rate density
#python measure_SFRD.py
#python plot_SFRD.py

#########################################
# LIGHT CONES
#########################################
# create the shells of the light cone

sh lc_create_shells_run.sh
sh lc_create_shells_run_L6.sh
# remains a problem for L6, some logSAR are not written
# based on  lc_create_shells.py

# merges the shells into a single light cone file
python3 lc_merge_shells.py L3
python3 lc_merge_shells.py L6

# Adds ra, dec, z 
python3 lc_add_sky_coordinates.py remaped_position_L6
python3 lc_add_sky_coordinates.py remaped_position_L3
 
L3 characteristics :
z< 1.0889947373832305 |ra [deg]|< 6.7529257176359 |dec [deg]|< 8.269819492449505
N points: 8037075 


# TO WRITE NEXT : CODE IS A PLACE HOLDER
# Adds AGN related columns: ADD DUTY CYCLE and ACTIVE FRACTION.
# CHECK LOGNLOGS IS CORRECT
python3 lc_add_agns.py  
# Adds cluster related columns 
python3 lc_add_clusters.py  
# Adds galaxy related columns 
python3 lc_add_galaxies.py  




