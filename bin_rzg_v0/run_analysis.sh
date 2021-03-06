# Running sequence to obtain EMERGE galaxies in the MultiDark simulation

module load git
module load anaconda
export PYTHONPATH="${PYTHONPATH}:/u/joco/software/pyEmerge/python/"

#preliminary pipeline
#--------------------

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
python emerge_init.py # first snapshot for each simulation
# all following snapshots on ds52 this way :
sh run_emerge_iterate_multiCore_run_md04.sh  
sh run_emerge_iterate_multiCore_run_md10.sh
# using this script
run_emerge_iterate_multiCore.py

# alternate implementation that does note really work fine
# all following snapshots on hydra this way :
python generate_and_submit_emerge_iterate.py
cd /u/joco/batch_emerge
llsubmit emerge_iterate_batch.sh
# it uses the command :
python emerge_iterate.py $ID

# remap the coordinates to have a cuboid
sh run_remap_MD04.s
sh run_remap_MD10.s
# executes the remapping for the MD04 box
# 12 cores is the maximum on ds52
python3 test_remap.py ARGS

#########################################
# DATA MODEL
#########################################

python3 print_data_structure.py 22 MD10
python3 print_data_structure.py 22 MD04


# Add clusters, agns, galaxies, from 4MOST.




#########################################
#########################################
#########################################
#########################################

read_h5_test.sh  
submit_test.sh


# machines for dev :
connect to rzgi (interactive) from mpe
module load anaconda
ipython
