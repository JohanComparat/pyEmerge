"""
Routine to create the light cones shells

L1 L2 L3   u11 u12 u13   u21 u22 u23   u31 u32 u33   (periodicity)
'5.9161', '0.4140', '0.4082', '5', '3', '1', '1', '1', '0', '0', '1', '0', '(1)'
'2.4495', '0.7071', '0.5774', '2', '1', '1', '1', '1', '0', '0', '1', '0', '(1)'

python3 create_light_cone_shells.py 10 MD10 1000


python3 create_light_cone_shells.py 10 MD10 1000

import numpy as n
import os
for ii in n.arange(50,115,1)[::-1]:
	comm="python3 create_light_cone_shells.py "+str(ii)+" MD10 1000"	
	print(comm)
	os.system(comm)

"""

import sys
env = 'MD10'# sys.argv[1] # 'MD10'
L_box = 1000./0.6777 # float(sys.argv[2]) / 0.6777
remap = sys.argv[1]

if remap=="L3":
	positions_group_name = 'remaped_position_L3'
	#x_obs, y_obs, z_obs = 0., 0.7071/2.*L_box, 0.5774/2.*L_box

if remap=="L3_z1":
	positions_group_name = 'remaped_position_L3_z1'
	#x_obs, y_obs, z_obs = -2.4495*L_box, 0.7071/2.*L_box, 0.5774/2.*L_box

if remap=="L6":
        positions_group_name = 'remaped_position_L6'
        #x_obs, y_obs, z_obs = 0., 0.4140/2.*L_box, 0.4082/2.*L_box

if remap=="L15":
        positions_group_name = 'remaped_position_L15'
        #x_obs, y_obs, z_obs = 0., 0.4140/2.*L_box, 0.4082/2.*L_box


import h5py    # HDF5 support
import os
import glob
import numpy as n

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)

h5_lc_dir = os.path.join(os.environ[env], 'h5_lc', 'cluster_shells_'+positions_group_name )
input_list = n.array(glob.glob(os.path.join(h5_lc_dir, "shell_hlist_?.?????.hdf5")))
input_list.sort()

file_out = os.path.join(os.environ[env], 'h5_lc', 'lc_cluster_'+positions_group_name+'.hdf5')
if os.path.isfile(file_out):
	os.system('rm '+file_out)

print("creates redshift list")
# creates the redshift list 
aexp = []
for file_1 in input_list : 
	aexp.append(float(os.path.basename(file_1)[:-5].split('_')[-1]))

aexp = n.array(aexp)

# loop over the shell files and concatenates everything ? + shell redshift.

def get_data(ii):
	f = h5py.File(input_list[ii],  "r")
	DATA = n.array([
		n.ones_like( 
		f['/halo_position/x'].value)*(1./aexp[ii]-1.), 
		f['/halo_position/x'].value, 
		f['/halo_position/y'].value, 
		f['/halo_position/z'].value, 
		f['/halo_position/vx'].value, 
		f['/halo_position/vy'].value, 
		f['/halo_position/vz'].value, 
		# DERIVED PROPERTIES
		f['/moster_2013_data/stellar_mass'].value, 
		f['/cluster_data/cool_class'].value,
		f['/cluster_data/kT'].value,
		f['/cluster_data/log_LX_05_24'].value,
		f['/cluster_data/log_LceX_05_24'].value,
		f['/cluster_data/log_Mgas'].value,
		#  HALO PROPERTIES
		f['halo_properties/scale']                          .value,
		f['halo_properties/id']                             .value,
		f['halo_properties/desc_scale']                     .value,
		f['halo_properties/desc_id']                        .value,
		f['halo_properties/num_prog']                       .value,
		f['halo_properties/pid']                            .value,
		f['halo_properties/upid']                           .value,
		f['halo_properties/desc_pid']                       .value,
		f['halo_properties/mvir']                           .value,
		f['halo_properties/rvir']                           .value,
		f['halo_properties/rs']                             .value,
		f['halo_properties/vrms']                           .value,
		f['halo_properties/mmp']                            .value,
		f['halo_properties/scale_of_last_MM']               .value,
		f['halo_properties/vmax']                           .value,
		f['halo_properties/Jx']                             .value,
		f['halo_properties/Jy']                             .value,
		f['halo_properties/Jz']                             .value,
		f['halo_properties/Spin']                           .value,
		f['halo_properties/Breadth_first_ID']               .value,
		f['halo_properties/Depth_first_ID']                 .value,
		f['halo_properties/Tree_root_ID']                   .value,
		f['halo_properties/Orig_halo_ID']                   .value,
		f['halo_properties/Next_coprogenitor_depthfirst_ID'].value,
		f['halo_properties/Last_progenitor_depthfirst_ID']  .value,
		f['halo_properties/Last_mainleaf_depthfirst_ID']    .value,
		f['halo_properties/Tidal_Force']                    .value,
		f['halo_properties/Tidal_ID']                       .value,
		f['halo_properties/Rs_Klypin']                      .value,
		f['halo_properties/Mmvir_all']                      .value,
		f['halo_properties/M200b']                          .value,
		f['halo_properties/M200c']                          .value,
		f['halo_properties/M500c']                          .value,
		f['halo_properties/M2500c']                         .value,
		f['halo_properties/Xoff']                           .value,
		f['halo_properties/Voff']                           .value,
		f['halo_properties/Spin_Bullock']                   .value,
		f['halo_properties/b_to_a']                         .value,
		f['halo_properties/c_to_a']                         .value,
		f['halo_properties/Ax']                             .value,
		f['halo_properties/Ay']                             .value,
		f['halo_properties/Az']                             .value,
		f['halo_properties/b_to_a_500c']                    .value,
		f['halo_properties/c_to_a_500c']                    .value,
		f['halo_properties/Ax_500c']                        .value,
		f['halo_properties/Ay_500c']                        .value,
		f['halo_properties/Az_500c']                        .value,
		f['halo_properties/TU']                             .value,
		f['halo_properties/M_pe_Behroozi']                  .value,
		f['halo_properties/M_pe_Diemer']                    .value,
		f['halo_properties/Macc']                           .value,
		f['halo_properties/Mpeak']                          .value,
		f['halo_properties/Vacc']                           .value,
		f['halo_properties/Vpeak']                          .value,
		f['halo_properties/Halfmass_Scale']                 .value,
		f['halo_properties/Acc_Rate_Inst']                  .value,
		f['halo_properties/Acc_Rate_100Myr']                .value,
		f['halo_properties/Acc_Rate_1Tdyn']                 .value,
		f['halo_properties/Acc_Rate_2Tdyn']                 .value,
		f['halo_properties/Acc_Rate_Mpeak']                 .value,
		f['halo_properties/Mpeak_Scale']                    .value,
		f['halo_properties/Acc_Scale']                      .value,
		f['halo_properties/First_Acc_Scale']                .value,
		f['halo_properties/First_Acc_Mvir']                 .value,
		f['halo_properties/First_Acc_Vmax']                 .value,
		f['halo_properties/VmaxAtMpeak']                    .value,
		f['halo_properties/Tidal_Force_Tdyn']               .value,
		f['halo_properties/logVmaxVmaxmaxTdynTmpeak']       .value,
		f['halo_properties/Time_to_future_merger']          .value,
		f['halo_properties/Future_merger_MMP_ID']           .value
		])
	f.close()
	return DATA

print(input_list)
DATA = get_data(0)
print(DATA.shape)
for ii in range(1,len(input_list),1):
	print(input_list[ii])
	DATA = n.hstack((DATA, get_data(ii)))


f = h5py.File(file_out, "a")
f.attrs['file_name'] = os.path.basename(file_out)
f.attrs['HDF5_Version']     = h5py.version.hdf5_version
f.attrs['h5py_version']     = h5py.version.version

halo_data = f.create_group('halo_position')

ds = halo_data.create_dataset('z_snap', data = DATA[0] )
ds = halo_data.create_dataset('x', data = DATA[1] )
ds = halo_data.create_dataset('y', data = DATA[2] )
ds = halo_data.create_dataset('z', data = DATA[3] )
ds = halo_data.create_dataset('vx', data = DATA[4] )
ds = halo_data.create_dataset('vy', data = DATA[5] )
ds = halo_data.create_dataset('vz', data = DATA[6] )

halo_data = f.create_group('moster_2013_data')
ds = halo_data.create_dataset('stellar_mass', data = DATA[7] )
	

halo_data = f.create_group('cluster_data')
ds = halo_data.create_dataset('cool_class',     data = DATA[8]  )
ds = halo_data.create_dataset('kT',             data = DATA[9]  )
ds = halo_data.create_dataset('log_LX_05_24',   data = DATA[10] )
ds = halo_data.create_dataset('log_LceX_05_24', data = DATA[11] )
ds = halo_data.create_dataset('log_Mgas',       data = DATA[12] )

halo_data = f.create_group('halo_properties')

ds = halo_data.create_dataset('scale'                              , data=DATA[13])
ds = halo_data.create_dataset('id'                                 , data=DATA[14])
ds = halo_data.create_dataset('desc_scale'                         , data=DATA[15])
ds = halo_data.create_dataset('desc_id'                            , data=DATA[16])
ds = halo_data.create_dataset('num_prog'                           , data=DATA[17])
ds = halo_data.create_dataset('pid'                                , data=DATA[18])
ds = halo_data.create_dataset('upid'                               , data=DATA[19])
ds = halo_data.create_dataset('desc_pid'                           , data=DATA[20])
ds = halo_data.create_dataset('mvir'                               , data=DATA[21])
ds = halo_data.create_dataset('rvir'                               , data=DATA[22])
ds = halo_data.create_dataset('rs'                                 , data=DATA[23])
ds = halo_data.create_dataset('vrms'                               , data=DATA[24])
ds = halo_data.create_dataset('mmp'                                , data=DATA[25])
ds = halo_data.create_dataset('scale_of_last_MM'                   , data=DATA[26])
ds = halo_data.create_dataset('vmax'                               , data=DATA[27])
ds = halo_data.create_dataset('Jx'                                 , data=DATA[28])
ds = halo_data.create_dataset('Jy'                                 , data=DATA[29])
ds = halo_data.create_dataset('Jz'                                 , data=DATA[30])
ds = halo_data.create_dataset('Spin'                               , data=DATA[31])
ds = halo_data.create_dataset('Breadth_first_ID'                   , data=DATA[32])
ds = halo_data.create_dataset('Depth_first_ID'                     , data=DATA[33])
ds = halo_data.create_dataset('Tree_root_ID'                       , data=DATA[34])
ds = halo_data.create_dataset('Orig_halo_ID'                       , data=DATA[35])
ds = halo_data.create_dataset('Next_coprogenitor_depthfirst_ID'    , data=DATA[36])
ds = halo_data.create_dataset('Last_progenitor_depthfirst_ID'      , data=DATA[37])
ds = halo_data.create_dataset('Last_mainleaf_depthfirst_ID'        , data=DATA[38])
ds = halo_data.create_dataset('Tidal_Force'                        , data=DATA[39])
ds = halo_data.create_dataset('Tidal_ID'                           , data=DATA[40])
ds = halo_data.create_dataset('Rs_Klypin'                          , data=DATA[41])
ds = halo_data.create_dataset('Mmvir_all'                          , data=DATA[42])
ds = halo_data.create_dataset('M200b'                              , data=DATA[43])
ds = halo_data.create_dataset('M200c'                              , data=DATA[44])
ds = halo_data.create_dataset('M500c'                              , data=DATA[45])
ds = halo_data.create_dataset('M2500c'                             , data=DATA[46])
ds = halo_data.create_dataset('Xoff'                               , data=DATA[47])
ds = halo_data.create_dataset('Voff'                               , data=DATA[48])
ds = halo_data.create_dataset('Spin_Bullock'                       , data=DATA[49])
ds = halo_data.create_dataset('b_to_a'                             , data=DATA[50])
ds = halo_data.create_dataset('c_to_a'                             , data=DATA[51])
ds = halo_data.create_dataset('Ax'                                 , data=DATA[52])
ds = halo_data.create_dataset('Ay'                                 , data=DATA[53])
ds = halo_data.create_dataset('Az'                                 , data=DATA[54])
ds = halo_data.create_dataset('b_to_a_500c'                        , data=DATA[55])
ds = halo_data.create_dataset('c_to_a_500c'                        , data=DATA[56])
ds = halo_data.create_dataset('Ax_500c'                            , data=DATA[57])
ds = halo_data.create_dataset('Ay_500c'                            , data=DATA[58])
ds = halo_data.create_dataset('Az_500c'                            , data=DATA[59])
ds = halo_data.create_dataset('TU'                                 , data=DATA[60])
ds = halo_data.create_dataset('M_pe_Behroozi'                      , data=DATA[61])
ds = halo_data.create_dataset('M_pe_Diemer'                        , data=DATA[62])
ds = halo_data.create_dataset('Macc'                               , data=DATA[63])
ds = halo_data.create_dataset('Mpeak'                              , data=DATA[64])
ds = halo_data.create_dataset('Vacc'                               , data=DATA[65])
ds = halo_data.create_dataset('Vpeak'                              , data=DATA[66])
ds = halo_data.create_dataset('Halfmass_Scale'                     , data=DATA[67])
ds = halo_data.create_dataset('Acc_Rate_Inst'                      , data=DATA[68])
ds = halo_data.create_dataset('Acc_Rate_100Myr'                    , data=DATA[69])
ds = halo_data.create_dataset('Acc_Rate_1Tdyn'                     , data=DATA[70])
ds = halo_data.create_dataset('Acc_Rate_2Tdyn'                     , data=DATA[71])
ds = halo_data.create_dataset('Acc_Rate_Mpeak'                     , data=DATA[72])
ds = halo_data.create_dataset('Mpeak_Scale'                        , data=DATA[73])
ds = halo_data.create_dataset('Acc_Scale'                          , data=DATA[74])
ds = halo_data.create_dataset('First_Acc_Scale'                    , data=DATA[75])
ds = halo_data.create_dataset('First_Acc_Mvir'                     , data=DATA[76])
ds = halo_data.create_dataset('First_Acc_Vmax'                     , data=DATA[77])
ds = halo_data.create_dataset('VmaxAtMpeak'                        , data=DATA[78])
ds = halo_data.create_dataset('Tidal_Force_Tdyn'                   , data=DATA[79])
ds = halo_data.create_dataset('logVmaxVmaxmaxTdynTmpeak'           , data=DATA[80])
ds = halo_data.create_dataset('Time_to_future_merger'              , data=DATA[81])
ds = halo_data.create_dataset('Future_merger_MMP_ID'               , data=DATA[82])

f.close()
