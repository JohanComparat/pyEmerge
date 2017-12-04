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
ii = int(sys.argv[1])
env = sys.argv[2] # 'MD10'
L_box = float(sys.argv[3]) / 0.6777
positions_group_name = sys.argv[4] # 'remaped_position_L3'

if positions_group_name == 'remaped_position_L3' :
	positions_group = 'remaped_position_L3'
	x_obs, y_obs, z_obs = 0., 0.7071/2.*L_box, 0.5774/2.*L_box

if positions_group_name == 'remaped_position_L3_z1' :
	positions_group = 'remaped_position_L3'
	x_obs, y_obs, z_obs = -2.4495*L_box, 0.7071/2.*L_box, 0.5774/2.*L_box

if positions_group_name == 'remaped_position_L6' :
	positions_group = 'remaped_position_L6'
	x_obs, y_obs, z_obs = 0., 0.4140/2.*L_box, 0.4082/2.*L_box

if positions_group_name == 'remaped_position_L15' :
	positions_group = 'remaped_position_L15'
	#1.4142', '1.0000', '0.7071
	x_obs, y_obs, z_obs = 0., 1.0000/2.*L_box, 0.7071/2.*L_box


import h5py    # HDF5 support
import os
import glob
import numpy as n

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)

h5_lc_dir = os.path.join(os.environ[env], 'h5_lc', 'cluster_shells_'+positions_group_name )
if os.path.isdir(h5_lc_dir)==False:
	os.mkdir(h5_lc_dir)

h5_dir = os.path.join(os.environ[env], 'cluster_h5' )

input_list_i = n.array(glob.glob(os.path.join(h5_dir, "hlist_?.?????.hdf5")))
input_list_i.sort()

# removing snapshots that cannote be remapped ...
input_list = n.delete(input_list_i,n.array([
  #n.argwhere(input_list_i== os.path.join(h5_dir, "hlist_0.16620_emerge.hdf5")), # LSAR  issue
  #n.argwhere(input_list_i== os.path.join(h5_dir, "hlist_0.17770_emerge.hdf5")), # LSAR  issue
  #n.argwhere(input_list_i== os.path.join(h5_dir, "hlist_0.18990_emerge.hdf5")), # LSAR  issue
  #n.argwhere(input_list_i== os.path.join(h5_dir, "hlist_0.19410_emerge.hdf5")), # LSAR  issue
  n.argwhere(input_list_i== os.path.join(h5_dir, "hlist_0.21210_emerge.hdf5")), # LSAR  issue
  n.argwhere(input_list_i== os.path.join(h5_dir, "hlist_0.24230_emerge.hdf5")), # LSAR  issue
  n.argwhere(input_list_i== os.path.join(h5_dir, "hlist_0.28920_emerge.hdf5")), # LSAR  issue
  n.argwhere(input_list_i== os.path.join(h5_dir, "hlist_0.27060_emerge.hdf5")), # remap issue
  n.argwhere(input_list_i== os.path.join(h5_dir, "hlist_0.43090_emerge.hdf5")), # remap issue 
  n.argwhere(input_list_i== os.path.join(h5_dir, "hlist_0.71730_emerge.hdf5")), # remap issue
  n.argwhere(input_list_i== os.path.join(h5_dir, "hlist_0.93570_emerge.hdf5"))  # remap issue
  ]) )


# creates the redshift list 
redshifts = []
for file_1 in input_list : 
	f1 = h5py.File(file_1,  "r")
	redshifts.append(f1.attrs['redshift'])
	f1.close()

redshifts = n.array(redshifts)

# creates the shell list
Dcom = cosmoMD.comoving_distance(redshifts).value
Dmax = n.hstack((Dcom[0],(Dcom[1:]+Dcom[:-1])/2.))
Dmin = n.hstack(((Dcom[1:]+Dcom[:-1])/2., Dcom[-1]))

def copylc_data(ii, option=False):
	"""
	Creates the selection array to obtain the shell in a snapshot to be added in the light cone
	
	Writes a lightcone shell for each snapshot
	"""
	file_1 = input_list[ii]
	file_out = os.path.join(h5_lc_dir, 'shell_'+os.path.basename( input_list[ii] ) )
	print(file_1, "==>>", file_out)
	f1 = h5py.File(file_1,  "r")
	print( "n halos=",f1['/halo_properties/'].attrs['N_halos'])
	x,y,z=f1[positions_group + '/xyz_Lbox'].value.T*L_box

	distance = ((x-x_obs)**2.+(y-y_obs)**2.+(z-z_obs)**2.)**0.5
	selection = (distance>=Dmin[ii])&(distance<Dmax[ii])
	print( len(distance[selection])," halos in shell ", Dmin[ii], "<d comoving<",Dmax[ii])
	
	if len(distance[selection])>1:
		f = h5py.File(file_out, "a")
		f.attrs['file_name'] = os.path.basename(file_out)
		f.attrs['HDF5_Version']     = h5py.version.hdf5_version
		f.attrs['h5py_version']     = h5py.version.version

		halo_data = f.create_group('halo_position')

		ds = halo_data.create_dataset('x', data = x[selection] )
		ds = halo_data.create_dataset('y', data = y[selection] )
		ds = halo_data.create_dataset('z', data = z[selection] )
		
		ds = halo_data.create_dataset('vx', data = f1['/halo_position/vx'].value[selection] )
		ds = halo_data.create_dataset('vy', data = f1['/halo_position/vy'].value[selection] )
		ds = halo_data.create_dataset('vz', data = f1['/halo_position/vz'].value[selection] )
		
		halo_data = f.create_group('moster_2013_data')
		ds = halo_data.create_dataset('stellar_mass', data = f1['/moster_2013_data/stellar_mass'].value[selection] )
		
		halo_data = f.create_group('remaped_position_L15')
		ds = halo_data.create_dataset('xyz_Lbox', data = f1['/remaped_position_L15/xyz_Lbox'].value[selection] )
		
		halo_data = f.create_group('remaped_position_L3')
		ds = halo_data.create_dataset('xyz_Lbox', data = f1['/remaped_position_L15/xyz_Lbox'].value[selection] )
		
		halo_data = f.create_group('remaped_position_L6')
		ds = halo_data.create_dataset('xyz_Lbox', data = f1['/remaped_position_L15/xyz_Lbox'].value[selection] )
		
		halo_data = f.create_group('cluster_data')
		ds = halo_data.create_dataset('cool_class', data = f1['/cluster_data/cool_class'].value[selection] )
		ds = halo_data.create_dataset('kT', data = f1['/cluster_data/kT'].value[selection] )
		ds = halo_data.create_dataset('log_LX_05_24', data = f1['/cluster_data/log_LX_05_24'].value[selection] )
		ds = halo_data.create_dataset('log_LceX_05_24', data = f1['/cluster_data/log_LceX_05_24'].value[selection] )
		ds = halo_data.create_dataset('log_Mgas', data = f1['/cluster_data/log_Mgas'].value[selection] )
		
		halo_data = f.create_group('halo_properties')
		halo_data.attrs['N_halos'] =  N_halo

		ds = halo_data.create_dataset('scale'                              , data=scale                            [selection])
		ds = halo_data.create_dataset('id'                                 , data=id                               [selection])
		ds = halo_data.create_dataset('desc_scale'                         , data=desc_scale                       [selection])
		ds = halo_data.create_dataset('desc_id'                            , data=desc_id                          [selection])
		ds = halo_data.create_dataset('num_prog'                           , data=num_prog                         [selection])
		ds = halo_data.create_dataset('pid'                                , data=pid                              [selection])
		ds = halo_data.create_dataset('upid'                               , data=upid                             [selection])
		ds = halo_data.create_dataset('desc_pid'                           , data=desc_pid                         [selection])
		ds = halo_data.create_dataset('mvir'                               , data=mvir                             [selection])
		ds = halo_data.create_dataset('rvir'                               , data=rvir                             [selection])
		ds = halo_data.create_dataset('rs'                                 , data=rs                               [selection])
		ds = halo_data.create_dataset('vrms'                               , data=vrms                             [selection])
		ds = halo_data.create_dataset('mmp'                                , data=mmp                              [selection])
		ds = halo_data.create_dataset('scale_of_last_MM'                   , data=scale_of_last_MM                 [selection])
		ds = halo_data.create_dataset('vmax'                               , data=vmax                             [selection])
		ds = halo_data.create_dataset('Jx'                                 , data=Jx                               [selection])
		ds = halo_data.create_dataset('Jy'                                 , data=Jy                               [selection])
		ds = halo_data.create_dataset('Jz'                                 , data=Jz                               [selection])
		ds = halo_data.create_dataset('Spin'                               , data=Spin                             [selection])
		ds = halo_data.create_dataset('Breadth_first_ID'                   , data=Breadth_first_ID                 [selection])
		ds = halo_data.create_dataset('Depth_first_ID'                     , data=Depth_first_ID                   [selection])
		ds = halo_data.create_dataset('Tree_root_ID'                       , data=Tree_root_ID                     [selection])
		ds = halo_data.create_dataset('Orig_halo_ID'                       , data=Orig_halo_ID                     [selection])
		ds = halo_data.create_dataset('Next_coprogenitor_depthfirst_ID'    , data=Next_coprogenitor_depthfirst_ID  [selection])
		ds = halo_data.create_dataset('Last_progenitor_depthfirst_ID'      , data=Last_progenitor_depthfirst_ID    [selection])
		ds = halo_data.create_dataset('Last_mainleaf_depthfirst_ID'        , data=Last_mainleaf_depthfirst_ID      [selection])
		ds = halo_data.create_dataset('Tidal_Force'                        , data=Tidal_Force                      [selection])
		ds = halo_data.create_dataset('Tidal_ID'                           , data=Tidal_ID                         [selection])
		ds = halo_data.create_dataset('Rs_Klypin'                          , data=Rs_Klypin                        [selection])
		ds = halo_data.create_dataset('Mmvir_all'                          , data=Mmvir_all                        [selection])
		ds = halo_data.create_dataset('M200b'                              , data=M200b                            [selection])
		ds = halo_data.create_dataset('M200c'                              , data=M200c                            [selection])
		ds = halo_data.create_dataset('M500c'                              , data=M500c                            [selection])
		ds = halo_data.create_dataset('M2500c'                             , data=M2500c                           [selection])
		ds = halo_data.create_dataset('Xoff'                               , data=Xoff                             [selection])
		ds = halo_data.create_dataset('Voff'                               , data=Voff                             [selection])
		ds = halo_data.create_dataset('Spin_Bullock'                       , data=Spin_Bullock                     [selection])
		ds = halo_data.create_dataset('b_to_a'                             , data=b_to_a                           [selection])
		ds = halo_data.create_dataset('c_to_a'                             , data=c_to_a                           [selection])
		ds = halo_data.create_dataset('Ax'                                 , data=Ax                               [selection])
		ds = halo_data.create_dataset('Ay'                                 , data=Ay                               [selection])
		ds = halo_data.create_dataset('Az'                                 , data=Az                               [selection])
		ds = halo_data.create_dataset('b_to_a_500c'                        , data=b_to_a_500c                      [selection])
		ds = halo_data.create_dataset('c_to_a_500c'                        , data=c_to_a_500c                      [selection])
		ds = halo_data.create_dataset('Ax_500c'                            , data=Ax_500c                          [selection])
		ds = halo_data.create_dataset('Ay_500c'                            , data=Ay_500c                          [selection])
		ds = halo_data.create_dataset('Az_500c'                            , data=Az_500c                          [selection])
		ds = halo_data.create_dataset('TU'                                 , data=TU                               [selection])
		ds = halo_data.create_dataset('M_pe_Behroozi'                      , data=M_pe_Behroozi                    [selection])
		ds = halo_data.create_dataset('M_pe_Diemer'                        , data=M_pe_Diemer                      [selection])
		ds = halo_data.create_dataset('Macc'                               , data=Macc                             [selection])
		ds = halo_data.create_dataset('Mpeak'                              , data=Mpeak                            [selection])
		ds = halo_data.create_dataset('Vacc'                               , data=Vacc                             [selection])
		ds = halo_data.create_dataset('Vpeak'                              , data=Vpeak                            [selection])
		ds = halo_data.create_dataset('Halfmass_Scale'                     , data=Halfmass_Scale                   [selection])
		ds = halo_data.create_dataset('Acc_Rate_Inst'                      , data=Acc_Rate_Inst                    [selection])
		ds = halo_data.create_dataset('Acc_Rate_100Myr'                    , data=Acc_Rate_100Myr                  [selection])
		ds = halo_data.create_dataset('Acc_Rate_1Tdyn'                     , data=Acc_Rate_1Tdyn                   [selection])
		ds = halo_data.create_dataset('Acc_Rate_2Tdyn'                     , data=Acc_Rate_2Tdyn                   [selection])
		ds = halo_data.create_dataset('Acc_Rate_Mpeak'                     , data=Acc_Rate_Mpeak                   [selection])
		ds = halo_data.create_dataset('Mpeak_Scale'                        , data=Mpeak_Scale                      [selection])
		ds = halo_data.create_dataset('Acc_Scale'                          , data=Acc_Scale                        [selection])
		ds = halo_data.create_dataset('First_Acc_Scale'                    , data=First_Acc_Scale                  [selection])
		ds = halo_data.create_dataset('First_Acc_Mvir'                     , data=First_Acc_Mvir                   [selection])
		ds = halo_data.create_dataset('First_Acc_Vmax'                     , data=First_Acc_Vmax                   [selection])
		ds = halo_data.create_dataset('VmaxAtMpeak'                        , data=VmaxAtMpeak                      [selection])
		ds = halo_data.create_dataset('Tidal_Force_Tdyn'                   , data=Tidal_Force_Tdyn                 [selection])
		ds = halo_data.create_dataset('logVmaxVmaxmaxTdynTmpeak'           , data=logVmaxVmaxmaxTdynTmpeak         [selection])
		ds = halo_data.create_dataset('Time_to_future_merger'              , data=Time_to_future_merger            [selection])
		ds = halo_data.create_dataset('Future_merger_MMP_ID'               , data=Future_merger_MMP_ID             [selection])

		f.close()

	f1.close()


copylc_data(ii)
