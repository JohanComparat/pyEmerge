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


import h5py    # HDF5 support
import os
import glob
import numpy as n

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)

h5_lc_dir = os.path.join(os.environ[env], 'h5_lc', 'shells_'+positions_group_name )
if os.path.isdir(h5_lc_dir)==False:
	os.mkdir(h5_lc_dir)

h5_dir = os.path.join(os.environ[env], 'h5' )

input_list_i = n.array(glob.glob(os.path.join(h5_dir, "hlist_?.?????_emerge.hdf5")))
input_list_i.sort()

# removing snapshots that cannote be remapped ...
input_list = n.delete(input_list_i,n.array([
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
	x,y,z=f1[positions_group + '/xyx_Lbox'].value.T*L_box

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
		ds.attrs['units'] = 'Mpc/h'
		ds.attrs['long_name'] = 'x' 
		ds = halo_data.create_dataset('y', data = y[selection] )
		ds.attrs['units'] = 'Mpc/h'
		ds.attrs['long_name'] = 'y' 
		ds = halo_data.create_dataset('z', data = z[selection] )
		ds.attrs['units'] = 'Mpc/h'
		ds.attrs['long_name'] = 'z' 

		ds = halo_data.create_dataset('vx', data = f1['/halo_position/vx'].value[selection] )
		ds.attrs['units'] = 'km/s'
		ds.attrs['long_name'] = 'vx' 
		ds = halo_data.create_dataset('vy', data = f1['/halo_position/vy'].value[selection] )
		ds.attrs['units'] = 'km/s'
		ds.attrs['long_name'] = 'vy' 
		ds = halo_data.create_dataset('vz', data = f1['/halo_position/vz'].value[selection] )
		ds.attrs['units'] = 'km/s'
		ds.attrs['long_name'] = 'vz' 

		halo_data = f.create_group('halo_properties')
		
		ds = halo_data.create_dataset('id', data = f1['/halo_properties/id'].value[selection] )
		ds.attrs['units'] = '-'
		ds.attrs['long_name'] = 'halo identifier' 

		ds = halo_data.create_dataset('pid', data = f1['/halo_properties/pid'].value[selection] )
		ds.attrs['units'] = '-'
		ds.attrs['long_name'] = 'parent identifier, -1 if distinct halo' 

		ds = halo_data.create_dataset('mvir', data = f1['/halo_properties/mvir'].value[selection] )
		ds.attrs['units'] = r'$h^{-1} M_\odot$'
		ds.attrs['long_name'] = r'$M_{vir}$' 

		ds = halo_data.create_dataset('rvir', data = f1['/halo_properties/rvir'].value[selection] )
		ds.attrs['units'] = r'$h^{-1} kpc$'
		ds.attrs['long_name'] = r'$r_{vir}$' 

		ds = halo_data.create_dataset('rs', data = f1['/halo_properties/rs'].value[selection] )
		ds.attrs['units'] = r'$h^{-1} kpc$'
		ds.attrs['long_name'] = r'$r_{s}$' 

		ds = halo_data.create_dataset('Vmax' , data = f1['/halo_properties/Vmax'].value[selection])
		ds.attrs['units'] = 'km/s'
		ds.attrs['long_name'] = r'$V_{max}$' 

		ds = halo_data.create_dataset('Mpeak' , data = f1['/halo_properties/Mpeak'].value[selection])
		ds.attrs['units'] = r'$h^{-1} M_\odot$'
		ds.attrs['long_name'] = r'$M_{peak}$' 

		moster_2013_data = f.create_group('moster_2013_data')
		
		ds = moster_2013_data.create_dataset('stellar_mass', data = f1['/moster_2013_data/stellar_mass'].value[selection])
		ds.attrs['units'] = r'$ M_\odot$'
		ds.attrs['long_name'] = 'stellar mass' 

		agn_properties = f.create_group('agn_properties')
		
		ds = agn_properties.create_dataset('log_lambda_sar', data = f1['/agn_properties/log_lambda_sar'].value[selection])
		ds.attrs['units'] = r'log lambda SAR'
		ds.attrs['long_name'] = 'log lambda SAR' 
		
		ds = agn_properties.create_dataset('agn_activity', data = f1['/agn_properties/agn_activity'].value[selection])

		emerge_data = f.create_group('emerge_data')
		
		ds = emerge_data.create_dataset('dMdt', data = f1['/emerge_data/dMdt'].value[selection])
		ds.attrs['units'] = r'$ M_\odot/yr$'
		ds.attrs['long_name'] = 'halo growth rate' 

		ds = emerge_data.create_dataset('mvir_dot', data = f1['/emerge_data/mvir_dot'].value[selection] )
		ds.attrs['units'] = r'$ M_\odot/yr$'
		ds.attrs['long_name'] = 'mvir variation with respect to last snapshot' 

		ds = emerge_data.create_dataset('rvir_dot', data = f1['/emerge_data/rvir_dot'].value[selection]                     )
		ds.attrs['units'] = r'$ kpc /yr $'
		ds.attrs['long_name'] = 'rvir variation with respect to last snapshot' 

		f.close()
	f1.close()


copylc_data(ii)
