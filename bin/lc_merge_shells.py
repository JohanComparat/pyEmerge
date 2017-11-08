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
	x_obs, y_obs, z_obs = 0., 0.7071/2.*L_box, 0.5774/2.*L_box

if remap=="L6":
        positions_group_name = 'remaped_position_L6'
        x_obs, y_obs, z_obs = 0., 0.4140/2.*L_box, 0.4082/2.*L_box


import h5py    # HDF5 support
import os
import glob
import numpy as n

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)

h5_lc_dir = os.path.join(os.environ[env], 'h5_lc', 'shells_'+positions_group_name )
input_list = n.array(glob.glob(os.path.join(h5_lc_dir, "shell_hlist_?.?????_emerge.hdf5")))
input_list.sort()

file_out = os.path.join(os.environ[env], 'h5_lc', 'lc_'+positions_group_name+'_.hdf5')

print("creates redshift list")
# creates the redshift list 
aexp = []
for file_1 in input_list : 
	aexp.append(float(os.path.basename(file_1).split('_')[-2]))

aexp = n.array(aexp)

# loop over the shell files and concatenates everything ? + shell redshift.

def get_data(ii):
	f = h5py.File(input_list[ii],  "r")
	DATA = n.array([n.ones_like( f['/halo_position/x'].value)*(1./aexp[ii]-1.), f['/halo_position/x'].value, f['/halo_position/y'].value, f['/halo_position/z'].value, f['/halo_position/vx'].value, f['/halo_position/vy'].value, f['/halo_position/vz'].value, f['/halo_properties/id'].value, f['/halo_properties/pid'].value, f['/halo_properties/mvir'].value, f['/halo_properties/rvir'].value, f['/halo_properties/rs'].value, f['/halo_properties/Vmax'].value, f['/halo_properties/Mpeak'].value, f['/moster_2013_data/stellar_mass'].value, f['/agn_properties/log_lambda_sar'].value,f['/agn_properties/agn_activity'].value, f['/emerge_data/dMdt'].value, f['/emerge_data/mvir_dot'].value, f['/emerge_data/rvir_dot'].value ])
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
ds.attrs['units'] = 'Mpc'
ds.attrs['long_name'] = 'x' 
ds = halo_data.create_dataset('y', data = DATA[2] )
ds.attrs['units'] = 'Mpc'
ds.attrs['long_name'] = 'y' 
ds = halo_data.create_dataset('z', data = DATA[3] )
ds.attrs['units'] = 'Mpc'
ds.attrs['long_name'] = 'z' 

ds = halo_data.create_dataset('vx', data = DATA[4] )
ds.attrs['units'] = 'km/s'
ds.attrs['long_name'] = 'vx' 
ds = halo_data.create_dataset('vy', data = DATA[5] )
ds.attrs['units'] = 'km/s'
ds.attrs['long_name'] = 'vy' 
ds = halo_data.create_dataset('vz', data = DATA[6] )
ds.attrs['units'] = 'km/s'
ds.attrs['long_name'] = 'vz' 

halo_data = f.create_group('halo_properties')

ds = halo_data.create_dataset('id'                   , data = DATA[7] )
ds.attrs['units'] = '-'
ds.attrs['long_name'] = 'halo identifier' 

ds = halo_data.create_dataset('pid'                  , data = DATA[8] )
ds.attrs['units'] = '-'
ds.attrs['long_name'] = 'parent identifier, -1 if distinct halo' 

ds = halo_data.create_dataset('mvir'                 , data = DATA[9] )
ds.attrs['units'] = r'$h^{-1} M_\odot$'
ds.attrs['long_name'] = r'$M_{vir}$' 

ds = halo_data.create_dataset('rvir'                 , data = DATA[10] )
ds.attrs['units'] = r'$h^{-1} kpc$'                          
ds.attrs['long_name'] = r'$r_{vir}$'                         
                                                             
ds = halo_data.create_dataset('rs'                   , data = DATA[11] )
ds.attrs['units'] = r'$h^{-1} kpc$'                          
ds.attrs['long_name'] = r'$r_{s}$'                           
                                                             
ds = halo_data.create_dataset('Vmax'                , data = DATA[12] )
ds.attrs['units'] = 'km/s'                                   
ds.attrs['long_name'] = r'$V_{max}$'                         
                                                             
ds = halo_data.create_dataset('Mpeak'                , data = DATA[13] )
ds.attrs['units'] = r'$h^{-1} M_\odot$'
ds.attrs['long_name'] = r'$M_{peak}$' 

agn_data = f.create_group('agn_properties')

ds= agn_data.create_dataset('log_lambda_sar'         , data =  DATA[15] )
ds= agn_data.create_dataset('agn_activity'         , data =  DATA[16] )

emerge_data = f.create_group('moster_2013_data')

ds = emerge_data.create_dataset('stellar_mass'          , data = DATA[14]  )
ds.attrs['units'] = r'$ M_\odot$'
ds.attrs['long_name'] = 'stellar mass' 

emerge_data = f.create_group('emerge_data')

ds = emerge_data.create_dataset('dMdt'          , data = DATA[17]  )
ds = emerge_data.create_dataset('mvir_dot'          , data = DATA[18]  )
ds = emerge_data.create_dataset('rvir_dot'          , data = DATA[19]  )

f.close()
