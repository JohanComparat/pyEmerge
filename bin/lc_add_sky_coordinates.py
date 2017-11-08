"""
Routine to add sky coordinates

L1 L2 L3   u11 u12 u13   u21 u22 u23   u31 u32 u33   (periodicity)
'5.9161', '0.4140', '0.4082', '5', '3', '1', '1', '1', '0', '0', '1', '0', '(1)'
'2.4495', '0.7071', '0.5774', '2', '1', '1', '1', '1', '0', '0', '1', '0', '(1)'

python3 lc_add_sky_coordinates.py 1.3

"""

import sys
env = 'MD10'# sys.argv[1] # 'MD10'
L_box = 1000./0.6777 # float(sys.argv[2]) / 0.6777

positions_group_name = sys.argv[1] # 'remaped_position_L3'

if positions_group_name=='remaped_position_L3':
  z_max = 1.3
  x_obs, y_obs, z_obs = 0., 0.7071/2.*L_box, 0.5774/2.*L_box
  strech_factor_los = 2.4495
  strech_factor_y = 0.7071
  strech_factor_z = 0.5774


if positions_group_name=='remaped_position_L6':
  z_max = 6.3
  x_obs, y_obs, z_obs = 0., 0.4140/2.*L_box, 0.4082/2.*L_box
  strech_factor_los = 5.9161
  strech_factor_y = 0.4140
  strech_factor_z = 0.4082


import h5py    # HDF5 support
import os
import glob
import numpy as n

from scipy.interpolate import interp1d
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)

print("reads the light cone")
file_lc = os.path.join(os.environ[env], 'h5_lc', 'lc_'+positions_group_name+'_.hdf5')
f = h5py.File(file_lc,  "r+")

a = 1. / ( 1. + f['/halo_position/z_snap'].value ) # scale factor
x = f['/halo_position/x'].value - x_obs
y = f['/halo_position/y'].value - y_obs
z = f['/halo_position/z'].value - z_obs
print('RA, DEC positions')
ra  = n.arctan( z / x ) * 180/n.pi
dec = n.arctan( y / x ) * 180/n.pi

vx = f['/halo_position/vx'].value
vy = f['/halo_position/vy'].value
vz = f['/halo_position/vz'].value

rr=(x**2 + y**2 + z**2)**0.5

print("interpolates z d_comoving")
z_array = n.arange(0, z_max, 0.00001)
dc_to_z = interp1d(cosmoMD.comoving_distance(z_array), z_array)

print("dimensions of the light cone")
z_reach = dc_to_z(strech_factor_los * L_box)
dec_max = (strech_factor_y*L_box*u.Mpc/cosmoMD.kpc_comoving_per_arcmin(z_reach).to(u.Mpc/u.degree)/2.).value
ra_max = (strech_factor_z*L_box*u.Mpc/cosmoMD.kpc_comoving_per_arcmin(z_reach).to(u.Mpc/u.degree)/2.).value
print("z<",z_reach,"|ra [deg]|<",ra_max, "|dec [deg]|<", dec_max)

print("evaluates redshifts")
redshift_R = dc_to_z(rr)
vPara = (vx * x + vy * y + vz * z )/rr 
rr_s = rr + vPara / (a * 67.77)
redshift_S = dc_to_z(rr_s)

selection = (redshift_R < z_reach) & (abs(ra)<ra_max) & (abs(dec)<dec_max)

print("select data in the LC",  len(redshift_R[selection]), " out of ",  len(redshift_R) )

halo_data = f.create_group('sky_position')

ds = halo_data.create_dataset('redshift_R', data = redshift_R )
ds = halo_data.create_dataset('redshift_S', data = redshift_S )
ds = halo_data.create_dataset('RA', data = ra )
ds.attrs['units'] = 'deg'
ds.attrs['long_name'] = 'right ascension' 
ds = halo_data.create_dataset('DEC', data = dec )
ds.attrs['units'] = 'deg'
ds.attrs['long_name'] = 'declination' 
ds = halo_data.create_dataset('selection', data = selection )
ds.attrs['long_name'] = 'True: inside the light cone, False: border effect' 

f.close()
