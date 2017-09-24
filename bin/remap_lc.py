"""
Script to remap coordinates into a 3 and 6 Gpc cuboid

here is how to get the transofrmation coefficients :

data = n.loadtxt("/data17s/darksim/software/cuboidremap-1.0/genremap/list7.txt", unpack=True, dtype='str')
lx = data[0].astype('float')
ly = data[1].astype('float')
lz = data[2].astype('float')
sel = (lx>5.9)&(lx<6.2)&(ly<0.5)&(lz<0.5)
data.T[sel]
sel = (lx>2.2)&(lx<3.5)&(ly<0.8)&(lz<0.8)
data.T[sel]
L1 L2 L3   u11 u12 u13   u21 u22 u23   u31 u32 u33   (periodicity)
'5.9161', '0.4140', '0.4082', '5', '3', '1', '1', '1', '0', '0', '1', '0', '(1)'
'2.4495', '0.7071', '0.5774', '2', '1', '1', '1', '1', '0', '0', '1', '0', '(1)'

writes in the h5 files


"""
import h5py    # HDF5 support
import os
import glob
import numpy as n
import sys

from multiprocessing import Pool
from remap import *
C6 = Cuboid(u1=(5, 3, 1), u2=(1, 1, 0), u3=(0, 1, 0))
C3 = Cuboid(u1=(2, 1, 1), u2=(1, 1, 0), u3=(0, 1, 0))

def f(x,y,z,L_box=1000.):
	return n.transpose([C6.Transform(aa,bb,cc) for aa,bb,cc in zip(x/L_box, y/L_box, z/L_box)])*L_box
    





ii = int(sys.argv[1])

h5_dir = os.path.join(os.environ['MD04'], 'h5' )
L_box = 400.

input_list = n.array(glob.glob(os.path.join(h5_dir, "hlist_?.?????_emerge.hdf5")))
input_list.sort()

file_1 = input_list[ii]
f1 = h5py.File(file_1,  "r+")

print "n halos=",f1['/halo_properties/'].attrs['N_halos']

if f1['/halo_properties/'].attrs['N_halos'] > 0:
	#print f1['/halo_position/x'].value, f1['/halo_position/y'].value, f1['/halo_position/z'].value
	x,y,z = n.transpose([C6.Transform(aa,bb,cc) for aa,bb,cc in zip(f1['/halo_position/x'].value/L_box, f1['/halo_position/y'].value/L_box, f1['/halo_position/z'].value/L_box)])*L_box
	#print x,y,z
	halo_data = f1.create_group('remaped_position_L6')
	ds = halo_data.create_dataset('x', data = x )
	ds.attrs['units'] = 'Mpc/h'
	ds.attrs['long_name'] = 'x' 
	ds = halo_data.create_dataset('y', data = y )
	ds.attrs['units'] = 'Mpc/h'
	ds.attrs['long_name'] = 'y' 
	ds = halo_data.create_dataset('z', data = z )
	ds.attrs['units'] = 'Mpc/h'
	ds.attrs['long_name'] = 'z' 

	x,y,z = n.transpose([C3.Transform(aa,bb,cc) for aa,bb,cc in zip(f1['/halo_position/x'].value/L_box, f1['/halo_position/y'].value/L_box, f1['/halo_position/z'].value/L_box)])*L_box
	#print x,y,z
	halo_data = f1.create_group('remaped_position_L3')
	ds = halo_data.create_dataset('x', data = x )
	ds.attrs['units'] = 'Mpc/h'
	ds.attrs['long_name'] = 'x' 
	ds = halo_data.create_dataset('y', data = y )
	ds.attrs['units'] = 'Mpc/h'
	ds.attrs['long_name'] = 'y' 
	ds = halo_data.create_dataset('z', data = z )
	ds.attrs['units'] = 'Mpc/h'
	ds.attrs['long_name'] = 'z' 

	f1.close()

else:
	f1.close()

if __name__ == '__main__':
    p = Pool(5)
    print(p.map(f, [1, 2, 3]))
