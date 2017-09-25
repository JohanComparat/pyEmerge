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
import time
print("start", time.time())
import sys
ii = int(sys.argv[1])
print("snapshot", ii)
import h5py    # HDF5 support
import os
import glob
import numpy as n
from multiprocessing import Pool
# imports the remapping library
from remap import Cuboid
C6 = Cuboid(u1=(5, 3, 1), u2=(1, 1, 0), u3=(0, 1, 0))
C3 = Cuboid(u1=(2, 1, 1), u2=(1, 1, 0), u3=(0, 1, 0))

def f6(aa,bb,cc):
	return C6.Transform(aa,bb,cc)

def f3(aa,bb,cc):
	return C3.Transform(aa,bb,cc)
    
def read_data(ii, L_box = 400., env= 'MD04'):
	"""
	Read all input data and returns 
	 - the h5 file: f1
	 - the coordinates to be mapped: x, y, z
	"""
h5_dir = os.path.join(os.environ[env], 'h5' )
input_list = n.array(glob.glob(os.path.join(h5_dir, "hlist_?.?????_emerge.hdf5")))
input_list.sort()
file_1 = input_list[ii]
f1 = h5py.File(file_1,  "r+")
	print( "n halos=",f1['/halo_properties/'].attrs['N_halos'])
	return f1, f1['/halo_position/x'].value/L_box, f1['/halo_position/y'].value/L_box, f1['/halo_position/z'].value/L_box

def write_mapped_coordinates(f1, out, L_box, group_name = 'remaped_position_L6'):
	"""
	Writes the new coordinates to file
	:param f1: h5 file
	:param x1,y1,z1: new coordinates
	:param group_name: name of the new group containing the new data in the h5 file. Example 'remaped_position_L6'
	"""
	print("writes")
	halo_data = f1.create_group(group_name)
	halo_data.attrs['L_box'] = L_box
	ds = halo_data.create_dataset('xyx_Lbox', data = out )
	ds.attrs['units'] = 'L box'
	ds.attrs['long_name'] = 'x,y,z' 

if __name__ == '__main__':
	p = Pool(12)
	# reads the data
	L_box = 400.
	env= 'MD04'
	f1, x0, y0, z0 = read_data(ii, L_box, env)
	# maps coordinates to L6
	out3 = p.starmap(f3, n.transpose([x0, y0, z0]))
	out6 = p.starmap(f6, n.transpose([x0, y0, z0]))
	# writes the results
	write_mapped_coordinates(f1, out6, L_box,  group_name = 'remaped_position_L6')
	write_mapped_coordinates(f1, out3, L_box,  group_name = 'remaped_position_L3')
	f1.close()
