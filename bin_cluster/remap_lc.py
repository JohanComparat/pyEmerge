"""
Script to remap coordinates into a 3 and 6 Gpc cuboid

here is how to get the transofrmation coefficients :
import numpy as n
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115)#, Ob0=0.048206)

data = n.loadtxt("/data17s/darksim/software/cuboidremap-1.0/genremap/list7.txt", unpack=True, dtype='str')
lx = data[0].astype('float')
ly = data[1].astype('float')
lz = data[2].astype('float')
sel = (ly>1.085)&(ly<1.1)#&(ly<1.)&(lz<1.)
data.T[sel]
sel = (lx>5.9)&(lx<6.2)&(ly<0.5)&(lz<0.5)
data.T[sel]
sel = (lx>2.2)&(lx<3.5)&(ly<0.8)&(lz<0.8)
data.T[sel]
L1 L2 L3   u11 u12 u13   u21 u22 u23   u31 u32 u33   (periodicity)
#'2.2361', '1.0954', '0.4082', '2', '1', '0', '1', '0', '1', '1', '0', '0', '(1)'
#'5.9161', '0.4140', '0.4082', '5', '3', '1', '1', '1', '0', '0', '1', '0', '(1)'
#'2.4495', '0.7071', '0.5774', '2', '1', '1', '1', '1', '0', '0', '1', '0', '(1)'

lz, along x : '1.4142', '1.0000', '0.7071', '1', '1', '0', '0', '0', '1', '1', '0', '0', '(12)'
mz, along y : '1.4142', '1.2247', '0.5774', '1', '1', '0', '1', '0', '1', '1', '0', '0', '(1)'
hz, along x : '1.7321', '0.8165', '0.7071', '1', '1', '1', '1', '0', '0', '0', '1', '0', '(1)'

writes in the h5 files


"""
import time
print("start", time.time())
import sys
ii = int(sys.argv[1])
env = sys.argv[2]
L_box = float(sys.argv[3])
print("snapshot", ii, env)
import h5py    # HDF5 support
import os
import glob
import numpy as n
from multiprocessing import Pool
# imports the remapping library
from remap import Cuboid
#C6 = Cuboid(u1=(5, 3, 1), u2=(1, 1, 0), u3=(0, 1, 0))
#C3 = Cuboid(u1=(2, 1, 1), u2=(1, 1, 0), u3=(0, 1, 0))
#C15 = Cuboid(u1=(1, 1, 0), u2=(0, 0, 1), u3=(1, 0, 0))
C15 = Cuboid(u1=(1, 1, 0), u2=(0, 0, 1), u3=(1, 0, 0)) 
C3 = Cuboid(u1=(1, 1, 0), u2=(1, 0, 1), u3=(1, 0, 0)) 
C6 = Cuboid(u1=(1, 1, 1), u2=(1, 0, 0), u3=(0, 1, 0))

def f6(aa,bb,cc):
	return C6.Transform(aa,bb,cc)

def f3(aa,bb,cc):
	return C3.Transform(aa,bb,cc)

def f15(aa,bb,cc):
	return C15.Transform(aa,bb,cc)
    
def read_data(ii, L_box = 1000., env= 'MD10'):
	"""
	Read all input data and returns 
	 - the h5 file: f1
	 - the coordinates to be mapped: x, y, z
	"""
	h5_dir = os.path.join(os.environ[env], 'cluster_h5' )
	input_list = n.array(glob.glob(os.path.join(h5_dir, "hlist_?.?????.hdf5")))
	input_list.sort()
	file_1 = input_list[ii]
	print("opens ",file_1)
	f1 = h5py.File(file_1,  "r+")
	print( "n halos=",f1['/halo_properties/'].attrs['N_halos'])
	return f1, f1['/halo_position/x'].value/L_box, f1['/halo_position/y'].value/L_box, f1['/halo_position/z'].value/L_box

def write_mapped_coordinates(f1, out, L_box, group_name = 'remaped_position_L6',status='create'):
	"""
	Writes the new coordinates to file
	:param f1: h5 file
	:param x1,y1,z1: new coordinates
	:param group_name: name of the new group containing the new data in the h5 file. Example 'remaped_position_L6'
	"""
	if status=='create':
		print("writes new group "+group_name)
		halo_data = f1.create_group(group_name)
		halo_data.attrs['L_box'] = L_box
		ds = halo_data.create_dataset('xyz_Lbox', data = out )
		ds.attrs['units'] = 'L box'
		ds.attrs['long_name'] = 'x,y,z' 
	if status=='update':
		print('writes update '+group_name)
		f1['/'+group_name+'/xyz_Lbox'][:] = out


if __name__ == '__main__':
	p = Pool(12)
	# reads the data
	#L_box = 400.
	#env= 'MD04'
	f1, x0, y0, z0 = read_data(ii, L_box, env)
	# map to L3
	out3 = p.starmap(f3, n.transpose([x0, y0, z0]))
	write_mapped_coordinates(f1, out3, L_box,  group_name = 'remaped_position_L3', status='update')
	# map to L6
	out6 = p.starmap(f6, n.transpose([x0, y0, z0]))
	write_mapped_coordinates(f1, out6, L_box,  group_name = 'remaped_position_L6', status='update')
	# map to L15
	out15 = p.starmap(f15, n.transpose([x0, y0, z0]))
	write_mapped_coordinates(f1, out15, L_box,  group_name = 'remaped_position_L15', status='update')
	f1.close()

