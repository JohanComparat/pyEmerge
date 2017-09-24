#import h5py    # HDF5 support
import os
import glob
import numpy as n
import sys
import time
from multiprocessing import Pool
from remap import *
C6 = Cuboid(u1=(5, 3, 1), u2=(1, 1, 0), u3=(0, 1, 0))
C3 = Cuboid(u1=(2, 1, 1), u2=(1, 1, 0), u3=(0, 1, 0))

def f(DATA):
	x,y,z = DATA
	return n.transpose([C6.Transform(aa,bb,cc) for aa,bb,cc in zip(x, y, z)])

def f(aa,bb,cc):
	return C6.Transform(aa,bb,cc)
    
L_box = 1000.
dx=0.001
x, y, z = n.arange(1,1000.,dx)/L_box, n.arange(1,1000.,dx)/L_box, n.arange(1,1000.,dx)/L_box

if __name__ == '__main__':
	t0=time.time()
	p = Pool(12)
	out = p.starmap(f, n.transpose([x, y, z]))
	print( time.time()-t0 )

