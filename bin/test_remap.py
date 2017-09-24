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

def f(x,y,z,L_box=1000.):
	return n.transpose([C6.Transform(aa,bb,cc) for aa,bb,cc in zip(x/L_box, y/L_box, z/L_box)])*L_box
    
    
x, y, z = n.arange(1000.), n.arange(1000.), n.arange(1000.)

if __name__ == '__main__':
	t0=time.time()
	p = Pool(2)
	out = p.starmap(f, zip(x, y, z))
	print( time.time()-t0 )

