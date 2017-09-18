import h5py    # HDF5 support
import os
import glob
import numpy as n

#data = n.loadtxt("/data17s/darksim/software/cuboidremap-1.0/genremap/list7.txt", unpack=True, dtype='str')
#lx = data[0].astype('float')
#ly = data[1].astype('float')
#lz = data[2].astype('float')
#sel = (lx>5.9)&(lx<6.2)&(ly<0.5)&(lz<0.5)
#data.T[sel]
#sel = (lx>2.2)&(lx<3.5)&(ly<0.8)&(lz<0.8)
#data.T[sel]
# L1 L2 L3   u11 u12 u13   u21 u22 u23   u31 u32 u33   (periodicity)
# '5.9161', '0.4140', '0.4082', '5', '3', '1', '1', '1', '0', '0', '1', '0', '(1)'
# '2.4495', '0.7071', '0.5774', '2', '1', '1', '1', '1', '0', '0', '1', '0', '(1)'



from remap import *
C6 = Cuboid(u1=(5, 3, 1), u2=(1, 1, 0), u3=(0, 1, 0))
C3 = Cuboid(u1=(2, 1, 1), u2=(1, 1, 0), u3=(0, 1, 0))

# x,y,z
x,y,z = test,test,test
xr, yr, zr = n.zeros_like(x), n.zeros_like(y), n.zeros_like(z)
Ngal = len(x)
ids = n.arange(len(x))

# optimum speed is per 1e6 arrays: 27k/s, 40s for 1e6 points.
import time
t0=time.time()
test = n.arange(0.1,1.00000001,0.0000001)
print len(test)

x,y,z = n.transpose([C6.Transform(aa,bb,cc) for aa,bb,cc in zip(test,test,test)])

dt = time.time()-t0
print dt, len(test)/(dt)