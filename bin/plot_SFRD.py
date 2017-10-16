import numpy as n
import glob
import h5py
import os
import time
import sys

out_dir = os.path.join(os.path.join("/afs/mpe/www/people/comparat/", "eRoMok", "h5", "star_formation_rate_density" ))

h5_files = n.array(glob.glob(os.path.join(os.environ['MD10'], "h5", "hlist_?.?????_emerge.hdf5")))
h5_files.sort()


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p


def get_sfrd(h5_file):
  f1 = h5py.File(h5_file,  "r")
  print(h5_file)
  redshift = f1.attrs['redshift']
  sfrd    = f1['star_formation_rate_density/sfrd'].value
  f1.close()
  return redshift, sfrd

sfrds,zs=[],[]
for h5_file in h5_files:
  try:
    z,sfrd=get_sfrd(h5_file)
    sfrds.append(sfrd)
    zs.append(z)
  except( ValueError, KeyError ):
    pass

sfrds = n.array([sfrds])
zs = n.array([zs])

psi = lambda z : 0.015*(1+z)**2.7/(1+((1+z)/2.9)**5.6)

p.figure(1, (6,6))
p.plot(zs, n.log10(sfrds), label='MD10')

p.plot(zs, n.log10(psi(zs)), label='Madau 14')

p.xlabel('redshift')
p.ylabel('SFRD')
p.xlim((0., 10))
p.ylim((-2.5,-0.4))
p.grid()
p.legend(loc=0, frameon=False)
p.savefig(os.path.join(out_dir, "MD10_SFRD.png"))
p.clf()
