"""
Routine to add AGN activity

python3 snapshots_add_AGN_activity_Bo16.py 113 MD10 
"""
import time
t0=time.time()

import sys
ii = int(sys.argv[1])
env = sys.argv[2] 
status = 'create'

import h5py    # HDF5 support
import os
import glob
import numpy as n
from scipy.interpolate import interp1d

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)

from multiprocessing import Pool
n_proc=12
pool = Pool(n_proc)


h5_dir = os.path.join(os.environ[env], 'h5' )
input_list = n.array(glob.glob(os.path.join(h5_dir, "hlist_?.?????_emerge.hdf5")))
input_list.sort()

file_1 = input_list[ii]
print(file_1)
f1 = h5py.File(file_1,  "r+")

z = f1.attrs['redshift']
logM = n.log10(f1['/moster_2013_data/stellar_mass'].value)
agn_random_number = n.random.random(len(logM))
activity = n.zeros(len(logM))
 
yi = f1['stellar_mass_function_moster_2013/duty_cycle'].value
xi = 0.5*(f1['stellar_mass_function_moster_2013/stellar_mass_low'].value+f1['stellar_mass_function_moster_2013/stellar_mass_up'].value)

itp = interp1d(xi,yi)

proba = itp(logM)

activity[agn_random_number < proba] = n.ones_like(activity[agn_random_number < proba])

if status=='create':
  f1['/agn_properties'].create_dataset('random_number_activity', data= agn_random_number)
  f1['/agn_properties'].create_dataset('agn_activity', data = activity)
  
if status=='update':
  f1['/agn_properties/activity'][:] = activity
  f1['/agn_properties/random_number_actiity'][:] = agn_random_number

f1.close()


