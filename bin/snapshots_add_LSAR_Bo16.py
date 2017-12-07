"""
Routine to add Moster et al. 2013 stellar masses

python3 lc_add_Ms_Mo13.py 115 MD10
"""
import time
t0=time.time()

import sys
ii = int(sys.argv[1])
env = sys.argv[2] # 'MD10'
status = sys.argv[3]

import h5py    # HDF5 support
import os
import glob
import numpy as n

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)

#from multiprocessing import Pool
#n_proc=12
#pool = Pool(n_proc)


h5_dir = os.path.join(os.environ[env], 'h5' )
input_list = n.array(glob.glob(os.path.join(h5_dir, "hlist_?.?????_emerge.hdf5")))
input_list.sort()

file_1 = input_list[ii]
print(file_1)
f1 = h5py.File(file_1,  "r+")

z = f1.attrs['redshift']
logM = n.log10(f1['/moster_2013_data/stellar_mass'].value)
agn_random_number = n.random.random(len(logM))
log_lSAR = n.zeros(len(logM))

dl=0.01
log_lambda_SAR_values = n.arange(32-dl,36+2*dl,dl)

#f_z = lambda z : n.piecewise(z, [z <= 1.1, z > 1.1], [ lambda z : (1.+z)**(5.82), lambda z : (1. + 1.1)**(5.82) * ((1.+z)/(1.+1.1))**(2.36)])
## single value
#fz = f_z( z )
#f_Mstar =lambda logM : (10**(logM - 10.99) )**(0.24) * n.e**( - 10**(logM - 10.99) )
##array of values
#fMstar = f_Mstar( logM )
#pre_factor = 10**(- 6.86) * fMstar * fz  

def f_lambda_sar( DATA ):
  logM, log_lambda_SAR = DATA
  log_lambda_SAR_var = 10**( log_lambda_SAR - 33.8 + 0.48 * (logM - 11.) )
  #return 1. / ( log_lambda_SAR_var**(1.01 - 0.58 * (z - 1.1)) + log_lambda_SAR_var**(3.72) )
  return 1. / ( log_lambda_SAR_var**(1.01 - 0.58 * (z - 1.1)) + log_lambda_SAR_var**(2.72) )


t0 = time.time()
ii0=0
ii_step=10000
for ii0 in n.arange(0, len(logM), ii_step):
  ii1=ii0+ii_step
  X,Y = n.meshgrid(logM[ii0:ii1], log_lambda_SAR_values)
  #Z = n.ones_like(X)*z
  probas_un = f_lambda_sar([ X, Y])#, Z ])
  norm = n.sum(probas_un, axis=0)
  probas = probas_un / norm
  cmat = n.array([ agn_random_number[ii0:ii1] > n.sum(probas.T[:,jj:], axis=1) for jj in n.arange(len(log_lambda_SAR_values)) ])
  #print(cmat.shape, cmat[0])
  #print(cmat.T[1])
  values = log_lambda_SAR_values[n.array([n.min(n.where(cmat.T[jj]==True)) for jj in n.arange(len(cmat.T)) ])]
  #print(values.shape, values[:10])
  log_lSAR[ii0:ii1] = values
  print(ii0, len(logM), n.max(log_lSAR), time.time()-t0)

print('max lsar',n.max(log_lSAR))
if status=='create':
  
  halo_data = f1.create_group('agn_properties')
  
  ds = halo_data.create_dataset('log_lambda_sar', data = log_lSAR )
  ds.attrs['units'] = '\log_{10}(\lambda_{SAR}/[])'
  ds.attrs['long_name'] = 'log10 of the specific accretion rate' 
  
  ds = halo_data.create_dataset('random_number_lambda_sar', data = agn_random_number )
  
if status=='update':
  f1['/agn_properties/log_lambda_sar'][:] = log_lSAR
  f1['/agn_properties/random_number_lambda_sar'][:] = agn_random_number

f1.close()


