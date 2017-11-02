import numpy as n
import glob
import h5py
import os
import time
import sys

out_dir = os.path.join(os.path.join("/afs/mpe/www/people/comparat/", "eRoMok", "h5", "stellar_mass_function" ))

h5_files = n.array(glob.glob(os.path.join(os.environ['MD10'], "h5", "hlist_?.?????_emerge.hdf5")))
h5_files.sort()

bins = n.arange(6,13,0.1)
xb = (bins[1:] + bins[:-1]) / 2.


smf_ilbert13 = lambda M, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s : ( phi_1s * (M/M_star) ** alpha_1s + phi_2s * (M/M_star) ** alpha_2s ) * n.e ** (-M/M_star) * (M/ M_star)
ll_dir = os.path.join(os.environ['GIT_NBODY_NPT'], 'data', 'stellar_mass_function')
path_ilbert13_SMF = os.path.join(ll_dir, "ilbert_2013_mass_function_params.txt")
zmin, zmax, N, M_comp, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s, log_rho_s = n.loadtxt(path_ilbert13_SMF, unpack=True)

smf_ilbert_fun = n.array([
lambda mass : smf_ilbert13( mass , 10**M_star[0], phi_1s[0]*10**(-3), alpha_1s[0], phi_2s[0]*10**(-3), alpha_2s[0] )
, lambda mass : smf_ilbert13( mass , 10**M_star[1], phi_1s[1]*10**(-3), alpha_1s[1], phi_2s[1]*10**(-3), alpha_2s[1] )
, lambda mass : smf_ilbert13( mass , 10**M_star[2], phi_1s[2]*10**(-3), alpha_1s[2], phi_2s[2]*10**(-3), alpha_2s[2] )
, lambda mass : smf_ilbert13( mass , 10**M_star[3], phi_1s[3]*10**(-3), alpha_1s[3], phi_2s[3]*10**(-3), alpha_2s[3] )
, lambda mass : smf_ilbert13( mass , 10**M_star[4], phi_1s[4]*10**(-3), alpha_1s[4], phi_2s[4]*10**(-3), alpha_2s[4] )
, lambda mass : smf_ilbert13( mass , 10**M_star[5], phi_1s[5]*10**(-3), alpha_1s[5], phi_2s[5]*10**(-3), alpha_2s[5] )
, lambda mass : smf_ilbert13( mass , 10**M_star[6], phi_1s[6]*10**(-3), alpha_1s[6], phi_2s[6]*10**(-3), alpha_2s[6] )
, lambda mass : smf_ilbert13( mass , 10**M_star[7], phi_1s[7]*10**(-3), alpha_1s[7], phi_2s[7]*10**(-3), alpha_2s[7] )
])

mbins = n.arange(8,12.5,0.25)

smf_ilbert_zmin = n.array([ 
0.2
, 0.5
, 0.8
, 1.1
, 1.5
, 2.0
, 2.5
, 3.0 ])

smf_ilbert_zmax = n.array([ 
0.5
, 0.8
, 1.1
, 1.5
, 2.0
, 2.5
, 3.0
, 4.0 ])

smf_ilbert_name = n.array([ "Il13 "+str(zmin)+"<z<"+str(zmax) for zmin, zmax in zip(smf_ilbert_zmin,smf_ilbert_zmax) ])


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p


def plot_SMF(h5_file):
  f1 = h5py.File(h5_file,  "r")
  print(h5_file)
  redshift = f1.attrs['redshift']
  p.figure(1, (6,6))
  for fun, name in zip(smf_ilbert_fun, smf_ilbert_name):
    p.plot(mbins, n.log10(fun(10**mbins)), label=name, ls='dashed', lw=0.5)

  logMs_low    = f1['stellar_mass_function/stellar_mass_low'].value
  logMs_up     = f1['stellar_mass_function/stellar_mass_up'].value
  counts       = f1['stellar_mass_function/counts'].value
  dN_dVdlogM_g = f1['stellar_mass_function/dN_dVdlogM'].value 
  AGN_HGMF = f1['/stellar_mass_function/AGN_HGMF'].value
  
  ok = (dN_dVdlogM_g>0)
  print( "SMF", n.min(logMs_low[ok]), n.max(logMs_up[ok]) )
  p.plot((logMs_low[ok] + logMs_up[ok])/2.-n.log10(0.6777), n.log10(dN_dVdlogM_g[ok])-2*n.log10(0.6777), label='SIM GAL')#, lw=2)
  p.plot((logMs_low[ok] + logMs_up[ok])/2., n.log10(AGN_HGMF[ok]), label='model HGMF')#, lw=2)
  #print(f1['/agn_model/stellar_mass'].value)
  #p.plot(f1['/agn_model/stellar_mass'].value-n.log10(0.6777), n.log10(f1['/agn_model/HGMF'].value), label='MODEL AGN')
  p.xlabel('stellar mass')
  p.ylabel('log Phi stellar mass')
  p.xlim((9., 12.2))
  p.ylim((-9,-2))
  p.title('z='+str(n.round(redshift,3)))
  p.grid()
  p.legend(loc=0, frameon=False)
  p.savefig(os.path.join(out_dir, "MD10_"+f1.attrs['file_name'].split('_')[1]+"_SMF.png"))
  p.clf()
  f1.close()

for h5_file in h5_files:
  try:
    plot_SMF(h5_file)
  except( ValueError, KeyError ):
    pass
