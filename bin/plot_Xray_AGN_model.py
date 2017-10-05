import numpy as n
import glob
import h5py
import os
import time
import sys

import XrayLuminosity
xr = XrayLuminosity.XrayLuminosity()

out_dir = os.path.join(os.path.join("/afs/mpe/www/people/comparat/", "eRoMok", "h5", "xray_agn_model" ))


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

all_z = n.arange(0, 6, 0.1)

p.figure(1, (6,6))
p.axes([0.17, 0.17, 0.78, 0.78])
p.plot(all_z, xr.f_z(all_z))
p.xlabel('redshift')
p.ylabel(r'$f(z)$')
#p.xlim((9., 12.2))
#p.ylim((-9,-2))
p.grid()
#p.legend(loc=0, frameon=False)
p.savefig(os.path.join(out_dir, "B016_fz.png"))
p.clf()

few_z = n.arange(0,3.6,0.5)
many_MS = n.arange(6,13,0.1)

p.figure(1, (6,6))
p.axes([0.17, 0.17, 0.78, 0.78])
#for zz in few_z:
p.plot(many_MS, xr.f_Mstar(many_MS, 0))#, label="z="str(0))

p.xlabel(r'$log(M*/M_\odot)$')
p.ylabel(r'$f_*(M)$')
#p.xlim((9., 12.2))
#p.ylim((-9,-2))
p.grid()
#p.legend(loc=0, frameon=False)
p.savefig(os.path.join(out_dir, "B016_fmstar.png"))
p.clf()

few_z = n.array([0.55, 1.15, 2.])
few_MS = n.arange(9,12.1,0.5)
log_lambda_SAR = n.arange(32,36.1,0.1)

for zz in few_z:
  p.figure(1, (6,6))
  p.axes([0.17, 0.17, 0.78, 0.78])
  p.plot(many_MS, n.array([xr.Phi_stellar_mass(MS, zz) for MS in many_MS]), label='all')
  p.plot(many_MS, n.array([xr.Phi_stellar_mass_mmin(MS, zz, 42-MS) for MS in many_MS]), label='>42')
  p.plot(many_MS, n.array([xr.Phi_stellar_mass_mmin(MS, zz, 43-MS) for MS in many_MS]), label='>43')
  p.plot(many_MS, n.array([xr.Phi_stellar_mass_mmin(MS, zz, 43.5-MS) for MS in many_MS]), label='>43.5')
  p.plot(many_MS, n.array([xr.Phi_stellar_mass_mmin(MS, zz, 44-MS) for MS in many_MS]), label='>44')
  p.plot(many_MS, n.array([xr.Phi_stellar_mass_mmin(MS, zz, 44.5-MS) for MS in many_MS]), label='>44.5')
  p.xlabel(r'$log(M*/M_\odot)$')
  p.ylabel(r'$\Phi_*(M)$')
  p.xlim((9., 12.5))
  p.ylim((1e-8,1e-2))
  p.yscale('log')
  p.grid()
  p.title('z='+str(zz))
  p.legend(loc=0, frameon=False)
  p.savefig(os.path.join(out_dir, "B016_psi_star_LX_"+str(zz)+".png"))
  p.clf()


for zz in few_z:
  p.figure(1, (6,6))
  p.axes([0.17, 0.17, 0.78, 0.78])
  for logM in few_MS:
    p.plot(log_lambda_SAR, xr.f_lambda_sar( logM, zz, log_lambda_SAR ), label=str(logM))
  p.xlabel(r'$log(\lambda_{SAR})$')
  p.ylabel(r'$f((\lambda_{SAR}),M)$')
  #p.xlim((9., 12.2))
  #p.ylim((-9,-2))
  p.yscale('log')
  p.title('z='+str(zz))
  p.grid()
  p.legend(loc=0, frameon=False)
  p.savefig(os.path.join(out_dir, "B016_flsar_"+str(zz)+".png"))
  p.clf()

for zz in few_z:
  p.figure(1, (6,6))
  p.axes([0.17, 0.17, 0.78, 0.78])
  for logM in few_MS:
    p.plot(log_lambda_SAR, xr.psi_log(log_lambda_SAR, logM, zz), label=str(logM))
  p.xlabel(r'$log(\lambda_{SAR})$')
  p.ylabel(r'$\Psi((\lambda_{SAR}),M,z)$')
  #p.xlim((9., 12.2))
  p.ylim((1e-10,1e-2))
  p.yscale('log')
  p.title('z='+str(zz))
  p.grid()
  p.legend(loc=0, frameon=False)
  p.savefig(os.path.join(out_dir, "B016_psi_"+str(zz)+".png"))
  p.clf()

for zz in few_z:
  p.figure(1, (6,6))
  p.axes([0.17, 0.17, 0.78, 0.78])
  p.plot(many_MS, n.array([xr.Phi_stellar_mass(MS, zz) for MS in many_MS]), label='>32')
  p.plot(many_MS, n.array([xr.Phi_stellar_mass_33(MS, zz) for MS in many_MS]), label='>33')
  p.plot(many_MS, n.array([xr.Phi_stellar_mass_34(MS, zz) for MS in many_MS]), label='>34')
  p.xlabel(r'$log(M*/M_\odot)$')
  p.ylabel(r'$\Phi_*(M)$')
  p.xlim((9., 12.5))
  p.ylim((1e-8,1e-2))
  p.yscale('log')
  p.grid()
  p.title('z='+str(zz))
  p.legend(loc=0, frameon=False)
  p.savefig(os.path.join(out_dir, "B016_psi_star_"+str(zz)+".png"))
  p.clf()
  
