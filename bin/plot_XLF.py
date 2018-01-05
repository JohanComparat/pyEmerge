import numpy as n
import glob
import h5py
import os
import time
import sys

out_dir = os.path.join(os.path.join("/afs/mpe/www/people/comparat/", "eRoMok", "h5", "LX_function" ))

h5_files = n.array(glob.glob(os.path.join(os.environ['MD10'], "h5", "hlist_?.?????_emerge.hdf5")))
h5_files.sort()

bins = n.arange(38,48,0.25)
xb = (bins[1:] + bins[:-1]) / 2.
hh = 0.6777
volume=1000.**3./hh**3.

zmin, zmax, z_center, Lxmin, Lxmax, Lx_c,   Nobj, phi, phierr = n.loadtxt(os.path.join(os.environ['DARKSIM_DIR'], 'observations', 'LXFunction', 'miyaji_2015.ascii'), unpack=True)


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p


def plot_XLF(h5_file):
  f1 = h5py.File(h5_file,  "r")
  print(h5_file)
  redshift = f1.attrs['redshift']
  p.figure(1, (6,6))
  sel = (redshift>zmin)&(redshift<zmax)
  if len(sel.nonzero()[0])>0:
    p.errorbar(Lx_c[sel], phi[sel], xerr=[Lx_c[sel]-Lxmin[sel], Lxmax[sel]-Lx_c[sel]], yerr=phierr[sel], label='Mi15 '+str(z_center[sel][0]), ls='dashed')
  mass = n.log10(f1['/moster_2013_data/stellar_mass'].value )
  lsar = f1['/agn_properties/log_lambda_sar'].value 
  active = (f1['/agn_properties/agn_activity'].value ==1 )
  sel = (mass>0) & (mass!=n.inf) & (n.isnan(mass)==False) & (lsar>0) & (active)
  print( h5_file, len(mass), len(mass[sel]), len(mass[sel])>0 )
  LX = mass[sel]+lsar[sel]
  counts, bb = n.histogram(LX, bins=bins)
  dN_dVdlogM = counts/(bins[1:]-bins[:-1])/volume/n.log(10)
  ok = (dN_dVdlogM>0)
  p.plot(xb[ok], dN_dVdlogM[ok], label='AGN sim Bo16', lw=2)#, ls='dashed')
  #g10(dN_dVdlogM_g_AGN[ok]), label='AGN simulated')#, lw=2)
  p.xlabel('LX [2-10 keV]')
  p.ylabel('Phi ')
  p.xlim((40., 46.))
  p.yscale('log')
  #p.ylim((-9,-2))
  p.title('z='+str(n.round(redshift,3)))
  p.grid()
  p.legend(loc=0, frameon=False)
  print(f1.attrs['file_name'])
  p.savefig(os.path.join(out_dir, "MD10_"+str(f1.attrs['aexp'])+"_XLF.png"))
  p.clf()
  f1.close()

#plot_SMF(h5_files[50])
#plot_SMF(h5_files[65])
for h5_file in h5_files[::-1]:
  ##try:
  plot_XLF(h5_file)
  ##except( ValueError, KeyError ):
  ##  pass
