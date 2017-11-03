import os
import glob
import numpy as n
from scipy.stats import norm

import EmergeStellarMass as sm
model = sm.StellarMass()

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p

mvir = 10**n.arange(8,15,0.01)
z=0.

meanSM= lambda Mh, z :  2. * 2. * ( 0.0351 - 0.0247 * z/(1.+z)) / ((Mh/ (10**(11.59 + 1.195 * z/(1.+z))) )**(- 1.376 + 0.826 * z/(1.+z)) + ( Mh /(10**(11.59 + 1.195 * z/(1.+z))) )**(0.608 + 0.329 *z/(1.+z)) ) 

meanSM_mix= lambda Mh, z :  2. * ( 0.0351 - 0.0247 * z/(1.+z)) / ((Mh/ (10**(11.79 + 1.195 * z/(1.+z))) )**(- 0.9 + 0.826  * z/(1.+z)) + ( Mh /(10**(11.79 + 1.195 * z/(1.+z))) )**(0.67 + 0.2 * z/(1.+z)) ) 

mean_SMHMR_13 = meanSM(mvir, z)
mean_SMHMR_mix = meanSM_mix(mvir, z)

mean_SMHMR_17 = model.epsilon(mvir, z) * 0.156 

out_dir = os.path.join(os.path.join("/afs/mpe/www/people/comparat/", "eRoMok" ))

p.figure(1, (6,6))
p.plot(mvir, mean_SMHMR_13, label='Mo13')
p.plot(mvir, mean_SMHMR_mix, label='mix')
p.plot(mvir, mean_SMHMR_17, label='Mo17')
p.xlabel('halo mass')
p.ylabel('stellar mass / halo mass')
#p.xlim((9., 12.2))
#p.ylim((-9,-2))
p.xscale('log')
p.yscale('log')
p.title('z='+str(n.round(z,3)))
p.grid()
p.legend(loc=0, frameon=False)
p.savefig(os.path.join(out_dir, "Mo13_Mo17_SMHMR.png"))
p.clf()