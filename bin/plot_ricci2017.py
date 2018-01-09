import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p
import numpy as n
from scipy.special import erf
ricci_ct_f = lambda z: 0.22 + 0.18 * z**0.4
fraction_ricci = lambda lsar, z : ricci_ct_f(z)+(0.8-ricci_ct_f(z))*(0.5+0.5*erf((-lsar+32.75)/0.4))

nhs = n.arange(32, 36, 0.1)

p.figure(1, (5,5))
for zz in n.arange(0.,3.,0.5):
	p.plot(nhs, fraction_ricci(nhs, zz), label='z='+str(zz))

p.axhline(0.22)
p.legend(frameon=False)
p.xlabel('lambda SAR')
p.ylabel('fraction')
p.ylim((0.,1.))
p.text(33,0.1,'thick 24<nh<26')
p.text(33,0.9,'unobscured 20<nH<22')
p.text(32,0.5,'thin')
p.text(32,0.4,'22<nH<24')
p.savefig('/home/comparat/Desktop/ricci2017.png')
p.clf()

