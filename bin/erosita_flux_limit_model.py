"""
Library of SSM to downsample each subsurveys according to science goals

"""

#import numpy as n
from scipy.special import erf


eRo_FL = {
  'd' : lambda x: 0.5+0.5*erf((x+14.22)/0.2),
  'm' : lambda x: 0.5+0.5*erf((x+14.09)/0.2),
  's' : lambda x: 0.5+0.5*erf((x+14.00)/0.3)
  }

import matplotlib.pyplot as p

log_F = n.arange(-15, -13, 0.01)

p.plot(10**log_F, eRo_FL['d'](log_F), 'r')
p.plot(10**log_F, eRo_FL['m'](log_F), 'g')
p.plot(10**log_F, eRo_FL['s'](log_F), 'b')
p.xscale('log')
p.grid()
p.savefig('test.png')
p.clf()

