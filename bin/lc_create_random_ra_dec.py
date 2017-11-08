import numpy as n

topdir = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_remaped_position_L3/'


size = 300000000
NN = 3 # or any positive integer
x = n.random.normal(size=(size, NN)) 
x /= n.linalg.norm(x, axis=1)[:, n.newaxis]

theta = n.arctan(x.T[1]/x.T[2])*180/n.pi
phi = n.arctan((x.T[0]**2+x.T[1]**2)**0.5/x.T[2])*180/n.pi

sel = (theta<6.7529257176359)&(theta>-6.7529257176359)&(phi>-8.269819492449505)&(phi<8.269819492449505)

n.savetxt(topdir+'random-ra-dec.txt', n.transpose([theta[sel], phi[sel]]))


