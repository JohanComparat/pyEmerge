import numpy as n

size = 300000000
NN = 3 # or any positive integer
x = n.random.normal(size=(size, NN)) 
x /= n.linalg.norm(x, axis=1)[:, n.newaxis]

theta = n.arctan(x.T[1]/x.T[2])*180/n.pi
phi = n.arctan((x.T[0]**2+x.T[1]**2)**0.5/x.T[2])*180/n.pi

#topdir = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_remaped_position_L3/'
#sel = (theta<6.7529257176359)&(theta>-6.7529257176359)&(phi>-8.269819492449505)&(phi<8.269819492449505)
#n.savetxt(topdir+'random-ra-dec.txt', n.transpose([theta[sel], phi[sel]]))


topdir = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_remaped_position_L6/'
sel = (theta<1.9766516114702513)&(theta>-1.9766516114702513)&(phi>-2.0047373031569915)&(phi<2.0047373031569915)
n.savetxt(topdir+'random-ra-dec.txt', n.transpose([theta[sel], phi[sel]]))


topdir = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_remaped_position_L15/'
sel = (theta<14.323944878104827)&(theta>-14.323944878104827)&(phi>-20.257311381848154)&(phi<20.257311381848154)
n.savetxt(topdir+'random-ra-dec.txt', n.transpose([theta[sel], phi[sel]]))

# L3 characteristics :
# z< 1.0889947373832305 |ra [deg]|< 6.7529257176359 |dec [deg]|< 8.269819492449505
# N points: 8037075 
# 
# L6 characteristics
# z< 6.697087333514605 |ra [deg]|< 1.9766516114702513 |dec [deg]|< 2.0047373031569915
# N points: 3287299
#
# L15 characteristics
# z< 0.5423857379098544 |ra [deg]|< 14.323944878104827 |dec [deg]|< 20.257311381848154
# N points: 8237518
