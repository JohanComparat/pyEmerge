
import h5py    # HDF5 support
import os
import glob
import numpy as n
import StellarMass as sm

f_loss = lambda t : 0.05*n.log( 1 + t / (1.4*10**6))

h5_dir = os.path.join(os.environ['HOME'], 'MD10', 'h5' )

input_list = n.array(glob.glob(os.path.join(h5_dir, "hlist_?.?????_emerge.hdf5")))
input_list.sort()

ii = 0 
file_0 = input_list[ii]
file_1 = input_list[ii+1]


f0 = h5py.File(file_0,  "r+")
for item in f0.attrs.keys():
    print( item + ":", f0.attrs[item])

f1 = h5py.File(file_1,  "r+")
for item in f1.attrs.keys():
    print( item + ":", f1.attrs[item])

    
pix=0
sel_pix =(f1['/halo_data/id'].value == f0['/halo_data/desc_id'].value[pix])

id_map = n.ravel(n.array([n.argwhere(f1['/halo_data/id'].value == pix_id) for pix_id in f0['/halo_data/desc_id'].value]) )

check_mapping_correct = (f1['/halo_data/id'].value[id_map] == f0['/halo_data/desc_id'].value)

#if check_mapping_correct.all(): 

# step 1 initiate on f0 arrays comparing with time=0.

# compute mvir_dot, rvir_dot, f_loss

emerge_data = f0.create_group('emerge_data')

# procedure for new halos, first snapshot :
mvir_dot = f0['/halo_data/mvir'].value / f0.attrs['age_yr']
rvir_dot = f0['/halo_data/rvir'].value / f0.attrs['age_yr']
f_lost = f_loss(f0.attrs['age_yr'])

c = f0['/halo_data/rvir'].value / f0['/halo_data/rs'].value
rho_nfw = f0['/halo_data/mvir'].value / (f0['/halo_data/rs'].value**3. * 4. * n.pi * c * (1+c)**2. * (n.log(1.+c)-c/(1.+c)))

dMdt = mvir_dot - 4.*n.pi*f0['/halo_data/rvir'].value *f0['/halo_data/rvir'].value * rvir_dot * rho_nfw

dmdt_star = sm.f_b * dMdt * sm.epsilon(f0['/halo_data/mvir'].value, f0.attrs['redshift']*n.ones_like(f0['/halo_data/mvir'].value))

# procedure for existing halos, after first snapshot is initiated :
mvir_dot = (f1['/halo_data/mvir'].value[id_map]-f0['/halo_data/mvir'].value) / (f1.attrs['age_yr'] - f0.attrs['age_yr'])
rvir_dot = (f1['/halo_data/rvir'].value[id_map]-f0['/halo_data/rvir'].value) / (f1.attrs['age_yr'] - f0.attrs['age_yr'])

f_lost = f_loss(f1.attrs['age_yr']-f0.attrs['age_yr'])

c = f1['/halo_data/rvir'].value[id_map] / f1['/halo_data/rs'].value[id_map]
rho_nfw = f1['/halo_data/mvir'].value[id_map] / (f1['/halo_data/rs'].value[id_map]**3. * 4. * n.pi * c * (1+c)**2. * (n.log(1.+c)-c/(1.+c)))

dMdt = mvir_dot - 4.*n.pi*f1['/halo_data/rvir'].value[id_map] *f1['/halo_data/rvir'].value[id_map] * rvir_dot * rho_nfw


# sqrt($12*$12*$12/("""+str(G)+"""*$11)), 


emerge_data.attrs['f_lost'] = f_lost

ds = emerge_data.create_dataset('mvir_dot', data = mvir_dot )
ds.attrs['units'] = r'$h^{-1} M_\odot / yr$'
ds.attrs['long_name'] = r'$d M_{vir} / dt$' 

ds = emerge_data.create_dataset('rvir_dot', data = rvir_dot )
ds.attrs['units'] = r'$h^{-1} kpc / yr$'
ds.attrs['long_name'] = r'$d r_{vir} / dt$' 




stellar_data = f0.create_group('stellar_data')
#stellar_data.attrs['N_halos'] =  N_halo

ds = halo_data.create_dataset('id', data = id )
ds.attrs['units'] = '-'
ds.attrs['long_name'] = 'halo identifier' 


f.close()
