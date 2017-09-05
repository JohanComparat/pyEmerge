
import h5py    # HDF5 support
import os
import glob
import numpy as n
import StellarMass as sm
model = sm.StellarMass()

f_loss = lambda t : 0.05*n.log( 1 + t / (1.4*10**6))

h5_dir = os.path.join(os.environ['HOME'], 'MD10', 'h5' )

input_list = n.array(glob.glob(os.path.join(h5_dir, "hlist_?.?????_emerge.hdf5")))
input_list.sort()

ii = 0 
file_0 = input_list[ii]


f0 = h5py.File(file_0,  "r+")
for item in f0.attrs.keys():
    print( item + ":", f0.attrs[item])

# step 1 initiate on f0 arrays comparing with time=0.


# procedure for new halos, first snapshot :
# evaluate equation (4)
mvir_dot = f0['/halo_data/mvir'].value / f0.attrs['age_yr']
rvir_dot = f0['/halo_data/rvir'].value / f0.attrs['age_yr']
c = f0['/halo_data/rvir'].value / f0['/halo_data/rs'].value
rho_nfw = f0['/halo_data/mvir'].value / (f0['/halo_data/rs'].value**3. * 4. * n.pi * c * (1+c)**2. * (n.log(1.+c)-c/(1.+c)))
# result
dMdt = mvir_dot - 4.*n.pi*f0['/halo_data/rvir'].value *f0['/halo_data/rvir'].value * rvir_dot * rho_nfw
# evaluate equation (1)
dmdt_star = model.f_b * dMdt * model.epsilon(f0['/halo_data/mvir'].value, f0.attrs['redshift']*n.ones_like(f0['/halo_data/mvir'].value))
# evaluate accretion: 0 in this first step
dmdt_star_accretion = n.zeros_like(dmdt_star)
# evaluate equation (11)
f_lost = f_loss(f0.attrs['age_yr']) # equation (12)
# evaluate stellar mass 
stellar_mass = (dmdt_star * (1. - f_lost) + dmdt_star_accretion )* f0.attrs['age_yr'] 
# intra-cluster mass is currently 0
m_icm = n.zeros_like(stellar_mass)

# save results

emerge_data = f0.create_group('emerge_data')

emerge_data.attrs['f_lost'] = f_lost

ds = emerge_data.create_dataset('emerge_data/mvir_dot', data = mvir_dot )
ds.attrs['units'] = r'$h^{-1} M_\odot / yr$'
ds.attrs['long_name'] = r'$d M_{vir} / dt$' 

ds = emerge_data.create_dataset('rvir_dot', data = rvir_dot )
ds.attrs['units'] = r'$h^{-1} kpc / yr$'
ds.attrs['long_name'] = r'$d r_{vir} / dt$' 

ds = emerge_data.create_dataset('dMdt', data = dMdt )
ds.attrs['units'] = r'$h^{-1} M_\odot / yr$'
ds.attrs['long_name'] = r'$\langle d M / dt \rangle$ (4)' 

ds = emerge_data.create_dataset('dmdt_star', data = dmdt_star )
ds.attrs['units'] = r'$h^{-1} M_\odot / yr$'
ds.attrs['long_name'] = r'$ d m_* / dt $ (1)' 

ds = emerge_data.create_dataset('dmdt_star_accretion', data = dmdt_star_accretion )
ds.attrs['units'] = r'$h^{-1} M_\odot / yr$'
ds.attrs['long_name'] = r'$ d m_{acc} / dt $ ' 

ds = emerge_data.create_dataset('stellar_mass', data = stellar_mass )
ds.attrs['units'] = r'$h^{-1} M_\odot $'
ds.attrs['long_name'] = r'$ m_* $ (11)' 

ds = emerge_data.create_dataset('m_icm', data = m_icm )
ds.attrs['units'] = r'$h^{-1} M_\odot $'
ds.attrs['long_name'] = r'$ m_{ICM}$ ' 


f0.close()

    
