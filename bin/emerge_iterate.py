# procedure for existing halos, after the first snapshot was initiated :
import sys
# only input parameter is the numbering of the snapshot. These have to be processed in sequence, cannot be done in parallel ... First 1, 2, 3,  ...
ii = int(sys.argv[1])
#ii = 10 
import time
t0 = time.time()

import h5py    
import os
import glob
import numpy as n
import EmergeStellarMass as sm
model = sm.StellarMass()
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)
import astropy.constants as constants

# generic functions
# =================

f_loss = lambda t : 0.05*n.log( 1 + t / (1.4*10**6))
t_dynamical = lambda rvir, mvir : (rvir**3./(9.797465327217671e-24*mvir))**0.5

def tau_quenching( m_star, tdyn, tau_0=4.282, tau_s=0.363):
	out = n.zeros_like(m_star)
	case_1 = (m_star < 1e10 )
	out[case_1] = tdyn[case_1] * tau_0
	case_2 = (case_1==False)
	out[case_2] = tdyn[case_2] * tau_0 * (m_star[case_2] * 10.**(-10.))**(tau_s)
	return out

# opens the files  
# ================

#h5_dir = os.path.join(os.environ['HOME'], 'MD10', 'h5' )
h5_dir = os.path.join(os.environ['MD10'], 'h5' )

input_list = n.array(glob.glob(os.path.join(h5_dir, "hlist_?.?????_emerge.hdf5")))
input_list.sort()

file_0 = input_list[ii-1]
file_1 = input_list[ii]

f0 = h5py.File(file_0,  "r")
for item in f0.attrs.keys():
    print( item + ":", f0.attrs[item])

f0_scale = os.path.basename(file_0).split('_')[1]

f1 = h5py.File(file_1,  "r+")
for item in f1.attrs.keys():
    print( item + ":", f1.attrs[item])

f1_scale = os.path.basename(file_1).split('_')[1]

# id mapping for halos present in the previous snapshot
# =====================================================

f0_desc_id_unique_list_all_descendents = n.unique(f0['/halo_properties/desc_id'].value)

f1_id_unique_list_descendents_detected_at_next_scale = n.intersect1d(f0_desc_id_unique_list_all_descendents, f1['/halo_properties/id'].value)

mask_f0_to_propagate = n.in1d(f0['/halo_properties/desc_id'].value, f1_id_unique_list_descendents_detected_at_next_scale)

mask_f0_lost = (mask_f0_to_propagate == False )

# evolving halos are given after applying this boolean mask to a f1 quantity :
mask_f1_evolved_from_previous = n.in1d( f1['/halo_properties/id'].value,  f1_id_unique_list_descendents_detected_at_next_scale )

# new halos are given after applying this boolean mask to a f1 quantity :
# new halos in f1, not present in f0
mask_f1_new_halos = (mask_f1_evolved_from_previous==False)
print('new halos', len(mask_f1_new_halos.nonzero()[0]))

# halos descending :
# mask_f0_to_propagate
# mask_f1_evolved_from_previous

s = pd.Series(f0['/halo_properties/desc_id'].value[mask_f0_to_propagate])
f1_id_with_multiple_progenitors = s[s.duplicated()].get_values()
# also = f0_desc_id merging into 1 halo in f1
# merging systems [many halos in f0 into a single f1 halo]
mask_f1_in_a_merging = n.in1d( f1['/halo_properties/id'].value,  f1_id_with_multiple_progenitors )
mask_f0_in_a_merging = n.in1d( f0['/halo_properties/desc_id'].value,  f1_id_with_multiple_progenitors )

# halos mapped 1::1 between snapshots	
mask_f0_evolving_11_halos = ( mask_f0_to_propagate ) & ( mask_f0_in_a_merging == False )
mask_f1_evolving_11_halos = ( mask_f1_evolved_from_previous ) & ( mask_f1_in_a_merging == False )

print('11 mapping', len(mask_f0_evolving_11_halos.nonzero()[0]), len(mask_f1_evolving_11_halos.nonzero()[0]))

print('merging systems', len(f1_id_with_multiple_progenitors))

#for dede in f1_id_with_multiple_progenitors :
	#sel = f0['/halo_properties/desc_id'].value == dede
	#print( dede )
	#print('desc id', f0['/halo_properties/desc_id'].value[sel])
	#print('id', f0['/halo_properties/id'].value[sel])
	#print('pid', f0['/halo_properties/pid'].value[sel])
	#print('mvir', f0['/halo_properties/mvir'].value[sel])
	#print('futur merger mmpid', f0['/halo_properties/Future_merger_MMP_ID'].value[sel])
	#print('time to future merger', f0['/halo_properties/Time_to_future_merger'].value[sel])
	#print('Ms', f0['/emerge_data/stellar_mass'].value[sel])
	#print('SFR',f0['/emerge_data/star_formation_rate'].value[sel])
	#print('mCIM',f0['/emerge_data/m_icm'].value[sel])
	#print('=================================')


# quantities computed for every halos 
# ===================================

# initializing outputs to 0

mvir_dot, rvir_dot, dMdt, dmdt_star, dmdt_star_accretion, f_lost, stellar_mass, star_formation_rate, m_icm = n.zeros_like(f1['/halo_properties/mvir'].value), n.zeros_like(f1['/halo_properties/mvir'].value), n.zeros_like(f1['/halo_properties/mvir'].value), n.zeros_like(f1['/halo_properties/mvir'].value), n.zeros_like(f1['/halo_properties/mvir'].value), 0., n.zeros_like(f1['/halo_properties/mvir'].value), n.zeros_like(f1['/halo_properties/mvir'].value), n.zeros_like(f1['/halo_properties/mvir'].value)

t_dynamical = t_dynamical( f1['/halo_properties/rvir'].value, f1['/halo_properties/mvir'].value )
# in years

# new halos are initiated :
def compute_qtys_new_halos(f_snap, selection):
	"""
	Creates a new galaxy along with the new halo.
	Integrates since the start of the Universe.
	"""
	# evaluate equation (4)
	mvir_dot = f_snap['/halo_properties/mvir'].value[selection] / f_snap.attrs['age_yr']
	rvir_dot = f_snap['/halo_properties/rvir'].value[selection] / f_snap.attrs['age_yr']
	c = f_snap['/halo_properties/rvir'].value[selection] / f_snap['/halo_properties/rs'].value[selection]
	rho_nfw = f_snap['/halo_properties/mvir'].value[selection] / (f_snap['/halo_properties/rs'].value[selection]**3. * 4. * n.pi * c * (1+c)**2. * (n.log(1.+c)-c/(1.+c)))
	# result: dMdt to intialize the halos (not corrected from pseudo evolution)
	#pseudo_evolution_correction = 4.*n.pi*f_snap['/halo_properties/rvir'].value[selection] *f_snap['/halo_properties/rvir'].value[selection] * rvir_dot * rho_nfw
	dMdt = mvir_dot #- pseudo_evolution_correction
	#dMdt[pseudo_evolution_correction > mvir_dot] = mvir_dot[pseudo_evolution_correction > mvir_dot] 
	# evaluate equation (1)
	dmdt_star = model.f_b * dMdt * model.epsilon(f_snap['/halo_properties/mvir'].value[selection], f_snap.attrs['redshift'] * n.ones_like(f_snap['/halo_properties/mvir'].value[selection]))
	# evaluate accretion: 0 in this first step
	dmdt_star_accretion = n.zeros_like(dmdt_star)
	# evaluate equation (11)
	f_lost = f_loss(f_snap.attrs['age_yr']) # equation (12)
	# evaluate stellar mass 
	star_formation_rate = dmdt_star * (1. - f_lost) + dmdt_star_accretion 
	stellar_mass = star_formation_rate * f_snap.attrs['age_yr'] 
	# intra-cluster mass is currently 0
	m_icm = n.zeros_like(stellar_mass)
	return mvir_dot, rvir_dot, dMdt, dmdt_star, dmdt_star_accretion, f_lost, stellar_mass, star_formation_rate, m_icm

def compute_qtys_evolving_halos(f_snap, f_previous, selection, selection_previous):
	# computing dMdt for the halo
	mvir_dot = (f_snap['/halo_properties/mvir'].value[selection]-f_previous['/halo_properties/mvir'].value[selection_previous]) / (f_snap.attrs['age_yr'] - f_previous.attrs['age_yr'])
	rvir_dot = (f_snap['/halo_properties/rvir'].value[selection]-f_previous['/halo_properties/rvir'].value[selection_previous]) / (f_snap.attrs['age_yr'] - f_previous.attrs['age_yr'])
	c = f_snap['/halo_properties/rvir'].value[selection] / f_snap['/halo_properties/rs'].value[selection]
	rho_nfw = f_snap['/halo_properties/mvir'].value[selection] / (f_snap['/halo_properties/rs'].value[selection]**3. * 4. * n.pi * c * (1+c)**2. * (n.log(1.+c)-c/(1.+c)))
	pseudo_evolution_correction = 4.*n.pi*f_snap['/halo_properties/rvir'].value[selection] *f_snap['/halo_properties/rvir'].value[selection] * rvir_dot * rho_nfw
	dMdt = mvir_dot - pseudo_evolution_correction
	
	# initialize the ICM mass to the previous value
	m_icm = f_previous['/emerge_data/m_icm'].value[selection_previous]
	
	# Direct estimates of stellar mass and SFR
	dmdt_star = model.f_b * dMdt * model.epsilon(f_snap['/halo_properties/mvir'].value[selection], f_snap.attrs['redshift'] * n.ones_like(f_snap['/halo_properties/mvir'].value[selection]))
	# evaluate accretion: 0 in this first step
	dmdt_star_accretion = n.zeros_like(dmdt_star)
	# evaluate equation (11)
	f_lost = f_loss(f_snap.attrs['age_yr']-f_previous.attrs['age_yr'])
	# evaluate stellar mass 
	star_formation_rate = dmdt_star * (1. - f_lost) + dmdt_star_accretion 
	stellar_mass = star_formation_rate * (f_snap.attrs['age_yr']-f_previous.attrs['age_yr']) + f_previous['/emerge_data/stellar_mass'].value[selection_previous]
	
	# Variations due to stripping, merging and quenching
	

	# quenching
	quenching = (f_snap['/halo_properties/mvir'].value[selection] < f_snap['/halo_properties/Mpeak'].value[selection]) & (f_snap['/halo_properties/Mpeak_scale'].value[selection] < float(f1_scale))
	t_quench = tau_quenching( f_previous['/emerge_data/stellar_mass'].value[selection_previous], t_dynamical[selection] )
	t_mpeak = cosmoMD.age(1./f_snap['/halo_properties/Mpeak_scale'].value[selection]-1.).to(u.yr).value
	
	# case 1. mdot = mdot at tpeak
	quench_1 = (quenching) & (f_snap.attrs['age_yr'] >= t_mpeak ) & ( f_snap.attrs['age_yr'] <= t_mpeak + t_quench) 
	if len(quench_1.nonzero()[0])>0:
		print("quenching1")
		star_formation_rate[quench_1] = n.ones_like(star_formation_rate[quench_1])*f_previous['/emerge_data/stellar_mass'].value[selection_previous][quench_1]
		stellar_mass[quench_1] = star_formation_rate[quench_1] * (f_snap.attrs['age_yr']-f_previous.attrs['age_yr']) + f_previous['/emerge_data/stellar_mass'].value[selection_previous][quench_1]
	

	# case 2. m dot =0
	quench_2 = (quenching) &(f_snap.attrs['age_yr'] >= t_mpeak + t_quench )
	if len(quench_2.nonzero()[0])>0:
		print("quenching2")
		star_formation_rate[quench_2] = n.zeros_like(star_formation_rate[quench_2])
		stellar_mass[quench_2] = f_previous['/emerge_data/stellar_mass'].value[selection_previous][quench_2]
		
	
	# stripping, case 1
	# negative growth value dMdt => 0
	stripping_1 = (dMdt < 0) 
	if len(stripping_1.nonzero()[0])>0:
		print("stripping1")
		m_icm[stripping_1] += f_previous['/emerge_data/stellar_mass'].value[selection_previous][stripping_1]
		stellar_mass[stripping_1] = n.zeros_like(stellar_mass[stripping_1])
		star_formation_rate[stripping_1] = n.zeros_like(star_formation_rate[stripping_1])
	
	# stripping, case 2
	# after reaching its peak mass,
	# if M < 0.122 * Mpeak, all mass goes to ICM, m=0, mdot=0
	stripping_2 = (f_snap['/halo_properties/mvir'].value[selection] < 0.122*f_snap['/halo_properties/Mpeak'].value[selection]) & (f_snap['/halo_properties/Mpeak_scale'].value[selection] < float(f1_scale))
	if len(stripping_2.nonzero()[0])>0:
		print("stripping2")
		m_icm[stripping_2] += f_previous['/emerge_data/stellar_mass'].value[selection_previous][stripping_2]
		stellar_mass[stripping_2] = n.zeros_like(stellar_mass[stripping_2])
		star_formation_rate[stripping_2] = n.zeros_like(star_formation_rate[stripping_2])

	return mvir_dot, rvir_dot, dMdt, dmdt_star, dmdt_star_accretion, f_lost, stellar_mass, star_formation_rate, m_icm

def merging_single_system(f_snap, f_previous, selection, selection_previous, selection_previous_merging):
	mvir_dot = (f_snap['/halo_properties/mvir'].value[selection]-f_previous['/halo_properties/mvir'].value[selection_previous]) / (f_snap.attrs['age_yr'] - f_previous.attrs['age_yr'])
	rvir_dot = (f_snap['/halo_properties/rvir'].value[selection]-f_previous['/halo_properties/rvir'].value[selection_previous]) / (f_snap.attrs['age_yr'] - f_previous.attrs['age_yr'])
	c = f_snap['/halo_properties/rvir'].value[selection] / f_snap['/halo_properties/rs'].value[selection]
	rho_nfw = f_snap['/halo_properties/mvir'].value[selection] / (f_snap['/halo_properties/rs'].value[selection]**3. * 4. * n.pi * c * (1+c)**2. * (n.log(1.+c)-c/(1.+c)))
	pseudo_evolution_correction = 4.*n.pi*f_snap['/halo_properties/rvir'].value[selection] *f_snap['/halo_properties/rvir'].value[selection] * rvir_dot * rho_nfw
	dMdt = mvir_dot - pseudo_evolution_correction
	
	# initialize the ICM mass to the previous value
	m_icm = f_previous['/emerge_data/m_icm'].value[selection_previous]
	
	# Direct estimates of stellar mass and SFR
	dmdt_star = model.f_b * dMdt * model.epsilon(f_snap['/halo_properties/mvir'].value[selection], f_snap.attrs['redshift'] * n.ones_like(f_snap['/halo_properties/mvir'].value[selection]))
	# evaluate accretion: 0 in this first step
	dmdt_star_accretion = n.zeros_like(dmdt_star)
	# evaluate equation (11)
	f_lost = f_loss(f_snap.attrs['age_yr']-f_previous.attrs['age_yr'])
	# evaluate stellar mass 
	star_formation_rate = dmdt_star * (1. - f_lost) + dmdt_star_accretion 
	stellar_mass = star_formation_rate * (f_snap.attrs['age_yr']-f_previous.attrs['age_yr']) + f_previous['/emerge_data/stellar_mass'].value[selection_previous]
	
	# merging
	# m_star_sat x f_esc => m_host_ICM
	# m_star_sat x (1-f_esc) => m_star_host
	# f_esc = 0.388
	
	#Time_to_future_merger: Time (in Gyr) until the given halo merges into a larger halo.  (-1 if no future merger happens)
	#Future_merger_MMP_ID: most-massive progenitor of the halo into which the given halo merges. (-1 if the main progenitor of the future merger halo does not exist at the given scale factor.)

	stellar_mass += (1.-0.388)*n.sum(f_previous['/emerge_data/stellar_mass'].value[selection_previous_merging])
	m_icm += 0.388*n.sum(f_previous['/emerge_data/stellar_mass'].value[selection_previous_merging])

	return mvir_dot, rvir_dot, dMdt, dmdt_star, dmdt_star_accretion, stellar_mass, star_formation_rate, m_icm

def compute_qtys_merging_halos(f_snap, f_previous, selection, selection_previous):
	# computing dMdt for the halo
	merger_ids = f_snap['/halo_properties/id'].value[selection]
	# merger_ids
	DATA = n.zeros((len(merger_ids),8))
	for jj, merger_id in enumerate(merger_ids):
		mask_f1_host = (selection) & (f_snap['/halo_properties/id'].value == merger_id)
		mask_f0_all = (selection_previous) & (f_previous['/halo_properties/desc_id'].value == merger_id)
		f0_host_id = f_previous['/halo_properties/id'].value[mask_f0_all][n.in1d(f_previous['/halo_properties/id'].value[mask_f0_all], n.unique(f_previous['/halo_properties/Future_merger_MMP_ID'].value[mask_f0_all]))][0]
		print "host id", f0_host_id
		mask_f0_host = (mask_f0_all) & (f_previous['/halo_properties/id'].value == f0_host_id)
		mask_f0_merging = (mask_f0_all) & (f_previous['/halo_properties/id'].value != f0_host_id)
		print merger_id, len(mask_f1_host.nonzero()[0]), len(mask_f0_host.nonzero()[0]), len(mask_f0_merging.nonzero()[0])
		out = n.hstack((merging_single_system(f_snap, f_previous, mask_f1_host, mask_f0_host, mask_f0_merging)))
		#print out, out.shape
		DATA[jj]=out
		#print DATA[jj]
	return DATA.T

if len(mask_f1_new_halos.nonzero()[0]) > 0 :
	mvir_dot[mask_f1_new_halos], rvir_dot[mask_f1_new_halos], dMdt[mask_f1_new_halos], dmdt_star[mask_f1_new_halos], dmdt_star_accretion[mask_f1_new_halos], f_lost, stellar_mass[mask_f1_new_halos], star_formation_rate[mask_f1_new_halos], m_icm[mask_f1_new_halos] = compute_qtys_new_halos(f1, mask_f1_new_halos)

if len(mask_f0_evolving_11_halos.nonzero()[0]) > 0 :
	mvir_dot[mask_f1_evolving_11_halos], rvir_dot[mask_f1_evolving_11_halos], dMdt[mask_f1_evolving_11_halos], dmdt_star[mask_f1_evolving_11_halos], dmdt_star_accretion[mask_f1_evolving_11_halos], f_lost, stellar_mass[mask_f1_evolving_11_halos], star_formation_rate[mask_f1_evolving_11_halos], m_icm[mask_f1_evolving_11_halos] = compute_qtys_evolving_halos(f1, f0, mask_f1_evolving_11_halos, mask_f0_evolving_11_halos)

if  len(f1_id_with_multiple_progenitors) > 0 :
	mvir_dot[mask_f1_in_a_merging], rvir_dot[mask_f1_in_a_merging], dMdt[mask_f1_in_a_merging], dmdt_star[mask_f1_in_a_merging], dmdt_star_accretion[mask_f1_in_a_merging], stellar_mass[mask_f1_in_a_merging], star_formation_rate[mask_f1_in_a_merging], m_icm[mask_f1_in_a_merging] = compute_qtys_merging_halos(f1, f0, mask_f1_in_a_merging, mask_f0_in_a_merging )

emerge_data = f1.create_group('emerge_data')

emerge_data.attrs['f_lost'] = f_lost

ds = emerge_data.create_dataset('mvir_dot', data = mvir_dot )
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

ds = emerge_data.create_dataset('star_formation_rate', data = star_formation_rate )
ds.attrs['units'] = r'$h^{-1} M_\odot / yr$'
ds.attrs['long_name'] = r'$ d m / dt $ ' 

ds = emerge_data.create_dataset('stellar_mass', data = stellar_mass )
ds.attrs['units'] = r'$h^{-1} M_\odot $'
ds.attrs['long_name'] = r'$ m_* $ (11)' 

ds = emerge_data.create_dataset('m_icm', data = m_icm )
ds.attrs['units'] = r'$h^{-1} M_\odot $'
ds.attrs['long_name'] = r'$ m_{ICM}$ ' 

f0.close()
f1.close()

print("time needed",time.time()-t0,"seconds")