
"""
.. class:: MultiDark

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class MultiDark is a wrapper to handle Multidark simulations results / outputs.

"""
#from os.path import join
import os
import h5py
import numpy as n
import os
import time

#import glob
#import time

#import cPickle
#import fileinput
#import astropy.io.fits as fits

#import numpy as n
#from scipy.interpolate import interp1d
#import scipy.spatial.ckdtree as t

#from astropy.cosmology import FlatLambdaCDM
#import astropy.units as u
#cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)
#cosmoDS = FlatLambdaCDM(H0=68.46*u.km/u.s/u.Mpc, Om0=0.298734, Ob0=0.046961)
#import astropy.constants as constants

#G =  constants.G.to(u.kpc**3/(u.solMass * u.yr**2)).value

#t_dynamical = lambda rvir, mvir : (rvir**3./(G*mvir))**0.5
#def tau_quenching(tdyn, tau_0, tau_s, m_star):
	#if m_star < 1e10 :
		#return tdyn * tau_0
	#else :
		#return tdyn * tau_0 * (m_star * 10.**(-10.))**(tau_s)

#f_loss = lambda t : 0.05*n.log( 1 + t / (1.4*10**6))


class MultiDarkSimulation :
	"""
	Loads the environement proper to the Multidark simulations. This is the fixed framework of the simulation.
			
	:param Lbox: length of the box in Mpc/h 
	:param wdir: Path to the multidark lightcone directory
	:param boxDir: box directory name
	:param snl: list of snapshots available
	:param zsl: list of redshift corresponding to the snapshots   
	:param zArray: redshift array to be considered to interpolate the redshift -- distance conversion
	:param Hbox: Hubble constant at redshift 0 of the box
	:param Melement: Mass of the resolution element in solar masses.   
	:param columnDict: dictionnary to convert column name into the index to find it in the snapshots
	"""

	def __init__(self, L_box = 1000. ):
		self.L_box = L_box
		self.h = 0.6777
		self.omega_lambda = 0.692885
		self.omega_matter = 0.307115
		self.omega_baryon = 0.048206
		self.ns = 0.96
		self.sigma8 = 0.8228
		self.G = 6.67428 * 10**(-9) # cm3 g-1 s-2
		self.Msun = 1.98892 * 10**(33.) # g
		#self.force_resolution = 5. # kpc /h
		self.Mpart_04Gpc = 9.63 * 10**7 # Msun
		self.Npart_04Gpc = 3840
		self.Mpart_10Gpc = 1.51 * 10**9. # Msun
		self.Npart_10Gpc = 3840
		self.Mpart_25Gpc = 2.359 * 10**10. # Msun
		self.Npart_25Gpc = 3840
		self.Mpart_40Gpc = 9.6 * 10**10. # Msun
		self.Npart_40Gpc = 4096
		self.columnDictHlist = {'scale': 0, 'id': 1, 'desc_scale': 2, 'desc_id': 3, 'num_prog': 4, 'pid': 5, 'upid': 6, 'desc_pid': 7, 'phantom': 8, 'sam_mvir': 9, 'mvir': 10, 'rvir': 11, 'rs': 12, 'vrms': 13, 'mmp?': 14, 'scale_of_last_MM': 15, 'vmax': 16, 'x': 17, 'y': 18, 'z': 19, 'vx': 20, 'vy': 21, 'vz': 22, 'Jx': 23, 'Jy': 24, 'Jz': 25, 'Spin': 26, 'Breadth_first_ID': 27, 'Depth_first_ID': 28, 'Tree_root_ID': 29, 'Orig_halo_ID': 30, 'Snap_num': 31, 'Next_coprogenitor_depthfirst_ID': 32, 'Last_progenitor_depthfirst_ID': 33, 'Last_mainleaf_depthfirst_ID': 34, 'Tidal_Force': 35, 'Tidal_ID': 36, 'Rs_Klypin': 37, 'Mmvir_all': 38, 'M200b': 39, 'M200c': 40, 'M500c': 41, 'M2500c': 42, 'Xoff': 43, 'Voff': 44, 'Spin_Bullock': 45, 'b_to_a': 46, 'c_to_a': 47, 'Ax': 48, 'Ay': 49, 'Az': 50, 'b_to_a_500c' : 51, 'c_to_a_500c' : 52, 'Ax_500c' : 53, 'Ay_500c' : 54, 'Az_500c' : 55, 'TU': 56, 'M_pe_Behroozi': 57, 'M_pe_Diemer': 58, 'Macc': 59, 'Mpeak': 60, 'Vacc': 61, 'Vpeak': 62, 'Halfmass_Scale': 63, 'Acc_Rate_Inst': 64, 'Acc_Rate_100Myr': 65, 'Acc_Rate_1Tdyn': 66, 'Acc_Rate_2Tdyn': 67, 'Acc_Rate_Mpeak': 68, 'Mpeak_Scale': 69, 'Acc_Scale': 70, 'First_Acc_Scale': 71, 'First_Acc_Mvir': 72, 'First_Acc_Vmax': 73, 'VmaxAtMpeak': 74, 'Tidal_Force_Tdyn': 75, 'logVmaxVmaxmaxTdynTmpeak': 76, 'Time_to_future_merger': 77, 'Future_merger_MMP_ID': 78 }
		self.columnDictHlist25 = {'scale': 0, 'id': 1, 'desc_scale': 2, 'desc_id': 3, 'num_prog': 4, 'pid': 5, 'upid': 6, 'desc_pid': 7, 'phantom': 8, 'sam_mvir': 9, 'mvir': 10, 'rvir': 11, 'rs': 12, 'vrms': 13, 'mmp?': 14, 'scale_of_last_MM': 15, 'vmax': 16, 'x': 17, 'y': 18, 'z': 19, 'vx': 20, 'vy': 21, 'vz': 22, 'Jx': 23, 'Jy': 24, 'Jz': 25, 'Spin': 26, 'Breadth_first_ID': 27, 'Depth_first_ID': 28, 'Tree_root_ID': 29, 'Orig_halo_ID': 30, 'Snap_num': 31, 'Next_coprogenitor_depthfirst_ID': 32, 'Last_progenitor_depthfirst_ID': 33, 'Last_mainleaf_depthfirst_ID': 34, 'Rs_Klypin': 35, 'Mmvir_all': 36, 'M200b': 37, 'M200c': 38, 'M500c': 39, 'M2500c': 40, 'Xoff': 41, 'Voff': 42, 'Spin_Bullock': 43, 'b_to_a': 44, 'c_to_a': 45, 'Ax': 46, 'Ay': 47, 'Az': 48, 'b_to_a_500c' : 49, 'c_to_a_500c' : 50, 'Ax_500c' : 51, 'Ay_500c' : 52, 'Az_500c' : 53, 'TU': 54, 'M_pe_Behroozi': 55, 'M_pe_Diemer': 56, 'Halfmass_Radius': 57, 'Macc': 58, 'Mpeak': 59, 'Vacc': 60, 'Vpeak': 61, 'Halfmass_Scale': 62, 'Acc_Rate_Inst': 63, 'Acc_Rate_100Myr': 64, 'Acc_Rate_1Tdyn': 65, 'Acc_Rate_2Tdyn': 66, 'Acc_Rate_Mpeak': 67, 'Mpeak_Scale': 68, 'Acc_Scale': 69, 'First_Acc_Scale': 70, 'First_Acc_Mvir': 71, 'First_Acc_Vmax': 72, 'VmaxAtMpeak': 73, 'Tidal_Force_Tdyn': 74, 'logVmaxVmaxmaxTdynTmpeak': 75, 'Time_to_future_merger': 76, 'Future_merger_MMP_ID': 77 }
		
	def transform_rockstar_hlist_catalog_into_emerge_input_catalog(self, path_2_input, path_2_output, mmin=1.51e11, option = 'ultra-light'):
		"""
		Extracts the columns of interest out of a snapshot of the Multidark simulation.        
		:param path_2_input: path to hlist input catalog
		:param path_2_output: path where the output is written
		:param mmin: minimum mass needed to be selected from the input and written in the output
		:param option: 'complete' copy of all columns or 'ultra-light' copy of 12 columns.
		#module load anaconda
		"""
		if option == 'complete':
			#self.header='scale id desc_scale desc_id': 3, 'num_prog': 4, 'pid': 5, 'upid': 6, 'desc_pid': 7, 'phantom': 8, 'sam_mvir': 9, 'mvir': 10, 'rvir': 11, 'rs': 12, 'vrms': 13, 'mmp?': 14, 'scale_of_last_MM': 15, 'vmax': 16, 'x': 17, 'y': 18, 'z': 19, 'vx': 20, 'vy': 21, 'vz': 22, 'Jx': 23, 'Jy': 24, 'Jz': 25, 'Spin': 26, 'Breadth_first_ID': 27, 'Depth_first_ID': 28, 'Tree_root_ID': 29, 'Orig_halo_ID': 30, 'Snap_num': 31, 'Next_coprogenitor_depthfirst_ID': 32, 'Last_progenitor_depthfirst_ID': 33, 'Last_mainleaf_depthfirst_ID': 34, 'Tidal_Force': 35, 'Tidal_ID': 36, 'Rs_Klypin': 37, 'Mmvir_all': 38, 'M200b': 39, 'M200c': 40, 'M500c': 41, 'M2500c': 42, 'Xoff': 43, 'Voff': 44, 'Spin_Bullock': 45, 'b_to_a': 46, 'c_to_a': 47, 'Ax': 48, 'Ay': 49, 'Az': 50, 'b_to_a_500c' : 51, 'c_to_a_500c' : 52, 'Ax_500c' : 53, 'Ay_500c' : 54, 'Az_500c' : 55, 'TU': 56, 'M_pe_Behroozi': 57, 'M_pe_Diemer': 58, 'Macc': 59, 'Mpeak': 60, 'Vacc': 61, 'Vpeak': 62, 'Halfmass_Scale': 63, 'Acc_Rate_Inst': 64, 'Acc_Rate_100Myr': 65, 'Acc_Rate_1Tdyn': 66, 'Acc_Rate_2Tdyn': 67, 'Acc_Rate_Mpeak': 68, 'Mpeak_Scale': 69, 'Acc_Scale': 70, 'First_Acc_Scale': 71, 'First_Acc_Mvir': 72, 'First_Acc_Vmax': 73, 'VmaxAtMpeak': 74, 'Tidal_Force_Tdyn': 75, 'logVmaxVmaxmaxTdynTmpeak': 76, 'Time_to_future_merger': 77, 'Future_merger_MMP_ID': 78 
			gawk_command = """gawk 'NR>63 {if ( $11 >= """ +str(mmin)+ """ ) print $1 - $79}' """ + path_2_input +" > " + path_2_output
			print("command to be executed: ",gawk_command)
			os.system(gawk_command)
		
		if option == 'ultra-light':
			self.header = """# id desc_scale desc_id pid mvir rvir rs Mpeak Mpeak_scale Acc_Rate_1Tdyn Time_to_future_merger Future_merger_MMP_ID"""
			print(self.header)
			gawk_command = """gawk 'NR>63 {if ( $11 >= """ +str(mmin)+ """ ) print $2, $3, $4, $6, $11, $12, $13, $61, $70, $67, $78, $79 }' """ + path_2_input +" > " + path_2_output
			print("command to be executed: ",gawk_command)
			os.system(gawk_command)
	
	def convert_to_emerge_input_catalog_to_h5_format(self, snap_name, aexp, redshift, age_yr, rho_crit, delta_vir):
		"""
		Converts these files to h5 format
		Reads into a numpy array
		:param path_2_input: path to hlist input catalog
		:param path_2_output: path where the output is written
		:param mmin: minimum mass needed to be selected from the input and written in the output
		:param option: 'complete' copy of all columns or 'ultra-light' copy of 12 columns.
		
		#module load anaconda
		
		or
		
		#module load python33/python/3.3    
		#module load python33/scipy/2015.10
		#module load python33/h5py/2.2      
		"""
		timestr = time.strftime("%Y%m%d-%H%M%S")
		
		path_2_snap = "/u/joco/data/MD/MD_1.0Gpc/emerge/hlist_" + snap_name + ".data"
		path_2_h5_file = "/u/joco/data/MD/MD_1.0Gpc/h5/hlist_"  + snap_name + "_emerge.hdf5"
		if os.path.isfile(path_2_h5_file):
			os.system("rm "+path_2_h5_file)
		
		dtype1 = n.dtype([
		('id'                     , 'i4' )
		,('desc_scale'             , 'f4' )
		,('desc_id'                , 'i4' )
		,('pid'                    , 'i4' )
		,('mvir'                   , 'f4' )
		,('rvir'                   , 'f4' )
		,('rs'                     , 'f4' )
		,('Mpeak'                  , 'f4' )
		,('Mpeak_scale'            , 'f4' )
		,('Acc_Rate_1Tdyn'         , 'f4' )
		,('Time_to_future_merger'  , 'f4' )
		,('Future_merger_MMP_ID'   , 'i4' )
		])

		id, desc_scale, desc_id, pid, mvir, rvir, rs, Mpeak, Mpeak_scale, Acc_Rate_1Tdyn, Time_to_future_merger, Future_merger_MMP_ID =  n.loadtxt(path_2_snap, unpack=True, dtype = dtype1)
		N_halo = len(id)

		# create the h5 container to host the data

		f = h5py.File(path_2_h5_file, "a")
		f.attrs['file_name'] = os.path.basename(path_2_h5_file)
		f.attrs['file_time'] = timestr
		f.attrs['creator']   = 'JC'
		f.attrs['HDF5_Version']     = h5py.version.hdf5_version
		f.attrs['h5py_version']     = h5py.version.version

		f.attrs['aexp']      = float(aexp)#[snap_id]
		f.attrs['redshift']  = float(redshift)#[snap_id]
		f.attrs['age_yr']    = float(age_yr)#[snap_id] 
		f.attrs['rho_crit']  = float(rho_crit)#[snap_id] 
		f.attrs['delta_vir'] = float(delta_vir)#[snap_id]

		halo_data = f.create_group('halo_data')
		halo_data.attrs['N_halos'] =  N_halo

		ds = halo_data.create_dataset('id'                   , data = id                       )
		ds.attrs['units'] = '-'
		ds.attrs['long_name'] = 'halo identifier' 

		ds = halo_data.create_dataset('desc_scale'           , data = desc_scale               )
		ds.attrs['units'] = '-'
		ds.attrs['long_name'] = 'scale factor of the descendant halo' 

		ds = halo_data.create_dataset('desc_id'              , data = desc_id                  )
		ds.attrs['units'] = '-'
		ds.attrs['long_name'] = 'identifier of the descendant halo in the snapshot at desc_scale' 

		ds = halo_data.create_dataset('pid'                  , data = pid                      )
		ds.attrs['units'] = '-'
		ds.attrs['long_name'] = 'parent identifier, -1 if distinct halo' 

		ds = halo_data.create_dataset('mvir'                 , data = mvir                     )
		ds.attrs['units'] = r'$h^{-1} M_\odot$'
		ds.attrs['long_name'] = r'$M_{vir}$' 

		ds = halo_data.create_dataset('rvir'                 , data = rvir                     )
		ds.attrs['units'] = r'$h^{-1} kpc$'
		ds.attrs['long_name'] = r'$r_{vir}$' 

		ds = halo_data.create_dataset('rs'                   , data = rs                       )
		ds.attrs['units'] = r'$h^{-1} kpc$'
		ds.attrs['long_name'] = r'$r_{s}$' 

		ds = halo_data.create_dataset('Mpeak'                , data = Mpeak                    )
		ds.attrs['units'] = r'$h^{-1} M_\odot$'
		ds.attrs['long_name'] = r'$M_{peak}$' 

		ds = halo_data.create_dataset('Mpeak_scale'          , data = Mpeak_scale              )
		ds.attrs['units'] = '-'
		ds.attrs['long_name'] = 'scale factor where M peak is reached' 

		ds = halo_data.create_dataset('Acc_Rate_1Tdyn'       , data = Acc_Rate_1Tdyn           )

		ds = halo_data.create_dataset('Time_to_future_merger', data = Time_to_future_merger    )

		ds = halo_data.create_dataset('Future_merger_MMP_ID' , data = Future_merger_MMP_ID     )

		f.close()

