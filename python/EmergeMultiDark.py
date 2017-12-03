
"""
.. class:: EmergeMultiDark

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class EmergeMultiDark is a wrapper to handle Multidark simulations results / outputs and link them to the EMERGE model.


Imports
-------

module load anaconda

 - import os
 - import h5py
 - import numpy as n
 - import os
 - import time


"""
import os
import h5py
import numpy as n
import os
import time

class MultiDarkSimulation :
	"""
	Loads the environement proper to the Multidark simulations. This is the fixed framework of the simulation.
			
	:param Lbox: length of the box in Mpc/h 
	
	Atributes
	---------
	
	 - h = 0.6777
	 - omega_lambda = 0.692885
	 - omega_matter = 0.307115
	 - omega_baryon = 0.048206
	 - ns = 0.96
	 - sigma8 = 0.8228
	 - G = 6.67428 * 10**(-9) # cm3 g-1 s-2
	 - Msun = 1.98892 * 10**(33.) # g
	 - Mpart_04Gpc = 9.63 * 10**7 # Msun
	 - Npart_04Gpc = 3840
	 - Mpart_10Gpc = 1.51 * 10**9. # Msun
	 - Npart_10Gpc = 3840
	 - Mpart_25Gpc = 2.359 * 10**10. # Msun
	 - Npart_25Gpc = 3840
	 - Mpart_40Gpc = 9.6 * 10**10. # Msun
	 - Npart_40Gpc = 4096
	 - columnDictHlist dictionary of the columns in the hlist files
	 - columnDictHlist25 same as above but for a specific set of files from MD25
	
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
		
	def transform_rockstar_hlist_catalog_into_emerge_input_catalog(self, path_2_input, path_2_output, mmin=1.51e11, option = 'light'):
		"""
		Extracts the columns of interest out of a snapshot of the Multidark simulation.
		Uses a gawk command to rewrite the files
		
		:param path_2_input: path to hlist input catalog
		:param path_2_output: path where the output is written
		:param mmin: minimum mass [mvir] needed to be selected from the input and written in the output
		:param option: 'complete': copy of all 79 columns or 'light': copy of 19 columns or 'ultra-light': copy of 12 columns.
		"""
		if option == 'complete':
			self.header = """# id"""
			gawk_command = """gawk 'NR>63 {if ( $11 >= """ +str(mmin)+ """ ) print $1 - $79}' """ + path_2_input +" > " + path_2_output
			print("command to be executed: ",gawk_command)
			os.system(gawk_command)
		
		if option == 'light':
			self.header = """# id desc_scale desc_id pid mvir rvir rs Mpeak Mpeak_scale Acc_Rate_1Tdyn Time_to_future_merger Future_merger_MMP_ID vmax x y z vx vy vz """
			print(self.header)
			gawk_command = """gawk 'NR>63 {if ( $11 >= """ +str(mmin)+ """ ) print $2, $3, $4, $6, $11, $12, $13, $61, $70, $67, $78, $79, $17, $18, $19, $20, $21, $22, $23 }' """ + path_2_input +" > " + path_2_output
			print("command to be executed: ",gawk_command)
			os.system(gawk_command)
		
		if option == 'ultra-light':
			self.header = """# id desc_scale desc_id pid mvir rvir rs Mpeak Mpeak_scale Acc_Rate_1Tdyn Time_to_future_merger Future_merger_MMP_ID"""
			print(self.header)
			gawk_command = """gawk 'NR>63 {if ( $11 >= """ +str(mmin)+ """ ) print $2, $3, $4, $6, $11, $12, $13, $61, $70, $67, $78, $79 }' """ + path_2_input +" > " + path_2_output
			print("command to be executed: ",gawk_command)
			os.system(gawk_command)

	def convert_to_emerge_input_catalog_to_h5_format_light(self, path_2_snap, path_2_h5_file, aexp, redshift, age_yr, rho_crit, delta_vir):
		"""
		Converts the ascii files create by the function 'transform_rockstar_hlist_catalog_into_emerge_input_catalog' to h5 format.
		
		Parameters
		..........
		
		:param path_2_snap: path to the ascii snapshot 
		:param path_2_h5_file: path to the h5 output file  
		:param aexp: expansion factor
		:param redshift: redshift
		:param age_yr: age of the Universe in years
		:param rho_crit: critical density at this redshift
		:param delta_vir: virial overdensity at this redshift (Bryan and Norman 98 approximation)
		
		Uses the summary file
		.....................
		
		aexp, redshift, age_yr, rho_crit, delta_vir = n.loadtxt("/u/joco/data/MD/MD_1.0Gpc/hlists_MD_1.0Gpc.ascii", unpack=True)
		to provide input parameters
		 - aexp      : expansion factor
		 - redshift  : redshift
		 - age_yr    : age of the Universe in years
		 - rho_crit  : critical density in Msun/Mpc3
		 - delta_vir : virial overdensity (Bryan and Norman 98)
		
		submitting
		..........
		
		pyEmerge/bin/generate_and_submit_convert-2-h5.py
		
		shows how to submit all jobs
		
		data model
		..........
		
		h5 file attributes
		
		 - f.attrs['file_name'] = h5 file name
		 - f.attrs['file_time'] = time stamp
		 - f.attrs['creator']   = 'JC'
		 - f.attrs['HDF5_Version']     = h5py.version.hdf5_version
		 - f.attrs['h5py_version']     = h5py.version.version
		 - f.attrs['aexp']      = float(aexp)#[snap_id]
		 - f.attrs['redshift']  = float(redshift)#[snap_id]
		 - f.attrs['age_yr']    = float(age_yr)#[snap_id] 
		 - f.attrs['rho_crit']  = float(rho_crit)#[snap_id] 
		 - f.attrs['delta_vir'] = float(delta_vir)#[snap_id]

		h5 file data 
		
		 - '/halo_position/' contains 
		   - '/halo_position/x',  position, Mpc/h
		   - '/halo_position/y',  position, Mpc/h
		   - '/halo_position/z',  position, Mpc/h
		   - '/halo_position/vx', velocity, km/s
		   - '/halo_position/vy', velocity, km/s
		   - '/halo_position/vz', velocity, km/s
		 - '/halo_properties/' contains       [attribute, number of halos: 'N_halos']
		   - '/halo_properties/id', halo identifier
		   - '/halo_properties/desc_scale', scale factor of the descendant halo
		   - '/halo_properties/desc_id', identifier of the descendant halo in the snapshot at desc_scale
		   - '/halo_properties/pid', parent identifier, -1 if distinct halo
		   - '/halo_properties/mvir', virial mass, r'$h^{-1} M_\odot$'
		   - '/halo_properties/rvir', virial radius, r$h^{-1} kpc$'
		   - '/halo_properties/rs', scale radius, r$h^{-1} kpc$'
		   - '/halo_properties/Vmax', maximum velocity, km/s
		   - '/halo_properties/Mpeak', r'$M_{peak}$', r'$h^{-1} M_\odot$'
		   - '/halo_properties/Mpeak_scale', scale factor where M peak is reached
		   - '/halo_properties/Acc_Rate_1Tdyn'
		   - '/halo_properties/Time_to_future_merger'
		   - '/halo_properties/Future_merger_MMP_ID'
		   
		"""
		timestr = time.strftime("%Y%m%d-%H%M%S")
		
		#path_2_snap = "/u/joco/data/MD/MD_1.0Gpc/emerge/hlist_" + snap_name + ".data"
		#path_2_h5_file = "/u/joco/data/MD/MD_1.0Gpc/h5/hlist_"  + snap_name + "_emerge.hdf5"
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
		,('Vmax'    	           , 'f4' )
		,('x'    	               , 'f4' )
		,('y'    	               , 'f4' )
		,('z'    	               , 'f4' )
		,('vy'    	               , 'f4' )
		,('vx'    	               , 'f4' )
		,('vz'    	               , 'f4' )
		])

		id, desc_scale, desc_id, pid, mvir, rvir, rs, Mpeak, Mpeak_scale, Acc_Rate_1Tdyn, Time_to_future_merger, Future_merger_MMP_ID, Vmax, x, y, z, vx, vy, vz =  n.loadtxt(path_2_snap, unpack=True, dtype = dtype1)
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

		halo_data = f.create_group('halo_position')
		
		ds = halo_data.create_dataset('x', data = x )
		ds.attrs['units'] = 'Mpc/h'
		ds.attrs['long_name'] = 'x' 
		ds = halo_data.create_dataset('y', data = y )
		ds.attrs['units'] = 'Mpc/h'
		ds.attrs['long_name'] = 'y' 
		ds = halo_data.create_dataset('z', data = z )
		ds.attrs['units'] = 'Mpc/h'
		ds.attrs['long_name'] = 'z' 

		ds = halo_data.create_dataset('vx', data = vx )
		ds.attrs['units'] = 'km/s'
		ds.attrs['long_name'] = 'vx' 
		ds = halo_data.create_dataset('vy', data = vy )
		ds.attrs['units'] = 'km/s'
		ds.attrs['long_name'] = 'vy' 
		ds = halo_data.create_dataset('vz', data = vz )
		ds.attrs['units'] = 'km/s'
		ds.attrs['long_name'] = 'vz' 
		
		halo_data = f.create_group('halo_properties')
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

		ds = halo_data.create_dataset('Vmax'                , data = Vmax                    )
		ds.attrs['units'] = 'km/s'
		ds.attrs['long_name'] = r'$V_{max}$' 
		
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
	
	def convert_to_emerge_input_catalog_to_h5_format_ultralight(self, snap_name, aexp, redshift, age_yr, rho_crit, delta_vir):
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

	def convert_to_h5_format_cluster(self, snap_name, path_2_h5_file, aexp, redshift):

		timestr = time.strftime("%Y%m%d-%H%M%S")
		
		#path_2_snap = "/u/joco/data/MD/MD_1.0Gpc/emerge/hlist_" + snap_name + ".data"
		#path_2_h5_file = "/u/joco/data/MD/MD_1.0Gpc/h5/hlist_"  + snap_name + "_emerge.hdf5"
		if os.path.isfile(path_2_h5_file):
			os.system("rm "+path_2_h5_file)

		dtype1 = n.dtype([
		( 'scale'                             , 'float64' )
		,('id'                                , 'int64' )
		,('desc_scale'                        , 'float64' )
		,('desc_id'                           , 'int64' )
		,('num_prog'                          , 'float64' )
		,('pid'                               , 'int64' )
		,('upid'                              , 'int64' )
		,('desc_pid'                          , 'int64' )
		,('phantom'                           , 'float64' )
		,('sam_mvir'                          , 'float64' )
		,('mvir'                              , 'float64' )
		,('rvir'                              , 'float64' )
		,('rs'                                , 'float64' )
		,('vrms'                              , 'float64' )
		,('mmp'                              , 'float64' )
		,('scale_of_last_MM'                  , 'float64' )
		,('vmax'                              , 'float64' )
		,('x'                                 , 'float64' )
		,('y'                                 , 'float64' )
		,('z'                                 , 'float64' )
		,('vx'                                , 'float64' )
		,('vy'                                , 'float64' )
		,('vz'                                , 'float64' )
		,('Jx'                                , 'float64' )
		,('Jy'                                , 'float64' )
		,('Jz'                                , 'float64' )
		,('Spin'                              , 'float64' )
		,('Breadth_first_ID'                  , 'int64' )
		,('Depth_first_ID'                    , 'int64' )
		,('Tree_root_ID'                      , 'int64' )
		,('Orig_halo_ID'                      , 'int64' )
		,('Snap_num'                          , 'int64' )
		,('Next_coprogenitor_depthfirst_ID'   , 'int64' )
		,('Last_progenitor_depthfirst_ID'     , 'int64' )
		,('Last_mainleaf_depthfirst_ID'       , 'int64' )
		,('Tidal_Force'                       , 'float64' )
		,('Tidal_ID'                          , 'int64' )
		,('Rs_Klypin'                         , 'float64' )
		,('Mmvir_all'                         , 'float64' )
		,('M200b'                             , 'float64' )
		,('M200c'                             , 'float64' )
		,('M500c'                             , 'float64' )
		,('M2500c'                            , 'float64' )
		,('Xoff'                              , 'float64' )
		,('Voff'                              , 'float64' )
		,('Spin_Bullock'                      , 'float64' )
		,('b_to_a'                            , 'float64' )
		,('c_to_a'                            , 'float64' )
		,('Ax'                                , 'float64' )
		,('Ay'                                , 'float64' )
		,('Az'                                , 'float64' )
		,('b_to_a_500c'                       , 'float64' )
		,('c_to_a_500c'                       , 'float64' )
		,('Ax_500c'                           , 'float64' )
		,('Ay_500c'                           , 'float64' )
		,('Az_500c'                           , 'float64' )
		,('TU'                                , 'float64' )
		,('M_pe_Behroozi'                     , 'float64' )
		,('M_pe_Diemer'                       , 'float64' )
		,('Macc'                              , 'float64' )
		,('Mpeak'                             , 'float64' )
		,('Vacc'                              , 'float64' )
		,('Vpeak'                             , 'float64' )
		,('Halfmass_Scale'                    , 'float64' )
		,('Acc_Rate_Inst'                     , 'float64' )
		,('Acc_Rate_100Myr'                   , 'float64' )
		,('Acc_Rate_1Tdyn'                    , 'float64' )
		,('Acc_Rate_2Tdyn'                    , 'float64' )
		,('Acc_Rate_Mpeak'                    , 'float64' )
		,('Mpeak_Scale'                       , 'float64' )
		,('Acc_Scale'                         , 'float64' )
		,('First_Acc_Scale'                   , 'float64' )
		,('First_Acc_Mvir'                    , 'float64' )
		,('First_Acc_Vmax'                    , 'float64' )
		,('VmaxAtMpeak'                       , 'float64' )
		,('Tidal_Force_Tdyn'                  , 'float64' )
		,('logVmaxVmaxmaxTdynTmpeak'          , 'float64' )
		,('Time_to_future_merger'             , 'float64' )
		,('Future_merger_MMP_ID'              , 'int64' )
		])
		
		scale                             ,
		id                                ,
		desc_scale                        ,
		desc_id                           ,
		num_prog                          ,
		pid                               ,
		upid                              ,
		desc_pid                          ,
		phantom                           ,
		sam_mvir                          ,
		mvir                              ,
		rvir                              ,
		rs                                ,
		vrms                              ,
		mmp                              ,
		scale_of_last_MM                  ,
		vmax                              ,
		x                                 ,
		y                                 ,
		z                                 ,
		vx                                ,
		vy                                ,
		vz                                ,
		Jx                                ,
		Jy                                ,
		Jz                                ,
		Spin                              ,
		Breadth_first_ID                  ,
		Depth_first_ID                    ,
		Tree_root_ID                      ,
		Orig_halo_ID                      ,
		Snap_num                          ,
		Next_coprogenitor_depthfirst_ID   ,
		Last_progenitor_depthfirst_ID     ,
		Last_mainleaf_depthfirst_ID       ,
		Tidal_Force                       ,
		Tidal_ID                          ,
		Rs_Klypin                         ,
		Mmvir_all                         ,
		M200b                             ,
		M200c                             ,
		M500c                             ,
		M2500c                            ,
		Xoff                              ,
		Voff                              ,
		Spin_Bullock                      ,
		b_to_a                            ,
		c_to_a                            ,
		Ax                                ,
		Ay                                ,
		Az                                ,
		b_to_a_500c                       ,
		c_to_a_500c                       ,
		Ax_500c                           ,
		Ay_500c                           ,
		Az_500c                           ,
		TU                                ,
		M_pe_Behroozi                     ,
		M_pe_Diemer                       ,
		Macc                              ,
		Mpeak                             ,
		Vacc                              ,
		Vpeak                             ,
		Halfmass_Scale                    ,
		Acc_Rate_Inst                     ,
		Acc_Rate_100Myr                   ,
		Acc_Rate_1Tdyn                    ,
		Acc_Rate_2Tdyn                    ,
		Acc_Rate_Mpeak                    ,
		Mpeak_Scale                       ,
		Acc_Scale                         ,
		First_Acc_Scale                   ,
		First_Acc_Mvir                    ,
		First_Acc_Vmax                    ,
		VmaxAtMpeak                       ,
		Tidal_Force_Tdyn                  ,
		logVmaxVmaxmaxTdynTmpeak          ,
		Time_to_future_merger             ,
		Future_merger_MMP_ID      =  n.loadtxt(path_2_snap, unpack=True, dtype = dtype1)
		
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
		#f.attrs['age_yr']    = float(age_yr)#[snap_id] 
		#f.attrs['rho_crit']  = float(rho_crit)#[snap_id] 
		#f.attrs['delta_vir'] = float(delta_vir)#[snap_id]

		halo_data = f.create_group('halo_position')
		
		ds = halo_data.create_dataset('x', data = x )
		ds.attrs['units'] = 'Mpc/h'
		ds.attrs['long_name'] = 'x' 
		ds = halo_data.create_dataset('y', data = y )
		ds.attrs['units'] = 'Mpc/h'
		ds.attrs['long_name'] = 'y' 
		ds = halo_data.create_dataset('z', data = z )
		ds.attrs['units'] = 'Mpc/h'
		ds.attrs['long_name'] = 'z' 

		ds = halo_data.create_dataset('vx', data = vx )
		ds.attrs['units'] = 'km/s'
		ds.attrs['long_name'] = 'vx' 
		ds = halo_data.create_dataset('vy', data = vy )
		ds.attrs['units'] = 'km/s'
		ds.attrs['long_name'] = 'vy' 
		ds = halo_data.create_dataset('vz', data = vz )
		ds.attrs['units'] = 'km/s'
		ds.attrs['long_name'] = 'vz' 
		
		halo_data = f.create_group('halo_properties')
		halo_data.attrs['N_halos'] =  N_halo

		ds = halo_data.create_dataset('scale'                              , data=scale                            )
		ds = halo_data.create_dataset('id'                                 , data=id                               )
		ds = halo_data.create_dataset('desc_scale'                         , data=desc_scale                       )
		ds = halo_data.create_dataset('desc_id'                            , data=desc_id                          )
		ds = halo_data.create_dataset('num_prog'                           , data=num_prog                         )
		ds = halo_data.create_dataset('pid'                                , data=pid                              )
		ds = halo_data.create_dataset('upid'                               , data=upid                             )
		ds = halo_data.create_dataset('desc_pid'                           , data=desc_pid                         )
		ds = halo_data.create_dataset('phantom'                            , data=phantom                          )
		ds = halo_data.create_dataset('sam_mvir'                           , data=sam_mvir                         )
		ds = halo_data.create_dataset('mvir'                               , data=mvir                             )
		ds = halo_data.create_dataset('rvir'                               , data=rvir                             )
		ds = halo_data.create_dataset('rs'                                 , data=rs                               )
		ds = halo_data.create_dataset('vrms'                               , data=vrms                             )
		ds = halo_data.create_dataset('mmp'                               , data=mmp                             )
		ds = halo_data.create_dataset('scale_of_last_MM'                   , data=scale_of_last_MM                 )
		ds = halo_data.create_dataset('vmax'                               , data=vmax                             )
		ds = halo_data.create_dataset('Jx'                                 , data=Jx                               )
		ds = halo_data.create_dataset('Jy'                                 , data=Jy                               )
		ds = halo_data.create_dataset('Jz'                                 , data=Jz                               )
		ds = halo_data.create_dataset('Spin'                               , data=Spin                             )
		ds = halo_data.create_dataset('Breadth_first_ID'                   , data=Breadth_first_ID                 )
		ds = halo_data.create_dataset('Depth_first_ID'                     , data=Depth_first_ID                   )
		ds = halo_data.create_dataset('Tree_root_ID'                       , data=Tree_root_ID                     )
		ds = halo_data.create_dataset('Orig_halo_ID'                       , data=Orig_halo_ID                     )
		ds = halo_data.create_dataset('Snap_num'                           , data=Snap_num                         )
		ds = halo_data.create_dataset('Next_coprogenitor_depthfirst_ID'    , data=Next_coprogenitor_depthfirst_ID  )
		ds = halo_data.create_dataset('Last_progenitor_depthfirst_ID'      , data=Last_progenitor_depthfirst_ID    )
		ds = halo_data.create_dataset('Last_mainleaf_depthfirst_ID'        , data=Last_mainleaf_depthfirst_ID      )
		ds = halo_data.create_dataset('Tidal_Force'                        , data=Tidal_Force                      )
		ds = halo_data.create_dataset('Tidal_ID'                           , data=Tidal_ID                         )
		ds = halo_data.create_dataset('Rs_Klypin'                          , data=Rs_Klypin                        )
		ds = halo_data.create_dataset('Mmvir_all'                          , data=Mmvir_all                        )
		ds = halo_data.create_dataset('M200b'                              , data=M200b                            )
		ds = halo_data.create_dataset('M200c'                              , data=M200c                            )
		ds = halo_data.create_dataset('M500c'                              , data=M500c                            )
		ds = halo_data.create_dataset('M2500c'                             , data=M2500c                           )
		ds = halo_data.create_dataset('Xoff'                               , data=Xoff                             )
		ds = halo_data.create_dataset('Voff'                               , data=Voff                             )
		ds = halo_data.create_dataset('Spin_Bullock'                       , data=Spin_Bullock                     )
		ds = halo_data.create_dataset('b_to_a'                             , data=b_to_a                           )
		ds = halo_data.create_dataset('c_to_a'                             , data=c_to_a                           )
		ds = halo_data.create_dataset('Ax'                                 , data=Ax                               )
		ds = halo_data.create_dataset('Ay'                                 , data=Ay                               )
		ds = halo_data.create_dataset('Az'                                 , data=Az                               )
		ds = halo_data.create_dataset('b_to_a_500c'                        , data=b_to_a_500c                      )
		ds = halo_data.create_dataset('c_to_a_500c'                        , data=c_to_a_500c                      )
		ds = halo_data.create_dataset('Ax_500c'                            , data=Ax_500c                          )
		ds = halo_data.create_dataset('Ay_500c'                            , data=Ay_500c                          )
		ds = halo_data.create_dataset('Az_500c'                            , data=Az_500c                          )
		ds = halo_data.create_dataset('TU'                                 , data=TU                               )
		ds = halo_data.create_dataset('M_pe_Behroozi'                      , data=M_pe_Behroozi                    )
		ds = halo_data.create_dataset('M_pe_Diemer'                        , data=M_pe_Diemer                      )
		ds = halo_data.create_dataset('Macc'                               , data=Macc                             )
		ds = halo_data.create_dataset('Mpeak'                              , data=Mpeak                            )
		ds = halo_data.create_dataset('Vacc'                               , data=Vacc                             )
		ds = halo_data.create_dataset('Vpeak'                              , data=Vpeak                            )
		ds = halo_data.create_dataset('Halfmass_Scale'                     , data=Halfmass_Scale                   )
		ds = halo_data.create_dataset('Acc_Rate_Inst'                      , data=Acc_Rate_Inst                    )
		ds = halo_data.create_dataset('Acc_Rate_100Myr'                    , data=Acc_Rate_100Myr                  )
		ds = halo_data.create_dataset('Acc_Rate_1Tdyn'                     , data=Acc_Rate_1Tdyn                   )
		ds = halo_data.create_dataset('Acc_Rate_2Tdyn'                     , data=Acc_Rate_2Tdyn                   )
		ds = halo_data.create_dataset('Acc_Rate_Mpeak'                     , data=Acc_Rate_Mpeak                   )
		ds = halo_data.create_dataset('Mpeak_Scale'                        , data=Mpeak_Scale                      )
		ds = halo_data.create_dataset('Acc_Scale'                          , data=Acc_Scale                        )
		ds = halo_data.create_dataset('First_Acc_Scale'                    , data=First_Acc_Scale                  )
		ds = halo_data.create_dataset('First_Acc_Mvir'                     , data=First_Acc_Mvir                   )
		ds = halo_data.create_dataset('First_Acc_Vmax'                     , data=First_Acc_Vmax                   )
		ds = halo_data.create_dataset('VmaxAtMpeak'                        , data=VmaxAtMpeak                      )
		ds = halo_data.create_dataset('Tidal_Force_Tdyn'                   , data=Tidal_Force_Tdyn                 )
		ds = halo_data.create_dataset('logVmaxVmaxmaxTdynTmpeak'           , data=logVmaxVmaxmaxTdynTmpeak         )
		ds = halo_data.create_dataset('Time_to_future_merger'              , data=Time_to_future_merger            )
		ds = halo_data.create_dataset('Future_merger_MMP_ID'               , data=Future_merger_MMP_ID             )

		f.close()
	
