"""
Writes clustering samples: ra, dec, z for a set of LX cuts
Based on the MDPL lightcones.
"""
import os

def write_param(Lname='L3', zz='03', LL='420'):
	topdir = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_remaped_position_'+Lname+'/'
	file_name = topdir + 'param_'+Lname+'_z_lt_'+zz+'_lx_gt_'+LL+'.ini'
	f=open(file_name,'w')
	f.write('data_filename= /data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_remaped_position_'+Lname+'/lc_'+Lname+'_z_lt_'+zz+'_lx_gt_'+LL+'.ascii \n')    
	f.write('random_filename= /data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_remaped_position_'+Lname+'/lc_'+Lname+'_z_lt_'+zz+'_lx_gt_'+LL+'.random \n')    
	f.write('input_format= 2 \n')
	f.write('output_filename= /data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_remaped_position_'+Lname+'/lc_'+Lname+'_z_lt_'+zz+'_lx_gt_'+LL+'.2pcf \n')
	f.write('corr_type= monopole \n')
	f.write('omega_M= 0.307 \n')
	f.write('omega_L= 0.693 \n')
	f.write('w= -1 \n')
	f.write('log_bin= 1 \n')
	f.write('dim1_max= 100. \n')
	f.write('dim1_min_logbin= 0.01 \n')
	f.write('dim1_nbin= 15 \n')
	f.write('dim2_max= 100. \n')
	f.write('dim2_nbin= 15 \n')
	f.write('dim3_min= 0.01 \n')
	f.write('dim3_max= 1.3 \n')
	f.write('dim3_nbin= 1 \n')
	f.close()

write_param("L6", '03', '415')
write_param("L6", '03', '420')
write_param("L6", '03', '425')
write_param("L6", '03', '430')
write_param("L6", '03', '435')
write_param("L6", '03', '440')
write_param("L6", '04', '415')
write_param("L6", '04', '420')
write_param("L6", '04', '425')
write_param("L6", '04', '430')
write_param("L6", '04', '435')
write_param("L6", '04', '440')

write_param("L15", '03', '415')
write_param("L15", '03', '420')
write_param("L15", '03', '425')
write_param("L15", '03', '430')
write_param("L15", '03', '435')
write_param("L15", '03', '440')
write_param("L15", '04', '415')
write_param("L15", '04', '420')
write_param("L15", '04', '425')
write_param("L15", '04', '430')
write_param("L15", '04', '435')
write_param("L15", '04', '440')

