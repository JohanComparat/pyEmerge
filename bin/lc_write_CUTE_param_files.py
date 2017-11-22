"""
Writes clustering samples: ra, dec, z for a set of LX cuts
Based on the MDPL lightcones.
"""
import os

Lname='L3'

topdir = '/data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_remaped_position_'+Lname+'/'

def write_param():
# input-output files and parameters
data_filename= /data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_remaped_position_L3/lc_L3_z_lt_04_lx_gt_430.ascii     
random_filename= /data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_remaped_position_L3/lc_L3_z_lt_04_lx_gt_430.random
input_format= 2
output_filename= /data17s/darksim/MD/MD_1.0Gpc/h5_lc/clustering_catalogs_remaped_position_L3/lc_L3_z_lt_04_lx_gt_430.2pcf

# estimation parameters
corr_type= monopole

# cosmological parameters
omega_M= 0.307
omega_L= 0.693
w= -1

# binning
log_bin= 1
dim1_max= 100.
dim1_min_logbin= 0.01
dim1_nbin= 15

dim2_max= 100.
dim2_nbin= 15

dim3_min= 0.01
dim3_max= 1.3
dim3_nbin= 1


write_samples("L3")
write_samples("L6")
write_samples("L15")
