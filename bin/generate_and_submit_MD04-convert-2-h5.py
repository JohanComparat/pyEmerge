import os

batch_dir = os.path.join(os.environ['HOME'], 'batch_emerge')
os.system("cd "+batch_dir)

import glob 
import numpy as n

aexp, redshift, age_yr, rho_crit, delta_vir = n.loadtxt("/u/joco/data/MD/MD_0.4Gpc/hlists_MD_0.4Gpc.ascii", unpack=True)
snap_names = n.array([ el.zfill(5)[0]+"."+el.zfill(5)[1:]+"0" for el in (aexp*10**4).astype('int').astype('str') ])
snap_ids = n.arange(len(aexp))
#snap_name = snap_names[snap_id]

import sys

main_text = lambda snap_name : """
# @ shell=/bin/bash
#
# Sample script for LoadLeveler
#
# @ error = /u/joco/run_monitor/smd_h5_"""+snap_name+""".err.$(jobid)
# @ output = /u/joco/run_monitor/smd_h5_"""+snap_name+""".out.$(jobid)
# @ job_type = parallel
# @ node_usage= not_shared
# @ node = 1
# @ tasks_per_node = 1
# @ resources = ConsumableCpus(1)
# @ wall_clock_limit = 01:00:00
# @ notification = error
# @ notify_user = comparat@mpe.mpg.de
# @ queue 

source /etc/profile.d/modules.sh
module load anaconda
export PYTHONPATH="${PYTHONPATH}:/u/joco/software/pyEmerge/python/"

cd /u/joco/software/pyEmerge/bin/
"""
print("snap ids considered", snap_ids )
for snap_id in snap_ids:
    print(snap_id)
    batch_file = os.path.join(batch_dir, "h5_"+snap_names[snap_id]+"_smd.sh")
    f=open( batch_file, 'w')
    f.write(main_text(snap_names[snap_id]))
    command = "python MD04-convert-2-h5.py " + snap_names[snap_id] +" "+str(aex[psnap_id])+" "+str( redshift[snap_id])+" "+str( age_yr[snap_id])+" "+str( rho_crit[snap_id])+" "+str( delta_vir[snap_id])
    print(command)
    f.write(command)
    f.close()
    os.system("llsubmit " + batch_file )

