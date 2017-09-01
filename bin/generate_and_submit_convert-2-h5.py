import os

batch_dir = os.path.join(os.environ['HOME'], 'batch_emerge')
os.system("cd "+batch_dir)

import glob 
import numpy as n

aexp, redshift, age_yr, rho_crit, delta_vir = n.loadtxt("/u/joco/data/MD/MD_1.0Gpc/hlists_MD_1.0Gpc.ascii", unpack=True)
snap_names = n.array([ el.zfill(6)[0]+"."+el.zfill(6)[1:] for el in (aexp*10**5).astype('int').astype('str') ])
#snap_id = -10
#snap_name = snap_names[snap_id]


main_text = lambda snap_name : """
# @ shell=/bin/bash
#
# Sample script for LoadLeveler
#
# @ error = /u/joco/run_monitor/h5_"""+snap_name+""".err.$(jobid)
# @ output = /u/joco/run_monitor/h5_"""+snap_name+""".out.$(jobid)
# @ job_type = parallel
# @ node_usage= not_shared
# @ node = 1
# @ tasks_per_node = 1
# @ resources = ConsumableCpus(1)
# @ wall_clock_limit = 00:30:00
# @ notification = complete
# @ notify_user = comparat@mpe.mpg.de
# @ queue 

source /etc/profile.d/modules.sh
module load anaconda
export PYTHONPATH="${PYTHONPATH}:/u/joco/software/pyEmerge/python/"

cd /u/joco/software/pyEmerge/bin/
"""

for snap_name in snapshot_names[:-10][::-1][:3]:
    print(snap_name)
    batch_file = os.path.join(batch_dir, "h5_"+snap_name+".sh")
    f=open( batch_file, 'w')
    f.write(main_text(snap_name))
    command = "python MD10-write-convert-2-h5.py "+snap_name
    print(command)
    f.write(command)
    f.close()
    
    #os.system("llsubmit " + batch_file )

