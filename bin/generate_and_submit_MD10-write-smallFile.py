import os

batch_dir = os.path.join(os.environ['HOME'], 'batch_emerge')
os.system("cd "+batch_dir)

import glob 
import numpy as n

input_list = n.array(glob.glob("/u/joco/data/MD/MD_1.0Gpc/hlists/hlist_*.list"))
input_list.sort()
snapshot_names = n.array([ os.path.basename(el)[:-5] for el in input_list ])
output_list = n.array([ os.path.join("/u/joco/data/MD/MD_1.0Gpc/emerge/", name+".data") for name in snapshot_names])

main_text = lambda snap_name : """
# @ shell=/bin/bash
#
# Sample script for LoadLeveler
#
# @ error = /u/joco/run_monitor/"""+snap_name+""".err.$(jobid)
# @ output = /u/joco/run_monitor/"""+snap_name+""".out.$(jobid)
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

for path_2_in, path_2_out, snap_name in zip(input_list,output_list, snapshot_names)[:13]:
    print(path_2_in, path_2_out, snap_name)
    batch_file = os.path.join(batch_dir, snap_name+".sh")
    f=open( batch_file, 'w')
    f.write(main_text(snap_name))
    command = "python MD10-write-smallFile.py "+path_2_in+" "+path_2_out
    print(command)
    f.write(command)
    f.close()
    
    os.system("llsubmit " + batch_file )

