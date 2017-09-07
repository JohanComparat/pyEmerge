import os
import numpy as n
import sys

batch_dir = os.path.join(os.environ['HOME'], 'batch_emerge')
os.system("cd "+batch_dir)
batch_file = os.path.join(batch_dir, "emerge_iterate_batch.sh")


main_text = """
# @ shell=/bin/bash
#
# Sample script for LoadLeveler
#
# @ error = /u/joco/run_monitor/emerge_iterate.err.$(jobid)
# @ output = /u/joco/run_monitor/emerge_iterate.out.$(jobid)
# @ job_type = parallel
# @ node_usage= not_shared
# @ node = 1
# @ tasks_per_node = 1
# @ resources = ConsumableCpus(1)
# @ wall_clock_limit = 20:00:00
# @ notification = error
# @ notify_user = comparat@mpe.mpg.de
# @ queue 

source /etc/profile.d/modules.sh
module load anaconda
export PYTHONPATH="${PYTHONPATH}:/u/joco/software/pyEmerge/python/"

cd /u/joco/software/pyEmerge/bin/
"""


f=open( batch_file, 'w')
f.write(main_text)

for num in n.arange(26, 109, 1):
	command = "python emerge_iterate.py " + str(num)+' \n'
	f.write(command)
	command2 = "touch ~/batch_emerge/emerge_iterate_" + str(num)+' \n'
	f.write(command2)

f.close()


