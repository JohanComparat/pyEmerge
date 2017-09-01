# @ shell=/bin/bash
#
# Sample script for LoadLeveler
#
# @ error = /u/joco/run_monitor/test.err.$(jobid)
# @ output = /u/joco/run_monitor/test.out.$(jobid)
# @ job_type = parallel
# @ node_usage= not_shared
# @ node = 1
# @ tasks_per_node = 1
# @ resources = ConsumableCpus(1)
# @ wall_clock_limit = 00:45:00
# @ notification = complete
# @ notify_user = comparat@mpe.mpg.de
# @ queue = single

source /etc/profile.d/modules.sh
module load anaconda
export PYTHONPATH="${PYTHONPATH}:/u/joco/software/pyEmerge/python/"

cd /u/joco/software/pyEmerge/bin/
python MD10-write-smallFile.py /ptmp/joco/MD/MD_1.0Gpc/hlists/hlist_0.06270.list /ptmp/joco/MD/MD_1.0Gpc/emerge/hlist_0.06270.data