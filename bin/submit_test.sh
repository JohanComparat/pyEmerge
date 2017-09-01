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
# @ queue 

source /etc/profile.d/modules.sh
module load anaconda
export PYTHONPATH="${PYTHONPATH}:/u/joco/software/pyEmerge/python/"

cd /u/joco/software/pyEmerge/bin/
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.06270.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.06270.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.06410.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.06410.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.06560.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.06560.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.06700.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.06700.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.06850.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.06850.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.07010.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.07010.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.07160.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.07160.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.07320.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.07320.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.07490.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.07490.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.07660.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.07660.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.07830.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.07830.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.08000.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.08000.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.08180.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.08180.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.08360.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.08360.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.08550.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.08550.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.08740.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.08740.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.08940.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.08940.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.09140.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.09140.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.09340.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.09340.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.09550.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.09550.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.09770.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.09770.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.09990.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.09990.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.10210.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.10210.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.10440.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.10440.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.10670.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.10670.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.10910.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.10910.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.11160.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.11160.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.11410.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.11410.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.11660.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.11660.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.11920.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.11920.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.12190.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.12190.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.12460.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.12460.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.12740.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.12740.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.13030.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.13030.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.13320.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.13320.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.13620.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.13620.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.13920.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.13920.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.14240.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.14240.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.14550.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.14550.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.14880.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.14880.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.15210.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.15210.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.15550.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.15550.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.15900.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.15900.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.16260.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.16260.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.16620.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.16620.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.17000.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.17000.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.17380.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.17380.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.17770.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.17770.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.18160.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.18160.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.18570.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.18570.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.18990.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.18990.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.19410.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.19410.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.19850.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.19850.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.20290.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.20290.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.20750.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.20750.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.21210.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.21210.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.21690.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.21690.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.22170.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.22170.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.22670.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.22670.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.23180.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.23180.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.23690.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.23690.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.24230.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.24230.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.24770.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.24770.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.25320.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.25320.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.25890.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.25890.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.26470.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.26470.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.27060.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.27060.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.27670.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.27670.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.28290.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.28290.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.28920.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.28920.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.29570.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.29570.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.30230.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.30230.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.30910.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.30910.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.31600.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.31600.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.32310.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.32310.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.33030.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.33030.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.33770.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.33770.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.34530.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.34530.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.35300.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.35300.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.36090.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.36090.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.36900.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.36900.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.37730.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.37730.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.38570.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.38570.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.39440.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.39440.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.40320.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.40320.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.41230.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.41230.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.42150.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.42150.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.43090.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.43090.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.44060.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.44060.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.45050.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.45050.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.46050.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.46050.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.47090.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.47090.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.48140.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.48140.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.49220.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.49220.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.50320.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.50320.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.51450.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.51450.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.52600.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.52600.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.53780.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.53780.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.54980.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.54980.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.56220.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.56220.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.57470.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.57470.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.58760.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.58760.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.60080.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.60080.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.61420.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.61420.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.62800.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.62800.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.64210.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.64210.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.65650.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.65650.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.67120.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.67120.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.68620.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.68620.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.70160.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.70160.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.71730.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.71730.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.73330.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.73330.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.74980.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.74980.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.76660.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.76660.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.78370.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.78370.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.80130.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.80130.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.81920.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.81920.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.83760.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.83760.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.85640.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.85640.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.87550.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.87550.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.89510.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.89510.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.91520.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.91520.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.93570.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.93570.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.95670.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.95670.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_0.97810.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_0.97810.data
python MD10-write-smallFile.py /u/joco/data/MD/MD_1.0Gpc/hlists/hlist_1.00000.list /u/joco/data/MD/MD_1.0Gpc/emerge/hlist_1.00000.data

