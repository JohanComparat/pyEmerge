import time
t0 = time.time()
import sys

sim_dir = sys.argv[1]
snap_id = int(sys.argv[2])
# python3 run_emerge_iterate_multiCore.py MD10 22

import numpy as n 
import pandas as pd

import EmergeIterate

from multiprocessing import Pool
t0=time.time()
n_proc=12
pool = Pool(n_proc)

iterate = EmergeIterate.EmergeIterate(snap_id, sim_dir)
iterate.open_snapshots()

iterate.init_new_quantities()
iterate.write_results()
