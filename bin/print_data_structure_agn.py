import sys
ii = int(sys.argv[1])
env = sys.argv[2]
# python3 print_data_structure.py 22 MD10
import glob
import os
import numpy as n
import EmergeIterate

iterate = EmergeIterate.EmergeIterate(ii, env)
iterate.open_snapshots()

print(ii,iterate.f0.attrs['aexp'])
print(iterate.f0['/agn_properties/log_lambda_sar'].value[:10])

