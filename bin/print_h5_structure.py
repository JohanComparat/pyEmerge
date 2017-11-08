import sys
file_name = sys.argv[1]
# python3 print_data_structure.py filename
import glob
import os
import numpy as n

import h5py    # HDF5 support

f0 = h5py.File(file_name,  "r")

def print_attr(h5item):
  for attr in h5item:
    print(attr, h5item[attr])
      
def print_all_key(h5item):
  for key in h5item.keys():
      print('========================================')
      print(key, h5item[key])
      print('- - - - - - - - - - - - - - - - - - - - ')
      print_attr(h5item[key])
    
def print_data_structure(h5item):
  print('+ + + + + + + HEADER + + + + + + + + +')
  print_attr(h5item.attrs)
  print('\n')
  print('+ + + + + + + DATA   + + + + + + + + + +')
  print_all_key(h5item)

print_data_structure(f0)

