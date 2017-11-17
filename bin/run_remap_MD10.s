#!/bin/bash
#python3 remap_lc.py 24  MD10 1000.
#python3 remap_lc.py 28  MD10 1000.
python3 remap_lc.py 55  MD10 1000.
python3 remap_lc.py 56  MD10 1000.
python3 remap_lc.py 99  MD10 1000.
python3 remap_lc.py 77  MD10 1000.
python3 remap_lc.py 100 MD10 1000.
python3 remap_lc.py 112 MD10 1000.

# failing because 
/*Traceback (most recent call last):
  File "remap_lc.py", line 86, in <module>
    out3 = p.starmap(f3, n.transpose([x0, y0, z0]))
  File "/usr/lib/python3.4/multiprocessing/pool.py", line 268, in starmap
    return self._map_async(func, iterable, starmapstar, chunksize).get()
  File "/usr/lib/python3.4/multiprocessing/pool.py", line 599, in get
    raise self._value
RuntimeError: (0.216542, 0.719859, 0.496683) not contained in any cell*/


python3 remap_lc.py 72  MD10 1000.
python3 remap_lc.py 76  MD10 1000.
python3 remap_lc.py 98 MD10 1000.
python3 remap_lc.py 99 MD10 1000.
python3 remap_lc.py 111 MD10 1000.

