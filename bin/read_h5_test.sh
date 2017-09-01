
import h5py    # HDF5 support

fileName = "hlist_0.08550_emerge.hdf5"
f = h5py.File(fileName,  "r")
for item in f.attrs.keys():
    print( item + ":", f.attrs[item])

f['/halo_data/id'].value

f.close()
