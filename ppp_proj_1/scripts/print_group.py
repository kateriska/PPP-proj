import h5py

filename = '/home/katerina/Documents/PPP-proj/ppp_proj_1/build/file_h_seq.h5'
data = h5py.File(filename, 'r')
for group in data.keys() :
    print (group)
    for dset in data[group].keys():
        print (dset)
        ds_data = data[group][dset] # returns HDF5 dataset object
        print (ds_data)
        print (ds_data.shape, ds_data.dtype)
        arr = data[group][dset][:] # adding [:] returns a numpy array
        print (arr.shape, arr.dtype)
        print (arr)
