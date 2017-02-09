import numpy as np

class obj_from_dict:
    def __init__(self, dictto):
        for key, value in dictto.iteritems():
            setattr(self, key, value)

def obj_to_dict(obj):
    dict_out={}
    members = dir(obj)
    for member in members:
        dict_out[member] = getattr(obj, member)
    return dict_out

def myloadmat(filename, squeeze = True):
    import scipy.io as sio
    dict_var=sio.loadmat(filename)
    if squeeze:
        for kk in dict_var.keys():
            try:
                dict_var[kk]=np.squeeze(dict_var[kk])
            except:
                pass
    return dict_var

def myloadmat_to_obj(filename, squeeze = True):
    return  obj_from_dict(myloadmat(filename, squeeze=squeeze))

def dict_of_arrays_and_scalar_from_h5(filename):
    import h5py
    with h5py.File(filename) as fid:
        f_dict = {}
        for kk in fid.keys():
            f_dict[kk] = np.array(fid[kk]).copy()
            if f_dict[kk].shape == ():
                f_dict[kk] = f_dict[kk].tolist()
    return  f_dict

def object_with_arrays_and_scalar_from_h5(filename):
    return  obj_from_dict(dict_of_arrays_and_scalar_from_h5(filename))

def bunchh5_to_dict(filename):
    import h5py
    with h5py.File(filename, 'r') as bunch_ev:
        bunch = bunch_ev['Bunch']
        bunch_dict = {}
        for kk in bunch.keys():
            bunch_dict[kk] = np.array(bunch[kk]).copy()

    return bunch_dict

def bunchh5_to_obj(filename):
    return  obj_from_dict(bunchh5_to_dict(filename))

def bunchh5list_to_dict(filename_list):
    bunch_dict = bunchh5_to_dict(filename_list[0])
    for i_file in xrange(1,len(filename_list)):
        bunch_dict_curr = bunchh5_to_dict(filename_list[i_file])
        for kk in bunch_dict.keys():
            bunch_dict[kk] = np.array(list(bunch_dict[kk])+list(bunch_dict_curr[kk]))

    return bunch_dict

def bunchh5list_to_obj(filename_list):
    return  obj_from_dict(bunchh5list_to_dict(filename_list))

# Only works for not nested h5 files
def h5_to_obj(filename):
    import h5py
    d = {}
    with h5py.File(filename, 'r') as f:
        for key in f:
            d[key] = np.array(f[key])
    return obj_from_dict(d)

# Only works for not nested attributes of object
def aligned_obj_to_h5(obj, h5):
    import h5py
    with h5py.File(h5, 'w') as h5_handle:
        h5_handle.create_dataset('timestamps', data=obj.timestamps)
        h5_handle.create_dataset('variables', data=obj.variables)
        h5_handle.create_dataset('data', data=obj.data)
