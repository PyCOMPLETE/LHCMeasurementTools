import os
import cPickle
import copy
import numpy as np

this_directory = os.path.dirname(os.path.abspath(__file__))
dict_file_2016 = this_directory + '/large_heat_load_dict_2016.pkl'
dict_file_2015 = this_directory + '/large_heat_load_dict_2015.pkl'

moment = 'stable_beams'
#moment = 'start_ramp'
#moment = 'sb+2_hrs'

with open(dict_file_2016, 'r') as f:
    main_dict_2016 = cPickle.load(f)
with open(dict_file_2015, 'r') as f:
    main_dict_2015 = cPickle.load(f)

def mask_dict(dictionary, mask):
    new_dict = copy.deepcopy(dictionary)
    _mask_recursively(new_dict,mask)
    return new_dict

def _mask_recursively(dictionary, mask):
    for key in dictionary:
        if type(dictionary[key]) is dict:
            _mask_recursively(dictionary[key],mask)
        else:
            dictionary[key] = dictionary[key][mask]

def merge_dicts(dict1,dict2):
    new_dict = copy.deepcopy(dict1)
    _merge_dicts_recursively(dict1,dict2,new_dict)
    return new_dict

def _merge_dicts_recursively(dict1, dict2, new_dict):
    for key in dict1:
        if type(dict1[key]) is dict:
            _merge_dicts_recursively(dict1[key],dict2[key], new_dict[key])
        elif type(dict1[key]) is np.ndarray:
            new_dict[key] = np.concatenate([dict1[key], dict2[key]])
        else:
            print('Unexpected type %s for key %s!' % (type(dict1[key]), key))

main_dict = merge_dicts(main_dict_2015, main_dict_2016)
