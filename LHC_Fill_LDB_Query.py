import pickle
import os
import lhc_log_db_query as lldb
import numpy as np


def save_variables_and_pickle(varlist, file_path_prefix, save_pkl, fills_dict, fill_sublist=None,save_to_pickle=True, save_pickle_every = 0):
    if os.path.isfile(save_pkl):
        with open(save_pkl, 'rb') as fid:
            saved_fills = pickle.load(fid)
    else:
        saved_fills = {}

    for i_fill, filln in enumerate(sorted(fills_dict.keys())):

        if fill_sublist is not None:
            if filln not in fill_sublist:
                continue

        t_start_fill = fills_dict[filln]['t_startfill']
        t_end_fill = fills_dict[filln]['t_endfill']


        if filln in saved_fills.keys() and (saved_fills[filln] == 'complete' or
                                            saved_fills[filln] == t_end_fill):
            continue

        fill_file = file_path_prefix + '_%d.csv'%filln
        print '\nSaving fill %d in file %s\n'%(filln, fill_file)

        lldb.dbquery(varlist, t_start_fill, t_end_fill, fill_file)

        if fills_dict[filln]['flag_complete'] is True:
            saved_fills[filln] = 'complete'
        else:
            saved_fills[filln] = t_end_fill
            
        if save_pickle_every>0:
            if int(np.mod(i_fill, save_pickle_every))==0:
                if save_to_pickle is True:
                    with open(save_pkl, 'wb') as fid:
                        pickle.dump(saved_fills, fid)
                    print '\nSaved pickle!\n'
                
    if save_to_pickle is True:
        with open(save_pkl, 'wb') as fid:
                pickle.dump(saved_fills, fid)
