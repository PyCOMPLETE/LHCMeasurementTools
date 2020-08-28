import pickle
import os
from . import lhc_log_db_query as lldb
import numpy as np

import LHCMeasurementTools.TimberManager as tm


def save_variables_and_pickle(varlist, file_path_prefix,
        save_pkl, fills_dict, fill_sublist=None,
        save_to_pickle=True, save_pickle_every = 0, db=None):

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


        if filln in list(saved_fills.keys()) and (saved_fills[filln] == 'complete' or
                                            saved_fills[filln] == t_end_fill):
            continue

        fill_file = file_path_prefix + '_%d.h5'%filln
        print('\nSaving fill %d in file %s\n'%(filln, fill_file))

        # Discontinued (save csv using java executable)
        # lldb.dbquery(varlist, t_start_fill, t_end_fill, fill_file)

        print('Start downloading...')
        if db is None:
            import pytimber
            db = pytimber.LoggingDB()
        data = {}
        for ii, vv in enumerate(varlist):
            print(f'{ii+1}/{len(varlist)} = {vv}')
            data.update(tm.CalsVariables_from_pytimber(
                            db.get([vv], t_start_fill, t_end_fill)))

        #data = tm.CalsVariables_from_pytimber(
        #                    db.get(varlist, t_start_fill, t_end_fill))
        print('Done downloading')

        print('Saving h5...')
        tm.CalsVariables_to_h5(data, fill_file)
        print('Done')

        if fills_dict[filln]['flag_complete'] is True:
            saved_fills[filln] = 'complete'
        else:
            saved_fills[filln] = t_end_fill

        if save_pickle_every>0:
            if int(np.mod(i_fill, save_pickle_every))==0:
                if save_to_pickle is True:
                    with open(save_pkl, 'wb') as fid:
                        pickle.dump(saved_fills, fid)
                    print('\nSaved pickle!\n')

    if save_to_pickle is True:
        with open(save_pkl, 'wb') as fid:
                pickle.dump(saved_fills, fid)

    return data
