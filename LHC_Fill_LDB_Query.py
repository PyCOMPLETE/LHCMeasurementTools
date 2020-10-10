import json
import os
# from . import lhc_log_db_query as lldb
import numpy as np

import LHCMeasurementTools.TimberManager as tm

def load_fill_dict_from_json(fname):
    with open(fname, 'r') as fid:
        ddd = json.load(fid)
        fill_dict = {int(kk): ddd[kk] for kk in ddd.keys()}
    return fill_dict

def save_variables_and_json(varlist, file_path_prefix,
        save_json, fills_dict, fill_sublist=None,
        save_to_json=True, save_json_every = 0, db=None,
        n_vars_per_extraction=1):

    if os.path.isfile(save_json):
        saved_fills =  load_fill_dict_from_json(save_json)
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

        try:
            print(f'pytimber source: {db._source}')
        except Exception:
            pass

        data = {}
        for ii in range(0, len(varlist), n_vars_per_extraction):
            thesevars = varlist[ii: ii + n_vars_per_extraction]
            print(f'{ii}/{len(varlist)}: {thesevars[0]} ... {thesevars[-1]}', end='\r', flush=True)
            data.update(tm.CalsVariables_from_pytimber(
                db.get(thesevars, t_start_fill, t_end_fill)))
        print('Done downloading')

        print('Saving h5...')
        tm.CalsVariables_to_h5(data, fill_file)
        print('Done')

        if fills_dict[filln]['flag_complete'] is True:
            saved_fills[filln] = 'complete'
        else:
            saved_fills[filln] = t_end_fill

        if save_json_every>0:
            if int(np.mod(i_fill, save_json_every))==0:
                if save_to_json is True:
                    with open(save_json, 'w') as fid:
                        json.dump(saved_fills, fid)
                    print('\nSaved json!\n')


    if save_to_json is True:
        with open(save_json, 'w') as fid:
                json.dump(saved_fills, fid)

