import h5py

import pytimber
import LHCMeasurementTools.TimberManager as tm

class data_period(object):
    def __init__(self, t_start, t_end, db):
        self.t_start = t_start
        self.t_end = t_end
        self.db = db

    def __getitem__(self, kk):
        var = tm.timber_variable_list()
        var.t_stamps, var.values = self.db.get(kk, self.t_start, self.t_end)[kk]
        return var

def save_to_h5(filename, fill_dict):
    with h5py.File(filename, 'w')as h5_handle:
        for varname, timber_variable_list in fill_dict.items():
            h5_handle.create_dataset(varname+'!t_stamps', data=timber_variable_list.t_stamps)
            h5_handle.create_dataset(varname+'!values', data=timber_variable_list.values)

def get_fill_dict(variable_list, t_start, t_end, db):
    data_period_obj = data_period(t_start, t_end, db)

    return {var: data_period_obj[var] for var in variable_list}

def save_variables_to_h5(filename, variable_list, t_start, t_end, db):
    fill_dict = get_fill_dict(variable_list, t_start, t_end, db)
    save_to_h5(filename, fill_dict)

