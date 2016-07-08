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
