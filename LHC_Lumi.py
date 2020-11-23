# -*- coding: utf-8 -*-
import numpy as np
from . import TimberManager as tm

class LUMI:
    def __init__(self, timber_variable, experiment='ATLAS'):

        if type(timber_variable) is str:            
            dict_timber = tm.parse_timber_file(timber_variable, verbose=True)
            timber_variable_LUMI = dict_timber[get_variable_dict()['LUMI_TOT_' + experiment]]
        elif type(timber_variable) is dict:
            timber_variable_LUMI = timber_variable[get_variable_dict()['LUMI_TOT_' + experiment]]            

        self.t_stamps = timber_variable_LUMI.t_stamps
        self.lumi_tot = timber_variable_LUMI.float_values()*1e30

        # self.bint = map(lambda x: np.array(map(float, x)), self.bint)
        # self.bint = np.array(self.bint)
        # self.t_stamps = np.array(self.t_stamps)
        # self.totint = np.array(map(sum, self.bint))

        self.lumi_tot = np.array(np.float_(self.lumi_tot))
        self.t_stamps = np.array(np.float_(self.t_stamps))        


    def uniform_time(self, t_inter=60.):
        t_stamps = self.t_stamps
        lumi = self.lumi_tot
        nslots = lumi_tot.shape[1]

        t_stamps_unif = np.arange(np.min(t_stamps), np.max(t_stamps), t_inter)
        lumi_tot_unif = 0.*np.zeros((len(t_stamps_unif), nslots))
        for ii in range(nslots):
            lumi_tot_unif[:,ii] = np.interp(t_stamps_unif, t_stamps, lumi_tot[:,ii])

        return t_stamps_unif, lumi_tot_unif


    def nearest_older_sample(self, t_obs, flag_return_time=False):
        ind_min = np.argmin(np.abs(self.t_stamps - t_obs))
        if self.t_stamps[ind_min] > t_obs:
            ind_min -= 1
        if flag_return_time:	
            if ind_min == -1:
                return 0.*self.lumi_tot[ind_min], -1
            else:	
                return self.lumi_tot[ind_min], self.t_stamps[ind_min]
        else:
            if ind_min == -1:
                return 0.*self.lumi_tot[ind_min]
            else:	
                return self.lumi_tot[ind_min]

        
def get_variable_dict():
    var_dict = {}
    var_dict['LUMI_TOT_ATLAS'] = 'ATLAS:LUMI_TOT_INST'
    var_dict['LUMI_TOT_CMS'] = 'CMS:LUMI_TOT_INST'

    return var_dict

def variable_list():
    var_list = []
    
    var_list = list(get_variable_dict().values())

    return var_list
