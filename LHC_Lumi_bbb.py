# -*- coding: utf-8 -*-
import numpy as np
from . import TimberManager as tm

class LUMI:
    def __init__(self, timber_variable, experiment='ATLAS'):

        if type(timber_variable) is str:            
            dict_timber = tm.parse_timber_file(timber_variable, verbose=True)
            timber_variable_LUMI = dict_timber[experiment+':BUNCH_LUMI_INST']
        elif type(timber_variable) is dict:
            timber_variable_LUMI = timber_variable[experiment+':BUNCH_LUMI_INST']            

        self.t_stamps = timber_variable_LUMI.t_stamps
        self.blumi = timber_variable_LUMI.values

        self.blumi = np.array(np.float_(self.blumi))
        self.t_stamps = np.array(np.float_(self.t_stamps))
        self.totlumi = np.sum(self.blumi, axis = 1)
	self.meanlumi = np.average(self.blumi, axis = 1)


    def uniform_time(self, t_inter=60.):
        t_stamps = self.t_stamps
        lumi = self.blumi
        nslots = blumi.shape[1]

        t_stamps_unif = np.arange(np.min(t_stamps), np.max(t_stamps), t_inter)
        blumi_unif = 0.*np.zeros((len(t_stamps_unif), nslots))
        for ii in range(nslots):
            blumi_unif[:,ii] = np.interp(t_stamps_unif, t_stamps, blumi[:,ii])

        return t_stamps_unif, blumi_unif


    def nearest_older_sample(self, t_obs, flag_return_time=False):
        ind_min = np.argmin(np.abs(self.t_stamps - t_obs))
        if self.t_stamps[ind_min] > t_obs:
            ind_min -= 1
        if flag_return_time:	
            if ind_min == -1:
                return 0.*self.blumi[ind_min], -1
            else:	
                return self.blumi[ind_min], self.t_stamps[ind_min]
        else:
            if ind_min == -1:
                return 0.*self.blumi[ind_min]
            else:	
                return self.blumi[ind_min]

		
def get_variable_dict():
    var_dict = {}
    var_dict['BUNCH_LUMI_INST_ATLAS'] = 'ATLAS:BUNCH_LUMI_INST'
    var_dict['BUNCH_LUMI_INST_CMS'] = 'CMS:BUNCH_LUMI_INST'

    return var_dict

def variable_list():
    var_list = []
    
    var_list = list(get_variable_dict().values())

    return var_list
