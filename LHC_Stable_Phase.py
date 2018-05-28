# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import TimberManager as tm
import numpy as np

class Phase:
    def __init__(self, timber_variable):

        if type(timber_variable) is str:
            dict_timber = parse_csv_file(timber_variable) 

        elif type(timber_variable) is dict:
            dict_timber = timber_variable

        timber_variable_phase = dict_timber['PHASE'] 
        timber_variable_bunches = dict_timber['BUNCHES']
        
        self.t_stamps = np.float_(np.array(timber_variable_phase.t_stamps)) 
        values_raw = np.float_(np.array(timber_variable_phase.values)) 
        bunches_raw = np.float_(np.array(timber_variable_bunches.values)) 

        # create variables with all bunch slots
        nslots = 3564
        self.bunches = np.arange(nslots)
        self.values = np.zeros((len(self.t_stamps), nslots))
        for ii in xrange(len(bunches_raw)):
            bunch = bunches_raw[ii]
            self.values[:,bunch-1] = values_raw[:,ii] 

    def stick_power_loss(self, fbct_obj, V_rf=6e6, bct_obj=None):
        qe = 1.602176565e-19
        T_rev = 88.9e-6 
        energy_loss_part = -V_rf * np.sin(self.values*np.pi/180.)
        energy_loss_bunch = 0*energy_loss_part
        n_traces = len(self.t_stamps)
        for ii in xrange(n_traces):
            fbct_raw, t_fbct_curr = fbct_obj.nearest_older_sample(self.t_stamps[ii], flag_return_time=True)
            if bct_obj is not None:
                bct_curr = bct_obj.nearest_older_sample(t_fbct_curr)
                fbct_curr = fbct_raw/np.float_(np.sum(fbct_raw))*bct_curr
            else:
                fbct_curr = fbct_raw

            energy_loss_bunch[ii,:] = energy_loss_part[ii,:] * fbct_curr 
        
        self.power_loss = energy_loss_bunch*qe/T_rev

    def nearest_older_sample_power_loss(self, t_obs, flag_return_time=False):
        ind_min = np.argmin(np.abs(self.t_stamps - t_obs))
        if self.t_stamps[ind_min] > t_obs:
            ind_min -= 1
        if flag_return_time:	
            if ind_min == -1:
                return 0.*self.power_loss[ind_min], -1
            else:	
                return self.power_loss[ind_min], self.t_stamps[ind_min]
        else:
            if ind_min == -1:
                return 0.*self.power_loss[ind_min]
            else:	
                return self.power_loss[ind_min]
	
	#~ @property #never tested, so I comment it out (Gianni)
	#~ def total_power_loss(self):
		#~ try:
			#~ tot = np.sum(self.power_loss, axis=1)
		#~ except NameError:
			#~ raise NameError ('power_loss not there!\nPlease run "stick_power_loss"')
		#~ return tot

def parse_csv_file(filename):

    with open(filename) as fid:
        csv_lines = fid.readlines() 

    bunches = tm.timber_variable_list()
    phase = tm.timber_variable_list()

    N_lines = len(csv_lines)
    i_ln = 1

    line = csv_lines[i_ln]
    line = line.split('\r\n')[0]
    bunch_list = line[1:].split(',')
    bunches.values = bunch_list
    i_ln += 1

    while i_ln < N_lines:
        line = csv_lines[i_ln]
        line = line.replace('--', '0.')
        line = line.split('\r\n')[0]
        line_obj = csv_data_line(line)
        (phase.t_stamps).append(line_obj.timestamp) 
        (phase.values).append(line_obj.data_strings)
        i_ln += 1

    variables = {}
    variables['PHASE'] = phase
    variables['BUNCHES'] = bunches

    return variables


class csv_data_line():
    def __init__(self, rline):
        list_values = rline.split(',')
        t_string = list_values[0]
        self.timestamp = tm.timb_timestamp2float(t_string.split('.')[0])
        self.data_strings = list_values[1:]
        
class PowerLoss:
    def __init__(self, timber_variable):

        if type(timber_variable) is str:
            dict_timber = parse_csv_file(timber_variable) 

        elif type(timber_variable) is dict:
            dict_timber = timber_variable
            
        if 'filln' in dict_timber.keys():
            
            filln = dict_timber['filln']
            beam = dict_timber['beam']
            
            print 'ObsBox mode'
            from scipy.constants import e as qe
            import pytimber
            
            ldb = pytimber.LoggingDB(source='ldb')
                       
            fillinfo = ldb.getLHCFillData(filln)
            t_start = fillinfo['startTime']
            t_stop = fillinfo['endTime']
            
            
            cav_list = ['ACSCA.UX45.L%dB%d:CAV_FIELD'%(icav, beam) for icav in range(1, 9)]

            print 'Downloading data. This might take a while...'
            data = {}
            data.update(ldb.getAligned([u'LONGDIAG.SR4.B%d:STABLE_PHASE_MEAN_REL'%beam, u'LHC.BCTFR.A6R4.B%d:BUNCH_INTENSITY'%beam]+cav_list, 
                master='LONGDIAG.SR4.B%d:STABLE_PHASE_MEAN_REL'%beam, t1=t_start, t2=t_stop))
            print 'Done downloading data.'
            
            timestamps = data['timestamps']
            intensity = data['LHC.BCTFR.A6R4.B%d:BUNCH_INTENSITY'%beam]
            phase = data['LONGDIAG.SR4.B%d:STABLE_PHASE_MEAN_REL'%beam]
            total_voltage = np.sum(np.array([data[var] for var in cav_list]), axis=0)*1e6
            
            T_rev = 88.9e-6 
            n_slots = len(phase[0])
            self.power_loss = -qe*np.array(n_slots*[total_voltage]).T*intensity*np.sin(phase*np.pi/360.)/T_rev
            self.t_stamps = timestamps
            self.bunches = np.arange(n_slots)
            
        else:
            timber_variable_ploss = dict_timber['PHASE']  #abusing a bit the csv_parser
            timber_variable_bunches = dict_timber['BUNCHES']
            
            self.t_stamps = np.float_(np.array(timber_variable_ploss.t_stamps)) 
            values_raw = np.float_(np.array(timber_variable_ploss.values)) 
            bunches_raw = np.float_(np.array(timber_variable_bunches.values)) 
            
            if len(bunches_raw) == (values_raw.shape[1]+1):
                print 'WARNING: first bunch seems not to be there! I remove it from the bunch list.'
                bunches_raw = bunches_raw[1:]

            # create variables with all bunch slots
            nslots = 3564
            self.bunches = np.arange(nslots)
            self.power_loss = np.zeros((len(self.t_stamps), nslots))
            for ii in xrange(len(bunches_raw)):
                bunch = bunches_raw[ii]
                self.power_loss[:,bunch-1] = values_raw[:,ii] 

   

    def nearest_older_sample_power_loss(self, t_obs, flag_return_time=False):
        ind_min = np.argmin(np.abs(self.t_stamps - t_obs))
        if self.t_stamps[ind_min] > t_obs:
            ind_min -= 1
        if flag_return_time:	
            if ind_min == -1:
                return 0.*self.power_loss[ind_min], -1
            else:	
                return self.power_loss[ind_min], self.t_stamps[ind_min]
        else:
            if ind_min == -1:
                return 0.*self.power_loss[ind_min]
            else:	
                return self.power_loss[ind_min]
	
    @property
    def total_power_loss(self):
		return np.sum(self.power_loss, axis=1)
		
