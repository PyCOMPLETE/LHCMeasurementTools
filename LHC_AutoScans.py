import numpy as np
import TimberManager as tm

class AutoScan:
	def __init__(self, timber_variable):

		if type(timber_variable) is str:
			print 'type var'
			# if not (IP == 1 or IP == 5):
			# 	raise ValueError('You need to specify which IP! (1 or 5)')
			dict_timber = tm.parse_timber_file(timber_variable, verbose=True)
			timber_variable_SCAN = dict_timber[get_variable_dict['AUTOSCAN']]
			# if IP == 1:
			# 	timber_variable_ENABLE = dict_timber[get_variable_dict['LUMILEVELING.IP1:ENABLE']]
			# elif IP == 5: 
			# 	timber_variable_ENABLE = dict_timber[get_variable_dict['LUMILEVELING.IP5:ENABLE']]

		elif hasattr(timber_variable, '__getitem__'):
			print 'type dic'
			try:
				timber_variable_SCAN = timber_variable['AUTOMATICSCAN:IP'] #get_variable_dict['AUTOSCAN']]
				# if IP == 1:
				# 	timber_variable_ENABLE = timber_variable['LUMILEVELING.IP1:ENABLE']
				# elif IP == 5: 
				# 	timber_variable_ENABLE = timber_variable['LUMILEVELING.IP5:ENABLE']
			except:
				print '# LHC_AutoScans : No Emittance Scan Information! Returning -1'
				return -1

		self.t_stamps        = timber_variable_SCAN.t_stamps
		self.autoscanIP      = timber_variable_SCAN.values


		self.autoscanIP     = np.array(np.float_(self.autoscanIP)).ravel()
		self.t_stamps            = np.array(np.float_(self.t_stamps))

	def nearest_older_sample(self, t_obs, flag_return_time=False):
		ind_min = np.argmin(np.abs(self.t_stamps - t_obs))
		if self.t_stamps[ind_min] > t_obs:
			ind_min -= 1
		if flag_return_time:    
			if ind_min == -1:
				return 0.*self.autoscanIP[ind_min], -1
			else:   
				return self.autoscanIP[ind_min], self.t_stamps[ind_min]
		else:
			if ind_min == -1:
				return 0.*self.autoscanIP[ind_min]
			else:   
				return self.autoscanIP[ind_min]
				

		
def get_variable_dict():
	var_dict = {}
	var_dict['AUTOSCAN'] = 'AUTOMATICSCAN:IP'
	#var_dict['LUMI_LEVEL_ENABLE_IP5'] = 'LUMILEVELING.IP5:ENABLE'

	return var_dict

def variable_list():
	var_list = []
	var_list += get_variable_dict().values()

	return var_list
