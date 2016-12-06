import numpy as np
import TimberManager as tm

class SetOfHomogeneousNumericVariables:
    def __init__(self, variable_list, timber_variables):

        self.variable_list = variable_list

        if type(timber_variables) is str:
            dict_timber = tm.parse_timber_file(timber_variables, verbose=False)
        elif hasattr(timber_variables, '__getitem__'):
            dict_timber = timber_variables

        self.timber_variables = {}
        for var_name in self.variable_list:
             self.timber_variables[var_name] = dict_timber[var_name]
             self.timber_variables[var_name].t_stamps = np.array(self.timber_variables[var_name].t_stamps,dtype=float)

             self.timber_variables[var_name].values = np.array(self.timber_variables[var_name].values, dtype=float)
             shape = self.timber_variables[var_name].values.shape
             if len(shape) > 1 and shape[1] == 1:
                 self.timber_variables[var_name].values = np.squeeze(self.timber_variables[var_name].values, axis=1)

    def aligned(self):
        aligned_list = []
        t_first = self.timber_variables[self.variable_list[0]].t_stamps
        for ii,kk in enumerate(self.variable_list):
            values = self.timber_variables[kk].values
            t_stamps = self.timber_variables[kk].t_stamps
            if len(values) == 0:
                aligned_list.append(np.zeros_like(t_first))
            else:
                aligned_list.append(np.interp(t_first, t_stamps, values))
        return t_first, np.array(aligned_list)

    def aligned_object(self, dt_seconds):
        # Time stamps
        t_min = np.inf
        t_max = -np.inf
        empty_vars = []
        for key, var in self.timber_variables.iteritems():
            try:
                t_min = min(t_min, var.t_stamps[0])
                t_max = max(t_max, var.t_stamps[-1])
            except IndexError:
                empty_vars.append(key)
        print('Variables without data:', empty_vars)
        tt = np.arange(t_min, t_max, dt_seconds)

        # Data
        aligned_data = np.zeros(shape=(len(tt), len(self.variable_list)))
        for ii, kk in enumerate(self.variable_list):
            values = self.timber_variables[kk].values
            t_stamps = self.timber_variables[kk].t_stamps
            if len(values) != 0:
                # This assumes smooth data, otherwise the nearest older sample should be used!
                aligned_data[:,ii] = np.interp(tt, t_stamps, values)
        return tm.AlignedTimberData(tt, aligned_data, self.variable_list)

    def mean(self):
        t_first, aligned = self.aligned()
        return t_first, np.mean(aligned, axis=0)

    def std(self):
        t_first, aligned = self.aligned()
        return t_first, np.std(aligned, axis=0)

    def average(self, weights=None, returned=False):
        t_first, aligned = self.aligned()
        average, weights = np.average(aligned, axis=0, weights=weights, returned=returned)
        return t_first, average, weights

    def correct_values(self, correction_factors):
        for ii,var_name in enumerate(self.variable_list):
            self.timber_variables[var_name].values = self.timber_variables[var_name].values*correction_factors[ii]
