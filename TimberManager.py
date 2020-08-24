import time
import numpy as np
import calendar
import os

timb_timestamp2float= lambda time_string: time.mktime(time.strptime(time_string,'%Y-%m-%d %H:%M:%S'))
timb_timestamp2float_ms= lambda time_string: time.mktime(time.strptime(time_string,'%Y-%m-%d %H:%M:%S.%f'))
timb_timestamp2float_UTC= lambda time_string: calendar.timegm(time.strptime(time_string,'%Y-%m-%d %H:%M:%S'))

def UnixTimeStamp2UTCTimberTimeString(t):
    return time.strftime('%Y-%m-%d %H:%M:%S.000', time.gmtime(t))

class timber_data_line(object):
    def __init__(self, rline, time_input_UTC = False):
        list_values = rline.split(',')
        t_string = list_values[0]
        if time_input_UTC:
            self.timestamp = timb_timestamp2float_UTC(t_string.split('.')[0])
            self.ms = float(t_string.split('.')[-1])
        else:
            self.timestamp = timb_timestamp2float(t_string.split('.')[0])
            self.ms = float(t_string.split('.')[-1])
        self.data_strings = list_values[1:]

class timber_variable_list(object):
    def __init__(self):
        self.t_stamps = []
        self.ms = []
        self.values = []

    def float_values(self):
        return np.squeeze(np.float_(self.values))

    def nearest_older_sample(self, t_stamp):
        index = np.argmin(np.abs(np.array(self.t_stamps) - t_stamp))
        if self.t_stamps[index] > t_stamp:
            index -= 1
        return self.values[index]

    def nearest_t_stamp(self, value):
        index = np.argmin(np.abs(np.array(self.values) - value))
        return self.t_stamps[index]

    def calc_avg(self, begin, end):
        return np.mean(self.selection(begin, end))

    def selection(self, begin, end):
        mask = np.logical_and(self.t_stamps > begin, self.t_stamps < end)
        try:
            return self.values[mask]
        except:
            print((self.values.shape, self.t_stamps.shape))
            raise
    def interp(self, t_stamps):
        return np.interp(t_stamps, self.t_stamps, self.values)


def make_timber_variable_list(t_stamps, values, ms=None):
    if ms == None:
        ms = np.zeros_like(t_stamps)
    assert len(t_stamps) == len(values) == len(ms)
    tvl = timber_variable_list()
    tvl.ms = ms
    tvl.t_stamps = t_stamps
    tvl.values = values
    return tvl

def parse_timber_file(timber_filename, verbose=True):

    with open(timber_filename) as fid:
        timber_lines = fid.readlines()

    time_input_UTC = False

    N_lines = len(timber_lines)

    i_ln = 0

    variables = {}
    while i_ln < N_lines:
        line = timber_lines[i_ln]
        line = line.split('\n')[0]
        i_ln = i_ln + 1

        if 'VARIABLE:' in line:
            vname = line.split(': ')[-1]
            if verbose:
                print('\n\nStarting variable: ' + vname)
            variables[vname] = timber_variable_list()
        else:
            try:
                currline_obj = timber_data_line(line, time_input_UTC=time_input_UTC)
                if currline_obj.data_strings == ['']:
                    raise ValueError
                variables[vname].t_stamps.append(currline_obj.timestamp)
                variables[vname].values.append(currline_obj.data_strings)
                variables[vname].ms.append(currline_obj.ms)
            except ValueError:
                if 'Timestamp (UTC_TIME)' in line:
                    time_input_UTC = True
                    if verbose: print('Set time to UTC')
                if verbose:
                    print('Skipped line: '+    line)
    return variables

def timber_variables_from_h5(filename):
    import h5py
    dict_data = {}
    with h5py.File(filename, 'r') as fid:
        for kk in list(fid.keys()):
            varname = kk.split('!')[0]
            part = kk.split('!')[1]
            if varname not in dict_data:
                dict_data[str(varname)] = timber_variable_list()
            if part=='t_stamps':
                dict_data[varname].t_stamps =  np.atleast_1d(fid[kk][:])
            elif part=='values':
                dict_data[varname].values =  list(np.atleast_2d(fid[kk][:]))
    return dict_data

def dbquery(varlist, t_start, t_stop, filename):

    if type(t_start) is not str:
        t_start_str_UTC = UnixTimeStamp2UTCTimberTimeString(t_start)
    else:
        t_start_str_UTC = t_start

    if type(t_stop) is not str:
        t_stop_str_UTC = UnixTimeStamp2UTCTimberTimeString(t_stop)
    else:
        t_stop_str_UTC = t_stop

    if type(varlist) is not list:
        raise TypeError

    varlist_str = ''
    for var in varlist:
        varlist_str += var +','
    varlist_str = varlist_str[:-1]

    execut = 'java -jar accsoft-cals-extractor-client-nodep.jar '
    config = ' -C ldb_UTC.conf '
    time_interval = ' -t1 "'+ t_start_str_UTC +'" -t2 "'+ t_stop_str_UTC +'"'
    variables = '-vs "%s"'%(varlist_str)
    outpfile = ' -N .//' + filename

    command = execut + config + variables + time_interval + outpfile
    print(command)
    os.system(command)

class AlignedTimberData(object):

    def __init__(self, timestamps, data, variables):
        dictionary = {}
        double = []
        for ctr, var in enumerate(variables):
            if var in dictionary:
                double.append(var)
                var += '_2'
            dictionary[var] = data[:,ctr]
        if double:
            print(('Duplicate variables: %s' % double))

        self.timestamps = timestamps
        self.data = data
        self.variables = variables
        self.dictionary = dictionary

    def nearest_older_index(self, t):
        index = np.argmin(np.abs(self.timestamps - t))
        if self.timestamps[index] > t:
            index -= 1
        return index

    def nearest_older_sample(self, t):
        return self.data[self.nearest_older_index(t),:]

def atd_from_h5(obj):
    return AlignedTimberData(obj.timestamps, obj.data, obj.variables)

def parse_aligned_csv_file(filename, empty_value=np.nan):
    timestamps = []
    np_data = []
    with open(filename,'r') as csv_file:
        lines = csv_file.read().splitlines()

    for line_n, line in enumerate(lines):
        csv = line.split(',')
        if line_n == 0:
            variables = csv[1:]
        else:
            timestamps.append(timb_timestamp2float_ms(csv[0]))
            values = []
            for string in csv[1:]:
                try:
                    values.append(float(string))
                except ValueError:
                    values.append(empty_value)
            np_data.append(values)

    np_data = np.array(np_data)
    timestamps = np.array(timestamps)
    return AlignedTimberData(timestamps, np_data, variables)
