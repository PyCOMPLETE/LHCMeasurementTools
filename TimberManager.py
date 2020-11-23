import time
import numpy as np
import calendar
import os
import six

import h5py

timb_timestamp2float= lambda time_string: time.mktime(time.strptime(time_string,'%Y-%m-%d %H:%M:%S'))
timb_timestamp2float_ms= lambda time_string: time.mktime(time.strptime(time_string,'%Y-%m-%d %H:%M:%S.%f'))
timb_timestamp2float_UTC= lambda time_string: calendar.timegm(time.strptime(time_string,'%Y-%m-%d %H:%M:%S'))

def UnixTimeStamp2UTCTimberTimeString(t):
    return time.strftime('%Y-%m-%d %H:%M:%S.000', time.gmtime(t))

class timber_data_line(object):
    '''
    Parses a single line of a Timber csv file.
    '''
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

class CalsVariable(object):
    '''
    Object containing a single CALS variable with multiple time stamps.
    '''
    def __init__(self, t_stamps=[], values=[], ms=[]):
        self.t_stamps = t_stamps
        self.ms = ms
        self.values = values

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

    raise ValueError('Function suppressed, you can simply use CalsVariable')
    # if ms == None:
    #     ms = np.zeros_like(t_stamps)
    # assert len(t_stamps) == len(values) == len(ms)
    # tvl = timber_variable_list()
    # tvl.ms = ms
    # tvl.t_stamps = t_stamps
    # tvl.values = values
    # return tvl

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
            variables[vname] = CalsVariable()
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

def CalsVariables_to_h5(data, filename, varlist=None):

    if varlist is not None:
        varnames = varlist
    else:
        varnames = data.keys()

    dict_to_h5 = {}

    for varname in varnames:
        values = data[varname].values
        np_vals = list(map(np.atleast_1d, values))

        lngts = list(map(len, np_vals))
        if len(lngts) > 0:
            minlen = np.min(lngts)
            maxlen = np.max(lngts)
        else:
            minlen = 0
            maxlen = 0


        if minlen < maxlen:
            n_entries = len(data[varname].t_stamps)
            vals_for_h5 = np.zeros((n_entries, maxlen))
            vals_for_h5[:,:] = np.nan
            for ii in range(n_entries):
                vals_for_h5[ii, :len(np_vals[ii])] = np_vals[ii]
            #vals_for_h5 = np.concatenate(np_vals)
        else:
            vals_for_h5 = np.array(np_vals)

        dict_to_h5[varname+'!t_stamps'] = np.atleast_1d(
                    np.float_(data[varname].t_stamps))
        dict_to_h5[varname+'!values'] = vals_for_h5

    with h5py.File(filename, 'w') as fid:
        for kk in list(dict_to_h5.keys()):
            #fid[kk] = dict_to_h5[kk]
            fid.create_dataset(kk, data=dict_to_h5[kk],
                    compression='lzf')
                    # 'lzf' filter used to have good speed
                    # https://docs.h5py.org/en/stable/high/dataset.html#lossless-compression-filters

def CalsVariables_from_h5(filename, remove_nans=True):
    dict_data = {}
    with h5py.File(filename, 'r') as fid:
        for kk in list(fid.keys()):
            varname = kk.split('!')[0]
            part = kk.split('!')[1]
            if varname not in dict_data:
                dict_data[str(varname)] = CalsVariable()
            if part=='t_stamps':
                dict_data[varname].t_stamps =  np.atleast_1d(fid[kk][:])
            elif part=='values':
                dict_data[varname].values =  list(np.atleast_2d(fid[kk][:]))

    if remove_nans:
        for kk in dict_data.keys():
            for ii, vv in enumerate(dict_data[kk].values):
                mask = ~np.isnan(vv)
                if np.any(mask):
                    dict_data[kk].values[ii] = vv[mask]

    return dict_data

def CalsVariables_from_pytimber(pt_variables):
    dict_data = {}
    for varname in list(pt_variables.keys()):
        dict_data[varname] = CalsVariable()
        dict_data[varname].t_stamps = np.atleast_1d(pt_variables[varname][0])
        dict_data[varname].values = list(map(np.atleast_1d, pt_variables[varname][1]))
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


# For backward compatibility --> To be removed!
timber_variable_list = CalsVariable
timber_variables_from_h5 = CalsVariables_from_h5

try:

    from pytimber.nxcals import NXCals

    # STILL TO BE COMPLETED, DOES NOT WORK WITH UNIX TIMESTAMPS
    class NXCalsFastQuery(NXCals):
        '''
        This class can replace pytimber to have fast extraction of many
        scalar variables.
        '''
        def __init__(self, *args, **kwargs):
            if 'system' in kwargs.keys():
                self.system = kwargs['system']
                kwargs.pop('system')
            super().__init__(*args, **kwargs)

        def toTimestring(self, t):
            if isinstance(t, six.string_types):
               return t
            else: #We assume linux timestamp
               return time.strftime("%Y-%m-%d %H:%M:%S",
                                            time.gmtime(t))+'.000'

        def get(self, variables, t1, t2, system=None):
            '''
                system should be either 'CMW' or 'WINCCOA'
            '''

            if system is None:
                system = self.system

            assert(system is not None)

            query = self.DataQuery.byVariables()\
                .system(system)\
                .startTime(self.toTimestring(t1))\
                .endTime(self.toTimestring(t2))

            for vv in variables:
                query = query.variable(vv)

            dfp = query.build()\
                    .sort("nxcals_variable_name","nxcals_timestamp")\
                    .na().drop()\
                    .select("nxcals_timestamp",
                            "nxcals_value", "nxcals_variable_name")

            data1=np.fromiter(
                    (tuple(dd.values()) for dd in dfp.collect()),
                    dtype=[('ts',int),('val',float),('var','U32')] )

            out={}
            for var in set(data1['var']):
              sel=data1['var']==var
              out[var]=(data1['ts'][sel]/1e9,data1['val'][sel])

            for vv in variables:
                if vv not in out.keys():
                    print(f'{vv} not found!')
                    out[vv] = (
                            np.array([], dtype=np.float64),
                            np.array([], dtype=np.float64))

            return out

except Exception:
    print('NXCALS not possible!!!')
