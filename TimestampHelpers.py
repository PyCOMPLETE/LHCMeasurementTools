import time
import numpy as np
import datetime
import matplotlib


def unixstamp2datestring(t_stamp):
    return time.strftime("%Y_%m_%d", time.localtime(t_stamp))

def unixstamp_of_midnight_before(t_stamp):
    day_string = unixstamp2datestring(t_stamp)
    return time.mktime(time.strptime(day_string+' 00:00:00', '%Y_%m_%d %H:%M:%S'))

def date_strings_interval(start_tstamp_unix, end_tstamp_unix):
    list_date_strings = []
    for t_day in np.arange(unixstamp_of_midnight_before(start_tstamp_unix), end_tstamp_unix+1., 24*3600.):
        list_date_strings.append(unixstamp2datestring(t_day))
    return list_date_strings

def localtime2unixstamp(local_time_str, strformat='%Y_%m_%d %H:%M:%S'):
    return time.mktime(time.strptime(local_time_str, strformat))

def unixstamp2localtime(t_stamp, strformat='%Y_%m_%d %H:%M:%S'):
    return time.strftime(strformat, time.localtime(t_stamp))

def unixstamp2localtimestamp(t_stamps):
    tlocal = map(datetime.datetime.fromtimestamp, t_stamps)
    return matplotlib.dates.date2num(tlocal)

def unixstampNow():
        return time.mktime(time.localtime())

class TimeConverter(object):
    def __init__(self, time_in, t_ref=0, t_plot_tick_h=None):
        self.time_in = time_in
        self.t_ref = t_ref
        self.t_plot_tick_h = t_plot_tick_h

    def from_unix(self, t_stamps):
        if self.time_in == 'h':
            ret = ((t_stamps-self.t_ref)/3600.)
        if self.time_in == 'd':
            ret = ((t_stamps-self.t_ref)/3600./24.)
        elif self.time_in == 'datetime' or self.time_in == 'hourtime':
            ret = matplotlib.dates.date2num(map(datetime.datetime.fromtimestamp, np.atleast_1d(t_stamps)))

        return ret

    def set_x_for_plot(self, fig, ax):
        if self.time_in == 'hourtime':
            hfmt = matplotlib.dates.DateFormatter('%H:%M')
            ax.xaxis.set_major_formatter(hfmt)
            if self.t_plot_tick_h is not None:
                ax.xaxis.set_major_locator(matplotlib.dates.HourLocator(np.arange(0, 24, self.t_plot_tick_h)))

        if self.time_in == 'datetime':
            hfmt = matplotlib.dates.DateFormatter('%a %d-%m %H:%M')
            ax.xaxis.set_major_formatter(hfmt)
            if self.t_plot_tick_h is not None:
                if self.t_plot_tick_h<24:
                    ax.xaxis.set_major_locator(matplotlib.dates.HourLocator(np.arange(0, 24, self.t_plot_tick_h)))
                elif self.t_plot_tick_h=='week':
                    ax.xaxis.set_major_locator(matplotlib.dates.WeekdayLocator(0))
                elif self.t_plot_tick_h=='2weeks':
                    ax.xaxis.set_major_locator(matplotlib.dates.WeekdayLocator(0, interval=2))
                elif self.t_plot_tick_h=='4weeks':
                    ax.xaxis.set_major_locator(matplotlib.dates.WeekdayLocator(0, interval=4))
            fig.autofmt_xdate()
