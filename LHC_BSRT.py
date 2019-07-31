import numpy as np
import TimberManager as tm

na = np.array


class BSRT:
    def __init__(self, timber_variable_bsrt, beam=0, calib_dict=None,
                 average_repeated_meas=False, filter_FESA_from=False):

        # assume timber_variable_bsrt is filename string for now
        if not (beam == 1 or beam == 2):
            raise ValueError('You need to specify which beam! (1 or 2)')

        self.beam = beam
        self.e_dict = calib_dict
        if hasattr(timber_variable_bsrt, '__getitem__'):
            dict_timber = timber_variable_bsrt
        else:
            dict_timber = tm.parse_timber_file(timber_variable_bsrt, verbose=True)

        sigma_h = dict_timber[get_variable_dict(beam)['SIGMA_H']]
        sigma_v = dict_timber[get_variable_dict(beam)['SIGMA_V']]
        gate = dict_timber[get_variable_dict(beam)['GATE_DELAY']]

        if (np.sum(np.abs(na(sigma_h.t_stamps) - na(sigma_v.t_stamps)))>1e-6 or
            np.sum(np.abs(na(sigma_h.t_stamps) - na(gate.t_stamps)))>1e-6):
            raise ValueError('Timestamps for the two channels (H and V) not equal!')
        #~ if (sigma_v.t_stamps != gate.t_stamps):
            #~ bunches_effectively_recorded = self.get_bunches_effectively_recorded(gate, sigma_v)
        else:
            bunches_effectively_recorded = gate

        def filter_FESA_bug(data, indx):
            for d in data:
                for v in d.values:
                    del v[np.s_[indx:]]
        if filter_FESA_from:
            filter_FESA_bug((sigma_h, sigma_v, gate), filter_FESA_from)

        self.t_stamps = []
        self.bunch_n = []
        self.sigma_h = []
        self.sigma_v = []

        for ii in xrange(len(bunches_effectively_recorded.t_stamps)):
            if np.mod(ii,10000) == 0:
                print 'expanding %.1f' % (
                    float(ii)/len(bunches_effectively_recorded.t_stamps)*100) + """%"""

            N_meas = len(bunches_effectively_recorded.values[ii])

            if average_repeated_meas:
                bunch_curr_aver = int(float(bunches_effectively_recorded.values[ii][0]))
                n_aver = 0
                sigma_h_sum = 0.
                sigma_v_sum = 0.
                for jj in xrange(N_meas):
                    sigma_h_sum+= float(sigma_h.values[ii][jj])
                    sigma_v_sum+= float(sigma_v.values[ii][jj])
                    n_aver+=1

                    if jj == N_meas-1:
                        self.bunch_n.append(bunch_curr_aver)
                        self.sigma_h.append(sigma_h_sum/float(n_aver))
                        self.sigma_v.append(sigma_v_sum/float(n_aver))
                        self.t_stamps.append(np.float_(sigma_v.t_stamps[ii]) + float(jj)*1e-3/float(N_meas))
                    elif int(float(bunches_effectively_recorded.values[ii][jj+1]))!=bunch_curr_aver:
                        self.bunch_n.append(bunch_curr_aver)
                        self.sigma_h.append(sigma_h_sum/float(n_aver))
                        self.sigma_v.append(sigma_v_sum/float(n_aver))
                        self.t_stamps.append(np.float_(sigma_v.t_stamps[ii]) + float(jj)*1e-3/float(N_meas))

                        bunch_curr_aver = int(float(bunches_effectively_recorded.values[ii][jj+1]))
                        n_aver = 0
                        sigma_h_sum = 0.
                        sigma_v_sum = 0.
            else:
                for jj in xrange(N_meas):
                    self.bunch_n.append(np.float_(bunches_effectively_recorded.values[ii][jj]))
                    self.sigma_h.append(np.float_(sigma_h.values[ii][jj]))
                    self.sigma_v.append(np.float_(sigma_v.values[ii][jj]))
                    self.t_stamps.append(np.float_(sigma_v.t_stamps[ii]) + float(jj)*1e-3/float(N_meas))

        self.t_stamps = np.array(self.t_stamps)
        self.bunch_n = np.array(self.bunch_n)
        self.sigma_h = np.array(self.sigma_h)
        self.sigma_v = np.array(self.sigma_v)


    # def calculate_emittances_slow(self, energy_ob):
    #     if self.e_dict is None:
    #         e_dict = emittance_dictionary()
    #     else:
    #         e_dict = self.e_dict
    #     self.norm_emit_h = []
    #     self.norm_emit_v = []
    #     for ii in xrange(len(self.t_stamps)):
    #         if np.mod(ii,10000)==0:
    #             print 'calc. emitt. %.1f'%(float(ii)/len(self.t_stamps)*100)+"""%"""

    #         norm_emit_h  = 0.
    #         norm_emit_v  = 0.

    #         energy = energy_ob.nearest_older_sample(self.t_stamps[ii])
    #         if energy > 400. and energy < 500.:
    #             energy = 450.
    #         elif energy > 6400. and energy < 6600.:
    #             energy = 6500.
    #         else:
    #             self.norm_emit_h.append(norm_emit_h)
    #             self.norm_emit_v.append(norm_emit_v)
    #             continue

    #         sigma_h_corr_sq = self.sigma_h[ii]**2 - e_dict['sigma_corr_h'][energy][self.beam]**2
    #         sigma_v_corr_sq = self.sigma_v[ii]**2 - e_dict['sigma_corr_v'][energy][self.beam]**2

    #         phys_emit_h = sigma_h_corr_sq/e_dict['betaf_h'][energy][self.beam]
    #         phys_emit_v = sigma_v_corr_sq/e_dict['betaf_v'][energy][self.beam]

    #         norm_emit_h  = phys_emit_h*e_dict['gamma'][energy]
    #         norm_emit_v  = phys_emit_v*e_dict['gamma'][energy]

    #         self.norm_emit_h.append(norm_emit_h)
    #         self.norm_emit_v.append(norm_emit_v)

    #     self.norm_emit_h = np.array(self.norm_emit_h)
    #     self.norm_emit_v = np.array(self.norm_emit_v)


    def calculate_emittances(self, energy_ob):
        if self.e_dict is None:
            e_dict = emittance_dictionary()
        else:
            e_dict = self.e_dict
        self.norm_emit_h = []
        self.norm_emit_v = []

        sigma_corr_h = 0.*self.sigma_h
        sigma_corr_v = 0.*self.sigma_v
        rescale_sigma_h = 0.*self.sigma_h
        rescale_sigma_v = 0.*self.sigma_v
        betaf_h = 0.*self.sigma_h+1.
        betaf_v = 0.*self.sigma_v+1.
        gamma = 0.*self.sigma_h

        energy = np.array(map(lambda x: energy_ob.nearest_older_sample(x), self.t_stamps))

        mask_450 = np.logical_and(energy > 400., energy < 500.)
        sigma_corr_h[mask_450] = e_dict['sigma_corr_h'][450.][self.beam]
        sigma_corr_v[mask_450] = e_dict['sigma_corr_v'][450.][self.beam]
        rescale_sigma_h[mask_450] = e_dict['rescale_sigma_h'][450.][self.beam]
        rescale_sigma_v[mask_450] = e_dict['rescale_sigma_v'][450.][self.beam]
        betaf_h[mask_450] = e_dict['betaf_h'][450.][self.beam]
        betaf_v[mask_450] = e_dict['betaf_v'][450.][self.beam]
        gamma[mask_450] = e_dict['gamma'][450.]

        mask_6500 = np.logical_and(energy > 6400., energy < 6600.)
        sigma_corr_h[mask_6500] = e_dict['sigma_corr_h'][6500.][self.beam]
        sigma_corr_v[mask_6500] = e_dict['sigma_corr_v'][6500.][self.beam]
        rescale_sigma_h[mask_6500] = e_dict['rescale_sigma_h'][6500.][self.beam]
        rescale_sigma_v[mask_6500] = e_dict['rescale_sigma_v'][6500.][self.beam]
        betaf_h[mask_6500] = e_dict['betaf_h'][6500.][self.beam]
        betaf_v[mask_6500] = e_dict['betaf_v'][6500.][self.beam]
        gamma[mask_6500] = e_dict['gamma'][6500.]

        sigma_h_corr_sq = (self.sigma_h*rescale_sigma_h)**2 - sigma_corr_h**2
        sigma_v_corr_sq = (self.sigma_v*rescale_sigma_v)**2 - sigma_corr_v**2

        phys_emit_h = sigma_h_corr_sq/betaf_h
        phys_emit_v = sigma_v_corr_sq/betaf_v

        self.norm_emit_h = phys_emit_h*gamma
        self.norm_emit_v = phys_emit_v*gamma


    def find_start_scans(self, scan_thresh):
        diff_bunch = np.diff(self.bunch_n)
        ind_start_scan_all = np.where(diff_bunch < -scan_thresh)[0]
        ind_start_scan = ind_start_scan_all[:-1][np.diff(ind_start_scan_all) > 10]
        self.t_start_scans = self.t_stamps[ind_start_scan]
        self.t_start_scans = np.array(list(self.t_start_scans)+[self.t_stamps[-1]])

        # return self.t_start_scans


    def find_closest_scan(self, t_start_requested, scan_thresh):
        self.find_start_scans(scan_thresh)
        ind_closest_scan = np.argmin(np.abs(t_start_requested - self.t_start_scans))
        t_start = self.t_start_scans[ind_closest_scan]
        if ind_closest_scan + 1 >= len(self.t_start_scans):
            raise IndexError('Index ind_closest_scan + 1 is out of bounds.\nYour requested scan times might be outside of the fill.')
        t_stop = self.t_start_scans[ind_closest_scan + 1]
        return Masked(self, t_start, t_stop)


    def get_bunches_effectively_recorded(self, gate_timber, sigma_timber):
        recorded_bunches = tm.timber_variable_list()
        i1 = 0
        for i2 in xrange(len(gate_timber.t_stamps)):
            if np.mod(i2,10000) == 0:
                print 'Cleaning %.1f'%(float(i2)/len(gate_timber.t_stamps)*100)+"""%"""
            if gate_timber.t_stamps[i2] == sigma_timber.t_stamps[i1]:
                recorded_bunches.t_stamps.append(gate_timber.t_stamps[i2])
                recorded_bunches.values.append(gate_timber.values[i2])
                i1 += 1

        return recorded_bunches


    def get_bbb_emit_evolution(self):
        bunch_n_un= np.sort(np.unique(self.bunch_n))
        emit_h_bbb=[]
        emit_v_bbb=[]
        t_bbb=[]

        dict_bunches = {}

        for i_line, i_bunch in enumerate(bunch_n_un):
                if np.mod(i_line,10)==0:
                    print 're-shuffle %.1f'%(float(i_line)/len(bunch_n_un)*100)+"""%"""
                inds=np.nonzero(self.bunch_n==i_bunch)
                x=self.norm_emit_h[inds]
                y=self.norm_emit_v[inds]
                t=self.t_stamps[inds]

                emit_h_bbb.append(x)
                emit_v_bbb.append(y)
                t_bbb.append(t)

                dict_bunches[i_bunch] = {}
                dict_bunches[i_bunch]['norm_emit_h'] = x
                dict_bunches[i_bunch]['norm_emit_v'] = y
                dict_bunches[i_bunch]['t_stamp'] = t


        return dict_bunches, t_bbb, emit_h_bbb, emit_v_bbb, bunch_n_un


class Masked:
    def __init__(self, bsrt, t_start, t_stop):
        self.t_start = t_start
        self.t_stop = t_stop
        #mask_bsrt = np.logical_and(bsrt.t_stamps >= self.t_start, bsrt.t_stamps < self.t_stop)
        mask_bsrt = np.logical_and(bsrt.t_stamps > self.t_start, bsrt.t_stamps <= self.t_stop)

        self.beam = bsrt.beam
        self.t_stamps = bsrt.t_stamps[mask_bsrt]
        self.bunch_n = bsrt.bunch_n[mask_bsrt]
        self.sigma_h = bsrt.sigma_h[mask_bsrt]
        self.sigma_v = bsrt.sigma_v[mask_bsrt]
        if hasattr(bsrt, 'norm_emit_h'):
            self.norm_emit_h = bsrt.norm_emit_h[mask_bsrt]
            self.norm_emit_v = bsrt.norm_emit_v[mask_bsrt]



def emittance_dictionary():
    e_dict = {'betaf_h':{}, 'betaf_v':{}, 'gamma':{},
              'sigma_corr_h':{}, 'sigma_corr_v':{}}
    e_dict['betaf_h'][450] = {1:205.5, 2:191.5}
    e_dict['betaf_v'][450] = {1:320., 2:387.8}
    e_dict['betaf_h'][6500] = {1:204.1, 2:191.5}
    e_dict['betaf_v'][6500] = {1:322.7, 2:395.}
    e_dict['gamma'][450] = 479.6
    e_dict['gamma'][6500] = 6927.6
    e_dict['sigma_corr_h'][450] = {1:0.,2:0.}#0.85
    e_dict['sigma_corr_v'][450] = {1:0., 2:0.}#0.87
    e_dict['sigma_corr_h'][6500] = {1:0.32, 2:0.34}#0.2 #0.35
    e_dict['sigma_corr_v'][6500] = {1:0.23, 2:0.28}#0.#0.2 #0.33

    return e_dict


def get_variable_dict(beam):
    beam_device_list = ['R','L']
    var_dict = {}
    var_dict['GATE_DELAY'] = 'LHC.BSRT.5%s4.B%d:GATE_DELAY'%(beam_device_list[beam-1],beam)
    var_dict['SIGMA_H'] = 'LHC.BSRT.5%s4.B%d:FIT_SIGMA_H'%(beam_device_list[beam-1],beam)
    var_dict['SIGMA_V'] = 'LHC.BSRT.5%s4.B%d:FIT_SIGMA_V'%(beam_device_list[beam-1],beam)

    return var_dict


def variable_list(beams = [1,2]):
    var_list = []
    for beam in beams:
                var_list += get_variable_dict(beam).values()

    return var_list
