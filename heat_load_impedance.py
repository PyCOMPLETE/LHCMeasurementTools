# All formula from Elias Metral presentation 2010 'Beam Screen Issues'
import os
import numpy as np
from scipy.constants import e as qe, c, m_p, epsilon_0
from scipy.interpolate import interp1d
from scipy.special import gamma

lhc_circumference = 26658.883 # design report
lhc_bending_radius = 2803.95 # design report
default_chamb_radius = 18.4e-3 # Elias
default_temperature = 20
default_weld_thickness = 1./60.* 2.*np.pi*default_chamb_radius # Elias
rho_stainless_steel = 6e-7 # Elias
Z0 = 376.73031
lhc_arc_length = 171.4*2 + 23*2*53.45

lhc_measurement_tools_path = os.path.dirname(os.path.realpath(__file__))
default_rho_vs_T_file = lhc_measurement_tools_path + '/copper_rho_vs_T.txt'

class HeatLoadImpedance(object):
    def __init__(self, rho_vs_T_file = default_rho_vs_T_file):
        # load resistivity
        data = np.loadtxt(rho_vs_T_file)
        self.copper_rho_Ohm_m = interp1d(data[:,0], data[:,1]*1e-8)
        
        
    def P_RW_Wm_1beam(self, ppb, sigma_t, temp, circ, chamb_radius, b_field=None,
             weld_Thickness=None, weld_Rho=None, n_bunches=None):
        """
        If n_bunches is not specified, ppb has to be either a numpy array with dimensions 
        time_steps x n_bunches or n_bunches.
        If n_bunches is specified, ppb can be either a number or an arry with dimensions of timesteps
        sigma_t is either a float or a np array of dimensions time_steps or timp_steps * n_bunches
        Do not specify n_bunches explicitly and implicitly at the same time!
        """
        
        if n_bunches is None:
            if type(ppb) is float or type(ppb) is int:
                raise ValueError('N_bunches has to be specified somehow!')
            else:
                n_bunches = 1.

        rho_B0_T = self.copper_rho_Ohm_m(temp) # 0.014 *1e-8 for 20 K
        
        if b_field is None:
            rho_B_T = rho_B0_T
        else:
            rho_B0_273 = self.copper_rho_Ohm_m(273.)
            try:
                rho_B0_4 = self.copper_rho_Ohm_m(4.)
            except ValueError:
                rho_B0_4 = 3.*self.copper_rho_Ohm_m(4.2) - 2.*self.copper_rho_Ohm_m(4.3)
            rho_B_T = rho_B0_T * (1. + 10.**(-2.65) * (b_field*rho_B0_273/rho_B0_4)**1.055)

        per_bunch_factor = ppb**2 * sigma_t**(-1.5)
        try:
            if len(ppb.shape) == 2:
                per_bunch_factor = np.sum(per_bunch_factor, axis=1)
        except AttributeError: 
            # ppb is not a numpy array
            pass

        P_no_weld = 1./circ*gamma(0.75)*n_bunches/chamb_radius * (qe/2./np.pi)**2* np.sqrt(c*rho_B_T*Z0/2.) * per_bunch_factor

        if weld_Thickness is not None:
            weld_Factor = 1. + np.sqrt(weld_Rho/rho_B_T)*weld_Thickness/(2.*np.pi*chamb_radius)
        else:
            weld_Factor = 1.

        return P_no_weld*weld_Factor

    def P_RW_Wm_1beam_fill_arc_half_cell(self, fill_dict, temp=default_temperature, chamb_radius=default_chamb_radius):
        """
        Returns the half cell heat load for a fill dict, which has to consist of basic and bunchbybunch data
        """
        from LHCMeasurementTools.LHC_BCT import BCT
        from LHCMeasurementTools.LHC_FBCT import FBCT
        from LHCMeasurementTools.LHC_BQM import filled_buckets, blength
        from LHCMeasurementTools.LHC_Heatloads import magnet_length
        from LHCMeasurementTools.LHC_Energy import energy

        l_halfcell = magnet_length['AVG_ARC'][0]
        l_dip_halfcell = 3. * magnet_length['special_HC_D2'][0]
 
        bunch_length_seconds = {}
        bunch_charge = {}

        # Different for both beams
        for beam_ctr in (1,2):

            bct = BCT(fill_dict, beam=beam_ctr)
            fbct = FBCT(fill_dict, beam_ctr)

            fbct_correction_factors = np.ones_like(fbct.bint)
            fbct_bunch_length = np.zeros_like(fbct.bint)
            bunch_length = blength(fill_dict, beam = beam_ctr)

            # rescale bct and bunch length to fbct time stamps
            for index_fbct in xrange(len(fbct.totint)):
                index_bunch_length = np.abs(bunch_length.t_stamps-fbct.t_stamps[index_fbct]).argmin()
                fbct_bunch_length[index_fbct,:] = bunch_length.blen[index_bunch_length,:]

                if fbct.totint[index_fbct] != 0:
                    index_bct = np.abs(bct.t_stamps - fbct.t_stamps[index_fbct]).argmin()
                    fbct_correction_factors[index_fbct,:] = bct.values[index_bct]/fbct.totint[index_fbct]
                    
            # correct empty bunches to avoid divisions by 0
            mask_bunch_length = fbct_bunch_length < 1e-15
            fbct_bunch_length[mask_bunch_length] = 1.
            fbct.bint[mask_bunch_length] = 0.

            bunch_length_seconds[beam_ctr] = fbct_bunch_length/4. # 2 sigmas*2 are stored. 1 sigma is required
            bunch_charge[beam_ctr] = fbct.bint * fbct_correction_factors

        # Identical for both beams
        beam_energy_GeV = np.zeros_like(fbct.totint)
        energy_ob = energy(fill_dict, beam=1)

        for index_fbct in xrange(len(fbct.totint)):
            index_energy = np.abs(energy_ob.t_stamps - fbct.t_stamps[index_fbct]).argmin()
            beam_energy_GeV[index_fbct] = energy_ob.energy[index_energy]
        dipole_B_Field = beam_energy_GeV / 7000. * 8.33

        # Add time steps to allow for easier plots
        self.t_stamps = fbct.t_stamps

        # Make one long beam out of the two LHC beams
        bunch_charge_both_beams = np.concatenate((bunch_charge[1], bunch_charge[2]), axis=1)
        bunch_length_seconds_both_beams = np.concatenate((bunch_length_seconds[1], bunch_length_seconds[2]), axis=1)

        # Lambda function to minimize copy pasting of code
        standard_lhc_p_rw_wm_1beam = lambda b_field: \
                self.P_RW_Wm_1beam(bunch_charge_both_beams, bunch_length_seconds_both_beams, temp, lhc_circumference, default_chamb_radius,
                b_field=b_field, weld_Thickness=default_weld_thickness, weld_Rho=rho_stainless_steel)

        p_dipole = standard_lhc_p_rw_wm_1beam(dipole_B_Field)

        p_drift = standard_lhc_p_rw_wm_1beam(None)

        return p_dipole*l_dip_halfcell + p_drift*(l_halfcell-l_dip_halfcell)




class HeatLoadSynchRad(object):
    def __init__(self, fill_dict=None, tot_int=None, bend_rad=None, energy_GeV=None, circ=None):
        if fill_dict is not None:
            self.init_fill(fill_dict)
        else:
            self.power_bend_m = synch_rad_power_bend_m(tot_int, bend_rad, energy_GeV, circ)

    def init_fill(self, fill_dict):
        from LHCMeasurementTools.LHC_BCT import BCT
        from LHCMeasurementTools.LHC_Energy import energy
        
        total_intensity = {}
        for beam_ctr in (1,2):
            bct = BCT(fill_dict, beam=beam_ctr)
            total_intensity[beam_ctr] = bct.values

        total_intensity_both_beams = total_intensity[1] + total_intensity[2]

        energy_ob = energy(fill_dict, beam=beam_ctr)
        energy_beam_GeV = energy_ob.interp(bct.t_stamps)
        
        self.t_stamps = bct.t_stamps
        self.power_bend_m = synch_rad_power_bend_m(total_intensity_both_beams, lhc_bending_radius, energy_beam_GeV, lhc_circumference)
        self.avg_arc_heat_load_m =  self.power_bend_m * 2*np.pi*lhc_bending_radius / 8 / lhc_arc_length
    
def synch_rad_power_bend_m(tot_int, bend_rad, energy_GeV, circ):
    beam_gamma = energy_GeV * qe * 1e9 / (m_p*c**2)
    energy_loss_turn_beam = (qe**2 * beam_gamma**4) / (3. * epsilon_0*bend_rad) * tot_int
    energy_loss_bend_m = energy_loss_turn_beam / (2*np.pi*bend_rad)
    power_bend_m = energy_loss_bend_m * (c/circ)

    return power_bend_m



    

            


