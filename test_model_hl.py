import sys, os
BIN = os.path.expanduser("../../../")
sys.path.append(BIN)

import numpy as np
from scipy.constants import e as const_e, c as const_c, m_p, epsilon_0
from scipy.special import gamma
from scipy.interpolate import interp1d

from LHCMeasurementTools.LHC_BCT import BCT
from LHCMeasurementTools.LHC_FBCT import FBCT
from LHCMeasurementTools.LHC_BQM import filled_buckets, blength
from LHCMeasurementTools.LHC_Energy import energy
from LHCMeasurementTools.LHC_Heatloads import magnet_length
import LHCMeasurementTools.TimberManager as tm

l_halfcell = magnet_length['AVG_ARC'][0]
l_dip_halfcell = 3. * magnet_length['special_HC_D2'][0]
aperture_radius = 18.4e-3 # Elias Metral presentation 2010 'Beam Screen Issues'
rho_stainless_steel_20K = 6e-7 # Elias 
weld_surface_factor = 1./60. # Elias

lhc_circumference = 26658.883 # design report
lhc_bending_radius = 2803.95  # design report
lhc_f_rev = const_c / lhc_circumference

const_Z0 = 376.73031
impedance_constant = gamma(0.75)/lhc_circumference * 1./aperture_radius * (const_e / (2.*np.pi))**2. * np.sqrt(const_c*const_Z0/2.)

class HeatLoadImpedanceSynchRad(object):

    def __init__(self, fill_dict, rho_vs_T_file='./copper_rho_vs_T.txt', temp_beam_screen_K=20, aperture_radius=aperture_radius):
        data = np.loadtxt(rho_vs_T_file)
        copper_rho_Ohm_m = interp1d(data[:,0], data[:,1]*1e-8)
        #self.filln = filln
        #
        #fill_dict = {}
        #fill_dict.update(tm.parse_timber_file('../fill_basic_data_csvs/basic_data_fill_%d.csv'%filln, verbose=False))
        #fill_dict.update(tm.parse_timber_file('../fill_bunchbybunch_data_csvs/bunchbybunch_data_fill_%d.csv'%filln, verbose=False))

        bunch_length_seconds = {}
        bunch_charge_squared = {}
        total_intensity_beams = {}

        for beam_ctr in (1,2):

            bct = BCT(fill_dict, beam=beam_ctr)
            fbct = FBCT(fill_dict, beam_ctr)

            bunch_int_sq = fbct.bint**2

            fbct_correction_factors = np.ones_like(fbct.bint)
            fbct_bunch_length = np.zeros_like(fbct.bint)
            bunch_length = blength(fill_dict, beam = beam_ctr)
            total_intensity = np.zeros_like(fbct.totint)

            for index_fbct in xrange(len(fbct.totint)):

                # rescale bct and bunch length to fbct time stamps
                index_bunch_length = np.abs(bunch_length.t_stamps-fbct.t_stamps[index_fbct]).argmin()
                fbct_bunch_length[index_fbct,:] = bunch_length.blen[index_bunch_length,:]

                if fbct.totint[index_fbct] != 0:
                    index_bct = np.abs(bct.t_stamps - fbct.t_stamps[index_fbct]).argmin()
                    total_intensity[index_fbct] = bct.values[index_bct]
                    fbct_correction_factors[index_fbct,:] = bct.values[index_bct]/fbct.totint[index_fbct]
                    
            # correct empty bunches to avoid divisions by 0
            mask_bunch_length = fbct_bunch_length == 0
            fbct_bunch_length[mask_bunch_length] = 1e-20
            bunch_int_sq[mask_bunch_length] = 0

            bunch_length_seconds[beam_ctr] = fbct_bunch_length
            bunch_int_sq_corrected = bunch_int_sq * fbct_correction_factors
            bunch_charge_squared[beam_ctr] = bunch_int_sq_corrected
            total_intensity_beams[beam_ctr] = total_intensity

        # Identical for both beams
        beam_energy_GeV = np.zeros_like(fbct.totint)
        energy_ob = energy(fill_dict, beam=1)
        for index_fbct in xrange(len(fbct.totint)):
            index_energy = np.abs(energy_ob.t_stamps - fbct.t_stamps[index_fbct]).argmin()
            beam_energy_GeV[index_fbct] = energy_ob.energy[index_energy]
        dipole_bField_T = beam_energy_GeV / 7000. * 8.33
        self.t_stamps = fbct.t_stamps
        self.t_stamps_hrs = (self.t_stamps - self.t_stamps[0]) / 3600.


        ## Impedance
        rho_0 = copper_rho_Ohm_m(temp_beam_screen_K)*np.ones_like(self.t_stamps)

        #extrapolate to 4.0
        f_rrr = copper_rho_Ohm_m(273)/(copper_rho_Ohm_m(4.2) - (copper_rho_Ohm_m(4.3) - copper_rho_Ohm_m(4.2))*2)
        delta_rho = rho_0*10.**-2.69*(f_rrr*dipole_bField_T)**1.055
        rho_tot_half_cell = rho_0 * l_halfcell + delta_rho * l_dip_halfcell
        bunch_contribution = np.sum((bunch_charge_squared[1]*bunch_length_seconds[1]**(-1.5)) + (bunch_charge_squared[2]*bunch_length_seconds[2]**(1.5)),axis=1) 
        weld_factor = np.sqrt(rho_stainless_steel_20K/copper_rho_Ohm_m(temp_beam_screen_K)) * weld_surface_factor 

        # neglect magnetic field in quadrupoles
        self.impedance_half_cell_heat_load = impedance_constant * (np.sqrt(rho_0)*(l_halfcell-l_dip_halfcell) + np.sqrt(rho_0+delta_rho)*l_dip_halfcell) * bunch_contribution * (1+weld_factor)

        
        ## Synchrotron radiation
        beam_gamma = beam_energy_GeV / (m_p*const_c**2/const_e/1e9)
        energy_loss_turn_single_proton = (const_e**2 * beam_gamma**4) / (3.*epsilon_0*lhc_bending_radius)
        beam_power_loss = energy_loss_turn_single_proton * (total_intensity_beams[1] + total_intensity_beams[2]) * lhc_f_rev

        self.synchRad_power_half_cell = beam_power_loss / (2.*np.pi*lhc_bending_radius) * l_dip_halfcell

    def show_a_plot(self):
        import matplotlib.pyplot as plt

        fig = plt.figure()
        title = 'Impedance and SynchRad heat load for fill %i' % self.filln
        fig.canvas.set_window_title(title)
        plt.suptitle(title, fontsize=25)
        sp = plt.subplot(1,1,1)
        sp.set_xlabel('Time [h]')
        sp.set_ylabel('Heat load [W/m]')
            
        sp.plot(self.t_stamps_hrs, self.synchRad_power_half_cell/l_halfcell, label='Synchrotron Radiation')
        sp.plot(self.t_stamps_hrs, self.impedance_half_cell_heat_load/l_halfcell, label='Impedance')

        plt.legend()
        plt.show()

# Rough estimate for Impedance
n_bunch = 2808.
sq_b_charge = (1.1e11)**2
rho = 1.7e-10
blen = 1.2e-9
n_beams = 2.
const = 1./lhc_circumference * gamma(0.75) / aperture_radius * (const_e/2./np.pi)**2 * np.sqrt(const_c*const_Z0/2.) 
weld_factor = np.sqrt(rho_stainless_steel_20K/rho) * weld_surface_factor
rough_estimate = const * (blen)**(-1.5) *n_bunch * sq_b_charge * np.sqrt(rho) * n_beams * (1 + weld_factor)

print('Rough estimate of the impedance power loss per m for %e bunch charge and 2808 bunches, with weld: %.3f mW/m' % (np.sqrt(sq_b_charge), rough_estimate*1e3))

# To test the module
hli = HeatLoadImpedanceSynchRad(5219)
hli.show_a_plot()
