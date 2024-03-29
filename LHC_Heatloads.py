import numpy as np
from . import TimberManager as tm

class Heatload:
    def __init__(self, timber_variable, sector):

        if type(timber_variable) is str:
            dict_timber = tm.parse_timber_file(timber_variable, verbose=True)
            timber_variable_hl = dict_timber[get_variable_dict(sector)['ARC_AVG']]

        elif type(timber_variable) is dict:
            timber_variable_hl = timber_variable[get_variable_dict(sector)['ARC_AVG']]

        # print np.squeeze(np.array(timber_variable_hl.values)).shape
        # print np.squeeze(np.array(timber_variable_hl.values))

        self.t_stamps = np.float_(np.array(timber_variable_hl.t_stamps))
        self.hl = np.squeeze(np.float_(np.array(timber_variable_hl.values)))



def get_variable_dict(sector):
    var_dict = {}
    var_dict['ARC_AVG'] = 'S%d_QBS_AVG_ARC.POSST'%sector

    return var_dict

def variable_list(sectors=[12,23,34,45,56,67,78,81]):
    var_list = []
    for sector in sectors:
        var_list += list(get_variable_dict(sector).values())

    return var_list

def sector_list():
    sector_list = [12,23,34,45,56,67,78,81]

    return sector_list


def arc_average_correction_factors():
    corr_factors = [1.3, 1.24, 1.22, 1.28, 1.26, 1.22, 1.24, 1.3]

    return corr_factors





average_arcs_variable_list = variable_list
variable_lists_heatloads = {}

variable_lists_heatloads['AVG_ARC'] = average_arcs_variable_list()

variable_lists_heatloads['Q4D2s_IR1'] = 'QRLFD_04L1_QBS947.POSST QRLFC_04R1_QBS947.POSST'.split()
variable_lists_heatloads['Q4D2s_IR5'] = 'QRLFC_04L5_QBS947.POSST QRLFD_04R5_QBS947.POSST'.split()
variable_lists_heatloads['Q4D2s_IR2'] = 'QRLFE_04L2_QBS947.POSST QRLFF_04R2_QBS947.POSST'.split()
variable_lists_heatloads['Q4D2s_IR8'] = 'QRLFE_04L8_QBS947.POSST QRLFF_04R8_QBS947.POSST'.split()

variable_lists_heatloads['Q6s_IR1'] = 'QRLEC_06L1_QBS947.POSST QRLEC_06R1_QBS947.POSST'.split()
variable_lists_heatloads['Q6s_IR5'] = 'QRLEC_06L5_QBS947.POSST QRLEC_06R5_QBS947.POSST'.split()
variable_lists_heatloads['Q6s_IR2'] = 'QRLEA_06L2_QBS947.POSST QRLEA_06R2_QBS947.POSST'.split()
variable_lists_heatloads['Q6s_IR8'] = 'QRLEA_06L8_QBS947.POSST QRLDE_06R8_QBS947.POSST'.split()

variable_lists_heatloads['Q5s_IR1'] = 'QRLEC_05L1_QBS947.POSST QRLEC_05R1_QBS947.POSST'.split()
variable_lists_heatloads['Q5s_IR5'] = 'QRLEC_05L5_QBS947.POSST QRLEC_05R5_QBS947.POSST'.split()
variable_lists_heatloads['Q5s_IR2'] = 'QRLEA_05L2_QBS947.POSST QRLEA_05R2_QBS947.POSST'.split()
variable_lists_heatloads['Q5s_IR8'] = 'QRLEA_05L8_QBS947.POSST QRLEA_05R8_QBS947.POSST'.split()

variable_lists_heatloads['IT_IR1'] = 'QRLGA_03L1_QBS947.POSST QRLGC_03R1_QBS947.POSST'.split()
variable_lists_heatloads['IT_IR5'] = 'QRLGD_03L5_QBS947.POSST QRLGB_03R5_QBS947.POSST'.split()
variable_lists_heatloads['IT_IR2'] = 'QRLGF_03L2_QBS947.POSST QRLGE_03R2_QBS947.POSST'.split()
variable_lists_heatloads['IT_IR8'] = 'QRLGF_03L8_QBS947.POSST QRLGE_03R8_QBS947.POSST'.split()

variable_lists_heatloads['special_HC_Q1'] = 'QRLAA_13L5_QBS943_Q1.POSST QRLAA_13R4_QBS947_Q1.POSST QRLAA_33L5_QBS947_Q1.POSST QRLAB_31L2_QBS943_Q1.POSST QRLAA_17L6_QBS943_Q1.POSST QRLAB_15R2_QBS943_Q1.POSST QRLAB_27L8_QBS947_Q1.POSST QRLAD_33R2_QBS947_Q1.POSST'.split()
variable_lists_heatloads['special_HC_D2'] = 'QRLAA_13L5_QBS943_D2.POSST QRLAA_13R4_QBS947_D2.POSST QRLAA_33L5_QBS947_D2.POSST QRLAB_31L2_QBS943_D2.POSST QRLAA_17L6_QBS943_D2.POSST QRLAB_15R2_QBS943_D2.POSST QRLAB_27L8_QBS947_D2.POSST QRLAD_33R2_QBS947_D2.POSST'.split()
variable_lists_heatloads['special_HC_D3'] = 'QRLAA_13L5_QBS943_D3.POSST QRLAA_13R4_QBS947_D3.POSST QRLAA_33L5_QBS947_D3.POSST QRLAB_31L2_QBS943_D3.POSST QRLAA_17L6_QBS943_D3.POSST QRLAB_15R2_QBS943_D3.POSST QRLAB_27L8_QBS947_D3.POSST QRLAD_33R2_QBS947_D3.POSST'.split()
variable_lists_heatloads['special_HC_D4'] = 'QRLAA_13L5_QBS943_D4.POSST QRLAA_13R4_QBS947_D4.POSST QRLAA_33L5_QBS947_D4.POSST QRLAB_31L2_QBS943_D4.POSST QRLAA_17L6_QBS943_D4.POSST QRLAB_15R2_QBS943_D4.POSST QRLAB_27L8_QBS947_D4.POSST QRLAD_33R2_QBS947_D4.POSST'.split()
variable_lists_heatloads['special_total'] = 'QRLAA_13L5_QBS943.POSST QRLAA_13R4_QBS947.POSST QRLAA_33L5_QBS947.POSST QRLAB_31L2_QBS943.POSST QRLAA_17L6_QBS943.POSST QRLAB_15R2_QBS943.POSST QRLAB_27L8_QBS947.POSST QRLAD_33R2_QBS947.POSST'.split()

variable_lists_heatloads['special_DS'] = 'QRLAB_11R2_QBS947_BP.POSST QRLAB_11R2_QBS947_CC1.POSST QRLAB_11R2_QBS947_CC2.POSST'.split()
variable_lists_heatloads['special_DS_total'] = ['QRLAB_11R2_QBS947.POSST']


for kk in ['special_' + m for m in ['HC_Q1', 'HC_D2', 'HC_D3', 'HC_D4', 'DS']]:
    for ll in ['B1', 'B2', '_COR', 'B1_COR', 'B2_COR']:
        var_list = []
        for mm in variable_lists_heatloads[kk]:
            name = mm.split('.POSST')[0]
            var_list.append(name + ll + '.POSST')
        variable_lists_heatloads[kk + ll] = var_list

for kk in ['special_total', 'special_DS_total']:
    var_list = []
    for mm in variable_lists_heatloads[kk]:
        name = mm.split('.POSST')[0]
        var_list.append(name + '_COR.POSST')
    variable_lists_heatloads[kk + '_COR'] = var_list


# Deprecated! Old variables are wrong!
variable_lists_heatloads['MODEL'] = ['LHC.QBS_CALCULATED_ARC_IMPED.B1', 'LHC.QBS_CALCULATED_ARC_IMPED.B2',
                                     'LHC.QBS_CALCULATED_ARC_SYNCH_RAD.B1', 'LHC.QBS_CALCULATED_ARC_SYNCH_RAD.B2',
                                     'LHC.QBS_CALCULATED_ARC.TOTAL']

heat_loads_plot_sets = {}
for kk in variable_lists_heatloads:
     heat_loads_plot_sets[kk] = variable_lists_heatloads[kk]


heat_loads_plot_sets['dipoles_31L2'] = 'QRLAB_31L2_QBS943_D2.POSST QRLAB_31L2_QBS943_D3.POSST QRLAB_31L2_QBS943_D4.POSST'.split()
heat_loads_plot_sets['dipoles_13L5'] = 'QRLAA_13L5_QBS943_D2.POSST QRLAA_13L5_QBS943_D3.POSST QRLAA_13L5_QBS943_D4.POSST'.split()
heat_loads_plot_sets['dipoles_33L5'] = 'QRLAA_33L5_QBS947_D2.POSST QRLAA_33L5_QBS947_D3.POSST QRLAA_33L5_QBS947_D4.POSST'.split()
heat_loads_plot_sets['dipoles_13R4'] = 'QRLAA_13R4_QBS947_D2.POSST QRLAA_13R4_QBS947_D3.POSST QRLAA_13R4_QBS947_D4.POSST'.split()
heat_loads_plot_sets['dipoles_17L6'] = 'QRLAA_17L6_QBS943_D2.POSST QRLAA_17L6_QBS943_D3.POSST QRLAA_17L6_QBS943_D4.POSST'.split()
heat_loads_plot_sets['dipoles_15R2'] = 'QRLAB_15R2_QBS943_D2.POSST QRLAB_15R2_QBS943_D3.POSST QRLAB_15R2_QBS943_D4.POSST'.split()
heat_loads_plot_sets['dipoles_27L8'] = 'QRLAB_27L8_QBS947_D2.POSST QRLAB_27L8_QBS947_D3.POSST QRLAB_27L8_QBS947_D4.POSST'.split() 
heat_loads_plot_sets['dipoles_33R2'] = 'QRLAD_33R2_QBS947_D2.POSST QRLAD_33R2_QBS947_D3.POSST QRLAD_33R2_QBS947_D4.POSST'.split()

heat_loads_plot_sets['quadrupole_31L2'] = 'QRLAB_31L2_QBS943_Q1.POSST'.split()
heat_loads_plot_sets['quadrupole_13L5'] = 'QRLAA_13L5_QBS943_Q1.POSST'.split()
heat_loads_plot_sets['quadrupole_33L5'] = 'QRLAA_33L5_QBS947_Q1.POSST'.split()
heat_loads_plot_sets['quadrupole_13R4'] = 'QRLAA_13R4_QBS947_Q1.POSST'.split()
heat_loads_plot_sets['quadrupole_17L6'] = 'QRLAA_17L6_QBS943_Q1.POSST'.split()
heat_loads_plot_sets['quadrupole_15R2'] = 'QRLAB_15R2_QBS943_Q1.POSST'.split()
heat_loads_plot_sets['quadrupole_27L8'] = 'QRLAB_27L8_QBS947_Q1.POSST'.split() 
heat_loads_plot_sets['quadrupole_33R2'] = 'QRLAD_33R2_QBS947_Q1.POSST'.split()

heat_loads_plot_sets['InnerTriplets_IR15'] = variable_lists_heatloads['IT_IR1']+variable_lists_heatloads['IT_IR5']
heat_loads_plot_sets['InnerTriplets_IR28'] = variable_lists_heatloads['IT_IR2']+variable_lists_heatloads['IT_IR8']
heat_loads_plot_sets['Arcs'] = variable_lists_heatloads['AVG_ARC']
heat_loads_plot_sets['Q5s_IR15'] = variable_lists_heatloads['Q5s_IR1']+variable_lists_heatloads['Q5s_IR5']
heat_loads_plot_sets['Q5s_IR28'] = variable_lists_heatloads['Q5s_IR2']+variable_lists_heatloads['Q5s_IR8']
heat_loads_plot_sets['Q6s_IR15'] = variable_lists_heatloads['Q6s_IR1']+variable_lists_heatloads['Q6s_IR5']
heat_loads_plot_sets['Q6s_IR28'] = variable_lists_heatloads['Q6s_IR2']+variable_lists_heatloads['Q6s_IR8']
heat_loads_plot_sets['special_HC_Q1']	= variable_lists_heatloads['special_HC_Q1']
heat_loads_plot_sets['special_HC_dipoles'] = heat_loads_plot_sets['dipoles_31L2'] + heat_loads_plot_sets['dipoles_13L5']\
                                             + heat_loads_plot_sets['dipoles_33L5'] +heat_loads_plot_sets['dipoles_13R4']\
                                             + heat_loads_plot_sets['dipoles_17L6'] + heat_loads_plot_sets['dipoles_15R2']\
                                             + heat_loads_plot_sets['dipoles_27L8'] + heat_loads_plot_sets['dipoles_33R2']

heat_loads_plot_sets['Q4D2s_IR15'] = variable_lists_heatloads['Q4D2s_IR1']+ variable_lists_heatloads['Q4D2s_IR5']
heat_loads_plot_sets['Q4D2s_IR28'] = variable_lists_heatloads['Q4D2s_IR2']+ variable_lists_heatloads['Q4D2s_IR8']

heat_loads_plot_sets['Q6s_IR37'] = 'QRLEA_06L3_QBS947.POSST QRLEA_06L7_QBS947.POSST QRLEA_06R3_QBS947.POSST QRLEA_06R7_QBS947.POSST'.split()

heat_loads_plot_sets['D3s_IR4'] = ['QRLEB_05L4_QBS947.POSST', 'QRLEB_05R4_QBS947.POSST']

heat_loads_plot_sets['dipoles_31L2_B1'] = 'QRLAB_31L2_QBS943_D2B1.POSST QRLAB_31L2_QBS943_D3B1.POSST QRLAB_31L2_QBS943_D4B1.POSST'.split()
heat_loads_plot_sets['dipoles_13L5_B1'] = 'QRLAA_13L5_QBS943_D2B1.POSST QRLAA_13L5_QBS943_D3B1.POSST QRLAA_13L5_QBS943_D4B1.POSST'.split()
heat_loads_plot_sets['dipoles_33L5_B1'] = 'QRLAA_33L5_QBS947_D2B1.POSST QRLAA_33L5_QBS947_D3B1.POSST QRLAA_33L5_QBS947_D4B1.POSST'.split()[:-1]
heat_loads_plot_sets['dipoles_13R4_B1'] = 'QRLAA_13R4_QBS947_D2B1.POSST QRLAA_13R4_QBS947_D3B1.POSST QRLAA_13R4_QBS947_D4B1.POSST'.split()[:-1]
heat_loads_plot_sets['dipoles_17L6_B1'] = 'QRLAA_17L6_QBS943_D2B1.POSST QRLAA_17L6_QBS943_D3B1.POSST QRLAA_17L6_QBS943_D4B1.POSST'.split()
heat_loads_plot_sets['dipoles_15R2_B1'] = 'QRLAB_15R2_QBS943_D2B1.POSST QRLAB_15R2_QBS943_D3B1.POSST QRLAB_15R2_QBS943_D4B1.POSST'.split()
heat_loads_plot_sets['dipoles_27L8_B1'] = 'QRLAB_27L8_QBS947_D2B1.POSST QRLAB_27L8_QBS947_D3B1.POSST QRLAB_27L8_QBS947_D4B1.POSST'.split() 
heat_loads_plot_sets['dipoles_33R2_B1'] = 'QRLAD_33R2_QBS947_D2B1.POSST QRLAD_33R2_QBS947_D3B1.POSST QRLAD_33R2_QBS947_D4B1.POSST'.split()
heat_loads_plot_sets['special_HC_dipoles_B1'] = heat_loads_plot_sets['dipoles_31L2_B1'] + heat_loads_plot_sets['dipoles_13L5_B1']\
                                                + heat_loads_plot_sets['dipoles_33L5_B1'] +heat_loads_plot_sets['dipoles_13R4_B1']\
                                                + heat_loads_plot_sets['dipoles_17L6_B1'] + heat_loads_plot_sets['dipoles_15R2_B1']\
                                                + heat_loads_plot_sets['dipoles_27L8_B1'] + heat_loads_plot_sets['dipoles_33R2_B1']


heat_loads_plot_sets['dipoles_31L2_B2'] = 'QRLAB_31L2_QBS943_D2B2.POSST QRLAB_31L2_QBS943_D3B2.POSST QRLAB_31L2_QBS943_D4B2.POSST'.split()
heat_loads_plot_sets['dipoles_13L5_B2'] = 'QRLAA_13L5_QBS943_D2B2.POSST QRLAA_13L5_QBS943_D3B2.POSST QRLAA_13L5_QBS943_D4B2.POSST'.split()
heat_loads_plot_sets['dipoles_33L5_B2'] = 'QRLAA_33L5_QBS947_D2B2.POSST QRLAA_33L5_QBS947_D3B2.POSST QRLAA_33L5_QBS947_D4B2.POSST'.split()[:-1]
heat_loads_plot_sets['dipoles_13R4_B2'] = 'QRLAA_13R4_QBS947_D2B2.POSST QRLAA_13R4_QBS947_D3B2.POSST QRLAA_13R4_QBS947_D4B2.POSST'.split()[:-1]
heat_loads_plot_sets['dipoles_17L6_B2'] = 'QRLAA_17L6_QBS943_D2B2.POSST QRLAA_17L6_QBS943_D3B2.POSST QRLAA_17L6_QBS943_D4B2.POSST'.split()
heat_loads_plot_sets['dipoles_15R2_B2'] = 'QRLAB_15R2_QBS943_D2B2.POSST QRLAB_15R2_QBS943_D3B2.POSST QRLAB_15R2_QBS943_D4B2.POSST'.split()
heat_loads_plot_sets['dipoles_27L8_B2'] = 'QRLAB_27L8_QBS947_D2B2.POSST QRLAB_27L8_QBS947_D3B2.POSST QRLAB_27L8_QBS947_D4B2.POSST'.split() 
heat_loads_plot_sets['dipoles_33R2_B2'] = 'QRLAD_33R2_QBS947_D2B2.POSST QRLAD_33R2_QBS947_D3B2.POSST QRLAD_33R2_QBS947_D4B2.POSST'.split()
heat_loads_plot_sets['special_HC_dipoles_B2'] = heat_loads_plot_sets['dipoles_31L2_B2'] + heat_loads_plot_sets['dipoles_13L5_B2']\
                                                + heat_loads_plot_sets['dipoles_33L5_B2'] +heat_loads_plot_sets['dipoles_13R4_B2']\
                                                + heat_loads_plot_sets['dipoles_17L6_B2'] + heat_loads_plot_sets['dipoles_15R2_B2']\
                                                + heat_loads_plot_sets['dipoles_27L8_B2'] + heat_loads_plot_sets['dipoles_33R2_B2']

instrum_prefixes = 'QRLAB_31L2_QBS943_ QRLAA_13L5_QBS943_ QRLAA_33L5_QBS947_ QRLAA_13R4_QBS947_ QRLAA_17L6_QBS943_ QRLAB_15R2_QBS943_ QRLAB_27L8_QBS947_ QRLAD_33R2_QBS947_'.split()
for pp in instrum_prefixes:
    cell_name = pp.split('_')[1]
    for nn in 'Q1 D2 D3 D4'.split():
        vlist = []
        for bb in ['B1', 'B2', '']:
            vlist.append(pp+nn+bb+'.POSST')
        heat_loads_plot_sets[cell_name+'_'+nn+bb] = vlist


def groups_dict():
    #~ dict_hl_groups = {}
    #~ dict_hl_groups['InnerTriplets'] = variable_lists_heatloads['IT_IR1']+variable_lists_heatloads['IT_IR5']+\
        #~ variable_lists_heatloads['IT_IR2']+variable_lists_heatloads['IT_IR8']
    #~ dict_hl_groups['Arcs'] = variable_lists_heatloads['AVG_ARC']
    #~ dict_hl_groups['Q5s'] = variable_lists_heatloads['Q5s_IR1']+variable_lists_heatloads['Q5s_IR5']+\
        #~ variable_lists_heatloads['Q5s_IR2']+variable_lists_heatloads['Q5s_IR8']
    #~ dict_hl_groups['Q6s'] = variable_lists_heatloads['Q6s_IR1']+variable_lists_heatloads['Q6s_IR5']+\
        #~ variable_lists_heatloads['Q6s_IR2']+variable_lists_heatloads['Q6s_IR8']
    #~ dict_hl_groups['Q4D2s'] =  variable_lists_heatloads['Q4D2s_IR1']+ variable_lists_heatloads['Q4D2s_IR5']+\
        #~ variable_lists_heatloads['Q4D2s_IR2']+ variable_lists_heatloads['Q4D2s_IR8']

    #~ dict_hl_groups['special_HC_Q1']	= variable_lists_heatloads['special_HC_Q1']
    #~ dict_hl_groups['special_HC_dipoles'] = variable_lists_heatloads['special_HC_D2']+\
        #~ variable_lists_heatloads['special_HC_D3']+variable_lists_heatloads['special_HC_D4']
    #~ dict_hl_groups['special_HC_total'] = variable_lists_heatloads['special_total']
    raise ValueError('Feature Discontinued!')

    #~ return heat_loads_plot_sets


cryogenic_length = {}

cryogenic_length['AVG_ARC'] = [53.45]
cryogenic_length['Arcs'] = [53.45]
cryogenic_length['MODEL'] = [53.45]

cryogenic_length['Q4D2s_IR1'] = [19.4]
cryogenic_length['Q4D2s_IR5'] = [19.4]
cryogenic_length['Q4D2s_IR2'] = [22.8]
cryogenic_length['Q4D2s_IR8'] = [22.8]

cryogenic_length['Q6s_IR1'] = [8.2]
cryogenic_length['Q6s_IR5'] = [8.2]
cryogenic_length['Q6s_IR2'] = [12.]
cryogenic_length['Q6s_IR8'] = [12.]

cryogenic_length['Q5s_IR1'] = [8.2]
cryogenic_length['Q5s_IR5'] = [8.2]
cryogenic_length['Q5s_IR2'] = [13.]
cryogenic_length['Q5s_IR8'] = [13.]

cryogenic_length['IT_IR1'] = [40.]
cryogenic_length['IT_IR5'] = [40.]
cryogenic_length['IT_IR2'] = [50.]
cryogenic_length['IT_IR8'] = [50.]

cryogenic_length['special_HC_Q1'] = [3.1]
cryogenic_length['special_HC_D2'] = [14.3]
cryogenic_length['special_HC_D3'] = [14.3]
cryogenic_length['special_HC_D4'] = [14.3]
cryogenic_length['special_total'] = [53.45]

cryogenic_length['dipoles_31L2'] = [14.3]

magnet_length = {}

magnet_length['AVG_ARC'] = [53.45]
magnet_length['Arcs'] = [53.45]
magnet_length['MODEL'] = [53.45]
magnet_length['Q4D2s_IR1'] = [18.08]
magnet_length['Q4D2s_IR5'] = [18.08]
magnet_length['Q4D2s_IR2'] = [21.39]
magnet_length['Q4D2s_IR8'] = [21.39]

magnet_length['Q6s_IR1'] = [4.8]
magnet_length['Q6s_IR5'] = [4.8]
magnet_length['Q6s_IR2'] = [8.567]
magnet_length['Q6s_IR8'] = [8.567]

magnet_length['Q5s_IR1'] = [4.8]
magnet_length['Q5s_IR5'] = [4.8]
magnet_length['Q5s_IR2'] = [7.181]
magnet_length['Q5s_IR8'] = [7.181]

magnet_length['IT_IR1'] = [36.98]
magnet_length['IT_IR5'] = [36.96]
magnet_length['IT_IR2'] = [44.91]
magnet_length['IT_IR8'] = [44.91]

magnet_length['special_HC_Q1'] = [3.1]
magnet_length['special_HC_D2'] = [14.3]
magnet_length['special_HC_D3'] = [14.3]
magnet_length['special_HC_D4'] = [14.3]
magnet_length['special_HC_dipoles'] = [14.3]
magnet_length['special_total'] = [53.45]

for kk in ['31L2', '13L5', '33L5', '13R4']:
    magnet_length['dipoles_'+kk] = [14.3]


def groups_length_dict(length='cryogenic_length'):

    name_dict = variable_lists_heatloads

    if length == 'magnet_length':
        len_dict = magnet_length
    elif length == 'cryogenic_length':
        len_dict = cryogenic_length

    dict_len_groups = {}

    dict_len_groups['InnerTriplets'] = []
    dict_len_groups['Arcs'] = []
    dict_len_groups['Q5s'] = []
    dict_len_groups['Q6s'] = []
    dict_len_groups['Q4D2s'] = []
    dict_len_groups['special_HC_Q1'] = []
    dict_len_groups['special_HC_dipoles'] = []
    dict_len_groups['special_HC_total'] = []


    dict_len_groups['InnerTriplets'].extend(len_dict['IT_IR1']*len(name_dict['IT_IR1']))
    dict_len_groups['InnerTriplets'].extend(len_dict['IT_IR5']*len(name_dict['IT_IR5']))
    dict_len_groups['InnerTriplets'].extend(len_dict['IT_IR2']*len(name_dict['IT_IR2']))
    dict_len_groups['InnerTriplets'].extend(len_dict['IT_IR8']*len(name_dict['IT_IR8']))

    dict_len_groups['Arcs'].extend(len_dict['AVG_ARC']*len(name_dict['AVG_ARC']))

    dict_len_groups['Q5s'].extend(len_dict['Q5s_IR1']*len(name_dict['Q5s_IR1']))
    dict_len_groups['Q5s'].extend(len_dict['Q5s_IR5']*len(name_dict['Q5s_IR5']))
    dict_len_groups['Q5s'].extend(len_dict['Q5s_IR2']*len(name_dict['Q5s_IR2']))
    dict_len_groups['Q5s'].extend(len_dict['Q5s_IR8']*len(name_dict['Q5s_IR8']))

    dict_len_groups['Q6s'].extend(len_dict['Q6s_IR1']*len(name_dict['Q6s_IR1']))
    dict_len_groups['Q6s'].extend(len_dict['Q6s_IR5']*len(name_dict['Q6s_IR5']))
    dict_len_groups['Q6s'].extend(len_dict['Q6s_IR2']*len(name_dict['Q6s_IR2']))
    dict_len_groups['Q6s'].extend(len_dict['Q6s_IR8']*len(name_dict['Q6s_IR8']))

    dict_len_groups['Q4D2s'].extend(len_dict['Q4D2s_IR1']*len(name_dict['Q4D2s_IR1']))
    dict_len_groups['Q4D2s'].extend(len_dict['Q4D2s_IR5']*len(name_dict['Q4D2s_IR5']))
    dict_len_groups['Q4D2s'].extend(len_dict['Q4D2s_IR2']*len(name_dict['Q4D2s_IR2']))
    dict_len_groups['Q4D2s'].extend(len_dict['Q4D2s_IR8']*len(name_dict['Q4D2s_IR8']))

    dict_len_groups['special_HC_Q1'].extend(len_dict['special_HC_Q1']*len(name_dict['special_HC_Q1']))
    dict_len_groups['special_HC_dipoles'].extend(len_dict['special_HC_D2']*len(name_dict['special_HC_D2']))
    dict_len_groups['special_HC_dipoles'].extend(len_dict['special_HC_D3']*len(name_dict['special_HC_D3']))
    dict_len_groups['special_HC_dipoles'].extend(len_dict['special_HC_D4']*len(name_dict['special_HC_D4']))

#    dict_len_groups['special_HC_total'].extend(len_dict['special_total']*len(name_dict['special_total']))

    return dict_len_groups


def sector_all_variables(sectors):
    sectors = np.array(sectors, ndmin=1)

    sector_variable_list = []
    for sector in sectors:
        sector_R = str(sector)[0]
        sector_L = str(sector)[1]

        variable_list_R = ['QRLBA_09R','QRLAB_11R','QRLAA_13R','QRLAB_15R','QRLAA_17R','QRLAB_19R','QRLAA_21R',
                           'QRLAB_23R','QRLAA_25R','QRLAB_27R','QRLAA_29R','QRLAC_31R','QRLAD_33R']
        variable_list_L = ['QRLAA_33L','QRLAB_31L','QRLAA_29L','QRLAB_27L','QRLAA_25L','QRLAB_23L','QRLAA_21L',
                           'QRLAB_19L','QRLAA_17L','QRLAB_15L','QRLAA_13L','QRLAB_11L','QRLBA_09L',
                           'QRLBB_09L','QRLAH_11L','QRLAG_13L','QRLAH_15L','QRLAG_17L','QRLAF_25L','QRLAE_25L']#this line is for special cases

        for variable in variable_list_R:
            for nqbs in [943,947]:
                curr_variable = variable+'%s_QBS%d.POSST'%(sector_R, nqbs)
                sector_variable_list.append(curr_variable)
        for variable in variable_list_L:
            for nqbs in [947,943]:
                curr_variable = variable+'%s_QBS%d.POSST'%(sector_L, nqbs)
                sector_variable_list.append(curr_variable)

    return sector_variable_list

def get_dict_magnet_lengths():
    dict_lengths = {}
    for kk in list(variable_lists_heatloads.keys()):
        for device in variable_lists_heatloads[kk]:
            try:
                dict_lengths[device] = magnet_length[kk][0]
            except KeyError:
                print(f"WARNING: {device}, {kk} not found in magnet_length!")
    return dict_lengths

def get_dict_cryostat_lengths():
    dict_lengths = {}
    for kk in list(variable_lists_heatloads.keys()):
        for device in variable_lists_heatloads[kk]:
            dict_lengths[device] = cryogenic_length[kk][0]
    return dict_lengths


arcs_varnames_static = [\
 'QRLBA_09L5_QBS943.POSST',
 'QRLAA_13L8_QBS947.POSST',
 'QRLAB_11L2_QBS943.POSST',
 'QRLAC_31R4_QBS943.POSST',
 'QRLAB_11L1_QBS947.POSST',
 'QRLAB_23L4_QBS943.POSST',
 'QRLAA_21L1_QBS943.POSST',
 'QRLAB_27L8_QBS943.POSST',
 'QRLAB_11L6_QBS943.POSST',
 'QRLAB_23L2_QBS947.POSST',
 'QRLAA_17R4_QBS947.POSST',
 'QRLAB_15R8_QBS947.POSST',
 'QRLAA_25R3_QBS943.POSST',
 'QRLAA_13R7_QBS943.POSST',
 'QRLAB_23L5_QBS943.POSST',
 'QRLAA_17L2_QBS943.POSST',
 'QRLBA_09R8_QBS943.POSST',
 'QRLAA_21L8_QBS947.POSST',
 'QRLAA_25R7_QBS943.POSST',
 'QRLAB_15R2_QBS947.POSST',
 'QRLAH_11L3_QBS943.POSST',
 'QRLAA_25R1_QBS943.POSST',
 'QRLAC_31R5_QBS943.POSST',
 'QRLAC_31R7_QBS943.POSST',
 'QRLAG_17L3_QBS943.POSST',
 'QRLAB_27R6_QBS943.POSST',
 'QRLAA_33L1_QBS947.POSST',
 'QRLAA_17R8_QBS943.POSST',
 'QRLBA_09R1_QBS947.POSST',
 'QRLAA_13L6_QBS943.POSST',
 'QRLAA_13L5_QBS947.POSST',
 'QRLAB_23R7_QBS943.POSST',
 'QRLAA_25R6_QBS943.POSST',
 'QRLAA_25L3_QBS947.POSST',
 'QRLBA_09R8_QBS947.POSST',
 'QRLAD_33R8_QBS947.POSST',
 'QRLAC_31R4_QBS947.POSST',
 'QRLAB_27L7_QBS947.POSST',
 'QRLAA_29R3_QBS943.POSST',
 'QRLAD_33R4_QBS943.POSST',
 'QRLAB_23R6_QBS943.POSST',
 'QRLAA_13R3_QBS947.POSST',
 'QRLAB_27R2_QBS943.POSST',
 'QRLAB_31L2_QBS947.POSST',
 'QRLAB_19L2_QBS947.POSST',
 'QRLAC_31R7_QBS947.POSST',
 'QRLAC_31R1_QBS947.POSST',
 'QRLAB_19L7_QBS943.POSST',
 'QRLAA_17R2_QBS947.POSST',
 'QRLAB_15R4_QBS943.POSST',
 'QRLAB_27R3_QBS943.POSST',
 'QRLAB_27R8_QBS943.POSST',
 'QRLAB_27L5_QBS947.POSST',
 'QRLAA_29L1_QBS947.POSST',
 'QRLAB_19L4_QBS943.POSST',
 'QRLAB_31L6_QBS943.POSST',
 'QRLAB_19R6_QBS943.POSST',
 'QRLAD_33R2_QBS943.POSST',
 'QRLAA_13R6_QBS943.POSST',
 'QRLAA_13R8_QBS947.POSST',
 'QRLAB_19R6_QBS947.POSST',
 'QRLAB_27L4_QBS943.POSST',
 'QRLAB_23R4_QBS943.POSST',
 'QRLAA_17R1_QBS943.POSST',
 'QRLAA_33L2_QBS947.POSST',
 'QRLAA_13L2_QBS943.POSST',
 'QRLAB_15L2_QBS947.POSST',
 'QRLAA_17R7_QBS947.POSST',
 'QRLAB_19R5_QBS947.POSST',
 'QRLAB_11R6_QBS943.POSST',
 'QRLAA_29R8_QBS943.POSST',
 'QRLAB_19L8_QBS947.POSST',
 'QRLAA_29R2_QBS943.POSST',
 'QRLAA_29R4_QBS947.POSST',
 'QRLAA_13R4_QBS943.POSST',
 'QRLBA_09L7_QBS943.POSST',
 'QRLAA_21L1_QBS947.POSST',
 'QRLAG_13L3_QBS943.POSST',
 'QRLAA_13L4_QBS943.POSST',
 'QRLAA_21L6_QBS943.POSST',
 'QRLAA_13R5_QBS947.POSST',
 'QRLAB_15L5_QBS943.POSST',
 'QRLAB_11L4_QBS943.POSST',
 'QRLAA_17L7_QBS943.POSST',
 'QRLAB_11R7_QBS943.POSST',
 'QRLAA_13R1_QBS947.POSST',
 'QRLAA_21L5_QBS947.POSST',
 'QRLAA_21R4_QBS943.POSST',
 'QRLAB_19R1_QBS943.POSST',
 'QRLAB_23R1_QBS947.POSST',
 'QRLAB_19L3_QBS943.POSST',
 'QRLAA_21L2_QBS943.POSST',
 'QRLAD_33R5_QBS947.POSST',
 'QRLAB_31L3_QBS943.POSST',
 'QRLAE_25L4_QBS947.POSST',
 'QRLAB_31L4_QBS943.POSST',
 'QRLAA_33L6_QBS947.POSST',
 'QRLAA_17L4_QBS947.POSST',
 'QRLAB_19R7_QBS943.POSST',
 'QRLAF_25L8_QBS947.POSST',
 'QRLAB_31L8_QBS947.POSST',
 'QRLAB_15L5_QBS947.POSST',
 'QRLAB_23R8_QBS947.POSST',
 'QRLAA_17R2_QBS943.POSST',
 'QRLAB_15L4_QBS943.POSST',
 'QRLAA_21R6_QBS947.POSST',
 'QRLAA_17R6_QBS943.POSST',
 'QRLAA_21R8_QBS947.POSST',
 'QRLAA_17L1_QBS947.POSST',
 'QRLAB_19R3_QBS947.POSST',
 'QRLAA_25R8_QBS947.POSST',
 'QRLAB_15R3_QBS947.POSST',
 'QRLAB_27L2_QBS943.POSST',
 'QRLAD_33R8_QBS943.POSST',
 'QRLAA_33L7_QBS947.POSST',
 'QRLAA_21L5_QBS943.POSST',
 'QRLAB_15L4_QBS947.POSST',
 'QRLAB_11R1_QBS947.POSST',
 'QRLAA_33L4_QBS943.POSST',
 'QRLAB_15R4_QBS947.POSST',
 'QRLAB_19L1_QBS947.POSST',
 'QRLAB_31L5_QBS943.POSST',
 'QRLAB_15R1_QBS947.POSST',
 'QRLAD_33R3_QBS947.POSST',
 'QRLAA_25L5_QBS943.POSST',
 'QRLAB_11L8_QBS947.POSST',
 'QRLAA_33L6_QBS943.POSST',
 'QRLAA_29L2_QBS943.POSST',
 'QRLAA_13L6_QBS947.POSST',
 'QRLAA_17R6_QBS947.POSST',
 'QRLAA_13L7_QBS947.POSST',
 'QRLBA_09R4_QBS947.POSST',
 'QRLAA_29R8_QBS947.POSST',
 'QRLAA_29R5_QBS947.POSST',
 'QRLAA_29R5_QBS943.POSST',
 'QRLAA_17L6_QBS943.POSST',
 'QRLAB_27R5_QBS947.POSST',
 'QRLBA_09L8_QBS943.POSST',
 'QRLAA_13L7_QBS943.POSST',
 'QRLAE_25L4_QBS943.POSST',
 'QRLAA_25R4_QBS947.POSST',
 'QRLAA_33L1_QBS943.POSST',
 'QRLAB_11L2_QBS947.POSST',
 'QRLAB_23L8_QBS943.POSST',
 'QRLAD_33R6_QBS947.POSST',
 'QRLAB_11R6_QBS947.POSST',
 'QRLBA_09L7_QBS947.POSST',
 'QRLAB_27R7_QBS947.POSST',
 'QRLAA_29R3_QBS947.POSST',
 'QRLAA_29R1_QBS943.POSST',
 'QRLAB_11R4_QBS947.POSST',
 'QRLBA_09L2_QBS947.POSST',
 'QRLAB_15R2_QBS943.POSST',
 'QRLAA_25R2_QBS947.POSST',
 'QRLAA_17L2_QBS947.POSST',
 'QRLAB_19L5_QBS943.POSST',
 'QRLAA_13R8_QBS943.POSST',
 'QRLAB_31L4_QBS947.POSST',
 'QRLBB_09L3_QBS947.POSST',
 'QRLAB_23R1_QBS943.POSST',
 'QRLAB_15R5_QBS947.POSST',
 'QRLAA_29L4_QBS947.POSST',
 'QRLAB_11L7_QBS947.POSST',
 'QRLAA_17R4_QBS943.POSST',
 'QRLAA_21R7_QBS947.POSST',
 'QRLAB_27R1_QBS947.POSST',
 'QRLAB_27L6_QBS943.POSST',
 'QRLAB_31L1_QBS947.POSST',
 'QRLAA_25L1_QBS947.POSST',
 'QRLBA_09L1_QBS947.POSST',
 'QRLAA_21R6_QBS943.POSST',
 'QRLAB_27L8_QBS947.POSST',
 'QRLAA_29L3_QBS947.POSST',
 'QRLBA_09L1_QBS943.POSST',
 'QRLAA_25L2_QBS947.POSST',
 'QRLAB_15R7_QBS947.POSST',
 'QRLAD_33R3_QBS943.POSST',
 'QRLAA_21L4_QBS943.POSST',
 'QRLAB_23L1_QBS943.POSST',
 'QRLAD_33R7_QBS943.POSST',
 'QRLAB_11L7_QBS943.POSST',
 'QRLAA_13R3_QBS943.POSST',
 'QRLBA_09R7_QBS947.POSST',
 'QRLAB_27R7_QBS943.POSST',
 'QRLBA_09R5_QBS943.POSST',
 'QRLAB_19R4_QBS943.POSST',
 'QRLBA_09R4_QBS943.POSST',
 'QRLAB_15R5_QBS943.POSST',
 'QRLAB_15L8_QBS943.POSST',
 'QRLAB_15R7_QBS943.POSST',
 'QRLAB_19L7_QBS947.POSST',
 'QRLAB_15L7_QBS943.POSST',
 'QRLAB_23R8_QBS943.POSST',
 'QRLAA_13R4_QBS947.POSST',
 'QRLAA_25R5_QBS943.POSST',
 'QRLAA_17L1_QBS943.POSST',
 'QRLAB_15L7_QBS947.POSST',
 'QRLAB_15R8_QBS943.POSST',
 'QRLAB_15L6_QBS947.POSST',
 'QRLBA_09L4_QBS947.POSST',
 'QRLAA_13L5_QBS943.POSST',
 'QRLAD_33R7_QBS947.POSST',
 'QRLAB_15L1_QBS947.POSST',
 'QRLAB_27R3_QBS947.POSST',
 'QRLAB_19R8_QBS947.POSST',
 'QRLBA_09R2_QBS947.POSST',
 'QRLAA_17L5_QBS943.POSST',
 'QRLAA_21L6_QBS947.POSST',
 'QRLAA_25R3_QBS947.POSST',
 'QRLAB_23L4_QBS947.POSST',
 'QRLAC_31R5_QBS947.POSST',
 'QRLAB_11L4_QBS947.POSST',
 'QRLAA_25L1_QBS943.POSST',
 'QRLAA_25R2_QBS943.POSST',
 'QRLAB_11R3_QBS947.POSST',
 'QRLBA_09R2_QBS943.POSST',
 'QRLBA_09R3_QBS943.POSST',
 'QRLAB_11R7_QBS947.POSST',
 'QRLAB_23R3_QBS943.POSST',
 'QRLAA_17R8_QBS947.POSST',
 'QRLAH_15L3_QBS943.POSST',
 'QRLAA_29L6_QBS947.POSST',
 'QRLAB_27L6_QBS947.POSST',
 'QRLAB_11L6_QBS947.POSST',
 'QRLAA_21L8_QBS943.POSST',
 'QRLAH_11L3_QBS947.POSST',
 'QRLAB_27L1_QBS943.POSST',
 'QRLAA_21L3_QBS947.POSST',
 'QRLAD_33R5_QBS943.POSST',
 'QRLAA_29R1_QBS947.POSST',
 'QRLAB_15L8_QBS947.POSST',
 'QRLAA_33L3_QBS947.POSST',
 'QRLAA_25L2_QBS943.POSST',
 'QRLAB_27L2_QBS947.POSST',
 'QRLAA_25L3_QBS943.POSST',
 'QRLAA_33L7_QBS943.POSST',
 'QRLAB_27L5_QBS943.POSST',
 'QRLAA_25L5_QBS947.POSST',
 'QRLAB_23L1_QBS947.POSST',
 'QRLAD_33R1_QBS947.POSST',
 'QRLAA_29R6_QBS943.POSST',
 'QRLAA_21L4_QBS947.POSST',
 'QRLAB_15L6_QBS943.POSST',
 'QRLAB_19L3_QBS947.POSST',
 'QRLAB_19L5_QBS947.POSST',
 'QRLAC_31R3_QBS943.POSST',
 'QRLAA_17R7_QBS943.POSST',
 'QRLAB_19L4_QBS947.POSST',
 'QRLAA_13L4_QBS947.POSST',
 'QRLAA_25R4_QBS943.POSST',
 'QRLAA_29R6_QBS947.POSST',
 'QRLAB_15R6_QBS947.POSST',
 'QRLAA_33L5_QBS947.POSST',
 'QRLAB_15L1_QBS943.POSST',
 'QRLAB_19L6_QBS947.POSST',
 'QRLAD_33R1_QBS943.POSST',
 'QRLAA_21R1_QBS943.POSST',
 'QRLAB_11R2_QBS943.POSST',
 'QRLAG_17L3_QBS947.POSST',
 'QRLAA_21R3_QBS947.POSST',
 'QRLAA_13R2_QBS943.POSST',
 'QRLAB_11R4_QBS943.POSST',
 'QRLAB_31L2_QBS943.POSST',
 'QRLAA_13R7_QBS947.POSST',
 'QRLAC_31R6_QBS947.POSST',
 'QRLAB_23L7_QBS943.POSST',
 'QRLAG_13L3_QBS947.POSST',
 'QRLAC_31R1_QBS943.POSST',
 'QRLAB_11R5_QBS943.POSST',
 'QRLAB_11R8_QBS947.POSST',
 'QRLAA_29L8_QBS947.POSST',
 'QRLAA_17R5_QBS947.POSST',
 'QRLAB_23L6_QBS947.POSST',
 'QRLAA_29L5_QBS947.POSST',
 'QRLAA_17L7_QBS947.POSST',
 'QRLAC_31R8_QBS947.POSST',
 'QRLAA_29R2_QBS947.POSST',
 'QRLAA_21R2_QBS947.POSST',
 'QRLAB_19R2_QBS947.POSST',
 'QRLAB_11L1_QBS943.POSST',
 'QRLAB_19L1_QBS943.POSST',
 'QRLAA_29R7_QBS943.POSST',
 'QRLAB_15R3_QBS943.POSST',
 'QRLAB_19R1_QBS947.POSST',
 'QRLAA_29L6_QBS943.POSST',
 'QRLAA_29L4_QBS943.POSST',
 'QRLAC_31R2_QBS947.POSST',
 'QRLBA_09R5_QBS947.POSST',
 'QRLBA_09L2_QBS943.POSST',
 'QRLAB_31L7_QBS943.POSST',
 'QRLAA_17R1_QBS947.POSST',
 'QRLAB_27R2_QBS947.POSST',
 'QRLAB_23R3_QBS947.POSST',
 'QRLAA_17L5_QBS947.POSST',
 'QRLAA_25L7_QBS947.POSST',
 'QRLAB_23R6_QBS947.POSST',
 'QRLAB_23L2_QBS943.POSST',
 'QRLAA_13R6_QBS947.POSST',
 'QRLAA_13L1_QBS947.POSST',
 'QRLAA_21L2_QBS947.POSST',
 'QRLAA_33L8_QBS947.POSST',
 'QRLAC_31R8_QBS943.POSST',
 'QRLAB_31L7_QBS947.POSST',
 'QRLAB_23L8_QBS947.POSST',
 'QRLAB_11R1_QBS943.POSST',
 'QRLAB_27L7_QBS943.POSST',
 'QRLAA_17L8_QBS947.POSST',
 'QRLAA_21L3_QBS943.POSST',
 'QRLAA_25R7_QBS947.POSST',
 'QRLAA_21R1_QBS947.POSST',
 'QRLAA_33L4_QBS947.POSST',
 'QRLAB_11R5_QBS947.POSST',
 'QRLAB_11R2_QBS947.POSST',
 'QRLBA_09R3_QBS947.POSST',
 'QRLBA_09R7_QBS943.POSST',
 'QRLAB_19R3_QBS943.POSST',
 'QRLAB_23R5_QBS947.POSST',
 'QRLAA_17R3_QBS947.POSST',
 'QRLAA_13R1_QBS943.POSST',
 'QRLAA_33L5_QBS943.POSST',
 'QRLAB_31L1_QBS943.POSST',
 'QRLAA_25L7_QBS943.POSST',
 'QRLAD_33R6_QBS943.POSST',
 'QRLBB_09L3_QBS943.POSST',
 'QRLAB_23R7_QBS947.POSST',
 'QRLAB_23R5_QBS943.POSST',
 'QRLAA_21L7_QBS943.POSST',
 'QRLAB_27L3_QBS943.POSST',
 'QRLAC_31R6_QBS943.POSST',
 'QRLAB_27R5_QBS943.POSST',
 'QRLAB_11R8_QBS943.POSST',
 'QRLAA_25R6_QBS947.POSST',
 'QRLAB_23L7_QBS947.POSST',
 'QRLAB_23R2_QBS947.POSST',
 'QRLAB_11L8_QBS943.POSST',
 'QRLAA_21R8_QBS943.POSST',
 'QRLAB_19R8_QBS943.POSST',
 'QRLAB_11L5_QBS943.POSST',
 'QRLBA_09L5_QBS947.POSST',
 'QRLAA_17L4_QBS943.POSST',
 'QRLAH_15L3_QBS947.POSST',
 'QRLAC_31R3_QBS947.POSST',
 'QRLAA_13R2_QBS947.POSST',
 'QRLAB_27R1_QBS943.POSST',
 'QRLAA_21R5_QBS947.POSST',
 'QRLAB_19L6_QBS943.POSST',
 'QRLAD_33R4_QBS947.POSST',
 'QRLAA_33L3_QBS943.POSST',
 'QRLAB_23L5_QBS947.POSST',
 'QRLAA_13L2_QBS947.POSST',
 'QRLAB_19R4_QBS947.POSST',
 'QRLAB_31L5_QBS947.POSST',
 'QRLAB_31L3_QBS947.POSST',
 'QRLAA_21R4_QBS947.POSST',
 'QRLAB_27R4_QBS943.POSST',
 'QRLAA_25L6_QBS947.POSST',
 'QRLAB_19L2_QBS943.POSST',
 'QRLAA_29L1_QBS943.POSST',
 'QRLAA_29L2_QBS947.POSST',
 'QRLAA_29R4_QBS943.POSST',
 'QRLAA_25R5_QBS947.POSST',
 'QRLAB_19R2_QBS943.POSST',
 'QRLAB_27L1_QBS947.POSST',
 'QRLAA_29L8_QBS943.POSST',
 'QRLAB_23R2_QBS943.POSST',
 'QRLAA_13L1_QBS943.POSST',
 'QRLAA_17R5_QBS943.POSST',
 'QRLAA_25R8_QBS943.POSST',
 'QRLAA_21R7_QBS943.POSST',
 'QRLAA_17L6_QBS947.POSST',
 'QRLAA_29L7_QBS943.POSST',
 'QRLAB_27R4_QBS947.POSST',
 'QRLAB_15R1_QBS943.POSST',
 'QRLAA_21R3_QBS943.POSST',
 'QRLAA_21L7_QBS947.POSST',
 'QRLAA_21R5_QBS943.POSST',
 'QRLAB_23L3_QBS947.POSST',
 'QRLAA_25L6_QBS943.POSST',
 'QRLAA_17R3_QBS943.POSST',
 'QRLAB_11R3_QBS943.POSST',
 'QRLAB_15R6_QBS943.POSST',
 'QRLAB_27L4_QBS947.POSST',
 'QRLAA_17L8_QBS943.POSST',
 'QRLAF_25L8_QBS943.POSST',
 'QRLAB_31L8_QBS943.POSST',
 'QRLAA_21R2_QBS943.POSST',
 'QRLAB_27R8_QBS947.POSST',
 'QRLAA_29L5_QBS943.POSST',
 'QRLBA_09L8_QBS947.POSST',
 'QRLAB_19R5_QBS943.POSST',
 'QRLAD_33R2_QBS947.POSST',
 'QRLAB_23L6_QBS943.POSST',
 'QRLAB_27L3_QBS947.POSST',
 'QRLAA_25R1_QBS947.POSST',
 'QRLAB_23L3_QBS943.POSST',
 'QRLBA_09R1_QBS943.POSST',
 'QRLAC_31R2_QBS943.POSST',
 'QRLAA_33L8_QBS943.POSST',
 'QRLAB_23R4_QBS947.POSST',
 'QRLAA_29R7_QBS947.POSST',
 'QRLAB_31L6_QBS947.POSST',
 'QRLAA_13R5_QBS943.POSST',
 'QRLAB_19L8_QBS943.POSST',
 'QRLAB_15L2_QBS943.POSST',
 'QRLAB_11L5_QBS947.POSST',
 'QRLAA_29L7_QBS947.POSST',
 'QRLAB_19R7_QBS947.POSST',
 'QRLAB_27R6_QBS947.POSST',
 'QRLAA_33L2_QBS943.POSST',
 'QRLAA_29L3_QBS943.POSST',
 'QRLBA_09L4_QBS943.POSST',
 'QRLAA_13L8_QBS943.POSST']

#traditionally not dowloaded
other_varnames_static=[
    'QRLEA_06L3_QBS947.POSST',
    'QRLEA_06L7_QBS947.POSST',
    'QRLEA_06R3_QBS947.POSST',
    'QRLEA_06R7_QBS947.POSST',
    'QRLEB_05L4_QBS947.POSST',
    'QRLEB_05R4_QBS947.POSST',
    'QRLFE_05R4_QBS947.POSST',
    'QRLFF_05L4_QBS947.POSST'
]

arc_cells_by_sector = {
 'S12': [
 'QRLBA_09R1_QBS947.POSST',
 'QRLBA_09R1_QBS943.POSST',
 'QRLAB_11R1_QBS947.POSST',
 'QRLAB_11R1_QBS943.POSST',
 'QRLAA_13R1_QBS947.POSST',
 'QRLAA_13R1_QBS943.POSST',
 'QRLAB_15R1_QBS947.POSST',
 'QRLAB_15R1_QBS943.POSST',
 'QRLAA_17R1_QBS947.POSST',
 'QRLAA_17R1_QBS943.POSST',
 'QRLAB_19R1_QBS947.POSST',
 'QRLAB_19R1_QBS943.POSST',
 'QRLAA_21R1_QBS947.POSST',
 'QRLAA_21R1_QBS943.POSST',
 'QRLAB_23R1_QBS947.POSST',
 'QRLAB_23R1_QBS943.POSST',
 'QRLAA_25R1_QBS947.POSST',
 'QRLAA_25R1_QBS943.POSST',
 'QRLAB_27R1_QBS947.POSST',
 'QRLAB_27R1_QBS943.POSST',
 'QRLAA_29R1_QBS947.POSST',
 'QRLAA_29R1_QBS943.POSST',
 'QRLAC_31R1_QBS947.POSST',
 'QRLAC_31R1_QBS943.POSST',
 'QRLAD_33R1_QBS947.POSST',
 'QRLAD_33R1_QBS943.POSST',
 'QRLAA_33L2_QBS947.POSST',
 'QRLAA_33L2_QBS943.POSST',
 'QRLAB_31L2_QBS947.POSST',
 'QRLAB_31L2_QBS943.POSST',
 'QRLAA_29L2_QBS947.POSST',
 'QRLAA_29L2_QBS943.POSST',
 'QRLAB_27L2_QBS947.POSST',
 'QRLAB_27L2_QBS943.POSST',
 'QRLAA_25L2_QBS947.POSST',
 'QRLAA_25L2_QBS943.POSST',
 'QRLAB_23L2_QBS947.POSST',
 'QRLAB_23L2_QBS943.POSST',
 'QRLAA_21L2_QBS947.POSST',
 'QRLAA_21L2_QBS943.POSST',
 'QRLAB_19L2_QBS947.POSST',
 'QRLAB_19L2_QBS943.POSST',
 'QRLAA_17L2_QBS947.POSST',
 'QRLAA_17L2_QBS943.POSST',
 'QRLAB_15L2_QBS947.POSST',
 'QRLAB_15L2_QBS943.POSST',
 'QRLAA_13L2_QBS947.POSST',
 'QRLAA_13L2_QBS943.POSST',
 'QRLAB_11L2_QBS947.POSST',
 'QRLAB_11L2_QBS943.POSST',
 'QRLBA_09L2_QBS947.POSST',
 'QRLBA_09L2_QBS943.POSST'],

'S23':[
 'QRLBA_09R2_QBS947.POSST',
 'QRLBA_09R2_QBS943.POSST',
 'QRLAB_11R2_QBS947.POSST',
 'QRLAB_11R2_QBS943.POSST',
 'QRLAA_13R2_QBS947.POSST',
 'QRLAA_13R2_QBS943.POSST',
 'QRLAB_15R2_QBS947.POSST',
 'QRLAB_15R2_QBS943.POSST',
 'QRLAA_17R2_QBS947.POSST',
 'QRLAA_17R2_QBS943.POSST',
 'QRLAB_19R2_QBS947.POSST',
 'QRLAB_19R2_QBS943.POSST',
 'QRLAA_21R2_QBS947.POSST',
 'QRLAA_21R2_QBS943.POSST',
 'QRLAB_23R2_QBS947.POSST',
 'QRLAB_23R2_QBS943.POSST',
 'QRLAA_25R2_QBS947.POSST',
 'QRLAA_25R2_QBS943.POSST',
 'QRLAB_27R2_QBS947.POSST',
 'QRLAB_27R2_QBS943.POSST',
 'QRLAA_29R2_QBS947.POSST',
 'QRLAA_29R2_QBS943.POSST',
 'QRLAC_31R2_QBS947.POSST',
 'QRLAC_31R2_QBS943.POSST',
 'QRLAD_33R2_QBS947.POSST',
 'QRLAD_33R2_QBS943.POSST',
 'QRLAA_33L3_QBS947.POSST',
 'QRLAA_33L3_QBS943.POSST',
 'QRLAB_31L3_QBS947.POSST',
 'QRLAB_31L3_QBS943.POSST',
 'QRLAA_29L3_QBS947.POSST',
 'QRLAA_29L3_QBS943.POSST',
 'QRLAB_27L3_QBS947.POSST',
 'QRLAB_27L3_QBS943.POSST',
 'QRLAA_25L3_QBS947.POSST',
 'QRLAA_25L3_QBS943.POSST',
 'QRLAB_23L3_QBS947.POSST',
 'QRLAB_23L3_QBS943.POSST',
 'QRLAA_21L3_QBS947.POSST',
 'QRLAA_21L3_QBS943.POSST',
 'QRLAB_19L3_QBS947.POSST',
 'QRLAB_19L3_QBS943.POSST',
 'QRLAG_17L3_QBS947.POSST',
 'QRLAG_17L3_QBS943.POSST',
 'QRLAH_15L3_QBS947.POSST',
 'QRLAH_15L3_QBS943.POSST',
 'QRLAG_13L3_QBS947.POSST',
 'QRLAG_13L3_QBS943.POSST',
 'QRLAH_11L3_QBS947.POSST',
 'QRLAH_11L3_QBS943.POSST',
 'QRLBB_09L3_QBS947.POSST',
 'QRLBB_09L3_QBS943.POSST'],
 'S34':[
 'QRLBA_09R3_QBS947.POSST',
 'QRLBA_09R3_QBS943.POSST',
 'QRLAB_11R3_QBS947.POSST',
 'QRLAB_11R3_QBS943.POSST',
 'QRLAA_13R3_QBS947.POSST',
 'QRLAA_13R3_QBS943.POSST',
 'QRLAB_15R3_QBS947.POSST',
 'QRLAB_15R3_QBS943.POSST',
 'QRLAA_17R3_QBS947.POSST',
 'QRLAA_17R3_QBS943.POSST',
 'QRLAB_19R3_QBS947.POSST',
 'QRLAB_19R3_QBS943.POSST',
 'QRLAA_21R3_QBS947.POSST',
 'QRLAA_21R3_QBS943.POSST',
 'QRLAB_23R3_QBS947.POSST',
 'QRLAB_23R3_QBS943.POSST',
 'QRLAA_25R3_QBS947.POSST',
 'QRLAA_25R3_QBS943.POSST',
 'QRLAB_27R3_QBS947.POSST',
 'QRLAB_27R3_QBS943.POSST',
 'QRLAA_29R3_QBS947.POSST',
 'QRLAA_29R3_QBS943.POSST',
 'QRLAC_31R3_QBS947.POSST',
 'QRLAC_31R3_QBS943.POSST',
 'QRLAD_33R3_QBS947.POSST',
 'QRLAD_33R3_QBS943.POSST',
 'QRLAA_33L4_QBS947.POSST',
 'QRLAA_33L4_QBS943.POSST',
 'QRLAB_31L4_QBS947.POSST',
 'QRLAB_31L4_QBS943.POSST',
 'QRLAA_29L4_QBS947.POSST',
 'QRLAA_29L4_QBS943.POSST',
 'QRLAB_27L4_QBS947.POSST',
 'QRLAB_27L4_QBS943.POSST',
 'QRLAE_25L4_QBS947.POSST',
 'QRLAE_25L4_QBS943.POSST',
 'QRLAB_23L4_QBS947.POSST',
 'QRLAB_23L4_QBS943.POSST',
 'QRLAA_21L4_QBS947.POSST',
 'QRLAA_21L4_QBS943.POSST',
 'QRLAB_19L4_QBS947.POSST',
 'QRLAB_19L4_QBS943.POSST',
 'QRLAA_17L4_QBS947.POSST',
 'QRLAA_17L4_QBS943.POSST',
 'QRLAB_15L4_QBS947.POSST',
 'QRLAB_15L4_QBS943.POSST',
 'QRLAA_13L4_QBS947.POSST',
 'QRLAA_13L4_QBS943.POSST',
 'QRLAB_11L4_QBS947.POSST',
 'QRLAB_11L4_QBS943.POSST',
 'QRLBA_09L4_QBS947.POSST',
 'QRLBA_09L4_QBS943.POSST'],
 'S45':[
 'QRLBA_09R4_QBS947.POSST',
 'QRLBA_09R4_QBS943.POSST',
 'QRLAB_11R4_QBS947.POSST',
 'QRLAB_11R4_QBS943.POSST',
 'QRLAA_13R4_QBS947.POSST',
 'QRLAA_13R4_QBS943.POSST',
 'QRLAB_15R4_QBS947.POSST',
 'QRLAB_15R4_QBS943.POSST',
 'QRLAA_17R4_QBS947.POSST',
 'QRLAA_17R4_QBS943.POSST',
 'QRLAB_19R4_QBS947.POSST',
 'QRLAB_19R4_QBS943.POSST',
 'QRLAA_21R4_QBS947.POSST',
 'QRLAA_21R4_QBS943.POSST',
 'QRLAB_23R4_QBS947.POSST',
 'QRLAB_23R4_QBS943.POSST',
 'QRLAA_25R4_QBS947.POSST',
 'QRLAA_25R4_QBS943.POSST',
 'QRLAB_27R4_QBS947.POSST',
 'QRLAB_27R4_QBS943.POSST',
 'QRLAA_29R4_QBS947.POSST',
 'QRLAA_29R4_QBS943.POSST',
 'QRLAC_31R4_QBS947.POSST',
 'QRLAC_31R4_QBS943.POSST',
 'QRLAD_33R4_QBS947.POSST',
 'QRLAD_33R4_QBS943.POSST',
 'QRLAA_33L5_QBS947.POSST',
 'QRLAA_33L5_QBS943.POSST',
 'QRLAB_31L5_QBS947.POSST',
 'QRLAB_31L5_QBS943.POSST',
 'QRLAA_29L5_QBS947.POSST',
 'QRLAA_29L5_QBS943.POSST',
 'QRLAB_27L5_QBS947.POSST',
 'QRLAB_27L5_QBS943.POSST',
 'QRLAA_25L5_QBS947.POSST',
 'QRLAA_25L5_QBS943.POSST',
 'QRLAB_23L5_QBS947.POSST',
 'QRLAB_23L5_QBS943.POSST',
 'QRLAA_21L5_QBS947.POSST',
 'QRLAA_21L5_QBS943.POSST',
 'QRLAB_19L5_QBS947.POSST',
 'QRLAB_19L5_QBS943.POSST',
 'QRLAA_17L5_QBS947.POSST',
 'QRLAA_17L5_QBS943.POSST',
 'QRLAB_15L5_QBS947.POSST',
 'QRLAB_15L5_QBS943.POSST',
 'QRLAA_13L5_QBS947.POSST',
 'QRLAA_13L5_QBS943.POSST',
 'QRLAB_11L5_QBS947.POSST',
 'QRLAB_11L5_QBS943.POSST',
 'QRLBA_09L5_QBS947.POSST',
 'QRLBA_09L5_QBS943.POSST'],
 'S56':[
 'QRLBA_09R5_QBS947.POSST',
 'QRLBA_09R5_QBS943.POSST',
 'QRLAB_11R5_QBS947.POSST',
 'QRLAB_11R5_QBS943.POSST',
 'QRLAA_13R5_QBS947.POSST',
 'QRLAA_13R5_QBS943.POSST',
 'QRLAB_15R5_QBS947.POSST',
 'QRLAB_15R5_QBS943.POSST',
 'QRLAA_17R5_QBS947.POSST',
 'QRLAA_17R5_QBS943.POSST',
 'QRLAB_19R5_QBS947.POSST',
 'QRLAB_19R5_QBS943.POSST',
 'QRLAA_21R5_QBS947.POSST',
 'QRLAA_21R5_QBS943.POSST',
 'QRLAB_23R5_QBS947.POSST',
 'QRLAB_23R5_QBS943.POSST',
 'QRLAA_25R5_QBS947.POSST',
 'QRLAA_25R5_QBS943.POSST',
 'QRLAB_27R5_QBS947.POSST',
 'QRLAB_27R5_QBS943.POSST',
 'QRLAA_29R5_QBS947.POSST',
 'QRLAA_29R5_QBS943.POSST',
 'QRLAC_31R5_QBS947.POSST',
 'QRLAC_31R5_QBS943.POSST',
 'QRLAD_33R5_QBS947.POSST',
 'QRLAD_33R5_QBS943.POSST',
 'QRLAA_33L6_QBS947.POSST',
 'QRLAA_33L6_QBS943.POSST',
 'QRLAB_31L6_QBS947.POSST',
 'QRLAB_31L6_QBS943.POSST',
 'QRLAA_29L6_QBS947.POSST',
 'QRLAA_29L6_QBS943.POSST',
 'QRLAB_27L6_QBS947.POSST',
 'QRLAB_27L6_QBS943.POSST',
 'QRLAA_25L6_QBS947.POSST',
 'QRLAA_25L6_QBS943.POSST',
 'QRLAB_23L6_QBS947.POSST',
 'QRLAB_23L6_QBS943.POSST',
 'QRLAA_21L6_QBS947.POSST',
 'QRLAA_21L6_QBS943.POSST',
 'QRLAB_19L6_QBS947.POSST',
 'QRLAB_19L6_QBS943.POSST',
 'QRLAA_17L6_QBS947.POSST',
 'QRLAA_17L6_QBS943.POSST',
 'QRLAB_15L6_QBS947.POSST',
 'QRLAB_15L6_QBS943.POSST',
 'QRLAA_13L6_QBS947.POSST',
 'QRLAA_13L6_QBS943.POSST',
 'QRLAB_11L6_QBS947.POSST',
 'QRLAB_11L6_QBS943.POSST',
 'QRLBA_10L6_QBS943.POSST',
 'QRLBA_10L6_QBS947.POSST'],
 'S67':[
 'QRLBA_10R6_QBS943.POSST',
 'QRLBA_10R6_QBS947.POSST',
 'QRLAB_11R6_QBS947.POSST',
 'QRLAB_11R6_QBS943.POSST',
 'QRLAA_13R6_QBS947.POSST',
 'QRLAA_13R6_QBS943.POSST',
 'QRLAB_15R6_QBS947.POSST',
 'QRLAB_15R6_QBS943.POSST',
 'QRLAA_17R6_QBS947.POSST',
 'QRLAA_17R6_QBS943.POSST',
 'QRLAB_19R6_QBS947.POSST',
 'QRLAB_19R6_QBS943.POSST',
 'QRLAA_21R6_QBS947.POSST',
 'QRLAA_21R6_QBS943.POSST',
 'QRLAB_23R6_QBS947.POSST',
 'QRLAB_23R6_QBS943.POSST',
 'QRLAA_25R6_QBS947.POSST',
 'QRLAA_25R6_QBS943.POSST',
 'QRLAB_27R6_QBS947.POSST',
 'QRLAB_27R6_QBS943.POSST',
 'QRLAA_29R6_QBS947.POSST',
 'QRLAA_29R6_QBS943.POSST',
 'QRLAC_31R6_QBS947.POSST',
 'QRLAC_31R6_QBS943.POSST',
 'QRLAD_33R6_QBS947.POSST',
 'QRLAD_33R6_QBS943.POSST',
 'QRLAA_33L7_QBS947.POSST',
 'QRLAA_33L7_QBS943.POSST',
 'QRLAB_31L7_QBS947.POSST',
 'QRLAB_31L7_QBS943.POSST',
 'QRLAA_29L7_QBS947.POSST',
 'QRLAA_29L7_QBS943.POSST',
 'QRLAB_27L7_QBS947.POSST',
 'QRLAB_27L7_QBS943.POSST',
 'QRLAA_25L7_QBS947.POSST',
 'QRLAA_25L7_QBS943.POSST',
 'QRLAB_23L7_QBS947.POSST',
 'QRLAB_23L7_QBS943.POSST',
 'QRLAA_21L7_QBS947.POSST',
 'QRLAA_21L7_QBS943.POSST',
 'QRLAB_19L7_QBS947.POSST',
 'QRLAB_19L7_QBS943.POSST',
 'QRLAA_17L7_QBS947.POSST',
 'QRLAA_17L7_QBS943.POSST',
 'QRLAB_15L7_QBS947.POSST',
 'QRLAB_15L7_QBS943.POSST',
 'QRLAA_13L7_QBS947.POSST',
 'QRLAA_13L7_QBS943.POSST',
 'QRLAB_11L7_QBS947.POSST',
 'QRLAB_11L7_QBS943.POSST',
 'QRLBA_09L7_QBS947.POSST',
 'QRLBA_09L7_QBS943.POSST'],
 'S78':[
 'QRLBA_09R7_QBS947.POSST',
 'QRLBA_09R7_QBS943.POSST',
 'QRLAB_11R7_QBS947.POSST',
 'QRLAB_11R7_QBS943.POSST',
 'QRLAA_13R7_QBS947.POSST',
 'QRLAA_13R7_QBS943.POSST',
 'QRLAB_15R7_QBS947.POSST',
 'QRLAB_15R7_QBS943.POSST',
 'QRLAA_17R7_QBS947.POSST',
 'QRLAA_17R7_QBS943.POSST',
 'QRLAB_19R7_QBS947.POSST',
 'QRLAB_19R7_QBS943.POSST',
 'QRLAA_21R7_QBS947.POSST',
 'QRLAA_21R7_QBS943.POSST',
 'QRLAB_23R7_QBS947.POSST',
 'QRLAB_23R7_QBS943.POSST',
 'QRLAA_25R7_QBS947.POSST',
 'QRLAA_25R7_QBS943.POSST',
 'QRLAB_27R7_QBS947.POSST',
 'QRLAB_27R7_QBS943.POSST',
 'QRLAA_29R7_QBS947.POSST',
 'QRLAA_29R7_QBS943.POSST',
 'QRLAC_31R7_QBS947.POSST',
 'QRLAC_31R7_QBS943.POSST',
 'QRLAD_33R7_QBS947.POSST',
 'QRLAD_33R7_QBS943.POSST',
 'QRLAA_33L8_QBS947.POSST',
 'QRLAA_33L8_QBS943.POSST',
 'QRLAB_31L8_QBS947.POSST',
 'QRLAB_31L8_QBS943.POSST',
 'QRLAA_29L8_QBS947.POSST',
 'QRLAA_29L8_QBS943.POSST',
 'QRLAB_27L8_QBS947.POSST',
 'QRLAB_27L8_QBS943.POSST',
 'QRLAF_25L8_QBS947.POSST',
 'QRLAF_25L8_QBS943.POSST',
 'QRLAB_23L8_QBS947.POSST',
 'QRLAB_23L8_QBS943.POSST',
 'QRLAA_21L8_QBS947.POSST',
 'QRLAA_21L8_QBS943.POSST',
 'QRLAB_19L8_QBS947.POSST',
 'QRLAB_19L8_QBS943.POSST',
 'QRLAA_17L8_QBS947.POSST',
 'QRLAA_17L8_QBS943.POSST',
 'QRLAB_15L8_QBS947.POSST',
 'QRLAB_15L8_QBS943.POSST',
 'QRLAA_13L8_QBS947.POSST',
 'QRLAA_13L8_QBS943.POSST',
 'QRLAB_11L8_QBS947.POSST',
 'QRLAB_11L8_QBS943.POSST',
 'QRLBA_09L8_QBS947.POSST',
 'QRLBA_09L8_QBS943.POSST'],
 'S81':[
 'QRLBA_09R8_QBS947.POSST',
 'QRLBA_09R8_QBS943.POSST',
 'QRLAB_11R8_QBS947.POSST',
 'QRLAB_11R8_QBS943.POSST',
 'QRLAA_13R8_QBS947.POSST',
 'QRLAA_13R8_QBS943.POSST',
 'QRLAB_15R8_QBS947.POSST',
 'QRLAB_15R8_QBS943.POSST',
 'QRLAA_17R8_QBS947.POSST',
 'QRLAA_17R8_QBS943.POSST',
 'QRLAB_19R8_QBS947.POSST',
 'QRLAB_19R8_QBS943.POSST',
 'QRLAA_21R8_QBS947.POSST',
 'QRLAA_21R8_QBS943.POSST',
 'QRLAB_23R8_QBS947.POSST',
 'QRLAB_23R8_QBS943.POSST',
 'QRLAA_25R8_QBS947.POSST',
 'QRLAA_25R8_QBS943.POSST',
 'QRLAB_27R8_QBS947.POSST',
 'QRLAB_27R8_QBS943.POSST',
 'QRLAA_29R8_QBS947.POSST',
 'QRLAA_29R8_QBS943.POSST',
 'QRLAC_31R8_QBS947.POSST',
 'QRLAC_31R8_QBS943.POSST',
 'QRLAD_33R8_QBS947.POSST',
 'QRLAD_33R8_QBS943.POSST',
 'QRLAA_33L1_QBS947.POSST',
 'QRLAA_33L1_QBS943.POSST',
 'QRLAB_31L1_QBS947.POSST',
 'QRLAB_31L1_QBS943.POSST',
 'QRLAA_29L1_QBS947.POSST',
 'QRLAA_29L1_QBS943.POSST',
 'QRLAB_27L1_QBS947.POSST',
 'QRLAB_27L1_QBS943.POSST',
 'QRLAA_25L1_QBS947.POSST',
 'QRLAA_25L1_QBS943.POSST',
 'QRLAB_23L1_QBS947.POSST',
 'QRLAB_23L1_QBS943.POSST',
 'QRLAA_21L1_QBS947.POSST',
 'QRLAA_21L1_QBS943.POSST',
 'QRLAB_19L1_QBS947.POSST',
 'QRLAB_19L1_QBS943.POSST',
 'QRLAA_17L1_QBS947.POSST',
 'QRLAA_17L1_QBS943.POSST',
 'QRLAB_15L1_QBS947.POSST',
 'QRLAB_15L1_QBS943.POSST',
 'QRLAA_13L1_QBS947.POSST',
 'QRLAA_13L1_QBS943.POSST',
 'QRLAB_11L1_QBS947.POSST',
 'QRLAB_11L1_QBS943.POSST',
 'QRLBA_09L1_QBS947.POSST',
 'QRLBA_09L1_QBS943.POSST']}


