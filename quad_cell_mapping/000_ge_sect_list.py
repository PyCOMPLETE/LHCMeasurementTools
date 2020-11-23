

# Q34s found in the machine
# MQ.34R1 MQ.34R2 MQ.34R3 MQ.34R4 MQ.34R5 MQ.34R6 MQ.34R7 MQ.34R8'''
# They are all on the R side!


ref_strings = \
'''Q12R1:QRLAA_13R1_CV947
Q12L2:QRLAB_11L2_CV947
Q12R2:QRLAA_13R2_CV947
Q12L3:QRLAH_11L3_CV947
Q12R3:QRLAA_13R3_CV947
Q12L4:QRLAB_11L4_CV947
Q12R4:QRLAA_13R4_CV947
Q12L5:QRLAB_11L5_CV947 
Q12R5:QRLAA_13R5_CV947
Q12L6:QRLAB_11L6_CV947
Q12R6:QRLAA_13R6_CV947
Q12L7:QRLAB_11L7_CV947
Q12R7:QRLAA_13R7_CV947
Q12L8:QRLAB_11L8_CV947
Q12R8:QRLAA_13R8_CV947
Q12L1:QRLAB_11L1_CV947'''.split()

#let's start

# I assume that 947 comes always before 943. Is that true???

Qmin = 9
Qmax = 33
check_limits = lambda q : (q>=Qmin and q<=Qmax)

for ref_str in ref_strings:
    side = ref_str.split(':')[0][-2:]
    Qref = int(ref_str[1:].split(side)[0])
    cell_ref = int(ref_str.split(':')[-1].split(side)[0].split('_')[-1])

    with open('cell_naming_%s.txt'%side, 'w') as fid:
        for ii in range(-(Qref-8), (34-Qref)+1, 2):
            if side[0]=='L':
                if check_limits(Qref+ii-1): fid.write('Q%d%s\t%d%s_CV943\n'%(Qref+ii-1, side, cell_ref+ii, side))
            if check_limits(Qref+ii): fid.write('Q%d%s\t%d%s_CV947\n'%(Qref+ii, side, cell_ref+ii, side))
            if side[0]=='R':
                if check_limits(Qref+ii+1): fid.write('Q%d%s\t%d%s_CV943\n'%(Qref+ii+1, side, cell_ref+ii, side))
        if side[0]=='L':
            # handle Q34
            tmp = int(side[1])-1
            if tmp==0: tmp=8
            fid.write('Q%dR%d\t%d%s_CV947\n'%(34, tmp,  33, side))
