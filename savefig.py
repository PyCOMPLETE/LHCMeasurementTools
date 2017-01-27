import sys
import os
import re

import matplotlib.pyplot as plt

pdijksta_dir = '/afs/cern.ch/work/l/lhcscrub/pdijksta_PyECLOUD_benchmark/plots/'
re_fill = re.compile('^\d{4}$')
re_script = re.compile('^(\d{3}[a-z]?_.{4})')

def get_file_title(fig, title=None):
    if title == None:
        title = fig.canvas.get_window_title()
    script_title = os.path.basename(sys.argv[0])
    info = re_script.match(script_title)
    if info != None:
        script_number = info.group(1)
    else:
        print('Warning! Script number could not be identified! %s' % script_title)
        script_number = script_title[:5]

    for arg in sys.argv[1:]:
        info = re_fill.match(arg)
        if info != None:
            filln = arg
            break
    else:
        print('No filln found')
        filln = '0'
    
    return '_'.join([script_number,filln,title])



def pdijksta(fig, title=None):
    file_title = get_file_title(fig, title)
    save_path = pdijksta_dir + file_title

    fig.savefig(save_path, dpi=200)
