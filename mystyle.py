from matplotlib import rc, rcdefaults
import pylab as pl
from colorsys import hsv_to_rgb
import numpy as np
import matplotlib.pyplot as plt


def mystyle(fontsz=16):
	font = {#'family' : 'normal',
			#'weight' : 'bold',
			'size'   : fontsz}
#	print fontsz
	rcdefaults()
	rc('font', **font)

def mystyle_arial(fontsz=16, dist_tick_lab=10):

	rcdefaults()
	rc('font',**{'family':'sans-serif','sans-serif':['arial'], 'size':fontsz})
	rc(('xtick.major','xtick.minor','ytick.major','ytick.minor'), pad=dist_tick_lab)

def sciy():
	pl.gca().ticklabel_format(style='sci', scilimits=(0,0),axis='y')

def scix():
	pl.gca().ticklabel_format(style='sci', scilimits=(0,0),axis='x')

def colorprog(i_prog, Nplots, v1 = .9, v2 = 1.):
    if hasattr(Nplots, '__len__'):
        Nplots = len(Nplots)
    return hsv_to_rgb(float(i_prog)/float(Nplots), v1, v2)
    #return [pl.cm.rainbow(k) for k in np.linspace(0, 1, Nplots)][i_prog]

def comb_legend(sp1, sp2, *args, **kwargs):
    """
    Combine legends for twinx()ed subplots
    """
    lines, labels = sp1.get_legend_handles_labels()
    lines2, labels2 = sp2.get_legend_handles_labels()
    sp2.legend(lines + lines2, labels + labels2, *args, **kwargs)

def figure(title, figs=None, figsize=(8*1.5, 6*1.5), **kwargs):
    fig = plt.figure(figsize=figsize, **kwargs)
    fig.canvas.set_window_title(title)
    if figs != None:
        figs.append(fig)
    fig.patch.set_facecolor('w')
    plt.suptitle(title, fontsize=20)
    return fig

