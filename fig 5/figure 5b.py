# ------IMPORTS-----
# Packages for ETA backend
import json
import etabackend.eta as eta  # Available at: https://github.com/timetag/ETA, https://eta.readthedocs.io/en/latest/
import etabackend.tk as etatk

import os
import time as t
from pathlib import Path

from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from scipy.constants import c

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

def time_to_dist(t):
    c =299792458
    d = t*1e-9*c/1.462
    return d

def dist_to_time(d):
    c = 299792458
    t = d*1.462/c*1e9
    return t



"""
Setup the plots
"""
#Plot the results
fig1, ax2= plt.subplots()

#---- Some variable ----
labelsize = 17
titlesize = 19
linethickness = 3
ticksize = 12
hfont = 'Arial'


#Data file
data_file = 'kista_13.5km_10km_10km_delay_det2_16.35uA_det1_13uA_30s_241030.txt'

#measurement Settings
delta_t, coin = np.loadtxt(data_file, unpack = True)
g2 = coin/np.mean(coin)

# ---- Plot the results ----
ax2.plot(delta_t/1e6, g2, c = '#58508d', label = 'Deployed fiber')

ax2.legend(loc='upper left')
ax2.set_xlabel('Time [$\mu$s]', fontsize = labelsize, fontname = hfont)
ax2.set_ylabel(r'$g^2$($\tau$)', fontsize = labelsize, fontname = hfont)
#ax2.grid()

xlim = (423.9, 424.6)
ax2.set_xlim(xlim[0], xlim[1])
ax2.set_ylim(0.85, 1.9)
ax2.minorticks_on()

#Create a secondary y-axis
ax12 = ax2.twinx()
ax12.minorticks_on()
ax12.set_ylim(860, 1860)
ax12.tick_params(labelsize = ticksize)
ax12.text(1.12, 0.35, 'Coincidences', transform=ax12.transAxes, fontsize=labelsize, fontname = hfont, rotation = 90, rotation_mode = 'anchor')

#ax12.plot(long_time, coin)

# Create a secondary x-axis on the top, which is a transformed version of the bottom
ax2t = ax2.secondary_xaxis('top', functions=(time_to_dist, dist_to_time))

#---- Make a zoomed in plot ----
ax2_zoom = ax2.inset_axes([0.4, 0.4, 0.5, 0.55])
ax2_zoom.plot((delta_t)/1e6, g2, c = '#58508d')
ax2_zoom.set_xlim(424.022, 424.032 )
ax2_zoom.set_ylim(0.9, 1.55)
ax2_zoom.tick_params(left = True, right = False , labelleft = True ,labelbottom = True, bottom = True)
#ax2_zoom.grid()
ax2.indicate_inset_zoom(ax2_zoom)

#Create a secondary y-axis
ax12 = ax2_zoom.twinx()
ax12.minorticks_on()
ax12.set_ylim(1000, 1620)
ax12.tick_params(labelsize = ticksize)
#ax12.text(1.17, 0.0, 'Coincidences', transform=ax12.transAxes, fontsize=labelsize, fontname = hfont, rotation = 90, rotation_mode = 'anchor')



#Add labels
x = -0
y = 1.12
ax2.text(x,y, '(a)', transform=ax2.transAxes, fontsize=labelsize)
ax2.text(0.5,y, 'Fiber length [km]', transform=ax2.transAxes,
         horizontalalignment='center', fontsize=titlesize, fontname = hfont)

plt.subplots_adjust(hspace= 0.75, bottom = 0.18, right = 0.69, top = 0.78)
plt.show()

