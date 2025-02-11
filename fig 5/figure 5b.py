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

# ------- ETA analysis ----------
def load_eta(recipe, **kwargs):
    print('Loading ETA')
    with open(recipe, 'r') as filehandle:
        recipe_obj = json.load(filehandle)

    eta_engine = eta.ETA()
    eta_engine.load_recipe(recipe_obj)

    # Set parameters in the recipe
    for arg in kwargs:
        eta_engine.recipe.set_parameter(arg, str(kwargs[arg]))

    eta_engine.load_recipe()

    return eta_engine

def eta_analysis(file, eta_engine):
    print('Starting ETA analysis')
    print(f'{file=}')
    cut = eta_engine.clips(Path(file), format=1)
    result = eta_engine.run({"timetagger1": cut}, group='quTAG')
    print('Finished ETA analysis')

    return result

def second_correleation(timetag_file, eta_engine, bins, binsisize):
    # ETA analys
    result = eta_analysis(timetag_file, eta_engine)

    # extract result
    hist1 = result["h3"]
    hist2 = result["h4"]
    hist0 = result["h4_zero"]
    hist1[0] += hist0[0]

    coin = np.concatenate((hist2[::-1], hist1))

    delta_t = np.arange(-bins, bins) * binsize * 1e-3

    g2 = coin / np.mean(coin)

    print(f'Peak coincidence: {np.max(coin)} \nTotal detected photon pairs: {np.sum(coin)}\n Noise lvl: {np.mean(coin)}\n')

    return coin, g2, delta_t


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


#Data analysis
eta_recipe = 'Correlation-swabian.eta'

#Data file
timetag_file5 = 'kista_13.5km_10km_10km_delay_det2_16.35uA_det1_13uA_30s_241030.timeres'

#ETA Settings
binsize = 20
bins = 350000
delay = int(1.95e8-2.8e6-0.7e3) + int(132.0e6) + int(99.8841e6 - 55.82e3 ) 


time_axis = np.arange(0,bins)*binsize
eta_engine = load_eta(eta_recipe, bins=bins, binsize=binsize, delay_1 = delay)
coin, g2, delta_t = second_correleation(timetag_file5, eta_engine,bins, binsize)


# ---- Plot the results ----
long_time = delta_t/1e3 + delay/1e6
ax2.plot(long_time, g2, c = '#58508d', label = 'Deployed fiber')

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
ax2_zoom.plot((delta_t+delay/1000)/1000, g2, c = '#58508d')
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

