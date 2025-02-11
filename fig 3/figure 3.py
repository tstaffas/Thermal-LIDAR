# ------IMPORTS-----
# Packages for ETA backend
import json
import etabackend.eta as eta  # Available at: https://github.com/timetag/ETA, https://eta.readthedocs.io/en/latest/
import etabackend.tk as etatk

import os
import time as t
from pathlib import Path

from matplotlib import pyplot as plt
import numpy as np
from scipy.constants import c

#Packages used for curve fitting
import lmfit as lm
from lmfit.models import GaussianModel, ConstantModel, SkewedGaussianModel
from lmfit import Model



#---- Curve fit ----
def siegert(t, amp, center, tc):
    r = 1 + amp*(np.e**(-np.abs((t-center))/tc))**2
    return r

def lm_curve_fit_siegert(xdata, ydata):
    gmodel = Model(siegert)

    #Make some initial values:
    amp = (np.max(ydata)-1)
    center = xdata[np.where(ydata == np.max(ydata))[0][0]]
    tc = 100

    params = gmodel.make_params(amp = amp, center = center, tc = tc)

    x_eval = np.linspace(xdata[0], xdata[-1], 2000)
    y_eval = gmodel.eval(params, t = x_eval)

    result = gmodel.fit(ydata, params, t = xdata)

    center = result.params['center'].value

    if False:
        fig, ax = plt.subplots()
        ax.plot(xdata, ydata, label = 'data')
        ax.plot(xdata, result.best_fit, label='fit')

        print(center)
        plt.show()

    return result.best_fit


def modified_siegert(t, amp, center, tc, fsr):

    return 1 + amp*(np.e**(-np.abs((t-center))/tc) * (1+np.cos(2*np.pi*fsr/1000*(t-center))))**2

def lm_curve_fit_mod_siegert(xdata, ydata):
    gmodel = Model(modified_siegert)

    #Make some initial values:
    amp = (np.max(ydata)-1)/4
    center = xdata[np.where(ydata == np.max(ydata))[0][0]]
    fsr = 2
    tc = 1000

    params = gmodel.make_params(amp = amp, center = center, tc = tc, fsr = fsr)

    x_eval = np.linspace(xdata[0], xdata[-1], 2000)
    y_eval = gmodel.eval(params, t = x_eval)

    result = gmodel.fit(ydata, params, t = xdata)

    center = result.params['center'].value
    fsr = result.params['fsr'].value
    #print(f'FSR: {fsr}')

    if False:
        fig, ax = plt.subplots()
        ax.plot(xdata, ydata, label = 'data')
        ax.plot(xdata, result.best_fit, label='fit')

        print(center)
        plt.show()

    return result.best_fit




#Data files
file1 = 'single_0_delay_det2_16.35uA_det1_13uA_40s_241030.txt'
file2 = 'multi_0_delay_det2_16.35uA_det1_13uA_20s_241030.txt'

DT_1, coin_1 = np.loadtxt(file1, unpack = True)
g2_1 = coin_1/np.mean(coin_1)

DT_2, coin_2 = np.loadtxt(file2, unpack=True)
g2_2 = coin_2/np.mean(coin_2)

#Curve fit
curve = lm_curve_fit_siegert(DT_1, g2_1)
curve2 = lm_curve_fit_mod_siegert(DT_2, g2_2)


#Plot the results
fig1, ((ax1,ax2))= plt.subplots(2,1)

#---- Some variable ----
labelsize = 17
titlesize = 18
linethickness = 3
ticksize = 12
hfont = 'Arial'


#Set the first column
ax1.plot((DT_1)*1e-3, g2_1, c = '#58508d', label = 'Measurement')
ax1.plot((DT_1)*1e-3, curve, c = '#ef8250', label = 'Siegert')
ax1.legend(loc='center', bbox_to_anchor = (0.2, 0.60))

ax1.set_xlabel('Time [ns]', fontsize = labelsize, fontname = hfont)
ax1.set_ylabel(r'$g^2$($\tau$)', fontsize = labelsize, fontname = hfont)
ax1.set_title('Single longitudinal mode', fontsize = titlesize, fontname = hfont)


#ax2.grid()
ax1.set_xlim(-5,5)
ax1.set_ylim(0.95, 2)
ax1.tick_params(labelsize = ticksize)
ax1.minorticks_on()

#Set the second column
ax2.plot((DT_2)*1e-3, g2_2, c = '#58508d', label = 'Measurement')
ax2.plot((DT_2)*1e-3, curve2, c = '#ef8250', label = 'Modified\nSiegert')
ax2.legend(loc='center', bbox_to_anchor = (0.2, 0.60))

ax2.set_xlabel('Time [ns]', fontsize = labelsize, fontname = hfont)
ax2.set_ylabel(r'$g^2$($\tau$)', fontsize = labelsize, fontname = hfont)
ax2.set_title('Two longitudinal modes', fontsize = titlesize, fontname = hfont)

#ax3.set_title(r'$\theta_{max}$', fontsize = titlesize)
#ax3.grid()
ax2.set_ylim(0.95, 2)
ax2.set_xlim(-5,5)
ax2.tick_params(labelsize = ticksize)
ax2.minorticks_on()

#---- Customize the plot ----

#Add labels
x = -0.1
y = 1.15
ax1.text(x,y, '(a)', transform=ax1.transAxes, fontsize=labelsize, fontname = hfont)
ax1.text(0.5,0.8, r'$\theta_{half}$', transform=ax1.transAxes,
         horizontalalignment='center', fontsize=labelsize, fontname = hfont)

ax2.text(x,y, '(b)', transform=ax2.transAxes, fontsize=labelsize, fontname = hfont)
ax2.text(0.5,0.8, r'$\theta_{max}$', transform=ax2.transAxes,
         horizontalalignment='center', fontsize=labelsize, fontname = hfont)

#--- Create second y axis ---
ax12 = ax1.twinx()
ax12.minorticks_on()
ax12.set_ylim(900, 2000)
ax12.tick_params(labelsize = ticksize)
ax12.text(1.17, 0.0, 'Coincidences', transform=ax12.transAxes, fontsize=labelsize, fontname = hfont, rotation = 90, rotation_mode = 'anchor')

ax22 = ax2.twinx()
ax22.minorticks_on()
ax22.set_ylim(200, 500)
ax22.tick_params(labelsize = ticksize)
ax22.text(1.17, 0.0, 'Coincidences', transform=ax22.transAxes, fontsize=labelsize, fontname = hfont, rotation = 90, rotation_mode = 'anchor')



plt.subplots_adjust(hspace=0.768, wspace=0.55, bottom = 0.124, right = 0.84)
plt.show()


