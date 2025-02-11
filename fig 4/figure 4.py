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

#Packages used for curve fitting
import lmfit as lm
from lmfit.models import GaussianModel, ConstantModel, SkewedGaussianModel
from lmfit import Model




def time_to_dist(t):
    c =299792458
    d = t*1e-9*c/1.462
    return d

def dist_to_time(d, mod = 1.462):
    c = 299792458
    t = d*mod/c*1e9
    return t
def modified_siegert(t, amp, center, tc, fsr):

    return 1 + amp*(np.e**(-np.abs((t-center))/tc) * (1+np.cos(2*np.pi*fsr/1000*(t-center))))**2

def lm_curve_fit(xdata, ydata):
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

    return center

def save_data(time, coin, file, folder_name = 'analysed files'):
    file = Path(file)
    save_folder = file.parent.joinpath(folder_name)
    if not os.path.isdir(save_folder):
        os.makedirs(save_folder)

    save_name = save_folder.joinpath(file.stem + '.txt')
    np.savetxt(save_name, np.transpose([time, coin]))

def time_to_distance(t):
    c = 299_792_458
    return 1000*(t*1e-12)*c/2 #returns distance in mm

def dist_to_time(d):
    c = 299792458
    d = 2*d/1000
    t = d/c*1e12
    return t


"""
Setup the plots
"""
#Plot the results
fig1, (ax0, ax1, ax2) = plt.subplots(3, figsize = (6, 5.5))

#---- Some plotting variables ----
labelsize = 15
titlesize = 18
linethickness = 3
ticksize = 12
hfont = 'Arial'



#---- Free-space measurements ----
#250205
labels = np.array([25.0, 23.75, 22.5, 21.25, 20.0, 18.75, 17.5, 16.25, 15.0 ,13.75, 12.5, 11.25, 10.0, 8.75, 7.5, 6.25, 5.0, 3.75, 2.5, 1.25])
files = [f'C:/Users/staff/Documents/Papers TBP/Thermal lidar/repository/fig 4/FS resolution/FS_target_{i}mm_72mA_drive_det2_14dB_det1_1kcts_det2_17uA_det1_13uA_300s_250205.txt' for i in labels]


x, y1, y2 = [], [], []
for i, file in enumerate(files):
    delta_t, coin = np.loadtxt(file, unpack = True)  #second_correleation(file, eta_engine,bins, binsize)
    g2 = coin / np.mean(coin)

    td = np.abs(delta_t[np.where(g2 == np.max(g2))[0][0]])
    td_fit = np.abs(lm_curve_fit(delta_t, g2))


    print(f'Dist: {labels[i]}: Peak time delay: {td}, Curve fit time delay: {td_fit}')

    if i == 0:
        td0 = td
        td0_fit = td_fit

    #dist = time_to_distance(np.abs(td))

    x.append(25-labels[i])
    y1.append(td - td0)
    y2.append(td_fit - td0_fit)


#Plot the results
ax2.scatter(x, y1, s = 50, c = '#58508d', label = 'Measurement')
ax2.scatter(x, y2, s = 50, c =  '#ef8250', label = 'Curve Fit', marker = 'x')

#---- Add theoretical line ----
distances = np.flip(labels-1.25)
times = dist_to_time(distances)
ax2.plot(distances, times, c='k', linestyle = 'dashed', label = 'Expected')
ax2.legend(loc='upper left')

#--- Work out some resolutions ---
print('\n---- Freespace measurements ----\n')
print(f'\nDifference between theoretical and measured')
diff_1 = np.array(y1)-np.array(times)
#print(diff_1)
print(f'RMS: {np.mean(np.sqrt(diff_1**2))}\n')

print('Differences between theoretical and curvefit')
diff_2 = np.array(y2)-np.array(times)
#print(diff_2)
print(f'RMS: {np.mean(np.sqrt(diff_2**2))}\n')

#-- Ax1 --
ax2.set_xlabel('Distance [mm]', fontsize = labelsize, fontname = hfont)
ax2.set_ylabel(r'Time delay [ps]', fontsize = labelsize, fontname = hfont)
ax2.set_title('Free-space delay', fontsize = titlesize, fontname = hfont)
ax2.set_ylim(-45,220)
ax2.set_xlim(-2,26)

ax2.tick_params(labelsize = ticksize)
ax2.minorticks_on()
ax2.minorticks_on()
ax2.legend(loc='upper left')



#----- Variable fiber delay measurements ----
#241030
labels = np.arange(0, 51, 2.5)

files = [f'C:/Users/staff/Documents/Papers TBP/Thermal lidar/repository/fig 4/Fiber resolution/{i}mm_delay_det2_16.35uA_det1_13uA_20s_241030.txt' for i in labels]

x, y1, y2 = [], [], []
for i, file in enumerate(files):
    delta_t, coin = np.loadtxt(file, unpack = True)  #second_correleation(file, eta_engine,bins, binsize)
    g2 = coin / np.mean(coin)

    td = np.abs(delta_t[np.where(g2 == np.max(g2))[0][0]])
    td_fit = np.abs(lm_curve_fit(delta_t, g2))


    print(f'Dist: {labels[i]}: Peak time delay: {td}, Curve fit time delay: {td_fit}')

    if i == 0:
        td0 = td
        td0_fit = td_fit


    x.append(labels[i])
    y1.append(td - td0)
    y2.append(td_fit - td0_fit)

ax1.scatter(x, y1, s = 50, c = '#58508d', label = 'Measurement')
ax1.scatter(x, y2, s = 50, c =  '#ef8250', label = 'Curve Fit', marker = 'x')

#---- Add theoretical line ----
distances = labels
times = dist_to_time(distances)
ax1.plot(distances, times, c='k', linestyle = 'dashed', label = 'Expected')
ax1.legend(loc='upper left')


#--- Work out some resolutions ---
print('\n---- Fiber measurements ----\n')
print(f'Difference between theoretical and measured')
diff_1 = np.array(y1)-np.array(times)
#(diff_1)
print(f'RMS: {np.mean(np.sqrt(diff_1**2))}\n')

print('Differences between theoretical and curvefit')
diff_2 = np.array(y2)-np.array(times)
#print(diff_2)
print(f'RMS: {np.mean(np.sqrt(diff_2**2))}\n')


#-- Ax2 --
ax1.set_xlabel('Distance [mm]', fontsize = labelsize, fontname = hfont)
ax1.set_ylabel(r'Time delay [ps]', fontsize = labelsize, fontname = hfont)
ax1.set_title('Variable fiber delay', fontsize = titlesize, fontname = hfont)
ax1.set_ylim(-45,410)
ax1.set_xlim(-2,54)

ax1.tick_params(labelsize = ticksize)
ax1.minorticks_on()
ax1.minorticks_on()
ax1.legend(loc='upper left')

#----- Long fiber delay measurements ----
#241009
timetag_file1 = 'C:/Users/staff/Documents/Papers TBP/Thermal lidar/repository/fig 4/fiber delays/Delay_10m_1548.93nm_80mA_-26dB_det2_16.35uA_det1_13uA_60s_241009.txt'
timetag_file2 = 'C:/Users/staff/Documents/Papers TBP/Thermal lidar/repository/fig 4/fiber delays/Delay_100m#2_1548.93nm_80mA_-23dB_det2_16.35uA_det1_13uA_60s_241009.txt'
timetag_file3 = 'C:/Users/staff/Documents/Papers TBP/Thermal lidar/repository/fig 4/fiber delays/Delay_1000m#2_1548.93nm_80mA_-26dB_det2_16.35uA_det1_13uA_60s_241009.txt'
labels = ['10 m','100 m', '1 km']
files = [timetag_file1, timetag_file2, timetag_file3]
z = [100, 200, 0]

col = [ '#58508d', '#ef8250', '#63b179']
for i, file in enumerate(files):
    delta_t, coin = np.loadtxt(file, unpack = True)

    g2 = coin/np.mean(coin)

    ax0.plot(delta_t, np.flip(g2), label = labels[i], zorder = z[i], c = col[i])



#-- Ax 0 --

ax0.set_xlabel('Time [ns]', fontsize = labelsize, fontname = hfont)
ax0.set_ylabel(r'$g^2$($\tau$)', fontsize = labelsize, fontname = hfont)
ax0.set_xlim(-10,7000)
ax0.legend(loc = 'upper right')

ax0.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}".replace(",", " ")))
ax0.tick_params(labelsize = ticksize)
ax0.minorticks_on()

# Create a secondary x-axis on the top, which is a transformed version of the bottom
ax1t = ax0.secondary_xaxis('top', functions=(time_to_dist, dist_to_time))
ax1t.set_xlabel('Fiber length [m]', fontsize = titlesize, labelpad = 10, fontname = hfont)


#---- Finish the plots ----
#Add labels
x = 0
y = 1.1
ax0.text(x,y+0.2, '(a)', transform=ax0.transAxes, fontsize=labelsize)
ax1.text(x,y, '(b)', transform=ax1.transAxes, fontsize=labelsize)
ax2.text(x,y, '(c)', transform=ax2.transAxes, fontsize=labelsize)



#---- Don't forget to plot ----
plt.subplots_adjust(bottom = 0.134, hspace= 0.82, right = 0.79)
plt.show()
