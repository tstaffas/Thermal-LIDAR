import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle

def gauss(t, c, b = 0.039, y0 =-70, a = 1):
    sigma = b/2.355
    return a*np.e**(-(t-c)**2/(2*sigma**2))

def exp_gauss(t, c, b = 0.029, y0 =-70, a = 1):
    sigma = b/2.355
    return y0 - a*(t-c)**2/(2*sigma**2)

#Backmatter
file_source = 'Laser.CSV'
file_etalon_high = 'Filtered.CSV'

WL1, source = np.loadtxt(file_source, skiprows=39, delimiter=',', unpack = True)
WL4, etalon_HP = np.loadtxt(file_etalon_high, skiprows=39, delimiter=',', unpack = True)


#---- Create a plot ----
fig, (ax1, ax2) = plt.subplots(2,1)
ax1.plot(WL1, source, c = '#58508d', label = 'Thermal source')
#ax1.plot(WL2, BPF, label = 'BPF')
ax2.plot(WL4, etalon_HP, c = '#58508d')

WL1 = np.linspace(np.min(WL1), np.max(WL1), 100000)
#-- Second y axis for filter transmission
ax12 = ax1.twinx()
#---- Plot BPF ----
b = gauss(WL1, c=1549.55, b=1.5, y0 = 0, a = 1)
ax12.plot(WL1, b, c = 'k', linestyle = 'solid', label = 'BPF')

#---- Plot the FSR ----
center = 1549.55
e1 = gauss(WL1, c =center, y0 = 0)
e2 = gauss(WL1, c = center-1.579, y0 = 0)
e3 = gauss(WL1, c = center+1.579,  y0 = 0)

ax12.plot(WL1, e1, linestyle = 'dashed', c = 'k', label = 'Etalon')
ax12.plot(WL1, e2, linestyle = 'dashed', c = 'k')
ax12.plot(WL1, e3, linestyle = 'dashed', c = 'k')


ax1.annotate(text='', xy = (center,-49), xytext = (center-1.579, -49), arrowprops=dict(arrowstyle='<->'))

#---- Some plotting variables ----
labelsize = 15
titlesize = 20
linethickness = 3
ticksize = 12


#---- Customize the plot ----
#ax1.legend(loc = 'upper left')

hfont = 'Arial'
#ax.set_title('Laser below threshold', fontsize = titlesize)
ax1.set_xlabel('Wavelength [nm]', fontsize = labelsize, fontname = hfont)
ax1.set_ylabel('Intensity [dBm]', fontsize = labelsize, fontname = hfont)
ax1.set_title('Unfiltered thermal spectrum', fontsize = titlesize, fontname = hfont)
ax1.minorticks_on()
ax1.tick_params(labelsize = ticksize)
ax1.set_ylim(-72,-43)
ax1.set_xlim(1546.9, 1552.1)

ax12.set_ylabel('Transmission [arb.u]', fontsize = labelsize, fontname = hfont)
ax12.minorticks_on()
ax12.tick_params(labelsize=0, color = 'w')

ax2.set_xlabel('Wavelength [nm]', fontsize = labelsize, fontname = hfont)
ax2.set_ylabel('Intensity [dBm]', fontsize = labelsize, fontname = hfont)
ax2.set_title('Filtered thermal spectrum', fontsize = titlesize, fontname = hfont)
ax2.minorticks_on()
ax2.tick_params(labelsize = ticksize)
ax2.set_ylim(-70,-43)
ax2.set_xlim(1546.9, 1552.1)
#ax2.grid()

#Add labels
x = 0.00
y = 1.05
ax1.text(x,y, '(a)', transform=ax1.transAxes, fontsize=labelsize, fontname = hfont)
ax1.text(0.217, 0.83, 'Etalon FSR', transform=ax1.transAxes, fontsize=labelsize, fontname = hfont)
ax2.text(x,y, '(b)', transform=ax2.transAxes, fontsize=labelsize, fontname = hfont)


plt.subplots_adjust(hspace=0.835, wspace=0.26, left = 0.15, top = 0.915, bottom = 0.141)
plt.show()
