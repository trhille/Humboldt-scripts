#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 13:43:07 2021

@author: trevorhillebrand
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

dirPath = '/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_1to10km/r04_Results/speedLimitTest/'
calv1kmyr = Dataset(dirPath + 'm7_MIROC5_VM170_shelfMelt20myr_1kmyrCalvLimit_globalStats.nc', 'r')
calv3kmyr = Dataset(dirPath + 'm7_MIROC5_VM170_shelfMelt20myr_3kmyrCalvLimit_globalStats.nc', 'r')
calv5kmyr = Dataset(dirPath + 'm7_MIROC5_VM170_shelfMelt20myr_5kmyrCalvLimit_globalStats.nc', 'r')
calv7kmyr = Dataset(dirPath + 'm7_MIROC5_VM170_shelfMelt20myr_7kmyrCalvLimit_globalStats.nc', 'r')
noCalvLimit = Dataset(dirPath + 'm7_MIROC5_VM170_shelfMelt20myr_noCalvLimit_globalStats.nc', 'r')
control = Dataset('/Users/trevorhillebrand/Documents/mpas/MALI_output/' +
                  'Humboldt_1to10km/r04_Results/m7/m7_control_MIROC5_globalStats.nc', 'r')

runsList = [calv1kmyr, calv3kmyr, calv5kmyr, calv7kmyr, noCalvLimit]
runsLabels = ['1 km yr$^{-1}$', '3 km yr$^{-1}$', 
              '5 km yr$^{-1}$', '7 km yr$^{-1}$', 'no limit']
fig, ax = plt.subplots(nrows=3, ncols=1, sharex=True)

#set colormap for plots, and remove red to be colorblind safe
colorList = ['tab:orange', 'tab:blue', 'tab:green', 'tab:olive', 'tab:pink']

# Hard-code these for plotting ease
VAF_orderOfMag = 12
floatingIceArea_orderOfMag = 8
rhoi = 910.
rhosw = 1028.

yrsVAF = control.variables["daysSinceStart"][:] / 365 + 2007.
controlVAF = control.variables["volumeAboveFloatation"][:]

def VAF2seaLevel(vol):
    return -vol / 3.62e14 * rhoi / rhosw * 1000.

def seaLevel2VAF(vol):
    return -vol * 3.62e14 * rhosw / rhoi / 1000.

def addSeaLevAx(axName):
    seaLevAx = axName.secondary_yaxis('right', functions=(VAF2seaLevel, seaLevel2VAF))
    seaLevAx.set_ylabel('Sea-level\nequivalent (mm)', fontsize=16)

# loop through runs and plot
for color, run, label in zip(colorList, runsList, runsLabels):
    yr = run.variables["daysSinceStart"][:] / 365 + 2007.

    controlInterp = np.interp(yr, yrsVAF, controlVAF)
    ax[0].plot(yr, run.variables["volumeAboveFloatation"][:] - controlInterp,
               color=color, label=label)
    ax[1].plot(yr, run.variables["floatingIceArea"][:] / 10**floatingIceArea_orderOfMag,
               color=color, label=label)
    ax[2].plot(yr, run.variables["surfaceSpeedMax"][:], color=color, label=label)
    
for axes in ax:
    axes.grid('on')

addSeaLevAx(ax[0])
    
ax[0].set_ylabel('Total change in volume above\nfloatation (10$^{' +
                 str(VAF_orderOfMag) +'}$ m$^3$)')
ax[1].set_ylabel('Floating ice area\n(10$^{' + 
                 str(floatingIceArea_orderOfMag) +'}$ m$^2$)')
ax[2].set_ylabel('Maximum ice\nspeed (m yr$^{-1}$)')

ax[2].set_ylim(5.e2, 5.e4)
ax[2].set_yscale('log')
ax[2].set_xlabel('Year')

ax[0].legend()

fig.set_size_inches(6, 15)
fig.subplots_adjust(hspace=0.05)
plt.show()
#fig.savefig('calvingRateTest', dpi=400, bbox_inches='tight', format='pdf')
