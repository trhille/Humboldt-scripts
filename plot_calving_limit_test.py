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

runsList = [calv1kmyr, calv3kmyr, calv5kmyr, calv7kmyr, noCalvLimit]
runsLabels = ['1 km yr$^{-1}$', '3 km yr$^{-1}$', 
              '5 km yr$^{-1}$', '7 km yr$^{-1}$', 'no limit']
fig, ax = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(4,8))

#set colormap for plots, and remove red to be colorblind safe
colorList = ['tab:blue', 'tab:orange', 'tab:green', 'tab:cyan', 'tab:purple']

# Hard-code these for plotting ease
VAF_orderOfMag = 13
floatingIceArea_orderOfMag = 8


# loop through runs and plot
for color, run, label in zip(colorList, runsList, runsLabels):
    ax[0].plot(run.variables["daysSinceStart"][:] / 365 + 2007., 
               run.variables["volumeAboveFloatation"][:] / 10**VAF_orderOfMag,
               color=color, label=label)
    ax[1].plot(run.variables["daysSinceStart"][:] / 365 + 2007., 
               run.variables["floatingIceArea"][:] / 10**floatingIceArea_orderOfMag,
               color=color, label=label)
    ax[2].plot(run.variables["daysSinceStart"][:] / 365 + 2007., 
               run.variables["surfaceSpeedMax"][:], color=color, label=label)
    
for axes in ax:
    axes.grid('on')
    
ax[0].set_ylabel('Volume above\nfloatation (10$^{' + 
                 str(VAF_orderOfMag) +'}$ m$^3$)', fontsize=12)
ax[1].set_ylabel('Floating ice area\n(10$^{' + 
                 str(floatingIceArea_orderOfMag) +'}$ m$^2$)', fontsize=12)
ax[2].set_ylabel('Maximum ice\nspeed (m yr$^{-1}$)', fontsize=12)

ax[2].set_ylim(5.e2, 5.e4)
ax[2].set_yscale('log')
ax[2].set_xlabel('Year', fontsize=12)

ax[0].legend()

fig.savefig('calvingRateTest', dpi=400, bbox_inches='tight', format='pdf')