#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 15:40:48 2021

@author: trevorhillebrand
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import os
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
from optparse import OptionParser
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.pyplot import cm
import matplotlib.tri as tri
import matplotlib.gridspec as gridspec
from matplotlib.colorbar import Colorbar
from matplotlib.colors import TwoSlopeNorm

rhoi = 910.0
rhosw = 1028.0
print("** Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = OptionParser(description=__doc__)
parser.add_option("-o", dest="obs", help="path to .nc file containing observations (strings separated by commas; no spaces)", default=None, metavar="FILENAME")
parser.add_option("-r", dest="runs", help="path to .nc file or dir containing output.nc file (strings separated by commas; no spaces)", default=None, metavar="FILENAME")
parser.add_option("-t", dest="timeLevels", help="integer time levels at which to plot (int separated by commas; no spaces)", default=-1)
parser.add_option("-b", dest="bedTopo", help="path to gridded bed topography data product for plotting",
                  default='/global/cfs/cdirs/piscees/GIS/BedMachineGreenland-2021-04-20_Humboldt.nc', metavar="FILENAME")
options, args = parser.parse_args()
runs = options.runs.split(',') # split run directories into list
runs.insert(0, options.obs)
bedTopoFile = options.bedTopo
timeLev = options.timeLevels.split(',')[0]  # split time levels into list
initialExtentValue = 1
dynamicValue = 2
floatValue = 4
groundingLineValue = 256

# Get fields from observations files 
obs = Dataset(options.obs, 'r')
xCell = obs.variables["xCell"][:] / 1000.
yCell = obs.variables["yCell"][:] / 1000.
obsSpeed = obs.variables["surfaceSpeed"][:] * 3.154e7  # m/yr

# Use rows for q values, columns for climate forcings.
# Currently this means plotting all time levels on the same
# axes, which may or may not work, depending on application.
# This could definitely be improved.
#fig, axs = plt.subplots(2, 3, sharex=True, sharey=True, 
#                        constrained_layout=True)
#nRows = 2
#nCols = len(runs) // nRows + 1  # in case you want to hard-code nRows
nCols = 5
nRows = len(runs) // nCols
# Add another row if necessary
#if nRows * (nCols) < len(runs)-1:
#   nRows += 1

figSize = (3 * nCols + 1, 3 * nRows + 1)
fig, axs_array = plt.subplots(nrows=nRows, ncols=nCols,
                       sharex=True, sharey=True,
                       figsize=figSize)
#                       constrained_layout=True)
#                       gridspec_kw={'height_ratios': (nRows - 1) * [6] + [1]})
axs = axs_array.ravel()
#fig = plt.figure(figsize=(3 * nCols + 1, 3 * nRows + 1))
# last column is for colorbars
gs = gridspec.GridSpec(nRows, nCols,
                       height_ratios=[1] * nRows,
                       width_ratios=(nCols-1)*[1] + [0.1]) 
#axs = []
#for row in np.arange(0, nRows):
#    for col in np.arange(0, nCols-1):
#        axs.append(plt.subplot(gs[row, col]))
##        axs[-1].sharex(axs[0])
##        axs[-1].sharey(axs[0])
#cbarAx1 = plt.subplot(gs[0,-1])
#cbarAx2 = plt.subplot(gs[1,-1])

# Plot bed topo and initial extent on all axes using first file
bedPlot = []
bedTopo = Dataset(bedTopoFile, 'r')
bed = bedTopo.variables["topg"][:]
bedx = bedTopo.variables["x1"][:] / 1000.
bedy = bedTopo.variables["y1"][:] / 1000.
bedX, bedY = np.meshgrid(bedx, bedy)
bedTopo.close()

initExtentPlot = []
if '.nc' not in runs[0]:
    runs[0] = runs[0] + '/output.nc'
f = Dataset(runs[0], 'r')
xCell = f.variables["xCell"][:] / 1000.
yCell = f.variables["yCell"][:] / 1000.
#cellMask = f.variables["cellMask"][:]
#initialExtentMask = (cellMask & initialExtentValue) // initialExtentValue

bedMin = -420.
bedMax = 780.
for ii,ax in enumerate(axs):
    if ii < len(runs):
        bedPlot.append(ax.pcolormesh(bedX, bedY, bed, cmap='BrBG',
                       norm=TwoSlopeNorm(
                       vmin=bedMin, vmax=bedMax, vcenter=0.)))
    else:
        ax.axis("off")
    #initExtentPlot.append(ax.tricontour(xCell, yCell, initialExtentMask[0,:], colors='black'))
f.close()
    
#Loop over runs
spdPlot = []
for ii,run in enumerate(runs):
    if '.nc' not in run:
        run = run + '/output.nc'
    f = Dataset(run, 'r')
    thk = f.variables["thickness"][:]
    spd = f.variables["surfaceSpeed"][:] * 3.154e7  # m/yr
    spdNan = spd  #add Nans below
    
    #obs are the first item in runs list
    if ii==0:
        timeLevPlot = 0
        obsSpdMask1 = spdNan[timeLevPlot, :] >= 100
        obsSpdMask2 = spdNan[timeLevPlot, :] >= 300
        obsSpdMask3 = spdNan[timeLevPlot, :] >= 600

    else:
        timeLevPlot = int(timeLev) #these are strings for some reason; make int to index
        
    if "cellMask" in f.variables.keys():
        cellMask = f.variables["cellMask"][:]
        floatMask = (cellMask & floatValue) // floatValue
        dynamicMask = (cellMask & dynamicValue) // dynamicValue
        # groundingLineMask = (cellMask & groundingLineValue) // groundingLineValue
        spdNan[timeLevPlot, 1 - dynamicMask] = np.nan

    thkMask = thk[timeLevPlot,:] <= 1.
    thk[timeLevPlot, thkMask] = np.nan
    spdNan[timeLevPlot, thkMask] = np.nan
    spdMask1 = spdNan[timeLevPlot, :] >= 100
    spdMask2 = spdNan[timeLevPlot, :] >= 300
    spdMask3 = spdNan[timeLevPlot, :] >= 600

    if ii==0:
        obsVel = spdNan[timeLevPlot,:]
        plotField = obsVel
        cmap='plasma'
        vmin=0.
        vmax=1000.
    else:
        plotField = spdNan[timeLevPlot,:] - obsVel
        cmap='RdBu'
        vmin=-350.
        vmax=350.

    spdPlot.append(axs[ii].tripcolor(xCell, yCell, 
                   plotField, cmap=cmap, shading='gouraud',
                   vmin=(vmin), vmax=(vmax)))
    for spdMask, obsSpdMask in zip([spdMask1, spdMask2, spdMask3],
                                   [obsSpdMask1, obsSpdMask2, obsSpdMask3]):
        axs[ii].tricontour(xCell, yCell, spdMask,
                           levels=[0.9999], colors='cyan', linewidths=[0.5])
        axs[ii].tricontour(xCell, yCell, obsSpdMask,
                           levels=[0.9999], colors='black', linewidths=[0.5])

    axs[ii].set_aspect('equal')

# Customize plots
for ax in axs:
    ax.set_ylim(top=-1020, bottom = -1120) # hard-code limits for now
    ax.set_xlim(left=-425, right=-300)

#cbar1 = Colorbar(ax=cbarAx1, mappable=spdPlot[0], orientation='vertical',
#                 label="Ice speed (m yr$^{-1}$)", extend='both')
#cbar2 = Colorbar(ax=cbarAx2, mappable=bedPlot[0], orientation='vertical',
#                 label="Bed elevation (m)")

cbar1 = plt.colorbar(ax=axs_array[-1,0], location='bottom',
                     mappable=spdPlot[0], orientation='horizontal',
                     label="Observed speed 2017-2018 (m yr$^{-1}$)")
cbar2 = plt.colorbar(ax=axs_array[-1,2], location='bottom',
                     mappable=spdPlot[1], orientation='horizontal',
                     label="Modeled - observed speed (m yr$^{-1}$)", extend='both')
cbar3 = plt.colorbar(ax=axs_array[-1,-1], mappable=bedPlot[0],
                     orientation='horizontal', location='bottom',
                     label="Bed elevation (m)", extend='both')

plt.show()
#fig.savefig('basalFrictionExpTuning', dpi=400, bbox_inches='tight')
