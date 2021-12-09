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

rhoi = 910.0
rhosw = 1028.0
print("** Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = OptionParser(description=__doc__)
parser.add_option("-r", dest="runs", help="path to .nc file or dir containing output.nc file (strings separated by commas; no spaces)", default=None, metavar="FILENAME")
parser.add_option("-t", dest="timeLevels", help="integer time levels at which to plot (int separated by commas; no spaces)", default=-1, metavar="FILENAME")
parser.add_option("-b", dest="bedTopo", help="path to gridded bed topography data product for plotting", default='/global/cfs/cdirs/piscees/GIS/BedMachineGreenland-2021-04-20_Humboldt.nc', metavar="FILENAME")
options, args = parser.parse_args()
runs = options.runs.split(',') # split run directories into list
bedTopoFile = options.bedTopo
timeLev = options.timeLevels.split(',')[0]  # split time levels into list
initialExtentValue = 1
dynamicValue = 2
floatValue = 4
groundingLineValue = 256

# Hack for asymmetric colorbap. Matplotlib is dumb.
# set the colormap and centre the colorbar
class MidpointNormalize(mpl.colors.Normalize):
    """Normalise the colorbar."""
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

# Use rows for q values, columns for climate forcings.
# Currently this means plotting all time levels on the same
# axes, which may or may not work, depending on application.
# This could definitely be improved.
#fig, axs = plt.subplots(2, 3, sharex=True, sharey=True, 
#                        constrained_layout=True)
fig = plt.figure(figsize=(15,7))
nRows = 2
nCols = 4
# last column is for colorbars
gs = gridspec.GridSpec(nRows, nCols,
                       height_ratios=[1,1],
                       width_ratios=[1,1,1,0.1]) 
axs = []
for row in np.arange(0, nRows):
    for col in np.arange(0, nCols-1):
        axs.append(plt.subplot(gs[row, col]))
#        axs[-1].sharex(axs[0])
#        axs[-1].sharey(axs[0])
cbarAx1 = plt.subplot(gs[0,-1])
cbarAx2 = plt.subplot(gs[1,-1])

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
cellMask = f.variables["cellMask"][:]
initialExtentMask = (cellMask & initialExtentValue) // initialExtentValue

bedMin = -420.
bedMax = 780.
for ax in axs: 
    bedPlot.append(ax.pcolormesh(bedX, bedY, bed, cmap='BrBG',
                   vmin=bedMin, vmax=bedMax,
                   norm=MidpointNormalize(bedMin, bedMax, 0.)))
    initExtentPlot.append(ax.tricontour(xCell, yCell, initialExtentMask[0,:], colors='black'))
f.close()

# Functions to simplify plotting
def get_line_color(run):
    if 'm5_' in run or 'm5/' in run:
        lowCalving = 'VM180'
        medCalving = 'VM170'
        highCalving = 'VM160'
    elif 'm7_' in run or 'm7/' in run:
        if 'HadGEM2' in run or 'CNRM' in run:
            lowCalving = 'VM190'
            medCalving = 'VM180'
            highCalving = 'VM170'
        elif 'MIROC5' in run:
            lowCalving = 'VM180'
            medCalving = 'VM170'
            highCalving = 'VM160'
    
    if '2017calvingFront' in run or 'calvingVelocityData' in run:
        lineColor = 'lightgrey'
    elif highCalving in run:
        lineColor = 'tab:purple'
    elif medCalving in run:
        lineColor = 'tab:blue'
    elif lowCalving in run:
        lineColor = 'tab:cyan'
    else:
        lineColor = 'white'
    
    return lineColor

def get_line_style(run):
    if 'm5_' in run or 'm5/' in run:
        lineStyle = 'solid'
    elif 'm7_' in run or 'm7/' in run:
        lineStyle = 'dashed'
    else:
        lineStyle = 'none'
        
    return lineStyle
    
#Loop over runs
for run in runs:
    if '.nc' not in run:
        run = run + '/output.nc'
    f = Dataset(run, 'r')
    thk = f.variables["thickness"][:]
    spd = f.variables["surfaceSpeed"][:] * 3.154e7  # m/yr
    spdNan = spd  #add Nans below
    cellMask = f.variables["cellMask"][:]
    floatMask = (cellMask & floatValue) // floatValue
    dynamicMask = (cellMask & dynamicValue) // dynamicValue
    groundingLineMask = (cellMask & groundingLineValue) // groundingLineValue
    
    # Assign to subplot based on climate forcing and value of m 
    if 'MIROC5' in run:
        col = 0
    elif 'HadGEM2' in run:
        col = 1
    elif 'CNRM' in run:
        col = 2
    else:
        col = int(input('Climate forcing not found in ' + run + '. Enter column number between 0 and 2: '))
    
    if 'm5_' in run or 'm5/' in run:
        row = 0
    elif 'm7_' in run or 'm7/' in run:
        row = 1
    else:
        row = int(input('Basal friction exponent not found in ' + run + '. Enter row number between 0 and 1: '))

    index = row * (nCols - 1) + col
    spdPlot = []
    # Loop over time levels
    timeLev = int(timeLev) #these are strings for some reason; make int to index
    thkMask = thk[timeLev,:] <= 1.
    thk[timeLev, thkMask] = np.nan
    spdNan[timeLev, thkMask] = np.nan
    spdMask = spd[timeLev, :] >= 3.e3
    triang = tri.Triangulation(xCell, yCell)
    faceColors = thk[timeLev,:] * thkMask
    axs[index].tricontour(xCell, yCell, groundingLineMask[timeLev, :], 
                          levels=[0.9999], colors=get_line_color(run), linestyles='solid')
    spdPlot.append(axs[index].tripcolor(xCell, yCell, 
                   np.log10(spdNan[timeLev,:]), cmap='plasma', shading='flat',
                   vmin=np.log10(20.), vmax=np.log10(2.e3)))
    axs[index].tricontour(xCell, yCell, spdMask, levels=[0.9999], colors='cyan')
    axs[index].set_aspect('equal')

        
# Customize plots
axs[0].set_title('MIROC5')
axs[1].set_title('HadGEM2')
axs[2].set_title('CNRM-CM6')
for ax in axs:
    ax.set_ylim(top=-1020, bottom = -1120) # hard-code limits for now
    ax.set_xlim(left=-425, right=-300)

for ind in [0,3]:
    axs[ind].set_ylabel('km')
for ind in [0, 1, 2]:
    axs[ind].set_xticklabels(['','', '', '', '', ''])
for ind in [1, 2, 4, 5]:
    axs[ind].set_yticklabels(['','', '', '', '', ''])
for ind in [3,4,5]:
    axs[ind].set_xlabel('km')

cbar1 = Colorbar(ax=cbarAx1, mappable=spdPlot[0], orientation='vertical',
                 label="log ice speed\n(10$^x$ m yr$^{-1}$)")
cbar2 = Colorbar(ax=cbarAx2, mappable=bedPlot[0], orientation='vertical',
                 label="Bed elevation (m)")

plt.show()
#fig.savefig('speed2100_highCalving', dpi=400, bbox_inches='tight')
