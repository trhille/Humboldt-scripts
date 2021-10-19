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

rhoi = 910.0
rhosw = 1028.0
print("** Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = OptionParser(description=__doc__)
parser.add_option("-r", dest="runs", help="path to .nc file or dir containing output.nc file (strings separated by commas; no spaces)", default=None, metavar="FILENAME")
parser.add_option("-t", dest="timeLevels", help="integer time levels at which to plot (int separated by commas; no spaces)", default=-1, metavar="FILENAME")

options, args = parser.parse_args()
#if type(options.runs) is list: 
runs = options.runs.split(',') # split run directories into list
#else:
#    runs = [options.runs] # Kind of cludgey, but allows for consistent indexing in loop
    
#if type(options.timeLevels) is list:
timeLevs = options.timeLevels.split(',') # split time levels into list
#else:
#    timeLevs = [options.timeLevels]

# Define cellMask bit values so we don't have to use bitmask conversion script
initialExtentValue = 1
dynamicValue = 2
floatValue = 4
groundingLineValue = 256

# Use rows for q values, columns for climate forcings.
# Currently this means plotting all time levels on the same
# axes, which may or may not work, depending on application.
# This could definitely be improved.
fig, axs = plt.subplots(2, 3, sharex=True, sharey=True) 


# Plot bed topo and initial extent on all axes using first file
bedPlot = []
initExtentPlot = []
f = Dataset(runs[0], 'r')
xCell = f.variables["xCell"][:]
yCell = f.variables["yCell"][:]
bed = f.variables["bedTopography"][:]
cellMask = f.variables["cellMask"][:]
initialExtentMask = (cellMask & initialExtentValue) // initialExtentValue

for ax in axs.ravel(): 
    bedPlot.append(ax.tricontourf(xCell, yCell, bed[0,:], levels=100, cmap='BrBG'))
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
        lineColor = 'tab:grey'
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

    
    # Loop over time levels
    for timeLev in timeLevs:
        timeLev = int(timeLev) #these are strings for some reason; make int to index
        axs[row,col].tricontour(xCell, yCell, groundingLineMask[timeLev, :], 
           levels=[0.9999], colors=get_line_color(run), linestyles='solid')
        

        
# Customize plots
axs[0,0].set_ylim(top=-1020000, bottom = -1120000) # hard-code limits for now
axs[0,0].set_xlim(left=-425000, right=-300000)

plt.show()