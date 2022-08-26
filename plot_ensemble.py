#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 21:10:36 2020
Plot ensemble volume and sea level. Plots a desired variable from globalStats.nc
files within a given set of directories. Runs in the same directory are plotted
in different colors (up to 9; change colormap if more are needed) with the same
linestyle. Additional directories will loop through the same colors with a 
different linestyle. Thus, this script is useful for separating out runs by a 
single characteristic (linestlye). For example, RCP2.6 runs with dashed lines,
RCP8.5 runs with solid lines. Also plots globa mean sea level equivalent 
for volumeAboveFloatation.
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
parser.add_option("-d", dest="ensembleDirs", help="directory containing ensemble members (strings separated by commas; no spaces)", metavar="FILENAME")
parser.add_option("-b", dest="boundsDirs", help="directory containing ensemble members (strings separated by commas; no spaces)", metavar="FILENAME")
parser.add_option("-v", dest="variableName", help="variable(s) to plot, separated by commas", default = "volumeAboveFloatation")
parser.add_option("-c", dest="subtractControl", help="subtract control run values from ensemble members", default = 'False')

parser.add_option("-u", dest="units", help="units for mass/volume: m3, kg, Gt", default="m3")
options, args = parser.parse_args()

ensembleDirs = options.ensembleDirs.split(',') # split ensemble directories into list
variableName = options.variableName.split(',')

if options.boundsDirs:
    boundsDirs = options.boundsDirs.split(',') # split control files into list
else:
    boundsDirs = None

print("Using ice density of {} kg/m3".format(rhoi))
print("Using seawater density of {} kg/m3".format(rhosw))

#set colormap for plots
colormap = mpl.colors.TABLEAU_COLORS
#colormap = mpl.colors.XKCD_COLORS
colorlist = list(colormap.items())
#Get units
units = options.units

#set linestyles to loop through (e.g., to separate out RCP scenarios)
linestyleList = ['solid', 'dashed', 'dotted', 'dashdot']
linestyleIndex = 0 # initialize for loop

# create axes to plot into
nCols = 3
varFig, varAx = plt.subplots(len(variableName), nCols, sharey='row', sharex=True)
varAx = varAx.ravel()
#ratioFig, ratioAx = plt.subplots(1,1)
for ax in varAx:
    ax.grid('on')

plotLines = [] #empty list to fill with lines
plotLineNames = [] #empty list to fill with filenames
plotBounds = []
plotBoundNames = []

def VAF2seaLevel(vol):
    return -vol / 3.62e14 * rhoi / rhosw * 1000.

def seaLevel2VAF(vol):
    return -vol * 3.62e14 * rhosw / rhoi / 1000.

def plotEnsembleBounds(boundsDir):
    boundsFiles = sorted(os.listdir(boundsDir)) # get filenames in directory
    
    for boundsFile in boundsFiles:
        if 'globalStats' not in boundsFile:
            boundsFiles.remove(boundsFile)
            
            
    f1 = Dataset(boundsDir + boundsFiles[0], 'r')
    yr = f1.variables['daysSinceStart'][:] / 365.0
    f2 = Dataset(boundsDir + boundsFiles[1], 'r')
    if controlFile:
        #interpolate control run onto ensemble member time vector
        controlData = Dataset(controlFile, 'r')
        controlInterp = np.interp(yr, controlData.variables['daysSinceStart'][:]/365.0, 
                       controlData.variables[options.variableName][:])
    
    var2plot1 = f1.variables[options.variableName][:] \
                 - f1.variables[options.variableName][0] \
                 - controlInterp + controlInterp[0]
                
    var2plot2 = np.interp(yr, f2.variables['daysSinceStart'][:]/365.0,
                          f2.variables[options.variableName][:] 
                          - f2.variables[options.variableName][0]) \
                          - controlInterp + controlInterp[0]
    if 'HadGEM2' in boundsFiles[0]:
        plotAx = varAx[1]
    elif 'MIROC5' in boundsFiles[0]:
        plotAx = varAx[0]
        
    tmpFill = plotAx.fill_between(yr+2007., var2plot1, var2plot2, facecolor='tab:grey', alpha = 0.6)
    plotBounds.append(tmpFill)
    plotBoundNames.append(boundsFiles[0])
    
    return plotBounds, plotBoundNames

def plotEnsemble(ensDir, row, variable):
#    print("Reading and plotting file: {}".format(fname))
    colorIndex = 0 #initialize index to loop through color list
    ensembleFiles = sorted(os.listdir(ensDir)) # get filenames in directory
    for ensembleMember in ensembleFiles:
        if 'globalStats' in ensembleMember and 'control' not in ensembleMember:           
            f = Dataset(ensDir+ensembleMember,'r')
            yr = f.variables['daysSinceStart'][:]/365.0
            deltat = f.variables['deltat'][:] / 3.154e7 # convert deltat to years
            
            # get units
            try: 
                units = f.variables[variable].units
            except:
                units = 'm^3'
        
            var2plot = f.variables[variable][:] \
                         #- f.variables[options.variableName][0]

            # find index for plotting and get control file
            if 'CNRM' in ensembleMember:
                col = 2
                controlFile = [i for i in ensembleFiles if 'control_CNRM' in i][0]
            elif 'HadGEM2' in ensembleMember:
                col = 1 
                controlFile = [i for i in ensembleFiles if 'control_HadGEM2' in i][0]
            elif 'MIROC5' in ensembleMember:
                col = 0                     
                controlFile = [i for i in ensembleFiles if 'control_MIROC5' in i][0]

            # subtract off variables from control run
            if options.subtractControl == 'True':
                #interpolate control run onto ensemble member time vector
                controlData = Dataset(ensDir+controlFile, 'r')
                controlInterp = np.interp(yr, controlData.variables['daysSinceStart'][:]/365.0, 
                               controlData.variables[variable][:])

                var2plot = var2plot - controlInterp #+ controlInterp[0]
                
            plotInd = row * nCols + col
            
            plotAx = varAx[plotInd]
                
            tmpLine, = plotAx.plot(yr+2007., var2plot, 
                                   label=ensembleMember)
            plotLines.append(tmpLine)
            plotLineNames.append(ensembleMember)
            
            colorIndex += 1 # go to next color
            f.close()
            if variable == 'volumeAboveFloatation':
                print(ensembleMember + ': {:10.1f} mm SLR at 2100'.format(VAF2seaLevel(var2plot)[-1]))

    return plotLines, plotLineNames
    

def addSeaLevAx(axName):
    seaLevAx = axName.secondary_yaxis('right', functions=(VAF2seaLevel, seaLevel2VAF))
    seaLevAx.set_ylabel('Sea-level\nequivalent (mm)', fontsize=16)

controlIndex=0
boundsIndex=0
controlFile=None

for directory in ensembleDirs:
    print("Ensemble {}".format(directory))
    for row,variable in enumerate(variableName):

        # This makes labeling flexible for one or more rows.
        yLabelInd = row * nCols
        
        if boundsDirs and boundsIndex <= len(boundsDirs)-1:
            plotEnsembleBounds(boundsDirs[boundsIndex])
        plotLines, plotLineNames = plotEnsemble(directory, row, variable)
        if variable == "volumeAboveFloatation":
            addSeaLevAx(varAx[(row + 1) * nCols - 1])
        if variable == 'volumeAboveFloatation':            
            varAx[yLabelInd].set_ylabel('Total change in volume above\nfloatation (10$^{12}$ m$^3$)', fontsize=16)
        elif variable == 'floatingIceVolume':
            varAx[yLabelInd].set_ylabel('Total floating ice volume (10$^{11}$ m$^3$)', fontsize=16)
        else:
            varAx[yLabelInd].set_ylabel('Total {} (${}$)'.format(variable, units), fontsize=16)

    linestyleIndex += 1
    boundsIndex += 1

#make this flexible for one or more rows
    varAx[-3].set_xlabel('Year', fontsize=16)
    varAx[0].set_title('MIROC5')
    varAx[-2].set_xlabel('Year', fontsize=16)
    varAx[1].set_title('HadGEM2')
    varAx[-1].set_xlabel('Year', fontsize=16)
    varAx[2].set_title('CNRM-CM6')


#varAx.legend()
varFig.tight_layout()
#varAx.set_ylim(bottom=-7e12, top=0)
varAx[0].set_xlim(left=2007, right=2100.)
varAx[1].set_xlim(left=2007, right=2100.)
varAx[2].set_xlim(left=2007, right=2100.)
#ratioAx.set_xlim(left=0, right=100.)
#set a reasonable fontsize
plt.rcParams.update({'font.size': 16})

# Special plotting for humboldt ensemble:

for line, lineName in zip(plotLines, plotLineNames):
    if 'shelfMelt10myr' in lineName:
        line.set_linewidth(1)
    elif 'shelfMelt20myr' in lineName:
        line.set_linewidth(2)
    elif 'shelfMelt30myr' in lineName:
        line.set_linewidth(3)
    if 'm5_' in lineName or 'm10_' in lineName:
        lowCalving = 'VM180'
        medCalving = 'VM170'
        highCalving = 'VM160'
    elif 'm7_' in lineName or 'm25_' in lineName:
        if 'HadGEM2' in lineName or 'CNRM' in lineName:
            lowCalving = 'VM190'
            medCalving = 'VM180'
            highCalving = 'VM170'
        elif 'MIROC5' in lineName:
            lowCalving = 'VM180'
            medCalving = 'VM170'
            highCalving = 'VM160'
    elif 'm1_' in lineName:
        lowCalving = 'VM180'
        medCalving = 'VM170'
        highCalving = 'VM150'
    if 'smb_only' in lineName:
        line.set_color('tab:pink')
    elif '2017calvingFront' in lineName or 'calvingVelocityData' in lineName:
        line.set_color('tab:grey')
    elif highCalving in lineName:
        line.set_color('tab:purple')
    elif medCalving in lineName:
        line.set_color('tab:blue')
    elif lowCalving in lineName:
        line.set_color('tab:cyan')
    if 'm1_' in lineName:
        line.set_linestyle('dotted')
    if 'm3_' in lineName:
        line.set_linestyle('none')
    if 'm5_' in lineName or 'm10_' in lineName:
        line.set_linestyle('solid')
    elif 'm7_' in lineName or 'm25_' in lineName:
        line.set_linestyle('dashed')
    if 'noFaceMelt' in lineName:
        line.set_color('tab:pink')
    if 'm10_' in lineName or 'm25_' in lineName:
        line.set_color('tab:pink')
    if '5kmyrSpeedLimit' in lineName:
        line.set_alpha(0.5)

for bound, boundName in zip(plotBounds, plotBoundNames):
    if 'm3_' in boundName:
        bound.set_color('tab:grey')
        bound.set_edgecolor('none')
        bound.set_alpha(0.4)
        bound.set_zorder(0)
    elif 'm1_' in boundName:
        bound.set_edgecolor('tab:grey')
        bound.set_facecolor('tab:grey')
        #bound.set_hatch('xxxxxx')
        bound.set_alpha(0.6)
        
varFig.set_size_inches(15, 7 * len(variableName))
varFig.subplots_adjust(wspace=0.15, hspace=0.1)

plt.show()
#varFig.savefig('/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_melt_calv_ensemble/Humboldt_only_runs/projections/SL_contributions.png', dpi=300)
