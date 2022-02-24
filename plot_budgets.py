#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 19:15:45 2020
This creates a mass budget (melting, calving, SMB) from model output.nc file. 
Output fields must include thickness, calvingThickness, faceMeltRateApplied, sfcMassBalApplied,
groundedMarineMarginMask

@author: trevorhillebrand
"""
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from optparse import OptionParser
import time

print("** Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = OptionParser(description=__doc__)
parser.add_option("-f", dest="filename", help="filename for plotting", metavar="FILENAME")
parser.add_option("-s", dest="startYear", help="Calendar year at time 0", default=2007, metavar="TIME")
options, args = parser.parse_args()

filenames = options.filename.split(',') #create list of filenames to loop over
startYear = options.startYear
nFiles = len(filenames)
rhoi = 910.
s_per_day = 86400.
nRows = 2 #nFiles * 2
nCols = 3 #nFiles

fig, axs = plt.subplots(nrows=nRows, ncols=nCols, sharex=True, sharey='row')
fig.set_size_inches(15, 12)
insetAxs = []
for ii,ax in enumerate(axs.ravel()):
    if ii <= 2: insetLoc = 3; insetWidth="40%"; insetHeight="45%"
    else: insetLoc = 2; insetWidth="40%"; insetHeight="35%"
    insetAxs.append(inset_axes(ax, insetWidth, height=insetHeight, loc=insetLoc))
#if nRows > 1 or nCols > 1:
#    axs = axs.ravel() #easier to index with flattened axs array; remove last subplot if nFiles is odd

floatValue = 4 # value in cellMask denoting floating ice
dynamicValue = 2 # value in cellMask denoting dynamic ice.
iceValue = 32
def flood_fill(seedMask, growMask):
    """
    Flood fill algorithm stolen from subroutine in mpas_li_calving.F, used for
    adding subglacial lakes to groundedMask (and maybe other things).
    """
    localLoopCount = 0
    newMaskCountLocal = 1
    while (newMaskCountLocal > 0):
        localLoopCount += 1

        newMaskCountLocal = 0
        seedMaskOld = seedMask.copy()

        for iCell in range(nCells):
            if ( growMask[iCell] > 0 ):
                # If it has a marked neighbor, then add it to the mask
                for n in range(nEdgesOnCell[iCell]):
                        jCell = cellsOnCell[iCell, n] - 1
                        if ( (jCell > -1) and (seedMaskOld[jCell] == 1) and (seedMask[iCell] != 1) ):
                           seedMask[iCell] = 1
                           newMaskCountLocal += 1
                           # skip the rest of this loop - no need to
                           # check additional neighbors
                           continue

    return seedMask


for ii, filename in enumerate(filenames):
    print(filename)
    if 'massBudgets' not in filename:
        f = Dataset(filename, 'r+')
        f.set_auto_mask(False)
        print('Calculating mass budgets from {}'.format(filename))
        g = Dataset(filename.replace('output_all_timesteps', 'globalStats'))
        deltat = np.gradient(f.variables["daysSinceStart"][:]) * s_per_day
        daysSinceStart = f.variables["daysSinceStart"][:]
        yr = daysSinceStart / 365. + startYear

        nCells = f.dimensions["nCells"].size
        nEdgesOnCell = f.variables['nEdgesOnCell'][:]
        cellsOnCell = f.variables["cellsOnCell"][:]

        GLflux = g.variables["groundingLineFlux"][:]
        GLflux = np.interp(yr, g.variables["daysSinceStart"][:] / 365., GLflux)
        g.close()
        thk = f.variables["thickness"][:]
        sfcMassBal = f.variables["sfcMassBalApplied"][:]
        basalMassBal = f.variables["basalMassBal"][:]
        faceMeltingThickness = f.variables["faceMeltingThickness"][:] #m
        calvingThickness = f.variables["calvingThickness"][:]
        xCell = f.variables["xCell"][:]
        areaCell = f.variables["areaCell"][:]
        
        cellMask = f.variables["cellMask"][:]
        iceMask = (cellMask & iceValue) // iceValue
        floatMaskKeep = (cellMask & floatValue) // floatValue
        floatMask = floatMaskKeep.copy()
        dynamicMask = (cellMask & dynamicValue) // dynamicValue
        nonDynamicMask = (1 - dynamicMask) * iceMask
        groundedMaskKeep = (thk > 1) - floatMask
        groundedMask = groundedMaskKeep.copy()

        # Add subglacial lakes to the grounded mask. Do this by starting with
        # a seedMask of floating cells that have an ocean open or non-dynamic
        # floating neighbor and a growMask that is all floating ice.
        print('Adding subglacial lake cells to groundedMask')
        tic = time.time()
        floodFillMask = cellMask * 0
        for iTime in np.arange(0, len(deltat)):
            seedMask = thk[0, :] * 0
            for iCell in np.arange(0, nCells):
                neighbors = cellsOnCell[iCell] - 1
                neighbors = neighbors[neighbors > -1]  # cellsOnCell = 0 in netCDF do not exist
                if ( (floatMask[iTime, iCell] == 1) and
                ( (np.sum(nonDynamicMask[iTime, neighbors]) >= 1)
                or (np.min(thk[iTime, neighbors]) <= 1) ) ):
                    seedMask[iCell] = 1

            floodFillMask[iTime, :] = flood_fill(
                                           seedMask=seedMask,
                                           growMask=floatMaskKeep[iTime, :] )

        subglacialLakeMask = iceMask - groundedMask - floodFillMask
        groundedMask += subglacialLakeMask
        floatMask -= subglacialLakeMask
        toc = time.time()
        print('Finished adding {} subglacial lake cells to grounded mask in {} s'.format(np.sum(subglacialLakeMask), toc - tic))

        # add non-dynamic cells fringing grounded ice to groundedMask
        print('Adding non-dynamic cells neighboring grounded ice to groundedMask')
        tic = time.time()
        for iTime in np.arange(0, len(deltat)):
            for iCell in np.arange(0, nCells):
                neighbors = cellsOnCell[iCell] - 1
                neighbors = neighbors[neighbors > -1]  # cellsOnCell = 0 in netCDF do not exist
                if (nonDynamicMask[iTime, iCell] == 1) and (np.sum(groundedMaskKeep[iTime, neighbors]) >= 1):
                    groundedMask[iTime, iCell] = 1
                    floatMask[iTime, iCell] = 0
        toc = time.time()
        print('Finished adding non-dynamic cells to groundedMask in {} s'.format(toc - tic))

#        groundedMask = np.concatenate(([groundedMask[0, :]], groundedMask[0:-1, :]))  # must use the mask from the previous timestep
#        floatMask = np.concatenate(([floatMask[0, :]], floatMask[0:-1, :]))
        masks = [groundedMask, floatMask]
        # write out masks for debugging
        f.variables['groundedMask'][:] = groundedMask
        f.variables['floatMask'][:] = floatMask
        f.variables['SGLmask'][:] = subglacialLakeMask
        f.close()
    else:
        groundedMask ='groundedMask'
        floatMask = 'floatMask'
        masks = [groundedMask, floatMask]
   
    # Define columns by climate forcing
    if 'MIROC5' in filename:
        col = 0
        colTitle = 'MIROC5'
    elif 'HadGEM2' in filename:
        col = 1
        colTitle = 'HadGEM2'
    elif 'CNRM' in filename:
        col = 2
        colTitle = 'CNRM'

    if 'm5' in filename:
        lineStyle='solid'
    elif 'm7' in filename:
        lineStyle='dashed'
    # loop over floating and grounded masks
    # Define row by mask
    for row, mask in enumerate(masks):
        ax = axs[row, col]
        inset = insetAxs[nCols*row + col]
        #calculate mass budgets if using output_all_timesteps.nc
        if 'massBudgets' not in filename:
            cellAreaArray = np.tile(areaCell, (np.shape(calvingThickness)[0],1))

            totalVol = np.sum(thk * mask * cellAreaArray, axis=1)
            calvingVolFlux = np.sum(calvingThickness * np.concatenate(([mask[0,:]], mask[0:-1, :])) * cellAreaArray,axis=1) #m^3
            faceMeltVolFlux = np.sum(faceMeltingThickness * cellAreaArray,axis=1) # m^3
            sfcMassBalVolFlux = np.sum(sfcMassBal * mask * cellAreaArray, axis=1) / 910. * deltat
            basalMassBalVolFlux = np.sum(basalMassBal * mask * cellAreaArray, axis=1) / 910. * deltat
            GLvolFlux = GLflux * deltat / 3.154e7 / 910. #m^3

            if mask is groundedMask:
                title = 'Grounded Ice'
                maskName = 'groundedMask'
                GLvolFlux = GLvolFlux * -1. #mass loss for grounded ice
#                basalMassBalVolFlux *= 0. # I don't know why, but there are huge negative spikes for grounded ice, likely related to surges?
            elif mask is floatMask:
                title = 'Floating Ice'
                maskName = 'floatMask'
                GLvolFlux = 0. * GLvolFlux #np.cumsum(totalVol - totalVol[0] - sfcMassBalVolFlux + faceMeltVolFlux + calvingVolFlux)

            #Now save these timeseries so we don't have to calculate them from the output files every time
            outfile = filename.replace('output_all_timesteps', 'massBudgets_'+maskName)
            print('Writing flux timeseries to ' + outfile)
            o = Dataset(outfile, 'w')
            o.createDimension('Time', len(daysSinceStart))
            o.createVariable('daysSinceStart', 'f', 'Time')
            o.createVariable('GLvolFlux', 'f', 'Time')
            o.createVariable('sfcMassBalVolFlux', 'f', 'Time')
            o.createVariable('basalMassBalVolFlux', 'f', 'Time')
            o.createVariable('faceMeltVolFlux', 'f', 'Time')
            o.createVariable('calvingVolFlux', 'f', 'Time')
            o.createVariable('totalVol', 'f', 'Time')

            o.variables['daysSinceStart'][:] = daysSinceStart
            o.variables['GLvolFlux'][:] = GLvolFlux
            o.variables['sfcMassBalVolFlux'][:] = sfcMassBalVolFlux
            o.variables['basalMassBalVolFlux'][:] = basalMassBalVolFlux
            o.variables['faceMeltVolFlux'][:] = faceMeltVolFlux
            o.variables['calvingVolFlux'][:] = calvingVolFlux
            o.variables['totalVol'][:] = totalVol
            o.close()

        else: #load budgets

            if mask is groundedMask:
                title = 'grounded Ice'
            elif mask is floatMask:
                title = 'floating Ice'

            f = Dataset(filename + '_' + mask + '.nc', 'r')
            f.set_auto_mask(False)
            yr = f.variables['daysSinceStart'][:] / 365. + startYear
            GLvolFlux = f.variables['GLvolFlux'][:] / 1.e12
            sfcMassBalVolFlux = f.variables['sfcMassBalVolFlux'][:] / 1.e12
            basalMassBalVolFlux = f.variables['basalMassBalVolFlux'][:] / 1.e12
            faceMeltVolFlux = f.variables['faceMeltVolFlux'][:] /1.e12
            calvingVolFlux = f.variables['calvingVolFlux'][:] / 1.e12
            totalVol = f.variables['totalVol'][:] / 1.e12
            f.close()
        #Now plot the budgets!
        #budgetSumPlot, = ax.plot(yr, np.cumsum(massBudget) - massBudget[0], c='tab:blue');
        massBudget = sfcMassBalVolFlux + basalMassBalVolFlux - faceMeltVolFlux - calvingVolFlux + GLvolFlux
        for plotAx in [inset, ax]:
            if mask is groundedMask:
                GLfluxPlot, = plotAx.plot(yr, (totalVol - totalVol[0]) - np.cumsum(basalMassBalVolFlux) - np.cumsum(sfcMassBalVolFlux) - np.cumsum(-calvingVolFlux) - np.cumsum(-faceMeltVolFlux), c='tab:orange', linestyle=lineStyle)
                #GLfluxPlot, = plotAx.plot(yr, np.cumsum(GLvolFlux), c='tab:orange', linestyle=lineStyle)
                faceMeltPlot, = plotAx.plot(yr, np.cumsum(-faceMeltVolFlux), c='tab:purple', linestyle=lineStyle)
                if plotAx is inset: inset.set_ylim(top=0.075, bottom=-.25)
            elif mask is floatMask:
                GLfluxPlot, = plotAx.plot(yr, (totalVol - totalVol[0]) - np.cumsum(sfcMassBalVolFlux) - 
                                               np.cumsum(-calvingVolFlux) - np.cumsum(basalMassBalVolFlux),
                                               c='tab:orange', linestyle=lineStyle)
                basalMassBalPlot, = plotAx.plot(yr, np.cumsum(basalMassBalVolFlux), c='tab:cyan', linestyle=lineStyle)
                if plotAx is inset: inset.set_ylim(top=.07, bottom=-.03)
            basalMassBalPlot, = plotAx.plot(yr, np.cumsum(basalMassBalVolFlux), c='tab:cyan', linestyle=lineStyle)
            sfcMassBalPlot, = plotAx.plot(yr, np.cumsum(sfcMassBalVolFlux), c='tab:pink', linestyle=lineStyle)
            calvingPlot, = plotAx.plot(yr, np.cumsum(-calvingVolFlux), c='tab:green', linestyle=lineStyle)
            totalVolChangePlot, = plotAx.plot(yr, totalVol - totalVol[0], c='black', linestyle=lineStyle); 
        ax.grid('on')
        if row==0:
            ax.set_title(colTitle)
            inset.tick_params(labelbottom=False, labeltop=True, labelleft=False, labelright=True)
            inset.tick_params(axis='both', which='major', pad=-3)
        elif row==1:
            inset.set_yticks([0, 0.05])
            inset.tick_params(labelbottom=True, labeltop=False, labelright=True, labelleft=False)
            inset.tick_params(axis='both', which='major', pad=-3)
        inset.grid('on')
        inset.set_xlim(left=2007., right=2021.)
        inset.set_xticks([2007, 2017])
        inset.set_xticklabels(["'07", "'17"], ha='left', rotation=0)

for ax in axs[1,:]:
    ax.set_xlabel('Year')
    
axs[0,0].set_ylabel('Total grounded volume\nchange ($10^{12}$ m$^3$)')
axs[1,0].set_ylabel('Total floating volume\nchange ($10^{12}$ m$^3$)')
axs[1,1].legend([GLfluxPlot, faceMeltPlot, sfcMassBalPlot,  basalMassBalPlot, calvingPlot, totalVolChangePlot],
               ['GL flux', 'undercutting', 'SMB', 'BMB', 'calving', 'total'], loc='upper center', bbox_to_anchor=(0.5, -0.25),
               fancybox=True, shadow=True, ncol=3)

#remove undercutting from floating plots, BMB from grounded plots
#fig.savefig('plot_budgets', dpi=400, bbox_inches='tight')
plt.show()
