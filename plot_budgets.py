#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 19:15:45 2020
This creates a mass budget (melting, calving, SMB) from model output.nc file. 
Output fields must include thickness, calvingThickness, faceMeltingThickness
sfcMassBalApplied, basalMassBal, bedTopography, and cellMask

@author: trevorhillebrand
"""
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from optparse import OptionParser
import time
from os.path import exists

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
fig.set_size_inches(8, 6)
insetAxs = []
budgFig, budgAx = plt.subplots(1,1)
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
    globalBudget = {}
    if 'massBudgets' not in filename:
        f = Dataset(filename, 'r+')
        f.set_auto_mask(False)
        print('Calculating mass budgets from {}'.format(filename))

        # Get mesh variables
        nCells = f.dimensions["nCells"].size
        xCell = f.variables['xCell'][:]
        yCell = f.variables['yCell'][:]
        nEdgesOnCell = f.variables['nEdgesOnCell'][:]
        cellsOnCell = f.variables["cellsOnCell"][:]
        xCell = f.variables["xCell"][:]
        areaCell = f.variables["areaCell"][:]

        # Get time variables
        deltat = np.gradient(f.variables["daysSinceStart"][:]) * s_per_day
        deltatArray = np.tile(deltat, (nCells, 1)).transpose()
        daysSinceStart = f.variables["daysSinceStart"][:]
        yr = daysSinceStart / 365. + startYear

        # Get geometry and budget variables
        thk = f.variables["thickness"][:]
        bed = f.variables["bedTopography"][:]
        sfcMassBal = f.variables["sfcMassBalApplied"][:]
        basalMassBal = f.variables["basalMassBal"][:]
        faceMeltingThickness = f.variables["faceMeltingThickness"][:] #m
        calvingThickness = f.variables["calvingThickness"][:]

        # Get masks
        cellMask = f.variables["cellMask"][:]
        iceMask = (cellMask & iceValue) // iceValue
        floatMaskKeep = (cellMask & floatValue) // floatValue  # do not update this mask
        floatMask = floatMaskKeep.copy()  # copy of floating ice mask to be updated below
        dynamicMask = (cellMask & dynamicValue) // dynamicValue
        nonDynamicMask = np.logical_not(dynamicMask) * iceMask
        groundedMaskKeep = np.logical_not(floatMask) * iceMask  # do not update this mask
        groundedMask = groundedMaskKeep.copy()  # copy of grounded ice mask to be updated below

        # add non-dynamic cells fringing grounded ice to groundedMask. These are
        # defined as non-dynamic cells with at least one grounded neighbor.
        print('Adding non-dynamic cells neighboring grounded ice to groundedMask')
        tic = time.time()
        for iTime in np.arange(0, len(deltat)):
            for iCell in np.arange(0, nCells):
              if (nonDynamicMask[iTime, iCell] == 1) and (floatMaskKeep[iTime, iCell] == 1):
                neighbors = cellsOnCell[iCell] - 1
                neighbors = neighbors[neighbors > -1]  # cellsOnCell = 0 in netCDF do not exist
                if (np.sum(groundedMaskKeep[iTime, neighbors]) >= 1):
                    groundedMask[iTime, iCell] = 1  # add this cell to groundedMask
                    floatMask[iTime, iCell] = 0  # remove it from floatMask
        toc = time.time()
        print('Finished adding non-dynamic cells to groundedMask in {:0.2f} s'.format(toc - tic))

        # Add subglacial lakes to the grounded mask. Do this by starting with
        # a seedMask of floating cells that have an ocean open or non-dynamic
        # floating neighbor and a growMask that is all floating ice.
        print('Adding subglacial lake cells to groundedMask')
        tic = time.time()
        floodFillMask = cellMask * 0
        for iTime in np.arange(0, len(deltat)):
            if (iTime % 50 == 0):
                print("Flood filling time level", iTime)
            seedMask = thk[0, :] * 0
            for iCell in np.arange(0, nCells):
                if floatMask[iTime, iCell] == 1:
                    neighbors = cellsOnCell[iCell] - 1
                    neighbors = neighbors[neighbors > -1]  # cellsOnCell = 0 in netCDF do not exist

                    # Seed the mask with floating cells that have open ocean neighbors.
                    for n in neighbors:
                        if thk[iTime,n] == 0.0 and bed[iTime,n] < 0.0:
                           seedMask[iCell] = 1
                           break

            floodFillMask[iTime, :] = flood_fill(
                                          seedMask=seedMask,
                                          growMask=floatMask[iTime, :] )  # Need to use the updated floating mask for growing to avoid growing back into the floating fringe we removed in the previous step

        floatMask = floodFillMask
        SGLmask = np.logical_and(np.logical_and(np.logical_not(floatMask),
                                                np.logical_not(groundedMask)),
                                                iceMask)
        groundedMask = np.logical_and(np.logical_not(floodFillMask), iceMask)
        toc = time.time()
        print('Finished adding {} subglacial lake cells to grounded mask in {} s'.format(np.sum(SGLmask), toc - tic))

        # Check for any shared mask elements
        print('Found {} cells shared between the floating and grounded masks'.format(np.sum(np.logical_and(floatMask,groundedMask))))

        # These masks will be missing some cells because calving might remove a
        # cell in the middle of a timestep, for instance. The best we can do is
        # add this to the nearest mask.
        floatMaskKeep2 = floatMask.copy()
        groundedMaskKeep2 = groundedMask.copy()
        tic = time.time()
        strandedCellCount = 0  # count how many cells are accounted for in this loop
        for iTime in np.arange(0, len(deltat)):
            for iCell in np.arange(0, nCells):
                if (groundedMaskKeep2[iTime, iCell] + floatMaskKeep2[iTime, iCell]) == 0 \
                    and ( (np.abs(calvingThickness[iTime, iCell]) > np.finfo('float64').eps) or
                          (np.abs(sfcMassBal[iTime, iCell]) >  np.finfo('float64').eps) or
                          (np.abs(basalMassBal[iTime, iCell]) > np.finfo('float64').eps) ):

                        distToGroundedMask = np.min(np.sqrt(
                                (xCell[np.where(groundedMaskKeep2[iTime,:]>0)] - xCell[iCell])**2
                                + (yCell[np.where(groundedMaskKeep2[iTime,:]>0)] - yCell[iCell])**2))
                        distToFloatMask = np.min(np.sqrt(
                                (xCell[np.where(floatMaskKeep2[iTime,:]>0)] - xCell[iCell])**2
                                + (yCell[np.where(floatMaskKeep2[iTime,:]>0)] - yCell[iCell])**2))

                        if distToGroundedMask < distToFloatMask:
                            groundedMask[iTime, iCell] = 1
                        elif distToFloatMask < distToGroundedMask:
                            floatMask[iTime, iCell] = 1
                        else:
                            # This shouldn't happen except maybe in very rare circumstances
                            print('Weird, the distance to the masks is the same' +
                                  'for iTime={}, iCell={}?'.format(iTime, iCell))
                        strandedCellCount += 1
        toc = time.time()
        print('finished adding {} stranded '.format(strandedCellCount) +
              'cells to grounded and float masks in {:0.2f} s'.format(toc - tic))
        masks = [groundedMask, floatMask]

        # write out masks for debugging
        if exists(filename + '_masks'):
            masksOut = Dataset(filename + '_masks', 'r+')
        else:
            print('Masks output file does not exist.'
                  ' Writing to new file {}_masks.nc'.format(filename))
            masksOut = Dataset(filename + '_masks', 'w')
            masksOut.createDimension('nCells', nCells)
            masksOut.createDimension('Time', None)
            masksOut.createVariable('groundedMask', 'i', ('Time', 'nCells'))
            masksOut.createVariable('floatMask', 'i', ('Time', 'nCells'))
            masksOut.createVariable('SGLmask', 'i', ('Time', 'nCells'))
        masksOut.variables['groundedMask'][:] = groundedMask
        masksOut.variables['floatMask'][:] = floatMask
        masksOut.variables['SGLmask'][:] = SGLmask

        f.close()
        masksOut.close
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
        colTitle = 'CNRM-CM6'
    else:
        col = 0
        colTitle = 'No climate forcing found'

    if 'm5' in filename:
        lineStyle='solid'
    elif 'm7' in filename:
        lineStyle='dashed'
    else:
        lineStyle='solid'
    # loop over floating and grounded masks
    # Define row by mask
    for row, mask in enumerate(masks):
        ax = axs[row, col]
        inset = insetAxs[nCols*row + col]
        #calculate mass budgets if using output_all_timesteps.nc
        if 'massBudgets' not in filename:
            cellAreaArray = np.tile(areaCell, (np.shape(calvingThickness)[0],1))
            totalVol = np.sum(thk * mask * cellAreaArray, axis=1)
            calvingVolFlux = np.sum(calvingThickness * mask * cellAreaArray,axis=1) #m^3
            faceMeltVolFlux = np.sum(faceMeltingThickness * cellAreaArray,axis=1) # m^3
            sfcMassBalVolFlux = np.sum(sfcMassBal * mask * cellAreaArray, axis=1) / 910. * deltat

            # Limit BMB to melting ice that is actually there. Sometimes a numerical
            # instability on one thin marginal cell creates BMB that is many orders
            # of magnitude larger than the ice thickness in that cell. This should
            # be dealt with in the model code, but deal with it here for now.
            badBasalMassBalInd = np.where( np.abs(basalMassBal) > thk * 910. / deltatArray )
            basalMassBal[badBasalMassBalInd] = (thk * 910. / deltatArray)[badBasalMassBalInd]
            basalMassBalVolFlux = np.sum(basalMassBal * mask * cellAreaArray, axis=1) / 910. * deltat

            # Calculate grounding line flux as the residual in the mass budget because groundingLineFlux
            # in globalStats.nc is not what we want here.
            if mask is groundedMask:
                title = 'Grounded Ice'
                maskName = 'groundedMask'
                GLvolFlux = ( totalVol-totalVol[0] -
                              np.cumsum(basalMassBalVolFlux) -
                              np.cumsum(sfcMassBalVolFlux) -
                              np.cumsum(-calvingVolFlux) -
                              np.cumsum(-faceMeltVolFlux) )
            elif mask is floatMask:
                title = 'Floating Ice'
                maskName = 'floatMask'
                GLvolFlux = ( totalVol-totalVol[0] -
                              np.cumsum(basalMassBalVolFlux) -
                              np.cumsum(sfcMassBalVolFlux) -
                              np.cumsum(-calvingVolFlux) )
                faceMeltVolFlux *= 0.

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

            try:
                f = Dataset(filename + '_' + mask + '_with_initial_solve.nc.cleaned', 'r')
            except:
                f = Dataset(filename + '_' + mask + '_with_initial_solve.nc', 'r')
            f.set_auto_mask(False)
            yr = f.variables['daysSinceStart'][:] / 365. + startYear
            GLvolFlux = f.variables['GLvolFlux'][:] / 1.e12
            sfcMassBalVolFlux = f.variables['sfcMassBalVolFlux'][:] / 1.e12
            basalMassBalVolFlux = f.variables['basalMassBalVolFlux'][:] / 1.e12
            faceMeltVolFlux = f.variables['faceMeltVolFlux'][:] /1.e12
            calvingVolFlux = f.variables['calvingVolFlux'][:] / 1.e12
            totalVol = f.variables['totalVol'][:] / 1.e12

            f.close()

        # Check global mass budget
        if mask is groundedMask:
            maskName = 'groundedMask'
        elif mask is floatMask:
            maskName = 'floatMask'

        globalBudget[maskName] = np.cumsum(sfcMassBalVolFlux + basalMassBalVolFlux - faceMeltVolFlux - calvingVolFlux)
        globalBudget[maskName + 'TotalVol'] = totalVol

        #Now plot the budgets!
        for plotAx in [inset, ax]:
            if mask is groundedMask:
                GLfluxPlot, = plotAx.plot(yr, GLvolFlux, c='tab:orange', linestyle=lineStyle)
                faceMeltPlot, = plotAx.plot(yr, np.cumsum(-faceMeltVolFlux), c='tab:purple', linestyle=lineStyle)
                # Uncomment line below to plot grounding line flux calculated from
                # grounded ice mask on the floating mass budget axes, for comparison
                #axs[row+1,col].plot(yr, -GLvolFlux, c='magenta', linestyle='dotted')
                if plotAx is inset:
                    inset.set_ylim(top=0.075, bottom=-.25)
            elif mask is floatMask:
                GLfluxPlot, = plotAx.plot(yr, GLvolFlux, c='tab:orange', linestyle=lineStyle)
                # Uncomment line below to plot grounding line flux calculated from
                # floating ice mask on the grounded mass budget axes, for comparison
                #axs[row-1,col].plot(yr, -GLvolFlux, c='magenta', linestyle='dotted')
                if plotAx is inset:
                    inset.set_ylim(top=.085, bottom=-.055)
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

    # Check global budget:
    globalBudget_percent_imbalance = ( (globalBudget['groundedMask'] + globalBudget['floatMask'] -
                                        globalBudget['groundedMaskTotalVol'] - globalBudget['floatMaskTotalVol'] +
                                        globalBudget['groundedMaskTotalVol'][0] + globalBudget['floatMaskTotalVol'][0]) /
                                       (globalBudget['groundedMaskTotalVol'] + globalBudget['floatMaskTotalVol'] -
                                        globalBudget['groundedMaskTotalVol'][0] - globalBudget['floatMaskTotalVol'][0]) )

    budgAx.set_title('global mass budget')
    budgAx.set_ylabel('Fractional imbalance')
    budgAx.plot(yr[1::], globalBudget_percent_imbalance[1::])



for ax in axs[1,:]:
    ax.set_xlabel('Year')
axs[0,0].set_ylabel('Total grounded volume\nchange ($10^{12}$ m$^3$)')
axs[1,0].set_ylabel('Total floating volume\nchange ($10^{12}$ m$^3$)')
axs[1,1].legend([GLfluxPlot, faceMeltPlot, sfcMassBalPlot,  basalMassBalPlot, calvingPlot, totalVolChangePlot],
               ['GL flux', 'undercutting', 'SMB', 'BMB', 'calving', 'total'], loc='upper center', bbox_to_anchor=(0.5, -0.25),
               fancybox=True, shadow=True, ncol=3)

#fig.savefig('plot_budgets.pdf', format='pdf', bbox_inches='tight')
plt.show()
