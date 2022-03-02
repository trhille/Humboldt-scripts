#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 09:43:28 2022

@author: trevorhillebrand
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import time
import pickle


f = Dataset('/lustre/scratch4/turquoise/trhille/debug_groundingLineFlux/tmp.nc', 'r+')
g = Dataset('/lustre/scratch4/turquoise/trhille/debug_groundingLineFlux/tmpmasks.nc', 'w')
#f = Dataset('output_all_timesteps.nc', 'r+')

f.set_auto_mask(False)

fluxAcrossGroundingLine = f.variables["fluxAcrossGroundingLine"][:]
thk = f.variables['thickness'][:]
bed = f.variables['bedTopography'][:]
sfcMassBal = f.variables["sfcMassBalApplied"][:]
basalMassBal = f.variables["basalMassBal"][:]
faceMeltingThickness = f.variables["faceMeltingThickness"][:] #m
calvingThickness = f.variables["calvingThickness"][:]
areaCell = f.variables['areaCell'][:]
edgeMask = f.variables['edgeMask'][:]
cellMask = f.variables['cellMask'][:]
cellsOnCell = f.variables['cellsOnCell'][:]
cOnE = f.variables['cellsOnEdge'][:]
nEdgesOnCell = f.variables['nEdgesOnCell'][:]
edgesOnCell = f.variables['edgesOnCell'][:]
dvEdge = f.variables["dvEdge"][:]
deltat = f.variables['deltat'][:]
xEdge = f.variables['xEdge'][:]
yEdge = f.variables['yEdge'][:]
xCell = f.variables['xCell'][:]
yCell = f.variables['yCell'][:]
nCells = f.dimensions['nCells'].size
nEdges = f.dimensions['nEdges'].size
GLyr = f.variables['daysSinceStart'][:] / 365.
cellAreaArray = np.tile(areaCell, (np.shape(thk)[0],1))

#f.close()
# Find grounding line cells that are not dynamic margin. This is the true
# grounding line (i.e., ice goes afloat as it passes this edge).
# ***Alternative would be grounding line cells that have floating neighbors.
# Compare these calculations***
dynamicMarginValue = 16
GLvalue = 256
dynamicIceValue = 2
floatingValue = 4
iceValue = 32

def flood_fill(seedMask, growMask):
    """
    Flood fill algorithm stolen from subroutine in mpas_li_calving.F, used for
    adding subglacial lakes to cellMask_grounded (and maybe other things).
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
        #print("Finished flood fill loop ", localLoopCount, " found new cells: ", newMaskCountLocal)

    return seedMask



GLmask = (edgeMask & GLvalue) // GLvalue
dynamicMarginMask = (edgeMask & dynamicMarginValue) // dynamicMarginValue
floatMask = (edgeMask & floatingValue) // floatingValue
dynamicIceMask = (edgeMask & dynamicIceValue) // dynamicIceValue

GLfluxMask = GLmask * (1 - dynamicMarginMask)

cellMask_ice = (cellMask & iceValue) // iceValue
cellMask_floatingKeep = (cellMask & floatingValue) // floatingValue
cellMask_floating = cellMask_floatingKeep.copy()
cellMask_groundedKeep = np.logical_not(cellMask_floating) * cellMask_ice
cellMask_grounded = cellMask_groundedKeep.copy()
cellMask_dynamic = (cellMask & dynamicIceValue) // dynamicIceValue
cellMask_nondynamic = np.logical_not(cellMask_dynamic) * cellMask_ice

doCalc=True
if doCalc:
   print('Adding non-dynamic cells neighboring grounded ice to groundedMask')
   tic = time.time()
   for iTime in np.arange(0, len(deltat)):
       for iCell in np.arange(0, nCells):
         if (cellMask_nondynamic[iTime, iCell] == 1) and (cellMask_floatingKeep[iTime, iCell] == 1):
           neighbors = cellsOnCell[iCell] - 1
           neighbors = neighbors[neighbors > -1]  # cellsOnCell = 0 in netCDF do not exist
           if (np.sum(cellMask_groundedKeep[iTime, neighbors]) >= 1):
               cellMask_grounded[iTime, iCell] = 1  # add this cell to groundedMask
               cellMask_floating[iTime, iCell] = 0  # remove it from floatMask
   toc = time.time()
   print('Finished adding non-dynamic cells to groundedMask in {:0.2f} s'.format(toc - tic))

   # Add subglacial lakes to the grounded mask. Do this by starting with
   # a seedMask of floating cells that have an ocean open or non-dynamic
   # floating neighbor and a growMask that is all floating ice.
   print('Adding subglacial lake cells to cellMask_grounded')
   tic = time.time()
   floodFillMask = cellMask * 0
   for iTime in np.arange(0, len(deltat)):
       print("Flood filling time level", iTime)
       seedMask = thk[0, :] * 0
       for iCell in np.arange(0, nCells):
           if cellMask_floating[iTime, iCell] == 1:
               neighbors = cellsOnCell[iCell] - 1
               neighbors = neighbors[neighbors > -1]  # cellsOnCell = 0 in netCDF do not exist
               #if ( (np.sum(cellMask_nondynamic[iTime, neighbors]) >= 1) or (np.min(thk[iTime, neighbors]) <= 1) ):

               # Seed the mask with floating cells that have open ocean neighbors.
               for n in neighbors:
                   if thk[iTime,n] == 0.0 and bed[iTime,n] < 0.0:
                      seedMask[iCell] = 1
                      break

       floodFillMask[iTime, :] = flood_fill(
                                      seedMask=seedMask,
                                      growMask=cellMask_floating[iTime, :] )  # Need to use the updated floating mask for growing to avoid growing back into the floating fringe we removed in the previous step
   # Subglacial lakes are cells in the iceMask that are not grounded and
   # are not connected to the ice shelf, as determined by flood_fill().
   #subglacialLakeMask = cellMask_ice * np.logical_not(cellMask_grounded) * np.logical_not(floodFillMask)
   #cellMask_grounded += subglacialLakeMask  # add them to grounded ice mask
   #cellMask_floating -= subglacialLakeMask  # and remove them from floating ice mask
   cellMask_floating = floodFillMask
   cellMask_grounded = np.logical_and(np.logical_not(floodFillMask), cellMask_ice)
   toc = time.time()
   #print('Finished adding {} subglacial lake cells to grounded mask in {} s'.format(np.sum(subglacialLakeMask), toc - tic))

   # These masks will be missing some cells because calving might remove a
   # cell in the middle of a timestep, for instance. The best we can do is
   # add this to the nearest mask.
   cellMask_floatingKeep2 = cellMask_floating.copy()
   cellMask_groundedKeep2 = cellMask_grounded.copy()
   tic = time.time()
   strandedCellCount = 0  # count how many cells are accounted for in this loop
   for iTime in np.arange(0, len(deltat)):
       for iCell in np.arange(0, nCells):
           if (cellMask_groundedKeep2[iTime, iCell] + cellMask_floatingKeep2[iTime, iCell]) == 0 \
               and ( (np.abs(calvingThickness[iTime, iCell]) > np.finfo('float64').eps) or
                     (np.abs(sfcMassBal[iTime, iCell]) >  np.finfo('float64').eps) or
                     (np.abs(basalMassBal[iTime, iCell]) > np.finfo('float64').eps) ):

                   distToGroundedMask = np.min(np.sqrt(
                           (xCell[np.where(cellMask_groundedKeep2[iTime,:]>0)] - xCell[iCell])**2
                           + (yCell[np.where(cellMask_groundedKeep2[iTime,:]>0)] - yCell[iCell])**2))
                   distToFloatMask = np.min(np.sqrt(
                           (xCell[np.where(cellMask_floatingKeep2[iTime,:]>0)] - xCell[iCell])**2
                           + (yCell[np.where(cellMask_floatingKeep2[iTime,:]>0)] - yCell[iCell])**2))

                   if distToGroundedMask < distToFloatMask:
                       cellMask_grounded[iTime, iCell] = 1
                   elif distToFloatMask < distToGroundedMask:
                       cellMask_floating[iTime, iCell] = 1
                   else:
                       # This shouldn't happen except maybe in very rare circumstances
                       print('Weird, the distance to the masks is the same' +
                             'for iTime={}, iCell={}?'.format(iTime, iCell))
                   strandedCellCount += 1
   toc = time.time()
   print('finished adding {} stranded '.format(strandedCellCount) +
         'cells to grounded and float masks in {:0.2f} s'.format(toc - tic))

   tic = time.time()
   print('Writing masks to .nc file.')
   # add some variables to the output file so we can write out masks for debugging/viz
   if not 'nCells' in g.dimensions:
      g.createDimension('nCells', nCells)
   if not 'Time' in g.dimensions:
      g.createDimension('Time', len(deltat))
   if not 'cellMask_grounded' in g.variables:
      g.createVariable('cellMask_grounded', 'd', ('Time','nCells'))
   if not 'cellMask_floating' in g.variables:
      g.createVariable('cellMask_floating', 'd', ('Time','nCells'))
   # Save the final masks
   g.variables['cellMask_grounded'][:]=cellMask_grounded[:]
   g.variables['cellMask_floating'][:]=cellMask_floating[:]
   toc = time.time()
   print('Finished writing masks to .nc file in {:0.2f} s'.format(toc-tic))
else:
   # load the masks from the file where we should have written them previously.
   cellMask_grounded[:] = f.variables['cellMask_grounded'][:]
   cellMask_floating[:] = f.variables['cellMask_floating'][:]

tic = time.time()
print('Calculating grToFlt')
grToFlt = np.zeros(thk.shape)
for t in range(1,len(deltat)):
    grToFlt[t-1,:]  = np.logical_and(cellMask_grounded[t-1], cellMask_floating[t]) * thk[t,:]
    grToFlt[t-1,:] -= np.logical_and(cellMask_grounded[t], cellMask_floating[t-1]) * thk[t,:]
toc = time.time()
print('finished calcuting grToFLt in {:0.2f} s'.format(toc-tic))
if not 'gTOf' in g.variables:
    g.createVariable('gTOf', 'd', ('Time','nCells'))
g.variables['gTOf'][:]=grToFlt

g2f = (grToFlt * cellAreaArray).sum(axis=1) # m3



print('Calculating GLF from scratch')
tic = time.time()
flux_array = np.zeros((len(deltat), nEdges))
for t in range(1,len(deltat)):
    maskTime = t-1
    lnv =  f.variables['layerNormalVelocity'][t,:,:]
    lt  =  f.variables['layerThicknessEdge'][t,:,:]
    for e in range(nEdges):
        c1 = cOnE[e,0]-1
        c2 = cOnE[e,1]-1
        if (((cellMask_grounded[maskTime,c1]==1) and (cellMask_floating[maskTime,c2]==1)) or
            ((cellMask_grounded[maskTime,c2]==1) and (cellMask_floating[maskTime,c1]==1))):
            # find sign of flux.  positive velo goes from cell 1 to cell 2
            if cellMask_grounded[maskTime,c1]==1:
                GLFsign=1.0 # will yield positive flux for velocity from grounded to floating
            else:
                GLFsign=-1.0
            flux_array[t,e] = (lnv[e,:]*lt[e,:]).sum() * GLFsign
toc = time.time()
print('Finished GLF calc in {:0.2f} s'.format(toc - tic))


# remove GLF on masked cells
print("cleaning GLF")
#for t in range(1,len(deltat)):
#    lnv =  f.variables['layerNormalVelocity'][t,:,:]
#    lt  =  f.variables['layerThicknessEdge'][t,:,:]
#    for iCell in range(nCells):
#        if grToFlt[t-1,iCell] != 0.0:
#            allEdges = edgesOnCell[iCell, :nEdgesOnCell[iCell]]
#            prevGLedges = flux_array[t, allEdges] != 0
#            newGLedges = flux_array[t, allEdges] == 0
#            flux_array[t, prevGLedges] = 0.0
#            for m in newGLedges:
#                c1 = cOnE[m,0]-1
#                c2 = cOnE[m,1]-1
#                if (((cellMask_grounded[t,c1]==1) and (cellMask_floating[t,c2]==1)) or
#                    ((cellMask_grounded[t,c2]==1) and (cellMask_floating[t,c1]==1))):
#                    # find sign of flux.  positive velo goes from cell 1 to cell 2
#                    if cellMask_grounded[t,c1]==1:
#                        GLFsign=1.0 # will yield positive flux for velocity from grounded to floating
#                    else:
#                        GLFsign=-1.0
#                    flux_array[t,m] = (lnv[m,:]*lt[m,:]).sum() * GLFsign




if not 'nEdges' in g.dimensions:
   g.createDimension('nEdges', nEdges)
if not 'myGLF2d' in g.variables:
   g.createVariable('myGLF2d', 'd', ('Time','nEdges')) # Note this field is on edges - won't show up in Paraview

g.variables['myGLF2d'][:]=flux_array

# Calculate scalar GLF values
myGLF = np.zeros(deltat.shape)
for t in range(1,len(deltat)):
    myGLF[t] = (flux_array[t,:]*dvEdge).sum() #units of m3/s




totalGroundedVol = np.sum(thk * cellMask_grounded * cellAreaArray, axis=1)
totalFloatingVol = np.sum(thk * cellMask_floating * cellAreaArray, axis=1)
# convert to m^3 / yr
GLfluxOnEdges = fluxAcrossGroundingLine * GLfluxMask * dvEdge
totalGLflux = np.sum(GLfluxOnEdges, axis=1)

calvingVolFlux = np.sum(calvingThickness * cellMask_grounded * cellAreaArray,axis=1) #m^3
faceMeltVolFlux = np.sum(faceMeltingThickness * cellAreaArray,axis=1) # m^3
sfcMassBalVolFlux = np.sum(sfcMassBal * cellMask_grounded * cellAreaArray, axis=1) / 910. * deltat
basalMassBalVolFlux = np.sum(basalMassBal * cellMask_grounded * cellAreaArray, axis=1) / 910. * deltat
GLflux_as_residual = ( np.cumsum(np.gradient(totalGroundedVol)) -
                       np.cumsum(basalMassBalVolFlux) -
                       np.cumsum(sfcMassBalVolFlux) -
                       np.cumsum(-calvingVolFlux) -
                       np.cumsum(-faceMeltVolFlux) )

RMSE = np.sqrt(np.mean( GLflux_as_residual - (-np.cumsum(myGLF * deltat)-np.cumsum(g2f) )))

fig, axs = plt.subplots(2)
axs[0].plot(GLyr + 2007., -np.cumsum(g2f), 'g', label='G2F')
axs[0].plot(GLyr + 2007., -np.cumsum(totalGLflux * deltat), 'r', label='cumulative GL flux')
axs[0].plot(GLyr + 2007., -np.cumsum(totalGLflux * deltat)-np.cumsum(g2f), 'r--', label='cumulative GL flux + G2F')
axs[0].plot(GLyr + 2007., -np.cumsum(myGLF * deltat), 'b', label='myGLF')
axs[0].plot(GLyr + 2007., -np.cumsum(myGLF * deltat)-np.cumsum(g2f), 'b--', label='myGLF+G2F')
axs[0].plot(GLyr + 2007., np.cumsum(np.gradient(totalGroundedVol)), 'k', label='total grounded volume change')
basalMassBalPlot, = axs[0].plot(GLyr + 2007., np.cumsum(basalMassBalVolFlux), c='tab:cyan', label='BMB')
sfcMassBalPlot, = axs[0].plot(GLyr + 2007., np.cumsum(sfcMassBalVolFlux), c='tab:pink', label='SMB')
calvingPlot, = axs[0].plot(GLyr + 2007., np.cumsum(-calvingVolFlux), c='tab:green', label='grounded calving')
faceMeltPlot, = axs[0].plot(GLyr + 2007., np.cumsum(-faceMeltVolFlux), c='tab:purple', label='face-melt')
Glflux_as_residual_plot = axs[0].plot(GLyr + 2007., GLflux_as_residual, c='magenta', label='GL flux as residual')
plt.ylabel('cumulative volume change')
axs[0].legend(loc='best', fontsize=8)


axs[1].plot(GLyr + 2007., ( ( GLflux_as_residual - (-np.cumsum(myGLF * deltat)-np.cumsum(g2f)) )
            / (-np.cumsum(myGLF * deltat)-np.cumsum(g2f)) ), label='fractional difference') 
#axs[1].plot(GLyr + 2007., np.cumsum(totalGLflux * deltat), label='cumulative GL flux')
#axs[1].plot(GLyr + 2007., np.cumsum(myGLF * deltat), label='myGLF')
#axs[1].plot(GLyr + 2007., np.cumsum(np.gradient(totalFloatingVol)), 'k', label='total floating volume change')
plt.ylabel('fractional difference')
axs[1].legend()

print('RMSE between residual and myGLF +G2F: {:0.2f}'.format(RMSE))

#axs[1].plot(GLyr + 2007., (thk * cellMask_floating * cellAreaArray).sum(axis=1), '-b')
#plt.ylabel('floating vol')
#axs[2].plot(GLyr + 2007., (thk * cellMask_grounded * cellAreaArray).sum(axis=1), '-r')
#plt.ylabel('grounded vol')

f.close()
g.close()

plt.show()
