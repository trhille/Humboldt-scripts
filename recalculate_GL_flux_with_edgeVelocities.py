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


f = Dataset('/global/cscratch1/sd/trhille/Humboldt_1to10km_r04_20210615/m5/MIROC5/VM170_shelfMelt20myr/debug_groundingLineFlux_noForcing/output_editable.nc', 'r+')
#f = Dataset('output_all_timesteps.nc', 'r+')

f.set_auto_mask(False)

fluxAcrossGroundingLine = f.variables["fluxAcrossGroundingLine"][:]
thk = f.variables['thickness'][:]
bed = f.variables['bedTopography'][:]
areaCell = f.variables['areaCell'][:]
edgeMask = f.variables['edgeMask'][:]
cellMask = f.variables['cellMask'][:]
cellsOnCell = f.variables['cellsOnCell'][:]
cOnE = f.variables['cellsOnEdge'][:]
nEdgesOnCell = f.variables['nEdgesOnCell'][:]
dvEdge = f.variables["dvEdge"][:]
deltat = f.variables['deltat'][:]
xEdge = f.variables['xEdge'][:]
yEdge = f.variables['yEdge'][:]
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
print('Finished adding non-dynamic cells to groundedMask in {} s'.format(toc - tic))

doCalc=True
if doCalc:
   # Add subglacial lakes to the grounded mask. Do this by starting with
   # a seedMask of floating cells that have an ocean open or non-dynamic
   # floating neighbor and a growMask that is all floating ice.
   print('Adding subglacial lake cells to groundedMask')
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
else:
   # load the masks from the file where we should have written them previously.
   cellMask_grounded[:] = f.variables['cellMask_grounded'][:]
   cellMask_floating[:] = f.variables['cellMask_floating'][:]


# add some variables to the output file so we can write out masks for debugging/viz
if not 'cellMask_grounded' in f.variables:
   f.createVariable('cellMask_grounded', 'd', ('Time','nCells'))
if not 'cellMask_floating' in f.variables:
   f.createVariable('cellMask_floating', 'd', ('Time','nCells'))
# Save the final masks
f.variables['cellMask_grounded'][:]=cellMask_grounded[:]
f.variables['cellMask_floating'][:]=cellMask_floating[:]


if not 'myGLF' in f.variables:
      f.createVariable('myGLF', 'd', ('Time','nEdges')) # Note this field is on edges - won't show up in Paraview
tic = time.time()
myGLF = np.zeros(deltat.shape)
for t in range(len(deltat)):
    lnv =  f.variables['layerNormalVelocity'][t,:,:]
    lt  =  f.variables['layerThicknessEdge'][t,:,:]
    flux =  np.zeros( (nEdges,) )
    for e in range(nEdges):
        c1 = cOnE[e,0]-1
        c2 = cOnE[e,1]-1
        if (((cellMask_grounded[t,c1]==1) and (cellMask_floating[t,c2]==1)) or
            ((cellMask_grounded[t,c2]==1) and (cellMask_floating[t,c1]==1))):
            # find sign of flux.  positive velo goes from cell 1 to cell 2
            if cellMask_grounded[t,c1]==1:
                GLFsign=1.0 # will yield positive flux for velocity from grounded to floating
            else:
                GLFsign=-1.0
            flux[e] = (lnv[e,:]*lt[e,:]).sum() * GLFsign
    myGLF[t] = (flux*dvEdge).sum() #units of m3/s
    f.variables['myGLF'][t,:]=flux
toc = time.time()
print('Finished GLF calc in {} s'.format(toc - tic))


grToFlt = np.zeros(thk.shape)
for t in range(1,len(deltat)):
    grToFlt[t-1,:]  = np.logical_and(cellMask_grounded[t-1], cellMask_floating[t]) * thk[t,:]
    grToFlt[t-1,:] -= np.logical_and(cellMask_grounded[t], cellMask_floating[t-1]) * thk[t,:]
if not 'gTOf' in f.variables:
    f.createVariable('gTOf', 'd', ('Time','nCells'))
f.variables['gTOf'][:]=grToFlt


g2f = (grToFlt * cellAreaArray).sum(axis=1) # m3



totalGroundedVol = np.sum(thk * cellMask_grounded * cellAreaArray, axis=1)
totalFloatingVol = np.sum(thk * cellMask_floating * cellAreaArray, axis=1)
# convert to m^3 / yr
GLfluxOnEdges = fluxAcrossGroundingLine * GLfluxMask * dvEdge
totalGLflux = np.sum(GLfluxOnEdges, axis=1)

fig, axs = plt.subplots(2)
axs[0].plot(GLyr + 2007., -np.cumsum(g2f), 'g', label='G2F')
axs[0].plot(GLyr + 2007., -np.cumsum(totalGLflux * deltat), 'r', label='cumulative GL flux')
axs[0].plot(GLyr + 2007., -np.cumsum(totalGLflux * deltat)-np.cumsum(g2f), 'r--', label='cumulative GL flux + G2F')
axs[0].plot(GLyr + 2007., -np.cumsum(myGLF * deltat), 'b', label='myGLF')
axs[0].plot(GLyr + 2007., -np.cumsum(myGLF * deltat)-np.cumsum(g2f), 'b--', label='myGLF+G2F')
axs[0].plot(GLyr + 2007., np.cumsum(np.gradient(totalGroundedVol)), 'k', label='total grounded volume change')
plt.ylabel('cumulative volume change')
axs[0].legend()

axs[1].plot(GLyr + 2007., np.cumsum(totalGLflux * deltat), label='cumulative GL flux')
axs[1].plot(GLyr + 2007., np.cumsum(myGLF * deltat), label='myGLF')
axs[1].plot(GLyr + 2007., np.cumsum(np.gradient(totalFloatingVol)), 'k', label='total floating volume change')
plt.ylabel('cumulative volume change')
axs[1].legend()


#axs[1].plot(GLyr + 2007., (thk * cellMask_floating * cellAreaArray).sum(axis=1), '-b')
#plt.ylabel('floating vol')
#axs[2].plot(GLyr + 2007., (thk * cellMask_grounded * cellAreaArray).sum(axis=1), '-r')
#plt.ylabel('grounded vol')


plt.show()
