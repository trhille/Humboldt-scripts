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


f = Dataset('/global/cscratch1/sd/trhille/Humboldt_1to10km_r04_20210615/m5/MIROC5/VM170_shelfMelt20myr/debug_groundingLineFlux_noForcing/output_all_timesteps.nc', 'r')

f.set_auto_mask(False)

fluxAcrossGroundingLine = f.variables["fluxAcrossGroundingLine"][:]
thk = f.variables['thickness'][:]
areaCell = f.variables['areaCell'][:]
edgeMask = f.variables['edgeMask'][:]
cellMask = f.variables['cellMask'][:]
cellsOnCell = f.variables['cellsOnCell'][:]
dvEdge = f.variables["dvEdge"][:]
deltat = f.variables['deltat'][:]
nCells = f.dimensions['nCells'].size
GLyr = f.variables['daysSinceStart'][:] / 365.
cellAreaArray = np.tile(areaCell, (np.shape(thk)[0],1))

f.close()
# Find grounding line cells that are not dynamic margin. This is the true 
# grounding line (i.e., ice goes afloat as it passes this edge).
# ***Alternative would be grounding line cells that have floating neighbors. 
# Compare these calculations***
dynamicMarginValue = 16
GLvalue = 256
dynamicIceValue = 2
floatingValue = 4

GLmask = (edgeMask & GLvalue) // GLvalue
dynamicMarginMask = (edgeMask & dynamicMarginValue) // dynamicMarginValue
floatMask = (edgeMask & floatingValue) // floatingValue
dynamicIceMask = (edgeMask & dynamicIceValue) // dynamicIceValue

GLfluxMask = GLmask * (1 - dynamicMarginMask)

cellMask_floatingKeep = (cellMask & floatingValue) // floatingValue
cellMask_floating = cellMask_floatingKeep.copy()
cellMask_groundedKeep = 1 - cellMask_floating
cellMask_grounded = cellMask_groundedKeep.copy()
cellMask_dynamic = (cellMask & dynamicIceValue) // dynamicIceValue
cellMask_nondynamic = 1 - cellMask_dynamic

print('Adding non-dynamic cells neighboring grounded ice to groundedMask')
tic = time.time()
for iTime in np.arange(0, len(deltat)):
    for iCell in np.arange(0, nCells):
        neighbors = cellsOnCell[iCell] - 1
        neighbors = neighbors[neighbors > -1]  # cellsOnCell = 0 in netCDF do not exist
        if (cellMask_nondynamic[iTime, iCell] == 1) and (np.sum(cellMask_groundedKeep[iTime, neighbors]) >= 1):
            cellMask_grounded[iTime, iCell] = 1  # add this cell to groundedMask
            cellMask_floating[iTime, iCell] = 0  # remove it from floatMask
toc = time.time()
print('Finished adding non-dynamic cells to groundedMask in {} s'.format(toc - tic))


totalGroundedVol = np.sum(thk * cellMask_grounded * cellAreaArray, axis=1)
# convert to m^3 / yr
GLfluxOnEdges = fluxAcrossGroundingLine * GLfluxMask * dvEdge 
totalGLflux = np.sum(GLfluxOnEdges, axis=1)
plt.plot(GLyr + 2007., -np.cumsum(totalGLflux * deltat), label='cumulative GL flux')
plt.plot(GLyr + 2007., np.cumsum(np.gradient(totalGroundedVol)), label='total grounded volume change')
plt.legend()
plt.show()                     
