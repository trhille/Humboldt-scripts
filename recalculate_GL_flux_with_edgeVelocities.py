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


f = Dataset('/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_' + 
            '1to10km/r04_Results/debugMassBudgets/flux.nc.cleaned',
            'r')

f.set_auto_mask(False)
fluxAcrossGroundingLine = f.variables["fluxAcrossGroundingLine"][:]
edgeMask = f.variables['edgeMask'][:]
dvEdge = f.variables["dvEdge"][:]
deltat = f.variables['deltat'][:]
GLyr = f.variables['daysSinceStart'][:] / 365.

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

# convert to m^3 / yr
GLfluxOnEdges = fluxAcrossGroundingLine * GLfluxMask * dvEdge 
totalGLflux = np.sum(GLfluxOnEdges, axis=1)
plt.plot(GLyr + 2007., np.cumsum(totalGLflux * deltat))

                
    

