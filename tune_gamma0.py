#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 14:48:12 2021
This script is for tuning the gamma0 parameter for the ISMIP6 ice shelf melt parameterization used on Humboldt Glacier
@author: trevorhillebrand
"""
from netCDF4 import Dataset
import numpy as np

forcingFile = '/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_1to10km/Humboldt_only_mesh/Humboldt_only_1to10km_MIROC5-rcp85_3dOceanThermalForcing.nc'
meshFile = '/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_1to10km/Humboldt_only_mesh/Humboldt_1to10km_r04_20210615.nc'

desiredMeltRate = 30. #m/yr
#get necessary initial condition fields
mesh = Dataset(meshFile, 'r')
mesh.set_auto_mask(False)
thk = mesh.variables["thickness"][:]
bed = mesh.variables["bedTopography"][:]
nCells = mesh.dimensions["nCells"].size
xCell = mesh.variables["xCell"][:]
yCell = mesh.variables["yCell"][:]
areaCell = mesh.variables["areaCell"][:]

f = Dataset(forcingFile, 'r')
f.set_auto_mask(False)
layerTF = f.variables["ismip6shelfMelt_3dThermalForcing"][:]
zOcean = f.variables["ismip6shelfMelt_zOcean"][:]
timeLevStart = 57
timeLevEnd = 68
layerTFmean = np.mean(layerTF[timeLevStart:timeLevEnd, :, :], axis=0)
Time = f.dimensions["Time"].size

#constants
rhosw = 1028.
rhoi = 910.
cp_seawater = 3.974e3
latent_heat_ice = 335.0e3
scyr = 60. * 60. * 24. * 365.

cste = (rhosw*cp_seawater/(rhoi*latent_heat_ice))**2.

#get lowerSurface from thickness, bedTopography, and floatation thickness
hFloat = -bed * rhosw / rhoi #flotation thickness
hFloat[bed>=0.] = 0.0
floatMask = ( (thk <= hFloat) & (thk > 0.0) )[0,:]

lowerSurface = bed * 0.0
lowerSurface[0, (1 - floatMask)] = bed[0, (1 - floatMask)]
lowerSurface[0, floatMask] = -rhoi / rhosw * thk[0, floatMask]

TFdraft = np.zeros(shape=(1, nCells))

for iCell in np.arange(0, nCells):
    TFdraft[0, iCell] = np.interp(bed[0, iCell], np.flip(zOcean), np.flip(layerTFmean[iCell,:])) #use TF at base of floating ice
    #TFdraft[0, iCell] = layerTFmean[iCell,2] #use TF at base of thermocline
    
# Mean thermal forcing at base of floating ice
meanTF = np.sum(TFdraft[0,floatMask] * areaCell[floatMask]) / np.sum(areaCell[floatMask])

# solve this for gamma0
# floatingBasalMassBal = ( -1.0 * gamma0 * cste / scyr * rhoi *
#                        TFdraft * np.abs(meanTF) )

desiredBMB = -desiredMeltRate * rhoi / scyr # kg / m^2 / s

gamma0 = -scyr * desiredBMB / (cste * rhoi * TFdraft[0,floatMask] * meanTF)
gamma0mean = np.mean(gamma0)


print("For desired melt rate of {} m/yr at base of floating ice, mean gamma0 = {}".
      format(desiredMeltRate,gamma0mean))
