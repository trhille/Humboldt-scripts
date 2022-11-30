#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 14:48:12 2021
This script is for tuning the gamma0 parameter for the ISMIP6 ice shelf melt parameterization used on Humboldt Glacier
@author: trevorhillebrand
"""
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt


forcingFile = '/global/cscratch1/sd/trhille/Thwaites_1to8km_r02_20220427/forcings/ocean_thermal_forcing/obs/processed_obs_TF_1995-2017_8km_x_60m.nc'
meshFile = '/global/cfs/cdirs/piscees/MALI_input_files/Thwaites_1to8km_r02/Thwaites_1to8km_r02_20220427.nc'

desiredMeanMeltRate = 27. #m/yr
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
timeLevStart = 0
timeLevEnd = 1
#layerTFmean = np.mean(layerTF[timeLevStart:timeLevEnd, :, :], axis=0)
layerTFmean = np.mean(layerTF[:, :, :], axis=0)
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
floatMask = np.logical_and( (thk <= hFloat), (thk > 0.0) )[0,:]

lowerSurface = bed.copy()
lowerSurface[0, floatMask] = -rhoi / rhosw * thk[0, floatMask]

TFdraft = np.zeros(shape=(1, nCells))

for iCell in np.arange(0, nCells):
    TFdraft[0, iCell] = np.interp(lowerSurface[0, iCell], np.flip(zOcean), np.flip(layerTFmean[iCell,:])) #use TF at base of floating ice
    #TFdraft[0, iCell] = layerTFmean[iCell,2] #use TF at base of thermocline
    
# Mean thermal forcing at base of floating ice
meanTF = np.sum(TFdraft[0,floatMask] * areaCell[floatMask]) / np.sum(areaCell[floatMask])

# solve this for gamma0
# floatingBasalMassBal = ( -1.0 * gamma0 * cste / scyr * rhoi *
#                        TFdraft * np.abs(meanTF) )

gamma0 = ( np.sum(desiredMeltRate * areaCell[floatMask]) / 
         np.sum(TFdraft[0,floatMask] * areaCell[floatMask]) ) / (cste * meanTF)

print("For desired melt rate of {} m/yr at base of floating ice, mean gamma0 = {}".
      format(desiredMeltRate,gamma0))
