#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 10:50:17 2021

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

MIROC5_smb = Dataset('/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_1to10km/Humboldt_only_mesh/Humboldt_only_1to10km_MIROC5-rcp85_sfcMassBal.nc', 'r')
MIROC5_runoff = Dataset('/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_1to10km/Humboldt_only_mesh/Humboldt_only_1to10km_MIROC5-rcp85_Runoff.nc', 'r')
MIROC5_TF = Dataset('/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_1to10km/Humboldt_only_mesh/Humboldt_only_1to10km_MIROC5-rcp85_oceanThermalForcing.nc', 'r')

HadGEM2_smb = Dataset('/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_1to10km/Humboldt_only_mesh/Humboldt_only_1to10km_HadGEM2-rcp85_sfcMassBal.nc', 'r')
HadGEM2_runoff = Dataset('/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_1to10km/Humboldt_only_mesh/Humboldt_only_1to10km_HadGEM2-rcp85_Runoff.nc', 'r')
HadGEM2_TF = Dataset('/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_1to10km/Humboldt_only_mesh/Humboldt_only_1to10km_HadGEM2-rcp85_oceanThermalForcing.nc', 'r')

CNRM_smb = Dataset('/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_1to10km/Humboldt_only_mesh/Humboldt_1to10km_CNRM-SSP585_sfcMassBal.nc', 'r')
CNRM_runoff = Dataset('/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_1to10km/Humboldt_only_mesh/Humboldt_1to10km_CNRM-SSP585_Runoff.nc', 'r')
CNRM_TF = Dataset('/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_1to10km/Humboldt_only_mesh/Humboldt_1to10km_CNRM-SSP585_oceanThermalForcing.nc', 'r')

CTD_dir = '/Users/trevorhillebrand/Documents/Greenland/OMG_CTD/Humboldt/'

mesh = Dataset('/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_1to10km/Humboldt_only_mesh/Humboldt_1to10km_r04_20210615.nc', 'r')
bed = mesh.variables["bedTopography"][:]
areaCell = mesh.variables["areaCell"][:]
nCells = mesh.dimensions['nCells'].size
thk = mesh.variables["thickness"][:]
yr = np.arange(0,np.shape(MIROC5_smb.variables["sfcMassBal"][:])[0]) + 1950.

def meanForcing(forcing, mask):
    #meanForcing = np.sum(areaCell * forcing * mask, axis = 1) / np.sum(areaCell)
    meanForcing = np.average(forcing*mask, axis=1, weights=areaCell)
    stdForcing = np.sqrt(np.average((forcing*mask-np.tile(meanForcing, (nCells,1)).transpose())**2, axis=1, weights=areaCell))
    return meanForcing, stdForcing

mean_smb = {'MIROC5':meanForcing(MIROC5_smb.variables["sfcMassBal"][:] / 910. * 3.154e7, 
                                 (thk > 0.) * (MIROC5_smb.variables["sfcMassBal"][:] > -0.5) / 1.e13),
       'HadGEM2':meanForcing(HadGEM2_smb.variables["sfcMassBal"][:] / 910. * 3.154e7, (
               thk > 0.) * (HadGEM2_smb.variables["sfcMassBal"][:] > -0.5) / 1.e13),
       'CNRM-CM6':meanForcing(CNRM_smb.variables["sfcMassBal"][:] / 910. * 3.154e7, 
                          (thk > 0.) * (CNRM_smb.variables["sfcMassBal"][:] > -0.5) / 1.e13)}

total_smb = {'MIROC5':np.sum(MIROC5_smb.variables["sfcMassBal"][:] * areaCell * 3.154e7 *
                             (thk > 0.) * (MIROC5_smb.variables["sfcMassBal"][:] > -0.5) / 1.e13 , axis=1),
       'HadGEM2':np.sum(HadGEM2_smb.variables["sfcMassBal"][:] * areaCell * 3.154e7 *
                        (thk > 0.) * (HadGEM2_smb.variables["sfcMassBal"][:] > -0.5) / 1.e13, axis=1),
       'CNRM-CM6':np.sum(CNRM_smb.variables["sfcMassBal"][:] * areaCell * 3.154e7 *
                     (thk > 0.) * (CNRM_smb.variables["sfcMassBal"][:] > -0.5) / 1.e13, axis=1)}

runoff = {'MIROC5':meanForcing(MIROC5_runoff.variables["ismip6Runoff"][:], thk > 0.)[0],
       'HadGEM2':meanForcing(HadGEM2_runoff.variables["ismip6Runoff"][:], thk > 0.)[0],
       'CNRM-CM6':meanForcing(CNRM_runoff.variables["ismip6Runoff"][:], thk > 0.)[0]}

mean_TF = {'MIROC5':meanForcing(MIROC5_TF.variables["ismip6_2dThermalForcing"][:], (bed < 0.) * (thk > 0.))[0],
       'HadGEM2':meanForcing(HadGEM2_TF.variables["ismip6_2dThermalForcing"][:], (bed < 0.) * (thk > 0.))[0],
       'CNRM-CM6':meanForcing(CNRM_TF.variables["ismip6_2dThermalForcing"][:], (bed < 0.) * (thk > 0.))[0]}

std_TF = {'MIROC5':meanForcing(MIROC5_TF.variables["ismip6_2dThermalForcing"][:], (bed < 0.) * (thk > 0.))[1],
       'HadGEM2':meanForcing(HadGEM2_TF.variables["ismip6_2dThermalForcing"][:], (bed < 0.) * (thk > 0.))[1],
       'CNRM-CM6':meanForcing(CNRM_TF.variables["ismip6_2dThermalForcing"][:], (bed < 0.) * (thk > 0.))[1]}

fig, ax = plt.subplots(2,2, figsize=(8,8), sharex=False)
ax = ax.ravel() #easier indexing
ind = 0
climateColors = ['tab:blue', 'tab:pink', 'tab:orange']
for forcing in [total_smb, runoff, mean_TF]:
    colorInd = 0
    for climate in ['MIROC5', 'HadGEM2', 'CNRM-CM6']:
        ax[ind].plot(yr, forcing[climate], c = climateColors[colorInd], label=climate)
        colorInd += 1
    ax[ind].grid('on')
    ax[ind].set_xlabel('Year', fontsize=12)
    ax[ind].set_xlim(left=2000, right=2100)
    ind += 1

## Plot CTD data
CTDfiles = sorted(os.listdir(CTD_dir))
# set bounds for CTD data to use
boundLatMin = 79.0
boundLatMax = 79.9
boundLonMin = -70.0
boundLonMax = -64.0

resultsList = []
catTemp = np.empty(0,)
catTF = np.empty(0,)
catDepth = np.empty(0,)
plotDepth = np.linspace(-400,0, 400)
SRef = 34.4 # reference salinity
def calc_Tfreeze(z):
    Tfreeze = 0.0901 - 0.057 * SRef + 7.61e-4 * z
    return Tfreeze

for CTDfile in CTDfiles:
   if '.nc' in CTDfile:
      ctd = Dataset(CTD_dir + CTDfile, 'r')
      ctd.set_auto_mask(False)
      lat = ctd.variables["lat"][:]
      lon = ctd.variables["lon"][:]
      if (boundLatMin <= lat) and (boundLatMax >= lat) and (boundLonMin <= lon) and (boundLonMax >= lon):
         resultsList.append(CTDfile)
         #shutil.copy(dataPath+CTDfile, destPath)
         temperature = ctd.variables["temperature"][0,:]
         time = ctd.variables["time"][:]
         #profile_time = ctd.variables["profile_time"][0,:]
         depth = -1.0 * ctd.variables["depth"][0,:]
         keepIndices = np.where(temperature>-25.) # trim out bad temps
         temperature = temperature[keepIndices]
         depth = depth[keepIndices]
         Tfreeze = calc_Tfreeze(depth) #depth-dependent freezing temperature
         TF = temperature - Tfreeze #thermal forcing
         catTF = np.concatenate((catTF,TF), axis=0)
         catTemp = np.concatenate((catTemp, temperature), axis=0) #concatenate for polynomial fit
         catDepth = np.concatenate((catDepth, depth), axis=0)
         ax[3].scatter(temperature, depth, c='tab:blue', s=1)

      ctd.close()

pz = np.polyfit(catDepth, catTemp, 5)
ax[3].plot(np.polyval(pz, plotDepth), plotDepth, linewidth=2, color='black')
ax[3].set_ylabel('Depth (m)', fontsize=12)
ax[3].set_xlabel('Temperature (°C)', fontsize=12)
ax[3].set_ylim(top=25, bottom=-400)
ax[3].set_xlim(left=-1.75, right=1.)

ax[3].grid('on')

## Plot std_TF as bands. This doesn't work well yet.
#colorInd = 0
#for climate in ['MIROC5', 'HadGEM2', 'CNRM-CM6']:
#   ax[-1].fill_between(yr, mean_TF[climate] + std_TF[climate], mean_TF[climate] 
#                   - std_TF[climate], color = climateColors[colorInd], alpha=0)
#   colorInd += 1
   
ax[0].legend(fontsize='medium')

fig.subplots_adjust(wspace=0.45)
ax[0].set_ylabel('Total surface mass\nbalance (10$^{13}$ kg yr$^{-1}$)', fontsize=12)
ax[1].set_ylabel('Mean subglacial\ndischarge (kg m$^{-2}$ s$^{-1}$)', fontsize=12)
ax[2].set_ylabel('Mean ocean\nthermal forcing (°C)', fontsize=12)
