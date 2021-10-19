#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 25 14:48:53 2021

@author: trevorhillebrand
"""
import os
import sys
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import argparse

#dataPath = '/Users/trevorhillebrand/Documents/Greenland/OMG_CTD/Humboldt/'
#forcingFile = '/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_1to10km/Humboldt_1to10km_cull5km/tmp.nc'
#meshFile = '/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_1to10km/Humboldt_1to10km_cull5km/Humboldt_1to10km_r04_20210617.nc'

print("== Gathering information.  (Invoke with --help for more details. All arguments are optional)\n")
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.description = __doc__
parser.add_argument("-d", "--dataPath", dest="dataPath", help="Path to directory containing CTD data", default="./", metavar="DIRNAME")
parser.add_argument("-f", "--forcingFile", dest="forcingFile", help="Path to file that contains 2D thermal forcing and will have 3D thermal forcing written to it.", metavar="FILENAME")
parser.add_argument("-m", "--meshFile", dest="meshFile", help="Path to file containing base mesh. Must contain bedTopography.", metavar="FILENAME")

args = parser.parse_args()

dataPath = args.dataPath
forcingFile = args.forcingFile
meshFile = args.meshFile


mesh = Dataset(meshFile, 'r')
f = Dataset(forcingFile, 'r+')
nCells = f.dimensions["nCells"].size
nTime = f.dimensions["Time"].size
nISMIP6OceanLayers = 4
oceanTF2D = f.variables["ismip6_2dThermalForcing"][:]
bed = mesh.variables["bedTopography"][0,:]

CTDfiles = sorted(os.listdir(dataPath))

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

fig, ax = plt.subplots(1,1)
ax.grid()

SRef = 34.4 # reference salinity
def calc_Tfreeze(z):
    Tfreeze = 0.0901 - 0.057 * SRef + 7.61e-4 * z
    return Tfreeze

for CTDfile in CTDfiles:
   if '.nc' in CTDfile:
      ctd = Dataset(dataPath + CTDfile, 'r')
      ctd.set_auto_mask(False)
      lat = ctd.variables["lat"][:]
      lon = ctd.variables["lon"][:]
      if (boundLatMin <= lat) and (boundLatMax >= lat) and (boundLonMin <= lon) and (boundLonMax >= lon):
         resultsList.append(CTDfile)
         #shutil.copy(dataPath+CTDfile, destPath)
         temperature = ctd.variables["temperature"][0,:]
         depth = -1.0 * ctd.variables["depth"][0,:]
         keepIndices = np.where(temperature>-25.) # trim out bad temps
         temperature = temperature[keepIndices]
         depth = depth[keepIndices]
         Tfreeze = calc_Tfreeze(depth) #depth-dependent freezing temperature
         TF = temperature - Tfreeze #thermal forcing
         catTF = np.concatenate((catTF,TF), axis=0)
         catTemp = np.concatenate((catTemp, temperature), axis=0) #concatenate for polynomial fit
         catDepth = np.concatenate((catDepth, depth), axis=0)
         ax.scatter(temperature, depth)

      ctd.close()

pz = np.polyfit(catDepth, catTemp, 5)
ax.plot(np.polyval(pz, plotDepth), plotDepth, linewidth=2, color='black')
ax.set_ylabel('Depth (m)')
ax.set_xlabel('Temperature (deg. C)')

# Find top and bottom of thermocline
thermoclineSlope = np.gradient(np.polyval(pz, plotDepth)) / np.gradient(plotDepth)
zThermoclineTop = plotDepth[np.argmin(np.polyval(pz, plotDepth))]
zThermoclineBot = -250.
ax.scatter(np.polyval(pz, [zThermoclineTop, zThermoclineBot]), 
           [zThermoclineTop, zThermoclineBot], color='black', marker='*')

layerDepths = np.array([0., zThermoclineTop, zThermoclineBot, -1.e3])
TfreezeLayerDepths = calc_Tfreeze(layerDepths)

layerTF = np.ndarray((nTime, nCells, nISMIP6OceanLayers))

for Time in np.arange(0, nTime):
    for iCell in np.arange(0, nCells):
        layerTemp = np.polyval(pz, [zThermoclineTop, zThermoclineTop, 
                             zThermoclineBot, zThermoclineBot])
        layerTF[Time, iCell, :] = layerTemp - TfreezeLayerDepths
        # Now adjust layerTF to match ismip6_2dThermalForcing at bottom of thermocline
        layerTF[Time, iCell, :] = ( layerTF[Time, iCell, :] + oceanTF2D[Time, iCell] - 
               (np.polyval(pz, zThermoclineBot) - TfreezeLayerDepths[2]) )

#get rid of negative thermal forcing
layerTF[layerTF<0.0] = 0.0
# Create necessary dimensions and variables
try:
    f.createDimension('nISMIP6OceanLayers', size=4)
except: 
    print('Dimension nISMIP6OceanLayers already exists')
try: 
    f.createVariable('ismip6shelfMelt_3dThermalForcing', datatype='d', dimensions=(
        'Time', 'nCells', 'nISMIP6OceanLayers'))
except: 
    print('Variable ismip6shelfMelt_3dThermalForcing already exists')
try:
    f.createVariable('ismip6shelfMelt_zOcean', datatype='d', dimensions=('nISMIP6OceanLayers'))
except:
    print('Variable ismip6shelfMelt_zOcean already exists')

# Populate variables
f.variables["ismip6shelfMelt_3dThermalForcing"][:] = layerTF
f.variables["ismip6shelfMelt_zOcean"][:] = layerDepths

print(resultsList) 
print('polynomial coefficients: {}'.format(pz)) 
plt.show() 

# Update history attribute of netCDF file
thiscommand = datetime.now().strftime("%a %b %d %H:%M:%S %Y") + ": " + " ".join(sys.argv[:])
thiscommand = thiscommand+"; Converted 2D thermal forcing to 3D thermal forcing using CTD files: " + ",".join(resultsList)
if hasattr(f, 'history'):
   newhist = '\n'.join([thiscommand, getattr(f, 'history')])
else:
   newhist = thiscommand
setattr(f, 'history', newhist )


f.close()