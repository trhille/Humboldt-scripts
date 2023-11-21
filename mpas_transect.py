#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 15:16:04 2020
Plot  bed and surface topography from MPAS netCDF along a transect defined in csv
@author: trevorhillebrand
"""

import numpy as np
import csv
from netCDF4 import Dataset
from optparse import OptionParser
from scipy.interpolate import LinearNDInterpolator
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

parser = OptionParser(description='Plot transect from MPAS netCDF')
parser.add_option("-d", "--data", dest="data_file", help="the MPAS netCDF file")
parser.add_option("-c", "--coords", dest="coords_file", help="csv file defining transect")
parser.add_option("--temperature", dest="interp_temp", action="store_true", help="interpolate temperature")
for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = parser.parse_args()

dataset = Dataset(options.data_file, 'r')
dataset.set_always_mask(False)
xCell = dataset.variables["xCell"][:]
yCell = dataset.variables["yCell"][:]
nCells = dataset.dimensions['nCells'].size
areaCell = dataset.variables["areaCell"][:]
layerThicknessFractions = dataset.variables["layerThicknessFractions"][:]
nVertLevels = dataset.dimensions['nVertLevels'].size

li_mask_ValueDynamicIce = 2
cellMask = dataset.variables['cellMask'][:]
cellMask_dynamicIce = (cellMask & li_mask_ValueDynamicIce) // li_mask_ValueDynamicIce
# only take thickness, speed, and temperature of dynamic ice
thk = dataset.variables["thickness"][:] * cellMask_dynamicIce
if "surfaceSpeed" in dataset.variables.keys():
    speed = dataset.variables["surfaceSpeed"][:] * cellMask_dynamicIce * 3600. * 24. * 365.
else:
    speed = np.sqrt(dataset.variables["uReconstructX"][:,:,0]**2. + 
                    dataset.variables["uReconstructY"][:,:,0]**2.) * cellMask_dynamicIce * 3600. * 24. * 365.

temperature = dataset.variables['temperature'][:]
    
bedTopo = dataset.variables["bedTopography"][0,:]

startYear = 2007
if "daysSinceStart" in dataset.variables.keys():
    time = dataset.variables["daysSinceStart"][:] / 365. + startYear
else:
    time = range(dataset.dimensions['Time'].size)

coords = []
x = []
y = []

with open(options.coords_file, newline='') as csvfile:
     reader = csv.reader(csvfile, delimiter=',')
     
     for row in reader:
         coords.append(row)
         x.append(float(row[0]))
         y.append(float(row[1]))

# increase sampling by a factor of 100
x_interp = np.interp(np.linspace(0, len(x)-1, 100*len(x)), np.linspace(0, len(x)-1, len(x)), x)
y_interp = np.interp(np.linspace(0, len(y)-1, 100*len(y)), np.linspace(0, len(y)-1, len(y)), y)

xArray = np.array(x_interp)
yArray= np.array(y_interp)

d_distance = np.zeros(len(xArray))
for ii in np.arange(1, len(xArray)):
    d_distance[ii] = np.sqrt( (xArray[ii] - xArray[ii-1])**2 + (yArray[ii] - yArray[ii-1])**2 )

distance = np.cumsum(d_distance)

transectFig, transectAx = plt.subplots(2,1, sharex=True, layout='constrained')
thickAx = transectAx[0]
thickAx.grid()
speedAx = transectAx[1]
speedAx.grid()
timeColors = cm.plasma(np.linspace(0,1,len(time)))

plt.rcParams.update({'font.size': 16})

bed_interpolant = LinearNDInterpolator(np.vstack((xCell, yCell)).T, bedTopo)
bed_transect = bed_interpolant(np.vstack((xArray, yArray)).T)

for timeSlice in np.arange(0, len(time)):

    thk_interpolant = LinearNDInterpolator(np.vstack((xCell, yCell)).T, thk[timeSlice,:])
    thk_transect = thk_interpolant(np.vstack((xArray, yArray)).T)
    lower_surf = np.maximum( -910. / 1028. * thk_transect, bed_transect)
    lower_surf_nan = lower_surf.copy()  # for plotting
    lower_surf_nan[thk_transect==0.] = np.nan
    upper_surf = lower_surf + thk_transect
    upper_surf_nan = upper_surf.copy()  # for plotting
    upper_surf_nan[thk_transect==0.] = np.nan
    thickAx.plot( (distance[-1]-distance)/1000., lower_surf_nan, color=timeColors[timeSlice])
    thickAx.plot( (distance[-1]-distance)/1000., upper_surf_nan, color=timeColors[timeSlice])
    
    speed_interpolant = LinearNDInterpolator(np.vstack((xCell, yCell)).T, speed[timeSlice,:])
    speed_transect = speed_interpolant(np.vstack((xArray, yArray)).T)
    speed_transect[thk_transect == 0.] = np.nan
    speedAx.plot( (distance[-1]-distance)/1000. , speed_transect, color=timeColors[timeSlice])

    if options.interp_temp:
        layer_thk = np.zeros((len(thk_transect), nVertLevels + 1))
        layer_midpoints = np.zeros((len(thk_transect), nVertLevels))
        layer_thk[:,0] = 0.
        for i in range(len(thk_transect)):
            layer_thk[i,1:] = np.cumsum(layerThicknessFractions * thk_transect[i])
            layer_midpoints[i,:] = upper_surf[i] - (layer_thk[i,1:] + layer_thk[i,0:-1]) / 2.

        temp_transect = np.zeros((len(xArray), nVertLevels))
        for lev in range(nVertLevels):
            print(f'Interpolating temperature for level {lev}')
            temp_interpolant = LinearNDInterpolator(np.vstack((xCell, yCell)).T, temperature[timeSlice,:,lev])
            temp_transect[:, lev] = temp_interpolant(np.vstack((xArray, yArray)).T)

thickAx.plot( (distance[-1]-distance)/1000., bed_transect, color='black')

if options.interp_temp:
    temp_transect[temp_transect == 0.] = np.nan
    temp_plot = thickAx.pcolormesh( np.tile( (distance[-1]-distance)/1000., (nVertLevels,1)).T,
                                    layer_midpoints[:,:], temp_transect[:,:], cmap='YlGnBu_r', vmin=240., vmax=273.15)

speedAx.set_xlabel('Distance (km)')
speedAx.set_ylabel('Surface\nspeed (m/yr)')
thickAx.set_ylabel('Elevation\n(m asl)')
transectFig.subplots_adjust(hspace=0.4)

cbar = plt.colorbar(temp_plot)
cbar.set_label('Temperature (K)')
cbar.ax.tick_params(labelsize=12)

plt.show()

dataset.close()
