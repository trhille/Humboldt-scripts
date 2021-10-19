#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import shutil
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
from optparse import OptionParser
import matplotlib.pyplot as plt
import matplotlib as mpl

dataPath =  '/Volumes/Storage2/OMG_CTD/podaac-tools.jpl.nasa.gov/drive/files/allData/omg/L2/CTD/CTD/2019/'
destPath = '/Users/trevorhillebrand/Documents/Greenland/OMG_CTD/Humboldt/'
boundLatMin = 79.0
boundLatMax = 79.9
boundLonMin = -70.0
boundLonMax = -64.0

CTDfiles = sorted(os.listdir(dataPath))

resultsList = []
catTemp = np.empty(0,)
catDepth = np.empty(0,)
plotDepth = np.linspace(-400,0, 400)

fig, ax = plt.subplots(1,1,)
ax.grid()
for CTDfile in CTDfiles:
   if '.nc' in CTDfile and '.md5' not in CTDfile:
      f = Dataset(dataPath + CTDfile, 'r')
      f.set_auto_mask(False)
      lat = f.variables["lat"][:]
      lon = f.variables["lon"][:]
      if (boundLatMin <= lat) and (boundLatMax >= lat) and (boundLonMin <= lon) and (boundLonMax >= lon):
         resultsList.append(CTDfile)
         shutil.copy(dataPath+CTDfile, destPath)
         temperature = f.variables["temperature"][0,:]
         depth = -1.0 * f.variables["depth"][0,:]
         keepIndices = np.where(temperature>-25.) # trim out bad temps
         temperature = temperature[keepIndices]
         depth = depth[keepIndices]
         catTemp = np.concatenate((catTemp, temperature), axis=0) #concatenate for polynomial fit
         catDepth = np.concatenate((catDepth, depth), axis=0)
         ax.scatter(temperature, depth)

      f.close()

pz = np.polyfit(catDepth, catTemp, 4)
ax.plot(np.polyval(pz, plotDepth), plotDepth, linewidth=2, color='black')
ax.set_ylabel('Depth (m)')
ax.set_xlabel('Temperature (deg. C)')

print(resultsList) 
print('polynomial coefficients: {}'.format(pz)) 
plt.show() 
