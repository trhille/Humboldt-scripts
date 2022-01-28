#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
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

parser = OptionParser(description=__doc__)
parser.add_option("-f", dest="filenames", help="path to globalStats.nc files (strings separated by commas; no spaces)", metavar="FILENAME")
parser.add_option("-v", dest="variables", help="variable(s) to plot, separated by commas", default = "volumeAboveFloatation,groundedIceArea")
parser.add_option("-l", dest="labels", help="labels for legend; one per file", default = None)

options, args = parser.parse_args()
filenames = options.filenames.split(',') # split ensemble directories into list
labels = options.labels.split(',')
variables = options.variables.split(',')

fig, ax = plt.subplots(1, len(variables), sharex=True, figsize=(8,4))

# make list indexing work if there is only one variable
if len(variables) == 1:
    ax = [ax]

changeOrderOfMag = []
colors = ['tab:pink', 'tab:green', 'tab:blue', 'tab:purple']
nFile = 0
for label,filename in zip(labels,filenames):
    f = Dataset(filename, 'r')
    f.set_auto_mask(False)

    for ii,variable in enumerate(variables):
        
        change = f.variables[variable][:] - f.variables[variable][0]
        units = f.variables[variable].units
        
        if filename == filenames[0]: #get order of magnitude of change from first file only
            changeOrderOfMag.append(int(np.floor(np.log10(np.max(np.abs(change)))))) # Get order of magnitude for plotting  

        ax[ii].plot(f.variables["daysSinceStart"][:] / 365. + 2007.,
                    change / 10**changeOrderOfMag[ii], color=colors[nFile], label=label)
        ax[ii].grid('on')
        ax[ii].set_xlim(left=2007, right=2050)
        ax[ii].set_xlabel('Year')
        if variable == "volumeAboveFloatation":
            ax[ii].set_ylabel('Total change in volume above floatation ($10^{' + str(changeOrderOfMag[ii]) + '}$ m$^3$)')
            ax[ii].set_ylim(top=0.5, bottom=-6)
        elif variable == "groundedIceArea":
            ax[ii].set_ylabel('Total change in grounded ice area ($10^{' + str(changeOrderOfMag[ii]) + '}$ m$^2$)')
            ax[ii].set_ylim(top=12, bottom=-8)
        else:
            ax[ii].set_ylabel(variable + ' ($10^{' + str(changeOrderOfMag[ii]) + '}$ $' + units + '$)') 

    f.close()
    nFile += 1

ax[-1].legend()
fig.subplots_adjust(wspace=0.4)
#fig.savefig('convergence_test', dpi=400, bbox_inches='tight')
plt.show()
