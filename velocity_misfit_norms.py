#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 16:11:59 2021
Use this script to compute 2-norm and inf-norm for velocity mismatch between
 MALI diagnostic solve and Albany optimization
@author: trevorhillebrand
"""
from netCDF4 import Dataset
import numpy as np
import argparse
import sys
import matplotlib.pyplot as plt
import scipy.stats

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.description = __doc__
parser.add_argument("-r", "--run", dest="runFiles", help="name of file with result from forward run", metavar="FILENAME")
parser.add_argument("-o", "--obs", dest="observationsFile", help="name of file with observed surfaceSpeed", metavar="FILENAME")
parser.add_argument("-t", "--timeLev", dest="timeLev", help="time level to compare with observations")
args = parser.parse_args()

timeLev = args.timeLev
scyr = 3600.0*24.0*365.0
runFilesList = args.runFiles.split(',') #parse variable names into a list

runNamesList = []

for runFile in runFilesList:
    runNamesList.append(runFile.replace('/output.nc', ''))

nFiles = len(runFilesList)
obs = Dataset(args.observationsFile, 'r')
obsSpeed = obs.variables['surfaceSpeed'][0,:]*scyr
obsSpeedMax = np.nanmax(obsSpeed)
obsMask = (obsSpeed > 0.0)
areaCell = obs.variables["areaCell"][:]
nEdges = obs.dimensions["nEdges"].size
nCells = obs.dimensions["nCells"].size
xCell = obs.variables["xCell"][:]
yCell = obs.variables["yCell"][:]
xEdge = obs.variables["xEdge"][:]
yEdge = obs.variables["yEdge"][:]
edgesOnCell = obs.variables["edgesOnCell"][:]
cellsOnCell = obs.variables["cellsOnCell"][:]
cellsOnEdge = obs.variables["cellsOnEdge"][:]

obsMarginEdgeMask = xEdge * 0
obsMarginEdgeMask = obsMarginEdgeMask.astype(bool)
modelMarginEdgeMask = xEdge * 0
modelMarginEdgeMask = modelMarginEdgeMask.astype(bool)

for iEdge in np.arange(0, nEdges):
    if np.sum(obsMask[cellsOnEdge[iEdge]-1]) == 1:  #that is, if one cell on the edge is in the ice mask and the other isn't, then this edge is the margin
        obsMarginEdgeMask[iEdge] = True

# Customize plotting
plotMarkers = {'MIROC5':'o', 'HadGEM2':'+', 'CNRM':'*'}

speedColors = ['tab:blue', 'tab:green', 'tab:purple', 'black']
bedFile = Dataset(runFilesList[0], 'r')
bed = bedFile.variables["bedTopography"][:]
bedFile.close()

normFig,normAx = plt.subplots(2,1, sharex=True)
speedFig,speedAx = plt.subplots(1,1)
runFileCount = 0
scores = {} #empty dictionary to fill with scores for each run

#define speed bands
speedBandsLower = [0., 300., 600., 0.]
speedBandsUpper = [300., 600., round(obsSpeedMax), round(obsSpeedMax)]
nBands = len(speedBandsLower)
marginScores = [] #fill with RMSdistanceToObsMargin below

for runFile in runFilesList:
    run = Dataset(runFile, 'r')
    modelSpeed = run.variables['surfaceSpeed'][timeLev,:]*scyr
    modelMarginEdgeMask[:] = False #start with empty mask each time
    distanceToObsMargin = xEdge * 0.0 #start with all zeros 
    modelMask = (modelSpeed > 0.0)
    colorIndex = 0
    scores[runFile] = {} #fill these with all the scores below
    scores[runFile]['mismatchWeighted'] = []
    scores[runFile]['mismatch'] = []
    scores[runFile]['twoNorm'] = []
    scores[runFile]['infNorm'] = []
    scores[runFile]['RMS'] = []

    for speed1, speed2 in zip(speedBandsLower, speedBandsUpper):
        mask = (obsSpeed>speed1) * (obsSpeed<speed2) * (modelSpeed>0.0) #mask using observed speeds and modeled extent
        #mismatchWeighted = (modelSpeed - obsSpeed) * mask * areaCell/np.max(areaCell) # do not penalize here for wrong ice extent. That will be scored separately
        mismatch =  (modelSpeed - obsSpeed) * mask
        mismatchWeighted = mismatch / np.median(obsSpeed[mask])
        twoNorm = np.sqrt(np.nansum(mismatchWeighted**2))
        infNorm = np.nanmax(np.abs(mismatch))
        RMS = ((mismatch**2).sum()/mask.sum())**0.5

        scores[runFile]['mismatchWeighted'].append(mismatchWeighted)
        scores[runFile]['mismatch'].append(mismatch)
        scores[runFile]['twoNorm'].append(twoNorm)
        scores[runFile]['infNorm'].append(infNorm)
        scores[runFile]['RMS'].append(RMS)
        #print("For speeds > {} m/yr".format(speed))
        #print("Speed misfit 2-Norm = {}".format(twoNorm))
        #print("MnormAx abs speed misfit = {}".format(infNorm))
        #print("RMS speed misfit = {}".format(RMS))
        #normAx[0].scatter(runFileCount, RMS, c=speedColors[colorIndex])
        #normAx[1].scatter(runFileCount, infNorm, c=speedColors[colorIndex])
        colorIndex += 1
     #loop over edges to find margins
    for iEdge in np.arange(0, nEdges):
        #if one cell on the edge is in the ice mask and the other isn't, and bed is below sea sea level. then this edge is the margin
         if (np.sum(modelMask[cellsOnEdge[iEdge]-1]) == 1) and (bed[0,cellsOnEdge[iEdge]-1] < 0.0).any():  
                 modelMarginEdgeMask[iEdge] = True 
                 distanceToObsMargin[iEdge] = np.min(np.sqrt( (xEdge[iEdge] - xEdge[obsMarginEdgeMask])**2 + (yEdge[iEdge] - yEdge[obsMarginEdgeMask])**2 ))
    distanceToObsMargin = distanceToObsMargin[modelMarginEdgeMask] # only use distances corresponding to actual model margin edges to avoid adding lots of extra zeros
    RMSdistanceToObsMargin = ((distanceToObsMargin**2).sum()/ modelMarginEdgeMask.sum())**0.5
    maxDistanceToObsMargin = np.nanmax(distanceToObsMargin)
    scores[runFile]['RMSdistToMargin'] = RMSdistanceToObsMargin
    marginScores.append(RMSdistanceToObsMargin) # for ranking below
    scores[runFile]['maxDistToMargin'] = maxDistanceToObsMargin
    normAx[1].scatter(runFileCount, RMSdistanceToObsMargin, c='black')
    #marginAx[1].scatter(runFileCount, maxDistanceToObsMargin, c='black')
    runFileCount += 1
    run.close()    

# Now loop through runs again and rank scores
bandScores = []
bandArea = []
bandRankings = []
fileScores = []
marginRankings = scipy.stats.rankdata(marginScores) 
for band in np.arange(0, nBands):
    for runFile in runFilesList:
        bandScores.append(scores[runFile]['RMS'][band])
    bandArea.append(np.sum(areaCell * (obsSpeed>speedBandsLower[band]) * (obsSpeed<speedBandsUpper[band])))
for band in np.arange(0, nBands):
    bandScoresTmp = bandScores[band*nFiles:((band+1) * nFiles)] # just the scores in this band
    bandScoresNormalized = bandScoresTmp / np.max(bandScoresTmp)
    #plot normalized RMS scores
    normAx[0].scatter(np.arange(0,nFiles)+band/25, bandScoresNormalized, c=speedColors[band],
                      marker='.', label=(str(speedBandsLower[band]) 
                      + ' - ' + str(speedBandsUpper[band]) + ' m yr$^{-1}$'))
    speedAx.scatter(np.arange(0,nFiles)+band/25, bandScoresTmp / np.mean([speedBandsLower[band], speedBandsUpper[band]]), c=speedColors[band],
                    marker='.', label=(str(speedBandsLower[band])
                    + ' - ' + str(speedBandsUpper[band]) + ' m yr$^{-1}$'))
    # Uncomment if you want to normalize by band area
    #speedAx.scatter(np.arange(0,nFiles)+band/25, bandScoresTmp / bandArea[band], c=speedColors[band],
    #                marker='.', label=(str(speedBandsLower[band])
    #                + ' - ' + str(speedBandsUpper[band]) + ' m yr$^{-1}$'))
    bandRankings.append(scipy.stats.rankdata(bandScores[band*nFiles:((band+1) * nFiles)]))

for iFile in np.arange(0, nFiles):
    fileScoreTmp = 0 
    for band in np.arange(0, nBands):
        fileScoreTmp += bandRankings[band][iFile]
    fileScoreTmp += marginRankings[iFile]
    fileScores.append(fileScoreTmp)

normLeg = normAx[0].legend()
normLeg.get_texts()[-1].set_text('>' + str(speedBandsLower[-1]) + ' m yr$^{-1}$')
normAx[0].set_title('MIROC5, q=1/5')
normAx[1].set_xlabel("$\sigma_{max}$ (kPa)")
#normAx[1].set_xlabel('Basal friction law exponent')
#normAx[0].set_ylabel('RMS velocity misfit\n(weighted by cell area)')
normAx[0].set_ylabel('Normalized RMS\nvelocity misfit')
normAx[1].set_xticks(np.arange(0,runFileCount))
#normAx[0].set_xticks(np.arange(0,runFileCount))
#normAx[0].set_xticklabels(['1', '1/3', '1/10', '1/20', '1/25', '1/30', '1/40', '1/50'])
#normAx[1].set_xticklabels(['150', '160', '170', '180', '190', '200'])
#normAx[1].set_xticklabels(runNamesList, rotation=15)  # if you want to label by the file path
normAx[1].set_ylabel('RMS distance to\nobserved margin (m)')
#normAx[1].set_ylabel('Max velocity misfit (m yr$^{-1}$)')
normFig.set_size_inches(4,7)
normFig.tight_layout()

#marginAx[0].set_ylabel('RMS distance to\nobserved margin (m)')
#marginAx[1].set_ylabel('Max distance to\nobserved margin (m)')
#marginAx[1].set_xticks(np.arange(0,runFileCount))
#marginAx[1].set_xticklabels(runNamesList, rotation=15)
#marginFig.tight_layout()
speedLeg = speedAx.legend()
speedLeg.get_texts()[-1].set_text('All velocities')
#speedAx.set_xlabel('basal friction law exponent')
speedAx.set_ylabel('RMS velocity misfit (m yr$^{-1}$)')
speedAx.set_xticklabels(['', '1', '1/3', '1/4', '1/5', '1/6', '1/7', '1/8', '1/9', '1/10'])
speedFig.set_size_inches(8, 4)
speedFig.tight_layout()
#speedFig.savefig('tuneExponent', dpi=300, bbox_inches='tight')
#normFig.savefig('tuneMIROC5m10', dpi=300, bbox_inches='tight')

print('file names: {}'.format(runFilesList))
print('file scores: {}'.format(fileScores))
plt.show()
