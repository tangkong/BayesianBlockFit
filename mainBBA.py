import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from os.path import basename

# Personal package import
from BlockData import BlockData
from peakFitResidIter import peakFit
from bumpFindFit import bumpFindFit
from reportFn import genOptParamCSV, genLitFWHMCSV

import time

##############################################################
##############INPUT HERE######################################

#path = os.path.expanduser('~/data/bl1-5/Nov2017/CoTaZr/data/J4/Processed/')
path = os.path.expanduser('~/Bayesianbumps/JaeProcPeakTst/')
savePath = path + '_peak_detection/'

peakShape = 'Voigt'
numCurves = 2 
fit_order = 2 
##############End Input#######################################
##############################################################


if not os.path.exists(savePath):
    os.makedirs(savePath)

files = glob.glob(os.path.join(path, '*.csv'))
if len(files) == 0:
    print('No files found')

plt.close('all')

fileGen = (x for x in files if basename(x)[-5]=='D')
loopTime = []
for file in fileGen:
    start = time.time()
    peakCnt = 0
    # File data into array
    print('============================================================')
    print file
    data = np.genfromtxt(file, delimiter = ',')
    Qlist = data[:,0]
    IntAve = data[:,1] 
    dataArray = np.array([Qlist, IntAve])

    #### Data Structure object instantiation (data, fit_order, ncp_prior)
    dataIn = BlockData(dataArray, fit_order, 0.5)
    #### has various functions
    
    # background subtracted with polynomial of order = fit_order, trims ends
    dataIn.bkgdSub(fit_order=fit_order) 
    #dataIn.trimSubData() # Trim off some values (taken from original script)
    
    # Plot bkgdSub Data
    plt.figure(figsize=(8,8))
    plt.plot(Qlist, IntAve, label='Raw data', marker='s', color='k')
    plt.plot(dataIn.subData[0], dataIn.subData[1],
                 label='Background subtracted', color='r')
    plt.plot(dataIn.subData[0], np.polyval(dataIn.fitCoeff, dataIn.subData[0]), 
                '--', label='Background', alpha=0.5, color='k')
    plt.legend()
    plt.savefig(savePath + basename(file)[:-7] + '_plot.png')
    plt.close()

    # Guess at noise level
    hld = dataIn.subData
    sigmaGuess = np.std(hld[1][hld[1] <= np.median(hld[1])])

    dataIn.cellData = sigmaGuess * np.ones(len(dataIn.subData[0]))
    
    # incorporate block information into data struct
    dataIn.blockFinder()
    
    
    # Get optimized parameters from fitting each block and plot
    optParams, litFWHM = bumpFindFit(dataIn, peakShape, numCurves, 
                            savePath, basename(file)[:-7])
    
    # Print information to terminal, print data to csv
    print('Fit ({0}) curve(s) per peak'.format(optParams['numCurves']))
    print('Using ({0}) peak shape'.format(optParams['peakShape']))

    genOptParamCSV(savePath, file, optParams)

    genLitFWHMCSV(savePath, file, litFWHM)
    end = time.time()
    
    loopTime += [(end-start)]

# Evaluate performance
avgTime = np.mean(loopTime)
maxTime = np.max(loopTime)
print('============================================================')
print('============================================================')
print('Files finished processing') 
print('-----Average {:.4f}s / file, max {:.4f}s / file'.format(avgTime, maxTime))
print('-----Total Time Elapsed {:.4f}s'.format(np.sum(loopTime)))
print('============================================================')
