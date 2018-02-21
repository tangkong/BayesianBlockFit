import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from os.path import basename

from BlockData import BlockData
from peakFit import peakFit
from bumpFindFit import bumpFindFit

##############################################################
##############INPUT HERE######################################

path = os.path.expanduser('~/Bayesianbumps/JaeProcPeakTst/')
save_path = path + '_peak_detection/'

peakShape = 'Voigt'
numCurves = 2

##############End Input#######################################
##############################################################


if not os.path.exists(save_path):
    os.makedirs(save_path)

files = glob.glob(os.path.join(path, '*.csv'))

plt.close('all')

fileGen = (x for x in files if basename(x)[-5]=='D')
for file in fileGen:
    peakCnt = 0
    # File data into array
    print file
    data = np.genfromtxt(file, delimiter = ',')
    Qlist = data[:,0]
    IntAve = data[:,1] 
    dataArray = np.array([Qlist, IntAve])

    #### Data Structure object instantiation
    dataIn = BlockData(dataArray, 2, 0.5)
    #### has various functions
    
    dataIn.bkgdSub() # background subtracted with polynomial of order = fit_order
    dataIn.trimSubData() # Trim off some values (taken from original script)

    # Guess at noise level
    hld = dataIn.subData
    sigmaGuess = np.std(hld[1][hld[1] <= np.median(hld[1])])

    dataIn.cellData = sigmaGuess * np.ones(len(dataIn.subData[0]))
    
    # incorporate block information into data struct
    dataIn.blockFinder()
    
    
    # Get optimized parameters from fitting each block
    optParams = bumpFindFit(dataIn, peakShape, numCurves)
    
    # bumpFindFit() plots, so can save plot after finishing 
    # plt.savefig(save_path + basename(file)[:-7] + '_fits.png')
    
    # Print information to terminal (crude, will rework later)
    for key in optParams:
        if key == 'numCurves':
            print('Fit ({0}) curve(s) per peak'.format(optParams[key]))
        elif key == 'peakShape':
            print('Using ({0}) peak shape'.format(optParams[key]))
            if optParams[key]=='Voigt':
                print('Array format: [x0, y0, I, alpha, gamma]')
            elif optParams[key] == 'Gaussian':
                print('Array format: [x0, y0, I, sigma]')
        else:
            for i in range(len(optParams[key])):
                print('Param array for peak {0}, curve {1}: {2}'.format(key, i+1, optParams[key][i]))
