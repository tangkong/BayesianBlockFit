import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from os.path import basename
import csv

# Personal package import
from BlockData import BlockData
from peakFit import peakFit
from bumpFindFit import bumpFindFit

##############################################################
##############INPUT HERE######################################

#path = os.path.expanduser('~/data/bl10-2/Jan2018/HiTp/data/Test/')
path = os.path.expanduser('~/Bayesianbumps/JaeProcPeakTst/')
savePath = path + '_peak_detection/'

peakShape = 'Gaussian'
numCurves = 2 

##############End Input#######################################
##############################################################


if not os.path.exists(savePath):
    os.makedirs(savePath)

files = glob.glob(os.path.join(path, '*.csv'))
if len(files) == 0:
    print('No files found')

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

    #### Data Structure object instantiation (data, fit_order, ncp_prior)
    dataIn = BlockData(dataArray, 2, 0.5)
    #### has various functions
    
    dataIn.bkgdSub() # background subtracted with polynomial of order = fit_order
    dataIn.trimSubData() # Trim off some values (taken from original script)
    
    # Plot bkgdSub Data
    plt.figure(figsize=(8,8))
    plt.plot(Qlist, IntAve, label='raw data')
    plt.plot(dataIn.subData[0], dataIn.subData[1], label='Background Subtracted')
    plt.savefig(savePath + basename(file)[:-7] + '_plot.png')
    plt.legend()
    plt.close()

    # Guess at noise level
    hld = dataIn.subData
    sigmaGuess = np.std(hld[1][hld[1] <= np.median(hld[1])])

    dataIn.cellData = sigmaGuess * np.ones(len(dataIn.subData[0]))
    
    # incorporate block information into data struct
    dataIn.blockFinder()
    
    
    # Get optimized parameters from fitting each block and plot
    optParams = bumpFindFit(dataIn, peakShape, numCurves, 
                            savePath, basename(file)[:-7])
    
    # Print information to terminal, print data to csv
    with open(savePath + basename(file)[:-7] + 'curveParams.csv', 'wb') as csv_file:
        writer = csv.writer(csv_file)
        print('Fit ({0}) curve(s) per peak'.format(optParams['numCurves']))
        print('Using ({0}) peak shape'.format(optParams['peakShape']))
        if optParams['peakShape'] =='Voigt':
            print('Array format: [x0, y0, I, alpha, gamma]')
            writer.writerow([ 'peak', 'curve', 'x0', 'y0', 'I', 'alpha', 'gamma'])
        elif optParams['peakShape'] == 'Gaussian':
            print('Array format: [x0, y0, I, sigma]')
            writer.writerow(['peak', 'curve', 'x0', 'y0', 'I', 'sigma'])
        for key, item in optParams.items():
            if key != 'numCurves' and key != 'peakShape':
                for i in range(len(item)):
                    print('Param array for peak {0}, curve {1}: {2}'.format(key, i+1, 
                                            np.array2string(item[i], precision=4)))

                    writer.writerow([key, i] + list(item[i]))
        
    break
