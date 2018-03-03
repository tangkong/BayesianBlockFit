import numpy as np
from peakFitResidIter import peakFit, calcFWHM

def bumpFindFit(dat, peakShape, numCurves, savePath = None, filename = None):
    '''
    Separate data into peaks and fit them.  Plot fitted peaks and plot output.
    *Use hill climbing to find local max
    *Incorporate blocks that climbed to bump as peaks
    *Fit bumps using desired curve and number per peak
    input: data array (x,y), desired peak shape, num curves per peak, save info
    output: peak data dictionary
            #: Array of arrays [numCurves x numParamsPerCurve]
    '''

    numBlocks = len(dat.rateVec)
    # add beginning and ending changePoints

    cpUse = np.concatenate(([1], dat.changePoints, [len(dat.subData[1])]))
    dat.changePoints = cpUse

    numCP = len(cpUse) - 1
    idLeftVec = np.zeros(numCP)
    idRightVec = np.zeros(numCP)

    for i in range(numCP):
        idLeftVec[i]  = cpUse[i]
        idRightVec[i] = cpUse[i+1]

    # Find maxima defining watersheds, scan for 
    # highest neighbor of each block
    idMax = np.zeros(numBlocks)
    for j in range(numBlocks):
        jL = (j-1)*(j>0) + 0*(j<=0) # prevent out of bounds
        jR = (j+1)*(j<(numBlocks-1)) + (numBlocks-1)*(j>=(numBlocks-1))

        rateL = dat.rateVec[jL]
        rateC = dat.rateVec[j]
        rateR = dat.rateVec[jR]
        rateList = [rateL, rateC, rateR]

        jMax = np.argmax(rateList) #get direction [0, 1, 2]
        idMax[j] = j + jMax - 1 # convert direction to delta

    idMax[ idMax > numBlocks] = numBlocks
    idMax[ idMax < 0] = 0

    # Implement hill climbing (HOP algorithm)
    hopIndex = np.array(range(numBlocks)) # init: all blocks point to self
    hopIndex = hopIndex.astype(int)  # cast all as int
    ctr = 0
    while ctr <= 10000000:
        newIndex = idMax[hopIndex]  # Point each to highest neighbor

        if np.array_equal(newIndex, hopIndex):
            break
        else:
            hopIndex = newIndex.astype(int)
            
        ctr += 1
        
        if ctr == 10000000:
            print('Hill climbing did not converge...?')

    idMax = np.unique(hopIndex)
    numMax = len(idMax)

    # collect data for each bump (peak).
    iCntMax = 0

    paramDict = {'numCurves': numCurves, 'peakShape': peakShape}
    litFWHM = {}
    for k in range(numMax):
        currMax = idMax[k] 

        # find info on current max block
        maxBumpInd = np.fix((idLeftVec[currMax] 
                             + idRightVec[currMax]) / 2) # center of block
        currVec = np.where(hopIndex == currMax)[0] # vector of blocks contributing to currMax
        currNumBlocks = len(currVec)
        
        # indices for start and end of bump
        leftDatum  = idLeftVec[currVec[0]] - 1
        rightDatum = idRightVec[currVec[currNumBlocks-1]] - 1
       
        ############################################################################## 
        # Get parameters from peak fit
        if (savePath != None) and (filename != None):
            optParam, FWHM = peakFit(dat.subData, leftDatum, rightDatum, 
                                peakShape, numCurves, savePath, filename)
        else:
            optParam, FWHM = peakFit(dat.subData, leftDatum, rightDatum, 
                                peakShape, numCurves)
        
        # Append FWHM to optimized parameters 
        for j in range(len(optParam)):
            optParam[j] = np.append(optParam[j], FWHM[j])

        ############################################################################## 
        # add to dictionary
        paramDict[k] = optParam

        
        ############################################################################## 
        # Calculate FWHM in a literal sense

        litFWHM[k] = calcFWHM(dat.subData, leftDatum, rightDatum)
        
    return paramDict, litFWHM
