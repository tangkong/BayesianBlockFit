from scipy.optimize import curve_fit, basinhopping
from scipy import integrate, signal
from scipy.special import wofz
from scipy.interpolate import UnivariateSpline, interp1d
from numpy.polynomial.chebyshev import chebval
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
from os.path import basename

#import personal modules
from peakShapes import voigtFn, gaussFn

class BlockData:
    '''
    Data wrapper object holding information for bayesian block analysis
    '''   
    
    def __init__(self, dataArray, fit_order, ncp_prior, peakShape):
        '''
        Initializes BlockData object with data array and basic parameters
        '''
        self.data = np.array(dataArray)
        self.fit_order = fit_order
        self.ncp_prior = ncp_prior
        self.peakShape = peakShape
        
    def bkgdSub(self, **kwarg):
        '''
        Perform modified chebyshev polynomial subtraction and penalize with cost fn
        Offset data to be positive
    
        Adapted from Fang Ren segmentation program
        '''
        
        # Truncate data 
        trimDataX = self.data[0, 50:-50]
        trimDataY = self.data[1, 50:-50]

        def downsamp(size, x, y):
            '''
            downsample data based on size parameter.  Return downsampled data
            '''
            yNew = y[::size] 
            xNew = x[::size] 
            
            return xNew, yNew
        
        def chebFunc(x, *params):
            '''
            Modified chebyshev function with 1/x 
            '''
            params = params[0]
            y = chebval(x, params[:4])
            E = params[4]
            y = y + E/x
            return y


        def objFunc(*params):
            '''
            implement cost function for modified chevbyshev function
            varies cost to emphasize low Q region
            Scoping pulls trimmed x, y data from surrounding scope

            return: cost of fit 
            '''
            
            params = params[0]
            J = 0
            fit = chebval(X, params[:4])
            E = params[4]
            fit = fit + E / X

            for i in range(len(Y)):
                if X[i] < 1:    # Treat low Q equally
                    J = J + (Y[i] - fit[i])**4
                else:
                    if Y[i] < fit[i]:
                        J = J + (Y[i] - fit[i])**4
                    if Y[i] >= fit[i]:
                        J = J + (Y[i] - fit[i])**2
            return J
        
        # Create a sparse data set for fitting
        X, Y = downsamp(20, trimDataX, trimDataY)

        x0 = [1,1,1,1,1]
        result = basinhopping(objFunc, x0)
        bkgd_sparse = chebFunc(X, result.x)
        # create function that interpolates sparse bkgd
        f = interp1d(X, bkgd_sparse, kind='cubic', bounds_error=False)

        # expressed background values
        bkgd = f(trimDataX)

        subDataY = trimDataY - bkgd
        if np.min(subDataY) < 0:
            subDataY = subDataY + np.absolute(np.min(subDataY))
    
        # Dump any nan values in all data
        finalDataY = subDataY[~np.isnan(subDataY)]
        finalDataX = trimDataX[~np.isnan(subDataY)]
        bkgd = bkgd[~np.isnan(subDataY)]

        # Save background subtracted data
        self.subData = np.array([finalDataX, finalDataY])
        self.bkgd = bkgd
        
        print 'background subtraction completed'

    def bkgdSubPoly(self, **kwarg):
        '''
        Perform crude polynomial background subtraction and 
        offset data to be positive
        '''
        if 'fitOrder' in kwarg:
            self.fit_order = fitOrder
            
        self.fitCoeff = np.polyfit(self.data[0], self.data[1], self.fit_order)
        fitDataX = self.data[0, 50:-50]
        fitDataY = self.data[1, 50:-50] - np.polyval(self.fitCoeff, 
                                                      self.data[0, 50:-50])
        if np.min(fitDataY) < 0:
            fitDataY = fitDataY + np.absolute(np.min(fitDataY))
        self.subData = np.array([fitDataX, fitDataY])
        self.bkgd = np.polyval(self.fitCoeff, self.data[0, 50:-50])
            
        print 'background subtraction completed'
            
    def trimSubData(self,trimLen=50):
        '''
        Trim 50 data points from both ends of data.
        '''
        self.subData = np.array([self.subData[0,trimLen:-trimLen], 
                                    self.subData[1,trimLen:-trimLen]])

    def blockFinder(self):
        '''
        Performs Bayesian Block finding on subData and records change point array
        Algorithm description is hard....
        '''
        data_mode = 3
        numPts = len(self.cellData)
        tt = np.arange(numPts)
        nnVec = []

        ## To implement: trimming/choosing of where to start/stop
        ## To implement: timing functions

        cp = []
        cntVec = []

        iterCnt = 0
        iterMax = 10

        while True:
            best = []
            last = []

            for r in range(numPts):
                # y-data background subtracted
                sumX1 = np.cumsum(self.subData[1,r::-1]) 
                sumX0 = np.cumsum(self.cellData[r::-1]) # sigma guess
                fitVec = (np.power(sumX1[r::-1],2) / (4 * sumX0[r::-1]))

                paddedBest = np.insert(best,0,0)
                best.append(np.amax( paddedBest + fitVec - self.ncp_prior ))
                last.append(np.argmax( paddedBest + fitVec - self.ncp_prior ) )

                # print('Best = {0},  Last = {1}'.format(best[r], last[r]))

            # Find change points by peeling off last block iteratively
            index = last[numPts-1]

            while index > 0:
                cp = np.concatenate( ([index], cp) )
                index = last[index-1]

            # Iterate if desired, to implement later
            break

        numCP = len(cp)
        numBlocks = numCP + 1

        rateVec  = np.zeros(numBlocks)
        numVec   = np.zeros(numBlocks)
        dtVec    = np.zeros(numBlocks)
        tt1Vec   = np.zeros(numBlocks)
        tt2Vec   = np.zeros(numBlocks)

        cptUse = np.insert(cp, 0, 0)

        #print('cptUse start: {0}, end: {1}, len: {2}'.format(cptUse[0],
                                                        #cptUse[-1], len(cptUse)))
        #print('lenCP: {0}'.format(len(cp)))
        
        # what does this stuff do I don't know...
        print('numBlocks: {0}'.format(numBlocks))
        for idBlock in range(numBlocks):
            # set ii1, ii2 as indexes.  Fix edge case at end of blocks
            ii1 = int(cptUse[idBlock])
            if idBlock < (numBlocks - 1):
                ii2 = int(cptUse[idBlock + 1] - 1)
            else:
                ii2 = int(numPts)
                    
            if data_mode == 3:
                # print('ii1 = {0}, ii2 = {1}'.format(ii1, ii2))
                subset = self.subData[1,ii1:ii2]
                weight = self.cellData[ii1:ii2]
                if ii1 == ii2:
                    subset = self.subData[1,ii1]
                    weight = self.cellData[ii1]
                rateVec[idBlock] = np.dot(weight, subset) / np.sum(weight)

                if np.sum(weight) == 0:
                    print('error, divide by zero at index: {0}'.format(idBlock))
                    print('-------ii1: {0}, ii2: {1}'.format(ii1, ii2))
            else:
                print('not yet implemented')
                
        self.changePoints = cp
        self.numVec       = numVec
        self.rateVec      = rateVec
        self.best         = best
        self.last         = last
        self.nn           = nnVec

        self.cptUse       = cptUse
        return None # end blockFinder method

    def genResidPlot(self, savePath, file):
        '''
        Generate full plot with peaks and residuals.
        Input: Save path, filename
        Input: optParams dictionary

        uses subData as data, 
        '''

        fitYData = np.array([])
        xData = self.subData[0]
        yData = self.subData[1]

        if self.peakShape == 'Voigt':
            func = voigtFn
        elif self.peakShape == 'Gaussian':
            func = gaussFn
        
        # Mega bandaid...to figure out indexing later
        temp = self.peakDomains[len(self.peakDomains)-1]
        self.peakDomains[len(self.peakDomains)-1] = (temp[0], temp[1]+1)
    
        # Stitch together fitted peaks
        for i in range(len(self.peakDomains)):
            LDat, RDat = self.peakDomains[i]

            domain = np.arange(int(LDat), int(RDat))
            xRange = xData[domain]

            # Get relevent fit params
            params = self.optParams[i]
            # Flatten / unroll fit params
            flatParams = np.concatenate(params).tolist()
            
            # concat data for each peak, all curves combined
            fitYData = np.append(fitYData, func(xData[domain], *flatParams))

        # Generate residual data array
        resid = yData - fitYData
        RSS = np.sum(resid**2)
        
        pctErr = [] 
        # Get average percent error in each block
        for k in range(len(self.peakDomains)):
            LDat, RDat = self.peakDomains[k]

            domain = np.arange(int(LDat), int(RDat))
            
            error = np.mean(np.absolute(resid[domain]) / (yData[domain]+1)) * 100
            pctErr += [error]

        # Store pctErr for each peak
        self.pctErr = pctErr
        # Plot all and save ######################################################
        # Plot data and fit
        fig = plt.figure(figsize=(10,8))
        frame1 = fig.add_axes((0.1, 0.3, 0.8, 0.6))
        plt.plot(xData, yData, 'sk', label='data')
        plt.plot(xData, fitYData, '-r', label='fit')

        # vertical lines at plot boundaries
        for j in range(len(self.peakDomains)):
            xLowerLoc = xData[int(self.peakDomains[j][0])]
            plt.axvline(x=xLowerLoc, linestyle='--', color='k')
            plt.text(xLowerLoc, np.max(yData), '{:.3f}%'.format(pctErr[j]))

        plt.legend()
        plt.grid()
        frame1.set_xticklabels([])

        # Residual plot
        frame2 = fig.add_axes((0.1,0.1,0.8,0.2))
        plt.plot(xData, resid, '-b', label='residual')
        # vertical lines at plot boundaries
        for j in range(len(self.peakDomains)):
            plt.axvline(x=xData[int(self.peakDomains[j][0])],
                            linestyle='--', color='k')
        plt.legend()
        plt.grid()
        plt.savefig(savePath + str(basename(file)[:-7]) + '_residPlot.png')
        
        plt.close()
        
        return pctErr 
    
