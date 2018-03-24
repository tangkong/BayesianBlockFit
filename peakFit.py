from scipy.optimize import curve_fit
from scipy import integrate, signal
from scipy.special import wofz
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

#import personal modules
from peakShapes import voigtFn, gaussFn

def peakFit(data, LDatum, RDatum, peakShape, numCurves,
             savePath = None, filename = None):
    '''
    Peak fitting function.  Fits with ?various? functions?  
    Implements custom functions for allowing multiple peaks of same shape
    Attempts to plot at the end
    
    Input: Data array [x, y], complete data set
    Input: Peak bound indices
    Input: Peak shape: {Gaussian, Voigt}
    Input: numCurves = number of curves to fit
    
    Output: ndarray of optimized parameters. 
        result[i] = opt. param of curve i
        FWHM calculated depending on curve type
    ''' 

    # Convert left and right bounds to integers for indexing
    LDat = int(LDatum)
    RDat = int(RDatum)
    domain = np.arange(LDat, RDat)
    # print('LDat: {0}, RDat: {1}'.format(LDat, RDat))
    
    # Separate x and y data series
    xData = data[0]
    yData = data[1]

    # Locate position of max in domain 
    maxInd = np.where(yData[domain] == np.max(yData[domain]))[0][0]
    # Shift according to LDat since data is complete 
    loc = int(LDat + maxInd)  # given all data
    # print('maxInd: {0}, loc: {1}'.format(maxInd, loc))
    
    #############################################################################
    ### Define various peak shape functions
    #############################################################################
    
    xRange = xData[RDat] - xData[LDat]

    if peakShape == 'Voigt':
        func = voigtFn
        # Initialize guess params and bounds for number of curves
        boundUpper = []
        boundLower = []
        guess = []
        
        # parameter array: [x0, y0, Intensity, alpha, gamma]
        for i in range(numCurves): # for each curve
            # Space peak locations evenly
            xPosGuess = (xData[LDat] + xRange * (i+1) / (numCurves+1))
            
            guess += [xPosGuess, np.min(yData[domain]), 
                        np.max(yData[domain]) - np.min(yData[domain]),
                         xRange / 10 / numCurves, xRange / 10 / numCurves]
            #guess += [xPosGuess, np.min(yData[domain]), 1, 0.05, 0.05] 

            # concatenate lists for bounds
            boundLower += [xData[0], np.min(yData), 0, 0, 0]
            boundUpper += [xData[-1], np.inf, np.inf, np.inf, np.inf]       
            
        # Combine bounds into tuple for input
        bounds = tuple([boundLower, boundUpper])
        print(bounds)
        
        try:
        # Curve fit function call using guess and bounds
            popt, pcov = curve_fit(voigtFn, xData[domain], yData[domain], 
                                    bounds=bounds, p0=guess)
        except RuntimeError as e:
            print(e) 
            popt = np.array(guess)
        
        # Calculate FWHM for each peak fit
        FWHM = []
        c0 = 2.0056
        c1 = 1.0593
        for i in range(0, len(popt), 5): # grab fwhm for each peak
            fg = 2*popt[i+3]*np.sqrt(2*np.log(2))
            fl = 2*popt[i+4]
            phi = fl / fg
            FWHM.append(fg * (1-c0*c1 + np.sqrt(phi**2 + 2*c1*phi + (c0*c1)**2)))

    
    elif peakShape == 'Gaussian':############################################
        func = gaussFn
        # Initialize guess params and bounds for number of curves
        boundUpper = []
        boundLower = []
        guess = []
        
        # parameter array: [x0, y0, Intensity, sigma]
        for i in range(numCurves): # for each curve
            # Space peak locations evenly
            xPosGuess = (xData[LDat] + xRange * (i+1) / (numCurves+1))
            
            guess += [xPosGuess, np.min(yData[domain]), 
                        np.max(yData[domain]) - np.min(yData[domain]), 
                        xRange/10/numCurves]
            
            # concatenate lists for bounds
            boundLower += [xData[0], np.min(yData), 0, 0]
            boundUpper += [xData[-1], np.inf, np.inf, np.inf]       
            
        # Combine bounds into tuple for input
        bounds = tuple([boundLower, boundUpper])

        try:
            # Curve fit function
            popt, pcov = curve_fit(gaussFn, xData[domain], yData[domain], 
                                    bounds=bounds, p0=guess)
        except RuntimeError as e:
            print(e)
            popt = np.array(guess)

        # Calculate FWHM for each peak fit
        FWHM = []
        for i in range(0, len(popt), 4): # grab fwhm for each peak
            FWHM.append(2*popt[i+3]*np.sqrt(2*np.log(2)))
   
    ###########################################################################
    ### Calculate full width half maximum of given peaks
    ###########################################################################

    ###########################################################################
    ### Plotting and saving
    ###########################################################################
    if (savePath != None) and (filename != None):
        plt.figure(figsize=(8,8))
    
    # Organize final parameters into array
    finalParams = []
    for j in range(numCurves):   # Plot each individual curve
        L = 0 + j * len(popt) / numCurves  # Sep popt array
        R = (j+1) * len(popt) / numCurves

        finalParams.append(popt[L:R])
        
        if (savePath != None) and (filename != None): # Plot setup
            plt.plot(xData[domain], func(xData[domain], *guess[L:R]), 
                    '--', alpha=0.5, label='guessed curve: ' + str(j))
            plt.plot(xData[domain], func(xData[domain], *popt[L:R]), '.', alpha=0.5, 
                label='opt. curve {:.0f} FWHM = {:.3f}'.format(j, FWHM[j]))
        
    # Finish plotting 
    if (savePath != None) and (filename != None):
        plt.plot(xData[domain], yData[domain], marker='s', color='k', label='data')
        plt.plot(xData[domain], func(xData[domain], *popt), 
                    color='r', label='combined data')
        #plt.plot(xData[loc], yData[loc], marker='o', label='max point') # Max position
        plt.legend()

        plt.savefig(savePath + str(filename) + 'peakAt_' + 
                        '{:.3f}'.format(xData[loc]) + '.png')

        #plt.text(np.max(xData[domain]), np.max(yData[domain]), r'FWHM = {0}'.format())
        
        plt.close()
    
        print('Plot generated for peak at: {:.3f}'.format(xData[loc]))
    return finalParams, FWHM
