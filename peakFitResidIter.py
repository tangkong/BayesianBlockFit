from scipy.optimize import curve_fit
from scipy import integrate, signal
from scipy.special import wofz
from scipy.interpolate import UnivariateSpline
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
    ### Fit with specified peak shape.
    ### Populate guess paramters with residual
    #############################################################################
    
    # Initialize guess and bound arrays with correct shape for each function
    xRange = xData[RDat] - xData[LDat]
    if peakShape == 'Voigt':
        func = voigtFn
        # x0 and I values to be replaced during loop
        guessTemp = [0, np.min(yData[domain]), 0,
                        xRange / 10, xRange / 10 ]
        # x bounds to be replaced
        boundLowTemp = [xData[LDat], np.min(yData), 0., 0., 0.]
        boundUppTemp = [xData[RDat], np.inf, np.inf, np.inf, np.inf]

        boundLowerPart = [xData[0], np.min(yData), 0, 0, 0]
        boundUpperPart = [xData[-1], np.inf, np.inf, np.inf, np.inf]       
    elif peakShape == 'Gaussian':
        func = gaussFn     
        # x0 and I values to be replaced during loop
        guessTemp = [0, np.min(yData[domain]), 0,
                        xRange / 10]
        # x bounds to be replaced
        boundLowTemp = [xData[LDat], np.min(yData), 0., 0.]
        boundUppTemp = [xData[RDat], np.inf, np.inf, np.inf]

        boundLowerPart = [xData[0], np.min(yData), 0, 0]
        boundUpperPart = [xData[-1], np.inf, np.inf, np.inf]       


    # Initialize guess params and bounds for number of curves
    boundUpper = []
    boundLower = []
    guess = []
        
    # init fit array to compare residual 
    fit = np.ones_like(yData[domain]) * np.min(yData[domain])
    
    # parameter array: [x0, y0, Intensity, alpha, gamma]
    # set up guesses
    for i in range(numCurves): # for each curve
        # Calculate residual
        resid = fit - yData[domain]
        
        # place peak at min residual
        xPosGuess = xData[domain][np.argmin(resid)]
        guessTemp[0] = xPosGuess
        guessTemp[2] = np.max(resid) - np.min(resid)
            
        # Deal with edge cases.. 
        xPosLow = xPosGuess - 0.05*xRange
        if xPosLow < xData[0]: xPosLow = xData[0]
        xPosHigh = xPosGuess + 0.05*xRange
        if xPosHigh > xData[-1]: xPosHigh = xData[-1]
            
        # Update temp bounds 
        boundLowTemp[0] = xPosLow
        boundUppTemp[0] = xPosHigh

        boundTemp = tuple([boundLowTemp, boundUppTemp])
        try: 
            poptTemp, pcovTemp = curve_fit(func, xData[domain], -resid, 
                                        bounds=boundTemp, p0=guessTemp)  
            #print('Fit to residual at {0}'.format(xPosGuess))
        except RuntimeError as e:
            print(e) 
            poptTemp = guessTemp

        # build guess array, update fit
        guess += list(poptTemp)
        fit = func(xData[domain], *guess)

        # concatenate lists for bounds for real fit
        boundLower += boundLowerPart 
        boundUpper += boundUpperPart

        # Combine bounds into tuple for input
        bounds = tuple([boundLower, boundUpper])
    try:
    # Curve fit function call using guess and bounds
        popt, pcov = curve_fit(func, xData[domain], yData[domain], 
                                    bounds=bounds, p0=guess)
    except RuntimeError as e:
        print(e) 
        popt = np.array(guess)
       

    # Calculate FWHM for each peak fit
    if peakShape == 'Voigt':
        FWHM = []
        c0 = 2.0056
        c1 = 1.0593
        for i in range(0, len(popt), 5): # grab fwhm for each peak
            fg = 2*popt[i+3]*np.sqrt(2*np.log(2))
            fl = 2*popt[i+4]
            phi = fl / fg
            FWHM.append(fg * (1-c0*c1 + np.sqrt(phi**2 + 2*c1*phi + (c0*c1)**2)))

    elif peakShape == 'Gaussian':
        FWHM = []
        for i in range(0, len(popt), 4): # grab fwhm for each peak
            FWHM.append(2*popt[i+3]*np.sqrt(2*np.log(2)))
   
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


def calcFWHM(data, LDatum, RDatum):
    '''
    Calculate full width half maximum literally
    '''
    roots = []
    # Convert left and right bounds to integers for indexing
    LDat = int(LDatum)
    RDat = int(RDatum)
    domain = np.arange(LDat, RDat)
    
    # Separate x and y data series
    xData = data[0]
    yData = data[1]
  
    # Locate position of max in domain 
    maxInd = np.where(yData[domain] == np.max(yData[domain]))[0][0]
    # Shift according to LDat since data is complete 
    loc = int(LDat + maxInd)  # given all data
 
    #offset to 0 
    yDataOffset = yData - np.min(yData[domain]) 
    yRange = np.max(yDataOffset[domain])
    spline = UnivariateSpline(xData[domain], yDataOffset[domain] - yRange/2, s=0)

    roots = spline.roots()

    if len(roots) == 2: 
        return [xData[loc], (roots[1] - roots[0])]
    else:
        print('number of roots ({0}) != 2'.format(len(roots)))
        return [xData[loc], 'N/A']
