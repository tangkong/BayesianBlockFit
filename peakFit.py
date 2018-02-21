def peakFit(data, LDatum, RDatum, peakShape, numCurves):
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
    ''' 
    from scipy.optimize import curve_fit
    from scipy import integrate, signal
    from scipy.special import wofz
    import numpy as np
    import matplotlib.pyplot as plt
    
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
    def voigtFn(x, *params):
        """
        Return Voigt line shape centered at x0 with intensity I
        alpha: Lorentzian comp HWHM
        gamma: Gaussian comp HWHM
        I: Intensity
        x0, y0: Offsets
        [x0, y0, I, alpha, gamma] * numCurves
        """
        result = 0
        for i in range(0,len(params),5):
            x0    = params[i]
            y0    = params[i+1]
            I     = params[i+2]
            alpha = params[i+3]
            gamma = params[i+4]
        
            sigma = alpha / np.sqrt(2 * np.log(2))

            result = result + I * np.real(wofz(((x-x0) + 1j*gamma)/sigma/np.sqrt(2))) / sigma / np.sqrt(2*np.pi) + y0
        
        return result


    def gaussFn(x, *params):
        """
        Return Gaussianb line shape centered at x0 with intensity I
        sigma: Gaussian comp stdev
        I: Intensity
        x0, y0: Offsets
        x0, y0, I, sigma
        """
        result = 0
        for i in range(0,len(params),4):
            x0    = params[i]
            y0    = params[i+1]
            I     = params[i+2]
            sigma = params[i+3]
            
            result += I * np.exp(-(x-x0)**2 / (2*sigma**2)) + y0
        
        return result
    #############################################################################
    ### End peak functions
    ### Begin fitting using desired function
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
            
            guess += [xPosGuess, np.mean(yData[domain]), 1, 0.05, 0.05]
            
            # concatenate lists for bounds
            boundLower += [xData[LDat], np.min(yData[domain]), 0, 0, 0]
            boundUpper += [xData[RDat], np.max(yData[domain]), np.inf, np.inf, np.inf]       
            
        # Combine bounds into tuple for input
        bounds = tuple([boundLower, boundUpper])
        
        # Curve fit function call using guess and bounds
        popt, pcov = curve_fit(voigtFn, xData[domain], yData[domain], 
                       bounds=bounds, p0=guess)
    elif peakshape == 'Gaussian':############################################
        func = gaussFn
        # Initialize guess params and bounds for number of curves
        boundUpper = []
        boundLower = []
        guess = []
        
        # parameter array: [x0, y0, Intensity, sigma]
        for i in range(numCurves): # for each curve
            # Space peak locations evenly
            xPosGuess = (xData[LDat] + xRange * (i+1) / (numCurves+1))
            
            guess += [xPosGuess, np.mean(yData[domain]), 1, 0.05]
            
            # concatenate lists for bounds
            boundLower += [xData[LDat], np.min(yData[domain]), 0, 0]
            boundUpper += [xData[RDat], np.max(yData[domain]), np.inf, np.inf]       
            
        # Combine bounds into tuple for input
        bounds = tuple([boundLower, boundUpper])
        popt, pcov = curve_fit(gaussFn, xData[domain], yData[domain], 
                       bounds=([0.9*xData[loc] ,np.min(yData),0,0],
                                [1.1*xData[loc],np.max(yData),np.inf,np.inf]))
   

    # Plot setup
    #plt.figure(1)
    #plt.plot(xData[domain], yData[domain], label='data')
    #plt.plot(xData[domain], func(xData[domain], *popt), label='combined data')
    
    # Organize final parameters into array
    finalParams = []
    for j in range(numCurves):   # Plot each individual curve
        L = 0 + j * len(popt) / numCurves  # Sep popt array
        R = (j+1) * len(popt) / numCurves
        #plt.plot(xData[domain], func(xData[domain], *guess[L:R]), 
                 #'--', label='guessed curve: ' + str(j))
        #plt.plot(xData[domain], func(xData[domain], *popt[L:R]), 
                 #'.', label='optimized curve: ' + str(j))
        
        finalParams.append(popt[L:R])
        
    #plt.plot(xData[loc], yData[loc], marker='o') # Max position
    
    return finalParams
