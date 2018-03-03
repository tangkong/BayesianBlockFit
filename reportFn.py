import csv
from os.path import basename
import numpy as np

def genOptParamCSV(savePath, file, optParams):
    '''
    Generate report CSV for given optParam dict and save path
    '''
    with open(savePath + basename(file)[:-7] + '_curveParams.csv', 'wb') as csv_file:
        writer = csv.writer(csv_file)
        if optParams['peakShape'] =='Voigt':
            #print('Array format: [x0, y0, I, alpha, gamma, FWHM]')
            writer.writerow([ 'peak', 'curve', 'x0', 'y0', 
                                'I', 'alpha', 'gamma', 'FWHM'])
        elif optParams['peakShape'] == 'Gaussian':
            #print('Array format: [x0, y0, I, sigma, FWHM]')
            writer.writerow(['peak', 'curve', 'x0', 'y0', 'I', 'sigma', 'FWHM'])
        for key, item in optParams.items():
            if key != 'numCurves' and key != 'peakShape':
                for i in range(len(item)):
                    #print('Param array for peak {0}, curve {1}: {2}'.format(key, i+1, 
                                            #np.array2string(item[i], precision=4)))

                    writer.writerow([key, i] + list(item[i]))
    
    print('OptParam report generated')

def genLitFWHMCSV(savePath, file, litFWHM):
    '''
    Generate report CSV for lit FWHM
    '''
    with open(savePath + basename(file)[:-7] + '_FWHM.csv', 'wb') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['peakNumber', 'peakLocation', 'FWHM'])

        for key, item in litFWHM.items():
            writer.writerow([key] + list(item))
    
    print('Literal FWHM report generated')
