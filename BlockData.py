import numpy as np
######Data wrapper object for BBA.  Contains parameters and data
### Parameters
### Data: np array.  depVar = data[0], indVar = data[1]
#######
class BlockData:
    
    def __init__(self, dataArray, fit_order, ncp_prior):
        self.data = np.array(dataArray)
        self.fit_order = fit_order
        self.ncp_prior = ncp_prior
        
    def bkgdSub(self, **kwarg):
        if 'fitOrder' in kwarg:
            self.fit_order = fitOrder
            
        if self.fit_order is 2:
            self.fitCoeff = np.polyfit(self.data[0], self.data[1], self.fit_order)
            fitDataX = self.data[0]
            fitDataY = self.data[1] - np.polyval(self.fitCoeff, self.data[0])
            self.subData = np.array([fitDataX, fitDataY])
            
            print 'Quadratic background subtraction completed'
            
    def trimSubData(self):
        self.subData = np.array([self.subData[0,11:-50], self.subData[1,11:-50]])

    def blockFinder(self):
    #### Bayesian Block finding algorithm
    ## Input: BlockData structure 
    ##      (Seems to be using point measurement type, code for that alone at first)
    ## Output: BlockData structure with new field denoting blocks and bumps.
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
                sumX1 = np.cumsum(self.subData[1,r::-1]) # y-data background subtracted
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

        #print('cptUse start: {0}, end: {1}, len: {2}'.format(cptUse[0],cptUse[-1], len(cptUse)))
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