import numpy as np
from numpy.linalg import inv
from scipy import stats
from _ast import Add
from subsetV import subsetV
# this code is translated from testMarkerUsingP3D, CompressedMLMusingDoubleMatrix in Tassel source code
def markerTest(add,dom,covariate,Y,nInd,inverseV ,nonMissI):
    x0 = [1 for i in range(len(Y))]  #this is for the intercept
    x1 = add
    x2 = dom
    dfMarker = 0
    #xMatrix = np.column_stack((x0,x1,x2))
    if(not np.array_equal(add,dom)) :
        dfMarker = 2
        xMatrix = np.column_stack((x0,x1,x2,covariate))
        try:
            print("done:numericX")
            xVx = np.dot(np.dot(xMatrix.transpose(),inverseV),xMatrix)
            invxVx = inv(xVx)
            xtV = np.dot(xMatrix.transpose(),inverseV)
            invxVxxV = np.dot(invxVx,xtV )
            beta  = np.dot(invxVxxV,Y)
            print("done: beta")
            print(beta)
            nParm = beta.shape[1]
            firstMarker = 1
            M = np.zeros((dfMarker,nParm),float)
              
            for i in range(dfMarker): 
                M[i, (i + firstMarker)] = 1;
            print(M)
            MB = np.dot(M,beta.transpose());
            invMiM = np.dot(M,np.dot(invxVx,M.transpose()))
            invMiM = inv(invMiM)
            preF = np.dot(MB.transpose(),np.dot(invMiM,MB))
            F= preF[0,0] # get the left top element
            F /= dfMarker
            print(F,dfMarker,nInd-nParm)
            
            p = 1-stats.f.cdf(F, dfMarker, nInd - nParm)
        except:
            p=0
            beta=0
    else:
        dfMarker = 1
        xMatrix = np.column_stack((x0,x1,covariate))
        try:
            print("done:numericX")
            xVx = np.dot(np.dot(xMatrix.transpose(),inverseV),xMatrix)
            invxVx = inv(xVx)
            xtV = np.dot(xMatrix.transpose(),inverseV)
            invxVxxV = np.dot(invxVx,xtV )
            beta  = np.dot(invxVxxV,Y)
            print("done: beta")
            print(beta)
            nParm = beta.shape[1]
            firstMarker = 1
            M = np.zeros((dfMarker,nParm),float)
              
            for i in range(dfMarker): 
                M[i, (i + firstMarker)] = 1;
            print(M)
            MB = np.dot(M,beta.transpose());
            invMiM = np.dot(M,np.dot(invxVx,M.transpose()))
            invMiM = inv(invMiM)
            preF = np.dot(MB.transpose(),np.dot(invMiM,MB))
            F= preF[0,0] # get the left top element
            F /= dfMarker
            print(F,dfMarker,nInd-nParm)
            
            p = 1-stats.f.cdf(F, dfMarker, nInd - nParm)
        except:
            p=0
            beta=0
        
    return [p,beta,dfMarker]