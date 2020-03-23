import numpy as np
from numpy.linalg import inv
from scipy import stats

def fixedMarkerTest(add,dom,covariate,Y,nInd):
    dfMarker=0
    try:
        x0 = [1 for i in range(len(Y))]  #this is for the intercept
        x1 = add
        x2 = dom
        if(not np.array_equal(add,dom)) :
            dfMarker=2
            xFull = np.column_stack((x0,x1,x2,covariate))
            xReduced = np.column_stack((x0,covariate))
            dfFull = xFull.shape[1]-1
            dfReduced = xReduced.shape[1]-1
            betaFull = np.dot(np.dot(inv(np.dot(xFull.transpose(),xFull)),xFull.transpose()),Y)
            betaReduced = np.dot(np.dot(inv(np.dot(xReduced.transpose(),xReduced)),xReduced.transpose()),Y)

            predFull = np.dot(xFull,betaFull)
            predReduced = np.dot(xReduced,betaReduced)

            SSEFull = sumOfSquare(Y-predFull)
            SSEReduced = sumOfSquare(Y-predReduced)
            errorDfFull = len(Y) - dfFull -1
            F = ((SSEReduced - SSEFull) / (dfFull-dfReduced))/(SSEFull/errorDfFull)

            p = 1-stats.f.cdf(F, dfFull-dfReduced, errorDfFull)
            print("F",F,"p",p)
        else:
            dfMarker=1
            xFull = np.column_stack((x0,x1,covariate))
            xReduced = np.column_stack((x0,covariate))
            dfFull = xFull.shape[1]-1
            dfReduced = xReduced.shape[1]-1
            betaFull = np.dot(np.dot(inv(np.dot(xFull.transpose(),xFull)),xFull.transpose()),Y)
            betaReduced = np.dot(np.dot(inv(np.dot(xReduced.transpose(),xReduced)),xReduced.transpose()),Y)
            
            predFull = np.dot(xFull,betaFull)
            predReduced = np.dot(xReduced,betaReduced)
            
            SSEFull = sumOfSquare(Y-predFull)
            SSEReduced = sumOfSquare(Y-predReduced)
            errorDfFull = len(Y) - dfFull -1 
            F = ((SSEReduced - SSEFull) / (dfFull-dfReduced))/(SSEFull/errorDfFull)
            
            p = 1-stats.f.cdf(F, dfFull-dfReduced, errorDfFull)
            print("F",F,"p",p)
    except:
        p=0
        betaFull = 0
        
    return(p,betaFull,dfMarker)
def fixedMarkerTest0(add,dom,Y,nInd):
    try:
        x0 = [1 for i in range(len(Y))]  #this is for the intercept
        x1 = add
        x2 = dom
        dfMarker=0
        if(not np.array_equal(add,dom)) :
            dfMarker=2
            xFull = np.column_stack((x0,x1,x2))
            dfFull = xFull.shape[1]-1 
            dfReduced = 0
            xTx = np.dot(xFull.transpose(),xFull)
            invxtx = inv(xTx)
            invxtxxt = np.dot(invxtx,xFull.transpose())
            betaFull = np.dot(invxtxxt,Y)
            predFull = np.dot(xFull,betaFull)
            
            predReduced = np.mean(Y)
            SSEFull = sumOfSquare(Y-predFull)
            SSEReduced = sumOfSquare(Y-predReduced)
            errorDfFull = len(Y) - dfFull -1 
            F = ((SSEReduced - SSEFull) / (dfFull-dfReduced))/(SSEFull/errorDfFull)
            
            p = 1-stats.f.cdf(F, (dfFull-dfReduced), errorDfFull)
            print("F",F,"p",p)  
        else:
            dfMarker=1
            xFull = np.column_stack((x0,x1))
            dfFull = xFull.shape[1]-1 
            dfReduced = 0
            xTx = np.dot(xFull.transpose(),xFull)
            invxtx = inv(xTx)
            invxtxxt = np.dot(invxtx,xFull.transpose())
            betaFull = np.dot(invxtxxt,Y)
            predFull = np.dot(xFull,betaFull)
            
            predReduced = np.mean(Y)
            SSEFull = sumOfSquare(Y-predFull)
            SSEReduced = sumOfSquare(Y-predReduced)
            errorDfFull = len(Y) - dfFull -1
            F = ((SSEReduced - SSEFull) / (dfFull-dfReduced))/(SSEFull/errorDfFull)
            
            p = 1-stats.f.cdf(F, (dfFull-dfReduced), errorDfFull)
            print("F",F,"p",p)   
    except:
        p=0
        betaFull = 0
        
    return(p,betaFull,dfMarker)     
def sumOfSquare(yourArray):
    SS = 0
    for i in yourArray:
        SS = SS + i * i
    return SS