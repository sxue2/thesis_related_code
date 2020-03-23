######################### important! dealing with NA
# because a and d columns has NA value, here will 
# make the matrix has the right dimension
# need to fix V
# nonMissI is the index for non missing 
import numpy as np
from numpy.linalg import inv

def subsetV(V,nonMissI):
    V1=[ (np.squeeze(np.asarray(V[i]))) for i in nonMissI]
    Vsub =np.asmatrix([V1[item][nonMissI] for item in range(len(V1))])
    
    print("Vshape",Vsub.shape)
    print("Do inverse OMG")
    inverseV = inv(Vsub)
    print("Inverse done")
    return(inverseV)