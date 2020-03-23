import numpy as np
from phenoGM import getPheno
from numpy.linalg import inv

def getSEofBeta(trait,Vg,Ve):
    covFolder ='/home/sxue2/thirdProject/Data_prep/'
    phenoFolder = '/home/sxue2/thirdProject/data/'
    phenoFileName = 'Traits_S0S1.txt'
    genoFolder='/home/sxue2/thirdProject/data/'
    outFolder= '/home/sxue2/thirdProject/outputs/'
    VCdir = '/home/sxue2/thirdProject/data/'
    
    covList = ['F.exp'] #a list of covariate names 
    phenoR = getPheno(phenoFolder, phenoFileName,trait, covList)
    covariate = phenoR[2]
    
    K =  np.loadtxt('/home/sxue2/thirdProject/data/Gmatrix_all_families.txt.gz') 
    #add = np.loadtxt(genoFolder+'ZeaSyn6_numeric_add_90trans.txt',skiprows=1,dtype=str) 
    #dom = np.loadtxt(genoFolder+'ZeaSyn6_numeric_dom_90trans.txt',skiprows=1,dtype=str)  
    add = np.loadtxt(genoFolder+'tryadd.txt',skiprows=1,dtype=str) 
    dom = np.loadtxt(genoFolder+'trydom.txt',skiprows=1,dtype=str)  
    varFile = open(outFolder +trait+'errorVar.txt','w',1)
    for index in range(add.shape[0]):
        thisAdd = add[index,:]
        thisDom = dom[index,:]
        
        keepAddI = [item for item in range(len(thisAdd)) if thisAdd[item] != 'NA']
        keepDomI = [item for item in range(len(thisDom)) if thisDom[item] != 'NA']
        nonMissI = list(set(keepAddI).intersection(keepDomI))
        # the following two lines confirmed that if add had NA then dom had NA too
        #if(set(keepAddI) != set(keepDomI)):
        #print("this one doesn't match",count)
        
        addFilter = [float(thisAdd[i]) for i in nonMissI]  
        domFilter = [float(thisDom[i]) for i in nonMissI]
        covFilter = [covariate[i] for i in nonMissI]
        indN = len(addFilter)
        x0 = [1 for i in range(indN)]  #this is for the intercept
        x1 = addFilter
        x2 = domFilter
        xMatrix = np.column_stack((x0,x1,x2,covFilter))
        
        xTx = np.dot(xMatrix.transpose(),xMatrix)
        upper = np.column_stack((xTx,xMatrix.transpose()))
        
        identity = np.zeros((indN,indN),float)
        np.fill_diagonal(identity, 1)
        lower = np.column_stack((xMatrix,identity+(Ve/Vg)*K))
        
        
        
        try:
            all = np.row_stack((upper,lower))
            allInv = inv(all)
            
            varFile.write(str(allInv[1,1]*Ve)+"\t"+str(allInv[2,2]*Ve)+"\n")
        except:
            varFile.write("NA" + "\t" + "NA" +"\n")
if __name__ =="__main__":
    trait="DTA"
    Vg = 3
    Ve= 5
    getSEofBeta(trait,Vg,Ve)    
        