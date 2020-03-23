import os as os
import numpy as np
from numpy.linalg import inv
from datetime import datetime
from markerTestGM import markerTest
from phenoGM import getPheno
from __builtin__ import str
from subsetV import subsetV
import sys
#1. linear model GWAS without inbreeding coefficient
#2. linear model GWAS with inbreeding coefficient, but no other adjustment
#3. linear model GWAS with inbreeding coefficient plus G.eigen as fixed covariates 
#4. mixed model GWAS with inbreeding coefficient plus G from all 10 chromosomes
#5. mixed model GWAS with inbreeding coefficient plus G from 9 chromosomes not including current one.
 
def fullScan(model,trait,DAratio,unitTime,rep):
    ### no magic number! here are constants
    nInd = 1846
    
    covFolder ='/home/sxue2/thirdProject/Data_prep/'
    phenoFolder = '/home/sxue2/thirdProject/simulation/data/'
    phenoFileName = 'allPhenoRep_'+trait+"_DAratio_" + DAratio+ "_unit_" + unitTime +'.txt'
    genoFolder='/home/sxue2/thirdProject/simulation/data/'
    outFolder= '/home/sxue2/thirdProject/simulation/outputs/'
    VCdir = '/home/sxue2/thirdProject/simulation/data/'
    
    covList = ['F.exp'] #a list of covariate names 
    repName = "rep"+str(rep)
    phenoR = getPheno(phenoFolder, phenoFileName,repName, covList)
    pheno = phenoR[1]
    covariate = phenoR[2]
    #print(type(covariate))
    #print(covariate)
    
    VarComps = [254.3738466, 1.851806e+01] #Vg and Verr for trait from asreml
    ########################## read in cov matrix  
    #read in the dense relationship matrix, it has no header or indicator columns
    #by default loadtxt makes a float matrix
    #current format of G matrix is square numeric with no row or column headings
    #sort order is enforced outside of this program to match row order of genotype scores
    K =  np.loadtxt('/home/sxue2/thirdProject/simulation/data/Gmatrix_all_families.txt.gz') 
    #type(K)   
    # read in corresponding parameters (from asreml output)
    
    Vg = float(VarComps[0])
    Ve = float(VarComps[1])
    ######################### error variance matrix
    error = np.zeros((nInd,nInd),float)
    np.fill_diagonal(error, Ve)
    ######################### Genetic covariance matrix
    # g-cov matrix
    ksigmaG = Vg*K 
    V = ksigmaG + error
    
    
    ######################### genotype
    add = np.loadtxt(genoFolder+'add_'+trait+"_DAratio_" + DAratio+ "_unit_" + unitTime + '_rep_'+rep+'.txt',skiprows=1,dtype=str) 
    dom = np.loadtxt(genoFolder+'dom_'+trait+"_DAratio_" + DAratio+ "_unit_" + unitTime + '_rep_'+rep+'.txt',skiprows=1,dtype=str)
    #add = np.loadtxt(genoFolder+'tryadd.txt',skiprows=1,dtype=str) 
    #dom = np.loadtxt(genoFolder+'tryadd.txt',skiprows=1,dtype=str)  
    
    
    f_pvalue = open(outFolder + model +trait+"_DAratio_" + DAratio+ "_unit_" + unitTime +'_rep_'+rep+'pvalues.txt','w',1)
    f_beta = open(outFolder + model + trait+"_DAratio_" + DAratio+ "_unit_" + unitTime +'_rep_'+rep+'beta.txt','w',1)
    f_dfMarker = open(outFolder + model + trait+"_DAratio_" + DAratio+ "_unit_" + unitTime +'_rep_'+rep+'dfMarker.txt','w',1)
    
    count = 1
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
        phenoFilter = [float(pheno[i][0]) for i in nonMissI]
        V1=[ (np.squeeze(np.asarray(V[i]))) for i in nonMissI]
        Vsub =np.asmatrix([V1[item][nonMissI] for item in range(len(V1))])
    
        print("Vshape",Vsub.shape)
        print("Do inverse OMG")
        inverseV = inv(Vsub)
        print("Inverse done")
        result = markerTest(addFilter,domFilter,covFilter, phenoFilter, nInd, inverseV,nonMissI)
        #if(model == "model5" ):
            #ksigmaG9 = 
            #V = ksigmaG9 + error
            #result = markerTest(addFilter,domFilter,covFilter, phenoFilter, nInd, V,nonMissI)
        p   = result[0]
        beta= result[1]
        print("p:",p)
        print("beta",beta)
        thisp =   [str(p)]
        if( type(beta) is int):
            thisBeta = [str(beta)]
        else:
            thisBeta =  beta.tolist()[0]
        pString = "\t".join(str(x) for x in thisp) 
        betaString = "\t".join(str(x) for x in thisBeta) 
        f_pvalue.write(pString+"\n")   
        f_beta.write(betaString+"\n")   
        dfMarker=result[2]
        f_dfMarker.write(str(dfMarker)+"\n")
        count=count+1
        print("count",count)      
    f_pvalue.close()   

if __name__ == '__main__':
#     for i in range(5):
#         p = multiprocessing.Process(target=fullScan())
#         p.start()    
#         p.join()
    fullScan("model4","EH",sys.argv[1],sys.argv[2],sys.argv[3])