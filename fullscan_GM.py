import os as os
import numpy as np
from numpy.linalg import inv
from datetime import datetime
from markerTestGM import markerTest
from phenoGM import getPheno
from __builtin__ import str
from subsetV import subsetV
#1. linear model GWAS without inbreeding coefficient
#2. linear model GWAS with inbreeding coefficient, but no other adjustment
#3. linear model GWAS with inbreeding coefficient plus G.eigen as fixed covariates 
#4. mixed model GWAS with inbreeding coefficient plus G from all 10 chromosomes
#5. mixed model GWAS with inbreeding coefficient plus G from 9 chromosomes not including current one.
 
def fullScan(model,trait):
    ### no magic number! here are constants
    nInd = 1846
    
    covFolder ='X:/DATA/Public/shang/Ginnie-20170118T190718Z/Ginnie/Data_prep/'
    phenoFolder = 'X:/DATA/Public/shang/Ginnie-20170118T190718Z/Ginnie/data/'
    phenoFileName = 'Traits_S0S1.txt'
    genoFolder='X:/DATA/Public/shang/Ginnie-20170118T190718Z/Ginnie/data/'
    outFolder= 'X:/DATA/Public/shang/Ginnie-20170118T190718Z/Ginnie/outputs/'
    VCdir = 'X:/DATA/Public/shang/Ginnie-20170118T190718Z/Ginnie/data/'
    
    covList = ['F.exp'] #a list of covariate names 
    phenoR = getPheno(phenoFolder, phenoFileName,trait, covList)
    pheno = phenoR[1]
    covariate = phenoR[2]
    #print(type(covariate))
    #print(covariate)
    
    ######################### genotype
    add = np.loadtxt(genoFolder+'ZeaSyn6_numeric_add_90trans.txt',skiprows=1,dtype=str) 
    dom = np.loadtxt(genoFolder+'ZeaSyn6_numeric_dom_90trans.txt',skiprows=1,dtype=str)  
    #add = np.loadtxt(genoFolder+'tryadd.txt',skiprows=1,dtype=str) 
    #dom = np.loadtxt(genoFolder+'tryadd.txt',skiprows=1,dtype=str)  
    
    totalMarker = add.shape[0] # shape[0] is marker number, shape[1] is individual number
    
    f_pvalue = open(outFolder + model +'pvalues.txt','w',1)
    f_beta = open(outFolder + model + 'beta.txt','w',1)
    
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
        #########################
        if(model == "model1" ):
            VeFixed = 500
            ######################### error variance matrix
            error = np.zeros((nInd,nInd),float)
            np.fill_diagonal(error, VeFixed)
            V =  error
            inverseV = subsetV(V,nonMissI)
            covFilter = []
            result = markerTest(addFilter,domFilter,covFilter, phenoFilter, nInd, V,nonMissI)
        if(model == "model2" ):
            V =  error
            result = markerTest(addFilter,domFilter,covFilter, phenoFilter, nInd, V,nonMissI)
        if(model == "model3" ):
            VeEigen = 500
            ######################### error variance matrix
            error = np.zeros((nInd,nInd),float)
            np.fill_diagonal(error, VeEigen)
            
            V =  error
            
            covList = ['F.exp','Geigen1','Geigen2','Geigen3',
                       'Geigen4','Geigen5','Geigen6','Geigen7',
                       'Geigen8','Geigen9','Geigen10'] #a list of covariate names 
            phenoR = getPheno(phenoFolder, phenoFileName,trait, covList)
            covariate = phenoR[2]
            result = markerTest(addFilter,domFilter,covFilter, phenoFilter, nInd, V,nonMissI)
        if(model == "model4" ):
            VarComps = [5492, 7053] #Vg and Verr for trait from asreml
            ########################## read in cov matrix  
            #read in the dense relationship matrix, it has no header or indicator columns
            #by default loadtxt makes a float matrix
            #current format of G matrix is square numeric with no row or column headings
            #sort order is enforced outside of this program to match row order of genotype scores
            K =  np.loadtxt('X:/DATA/Public/shang/Ginnie-20170118T190718Z/Ginnie/data/Gmatrix_all_families.txt.gz') 
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
            result = markerTest(addFilter,domFilter,covFilter, phenoFilter, nInd, V,nonMissI)
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
        
        count=count+1
        print("count",count)      
    f_pvalue.close()   

if __name__ == '__main__':
#     for i in range(5):
#         p = multiprocessing.Process(target=fullScan())
#         p.start()    
#         p.join()
    fullScan("model4","DTA")