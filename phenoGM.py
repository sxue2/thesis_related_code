import csv
import numpy as np

def getPheno(phenoFolder,phenoFileName,trait, covList):
    '''
    Function to take in a phenotype file formatted as csv
    Assume it has one header line
    First column is a unique identifier of strings of the line IDs
    One column has a name that matches the trait 'trait' argument
    covList is a list of strings that match one or more column names to be used as the covariates
    Take the column that matches trait
    Return a list of three numpy arrays, first is an array of line names as strings, 
    second is an array of trait values, 
    third is an array of possibly multiple covariate values 
    '''
    ######################### pheno file
    with open(phenoFolder + phenoFileName,'r') as file:
        file_content = file.readlines()
    
    headerLine = file_content[0].rstrip() #Top line has useful header information. This is one big long string, it needs to be split
    headerList = headerLine.split('\t') #this splits it into a list
    print(trait)
    traitCol = [index for index,value in enumerate(headerList) if value == trait] #gets the column that matches the trait
    covCols = [index for index,value in enumerate(headerList) if value in covList] #gets the column that matches the covariate
    
    
    phenoData = np.genfromtxt(phenoFolder + phenoFileName, dtype = 'str', delimiter='\t', skip_header=1, missing_values='NA') 
    print(headerList,traitCol,covCols)
    lineNames = phenoData[:,0]
    traitValues = phenoData[:,traitCol]
    covariate = phenoData[:,covCols].astype(float)
    #covariate = [item for item in covariate]
  
    
    #deal with missing phenotypes
    #print(type(traitValues))
    #print(traitValues[1])
    
    #keepIndex = [item for item in range(len(traitValues)) if traitValues[item] != 'NA']

    return [lineNames, traitValues, covariate]
    
#testing code
def main():
    phenoFolder = 'X:/DATA/Public/shang/Ginnie-20170118T190718Z/Ginnie/data/'
    phenoFileName = 'Traits_S0S1new.txt'
    trait="DTA"
    DAratio=0
    unitTime=1
    phenoFolder = 'X:/DATA/Public/shang/Ginnie-20170118T190718Z/Ginnie/simulation/dataset/'
    phenoFileName = 'allPhenoRep_'+trait+"_DAratio_" + str(DAratio)+ "_unit_" + str(unitTime) +'.txt'
   
    trait = "rep1"
    covList = ['F.exp'] #a list of covariate names 
    phenoR = getPheno(phenoFolder, phenoFileName,trait, covList)
    pheno = phenoR[1]
    covariate = phenoR[2]
if __name__ =="__main__":
    main()