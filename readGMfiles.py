import gzip
import numpy as np

#get the numeric genotype scores for additive effects from ZeaSyn Full Sib Families

with gzip.open('C:/Users/jholland/Google Drive/Ginnie/data/ZeaSyn6_S0_numeric_add.txt.gz', 'rb') as f:
    file_content = f.readlines()

headerLine = file_content[0] #this is one big long string, it needs to be split
headerList = headerLine.split() #this splits it into a list
print(headerList[:11])

#read it as character, because NA causes fail when read as numeric. Then we convert Na to NaN and convert all to numeric
#skip first column [0] with usecols=range(1,...)
#get complete data set with usecols = range(1,(len(headerList)+1))
genoA =  np.loadtxt('C:/Users/jholland/Google Drive/Ginnie/data/ZeaSyn6_S0_numeric_add.txt.gz', dtype= 'str', comments='#', delimiter=None, skiprows=1, usecols=range(1,11)) 
genoA[genoA == "NA"] = "NaN"
genoA = genoA.astype(np.float)
print(genoA[:10, :6])
print(genoA.dtype)

#get the numeric genotype scores for dominance effects from ZeaSyn Full Sib Families

with gzip.open('C:/Users/jholland/Google Drive/Ginnie/data/ZeaSyn6_S0_numeric_dom.txt.gz', 'rb') as f:
    file_contentD = f.readlines()

headerLineD = file_contentD[0] #this is one big long string, it needs to be split
headerListD = headerLineD.split() #this splits it into a list
print(headerListD[:11])

genoD =  np.loadtxt('C:/Users/jholland/Google Drive/Ginnie/data/ZeaSyn6_S0_numeric_dom.txt.gz', dtype= 'str', comments='#', delimiter=None, skiprows=1, usecols=range(1,11)) 
genoD[genoD == "NA"] = "NaN"
genoD = genoD.astype(np.float)
print(genoD[:10, :6])
print(genoD.dtype)

