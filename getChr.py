import numpy as np

def getChr(index):
    folder= "X:/DATA/Public/shang/Ginnie-20170118T190718Z/Ginnie/data/"
    fileName = "Sites_in_filtered_num_coeff.txt"
    markerName = np.genfromtxt(folder + fileName, dtype = 'str', delimiter=',', missing_values='NA') 
    thisMarker = markerName[index]
    Schr = thisMarker.split('_')[0]
    print(Schr[1:])
if __name__ == "__main__":
    getChr(0)
    getChr(1)
    getChr(2)
    getChr(3)
    getChr(4)
    getChr(10000)
    getChr(30184)# last index