## rawDN is the rawD number
## YN is the replication number
def extractVC(VCdir,rawDN,YN):
    fileName = "rawD"+str(rawDN)+"_Anl5Y"+str(YN)+".asr"
    f = open(VCdir+fileName,"r")
    i=0
    for line in f:
        #print line
        #if("Model_Term" in line):# this is for hpc asreml
        if("Source" in line and "Component" in line):    #this is for standalone asreml
            i=i+1
        if(i>=1 and i<=4):
            
            if("giv" in line):
                Vg = line.split()[4]
            #if("Residual" in line):#this is for hpc asreml
            if("Variance" in line):##this is for standalone asreml
                Ve = line.split()[4]
            if("env" in line):
                Venv = line.split()[4]
            i=i+1
    VClist = [Vg,Venv,Ve]
    return VClist
