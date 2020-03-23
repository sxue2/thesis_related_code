## 
def extractVC(VCdir,trait,chr):
    fileName = "asreml_"+trait+"_leave_"+"Chr"+str(chr)+".asr"
    f = open(VCdir+fileName,"r")
    i=0
    start="false"
    for line in f:
        #print line
        #if("Model_Term" in line):# this is for hpc asreml
        if("Source" in line and "Component" in line):    #this is for standalone asreml
            start="true"
               
        if(start == "true" and i>=0 and i<=2):
            print line
            if("giv(entry,1)" in line):
                Vg = line.split()[4]
            #if("Residual" in line):#this is for hpc asreml
            if("Variance" in line):##this is for standalone asreml
                Ve = line.split()[4]
        if(start=="true"):
            i=i+1
    VClist = [Vg,Ve]
    return VClist
if __name__ == "__main__":
    VCdir = 'X:/DATA/Public/shang/Ginnie-20170118T190718Z/Ginnie/Data_prep/leave_oneGmatrix_all_families.giv/asremlCode_estimate_VC/code_for_asreml/'
    
    print extractVC(VCdir,"KWt.20",1)
