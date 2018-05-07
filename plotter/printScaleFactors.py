#from dataMC import labels,frange
#
#import ROOT
#from ROOT import *
#from array import array
#from include import *
import numpy as np
import pickle as pkl

import matplotlib.pyplot as plt

#pklZZ = pkl.load( open( "plots_CR/factors/ControlRegionZZ.pkl", "rb" ))
#pklDY = pkl.load( open( "plots_CR/factors/ControlRegionDY.pkl", "rb" ))
#pklTT = pkl.load( open( "plots_CR/factors/ControlRegionTT.pkl", "rb" ))
#pklWZ = pkl.load( open( "plots_CR/factors/ControlRegionWZ.pkl", "rb" ))
pklZZ = pkl.load( open( "plots_CR_zz/factors/CRZZ.pkl", "rb" ))
pklDY = pkl.load( open( "plots_CR_dy/factors/CRDY.pkl", "rb" ))
pklTT = pkl.load( open( "plots_CR_tt/factors/CRTT.pkl", "rb" ))
pklWZ = pkl.load( open( "plots_CR_wz/factors/CRWZ.pkl", "rb" ))

zz=[]
zzErr=[]
dy=[]
dyEE=[]
dyEEErr=[]
dyMM=[]
dyMMErr=[]
dyLL=[]
dyLLErr=[]
tt=[]
ttErr=[]
wz=[]
wzErr=[]

print ""
print "------- ZZ --------"
#for key in pklZZ["LL"]:
for variable in pklZZ["LL"]:
    print variable, pklZZ["LL"][variable]
    zz.append(pklZZ["LL"][variable][0])
    zzErr.append(pklZZ["LL"][variable][1])
print ""
print "------- TT --------"
for key in pklTT:
    for variable in pklTT[key]:
        print variable, pklTT[key][variable]
        tt.append(pklTT[key][variable][0])
        ttErr.append(pklTT[key][variable][1])
print ""
print "------- DY --------"
#for key in pklDY:
    #for variable in pklDY[key]:
        #print variable, pklDY[key][variable]
        #dy.append(pklDY[key][variable][0])
#for key in pklDY:
for variable in pklDY["EE"]:
    print variable, pklDY["EE"][variable]
    dyEE.append(pklDY["EE"][variable][0])
    dyEEErr.append(pklDY["EE"][variable][1])
for variable in pklDY["MM"]:
    print variable, pklDY["MM"][variable]
    dyMM.append(pklDY["MM"][variable][0])
    dyMMErr.append(pklDY["MM"][variable][1])
for variable in pklDY["LL"]:
    print variable, pklDY["LL"][variable]
    dyLL.append(pklDY["LL"][variable][0])
    dyLLErr.append(pklDY["LL"][variable][1])
print ""
print "------- WZ --------"
#for key in pklWZ["LL"]:
for variable in pklWZ["LL"]:
    print variable, pklWZ["LL"][variable]
    wz.append(pklWZ["LL"][variable][0])
    wzErr.append(pklWZ["LL"][variable][1])


zz=np.array(zz)
dy=np.array(dy)
tt=np.array(tt)
wz=np.array(wz)
print ""
print "ZZ:"
print "mean: ", np.mean(zz), " (+-) ", np.std(zz), " +- ", np.mean(zzErr), " ", np.mean(zzErr)/np.mean(zz)*100.,"%"
print "TT:"
print "mean: ", np.mean(tt), " (+-) ", np.std(tt), " +- ", np.mean(ttErr), " ", np.mean(ttErr)/np.mean(tt)*100.,"%"
print "DYEE:"
print "mean: ", np.mean(dyEE), " (+-) ", np.std(dyEE), " +- ", np.mean(dyEEErr), " ", np.mean(dyEEErr)/np.mean(dyEE)*100.,"%"
print "DYMM:"
print "mean: ", np.mean(dyMM), " (+-) ", np.std(dyMM), " +- ", np.mean(dyMMErr), " ", np.mean(dyMMErr)/np.mean(dyMM)*100.,"%"
print "DYLL:"
print "mean: ", np.mean(dyLL), " (+-) ", np.std(dyLL), " +- ", np.mean(dyLLErr), " ", np.mean(dyLLErr)/np.mean(dyLL)*100.,"%"
print "WZ:"
print "mean: ", np.mean(wz), " (+-) ", np.std(wz), " +- ", np.mean(wzErr), " ", np.mean(wzErr)/np.mean(wz)*100.,"%"
