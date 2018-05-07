#from dataMC import labels,frange

import ROOT
from ROOT import *
from array import array
from include import *
import numpy as np
import pickle as pkl
#import CR_tt,CR_DY,CR_WZ,CR_ZZ

import style
style.defaultStyle()

pklZZ = pkl.load( open( "plots_CR_zz/factors/CRZZ.pkl", "rb" ) )
pklDY = pkl.load( open( "plots_CR_dy/factors/CRDY.pkl", "rb" ) )
pklTT = pkl.load( open( "plots_CR_tt/factors/CRTT.pkl", "rb" ) )
pklWZ = pkl.load( open( "plots_CR_wz/factors/CRWZ.pkl", "rb" ) )
pklZZchi = pkl.load( open( "plots_CR/chi/ZZ_chi.pkl", "rb" ) )
#pklDYEEchi = pkl.load( open( "plots_CR/chi/DY_chi_EE.pkl", "rb" ) )
#pklDYMMchi = pkl.load( open( "plots_CR/chi/DY_chi_MM.pkl", "rb" ) )
#pklDYLLchi = pkl.load( open( "plots_CR/chi/DY_chi_LL.pkl", "rb" ) )
pklTTchi = pkl.load( open( "plots_CR/chi/TT_chi.pkl", "rb" ) )
pklWZchi = pkl.load( open( "plots_CR/chi/WZ_chi.pkl", "rb" ) )

valueZZInt=pklZZ["LL"]["eta1"][0]
errZZInt=pklZZ["LL"]["eta1"][1]
arZZChi=[]
arZZChiUp=[]
arZZChiDn=[]
arZZChiName=[]

for comb in pklZZchi:
    #print comb
    for variable in pklZZchi[comb]:
        #print variable,pklZZchi[comb][variable]["value"]
        arZZChi.append(pklZZchi[comb][variable]["value"])
        arZZChiUp.append(pklZZchi[comb][variable]["erUp"])
        arZZChiDn.append(pklZZchi[comb][variable]["erDown"])
        arZZChiName.append(variable)

#print valueZZInt, errZZInt
print arZZChi
print arZZChiUp
print arZZChiDn
print arZZChiName

c = TCanvas("canvas","",1200,800)


#ge = TGraphErrors()
#ge.SetPoint(0,valueZZInt,4.5)
#ge.SetPointError(0,errZZInt,4.5)
#ge.SetFillColor(ROOT.kBlue)
#ge.SetFillStyle(1001)
#ge.Draw("2")
#c.Update()


gr = TGraphAsymmErrors()

for i in range(len(arZZChi)):
    y= float(i+1)
    gr.SetPoint(i,arZZChi[i],y)
    gr.SetPointError(i,arZZChiDn[i],arZZChiUp[i],0.,0.)

gr.SetTitle("; Scale Factor #alpha; ")
gr.SetMarkerColor(4)
gr.SetMarkerStyle(21)
gr.Draw("same AP")
#c.Update()

line=TLine(valueZZInt,0,valueZZInt,11)
line.SetLineColor(kBlue)
line.SetLineWidth(2)
line.Draw("SAME")
lineUp=TLine(valueZZInt+errZZInt,0,valueZZInt+errZZInt,11)
lineUp.SetLineColor(kBlue)
lineUp.SetLineWidth(1)
lineUp.Draw("SAME")
lineDn=TLine(valueZZInt-errZZInt,0,valueZZInt-errZZInt,11)
lineDn.SetLineColor(kBlue)
lineDn.SetLineWidth(1)
lineDn.Draw("SAME")






c.SaveAs('plots_CR/chi/ZZ_Compare.pdf')





