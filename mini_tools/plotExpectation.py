from ROOT import *
import pickle



#import matplotlib.pyplot as plt
#import numpy as np

lumi=35867.
xSec_tching = pickle.load(open("../SUSYxSections/xSec_SMS_C1C1_13TeV.pkl","rb"))
xSec_tching_neutr = pickle.load(open("../SUSYxSections/xSec_SMS_N2C1_13TeV.pkl","rb"))
xSec_t5zg = pickle.load(open("../SUSYxSections/xSec_SMS_Gluino_13TeV.pkl","rb"))

xSec_squark= pickle.load(open("../SUSYxSections/xSec_SMS_Squark_13TeV.pkl","rb"))

#brZG_tching=0.25
#brZG_t5zg=0.5
brZG_tching=1.
brZG_t5zg=1.
brZG_stealth=0.5

ratioSquark=1.0
#ratioSquark=0.1

#brZll=0.034
brZll=0.068

eE=0.8
eM=0.98
eP=0.9




weak = TH1F("hist","TChiNG",64,0.,1600.)
weakEffEE = TH1F("hist2","",64,0.,1600.)
weakEffMM = TH1F("hist3","",64,0.,1600.)

strong = TH1F("hist2","T5bbbbZG",500,0.,2500.)

squark = TH1F("hist4","Stealth Susy - T2 based",600,0.,3000.)

for key in xSec_tching:
    count= lumi * (xSec_tching[key][0]+xSec_tching_neutr[key][0])*brZG_tching*brZll
    weak.Fill(key,count)
    #weakEffEE.Fill(key,count*eE*eE*eP)
    #weakEffMM.Fill(key,count*eM*eM*eP)
    
    
for key in xSec_t5zg:
    count = lumi *xSec_t5zg[key][0]*brZG_t5zg*brZll
    strong.Fill(key,count)    

for key in xSec_squark:
    count = lumi *xSec_squark[key][0]*ratioSquark*brZll*brZG_stealth
    squark.Fill(key,count)    




gStyle.SetOptStat(0)
    
weak.GetXaxis().SetTitle("m_{NLSP}")
weak.GetYaxis().SetTitle("Expected Events")
    
c = TCanvas("canvas","",800,800)

text = TLatex( 600., 100000., "#splitline{#font[22]{Luminosity} #upoint #sigma(m_{#tilde{g}}) #upoint BF(#chi_{0}^{1}#rightarrow Z, #chi_{0}^{1}#rightarrow #gamma) #upoint BF(Z#rightarrow (ee,#mu#mu))}{= 35.867fb^{-1} #upoint #sigma #upoint 50% #upoint 6.8%}")
text.SetTextSize(0.025)

text2 = TLatex( 300., 500., "#splitline{#font[22]{Luminosity} #upoint #sigma(m_{NLSP}) #upoint BF(#chi_{0}^{1}#rightarrow Z, #chi_{0}^{1}#rightarrow #gamma) #upoint BF(Z#rightarrow (ee,#mu#mu))}{= 35.867fb^{-1} #upoint #sigma #upoint 25% #upoint 6.8%}")
text2.SetTextSize(0.025)

#text3 = TLatex( 600., 50000., "#splitline{#font[22]{Luminosity} #upoint #sigma(m_{#tilde{q}}) #upoint BF(#chi_{0}^{1}#rightarrow Z, #chi_{0}^{1}#rightarrow #gamma) #upoint BF(Z#rightarrow (ee,#mu#mu))}{= 35.867fb^{-1} #upoint #sigma #upoint 0.4 #upoint 50% #upoint 6.8%}")
text3 = TLatex( 600., 50000., "#splitline{#font[22]{Luminosity} #upoint #sigma(m_{#tilde{q}}) #upoint BF(#chi_{0}^{1}#rightarrow Z, #chi_{0}^{1}#rightarrow #gamma) #upoint BF(Z#rightarrow (ee,#mu#mu))}{= 35.867fb^{-1} #upoint #sigma #upoint 50% #upoint 6.8%}")
text3.SetTextSize(0.025)

line = TLine(0.,1.,1600.,1.)
line.SetLineColor(kBlack)
line.SetLineStyle(9)

jo = TLine(950.,0.,950.,200.)
jo.SetLineColor(kRed)


stealthBino = TLine(1100.,0.,1100.,1500.)
stealthBino.SetLineColor(kRed)
stealthWino = TLine(600.,0.,600.,1500.)
stealthWino.SetLineColor(kRed)

#knut = TLine()

#c.SetGrid()

weak.Draw("hist")
#weakEffEE.Draw("hist same")
#weakEffMM.Draw("hist same")
line.Draw()
jo.Draw()
c.SetLogy()
text2.Draw()
c.Update()
c.SaveAs('tching.pdf')

line = TLine(0.,1.,2500.,1.)
line.SetLineColor(kBlack)
line.SetLineStyle(9)


knut = TLine(1600.,0.,1600.,1500.)
knut2 = TLine(1900.,0.,1900.,1500.)
knut.SetLineColor(kRed)
knut2.SetLineColor(kRed)

strong.GetXaxis().SetTitle("m_{#tilde{g}}")
strong.GetYaxis().SetTitle("Expected Events")

strong.Draw("hist")
line.Draw()
knut.Draw()
knut2.Draw()
text.Draw()
#text2.Draw()
c.Update()
c.SaveAs("t5zg.pdf")

squark.GetXaxis().SetTitle("m_{#tilde{q}}")
squark.GetYaxis().SetTitle("Expected Events")

squark.Draw("hist")
#line.Draw()
#knut.Draw()
#knut2.Draw()
#text.Draw()
text3.Draw()
stealthBino.Draw()
stealthWino.Draw()
line.Draw()
c.Update()
c.SaveAs("squark.pdf")




