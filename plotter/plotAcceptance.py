from ROOT import *
import pickle
import ROOT
from array import array
from include import *
import numpy as np


def getPointFromDir(name):
    m = re.match("(.*)_(.*)_(.*)", name)
    combi, m1, m2 = m.groups()
    m1,m2 = int(m1), int(m2)
    return combi, m1, m2



#Weak SMS TChiNG

#path="../minimal/output_signalScan/"
#path="../minimal/output/"
path="../minimal/output_noVeto/"

lumi=35867.
xSec_tching = pickle.load(open("../SUSYxSections/xSec_SMS_C1C1_13TeV.pkl","rb"))
xSec_tching_neutr = pickle.load(open("../SUSYxSections/xSec_SMS_N2C1_13TeV.pkl","rb"))
xSec_tching_total = pickle.load(open("data/xSec_SMS_TChiNG_13TeV.pkl","rb"))
#xSec_t5zg = pickle.load(open("../SUSYxSections/xSec_SMS_Gluino_13TeV.pkl","rb"))
xSec_t5zg = pickle.load(open("data/xSec_SMS_Gluino_13TeV.pkl","rb"))

brZG_tching=0.25
brZG_t5zg=0.5


weakSampleName = "SMS-TChiNG_BF50N50G_hists.root"
strongSampleName = "SMS-T5bbbbZg_hists.root"
GGM12SampleName = "GGM_GravitinoLSP_M1-200to1500_M2-200to1500_hists.root"
GGM13SampleName = "GGM_GravitinoLSP_M1-50to1500_M3-1000to2500_hists.root"




file_weak = ROOT.TFile(path+weakSampleName)
dirs = [k.GetName() for k in file_weak.GetListOfKeys() if k.GetName().startswith("Ng")]

weak = TGraph();

i=0
for key in dirs:
    folder = "sig/"
    hName = "met"
    point = getPointFromDir(key)
    #print point
    mNLSP = point[1]
    mDummy = point[2]
    histo = file_weak.Get(key+"/"+folder+hName)
    
    avgTopPtWeightHisto = file_weak.Get(key+"/"+folder+"weight_topPt")
    avgNIsrWeightHisto = file_weak.Get(key+"/"+folder+"weight_nISR")
    avgEWKinoWeightHisto = file_weak.Get(key+"/"+folder+"weight_EWKinoPairPt")
    avgleptonWeightHisto = file_weak.Get(key+"/"+folder+"weight_leptonPairPt")
    
    if avgTopPtWeightHisto.Integral()>0.:
        avgTopPtWeight = avgTopPtWeightHisto.GetMean()
    else:
        avgTopPtWeight=1.
    if avgNIsrWeightHisto.Integral()>0.:
        avgNIsrWeight = avgNIsrWeightHisto.GetMean()
    else:
        avgNIsrWeight=1.
    if avgEWKinoWeightHisto.Integral()>0.:
        avgEWKinoWeight = avgEWKinoWeightHisto.GetMean()
    else: avgEWKinoWeight=1.
    if avgleptonWeightHisto.Integral()>0.:
        avgleptonWeight = avgleptonWeightHisto.GetMean()
    else:
        avgleptonWeight=1.
#
    histo.Scale(1./avgTopPtWeight)
    histo.Scale(1./avgNIsrWeight)
    histo.Scale(1./avgEWKinoWeight)
    histo.Scale(1./avgleptonWeight)
    
    #print key+"/"+folder+hName
    #acc = histo.Integral(150.,-1)
    #acc = histo.Integral(histo.FindFixBin(150.),-1)*100.
    acc = histo.Integral()*100.
    
    weak.SetPoint(i,mNLSP,acc)
    i+=1




file_strong = ROOT.TFile(path+strongSampleName)
dirs = [k.GetName() for k in file_strong.GetListOfKeys() if k.GetName().startswith("Zg")]

strong = TGraph2D();

i=0
for key in dirs:
    folder = "sig/"
    hName = "met"
    point = getPointFromDir(key)
    mGluino = point[1]
    mNeutralino = point[2]
    histo = file_strong.Get(key+"/"+folder+hName)
    #acc = histo.Integral()*100.
    avgTopPtWeightHisto = file_strong.Get(key+"/"+folder+"weight_topPt")
    avgNIsrWeightHisto = file_strong.Get(key+"/"+folder+"weight_nISR")
    avgEWKinoWeightHisto = file_strong.Get(key+"/"+folder+"weight_EWKinoPairPt")
    avgleptonWeightHisto = file_strong.Get(key+"/"+folder+"weight_leptonPairPt")
    
    if avgTopPtWeightHisto.Integral()>0.:
        avgTopPtWeight = avgTopPtWeightHisto.GetMean()
    else:
        avgTopPtWeight=1.
    if avgNIsrWeightHisto.Integral()>0.:
        avgNIsrWeight = avgNIsrWeightHisto.GetMean()
    else:
        avgNIsrWeight=1.
    if avgEWKinoWeightHisto.Integral()>0.:
        avgEWKinoWeight = avgEWKinoWeightHisto.GetMean()
    else: avgEWKinoWeight=1.
    if avgleptonWeightHisto.Integral()>0.:
        avgleptonWeight = avgleptonWeightHisto.GetMean()
    else:
        avgleptonWeight=1.
#
    histo.Scale(1./avgTopPtWeight)
    histo.Scale(1./avgNIsrWeight)
    histo.Scale(1./avgEWKinoWeight)
    histo.Scale(1./avgleptonWeight)
    #acc = histo.Integral(histo.FindFixBin(150.),-1)*100.
    acc = histo.Integral()*100.
    strong.SetPoint(i,mGluino,mNeutralino,acc)
    i+=1
    

file_GGM12 = ROOT.TFile(path+GGM12SampleName)
dirs = [k.GetName() for k in file_GGM12.GetListOfKeys() if k.GetName().startswith("GGM")]
#print dirs

ggm1m2 = TGraph2D();

i=0
for key in dirs:
    folder = "sig/"
    hName = "met"
    point = getPointFromDir(key)
    #print point
    m1 = point[1]
    m2 = point[2]
    histo = file_GGM12.Get(key+"/"+folder+hName)
    #acc = histo.Integral()*100.
    acc = histo.Integral(histo.FindFixBin(150.),-1)*100.
    ggm1m2.SetPoint(i,m1,m2,acc)
    i+=1
    
file_GGM13 = ROOT.TFile(path+GGM13SampleName)
dirs = [k.GetName() for k in file_GGM13.GetListOfKeys() if k.GetName().startswith("GGM")]

ggm1m3 = TGraph2D();

i=0
for key in dirs:
    folder = "sig/"
    hName = "met"
    point = getPointFromDir(key)
    m1 = point[1]
    m3 = point[2]
    histo = file_GGM13.Get(key+"/"+folder+hName)
    #acc = histo.Integral()*100.
    acc = histo.Integral(histo.FindFixBin(150.),-1)*100.
    ggm1m3.SetPoint(i,m1,m3,acc)
    i+=1


import style
style.defaultStyle()
#style.defaultStyle()
#s= style.style2d()
#s.SetPadLeftMargin(0.18)
#style.setPaletteRWB()
    


weak.SetTitle("; m_{NLSP} (GeV); Acceptance x Efficiency [%]")

c = TCanvas("canvas","",800,800)

#l.Draw()

weak.SetMarkerStyle(20)
weak.SetMarkerSize(1)
#weak.Draw("AL*")
weak.Draw("A*")
#c.SetLogy()
#l = aux.Label(sim=True)
l = ROOT.TLatex( 0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
l2 = ROOT.TLatex( 0.21, .88, "#scale[0.76]{#font[52]{Simulation}}")
l.SetNDC()
l2.SetNDC()
l.Draw()
l2.Draw()
lum = ROOT.TLatex( .62, .95, "%.1f fb^{-1} (%s TeV)"%(aux.intLumi/1000., aux.Label.cmsEnergy) )
lum.SetNDC()
lum.Draw()
c.Update()
c.SaveAs('tching.pdf')


#style.defaultStyle()
style.style2d()
s= style.style2d()
s.SetPadLeftMargin(0.18)
c = TCanvas("canvas","",800,800)
strong.SetTitle("; m_{#tilde{g}} (GeV);m_{#tilde{#chi_{0}^{1}}} (GeV); Acceptance x Efficiency [%]")
ggm1m2.SetTitle("; M1;M2; Acceptance x Efficiency [%]")

lum = ROOT.TLatex( .62, .95, "#scale[0.76]{%.1f fb^{-1} (%s TeV)}"%(aux.intLumi/1000., aux.Label.cmsEnergy) )
lum.SetNDC()
strong.Draw("COLZ")
#l = aux.Label(sim=True)
l.Draw()
l2.Draw()
lum.Draw()
c.Update()
c.SaveAs("t5zg.pdf")

#c.Clear()
ggm1m2.Draw("COLZ")
#l = aux.Label(sim=True)
l.Draw()
lum.Draw()
c.Update()
c.SaveAs("ggm1m2.pdf")


#c.Clear()
ggm1m3.Draw("COLZ")
#l = aux.Label(sim=True)
l.Draw()
lum.Draw()
c.Update()
c.SaveAs("ggm1m3.pdf")
