from ROOT import *
import pickle
import ROOT
from array import array
from include import *
import numpy as np


def getPointFromDir(name):
    m = re.match("(.*)_(.*)_(.*)", name)
    combi, m1, m2 = m.groups()
    m1, m2 = int(m1), int(m2)
    return combi, m1, m2


# Weak SMS TChiNG

# path="../minimal/output_signalScan/"
# path="../minimal/output/"
# path="../myAnalyzer/output_noTopPt_noNIsr/"
# path="../myAnalyzer/output/"
# path="../myAnalyzer/output_2/"
path = "../myAnalyzer/output_AN/"
# path="../minimal/output_noVeto/"

lumi = 35867.
xSec_tching = pickle.load(
    open("../SUSYxSections/xSec_SMS_C1C1_13TeV.pkl", "rb"))
xSec_tching_neutr = pickle.load(
    open("../SUSYxSections/xSec_SMS_N2C1_13TeV.pkl", "rb"))
xSec_tching_total = pickle.load(open("data/xSec_SMS_TChiNG_13TeV.pkl", "rb"))
#xSec_t5zg = pickle.load(open("../SUSYxSections/xSec_SMS_Gluino_13TeV.pkl","rb"))
xSec_t5zg = pickle.load(open("data/xSec_SMS_Gluino_13TeV.pkl", "rb"))

brZG_tching = 0.25
brZG_t5zg = 0.5


weakSampleName = "SMS-TChiNG_BF50N50G_hists.root"
strongSampleName = "SMS-T5bbbbZg_hists.root"
GGM12SampleName = "GGM_GravitinoLSP_M1-200to1500_M2-200to1500_hists.root"
GGM13SampleName = "GGM_GravitinoLSP_M1-50to1500_M3-1000to2500_hists.root"

gmsbSampleName = "GMSB_GravitinoLSP_N1decays_hists.root"


file_weak = ROOT.TFile(path + weakSampleName)
dirs = [k.GetName() for k in file_weak.GetListOfKeys()
        if k.GetName().startswith("Ng")]

#weak = TGraph();
weak = TGraphErrors()

i = 0
for key in dirs:
    folder = "sig/LL/nom/"
    hName = "met"
    point = getPointFromDir(key)
    #print point
    mNLSP = point[1]
    mDummy = point[2]
    histo = file_weak.Get(key + "/" + folder + hName)
    topPtWeight = aux.getWeightForWeights(
        path + weakSampleName, histoName="weightHisto_" + str(mNLSP), whichWeight="pu_mc_toppt")
    nIsrWeight = aux.getWeightForWeights(
        path + weakSampleName, histoName="weightHisto_" + str(mNLSP), whichWeight="pu_mc_nisr")
    ewkWeight = aux.getWeightForWeights(
        path + weakSampleName, histoName="weightHisto_" + str(mNLSP), whichWeight="pu_mc_ewk")
    #Weight = aux.getWeightForWeights(path+weakSampleName,histoName="weightHisto_"+str(mNLSP),whichWeight="pu_mc_toppt")

    histo.Scale(topPtWeight * nIsrWeight * ewkWeight)

    #print key+"/"+folder+hName
    #acc = histo.Integral(150.,-1)
    #acc = histo.Integral(histo.FindFixBin(150.),-1)*100.
    acc = histo.Integral() * 100.
    entries = histo.GetEntries()
    errNrel = np.sqrt(entries) / entries
    print errNrel
    weak.SetPoint(i, mNLSP, acc)
    weak.SetPointError(i, 0., acc * errNrel)
    i += 1


file_strong = ROOT.TFile(path + strongSampleName)
dirs = [k.GetName() for k in file_strong.GetListOfKeys()
        if k.GetName().startswith("Zg")]

strong = TGraph2D()

i = 0
for key in dirs:
    folder = "sig/LL/nom/"
    hName = "met"
    point = getPointFromDir(key)
    mGluino = point[1]
    mNeutralino = point[2]
    histo = file_strong.Get(key + "/" + folder + hName)
    topPtWeight = aux.getWeightForWeights(
        path + weakSampleName, histoName="weightHisto_" + str(mNLSP), whichWeight="pu_mc_toppt")
    nIsrWeight = aux.getWeightForWeights(
        path + weakSampleName, histoName="weightHisto_" + str(mNLSP), whichWeight="pu_mc_nisr")
    ewkWeight = aux.getWeightForWeights(
        path + weakSampleName, histoName="weightHisto_" + str(mNLSP), whichWeight="pu_mc_ewk")
    #Weight = aux.getWeightForWeights(path+weakSampleName,histoName="weightHisto_"+str(mNLSP),whichWeight="pu_mc_toppt")

    histo.Scale(topPtWeight * nIsrWeight * ewkWeight)
    #acc = histo.Integral()*100.

    #acc = histo.Integral(histo.FindFixBin(150.),-1)*100.
    acc = histo.Integral() * 100.
    strong.SetPoint(i, mGluino, mNeutralino, acc)
    i += 1

file_gmsb = ROOT.TFile(path + gmsbSampleName)
dirs = [k.GetName() for k in file_gmsb.GetListOfKeys()
        if k.GetName().startswith("GMSB")]

gmsb = TGraph2D()

k = 0
for key in dirs:
    folder = "sig/LL/nom/"
    hName = "met"
    point = getPointFromDir(key)
    mGluino = point[1]
    mNeutralino = point[2]
    histo = file_gmsb.Get(key + "/" + folder + hName)
    topPtWeight = aux.getWeightForWeights(
        path + weakSampleName, histoName="weightHisto_" + str(mNLSP), whichWeight="pu_mc_toppt")
    nIsrWeight = aux.getWeightForWeights(
        path + weakSampleName, histoName="weightHisto_" + str(mNLSP), whichWeight="pu_mc_nisr")
    ewkWeight = aux.getWeightForWeights(
        path + weakSampleName, histoName="weightHisto_" + str(mNLSP), whichWeight="pu_mc_ewk")
    #Weight = aux.getWeightForWeights(path+weakSampleName,histoName="weightHisto_"+str(mNLSP),whichWeight="pu_mc_toppt")

    histo.Scale(topPtWeight * nIsrWeight * ewkWeight)
    #acc = histo.Integral()*100.

    #acc = histo.Integral(histo.FindFixBin(150.),-1)*100.
    acc = histo.Integral() * 100.
    # gmsb.SetPoint(i,mGluino,mNeutralino,acc)
    gmsb.SetPoint(k, mNeutralino, mGluino, acc)
    k += 1


file_GGM12 = ROOT.TFile(path + GGM12SampleName)
dirs = [k.GetName() for k in file_GGM12.GetListOfKeys()
        if k.GetName().startswith("GGM")]
#print dirs

ggm1m2 = TGraph2D()

i = 0
for key in dirs:
    folder = "sig/LL/nom/"
    hName = "met"
    point = getPointFromDir(key)
    #print point
    m1 = point[1]
    m2 = point[2]
    histo = file_GGM12.Get(key + "/" + folder + hName)
    #acc = histo.Integral()*100.
    acc = histo.Integral(histo.FindFixBin(150.), -1) * 100.
    ggm1m2.SetPoint(i, m1, m2, acc)
    i += 1

file_GGM13 = ROOT.TFile(path + GGM13SampleName)
dirs = [k.GetName() for k in file_GGM13.GetListOfKeys()
        if k.GetName().startswith("GGM")]

ggm1m3 = TGraph2D()

i = 0
for key in dirs:
    folder = "sig/LL/nom/"
    hName = "met"
    point = getPointFromDir(key)
    m1 = point[1]
    m3 = point[2]
    histo = file_GGM13.Get(key + "/" + folder + hName)
    #acc = histo.Integral()*100.
    acc = histo.Integral(histo.FindFixBin(150.), -1) * 100.
    ggm1m3.SetPoint(i, m1, m3, acc)
    i += 1


import style
style.defaultStyle()
# style.defaultStyle()
#s= style.style2d()
# s.SetPadLeftMargin(0.18)
# style.setPaletteRWB()


weak.SetTitle("; m_{NLSP} (GeV); Acceptance x Efficiency [%]")

c = TCanvas("canvas", "", 800, 800)
# c = ROOT.TCanvas()

# l.Draw()

weak.SetMarkerStyle(20)
weak.SetMarkerSize(1)
# weak.Draw("AL*")
# weak.Draw("A*")
# weak.Draw("ALP")
# weak.Draw("A")

gs = ROOT.TGraphSmooth()
gs2 = ROOT.TGraphSmooth()
weak2 = gs.SmoothSuper(weak)
weak2.SetLineColor(ROOT.kBlack)
weak2.SetTitle("; m_{NLSP} (GeV); Acceptance x Efficiency [%]")
weak2.SetMarkerStyle(7)
weak2.SetMarkerSize(2)
# weak2.Draw("LP same")
weak2.Draw("ALP")
# weak3 = gs2.SmoothLowess(weak)
# weak3.SetLineColor(ROOT.kGreen)
# weak3.Draw("LP a4 same")
# weak3.Draw("a3 same")
# c.SetLogy()
#l = aux.Label(sim=True)
l = ROOT.TLatex(
    0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
l2 = ROOT.TLatex(0.17, .95, "#scale[0.76]{#font[52]{Simulation}}")
l.SetNDC()
l2.SetNDC()
# l.Draw()
# l2.Draw()
lum = ROOT.TLatex(.62, .95, "%.1f fb^{-1} (%s TeV)" %
                  (aux.intLumi / 1000., aux.Label.cmsEnergy))
lum.SetNDC()
# lum.Draw()
aux.Label(status="Simulation")
gm = ROOT.TLatex()
gm.DrawLatexNDC(0.3, 0.85, "#scale[0.76]{TChiZG}")
c.Update()
c.SaveAs('tching.pdf')


# style.defaultStyle()
style.style2d()
s = style.style2d()
# s.SetPadLeftMargin(0.18)
# c = TCanvas("canvas", "", 800, 800)
c = TCanvas()
strong.SetTitle(
    "; m_{#tilde{g}} (GeV);m_{NLSP} (GeV); Acceptance x Efficiency [%]")
# "; m_{#tilde{g}} (GeV);m_{#tilde{#chi_{0}^{1}}} (GeV); Acceptance x Efficiency [%]")
gmsb.SetTitle(
    "; m_{#tilde{B}} (GeV);m_{#tilde{W}} (GeV); Acceptance x Efficiency [%]")
ggm1m2.SetTitle("; M1;M2; Acceptance x Efficiency [%]")

lum = ROOT.TLatex(.62, .95, "#scale[0.76]{%.1f fb^{-1} (%s TeV)}" %
                  (aux.intLumi / 1000., aux.Label.cmsEnergy))
lum.SetNDC()
strong.GetZaxis().SetTitleOffset(1.42)
strong.GetYaxis().SetTitleOffset(gmsb.GetYaxis().GetTitleOffset() * .95)
strong.Draw("COLZ")
#l = aux.Label(sim=True)
# l.Draw()
# l2.Draw()
# lum.Draw()
aux.Label2D(status="Simulation")
gm = ROOT.TLatex()
gm.DrawLatexNDC(0.3, 0.85, "#scale[0.76]{T5bbbbZG}")
c.Update()
c.SaveAs("t5zg.pdf")

# c.Clear()
ggm1m2.Draw("COLZ")
#l = aux.Label(sim=True)
# l.Draw()
lum.Draw()
c.Update()
c.SaveAs("ggm1m2.pdf")

# c.Clear()
# gmsb.SetPoint(k+1,0,1300,.0)
gmsb.GetZaxis().SetTitleOffset(1.42)
gmsb.GetYaxis().SetTitleOffset(gmsb.GetYaxis().GetTitleOffset() * .95)
gmsb.Draw("COLZ")
# style.style2d()
c.Update()
c.Modified()
#l = aux.Label(sim=True)
# l.Draw()
# lum.Draw()
aux.Label2D(status="Simulation")
gm = ROOT.TLatex()
gm.DrawLatexNDC(0.3, 0.15, "#scale[0.76]{GMSB electroweak production}")

c.Update()
c.SaveAs("gmsb.pdf")


# c.Clear()
ggm1m3.Draw("COLZ")
#l = aux.Label(sim=True)
# l.Draw()
lum.Draw()
c.Update()
c.SaveAs("ggm1m3.pdf")
