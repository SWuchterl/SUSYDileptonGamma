from ROOT import *
import pickle as pkl
import ROOT
from array import array
from include import *
import numpy as np


def getPointFromDir(name):
    m = re.match("(.*)_(.*)_(.*)", name)
    combi, m1, m2 = m.groups()
    m1, m2 = int(m1), int(m2)
    return combi, m1, m2


# path="../myAnalyzer/output_mediumPOG/"
# path="../myAnalyzer/output/"
# path="../myAnalyzer/output_2/"
path = "../myAnalyzer/output_AN/"

lumi = 35867.
xSec_tching = pkl.load(open("../SUSYxSections/xSec_SMS_C1C1_13TeV.pkl", "rb"))
xSec_tching_neutr = pkl.load(
    open("../SUSYxSections/xSec_SMS_N2C1_13TeV.pkl", "rb"))
xSec_tching_total = pkl.load(open("data/xSec_SMS_TChiNG_13TeV.pkl", "rb"))
xSec_t5zg = pkl.load(open("data/xSec_SMS_Gluino_13TeV.pkl", "rb"))

xSec_gmsb = pkl.load(open("data/xSec_GMSB.pkl", "rb"))
xSec_t6 = pkl.load(open("data/xSec_SMS_Stop_13TeV.pkl"))

xSec_GGM_M1_M2 = pkl.load(open("data/xsec_GGM_M1_M2.pkl"))
xSec_GGM_M1_M3 = pkl.load(open("data/xsec_GGM_M1_M3.pkl"))

brZG_tching = 0.25
brZG_t5zg = 0.5


weakSampleName = "SMS-TChiNG_BF50N50G_hists.root"
strongSampleName = "SMS-T5bbbbZg_hists.root"
GGM12SampleName = "GGM_GravitinoLSP_M1-200to1500_M2-200to1500_hists.root"
GGM13SampleName = "GGM_GravitinoLSP_M1-50to1500_M3-1000to2500_hists.root"

gmsbSampleName = "GMSB_GravitinoLSP_N1decays_hists.root"

t6SampleName = "SMS-T6ttZg_hists.root"

folder = "CRZZ/LL/nom/"
hName = "met"

style.additionalPoissonUncertainty = False

zgHist = aux.stdHistWithWeights(
    zgamma, folder + hName, ["nISR", "topPt", "ewk"])
ttgHist = aux.stdHistWithWeights(
    ttgamma, folder + hName, ["nISR", "topPt", "ewk"])
zzHist = aux.stdHistWithWeights(zz, folder + hName, ["nISR", "topPt", "ewk"])
wwgHist = aux.stdHistWithWeights(
    wwgamma, folder + hName, ["nISR", "topPt", "ewk"])
wzgHist = aux.stdHistWithWeights(
    wzgamma, folder + hName, ["nISR", "topPt", "ewk"])
dyHist = aux.stdHistWithWeights(
    DYjetsNLO, folder + hName, ["nISR", "topPt", "ewk"])
wjetsHist = aux.stdHistWithWeights(
    wjets, folder + hName, ["nISR", "topPt", "ewk"])
ttHist = aux.stdHistWithWeights(tt, folder + hName, ["nISR", "topPt", "ewk"])
singletopHist = aux.stdHistWithWeights(
    singletop, folder + hName, ["nISR", "topPt", "ewk"])
wzHist = aux.stdHistWithWeights(wz, folder + hName, ["nISR", "topPt", "ewk"])
wwHist = aux.stdHistWithWeights(ww, folder + hName, ["nISR", "topPt", "ewk"])
zz4lHist = aux.stdHistWithWeights(
    zz4l, folder + hName, ["nISR", "topPt", "ewk"])
wgHist = aux.stdHistWithWeights(
    wgamma, folder + hName, ["nISR", "topPt", "ewk"])


histsToScale = [zgHist, ttgHist, zzHist, wwgHist, wzgHist, dyHist,
                wjetsHist, ttHist, singletopHist, wzHist, wwHist, zz4lHist, wgHist]


pklZZ = pkl.load(open("plots_CR_zz/factors/CRZZ.pkl", "rb"))
pklDY = pkl.load(open("plots_CR_dy/factors/CRDY.pkl", "rb"))
pklTT = pkl.load(open("plots_CR_tt/factors/CRTT.pkl", "rb"))
pklWZ = pkl.load(open("plots_CR_wz/factors/CRWZ.pkl", "rb"))
sfZZ, sfZZErr = pklZZ["LL"]["m_ll"]
sfDY, sfDYErr = pklDY["LL"]["eta1"]
sfTT = pklTT["EM"]["eta1"][0]
sfTT080Err = 0.04
sfTT80Err = 0.4
sfWZ, sfWZErr = pklWZ["LL"]["eta1"]


zzHist.Scale(sfZZ)
zz4lHist.Scale(sfZZ)
dyHist.Scale(sfDY)
zgHist.Scale(sfDY)
wzHist.Scale(sfWZ)
ttHist.Scale(sfTT)
ttgHist.Scale(sfTT)

dirHist = aux.addHists(*histsToScale)


yieldBKG = dirHist.Integral()


file_weak = ROOT.TFile(path + weakSampleName)
dirs = [k.GetName() for k in file_weak.GetListOfKeys()
        if k.GetName().startswith("Ng")]

weak = TGraph()


contGMSB = {}
contWEAK = {}
contSTRONG = {}
contM1M2 = {}
contM1M3 = {}
contT6 = {}

if(True):

    i = 0
    for key in dirs:
        #print key
        name, m1, m2 = key.split("_")
        folder = "CRTT/EM/nom/"
        hName = "met"
        point = getPointFromDir(key)
        mNLSP = point[1]
        mDummy = point[2]
        histo = file_weak.Get(key + "/" + folder + hName)

        #print point[1]
        sftoppt = aux.getWeightForWeights(
            path + weakSampleName, histoName="weightHisto" + "_" + str(mNLSP), whichWeight="pu_mc_toppt")
        sfnisr = aux.getWeightForWeights(
            path + weakSampleName, histoName="weightHisto" + "_" + str(mNLSP), whichWeight="pu_mc_nisr")
        sfewk = aux.getWeightForWeights(
            path + weakSampleName, histoName="weightHisto" + "_" + str(mNLSP), whichWeight="pu_mc_ewk")

        histo.Scale(sftoppt * sfnisr * sfewk)
        xSec = xSec_tching_total[point[1]][0]

        acc = histo.Integral()  # for no Higgs Decays
        yieldSignal = acc * xSec * lumi

        #style.additionalPoissonUncertainty = False

        cont = yieldSignal / yieldBKG * 100.

        # weak.SetPoint(i,mNLSP,acc)
        weak.SetPoint(i, mNLSP, cont)
        contWEAK[mNLSP] = cont / 100.
        i += 1


file_strong = ROOT.TFile(path + strongSampleName)
dirs = [k.GetName() for k in file_strong.GetListOfKeys()
        if k.GetName().startswith("Zg")]

strong = TGraph2D()

i = 0
for key in dirs:
    folder = "CRTT/EM/nom/"
    hName = "met"
    point = getPointFromDir(key)
    mGluino = point[1]
    mNeutralino = point[2]
    histo = file_strong.Get(key + "/" + folder + hName)

    xSec = xSec_t5zg[mGluino][0]

    sftoppt = aux.getWeightForWeights(path + strongSampleName, histoName="weightHisto" + "_" + str(
        mGluino) + "_" + str(mNeutralino), whichWeight="pu_mc_toppt")
    sfnisr = aux.getWeightForWeights(path + strongSampleName, histoName="weightHisto" + "_" + str(
        mGluino) + "_" + str(mNeutralino), whichWeight="pu_mc_nisr")
    sfewk = aux.getWeightForWeights(path + strongSampleName, histoName="weightHisto" +
                                    "_" + str(mGluino) + "_" + str(mNeutralino), whichWeight="pu_mc_ewk")

    histo.Scale(sftoppt * sfnisr * sfewk)

    acc = histo.Integral()
    yieldSignal = acc * xSec * lumi
    cont = yieldSignal / yieldBKG * 100.

    strong.SetPoint(i, mGluino, mNeutralino, cont)
    if(not bool(contSTRONG.get(mGluino))):
        contSTRONG[mGluino] = {}
    contSTRONG[mGluino][mNeutralino] = cont / 100.
    i += 1

file_t6 = ROOT.TFile(path + t6SampleName)
dirs = [k.GetName() for k in file_t6.GetListOfKeys()
        if k.GetName().startswith("Zg")]

t6 = TGraph2D()

i = 0
for key in dirs:
    folder = "CRTT/EM/nom/"
    hName = "met"
    point = getPointFromDir(key)
    mGluino = point[1]
    mNeutralino = point[2]
    histo = file_t6.Get(key + "/" + folder + hName)

    xSec = xSec_t6[mGluino][0]

    sftoppt = aux.getWeightForWeights(path + t6SampleName, histoName="weightHisto" + "_" + str(
        mGluino) + "_" + str(mNeutralino), whichWeight="pu_mc_toppt")
    sfnisr = aux.getWeightForWeights(path + t6SampleName, histoName="weightHisto" + "_" + str(
        mGluino) + "_" + str(mNeutralino), whichWeight="pu_mc_nisr")
    sfewk = aux.getWeightForWeights(path + t6SampleName, histoName="weightHisto" +
                                    "_" + str(mGluino) + "_" + str(mNeutralino), whichWeight="pu_mc_ewk")

    histo.Scale(sftoppt * sfnisr * sfewk)

    acc = histo.Integral()
    yieldSignal = acc * xSec * lumi
    cont = yieldSignal / yieldBKG * 100.

    t6.SetPoint(i, mGluino, mNeutralino, cont)
    if(not bool(contT6.get(mGluino))):
        contT6[mGluino] = {}
    contT6[mGluino][mNeutralino] = cont / 100.
    i += 1

file_gmsb = ROOT.TFile(path + gmsbSampleName)
dirs = [k.GetName() for k in file_gmsb.GetListOfKeys()
        if k.GetName().startswith("GMSB")]

gmsb = TGraph2D()

k = 0
for key in dirs:
    #print k,"/",len(dirs)
    folder = "CRZZ/LL/nom"
    hName = "/met"
    point = getPointFromDir(key)
    mGluino = point[1]
    mNeutralino = point[2]

    xSec = xSec_gmsb[point[1]][point[2]][0]

    sftoppt = aux.getWeightForWeights(path + gmsbSampleName, histoName="weightHisto" + "_" + str(
        mGluino) + "_" + str(mNeutralino), whichWeight="pu_mc_toppt")
    sfnisr = aux.getWeightForWeights(path + gmsbSampleName, histoName="weightHisto" + "_" + str(
        mGluino) + "_" + str(mNeutralino), whichWeight="pu_mc_nisr")
    sfewk = aux.getWeightForWeights(path + gmsbSampleName, histoName="weightHisto" + "_" + str(
        mGluino) + "_" + str(mNeutralino), whichWeight="pu_mc_ewk")
    histo = file_gmsb.Get(key + "/" + folder + hName)
    histo.Scale(sftoppt * sfnisr * sfewk)

    acc = histo.Integral()
    yieldSignal = acc * xSec * lumi

    #print acc

    cont = yieldSignal / yieldBKG * 100.

    #print point,xSec,yieldBKG,yieldSignal,cont,acc

    gmsb.SetPoint(k, mNeutralino, mGluino, cont)
    if(not bool(contGMSB.get(mGluino))):
        contGMSB[mGluino] = {}
    contGMSB[mGluino][mNeutralino] = cont / 100.
    k += 1


file_GGM12 = ROOT.TFile(path + GGM12SampleName)
dirs = [k.GetName() for k in file_GGM12.GetListOfKeys()
        if k.GetName().startswith("GGM")]

ggm1m2 = TGraph2D()

i = 0
for key in dirs:
    folder = "CRTT/EM/nom/"
    hName = "met"
    point = getPointFromDir(key)
    m1 = point[1]
    m2 = point[2]
    sftoppt = aux.getWeightForWeights(
        path + GGM12SampleName, histoName="weightHisto" + "_" + str(m1) + "_" + str(m2), whichWeight="pu_mc_toppt")
    sfnisr = aux.getWeightForWeights(
        path + GGM12SampleName, histoName="weightHisto" + "_" + str(m1) + "_" + str(m2), whichWeight="pu_mc_nisr")
    sfewk = aux.getWeightForWeights(
        path + GGM12SampleName, histoName="weightHisto" + "_" + str(m1) + "_" + str(m2), whichWeight="pu_mc_ewk")

    histo = file_GGM12.Get(key + "/" + folder + hName)
    histo.Scale(sftoppt * sfnisr * sfewk)

    #acc = histo.Integral(histo.FindFixBin(150.),-1)*100.
    acc = histo.Integral() * 100.
    yieldSignal = acc * xSec * lumi
    cont = yieldSignal / yieldBKG * 100.
    ggm1m2.SetPoint(i, m1, m2, cont)
    if(not bool(contM1M2.get(m1))):
        contM1M2[m1] = {}
    contM1M2[m1][m2] = cont / 100.
    i += 1

file_GGM13 = ROOT.TFile(path + GGM13SampleName)
dirs = [k.GetName() for k in file_GGM13.GetListOfKeys()
        if k.GetName().startswith("GGM")]

ggm1m3 = TGraph2D()

i = 0
for key in dirs:
    folder = "CRTT/EM/nom/"
    hName = "met"
    point = getPointFromDir(key)
    m1 = point[1]
    m3 = point[2]
    sftoppt = aux.getWeightForWeights(
        path + GGM13SampleName, histoName="weightHisto" + "_" + str(m1) + "_" + str(m3), whichWeight="pu_mc_toppt")
    sfnisr = aux.getWeightForWeights(
        path + GGM13SampleName, histoName="weightHisto" + "_" + str(m1) + "_" + str(m3), whichWeight="pu_mc_nisr")
    sfewk = aux.getWeightForWeights(
        path + GGM13SampleName, histoName="weightHisto" + "_" + str(m1) + "_" + str(m3), whichWeight="pu_mc_ewk")

    histo = file_GGM13.Get(key + "/" + folder + hName)
    histo.Scale(sftoppt * sfnisr * sfewk)

    acc = histo.Integral() * 100.
    yieldSignal = acc * xSec * lumi
    cont = yieldSignal / yieldBKG * 100.

    ggm1m3.SetPoint(i, m1, m3, cont)
    if(not bool(contM1M3.get(m1))):
        contM1M3[m1] = {}
    contM1M3[m1][m3] = cont / 100.
    i += 1


import style
style.defaultStyle()


weak.SetTitle("; m_{NLSP} (GeV); signal fraction [%]")

c = TCanvas("canvas", "", 800, 800)


weak.SetMarkerStyle(20)
weak.SetMarkerSize(1)

weak.Draw("A*")
l = ROOT.TLatex(
    0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Private Work}}")
# 0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
l2 = ROOT.TLatex(0.21, .88, "#scale[0.76]{#font[52]{Simulation}}")
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
gm.DrawLatexNDC(0.2, 0.8, "#scale[0.76]{TChiZG}")
gm.DrawLatexNDC(
    0.2, 0.85, "#scale[0.66]{#font[52]{ZZ Control Region}}")
c.Update()
c.SaveAs('contamination/tching_ZZ.pdf')


style.style2d()
s = style.style2d()
# s.SetPadLeftMargin(0.18)
# c = TCanvas("canvas", "", 800, 800)
c = TCanvas()
# strong.SetTitle("; m_{#tilde{g}} (GeV);m_{#tilde{#chi_{0}^{1}}} (GeV); Acceptance x Efficiency [%]")
gmsb.SetTitle("; m_{#tilde{B}} (GeV);m_{#tilde{W}} (GeV); signal fraction [%]")
ggm1m2.SetTitle(
    "; m_{#tilde{B}} (GeV);m_{#tilde{W}} (GeV); signal fraction [%]")
ggm1m3.SetTitle(
    "; m_{#tilde{B}} (GeV);m_{#tilde{W}} (GeV); signal fraction [%]")
strong.SetTitle(
    "; m_{#tilde{g}} (GeV);m_{NLSP} (GeV); signal fraction [%]")
t6.SetTitle(
    "; m_{#tilde{B}} (GeV);m_{#tilde{W}} (GeV); signal fraction [%]")
#ggm1m2.SetTitle("; M1;M2; Acceptance x Efficiency [%]")

lum = ROOT.TLatex(.62, .95, "#scale[0.76]{%.1f fb^{-1} (%s TeV)}" %
                  (aux.intLumi / 1000., aux.Label.cmsEnergy))
lum.SetNDC()
strong.Draw("COLZ")
# l.Draw()
# l2.Draw()
# lum.Draw()
aux.Label2D(status="Simulation")
gm = ROOT.TLatex()
gm.DrawLatexNDC(0.2, 0.85, "#scale[0.76]{T5bbbbZG}")
gm.DrawLatexNDC(
    0.2, 0.8, "#scale[0.66]{#font[52]{ZZ Control Region}}")
c.Update()
c.SaveAs("contamination/t5zg_ZZ.pdf")

ggm1m2.Draw("COLZ")
# l.Draw()
# lum.Draw()
# l2.Draw()
c.Update()
c.SaveAs("contamination/ggm1m2_ZZ.pdf")

t6.Draw("COLZ")
# l.Draw()
# lum.Draw()
# l2.Draw()
c.Update()
c.SaveAs("contamination/t6_ZZ.pdf")


gmsb.Draw("COLZ")
# l.Draw()
# lum.Draw()
# l2.Draw()
aux.Label2D(status="Simulation")
gm = ROOT.TLatex()
gm.DrawLatexNDC(0.3, 0.15, "#scale[0.76]{GMSB electroweak production}")
gm.DrawLatexNDC(
    0.3, 0.2, "#scale[0.66]{#font[52]{ZZ Control Region}}")
c.Update()
c.SaveAs("contamination/gmsb_ZZ.pdf")

import pickle as pkl
pkl.dump(contWEAK, open("contamination/ZZ_tching.pkl", "wb"))
print "Created contamination/ZZ_tching.pkl."
pkl.dump(contSTRONG, open("contamination/ZZ_t5zg.pkl", "wb"))
print "Created contamination/ZZ_t5zg.pkl."
pkl.dump(contGMSB, open("contamination/ZZ_gmsb.pkl", "wb"))
print "Created contamination/ZZ_gmsb.pkl."
pkl.dump(contM1M2, open("contamination/ZZ_m1m2.pkl", "wb"))
print "Created contamination/ZZ_m1m2.pkl."
pkl.dump(contM1M3, open("contamination/ZZ_m1m3.pkl", "wb"))
print "Created contamination/ZZ_m1m3.pkl."
pkl.dump(contT6, open("contamination/ZZ_t6.pkl", "wb"))
print "Created contamination/ZZ_t6.pkl."

ggm1m3.Draw("COLZ")
l.Draw()
lum.Draw()
l2.Draw()
c.Update()
c.SaveAs("contamination/ggm1m3_ZZ.pdf")
