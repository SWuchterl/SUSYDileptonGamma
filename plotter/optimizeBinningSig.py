#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
if sys.version_info[:2] == (2, 6):
    print "Initialize correct python version first!"
    sys.exit()


variable="mt2"
#variable="met"

from dataMC import labels, frange


import ConfigParser
import ROOT
import math
import argparse
import re
from random import randint
from sys import maxint

# private libs
import ratio
import style
import multiplot
from rwthColors import rwth
import limitTools

import auxiliary as aux
from datasets import *
import pickle as pkl


def getPoissonUnc(n):

    # calculate x = 1-alpha ( approx 68% )
    x = ROOT.TMath.Erf(1. / ROOT.TMath.Sqrt(2))
    alpha = 1 - x

    # for central confidence intervals, alpha_lo =alpha_up = alpha/2
    alpha_lo = alpha / 2
    alpha_up = alpha / 2

    # confidence interval is [ xlo, xup ]
    xlo = 0.5 * ROOT.TMath.ChisquareQuantile(alpha_lo, 2 * n)
    xup = 0.5 * ROOT.TMath.ChisquareQuantile(1 - alpha_up, 2 * (n + 1))
    return n - xlo, xup - n, n, xlo, xup


def writeSimplifiedDatacard(bkgInts, sigInts, name="tmpDataCard.txt"):
    # assuming no observation, 50% bkg uncert, 15% signal uncert
    rebinning = len(bkgInts)
    header = """
imax %d
jmax 1
kmax 2
bin %s
observation %s\n""" % (rebinning, " ".join(["b%d" % i for i in range(rebinning)]), " ".join([str(int(i)) for i in bkgInts]))
    infos = [["bin", "process", "process", "rate", "signal lnN", "bgUnc lnN"]]
    for iBin in range(len(bkgInts)):
        bkgInfos = ["b%d" % iBin, "bg", 1, bkgInts[iBin], "-", 1.5]
        sigInfos = ["b%d" % iBin, "sig", 0, sigInts[iBin], 1.15, "-"]
        bkgInfos = [str(i) for i in bkgInfos if i is not None]
        sigInfos = [str(i) for i in sigInfos if i is not None]
        infos.append(bkgInfos)
        infos.append(sigInfos)
    for l in zip(*infos):
        header += " ".join(l) + "\n"

    with open(name, "w") as f:
        f.write(header)


bkgSet = zgamma + ttgamma + zz + wwgamma + wzgamma + \
    DYjetsNLO + wjets + tt + singletop + wz + ww + zz4l + wjets
#sigSet = signal["T5Wg_1550_100"]
sigSet = tching_600
#sigSet = gmsb_240_230
#sigSet = gmsb_290_205
# sigSet3 = gmsb_290_205
sigSet3 = gmsb_415_355
# sigSet = t5bbbbzg_1500_600
# sigSet2 = t5bbbbzg_1500_1400
sigSet2 = t5bbbbzg_1500_600
sigSet4 = tching_400
sigSet5 = gmsb_240_230
# sigSet6 = t5bbbbzg_1800_1700
sigSet6 = t5bbbbzg_1500_1400

pklZZ = pkl.load(open("plots_CR_zz/factors/CRZZ.pkl", "rb"))
pklDY = pkl.load(open("plots_CR_dy/factors/CRDY.pkl", "rb"))
pklTT = pkl.load(open("plots_CR_tt/factors/CRTT.pkl", "rb"))
pklWZ = pkl.load(open("plots_CR_wz/factors/CRWZ.pkl", "rb"))
ZZsf = pklZZ["LL"]["m_ll"]
DYsf = pklDY["LL"]["eta1"]
TTsf = pklTT["EM"]["eta1"]
WZsf = pklWZ["LL"]["eta1"]

#rebinning = frange(150., 1100., 10.)
rebinning = frange(100., 700., 20.)


dirDir = "xx_0_0/sig/LL/nom"

zgHist = aux.stdHistWithoutNGenWithWeights(
    zgamma, dirDir + "/"+variable, ["nISR", "topPt", "ewk"], rebinning)
ttgHist = aux.stdHistWithoutNGenWithWeights(
    ttgamma, dirDir + "/"+variable, ["nISR", "topPt", "ewk"], rebinning)
zzHist = aux.stdHistWithoutNGenWithWeights(
    zz, dirDir + "/"+variable, ["nISR", "topPt", "ewk"], rebinning)
wwgHist = aux.stdHistWithoutNGenWithWeights(
    wwgamma, dirDir + "/"+variable, ["nISR", "topPt", "ewk"], rebinning)
wzgHist = aux.stdHistWithoutNGenWithWeights(
    wzgamma, dirDir + "/"+variable, ["nISR", "topPt", "ewk"], rebinning)
#dyHist = aux.stdHistWithoutNGenWithWeights(
    #DYjetsNLO, dirDir + "/met"+variable, ["nISR", "topPt", "ewk"], rebinning)
wjetsHist = aux.stdHistWithoutNGenWithWeights(
    wjets, dirDir + "/"+variable, ["nISR", "topPt", "ewk"], rebinning)
ttHist = aux.stdHistWithoutNGenWithWeights(
    tt, dirDir + "/"+variable, ["nISR", "topPt", "ewk"], rebinning)
singletopHist = aux.stdHistWithoutNGenWithWeights(
    singletop, dirDir + "/"+variable, ["nISR", "topPt", "ewk"], rebinning)
wzHist = aux.stdHistWithoutNGenWithWeights(
    wz, dirDir + "/"+variable, ["nISR", "topPt", "ewk"], rebinning)
wwHist = aux.stdHistWithoutNGenWithWeights(
    ww, dirDir + "/"+variable, ["nISR", "topPt", "ewk"], rebinning)
zz4lHist = aux.stdHistWithoutNGenWithWeights(
    zz4l, dirDir + "/"+variable, ["nISR", "topPt", "ewk"], rebinning)
wgHist = aux.stdHistWithoutNGenWithWeights(
    wgamma, dirDir + "/"+variable, ["nISR", "topPt", "ewk"], rebinning)

#final_dyHist = aux.addHists(dyHist, zgHist)
final_dyHist = aux.addHists( zgHist)
final_ttHist = aux.addHists(ttHist, ttgHist)
final_zzHist = aux.addHists(zzHist, zz4lHist)
final_wzHist = aux.addHists(wzHist)
final_otherHist = aux.addHists(
    wwgHist, wzgHist, wjetsHist, singletopHist, wwHist, wgHist)

# final_tt080Hist.Scale(sfTT[0])
# final_tt80Hist.Scale(sfTT[0])
final_ttHist.Scale(TTsf[0])
final_zzHist.Scale(ZZsf[0])
final_wzHist.Scale(WZsf[0])
final_dyHist.Scale(DYsf[0])


#bkgHist = bkgSet.getHist("tr/mt2")
#sigHist = sigSet.getHist("tr/mt2")
#bkgHist = bkgSet.getHist("onZMet150/LL/mt2")
# bkgHist = bkgSet.getHistWithoutNGen("xx_0_0/sig/LL/nom/mt2")
bkgHist = aux.addHists(final_ttHist, final_zzHist,
                       final_wzHist, final_otherHist)
print bkgHist.Integral(bkgHist.FindBin(150), bkgHist.FindBin(200) - 1)

#sigHist = sigSet.getHist("onZMet150/LL/mt2")
sigHist = sigSet.getHistWithoutNGen("Ng_0_0/sig/LL/nom/"+variable)
sigHist2 = sigSet2.getHistWithoutNGen("Zg_0_0/sig/LL/nom/"+variable)
sigHist3 = sigSet3.getHistWithoutNGen("GMSB_0_0/sig/LL/nom/"+variable)
sigHist4 = sigSet4.getHistWithoutNGen("Ng_0_0/sig/LL/nom/"+variable)
sigHist5 = sigSet5.getHistWithoutNGen("GMSB_0_0/sig/LL/nom/"+variable)
sigHist6 = sigSet6.getHistWithoutNGen("Zg_0_0/sig/LL/nom/"+variable)
#bkgHist = bkgSet.getHist("onZMet150/LL/m_llg")
#sigHist = sigSet.getHist("onZMet150/LL/m_llg")
#bkgHist = bkgSet.getHist("onZMet150/LL/zpt")
#sigHist = sigSet.getHist("onZMet150/LL/zpt")
#bkgHist = bkgSet.getHist("onZMet150/LL/deltaPhiLL")
#sigHist = sigSet.getHist("onZMet150/LL/deltaPhiLL")


# rebinning=frange(150.,5001.,5.)
# rebinning=frange(0.,1501.,25.)
# rebinning = frange(0., 1501., 50.)
# rebinning = frange(150., 1501., 50.)
# rebinning = frange(150., 1101., 50.)
# rebinning = frange(150., 1101., 10.)
# rebinning = frange(150., 260., 10.)
# rebinning=frange(0.,4.5,0.1)

bkgHist = aux.rebin(bkgHist, rebinning)
sigHist = aux.rebin(sigHist, rebinning)
sigHist2 = aux.rebin(sigHist2, rebinning)
sigHist3 = aux.rebin(sigHist3, rebinning)
sigHist4 = aux.rebin(sigHist4, rebinning)
sigHist5 = aux.rebin(sigHist5, rebinning)
sigHist6 = aux.rebin(sigHist6, rebinning)

binBoarders = [bkgHist.GetNbinsX() + 10]
#print binBoarders
oldSig = 0.  # maxint wouly be nicer


arCut = []
arSig = []
arSig2 = []
arSig3 = []
arSig4 = []
arSig5 = []
arSig6 = []
arBKG = []


#def returnSig(sA, bA):
    #arS = []
    #for i in range(len(sA)):
        #sig = ROOT.TMath.Sqrt(
            #2 * (sA[i] + bA[i]) * ROOT.TMath.Log(1 + sA[i] / bA[i]) - 2 * sA[i])
        #arS.append(sig)
    #return max(arS)
def returnSig(sA, bA):
    sig = ROOT.TMath.Sqrt(
        2 * (sA + bA) * ROOT.TMath.Log(1 + sA / bA) - 2 * sA)
    return sig


def returnSigAs(s, b, sigB):
    if s == 0:
        s = 0.01
    part = (s + b) * ROOT.TMath.Log((s + b) * (b + sigB * sigB) / (b * b + (s + b) * sigB * sigB)) - \
        (b * b) / (sigB * sigB) * ROOT.TMath.Log(1 +
                                                 (sigB * sigB * s) / (b * (b + sigB * sigB)))
    print part, s, b, sigB
    return (2. * part)**(0.5)

# def sig2(n,b,sigB):


def returnSigAs2(n, b, sigB):
    if n - b == 0:
        n = b + 0.01
    part = (n) * ROOT.TMath.Log((n) * (b + sigB * sigB) / (b * b + (n) * sigB * sigB)) - (b * b) / \
        (sigB * sigB) * ROOT.TMath.Log(1 +
                                       (sigB * sigB * (n - b)) / (b * (b + sigB * sigB)))
    return (2. * part)**(0.5)

# def returnSig2(sA, bA):
#     if bA == 0:
#         return 0.01
#     sig = ROOT.TMath.Sqrt(
#         2 * (sA + bA) * ROOT.TMath.Log(1 + sA / bA) - 2 * sA)
#     return sig


def returnSig2(sA, bA):
    if bA == 0:
        return 0.01
    sig = sA / ROOT.TMath.Sqrt(bA)
    return sig
def returnSig3(sA, bA):
    if bA == 0:
        return 0.01
    sig = sA / ROOT.TMath.Sqrt(sA+bA)
    return sig


# for bin in range( sigHist.GetNbinsX()+2, 0, -1 ):
# for bin in range(sigHist.GetNbinsX() + 2, 0, -1):
# rebinning2 = frange(150, 200, 10)
#rebinning2 = frange(150, 1000, 50)
#rebinning2 = frange(150, 1000, 10)
#rebinning2 = frange(100, 1000, 20)
rebinning2 = frange(100, 400, 20)
# rebinning2 = frange(150, 1000, 10)
# for cut in rebinning:
for cut in rebinning2:

    # bkgInts = [bkgHist.Integral(binBoarders[i + 1], binBoarders[i])
    #            for i in range(len(binBoarders) - 1)]
    # sigInts = [sigHist.Integral(binBoarders[i + 1], binBoarders[i])
    #            for i in range(len(binBoarders) - 1)]
    bkgInt = bkgHist.Integral(bkgHist.FindBin(cut), -1)
    sigInt = sigHist.Integral(sigHist.FindBin(cut), -1)
    sigInt2 = sigHist2.Integral(sigHist2.FindBin(cut), -1)
    sigInt3 = sigHist3.Integral(sigHist3.FindBin(cut), -1)
    sigInt4 = sigHist4.Integral(sigHist4.FindBin(cut), -1)
    # sigInt5 = sigHist5.Integral(sigHist5.FindBin(cut), -1)
    sigInt6 = sigHist6.Integral(sigHist6.FindBin(cut), -1)
    # bkgInt = bkgHist.Integral(bkgHist.FindBin(cut), bkgHist.FindBin(200) - 1)
    # sigInt = sigHist.Integral(sigHist.FindBin(cut), bkgHist.FindBin(200) - 1)
    # sigInt2 = sigHist2.Integral(
    #     sigHist2.FindBin(cut), bkgHist.FindBin(200) - 1)
    # sigInt3 = sigHist3.Integral(
    #     sigHist3.FindBin(cut), bkgHist.FindBin(200) - 1)
    # sigInt4 = sigHist4.Integral(
    #     sigHist4.FindBin(cut), bkgHist.FindBin(200) - 1)
    # sigInt5 = sigHist5.Integral(
    #     sigHist5.FindBin(cut), bkgHist.FindBin(200) - 1)
    # sigInt6 = sigHist6.Integral(
    #     sigHist6.FindBin(cut), bkgHist.FindBin(200) - 1)

    # bkgInts.append(bkgHist.Integral(bin, binBoarders[-1]))
    # sigInts.append(sigHist.Integral(bin, binBoarders[-1]))

    # if min(bkgInts) < 1e-6:
    # if min(bkgInts) < 0.5:
    # continue
    # if min(sigInts) < 1e-6:
    # continue
    # print "bkg", bkgInts, "signal", sigInts
    print "bkg", bkgInt, "signal", sigInt

    # sig = returnSig2(sigInt, bkgInt)
    # sig2 = returnSig2(sigInt2, bkgInt)
    # sig3 = returnSig2(sigInt3, bkgInt)
    # sig4 = returnSig2(sigInt4, bkgInt)
    # sig5 = returnSig2(sigInt5, bkgInt)
    # sig6 = returnSig2(sigInt6, bkgInt)
    sig = returnSig(sigInt, bkgInt)
    sig2 = returnSig(sigInt2, bkgInt)
    sig3 = returnSig(sigInt3, bkgInt)
    sig4 = returnSig(sigInt4, bkgInt)
    ####sig5 = returnSig2(sigInt5, bkgInt)
    sig6 = returnSig(sigInt6, bkgInt)
    #sig = returnSigAs(sigInt, bkgInt, getPoissonUnc(bkgInt)[0])
    #sig2 = returnSigAs(sigInt2, bkgInt, getPoissonUnc(bkgInt)[0])
    #sig3 = returnSigAs(sigInt3, bkgInt, getPoissonUnc(bkgInt)[0])
    #sig4 = returnSigAs(sigInt4, bkgInt, getPoissonUnc(bkgInt)[0])
    ## sig5 = returnSigAs(sigInt5, bkgInt, getPoissonUnc(bkgInt)[0])
    #sig6 = returnSigAs(sigInt6, bkgInt, getPoissonUnc(bkgInt)[0])
    #sig = returnSigAs2(sigInt + bkgInt, bkgInt, getPoissonUnc(bkgInt)[0])
    #sig2 = returnSigAs2(sigInt2 + bkgInt, bkgInt, getPoissonUnc(bkgInt)[0])
    #sig3 = returnSigAs2(sigInt3 + bkgInt, bkgInt, getPoissonUnc(bkgInt)[0])
    #sig4 = returnSigAs2(sigInt4 + bkgInt, bkgInt, getPoissonUnc(bkgInt)[0])
     ###sig5 = returnSigAs2(sigInt5+bkgInt, bkgInt, getPoissonUnc(bkgInt)[0])
    #sig6 = returnSigAs2(sigInt6 + bkgInt, bkgInt, getPoissonUnc(bkgInt)[0])

    # arCut.append(bkgHist.GetBinLowEdge(
    #     bin))
    arCut.append(cut)
    arSig.append(sig)
    arSig2.append(sig2)
    arSig3.append(sig3)
    arSig4.append(sig4)
    # arSig5.append(sig5)
    arSig6.append(sig6)
    # arBKG.append(bkgInts[len(binBoarders) - 1])
    # arBKG.append(bkgInts[0])
    arBKG.append(bkgInt)

    # writeSimplifiedDatacard(
    #     bkgInts, sigInts, name="tmpDataCard_" + sigSet.names[0] + "_.txt")
    # r = limitTools.infosFromDatacard(
    #     "tmpDataCard_" + sigSet.names[0] + "_.txt")["exp"]
    # print "bin", bin, "lowEdge", bkgHist.GetBinLowEdge(
    #     bin), "sig", sig, "oldSig", oldSig
    print "bin", bin, "lowEdge", cut, "sig", sig, "oldSig", oldSig
    if (sig - oldSig) / sig < -0.05:  # change must me larger than 5%
        print "append"
        # binBoarders.append(bin)
    oldSig = sig

print "final sig =", sig
# print binBoarders
# metBoarders = [sigHist.GetBinLowEdge(i) for i in binBoarders]
# metBoarders.reverse()
# print metBoarders[0:-1]


import style
style.defaultStyle()


# weak.SetTitle("; m_{NLSP} (GeV); signal fraction [%]")

c = ROOT.TCanvas("canvas", "", 800, 800)


# weak.Draw("A*")
# l = ROOT.TLatex(
#     # 0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
#     0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Private Work}}")
# l2 = ROOT.TLatex(0.21, .88, "#scale[0.76]{#font[52]{Simulation}}")
# l.SetNDC()
# l2.SetNDC()
# l.Draw()
# l2.Draw()
# lum = ROOT.TLatex(.62, .95, "%.1f fb^{-1} (%s TeV)" %
#                   (aux.intLumi / 1000., aux.Label.cmsEnergy))
# lum.SetNDC()
# lum.Draw()
import style
style.defaultStyle()
grSig = ROOT.TGraph()
grSig2 = ROOT.TGraph()
grSig3 = ROOT.TGraph()
grSig4 = ROOT.TGraph()
grSig5 = ROOT.TGraph()
grSig6 = ROOT.TGraph()
grBKG = ROOT.TGraph()
# grCut=TGraph()

for i in range(len(arSig)):
    grSig.SetPoint(i, arCut[i], arSig[i])
    grSig2.SetPoint(i, arCut[i], arSig2[i])
    grSig3.SetPoint(i, arCut[i], arSig3[i])
    grSig4.SetPoint(i, arCut[i], arSig4[i])
    # grSig5.SetPoint(i, arCut[i], arSig5[i])
    grSig6.SetPoint(i, arCut[i], arSig6[i])
    grBKG.SetPoint(i, arCut[i], arBKG[i])
    # grCut.SetPoint(arCut[i])
grBKG.SetLineWidth(2)
grSig.SetLineWidth(2)
grSig2.SetLineWidth(2)
grSig3.SetLineWidth(2)
grSig4.SetLineWidth(2)
grSig5.SetLineWidth(2)
grSig6.SetLineWidth(2)
grSig.SetMarkerStyle(20)
grSig2.SetMarkerStyle(20)
grSig3.SetMarkerStyle(20)
grSig4.SetMarkerStyle(20)
grSig5.SetMarkerStyle(20)
grSig6.SetMarkerStyle(20)
grSig.SetMarkerSize(0.5)
grSig2.SetMarkerSize(0.5)
grSig3.SetMarkerSize(0.5)
grSig4.SetMarkerSize(0.5)
grSig5.SetMarkerSize(0.5)
grSig6.SetMarkerSize(0.5)
grSig.SetMarkerColor(ROOT.kRed + 3)
grSig2.SetMarkerColor(ROOT.kGreen + 3)
grSig3.SetMarkerColor(ROOT.kBlue + 3)
grSig4.SetMarkerColor(ROOT.kRed + 2)
grSig5.SetMarkerColor(ROOT.kBlue - 2)
grSig6.SetMarkerColor(ROOT.kGreen + 2)
grSig.SetLineColor(ROOT.kRed + 3)
grSig2.SetLineColor(ROOT.kGreen + 3)
grSig3.SetLineColor(ROOT.kBlue + 3)
grSig4.SetLineColor(ROOT.kRed + 2)
grSig5.SetLineColor(ROOT.kBlue - 2)
grSig6.SetLineColor(ROOT.kGreen + 2)
grBKG.Draw("APL")
grSig.Draw("PL")
grSig2.Draw("PL")
grSig3.Draw("PL")
grSig4.Draw("PL")

leg = ROOT.TLegend(0.45, 0.6, 0.9, 0.9)
leg.AddEntry(grBKG, "Number of Background Events", "lp")
leg.AddEntry(grSig, "Z_{A} TChiZG m(NLSP)=600 GeV", "lp")
leg.AddEntry(grSig4, "Z_{A} TChiZG m(NLSP)=400 GeV", "lp")
leg.AddEntry(
    grSig3, "Z_{A} GMSB m(#tilde{W})=415 GeV m(#tilde{B})=355 GeV", "lp")
leg.AddEntry(
    grSig2, "Z_{A} T5bbbbZG m(#tilde{g})=1500 GeV m(#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}})=600 GeV", "lp")
leg.AddEntry(
    grSig6, "Z_{A} T5bbbbZG m(#tilde{g})=1500 GeV m(#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}})=1400 GeV", "lp")
# leg.AddEntry(grSig2,"NG 400")
# grSig5.Draw("P")
grSig6.Draw("PL")
# grBKG.Draw("PL")
line = ROOT.TLine(200, 0, 200, 6)
line2 = ROOT.TLine(150, 3, 1000, 3)
# line.Draw("same")
# line2.Draw("same")
grBKG.GetXaxis().SetRangeUser(100, 1000)
grBKG.SetTitle("; p_{T}^{miss} (GeV); #splitline{Z_{A}}{Bkg. Events}")
l2 = ROOT.TLatex(0.21, .88, "#scale[0.76]{#font[52]{Simulation}}")
l2.SetNDC()
l2.Draw()
lum = ROOT.TLatex(.62, .95, "%.1f fb^{-1} (%s TeV)" %
                  (aux.intLumi / 1000., aux.Label.cmsEnergy))
lum.SetNDC()
lum.Draw()
grBKG.GetYaxis().SetTitleOffset(1.)
leg.Draw()

c.Update()
c.SaveAs('optimizeBinning/'+variable+'.pdf')
