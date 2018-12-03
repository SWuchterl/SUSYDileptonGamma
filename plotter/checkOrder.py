from dataMC import labels, frange

import ROOT
from ROOT import *
from array import array
from include import *
import numpy as np
import pickle as pkl
import os


import matplotlib.pyplot as plt

from CR_tt import calculateSFAndError, drawTTCR
from CR_DY import drawDYCR
from CR_WZ import drawCRWZ
from CR_ZZ import drawCRZZ

binnings = {
    # 'pt1':              frange(20,100,10)+frange(100,200,25)+range(200,350,50),
    # 'pt1':              frange(25,160,5),
    # 'pt1':              frange(20,150,10),
    # 'pt1':              frange(20,150,10),
    'pt1':              frange(20., 200., 10.),
    # 'pt2':              frange(20,100,10)+frange(100,200,25),
    'pt2':              frange(20, 150, 10),
    'pt3':              range(0, 200, 10),
    'pt4':              range(0, 200, 10),
    'eta1':             frange(0., 2.45, 0.1),
    # 'eta1':             frange(0., 2.4, 0.05),
    'eta2':             frange(0., 2.45, 0.1),
    # 'eta2':             frange(0., 2.4, 0.05),
    'eta3':             frange(0., 2.45, 0.1),
    'eta4':             frange(0., 2.45, 0.1),
    'phi1':             frange(0., 3.15, 0.1),
    'phi2':             frange(0., 3.15, 0.1),
    'phi3':             frange(0., 3.15, 0.1),
    'phi4':             frange(0., 3.15, 0.1),
    'ht':               frange(0., 1000., 50),
    # 'met':               [0,25,50,75,100,150,190,230,500],
    'met':               [0, 25, 50, 75, 100, 150, 250, 450],
    'm_ll':             frange(50., 100., 10) + frange(100., 300., 20.),
    'm_ll2':             frange(50., 100., 10) + frange(100., 300., 20.),
    'm_llg':             frange(0, 100, 10) + frange(100, 200, 10) + frange(200, 500, 50),
    'pt_llg':             frange(0, 100, 10) + frange(100, 200, 10) + frange(200, 500, 50),
    'n_jets':           frange(0., 7., 1),
    'n_photons':         frange(0., 4., 1),
    'n_vtx':            frange(0., 40., 1),
    # 'pt_g1':            frange(20,100,10)+frange(100,150,25)+frange(150,250,50),
    "pt_g1":            frange(25, 80, 5) + frange(80, 140, 10) + frange(140, 200, 20),
    'eta_g1':           frange(0., 2.60, 0.1),
    'phi_g1':           frange(0., 3.50, 0.1),
    'sigmaIetaIeta_g1': frange(0., 0.04, 0.01),
    'sigmaIphiIphi_g1': frange(0., 0.2, 0.01),
    'deltaR1_g1':       frange(0., 1., 0.05),
    'deltaR2_g1':       frange(0., 1., 0.05),
    'r9_g1':            frange(0., 1.5, 0.05),
    'hOverE_g1':        frange(0., 0.1, 0.01),
    'deltaEtaLL':       frange(0., 6., 0.1),
    'deltaPhiLL':       frange(0., 6., 0.1),
    'deltaEtaLLG':      frange(0., 6., 0.1),
    'deltaPhiLLG':      frange(0., 6., 0.01),
    'deltaRLL':         frange(0., 1., 0.1),
    'deltaRLLG':        frange(0., 6., 0.1),
    'st':               frange(0., 1000., 100),
    'stmet':            frange(0, 1000, 100) + frange(1000, 4000, 500),
    'zpt':              frange(0., 2000., 100),
    'mtll':             frange(0, 200, 25) + frange(200, 500, 50) + frange(500, 1250, 250),
    'mtllg':            frange(0, 200, 25) + frange(200, 500, 50) + frange(500, 1750, 250),
    'mtgmet':            frange(0, 200, 25) + frange(200, 500, 50) + frange(500, 1750, 250),
    'mtl1met':            frange(0, 200, 25) + frange(200, 500, 50) + frange(500, 1750, 250),
    'mtl2met':            frange(0, 200, 25) + frange(200, 500, 50) + frange(500, 1750, 250),
    'mtllmet':            frange(0, 200, 25) + frange(200, 500, 50) + frange(500, 1750, 250),
    'mtllgmet':            frange(0, 200, 25) + frange(200, 500, 50) + frange(500, 1750, 250),
    'mt2':            frange(0., 425., 25.),
    'mzg_exo':            frange(0, 100, 50) + frange(100, 200, 50) + frange(200, 500, 50),
    'gammaMotherID':            frange(0., 200., 1.),
    'genPhotonPT':            np.concatenate((np.arange(0, 100, 10), np.arange(100, 350, 50)), axis=0),
    'genPhotonPT_Veto':            np.concatenate((np.arange(0, 100, 10), np.arange(100, 350, 50)), axis=0),
    'PhotonPT_Veto':            np.concatenate((np.arange(0, 100, 10), np.arange(100, 350, 50)), axis=0),
    'genPhotonPT_NoVeto':            np.concatenate((np.arange(0, 100, 10), np.arange(100, 350, 50)), axis=0),
    'PhotonPT_NoVeto':            np.concatenate((np.arange(0, 100, 10), np.arange(100, 350, 50)), axis=0),
    'VetoCompare':            np.arange(0., 200., 1.),
    'DeltaPhiLLMet':            frange(0., 6., 0.3),
    'DeltaEtaLLMet':            frange(0., 6., 0.3),
    'DeltaRLLMet':            frange(0., 6., 0.3),
    'nElectrons': frange(0, 10, 1),
    'nMuons': frange(0, 10, 1),
    'mTL3Met': frange(0., 300., 10),
    'Fakes': frange(0, 5, 1)
}

binnings_ = binnings.copy()

if(os.path.exists("plots_CR_zz/factors/CRZZ.pkl")):
    pklZZ = pkl.load(open("plots_CR_zz/factors/CRZZ.pkl", "rb"))
    ZZsf = pklZZ["LL"]["m_ll"][0]
else:
    ZZsf = 1.
if(os.path.exists("plots_CR_dy/factors/CRDY.pkl")):
    pklDY = pkl.load(open("plots_CR_dy/factors/CRDY.pkl", "rb"))
    DYsf = pklDY["LL"]["eta1"][0]
else:
    DYsf = 1.
if(os.path.exists("plots_CR_wz/factors/CRWZ.pkl")):
    pklWZ = pkl.load(open("plots_CR_wz/factors/CRWZ.pkl", "rb"))
    WZsf = pklWZ["LL"]["eta1"][0]
else:
    WZsf = 1.
if(os.path.exists("plots_CR_tt/factors/CRTT.pkl")):
    pklTT = pkl.load(open("plots_CR_tt/factors/CRTT.pkl", "rb"))
    TTsf = pklTT["EM"]["eta1"][0]
else:
    TTsf = 1.


def main():
    # bla
    # (["CRTT"], ["eta1"], "tt", "EM")
    # (groups, variables, region, flavorComb)

    arZZ = []
    arZZE = []
    arDY = []
    arDYE = []
    arTT = []
    arTTE = []
    arWZ = []
    arWZE = []
    arNames = []
    # zz,zzE,dy,dyE,tt,ttE,wz,wzE
    dummy1 = doScaling("ZZ", "DY", "TT", "WZ")
    dummy2 = doScaling("ZZ", "DY", "WZ", "TT")
    dummy3 = doScaling("ZZ", "TT", "WZ", "DY")
    dummy4 = doScaling("ZZ", "TT", "DY", "WZ")
    dummy5 = doScaling("ZZ", "WZ", "TT", "DY")
    dummy6 = doScaling("ZZ", "WZ", "DY", "TT")

    dummy7 = doScaling("DY", "ZZ", "TT", "WZ")
    dummy8 = doScaling("DY", "ZZ", "WZ", "TT")
    dummy9 = doScaling("DY", "TT", "WZ", "ZZ")
    dummy10 = doScaling("DY", "TT", "ZZ", "WZ")
    dummy11 = doScaling("DY", "WZ", "TT", "ZZ")
    dummy12 = doScaling("DY", "WZ", "ZZ", "TT")

    dummy13 = doScaling("WZ", "ZZ", "TT", "DY")
    dummy14 = doScaling("WZ", "ZZ", "DY", "TT")
    dummy15 = doScaling("WZ", "TT", "DY", "ZZ")
    dummy16 = doScaling("WZ", "TT", "ZZ", "DY")
    dummy17 = doScaling("WZ", "DY", "TT", "ZZ")
    dummy18 = doScaling("WZ", "DY", "ZZ", "TT")

    dummy19 = doScaling("TT", "ZZ", "WZ", "DY")
    dummy20 = doScaling("TT", "ZZ", "DY", "WZ")
    dummy21 = doScaling("TT", "WZ", "DY", "ZZ")
    dummy22 = doScaling("TT", "WZ", "ZZ", "DY")
    dummy23 = doScaling("TT", "DY", "WZ", "ZZ")
    dummy24 = doScaling("TT", "DY", "ZZ", "WZ")

    arNames.append("ZZ DY TT WZ")
    arNames.append("ZZ DY WZ TT")
    arNames.append("ZZ TT WZ DY")
    arNames.append("ZZ TT DY WZ")
    arNames.append("ZZ WZ TT DY")
    arNames.append("ZZ WZ DY TT")

    arNames.append("DY ZZ TT WZ")
    arNames.append("DY ZZ WZ TT")
    arNames.append("DY TT WZ ZZ")
    arNames.append("DY TT ZZ WZ")
    arNames.append("DY WZ TT ZZ")
    arNames.append("DY WZ ZZ TT")

    arNames.append("WZ ZZ TT DY")
    arNames.append("WZ ZZ DY TT")
    arNames.append("WZ TT DY ZZ")
    arNames.append("WZ TT ZZ DY")
    arNames.append("WZ DY TT ZZ")
    arNames.append("WZ DY ZZ TT")

    arNames.append("TT ZZ WZ DY")
    arNames.append("TT ZZ DY WZ")
    arNames.append("TT WZ DY ZZ")
    arNames.append("TT WZ ZZ DY")
    arNames.append("TT DY WZ ZZ")
    arNames.append("TT DY ZZ WZ")

    for d in [dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7, dummy8, dummy9, dummy10, dummy11, dummy12, dummy13, dummy14, dummy15, dummy16, dummy17, dummy18, dummy19, dummy20, dummy21, dummy22, dummy23, dummy24]:
        arZZ.append(d[0])
        arZZE.append(d[1])
        arDY.append(d[2])
        arDYE.append(d[3])
        arTT.append(d[4])
        arTTE.append(d[5])
        arWZ.append(d[6])
        arWZE.append(d[7])

    doPlot(arTT, arTTE, "TT", arNames)
    doPlot(arWZ, arWZE, "WZ", arNames)
    doPlot(arDY, arDYE, "DY", arNames)
    doPlot(arZZ, arZZE, "ZZ", arNames)

    # find outliers
    nomZZ = arZZ[0]
    nomZZE = arZZE[0]
    nomWZ = arWZ[0]
    nomWZE = arWZE[0]
    nomTT = arTT[0]
    nomTTE = arTTE[0]
    nomDY = arDY[0]
    nomDYE = arDYE[0]

    cloneZZ = arZZ
    cloneZZ = [abs(cloneZZ[i] - nomZZ) for i in range(len(cloneZZ))]
    outZZ = max(cloneZZ)
    cloneWZ = arWZ
    cloneWZ = [abs(cloneWZ[i] - nomWZ) for i in range(len(cloneWZ))]
    outWZ = max(cloneWZ)
    cloneDY = arDY
    cloneDY = [abs(cloneDY[i] - nomDY) for i in range(len(cloneDY))]
    outDY = max(cloneDY)
    cloneTT = arTT
    cloneTT = [abs(cloneTT[i] - nomTT) for i in range(len(cloneTT))]
    outTT = max(cloneTT)

    print "nominal"
    print "ZZ DY TT WZ"
    print "ZZ", nomZZ, nomZZE, "=", str(
        nomZZE / nomZZ * 100) + "%", "-", outZZ, "=", str(outZZ / nomZZ * 100) + "%"
    print "WZ", nomWZ, nomWZE, "=", str(
        nomWZE / nomWZ * 100) + "%", "-", outWZ, "=", str(outWZ / nomWZ * 100) + "%"
    print "DY", nomDY, nomDYE, "=", str(
        nomDYE / nomDY * 100) + "%", "-", outDY, "=", str(outDY / nomDY * 100) + "%"
    print "TT", nomTT, nomTTE, "=", str(
        nomTTE / nomTT * 100) + "%", "-", outTT, "=", str(outTT / nomTT * 100) + "%"

    toSave = {}
    toSave["ZZ"] = outZZ
    toSave["DY"] = outDY
    toSave["TT"] = outTT
    toSave["WZ"] = outWZ
    pkl.dump(toSave, open("plots_Order/orderUnc.pkl", "wb"))


def doPlot(ar, arE, name, arN):
    c = TCanvas("canvas", "", 800, 1000)
    x = [ar[0]]
    ex = [arE[0]]
    y = [len(ar) / 2.]
    eyU = [round(len(ar) / 2.) + 1]
    eyD = [len(ar) / 2.]
    ge = TGraphAsymmErrors()
    for i in range(len(x)):
        ge.SetPoint(i, x[i], y[i])
        ge.SetPointError(i, ex[i], ex[i], eyD[i], eyU[i])
    ge.SetFillStyle(3003)
    ge.SetFillColor(ROOT.kBlue)
    ge.Draw("same a2")

    gr = TGraphAsymmErrors()

    for i in range(len(ar)):
        y = float(i + 1)
        gr.SetPoint(i, ar[i], y)
        gr.SetPointError(i, arE[i], arE[i], 0., 0.)

    gr.SetTitle("; Scale Factor #alpha; ")
    ge.SetTitle("; Scale Factor #alpha; ")
    gr.SetMarkerColor(kRed)
    gr.SetMarkerStyle(21)
    gr.Draw("sameP")
    line = TLine(ar[0], 0., ar[0], len(ar) + 1.)
    line.SetLineColor(kBlue)
    line.SetLineWidth(2)
    line.Draw("SAME")

    # ge.GetXaxis().SetLimits(0.75, 0.85)
    ge.GetYaxis().SetLabelOffset(1)
    c.SetGridy()
    c.Update()

    t = TLatex()
    t.SetTextAlign(32)
    t.SetTextSize(0.02)
    # t.SetTextFont(72)
    for i in range(len(arN)):
        if name == "TT":
            t.DrawLatex(0.765, (i + 1), arN[i])
        if name == "WZ":
            t.DrawLatex(1.09, (i + 1), arN[i])
        if name == "DY":
            t.DrawLatex(1.0545, (i + 1), arN[i])
        if name == "ZZ":
            t.DrawLatex(1.03, (i + 1), arN[i])

    l = ROOT.TLatex(
        0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Private Work}}")
    l.SetNDC()
    l.Draw()

    lum = ROOT.TLatex(.62, .95, "%.1f fb^{-1} (%s TeV)" %
                      (aux.intLumi / 1000., aux.Label.cmsEnergy))
    lum.SetNDC()
    lum.Draw()

    # c.SaveAs('plots_Order/TT_Compare.pdf')
    c.SaveAs('plots_Order/' + name + '_Compare.pdf')


def doScaling(first, second, third, fourth):
    # remove all pkl first

    # if(os.path.exists("plots_CR_zz/factors/CRZZ.pkl")):
    #     os.remove("plots_CR_zz/factors/CRZZ.pkl")
    # if(os.path.exists("plots_CR_dy/factors/CRDY.pkl")):
    #     os.remove("plots_CR_dy/factors/CRDY.pkl")
    # if(os.path.exists("plots_CR_wz/factors/CRWZ.pkl")):
    #     os.remove("plots_CR_wz/factors/CRWZ.pkl")
    # if(os.path.exists("plots_CR_tt/factors/CRTT.pkl")):
    #     os.remove("plots_CR_tt/factors/CRTT.pkl")

    SF_ZZ, SF_ZZErr, SF_DY, SF_DYErr, SF_TT, SF_TTErr, SF_WZ, SF_WZErr = [
        1., 0., 1., 0., 1., 0., 1., 0.]

    for i, rank in enumerate([first, second, third, fourth], 1):
        if rank == "ZZ":
            SF_ZZ, SF_ZZErr = drawCRZZ("LL", "CRZZ/LL/nom/eta1", dataLL, binning=binnings_[
                                       "eta1"], xTitle=labels["eta1"][0], weightsToUse=["nISR", "topPt", "ewk"], SF_TT=SF_TT, SF_WZ=SF_WZ, SF_DY=SF_DY, SF_ZZ=SF_ZZ)
        elif rank == "DY":
            SF_DY, SF_DYErr = drawDYCR("LL", "CRDY" + "/" + "LL" + "/nom/" + "eta1", dataLL, binning=binnings_[
                "eta1"], xTitle=labels["eta1"][0], weightsToUse=["nISR", "topPt", "ewk"], SF_TT=SF_TT, SF_WZ=SF_WZ, SF_DY=SF_DY, SF_ZZ=SF_ZZ)
        elif rank == "TT":
            SF_TT, SF_TTErr = drawTTCR("EM", "CRTT" + "/" + "EM" + "/nom/" + "eta1", binning=binnings_[
                "eta1"], xTitle=labels["eta1"][0], weightsToUse=["nISR", "topPt", "ewk"], SF_TT=SF_TT, SF_WZ=SF_WZ, SF_DY=SF_DY, SF_ZZ=SF_ZZ)
        elif rank == "WZ":
            SF_WZ, SF_WZErr = drawCRWZ("LL", "CRWZ" + "/" + "LL" + "/nom/" + "eta1", dataLL, binning=binnings_[
                "eta1"], xTitle=labels["eta1"][0], weightsToUse=["nISR", "topPt", "ewk"], SF_TT=SF_TT, SF_WZ=SF_WZ, SF_DY=SF_DY, SF_ZZ=SF_ZZ)
        else:
            if i == 1:
                print "ERROR: First not found!", first
            if i == 2:
                print "ERROR: Second not found!", second
            if i == 3:
                print "ERROR: Third not found!", third
            if i == 4:
                print "ERROR: Fourth not found!", fourth

    # print "........................"
    # print first, second, third, fourth
    # print "ZZ", SF_ZZ, SF_ZZErr
    # print "DY", SF_DY, SF_DYErr
    # print "TT", SF_TT, SF_TTErr
    # print "WZ", SF_WZ, SF_WZErr
    # print ""

    return SF_ZZ, SF_ZZErr, SF_DY, SF_DYErr, SF_TT, SF_TTErr, SF_WZ, SF_WZErr

    # ZZ,DY,TT,WZ
    # ZZ,DY,WZ,TT
    # ZZ,TT,DY,WZ
    # ZZ,TT,Wz,DY
    # ZZ,WZ,TT,DY
    # ZZ,WZ,DY,TT


if __name__ == "__main__":
    main()
