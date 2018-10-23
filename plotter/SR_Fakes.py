from dataMC import labels, frange

import ROOT
from ROOT import *
from array import array
from include import *
import numpy as np
import pickle as pkl
import os

binnings = {
    # 'pt1':              frange(20,100,10)+frange(100,155,25),
    # 'pt1':              frange(20,100,10)+frange(100,155,10),
    # 'pt1':              frange(25,151,25),
    'pt1':              frange(20, 140, 20),
    'pt2':              frange(0, 151, 25),
    # 'pt2':              frange(0, 300, 10),
    'eta1':             frange(0., 2.7, 0.52),
    'eta2':             frange(0., 2.7, 0.52),
    'phi1':             frange(0., 3.51, 0.7),
    'phi2':             frange(0., 3.51, 0.7),
    'ht':               frange(0., 2000., 250),
    # 'met':               [0,25,50,75,100,150,190,230,500],
    # 'met':               frange(100,155,10),
    'met':               frange(100, 280, 20),
    # 'met':               frange(100,200,10)+frange(220,300,20),
    'm_ll':             frange(80., 111., 10),
    'm_llg':             frange(0, 100, 10) + frange(100, 200, 10) + frange(200, 500, 50),
    'pt_llg':             frange(0, 100, 10) + frange(100, 200, 10) + frange(200, 500, 50),
    'n_jets':           frange(0., 7., 1),
    'n_photons':         frange(0., 4., 1),
    'n_vtx':            frange(0., 40., 1),
    # 'pt_g1':            frange(20,100,10)+frange(100,210,50),
    'pt_g1':            frange(25, 126, 25),
    'eta_g1':           frange(0., 2.60, 0.1),
    'phi_g1':           frange(0., 3.50, 0.1),
    'sigmaIetaIeta_g1': frange(0., 0.04, 0.01),
    'sigmaIphiIphi_g1': frange(0., 0.2, 0.01),
    'deltaR1_g1':       frange(0., 6., 0.1),
    'deltaR2_g1':       frange(0., 6., 0.1),
    'r9_g1':            frange(0., 1.5, 0.05),
    'hOverE_g1':        frange(0., 0.1, 0.01),
    'deltaEtaLL':       frange(0., 6., 0.1),
    'deltaPhiLL':       frange(0., 6., 0.1),
    'deltaEtaLLG':      frange(0., 6., 0.1),
    'deltaPhiLLG':      frange(0., 6., 0.01),
    'deltaRLL':         frange(0., 6., 0.12),
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


def calculateSFAndError(numerator_data, denominator_toScale, additional_fix):
    num_dataErr = (ROOT.Double(0))
    num_data = numerator_data.IntegralAndError(0, -1, num_dataErr)

    den_toScaleErr = (ROOT.Double(0))
    den_toScale = denominator_toScale.IntegralAndError(0, -1, den_toScaleErr)

    add_fixError = ROOT.Double(0)
    add_fix = additional_fix.IntegralAndError(0, -1, add_fixError)

    alpha = (num_data - add_fix) / (den_toScale)
    alphaErr = np.sqrt((num_dataErr / den_toScale)**2. + (add_fixError / den_toScale)
                       ** 2. + (den_toScaleErr * (num_data - add_fix) / (den_toScale)**2.)**2.)

    return [alpha, alphaErr / alpha] if den_toScale else [1., 0.]


def drawVR(sampleNames, name, binning=None, binningName="", xTitle=None, yTitle=None, sfZZ=[], sfDY=[], sfTT=[], sfWZ=[]):
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()

    #folder= (name.split("/"))[0]
    #folder= (name.split("/"))[0]+"/"+(name.split("/"))[1]+"/"+(name.split("/"))[2]
    folder = (name.split("/"))[0] + "/" + (name.split("/"))[1]

    style.divideByBinWidth = False

    style.minimumOne = True
    #print name
    #dataHist = aux.stdHist(dataMuonEG, name, binning)
    dataHist = aux.stdHist(dataLL, name, binning)
    #dataHist = aux.stdHist(dataDoubleSF, name, binning)
    aux.drawOpt(dataHist, "data")

    #zgHist = aux.stdHist(zgamma, name, binning)
    #ttgHist = aux.stdHist(ttgamma, name, binning)
    #ttg080Hist = aux.stdHist(ttgamma, name.replace("VR","VR080"), binning)
    #ttg80Hist = aux.stdHist(ttgamma, name.replace("VR","VR80"), binning)
    #zzHist = aux.stdHist(zz, name, binning)
    #wwgHist = aux.stdHist(wwgamma, name, binning)
    #wzgHist = aux.stdHist(wzgamma, name, binning)
    #dyHist = aux.stdHist(DYjetsNLO, name, binning)
    #wjetsHist = aux.stdHist(wjets, name, binning)
    #ttHist = aux.stdHist(tt, name, binning)
    #tt080Hist = aux.stdHist(tt, name.replace("VR","VR080"), binning)
    #tt80Hist = aux.stdHist(tt, name.replace("VR","VR80"), binning)
    #singletopHist=aux.stdHist(singletop, name, binning)
    #wzHist=aux.stdHist(wz, name, binning)
    #wwHist=aux.stdHist(ww, name, binning)
    #zz4lHist=aux.stdHist(zz4l, name, binning)
    #wgHist=aux.stdHist(wgamma, name, binning)
    zgHist = aux.stdHistWithWeights(
        zgamma, name, ["nISR", "topPt", "ewk"], binning)
    ttgHist = aux.stdHistWithWeights(
        ttgamma, name, ["nISR", "topPt", "ewk"], binning)
    ttg080Hist = aux.stdHistWithWeights(ttgamma, name.replace(
        "VR", "VR080"), ["nISR", "topPt", "ewk"], binning)
    ttg80Hist = aux.stdHistWithWeights(ttgamma, name.replace(
        "VR", "VR80"), ["nISR", "topPt", "ewk"], binning)
    zzHist = aux.stdHistWithWeights(
        zz, name, ["nISR", "topPt", "ewk"], binning)
    wwgHist = aux.stdHistWithWeights(
        wwgamma, name, ["nISR", "topPt", "ewk"], binning)
    wzgHist = aux.stdHistWithWeights(
        wzgamma, name, ["nISR", "topPt", "ewk"], binning)
    dyHist = aux.stdHistWithWeights(
        DYjetsNLO, name, ["nISR", "topPt", "ewk"], binning)
    wjetsHist = aux.stdHistWithWeights(
        wjets, name, ["nISR", "topPt", "ewk"], binning)
    ttHist = aux.stdHistWithWeights(
        tt, name, ["nISR", "topPt", "ewk"], binning)
    tt080Hist = aux.stdHistWithWeights(tt, name.replace("VR", "VR080"), [
                                       "nISR", "topPt", "ewk"], binning)
    tt80Hist = aux.stdHistWithWeights(tt, name.replace("VR", "VR80"), [
                                      "nISR", "topPt", "ewk"], binning)
    singletopHist = aux.stdHistWithWeights(
        singletop, name, ["nISR", "topPt", "ewk"], binning)
    wzHist = aux.stdHistWithWeights(
        wz, name, ["nISR", "topPt", "ewk"], binning)
    wwHist = aux.stdHistWithWeights(
        ww, name, ["nISR", "topPt", "ewk"], binning)
    zz4lHist = aux.stdHistWithWeights(
        zz4l, name, ["nISR", "topPt", "ewk"], binning)
    wgHist = aux.stdHistWithWeights(
        wgamma, name, ["nISR", "topPt", "ewk"], binning)

    final_dyHist = aux.addHists(dyHist, zgHist)
    final_ttHist = aux.addHists(ttHist, ttgHist)
    final_tt080Hist = aux.addHists(tt080Hist, ttg080Hist)
    final_tt80Hist = aux.addHists(tt80Hist, ttg80Hist)
    final_zzHist = aux.addHists(zzHist, zz4lHist)
    final_wzHist = aux.addHists(wzHist)
    final_otherHist = aux.addHists(
        wwgHist, wzgHist, wjetsHist, singletopHist, wwHist, wgHist)

    final_tt080Hist.Scale(sfTT[0])
    final_tt80Hist.Scale(sfTT[0])
    final_ttHist.Scale(sfTT[0])
    final_zzHist.Scale(sfZZ[0])
    final_wzHist.Scale(sfWZ[0])
    final_dyHist.Scale(sfDY[0])

    # final_otherHist.SetLineColor(ROOT.kCyan+3)
    # final_ttHist.SetLineColor(ROOT.kOrange+8)
    # final_ttHist.SetLineColor(ROOT.kRed+1)
    final_ttHist.SetLineColor(ROOT.kOrange + 7)
    final_zzHist.SetLineColor(ROOT.kOrange - 2)
    # final_wzHist.SetLineColor(ROOT.kAzure-5)
    final_wzHist.SetLineColor(ROOT.kAzure - 6)
    # final_wzHist.SetLineColor(ROOT.kBlue+3)
    # final_dyHist.SetLineColor(ROOT.kGreen-3)
    final_dyHist.SetLineColor(ROOT.kGreen + 3)
    # final_dyHist.SetLineColor(ROOT.kGreen+2)
    final_otherHist.SetLineColor(ROOT.kGray + 2)

    zgHist.SetLineColor(ROOT.kGreen - 3)
    ttgHist.SetLineColor(ROOT.kRed + 1)
    zzHist.SetLineColor(ROOT.kYellow)
    wwgHist.SetLineColor(ROOT.kCyan + 3)
    wzgHist.SetLineColor(ROOT.kBlue + 3)
    dyHist.SetLineColor(ROOT.kGreen + 3)
    wjetsHist.SetLineColor(ROOT.kBlue - 9)
    ttHist.SetLineColor(ROOT.kOrange + 8)
    singletopHist.SetLineColor(ROOT.kOrange + 4)
    # wzHist.SetLineColor(ROOT.kAzure-5)
    wzHist.SetLineColor(ROOT.kAzure - 6)
    # wzHist.SetLineColor(ROOT.kBlue-3)
    wwHist.SetLineColor(ROOT.kCyan - 2)
    zz4lHist.SetLineColor(ROOT.kOrange - 2)
    wgHist.SetLineColor(ROOT.kRed + 3)

    #print sfTT

    ttHist.Scale(sfTT[0])
    ttgHist.Scale(sfTT[0])
    tt080Hist.Scale(sfTT[0])
    ttg080Hist.Scale(sfTT[0])
    tt80Hist.Scale(sfTT[0])
    ttg80Hist.Scale(sfTT[0])
    zgHist.Scale(sfDY[0])
    dyHist.Scale(sfDY[0])
    zzHist.Scale(sfZZ[0])
    zz4lHist.Scale(sfZZ[0])
    wzHist.Scale(sfWZ[0])

    if "Fakes" in name:
        normFactor = 1. / aux.addHists(final_ttHist, final_otherHist,
                                       final_dyHist, final_zzHist, wzHist).Integral(1, 1)
        final_ttHist.Scale(normFactor)
        final_otherHist.Scale(normFactor)
        final_dyHist.Scale(normFactor)
        wzHist.Scale(normFactor)
        final_zzHist.Scale(normFactor)
        dataHist.Scale(normFactor)

    mcSystUncert = 0.0  # SF, lumi, trigger
    #zgSyst = aux.getSysHisto(zgHist, sfDY[1])
    #ttgSyst = aux.getSysHisto(ttgHist, sfTT[1])
    #zzSyst = aux.getSysHisto(zzHist, sfZZ[1])
    #wwgSyst = aux.getSysHisto(wwgHist, mcSystUncert)
    #wzgSyst = aux.getSysHisto(wzgHist, mcSystUncert)
    #dySyst = aux.getSysHisto(dyHist, sfDY[1])
    #wjetsSyst = aux.getSysHisto(wjetsHist, mcSystUncert)
    #ttSyst = aux.getSysHisto(ttHist, sfTT[1])
    # singletopSyst=aux.getSysHisto(singletopHist,mcSystUncert)
    # wzSyst=aux.getSysHisto(wzHist,sfWZ[1])
    # wwSyst=aux.getSysHistoW(wwHist,mcSystUncert)
    # zz4lSyst=aux.getSysHisto(zz4lHist,sfZZ[1])
    # wgSyst=aux.getSysHisto(wgHist,mcSystUncert)

    # ttg080Syst=aux.getSysHisto(ttg080Hist,0.04)
    # ttg80Syst=aux.getSysHisto(ttg80Hist,0.4)
    # ttg80Syst=aux.getSysHisto(ttg80Hist,0.2)
    # tt080Syst=aux.getSysHisto(tt080Hist,0.04)
    # tt80Syst=aux.getSysHisto(tt80Hist,0.4)
    # tt80Syst=aux.getSysHisto(tt80Hist,0.2)

    otherSyst = aux.getSysHistoWithMeanWeight(
        final_otherHist, mcSystUncert, final_otherHist.Integral() / sum(otherBKG.ngens))

    wzSFSyst = aux.getSysHisto(wzHist, sfWZ[1])
    final_tt080SFSyst = aux.getSysHisto(final_tt080Hist, sfTT[1])
    final_tt80SFSyst = aux.getSysHisto(final_tt80Hist, 0.2)
    final_zzSFSyst = aux.getSysHisto(final_zzHist, sfZZ[1])
    final_dySFSyst = aux.getSysHistoWithMeanWeight(
        final_dyHist, sfDY[1], final_dyHist.Integral() / (sum(DYjetsNLO.ngens) + sum(zgamma.ngens)))

    final_ttTotSFSyst = aux.addHists(final_tt080SFSyst, final_tt80SFSyst)

    totStat = aux.addHists(final_zzHist, final_dyHist,
                           final_ttHist, wzHist, final_otherHist)
    totSyst = aux.addHists(final_dySFSyst, final_zzSFSyst,
                           final_ttTotSFSyst, wzSFSyst, otherSyst)

    # if "pt_g1" in name:
    # cutBin=ttHist.FindBin(80)
    # scaleSyst=aux.getSysHistoCut(scaleHist,0.04,0.4,cutBin)
    # ttgSyst=aux.getSysHistoCut(ttgHist,0.04,0.4,cutBin)
    # ttSyst=aux.getSysHistoCut(ttHist,0.04,0.4,cutBin)
    # else:
    # scaleSyst=aux.getSysHisto(scaleHist,sferr)
    # ttgSyst=aux.getSysHisto(ttgHist,0.04)
    # ttSyst=aux.getSysHisto(ttHist,0.04)

    #dataSyst= aux.getSysHisto(dataHist,mcSystUncert)

    #totStat = aux.addHists(zgHist, ttgHist, zzHist, wwgHist, wzgHist, dyHist, wjetsHist, ttHist, singletopHist, wzHist, wwHist, zz4lHist, wgHist)
    #totSyst = aux.addHists(zgSyst, ttgSyst, zzSyst, wwgSyst, wzgSyst, dySyst, wjetsSyst, ttSyst, singletopSyst, wzSyst, wwSyst, zz4lSyst, wgSyst)
    #totSyst = aux.addHists(zgSyst, ttg080Syst,ttg80Syst, zzSyst, wwgSyst, wzgSyst, dySyst, wjetsSyst, tt080Syst,tt80Syst, singletopSyst, wzSyst, wwSyst, zz4lSyst, wgSyst)
    #totSyst = aux.addHists(zgSyst, ttg080Syst,ttg80Syst, zzSyst, dySyst, tt080Syst,tt80Syst, wzSyst, zz4lSyst)

    totUnc = aux.addHistUncert(totStat, totSyst)

    # for bin in range(ttg80Syst.GetNbinsX()+1):
    #print "ttg80",ttg80Syst.GetBinLowEdge(bin),ttg80Syst.GetBinContent(bin),ttg80Syst.GetBinError(bin)
    #print "tt80",tt80Syst.GetBinLowEdge(bin),tt80Syst.GetBinContent(bin),tt80Syst.GetBinError(bin)
    #print "ttg080",ttg080Syst.GetBinLowEdge(bin),ttg080Syst.GetBinContent(bin),ttg080Syst.GetBinError(bin)
    #print "tt080",tt080Syst.GetBinLowEdge(bin),tt080Syst.GetBinContent(bin),tt080Syst.GetBinError(bin)
    #print "tot",totSyst.GetBinLowEdge(bin),totSyst.GetBinContent(bin),totSyst.GetBinError(bin)

    aux.drawOpt(totUnc, "totUnc")
    aux.drawOpt(totSyst, "sysUnc")

    #print "mc:",totStat.Integral(),"data",dataHist.Integral()
    #print "zg", zgHist.Integral(),zgHist.GetEntries()
    #print "ttg", ttgHist.Integral(),ttgHist.GetEntries()
    #print "zz", zzHist.Integral(),zzHist.GetEntries()
    #print "wwg", wwgHist.Integral(),wwgHist.GetEntries()
    #print "wzg", wzgHist.Integral(),wzgHist.GetEntries()
    #print "dy", dyHist.Integral(),dyHist.GetEntries()
    #print "wjets", wjetsHist.Integral(),wjetsHist.GetEntries()
    #print "tt", ttHist.Integral(),ttHist.GetEntries()
    #print "singletop", singletopHist.Integral(),singletopHist.GetEntries()
    #print "wz", wzHist.Integral(),wzHist.GetEntries()
    #print "ww", wwHist.Integral(),wwHist.GetEntries()
    #print "zz4l", zz4lHist.Integral(),zz4lHist.GetEntries()
    #print "wg", wgHist.Integral(),wgHist.GetEntries()

    #signal2 = aux.stdHistWithWeights(t5bbbbzg_1500_400, name, ["nISR","topPt","ewk"],binning)
    #signal1 = aux.stdHistWithWeights(tching_400, name,["nISR","topPt","ewk"], binning)
    signal2 = aux.stdHistWithWeights(
        gmsb_240_230, name, ["nISR", "topPt", "ewk"], binning)
    signal1 = aux.stdHistWithWeights(
        gmsb_290_205, name, ["nISR", "topPt", "ewk"], binning)
    #signal1 = aux.stdHist(tching_400, name, binning)

    #signal2_AvgTopPtWeightHisto = t5bbbbzg_1500_400.getHist(folder+"/weight_topPt")
    #signal2_AvgNIsrWeightHisto = t5bbbbzg_1500_400.getHist(folder+"/weight_nISR")
    #signal2_AvgEWKinoWeightHisto = t5bbbbzg_1500_400.getHist(folder+"/weight_EWKinoPairPt")
    #signal1_AvgTopPtWeightHisto = tching_600.getHist(folder+"/weight_topPt")
    #signal1_AvgNIsrWeightHisto = tching_600.getHist(folder+"/weight_nISR")
    #signal1_AvgEWKinoWeightHisto = tching_600.getHist(folder+"/weight_EWKinoPairPt")

    #signal2.Scale((1./(signal2_AvgTopPtWeightHisto.GetMean() if signal2_AvgTopPtWeightHisto.GetMean() >0. else 1.)))
    #signal2.Scale((1./(signal2_AvgNIsrWeightHisto.GetMean() if signal2_AvgNIsrWeightHisto.GetMean() >0. else 1.)))
    #signal2.Scale((1./(signal2_AvgEWKinoWeightHisto.GetMean() if signal2_AvgEWKinoWeightHisto.GetMean() >0. else 1.)))
    #signal1.Scale((1./(signal1_AvgTopPtWeightHisto.GetMean() if signal1_AvgTopPtWeightHisto.GetMean() >0. else 1.)))
    #signal1.Scale((1./(signal1_AvgNIsrWeightHisto.GetMean() if signal1_AvgNIsrWeightHisto.GetMean() >0. else 1.)))
    #signal1.Scale((1./(signal1_AvgEWKinoWeightHisto.GetMean() if signal1_AvgEWKinoWeightHisto.GetMean() >0. else 1.)))

    for h in signal1, signal2:
        aux.drawOpt(h, "signal")
    signal1.SetLineColor(ROOT.kBlue + 3)
    signal2.SetLineColor(ROOT.kBlue + 3)
    signal2.SetLineStyle(2)

    #m.add(signal2, "T5bbbbZg")
    #m.add(signal1, "TChiNg")
    if not "Fakes" in name:
        m.add(signal2, "gmsb_240_230")
        m.add(signal1, "gmsb_290_205")

    if "Fakes" in name:
        labels = ["all", "genPhoton",
                  "genElectron", "genJet", "noMatch"]
        for i in range(len(labels)):
            final_ttHist.GetXaxis().SetBinLabel(i + 1, labels[i])
            final_ttHist.GetXaxis().LabelsOption("d")
            final_ttHist.SetLabelOffset(10.5, "X")
            final_ttHist.GetXaxis().SetLabelOffset(10.5)
            final_ttHist.LabelsOption("d", "X")
            final_zzHist.GetXaxis().SetBinLabel(i + 1, labels[i])
            final_zzHist.GetXaxis().LabelsOption("d")
            final_zzHist.LabelsOption("d", "X")
            final_zzHist.SetLabelOffset(10.5, "X")
            final_zzHist.GetXaxis().SetLabelOffset(10.5)
            final_otherHist.GetXaxis().SetBinLabel(i + 1, labels[i])
            final_otherHist.GetXaxis().LabelsOption("d")
            final_otherHist.LabelsOption("d", "X")
            final_otherHist.SetLabelOffset(10.5, "X")
            final_otherHist.GetXaxis().SetLabelOffset(10.5)
            dataHist.GetXaxis().SetBinLabel(i + 1, labels[i])
            dataHist.GetXaxis().LabelsOption("d")
            dataHist.LabelsOption("d", "X")
            dataHist.GetXaxis().SetLabelOffset(10.5)
            final_ttHist.SetYTitle("arbitrary units")
            final_zzHist.SetYTitle("arbitrary units")
            final_otherHist.SetYTitle("arbitrary units")
            dataHist.SetYTitle("arbitrary units")

    # m.addStack(zgHist, "Z#gamma")
    # m.addStack(ttgHist, "t#bar{t}#gamma")
    #m.addStack(zzHist, "ZZ")
    # m.addStack(wwgHist, "WW#gamma")
    # m.addStack(wzgHist, "WZ#gamma")
    #m.addStack(dyHist, "Drell-Yan/Z")
    #m.addStack(wjetsHist, "W+jets")
    # m.addStack(ttHist, "t#bar{t}")
    #m.addStack(singletopHist, "single t")
    #m.addStack(wwHist, "WW")
    m.addStack(wzHist, "WZ")
    #m.addStack(zz4lHist, "ZZ(4l)")
    # m.addStack(wgHist, "W#gamma")

    # m.addStack(aux.addHists(wwgHist,wzgHist,wjetsHist,singletopHist,wwHist,wgHist),"other")
    # m.addStack(aux.addHists(ttHist,ttgHist),"t#bar{t}#gamma(+#gamma)")
    # m.addStack(aux.addHists(tt080Hist,tt80Hist,ttg080Hist,ttg80Hist),"t#bar{t}#gamma(+#gamma)")
    # m.addStack(aux.addHists(zgHist,dyHist),"Drell-Yan/Z(+#gamma)")
    # m.addStack(aux.addHists(zz4lHist,zzHist),"ZZ")

    m.addStack(final_otherHist, "other")
    m.addStack(final_ttHist, "t#bar{t}(+#gamma)")
    # m.addStack(aux.addHists(tt080Hist,tt80Hist,ttg080Hist,ttg80Hist),"t#bar{t}#gamma(+#gamma)")
    m.addStack(final_dyHist, "Drell-Yan/Z(+#gamma)")
    m.addStack(final_zzHist, "ZZ")

    # for bin in range(signal1.GetNbinsX()+1):
    #print signal1.GetBinLowEdge(bin),signal1.GetBinContent(bin),totStat.GetBinContent(bin)
    #print signal2.GetBinLowEdge(bin),signal2.GetBinContent(bin),totStat.GetBinContent(bin)
    #print signal1.Integral(),totStat.Integral()
    #print signal2.Integral(),totStat.Integral()

    m.sortStackByIntegral()

    if not "Fakes" in name:
        m.add(dataHist, "Data")
        m.add(totUnc, "Total uncertainty")
        m.add(totSyst, "syst. uncertainty")

    #legInfo = "Validation Region"
    # m.leg.SetHeader(legInfo)
    # m.leg.SetY1(.56)
    # m.leg.SetX1(.56)
    # m.leg.SetX1(.46)
    # m.leg.SetX2(.99)
    # m.leg.SetX2(.89)

    m.Draw()
    if not "Fakes" in name:
        KS = totStat.Clone().KolmogorovTest(dataHist.Clone(), "UO")
        ksText = ROOT.TLatex()
        ksText.SetTextSize(0.45 * ksText.GetTextSize())
        ksText.DrawLatexNDC(0.57, 0.65, "KS Value= " + str(np.round(KS, 5)))

        r = ratio.Ratio(
            "#scale[.9]{#lower[.24]{#splitline{Data/Pred.}{Bkg. frac.}}}", dataHist, totStat, sysHisto=totSyst)
        #rMax = 1.5
        rMax = 2.
        #r.draw(0., rMax, m.getStack(), True)
        r.draw(0., rMax, m.getStack())
    else:
        r = ratio.Ratio(
            "#scale[.9]{#lower[.24]{#splitline{Data/Pred.}{Bkg. frac.}}}", dataHist, totStat)
        #rMax = 1.5
        rMax = 2.
        #r.draw(0., rMax, m.getStack(), True)
        r.draw(0., rMax, m.getStack())

    aux.Label(sim=False, status="Work in Progress")
    #aux.save(name, normal=False, changeMinMax=False)
    directory = "plots/"
    if not os.path.exists(directory):
        os.makedirs(directory)
    if not "Fakes" in name:
        aux.save(name.replace("/", "_"), folder=directory)
    else:
        aux.save(name.replace("/", "_"), folder=directory, changeMinMax=False)
    # aux.save(name.replace("/","_"),folder=directory)

    # return sf,sferr


# variables=["mt2","eta1","eta2","pt1","pt2","n_jets","n_vtx","phi1","phi2","m_ll","ht","n_photons","pt_g1","met","nElectrons","nMuons","deltaR1_g1","deltaR2_g1","deltaRLL","Fakes"]
variables = ["Fakes"]
# variables=["pt_g1"]
bkgs = [DYjetsNLO, zgamma, tt, ttgamma, wwgamma, wzgamma,
        zz, wjets, wgamma, singletop, zz4l, wz, zz, ww]
# groups=["ValidationRegion"]
groups = ["onZMet150"]
binnings_ = binnings.copy()
# binnings_["pt_g1"]=frange(25,80,5)+frange(80,140,10)+frange(140,200,20)
# 'pt_g1':            frange(25,126,25),
# binnings_["pt_g1"]=frange(20,60,10)+frange(60,80,20)+frange(80,160,40)
# binnings_["pt_g1"]=frange(20,150,10)
binnings_["pt_g1"] = frange(20, 100, 10)
binnings_["mt2"] = frange(80, 200, 20)
#pklZZ = pkl.load( open( "plots_CR_zz/factors/ControlRegionZZ.pkl", "rb" ) )
#pklDY = pkl.load( open( "plots_CR_dy/factors/ControlRegionDY.pkl", "rb" ) )
#pklTT = pkl.load( open( "plots_CR_tt/factors/ControlRegionTT.pkl", "rb" ) )
#pklWZ = pkl.load( open( "plots_CR_wz/factors/ControlRegionWZ.pkl", "rb" ) )
pklZZ = pkl.load(open("plots_CR_zz/factors/CRZZ.pkl", "rb"))
pklDY = pkl.load(open("plots_CR_dy/factors/CRDY.pkl", "rb"))
pklTT = pkl.load(open("plots_CR_tt/factors/CRTT.pkl", "rb"))
pklWZ = pkl.load(open("plots_CR_wz/factors/CRWZ.pkl", "rb"))
ZZsf = pklZZ["LL"]["m_ll"]
DYsf = pklDY["LL"]["eta1"]
TTsf = pklTT["EM"]["eta1"]
WZsf = pklWZ["LL"]["eta1"]
for group in groups:
    for variable in variables:
        # drawVR("EE",group+"/EE/"+variable,binning=binnings_[variable],xTitle=labels[variable][0],sfZZ=ZZsf,sfDY=DYsf,sfTT=TTsf,sfWZ=WZsf)
        # drawVR("MM",group+"/MM/"+variable,binning=binnings_[variable],xTitle=labels[variable][0],sfZZ=ZZsf,sfDY=DYsf,sfTT=TTsf,sfWZ=WZsf)
        drawVR("LL", group + "/LL/" + variable, binning=binnings_[
               variable], xTitle=labels[variable][0], sfZZ=ZZsf, sfDY=DYsf, sfTT=TTsf, sfWZ=WZsf)
        # drawVR("EM",group+"/EM/"+variable,binning=binnings_[variable],xTitle=labels[variable][0],sfZZ=ZZsf,sfDY=DYsf,sfTT=TTsf,sfWZ=WZsf)
