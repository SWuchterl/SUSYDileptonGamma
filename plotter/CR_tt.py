from dataMC import labels, frange

import ROOT
from ROOT import *
from array import array
from include import *
import numpy as np
import pickle as pkl
import os

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

# pklTT = pkl.load(open("plots_CR_tt/factors/CRTT.pkl", "rb"))
# TTsf = pklTT["EM"]["eta1"]


binnings = {
    # 'pt1':              frange(20,100,10)+frange(100,200,25)+range(200,350,50),
    # 'pt1':              frange(25,160,5),
    # 'pt1':              frange(20,150,10),
    # 'pt1':              frange(20,150,10),
    # 'pt1':              frange(20.,200.,10.),
    # 'pt1':              frange(20.,200.,5.),
    # 'pt1':              frange(20.,220.,5.)+frange(230,250,10),
    # 'pt1':              frange(20., 250., 5.),
    'pt1':              frange(20., 300., 5.),
    # 'pt1':              frange(20., 300., 2.),
    # 'pt2':              frange(20,100,10)+frange(100,200,25),
    # 'pt2':              frange(20,150,10),
    # 'pt2':              frange(20,110,5)+frange(120,200,10),
    'pt2':              frange(20, 180, 5),
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
    # 'm_ll':             frange(50.,90.,10)+frange(100., 300.,20.),
    'm_ll':             frange(50., 400., 10),
    'm_ll2':             frange(50., 100., 10) + frange(100., 300., 20.),
    # 'm_llg':             frange(0,90,10)+frange(100,200,10)+frange(200,500,50),
    'm_llg':             frange(50, 400, 10),
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
    'deltaPhiLLG':      frange(0., 3., 0.15),
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
    'mt2':            frange(50., 550., 25.),
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
    'Fakes': frange(0, 5, 1),
    'DeltaPhiL1G': frange(-3., 3., 0.3),
    'DeltaPhiL2G': frange(-3., 3., 0.3),
    'm_l1g': frange(10., 200., 5.),
    'm_l2g': frange(10., 200., 5.)

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
    if __name__ == "__main__":
        print "data", num_data, "fix", add_fix, "toScale", den_toScale
        print "raw fix", additional_fix.GetEntries(
        ), "toScaleFix", denominator_toScale.GetEntries()

    return [alpha, alphaErr / alpha] if den_toScale else [1., 0.]


def drawTTCR(sampleNames, name, binning=None, binningName="", xTitle=None, yTitle=None, weightsToUse=["nISR", "topPt", "ewk"], SF_TT=1., SF_WZ=WZsf, SF_DY=DYsf, SF_ZZ=ZZsf):
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()

    #folder= (name.split("/"))[0]
    folder = (name.split("/"))[0] + "/" + \
        (name.split("/"))[1] + "/" + (name.split("/"))[2]

    style.divideByBinWidth = False

    style.minimumOne = True

    if __name__ == "__main__":
        dataHist = aux.stdHist(dataLL, name, binning)
    else:
        dataHist = aux.stdHist(dataLL, "CRTT/EM/nom/" +
                               (name.split("/"))[-1], binning)
    aux.drawOpt(dataHist, "data")

    zgHist = aux.stdHistWithWeights(zgamma, name, weightsToUse, binning)
    ttgHist = aux.stdHistWithWeights(ttgamma, name, weightsToUse, binning)
    ttg080Hist = aux.stdHistWithWeights(
        ttgamma, name.replace("CRTT", "CRTT080"), weightsToUse, binning)
    ttg80Hist = aux.stdHistWithWeights(
        ttgamma, name.replace("CRTT", "CRTT80"), weightsToUse, binning)
    zzHist = aux.stdHistWithWeights(zz, name, weightsToUse, binning)
    wwgHist = aux.stdHistWithWeights(wwgamma, name, weightsToUse, binning)
    wzgHist = aux.stdHistWithWeights(wzgamma, name, weightsToUse, binning)
    dyHist = aux.stdHistWithWeights(DYjetsNLO, name, weightsToUse, binning)
    wjetsHist = aux.stdHistWithWeights(wjets, name, weightsToUse, binning)
    ttHist = aux.stdHistWithWeights(tt, name, weightsToUse, binning)
    # tt080Hist = aux.stdHistWithWeights(tt, name.replace(
    #     "CRTT", "CRTT080"), weightsToUse, binning)
    # tt80Hist = aux.stdHistWithWeights(tt, name.replace(
    #     "CRTT", "CRTT80"), weightsToUse, binning)
    # if __name__=="__main__":
    #singletopHist=aux.stdHist(singletop, name, binning)
    # else:
    #singletopHist = aux.stdHist(singletop,"CRTT/EM/nom/"+(name.split("/"))[-1],binning)
    if __name__ == "__main__":
        singletopHist = aux.stdHistWithWeights(
            singletop, name, weightsToUse, binning)
    else:
        singletopHist = aux.stdHistWithWeights(
            singletop, "CRTT/EM/nom/" + (name.split("/"))[-1], weightsToUse, binning)
    #singletopHist=aux.stdHist(singletop, name, binning)

    wzHist = aux.stdHistWithWeights(wz, name, weightsToUse, binning)
    wwHist = aux.stdHistWithWeights(ww, name, weightsToUse, binning)
    zz4lHist = aux.stdHistWithWeights(zz4l, name, weightsToUse, binning)
    wgHist = aux.stdHistWithWeights(wgamma, name, weightsToUse, binning)

    zgHist.SetLineColor(ROOT.kGreen - 3)
    ttgHist.SetLineColor(ROOT.kRed + 1)
    # ttg080Hist.SetLineColor(ROOT.kRed+1)
    # ttg80Hist.SetLineColor(ROOT.kRed+1)
    zzHist.SetLineColor(ROOT.kYellow)
    wwgHist.SetLineColor(ROOT.kCyan + 3)
    wzgHist.SetLineColor(ROOT.kBlue + 3)
    dyHist.SetLineColor(ROOT.kGreen + 3)
    wjetsHist.SetLineColor(ROOT.kBlue - 9)
    # ttHist.SetLineColor(ROOT.kOrange+8)
    ttHist.SetLineColor(ROOT.kOrange + 7)
    # tt080Hist.SetLineColor(ROOT.kOrange+8)
    # tt80Hist.SetLineColor(ROOT.kOrange+8)
    singletopHist.SetLineColor(ROOT.kOrange + 4)
    # wzHist.SetLineColor(ROOT.kAzure-5)
    wzHist.SetLineColor(ROOT.kAzure - 6)
    wwHist.SetLineColor(ROOT.kCyan - 2)
    zz4lHist.SetLineColor(ROOT.kOrange - 2)
    wgHist.SetLineColor(ROOT.kRed + 3)

    final_dyHist = aux.addHists(dyHist, zgHist)
    final_ttHist = aux.addHists(ttHist, ttgHist)
    # final_tt080Hist = aux.addHists(tt080Hist, ttg080Hist)
    # final_tt80Hist = aux.addHists(tt80Hist, ttg80Hist)
    final_zzHist = aux.addHists(zzHist, zz4lHist)
    # final_wzHist = aux.addHists(wzHist)
    final_wzHist = wzHist.Clone()
    final_otherHist = aux.addHists(
        wwgHist, wzgHist, wjetsHist, singletopHist, wwHist, wgHist)

    # final_otherHist.SetLineColor(ROOT.kCyan+3)
    final_otherHist.SetLineColor(ROOT.kGray + 2)
    # final_ttHist.SetLineColor(ROOT.kOrange+8)
    final_ttHist.SetLineColor(ROOT.kOrange + 7)
    final_zzHist.SetLineColor(ROOT.kOrange - 2)
    final_wzHist.SetLineColor(ROOT.kAzure - 5)
    final_dyHist.SetLineColor(ROOT.kGreen - 3)

    # scaleHist=aux.addHists(ttHist,ttgHist)
    if __name__ == "__main__":
        print "main"
        # scaleHist=aux.addHists(tt080Hist,tt80Hist,ttg080Hist,ttg80Hist)
        # scaleHist=aux.addHists(ttg080Hist,ttg80Hist)
        # scaleHist=aux.addHists(ttgHist,ttHist)
        scaleHist = final_ttHist.Clone()
        # scaleHist=aux.addHists(ttgHist)
    else:
        # scaleHist=aux.addHists(ttHist,ttgHist)
        scaleHist = final_ttHist.Clone()
    # fixHist=aux.addHists(zgHist,zzHist,wwgHist,wzgHist,dyHist,wjetsHist,singletopHist,wzHist,wwHist,zz4lHist,wgHist)

    # try SCALING
    # print ZZsf, DYsf, WZsf
    final_zzHist.Scale(SF_ZZ)
    final_dyHist.Scale(SF_DY)
    final_wzHist.Scale(SF_WZ)
    # final_zzHist.Scale(ZZsf)
    # final_dyHist.Scale(DYsf)
    # final_wzHist.Scale(WZsf)

    fixHist = aux.addHists(final_otherHist, final_zzHist,
                           final_dyHist, final_wzHist)
    if __name__ == "__main__":
        print "ttRaw", ttHist.GetEntries(), "ttgRaw", ttgHist.GetEntries()
        print "tt", ttHist.Integral(), "ttg", ttgHist.Integral()

    sf, sferr = calculateSFAndError(dataHist, scaleHist, fixHist)
    #sf,sferr = 1.,0.2

    ttHist.Scale(sf)
    ttgHist.Scale(sf)
    # ttHist.Scale(1.89)
    # ttgHist.Scale(0.988)
    # ttHist.Scale(1.)
    # ttgHist.Scale(1.)
    # ttg080Hist.Scale(sf)
    # ttg80Hist.Scale(sf)
    # tt080Hist.Scale(sf)
    # tt80Hist.Scale(sf)

    # final_tt080Hist.Scale(sf)
    # final_tt80Hist.Scale(sf)
    final_ttHist.Scale(sf)
    # final_zzHist.Scale(sfZZ[0])
    # final_wzHist.Scale(sfWZ[0])
    # final_dyHist.Scale(sfDY[0])

    scaleHist.Scale(sf)
    # scaleHist.Scale(1.)

    if "Fakes" in name:
        normFactor = 1. / aux.addHists(ttHist, ttgHist, final_otherHist,
                                       final_dyHist, final_zzHist, final_wzHist).Integral(1, 1)
        ttHist.Scale(normFactor)
        ttgHist.Scale(normFactor)
        final_otherHist.Scale(normFactor)
        final_dyHist.Scale(normFactor)
        final_wzHist.Scale(normFactor)
        final_zzHist.Scale(normFactor)
        dataHist.Scale(normFactor)

    #print sf

    mcSystUncert = 0.0  # SF, lumi, trigger
    #zgSyst = aux.getSysHisto(zgHist, mcSystUncert)
    #ttgSyst = aux.getSysHisto(ttgHist, sferr)
    #zzSyst = aux.getSysHisto(zzHist, mcSystUncert)
    #wwgSyst = aux.getSysHisto(wwgHist, mcSystUncert)
    #wzgSyst = aux.getSysHisto(wzgHist, mcSystUncert)
    #dySyst = aux.getSysHisto(dyHist, mcSystUncert)
    #wjetsSyst = aux.getSysHisto(wjetsHist, mcSystUncert)
    #ttSyst = aux.getSysHisto(ttHist, sferr)
    # singletopSyst=aux.getSysHisto(singletopHist,mcSystUncert)
    # wzSyst=aux.getSysHisto(wzHist,mcSystUncert)
    # wwSyst=aux.getSysHisto(wwHist,mcSystUncert)
    # zz4lSyst=aux.getSysHisto(zz4lHist,mcSystUncert)
    # wgSyst=aux.getSysHisto(wgHist,mcSystUncert)

    # cutBin=ttHist.FindBin(80)

    # ttg080Syst=aux.getSysHisto(ttg080Hist,0.04)
    # ttg80Syst=aux.getSysHisto(ttg80Hist,0.4)
    # tt080Syst=aux.getSysHisto(tt080Hist,0.04)
    # tt80Syst=aux.getSysHisto(tt80Hist,0.4)

    # scaleSyst080=aux.getSysHisto(aux.addHists(ttg080Hist,tt080Hist),0.04)
    # scaleSyst80=aux.getSysHisto(aux.addHists(ttg80Hist,tt80Hist),0.4)
    # scaleSyst80=aux.getSysHisto(aux.addHists(ttg80Hist,tt80Hist),0.2)
    # scaleSyst080=aux.getSysHisto(aux.addHists(ttg080Hist),0.04)
    # scaleSyst80=aux.getSysHisto(aux.addHists(ttg80Hist),0.4)

    # scaleSyst=aux.addHists(tt080Syst,tt80Syst,ttg080Syst,ttg80Syst)

    # scaleSyst=aux.addHists(scaleSyst080,scaleSyst80)

    # scaleSyst=aux.getSysHisto(aux.addHists(ttgHist,ttHist),0.04) ##<<<<<

    # scaleSyst=aux.getSysHisto(aux.addHists(ttgHist,),0.04)

    # scaleSyst=aux.getSysHisto(aux.addHists(ttg80Hist,tt80Hist),0.4)

    # if "pt_g1" in name:
    # scaleSyst=aux.getSysHistoCut(scaleHist,0.04,0.4,cutBin)
    # else:
    # scaleSyst=aux.getSysHisto(scaleHist,sferr)
    # scaleSyst=aux.getSysHisto(scaleHist,0.04)

    otherSyst = aux.getSysHistoWithMeanWeight(
        final_otherHist, mcSystUncert, final_otherHist.Integral() / sum(otherBKG.ngens))

    wzSFSyst = aux.getSysHisto(wzHist, mcSystUncert)
    # final_tt080SFSyst = aux.getSysHisto(final_tt080Hist, sferr)
    # final_tt80SFSyst = aux.getSysHisto(final_tt80Hist, 0.2)
    final_zzSFSyst = aux.getSysHisto(final_zzHist, mcSystUncert)
    final_dySFSyst = aux.getSysHistoWithMeanWeight(
        final_dyHist, mcSystUncert, final_dyHist.Integral() / (sum(DYjetsNLO.ngens) + sum(zgamma.ngens)))

    # final_ttTotSFSyst=aux.addHists(final_tt080SFSyst,final_tt80SFSyst)
    final_ttTotSFSyst = aux.getSysHisto(final_ttHist, sferr)

    totStat = aux.addHists(final_zzHist, final_dyHist,
                           final_ttHist, final_wzHist, final_otherHist)
    totSyst = aux.addHists(final_dySFSyst, final_zzSFSyst,
                           final_ttTotSFSyst, wzSFSyst, otherSyst)

    #dataSyst= aux.getSysHisto(dataHist,mcSystUncert)

    #totStat = aux.addHists(zgHist, ttgHist, zzHist, wwgHist, wzgHist,dyHist,wjetsHist,ttHist,singletopHist,wzHist,wwHist,zz4lHist,wgHist)
    #totSyst = aux.addHists(zgSyst, zzSyst, wwgSyst, wzgSyst,dySyst,wjetsSyst,singletopSyst,wzSyst,wwSyst,zz4lSyst,wgSyst,scaleSyst)
    #totStat = aux.addHists(zgHist, ttgHist, zzHist, wwgHist, wzgHist,dyHist,wjetsHist,singletopHist,wzHist,wwHist,zz4lHist,wgHist)
    #totSyst = aux.addHists(zgSyst, zzSyst, wwgSyst, wzgSyst,dySyst,wjetsSyst,singletopSyst,wzSyst,wwSyst,zz4lSyst,wgSyst,scaleSyst)

    totUnc = aux.addHistUncert(totStat, totSyst)
    aux.drawOpt(totUnc, "totUnc")
    aux.drawOpt(totSyst, "sysUnc")

    #signal2 = aux.stdHist(t5bbbbzg_1500_400, name, binning)
    #signal1 = aux.stdHist(tching_400, name, binning)
    #signal2 = aux.stdHist(gmsb_240_230, name, binning)
    #signal1 = aux.stdHist(gmsb_290_205, name, binning)

    ##signal1 = aux.stdHist(tching_1200, name, binning)

    ##signal3 = aux.stdHist(t6ttZg_600_300, name, binning)
    # signal4 = aux.stdHist(t6ttZg_600_200, name, binning)
    # signal5 = aux.stdHist(t6ttZg_400_200, name, binning)

    # for h in signal1, signal2:
    # for h in signal1, signal2,signal3:
    # for h in signal4, signal5:
    #     aux.drawOpt(h, "signal")
    # signal1.SetLineColor(ROOT.kBlue+3)
    # signal2.SetLineColor(ROOT.kBlue+3)
    # signal2.SetLineStyle(2)
    # signal3.SetLineColor(ROOT.kRed+3)
    # signal3.SetLineStyle(2)

    #m.add(signal2, "T5bbbbZg")
    #m.add(signal1, "TChiNg")
    #m.add(signal3, "T6ttZg")

    # m.add(signal4, "T6ttZg")
    # m.add(signal5, "T6ttZg")
    addedHists = aux.addHists(
        final_otherHist, final_dyHist, final_zzHist, final_wzHist)

    dataHist.GetXaxis().SetTitle(xTitle)
    ttgHist.GetXaxis().SetTitle(xTitle)
    addedHists.GetXaxis().SetTitle(xTitle)

    if "Fakes" in name:
        labels = ["all", "genPhoton",
                  "genElectron", "genJet", "noMatch"]
        for i in range(len(labels)):
            ttgHist.GetXaxis().SetBinLabel(i + 1, labels[i])
            ttgHist.GetXaxis().LabelsOption("d")
            ttgHist.SetLabelOffset(10.5, "X")
            ttgHist.GetXaxis().SetLabelOffset(10.5)
            ttgHist.LabelsOption("d", "X")
            ttHist.GetXaxis().SetBinLabel(i + 1, labels[i])
            ttHist.GetXaxis().LabelsOption("d")
            ttHist.LabelsOption("d", "X")
            ttHist.SetLabelOffset(10.5, "X")
            ttHist.GetXaxis().SetLabelOffset(10.5)
            addedHists.GetXaxis().SetBinLabel(i + 1, labels[i])
            addedHists.GetXaxis().LabelsOption("d")
            addedHists.LabelsOption("d", "X")
            addedHists.SetLabelOffset(10.5, "X")
            addedHists.GetXaxis().SetLabelOffset(10.5)
            dataHist.GetXaxis().SetBinLabel(i + 1, labels[i])
            dataHist.GetXaxis().LabelsOption("d")
            dataHist.LabelsOption("d", "X")
            dataHist.GetXaxis().SetLabelOffset(10.5)
            ttgHist.SetYTitle("arbitrary units")
            ttHist.SetYTitle("arbitrary units")
            addedHists.SetYTitle("arbitrary units")
            dataHist.SetYTitle("arbitrary units")

    # m.addStack(zgHist, "Z#gamma")
    m.addStack(ttgHist, "t#bar{t}#gamma")
    #m.addStack(zzHist, "ZZ")
    # m.addStack(wwgHist, "WW#gamma")
    # m.addStack(wzgHist, "WZ#gamma")
    #m.addStack(dyHist, "Drell-Yan/Z")
    #m.addStack(wjetsHist, "W+jets")
    m.addStack(ttHist, "t#bar{t}")
    #m.addStack(singletopHist, "single t")
    #m.addStack(wwHist, "WW")
    #m.addStack(wzHist, "WZ")
    #m.addStack(zz4lHist, "ZZ(4l)")
    # m.addStack(wgHist, "W#gamma")

    m.addStack(addedHists, "Other")

    m.sortStackByIntegral()

    if not "Fakes" in name:
        m.add(dataHist, "Data")
        m.add(totUnc, "Total uncertainty")
        m.add(totSyst, "syst. uncertainty")

    legInfo = "t#bar{t}(#gamma) Control Region"
    m.leg.SetHeader(legInfo)
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
        ksText.DrawLatexNDC(0.57, 0.65, "KS Value= " + str(np.round(KS, 3)))

        purity = aux.addHists(ttHist, ttgHist).Integral() / aux.addHists(ttHist, ttgHist,
                                                                         final_otherHist, final_dyHist, final_zzHist, final_wzHist).Integral() * 100.
        purText = ROOT.TLatex()
        purText.SetTextSize(0.45 * purText.GetTextSize())
        purText.DrawLatexNDC(0.57, 0.6, "purity= " +
                             str(np.round(purity, 2)) + " %")

        r = ratio.Ratio(
            "#scale[.9]{#lower[.24]{#splitline{Data/Pred.}{Bkg. frac.}}}", dataHist, totStat, sysHisto=totSyst)
        #rMax = 1.5
        rMax = 2.
        #r.draw(0., rMax, m.getStack(), True)
        r.draw(0., rMax, m.getStack())
    else:
        # r = ratio.Ratio(
        #     "#scale[.9]{#lower[.24]{#splitline{Data/Pred.}{Bkg. frac.}}}", dataHist, totStat)
        r = ratio.Ratio(
            "#scale[.9]{#lower[.24]{#splitline{Data/Pred.}{Bkg. frac.}}}", totStat, totStat)
        #rMax = 1.5
        rMax = 2.
        #r.draw(0., rMax, m.getStack(), True)
        r.draw(0., rMax, m.getStack())

    #aux.Label(sim=False, status="Work in Progress")
    # aux.Label(sim=False, status="Private Work")
    aux.Label(status="")
    #aux.save(name, normal=False, changeMinMax=False)
    directory = "plots_CR_tt/"
    # directory = "plots_CR_tt_rishi/"
    if not os.path.exists(directory):
        os.makedirs(directory)
    if __name__ == "__main__":
        if not "Fakes" in name:
            aux.save(name.replace("/", "_"), folder=directory)
        else:
            aux.save(name.replace("/", "_"),
                     folder=directory, changeMinMax=False)
    if __name__ == "__main__":
        print "SF", sf, "SFErr", sferr

    return sf, sferr


def main():
    # variables=["eta1","eta2","pt1","pt2","n_jets","n_vtx","phi1","phi2","m_ll","ht","n_photons","pt_g1","met","nElectrons","nMuons","deltaR1_g1","deltaR2_g1","deltaRLL"]
    # variables=["mt2","eta1","eta2","pt1","pt2","n_jets","phi1","phi2","m_ll","ht","n_photons","pt_g1","met","Fakes","deltaPhiLLG","DeltaPhiL1G","DeltaPhiL2G","m_l1g","m_l2g","m_llg"]
    variables = ["mt2", "eta1", "eta2", "pt1", "pt2", "n_jets", "phi1", "phi2", "m_ll", "ht",
                 "n_photons", "pt_g1", "met", "deltaPhiLLG", "DeltaPhiL1G", "DeltaPhiL2G", "m_l1g", "m_l2g", "m_llg", "eta_g1"]
    # variables = ["Fakes"]
    # variables=["pt_g1"]
    # variables = ["pt1"]
    # variables=["sigmaIetaIeta_g1","r9_g1"]
    bkgs = [DYjetsNLO, zgamma, tt, ttgamma, wwgamma, wzgamma,
            zz, wjets, wgamma, singletop, zz4l, wz, zz, ww]
    # groups=["ControlRegionTT"]
    groups = ["CRTT"]
    binnings_ = binnings.copy()
    # binnings_["pt_g1"]=frange(25,80,5)+frange(80,140,10)+frange(140,200,20)
    # binnings_["pt_g1"]=frange(25,75,5)+frange(80,130,10)+frange(140,200,20)
    # binnings_["pt_g1"]=frange(20,200,10)
    # binnings_["pt_g1"]=frange(20,120,5)+frange(130,200,7)
    # binnings_["pt_g1"]=frange(20,250,5)
    # binnings_["pt_g1"]=frange(20,250,10)
    # binnings_["pt_g1"]=frange(20,150,10)+frange(165,225,15)
    # binnings_["pt_g1"]=frange(20,150,10)+frange(160,300,10)+frange(325,400,25)
    # binnings_["pt_g1"]=frange(20,150,10)+frange(160,250,10)+frange(275,400,25)
    binnings_["pt_g1"] = frange(20, 300, 10)
    # binnings_["met"]=frange(0,200,10)+frange(225,250,25)+frange(300,500,50)
    binnings_["met"] = frange(0, 300, 10)
    # binnings_["pt_g1"]=frange(20,180,5)+frange(190,250,10)
    toSave = {}
    toSave["EM"] = {}
    for group in groups:
        for variable in variables:
            toSave["EM"][variable] = drawTTCR("EM", group + "/EM/nom/" + variable, binning=binnings_[
                                              variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
        if not os.path.exists("plots_CR_tt/factors"):
            os.makedirs("plots_CR_tt/factors")
        pkl.dump(toSave, open("plots_CR_tt/factors/" + group + ".pkl", "wb"))


if __name__ == "__main__":
    main()
