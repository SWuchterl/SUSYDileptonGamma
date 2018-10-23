from dataMC import labels, frange

import ROOT
from ROOT import *
from array import array
from include import *
import numpy as np
import pickle as pkl
import os

binnings = {
    # 'pt1':              frange(0,100,10)+frange(100,200,25)+range(200,350,50),
    # 'pt1':              frange(20,100,10)+frange(100,200,25)+range(200,350,50),
    # 'pt1':              frange(25,160,5),
    'pt1':              frange(20, 350, 10),
    # 'pt2':              frange(0,100,10)+frange(100,200,25),
    # 'pt2':              frange(20,100,10)+frange(100,200,25),
    'pt2':              frange(20, 200, 10),
    'pt3':              range(0, 200, 10),
    'pt4':              range(0, 200, 10),
    'eta1':             frange(0., 2.4, 0.1),
    'eta2':             frange(0., 2.4, 0.1),
    'eta3':             frange(0., 2.4, 0.1),
    'eta4':             frange(0., 2.4, 0.1),
    'phi1':             frange(0., 3.10, 0.1),
    'phi2':             frange(0., 3.10, 0.1),
    'phi3':             frange(0., 3.10, 0.1),
    'phi4':             frange(0., 3.10, 0.1),
    'ht':               frange(0., 1000., 50),
    # 'met':               [0,25,50,75,100,150,190,230,500],
    # 'met':               frange(0,105,10),
    'met':               frange(0, 105, 1),
    # 'm_ll':             frange(50.,100.,10)+frange(100., 300.,20.),
    'm_ll':             frange(81., 101., 1),
    'm_ll2':             frange(50., 100., 10) + frange(100., 300., 20.),
    'm_llg':             frange(0, 100, 10) + frange(100, 200, 10) + frange(200, 500, 50),
    'pt_llg':             frange(0, 100, 10) + frange(100, 200, 10) + frange(200, 500, 50),
    'n_jets':           frange(0., 7., 1),
    'n_photons':         frange(0., 4., 1),
    'n_vtx':            frange(0., 40., 1),
    # 'pt_g1':            frange(20,100,10)+frange(125,150,25)+frange(200,250,50),
    'pt_g1':            frange(20, 350, 10),
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
    'mt2':            frange(80., 150., 5.),
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


def drawDYCR(sampleNames, name, datasetToUse, binning=None, binningName="", xTitle=None, yTitle=None, weightsToUse=["nISR", "topPt", "ewk"]):
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()

    #folder= (name.split("/"))[0]
    folder = (name.split("/"))[0] + "/" + \
        (name.split("/"))[1] + "/" + (name.split("/"))[2]

    style.divideByBinWidth = False

    style.minimumOne = True

    if __name__ == "__main__":
        dataHist = aux.stdHist(datasetToUse, name, binning)
    else:
        dataHist = aux.stdHist(
            datasetToUse, "CRDY/LL/nom/" + (name.split("/"))[-1], binning)
    #dataHist = aux.stdHist(datasetToUse, name, binning)
    aux.drawOpt(dataHist, "data")

    zgHist = aux.stdHistWithWeights(zgamma, name, weightsToUse, binning)
    ttgHist = aux.stdHistWithWeights(ttgamma, name, weightsToUse, binning)
    zzHist = aux.stdHistWithWeights(zz, name, weightsToUse, binning)
    wwgHist = aux.stdHistWithWeights(wwgamma, name, weightsToUse, binning)
    wzgHist = aux.stdHistWithWeights(wzgamma, name, weightsToUse, binning)
    dyHist = aux.stdHistWithWeights(DYjetsNLO, name, weightsToUse, binning)
    wjetsHist = aux.stdHistWithWeights(wjets, name, weightsToUse, binning)
    ttHist = aux.stdHistWithWeights(tt, name, weightsToUse, binning)
    #singletopHist=aux.stdHist(singletop, name, binning)
    if __name__ == "__main__":
        singletopHist = aux.stdHistWithWeights(
            singletop, name, weightsToUse, binning)
    else:
        singletopHist = aux.stdHistWithWeights(
            singletop, "CRDY/LL/nom/" + (name.split("/"))[-1], weightsToUse, binning)
    wzHist = aux.stdHistWithWeights(wz, name, weightsToUse, binning)
    wwHist = aux.stdHistWithWeights(ww, name, weightsToUse, binning)
    zz4lHist = aux.stdHistWithWeights(zz4l, name, weightsToUse, binning)
    wgHist = aux.stdHistWithWeights(wgamma, name, weightsToUse, binning)

    zgHist.SetLineColor(ROOT.kGreen - 3)
    ttgHist.SetLineColor(ROOT.kRed + 1)
    zzHist.SetLineColor(ROOT.kYellow)
    wwgHist.SetLineColor(ROOT.kCyan + 3)
    wzgHist.SetLineColor(ROOT.kBlue + 3)
    dyHist.SetLineColor(ROOT.kGreen + 3)
    wjetsHist.SetLineColor(ROOT.kBlue - 9)
    ttHist.SetLineColor(ROOT.kOrange + 8)
    singletopHist.SetLineColor(ROOT.kOrange + 4)
    wzHist.SetLineColor(ROOT.kAzure - 5)
    wwHist.SetLineColor(ROOT.kCyan - 2)
    zz4lHist.SetLineColor(ROOT.kOrange - 2)
    wgHist.SetLineColor(ROOT.kRed + 3)

    final_dyHist = aux.addHists(dyHist, zgHist)
    final_ttHist = aux.addHists(ttHist, ttgHist)
    final_zzHist = aux.addHists(zzHist, zz4lHist)
    final_wzHist = aux.addHists(wzHist)
    final_otherHist = aux.addHists(
        wwgHist, wzgHist, wjetsHist, singletopHist, wwHist, wgHist)

    # final_otherHist.SetLineColor(ROOT.kCyan+3)
    final_otherHist.SetLineColor(ROOT.kGray + 2)
    final_ttHist.SetLineColor(ROOT.kOrange + 8)
    final_zzHist.SetLineColor(ROOT.kOrange - 2)
    final_wzHist.SetLineColor(ROOT.kAzure - 5)
    final_dyHist.SetLineColor(ROOT.kGreen - 3)

    # scaleHist=aux.addHists(dyHist,zgHist)
    scaleHist = final_dyHist.Clone()
    # fixHist=aux.addHists(ttgHist,zzHist,wwgHist,wzgHist,ttHist,wjetsHist,singletopHist,wzHist,wwHist,zz4lHist,wgHist)
    fixHist = aux.addHists(final_otherHist, final_ttHist,
                           final_wzHist, final_zzHist)

    sf, sferr = calculateSFAndError(dataHist, scaleHist, fixHist)

    #print sf

    dyHist.Scale(sf)
    zgHist.Scale(sf)
    scaleHist.Scale(sf)

    final_dyHist.Scale(sf)

    if "Fakes" in name:
        normFactor = 1. / aux.addHists(final_ttHist, final_otherHist,
                                       dyHist, zgHist, final_zzHist, final_wzHist).Integral(1, 1)
        final_ttHist.Scale(normFactor)
        final_otherHist.Scale(normFactor)
        dyHist.Scale(normFactor)
        zgHist.Scale(normFactor)
        final_wzHist.Scale(normFactor)
        final_zzHist.Scale(normFactor)
        dataHist.Scale(normFactor)

    mcSystUncert = 0.0  # SF, lumi, trigger
    #zgSyst = aux.getSysHisto(zgHist, mcSystUncert)
    #ttgSyst = aux.getSysHisto(ttgHist, mcSystUncert)
    #zzSyst = aux.getSysHisto(zzHist, mcSystUncert)
    #wwgSyst = aux.getSysHisto(wwgHist, mcSystUncert)
    #wzgSyst = aux.getSysHisto(wzgHist, mcSystUncert)
    #dySyst = aux.getSysHisto(dyHist, mcSystUncert)
    #wjetsSyst = aux.getSysHisto(wjetsHist, mcSystUncert)
    #ttSyst = aux.getSysHisto(ttHist, mcSystUncert)
    # singletopSyst=aux.getSysHisto(singletopHist,mcSystUncert)
    # wzSyst=aux.getSysHisto(wzHist,mcSystUncert)
    # wwSyst=aux.getSysHisto(wwHist,mcSystUncert)
    # zz4lSyst=aux.getSysHisto(zz4lHist,mcSystUncert)
    # wgSyst=aux.getSysHisto(wgHist,mcSystUncert)

    # cutBin=ttHist.FindBin(80)
    #
    #
    # if "pt_g1" in name:
    # scaleSyst=aux.getSysHistoCut(scaleHist,0.04,0.4,cutBin)
    # else:
    # scaleSyst=aux.getSysHisto(scaleHist,sferr)

    #dataSyst= aux.getSysHisto(dataHist,mcSystUncert)

    #totStat = aux.addHists(zgHist, ttgHist, zzHist, wwgHist, wzgHist,dyHist,wjetsHist,ttHist,singletopHist,wzHist,wwHist,zz4lHist,wgHist)
    #totSyst = aux.addHists(ttgSyst, zzSyst, wwgSyst, wzgSyst,ttSyst,wjetsSyst,singletopSyst,wzSyst,wwSyst,zz4lSyst,wgSyst,scaleSyst)

    otherSyst = aux.getSysHistoWithMeanWeight(
        final_otherHist, mcSystUncert, final_otherHist.Integral() / sum(otherBKG.ngens))
    wzSFSyst = aux.getSysHisto(wzHist, mcSystUncert)
    final_zzSFSyst = aux.getSysHisto(final_zzHist, mcSystUncert)
    final_ttSFSyst = aux.getSysHisto(final_ttHist, mcSystUncert)
    #final_dySFSyst= aux.getSysHistoWithMeanWeight(final_dyHist,sferr,final_dyHist.Integral()/(sum(DYjetsNLO.ngens)+sum(zgamma.ngens)))
    final_dySFSyst = aux.getSysHisto(final_dyHist, sferr)
    print sferr

    totStat = aux.addHists(final_zzHist, final_dyHist,
                           final_ttHist, final_wzHist, final_otherHist)
    totSyst = aux.addHists(final_dySFSyst, final_zzSFSyst,
                           final_ttSFSyst, wzSFSyst, otherSyst)

    totUnc = aux.addHistUncert(totStat, totSyst)
    aux.drawOpt(totUnc, "totUnc")
    aux.drawOpt(totSyst, "sysUnc")

    #signal2 = aux.stdHist(t5bbbbzg_1500_400, name, binning)
    #signal1 = aux.stdHist(tching_400, name, binning)
    #signal1 = aux.stdHist(tching_1200, name, binning)

    #signal2_AvgTopPtWeightHisto = t5bbbbzg_1500_400.getHist(folder+"/weight_topPt")
    #signal2_AvgNIsrWeightHisto = t5bbbbzg_1500_400.getHist(folder+"/weight_nISR")
    #signal2_AvgEWKinoWeightHisto = t5bbbbzg_1500_400.getHist(folder+"/weight_EWKinoPairPt")
    #signal1_AvgTopPtWeightHisto = tching_600.getHist(folder+"/weight_topPt")
    #signal1_AvgNIsrWeightHisto = tching_600.getHist(folder+"/weight_nISR")
    #signal1_AvgEWKinoWeightHisto = tching_600.getHist(folder+"/weight_EWKinoPairPt")

    # signal2.Scale((1./signal2_AvgTopPtWeightHisto.GetMean()))
    # signal2.Scale((1./signal2_AvgNIsrWeightHisto.GetMean()))
    # signal2.Scale((1./signal2_AvgEWKinoWeightHisto.GetMean()))
    # signal1.Scale((1./signal1_AvgTopPtWeightHisto.GetMean()))
    # signal1.Scale((1./signal1_AvgNIsrWeightHisto.GetMean()))
    # signal1.Scale((1./signal1_AvgEWKinoWeightHisto.GetMean()))

    # for h in signal1, signal2:
    #aux.drawOpt(h, "signal")
    # signal1.SetLineColor(ROOT.kBlue+3)
    # signal2.SetLineColor(ROOT.kBlue+3)
    # signal2.SetLineStyle(2)

    #m.add(signal2, "T5bbbbZg")
    #m.add(signal1, "TChiNg")
    addedHists = aux.addHists(
        final_otherHist, final_zzHist, final_ttHist, final_wzHist)
    if "Fakes" in name:
        labels = ["all", "genPhoton",
                  "genElectron", "genJet", "noMatch"]
        for i in range(len(labels)):
            zgHist.GetXaxis().SetBinLabel(i + 1, labels[i])
            zgHist.GetXaxis().LabelsOption("d")
            zgHist.SetLabelOffset(10.5, "X")
            zgHist.GetXaxis().SetLabelOffset(10.5)
            zgHist.LabelsOption("d", "X")
            dyHist.GetXaxis().SetBinLabel(i + 1, labels[i])
            dyHist.GetXaxis().LabelsOption("d")
            dyHist.LabelsOption("d", "X")
            dyHist.SetLabelOffset(10.5, "X")
            dyHist.GetXaxis().SetLabelOffset(10.5)
            addedHists.GetXaxis().SetBinLabel(i + 1, labels[i])
            addedHists.GetXaxis().LabelsOption("d")
            addedHists.LabelsOption("d", "X")
            addedHists.SetLabelOffset(10.5, "X")
            addedHists.GetXaxis().SetLabelOffset(10.5)
            dataHist.GetXaxis().SetBinLabel(i + 1, labels[i])
            dataHist.GetXaxis().LabelsOption("d")
            dataHist.LabelsOption("d", "X")
            dataHist.GetXaxis().SetLabelOffset(10.5)
            zgHist.SetYTitle("arbitrary units")
            dyHist.SetYTitle("arbitrary units")
            addedHists.SetYTitle("arbitrary units")
            dataHist.SetYTitle("arbitrary units")

    dataHist.GetXaxis().SetTitle(xTitle)
    dyHist.GetXaxis().SetTitle(xTitle)
    zgHist.GetXaxis().SetTitle(xTitle)
    addedHists.GetXaxis().SetTitle(xTitle)
    wzHist.GetXaxis().SetTitle(xTitle)

    m.addStack(zgHist, "Z#gamma")
    # m.addStack(ttgHist, "t#bar{t}#gamma")
    #m.addStack(zzHist, "ZZ")
    # m.addStack(wwgHist, "WW#gamma")
    # m.addStack(wzgHist, "WZ#gamma")
    m.addStack(dyHist, "Drell-Yan/Z")
    #m.addStack(wjetsHist, "W+jets")
    # m.addStack(ttHist, "t#bar{t}")
    #m.addStack(singletopHist, "single t")
    #m.addStack(wwHist, "WW")
    #m.addStack(wzHist, "WZ")
    #m.addStack(zz4lHist, "ZZ(4l)")
    # m.addStack(wgHist, "W#gamma")
    # m.addStack(aux.addHists(wwgHist,ttgHist,zzHist,wzgHist,wjetsHist,ttHist,singletopHist,wwHist,wzHist,zz4lHist,wgHist),"other")
    # m.addStack(aux.addHists(final_otherHist,final_zzHist,final_ttHist,final_wzHist),"other")
    m.addStack(addedHists, "other")

    m.sortStackByIntegral()
    if not "Fakes" in name:
        if "EE" in sampleNames:
            m.add(dataHist, "Data (ee)")
        elif "MM" in sampleNames:
            m.add(dataHist, "Data (#mu#mu)")
        else:
            m.add(dataHist, "Data (ee+#mu#mu)")
        m.add(totUnc, "Total uncertainty")
        m.add(totSyst, "syst. uncertainty")

    legInfo = "DY/Z(#gamma) Control Region"
    m.leg.SetHeader(legInfo)
    # m.leg.SetY1(.56)
    # m.leg.SetX1(.56)
    # m.leg.SetX1(.46)
    # m.leg.SetX2(.99)
    # m.leg.SetX2(.89)

    m.Draw()
    if not "Fakes" in name:
        # KS=totStat.Clone().KolmogorovTest(dataHist.Clone(),"UO")
        KS = totStat.Clone().KolmogorovTest(dataHist.Clone(), "UOD")
        # KS=totStat.Clone().KolmogorovTest(dataHist.Clone())
        # KS=totStat.Clone().KolmogorovTest(dataHist.Clone(),"UOXN")
        ksText = ROOT.TLatex()
        ksText.SetTextSize(0.45 * ksText.GetTextSize())
        ksText.DrawLatexNDC(0.57, 0.65, "KS Value= " + str(np.round(KS, 5)))

        purity = aux.addHists(dyHist, zgHist).Integral() / aux.addHists(ttHist, ttgHist,
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
        r = ratio.Ratio(
            "#scale[.9]{#lower[.24]{#splitline{Data/Pred.}{Bkg. frac.}}}", dataHist, totStat)
        #rMax = 1.5
        rMax = 2.
        #r.draw(0., rMax, m.getStack(), True)
        r.draw(0., rMax, m.getStack())

    aux.Label(sim=False, status="Work in Progress")
    #aux.save(name, normal=False, changeMinMax=False)
    directory = "plots_CR_dy/"
    # directory = "plots_CR_dy_rishi/"
    if not os.path.exists(directory):
        os.makedirs(directory)
    if __name__ == "__main__":
        # aux.save(name.replace("/","_"),folder=directory)
        if not "Fakes" in name:
            aux.save(name.replace("/", "_"), folder=directory)
        else:
            aux.save(name.replace("/", "_"),
                     folder=directory, changeMinMax=False)

    return sf, sferr


def main():
    # variables=["eta1","eta2","pt1","pt2","n_jets","n_vtx","phi1","phi2","m_ll","ht","n_photons","pt_g1","met","nElectrons","nMuons","deltaR1_g1","deltaR2_g1","deltaRLL"]
    # variables=["mt2","eta1","eta2","pt1","pt2","n_jets","phi1","phi2","m_ll","ht","n_photons","pt_g1","met","Fakes"]
    variables = ["mt2", "eta1", "eta2", "pt1", "pt2", "n_jets",
                 "phi1", "phi2", "m_ll", "ht", "n_photons", "pt_g1", "met", "eta_g1"]
    # variables = ["Fakes"]
    # variables=["pt_g1"]
    bkgs = [DYjetsNLO, zgamma, tt, ttgamma, wwgamma, wzgamma,
            zz, wjets, wgamma, singletop, zz4l, wz, zz, ww]
    groups = ["CRDY"]
    binnings_ = binnings.copy()
    # binnings_["pt_g1"]=frange(25,80,5)+frange(80,140,10)+frange(140,200,20)
    # binnings_["met"]=frange(0,100,10)
    binnings_["met"] = frange(0, 100, 5)
    toSave = {}
    # toSave["EE"]={}
    # toSave["MM"]={}
    toSave["LL"] = {}
    for group in groups:
        for variable in variables:
            #toSave["EE"][variable] = drawTTCR("EE",group+"EE/"+variable,dataDoubleEG,binning=binnings_[variable],xTitle=labels[variable][0])
            #toSave["MM"][variable] = drawTTCR("MM",group+"MM/"+variable,dataDoubleMuon,binning=binnings_[variable],xTitle=labels[variable][0])
            #toSave["LL"][variable] = drawTTCR("LL",group+"/"+variable,dataDoubleSF,binning=binnings_[variable],xTitle=labels[variable][0])
            #toSave["EE"][variable] = drawTTCR("EE",group+"/EE/"+variable,dataLL,binning=binnings_[variable],xTitle=labels[variable][0])
            #toSave["MM"][variable] = drawTTCR("MM",group+"/MM/"+variable,dataLL,binning=binnings_[variable],xTitle=labels[variable][0])
            #toSave["LL"][variable] = drawTTCR("LL",group+"/LL/"+variable,dataLL,binning=binnings_[variable],xTitle=labels[variable][0])
            #toSave["EE"][variable] = drawDYCR("EE",group+"/EE/nom/"+variable,dataLL,binning=binnings_[variable],xTitle=labels[variable][0],weightsToUse=["nISR","topPt","ewk"])
            #toSave["MM"][variable] = drawDYCR("MM",group+"/MM/nom/"+variable,dataLL,binning=binnings_[variable],xTitle=labels[variable][0],weightsToUse=["nISR","topPt","ewk"])
            toSave["LL"][variable] = drawDYCR("LL", group + "/LL/nom/" + variable, dataLL, binning=binnings_[
                                              variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
        if not os.path.exists("plots_CR_dy/factors"):
            os.makedirs("plots_CR_dy/factors")
        pkl.dump(toSave, open("plots_CR_dy/factors/" + group + ".pkl", "wb"))


if __name__ == "__main__":
    main()
