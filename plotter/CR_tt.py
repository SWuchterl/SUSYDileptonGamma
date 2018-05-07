from dataMC import labels,frange

import ROOT
from ROOT import *
from array import array
from include import *
import numpy as np
import pickle as pkl
import os

binnings = {
    #'pt1':              frange(20,100,10)+frange(100,200,25)+range(200,350,50),
    #'pt1':              frange(25,160,5),
    'pt1':              frange(20,160,10),
    'pt2':              frange(20,100,10)+frange(100,200,25),
    'pt3':              range(0, 200, 10),
    'pt4':              range(0, 200, 10),
    'eta1':             frange(0., 2.45, 0.1),
    'eta2':             frange(0., 2.45, 0.1),
    'eta3':             frange(0., 2.45, 0.1),
    'eta4':             frange(0., 2.45, 0.1),
    'phi1':             frange(0., 3.15, 0.1),
    'phi2':             frange(0., 3.15, 0.1),
    'phi3':             frange(0., 3.15, 0.1),
    'phi4':             frange(0., 3.15, 0.1),
    'ht':               frange(0., 1000.,50),
    #'met':               [0,25,50,75,100,150,190,230,500],
    'met':               [0,25,50,75,100,150,250,450],
    'm_ll':             frange(50.,100.,10)+frange(100., 300.,20.),
    'm_ll2':             frange(50.,100.,10)+frange(100., 300.,20.),
    'm_llg':             frange(0,100,10)+frange(100,200,10)+frange(200,500,50),
    'pt_llg':             frange(0,100,10)+frange(100,200,10)+frange(200,500,50),
    'n_jets':           frange(0.,7.,1),
    'n_photons':         frange(0.,4.,1),
    'n_vtx':            frange(0.,40.,1),
    #'pt_g1':            frange(20,100,10)+frange(100,150,25)+frange(150,250,50),
    "pt_g1":            frange(25,80,5)+frange(80,140,10)+frange(140,200,20),
    'eta_g1':           frange(0., 2.60, 0.1),
    'phi_g1':           frange(0., 3.50, 0.1),
    'sigmaIetaIeta_g1': frange(0., 0.04,0.01),
    'sigmaIphiIphi_g1': frange(0., 0.2,0.01),
    'deltaR1_g1':       frange(0., 1.,0.05),
    'deltaR2_g1':       frange(0., 1.,0.05),
    'r9_g1':            frange(0., 1.5,0.05),
    'hOverE_g1':        frange(0., 0.1,0.01),
    'deltaEtaLL':       frange(0.,6.,0.1),
    'deltaPhiLL':       frange(0.,6.,0.1),
    'deltaEtaLLG':      frange(0.,6.,0.1),
    'deltaPhiLLG':      frange(0.,6.,0.01),
    'deltaRLL':         frange(0.,1.,0.1),
    'deltaRLLG':        frange(0.,6.,0.1),
    'st':               frange(0.,1000.,100),
    'stmet':            frange(0,1000,100)+frange(1000,4000,500),
    'zpt':              frange(0.,2000.,100),
    'mtll':             frange(0,200,25)+frange(200,500,50)+frange(500,1250,250),
    'mtllg':            frange(0,200,25)+frange(200,500,50)+frange(500,1750,250),
    'mtgmet':            frange(0,200,25)+frange(200,500,50)+frange(500,1750,250),
    'mtl1met':            frange(0,200,25)+frange(200,500,50)+frange(500,1750,250),
    'mtl2met':            frange(0,200,25)+frange(200,500,50)+frange(500,1750,250),
    'mtllmet':            frange(0,200,25)+frange(200,500,50)+frange(500,1750,250),
    'mtllgmet':            frange(0,200,25)+frange(200,500,50)+frange(500,1750,250),
    'mt2':            frange(0.,425.,25.),
    'mzg_exo':            frange(0,100,50)+frange(100,200,50)+frange(200,500,50),
    'gammaMotherID':            frange(0.,200.,1.),
    'genPhotonPT':            np.concatenate((np.arange(0,100,10),np.arange(100,350,50)),axis=0),
    'genPhotonPT_Veto':            np.concatenate((np.arange(0,100,10),np.arange(100,350,50)),axis=0),
    'PhotonPT_Veto':            np.concatenate((np.arange(0,100,10),np.arange(100,350,50)),axis=0),
    'genPhotonPT_NoVeto':            np.concatenate((np.arange(0,100,10),np.arange(100,350,50)),axis=0),
    'PhotonPT_NoVeto':            np.concatenate((np.arange(0,100,10),np.arange(100,350,50)),axis=0),
    'VetoCompare':            np.arange(0.,200.,1.),
    'DeltaPhiLLMet':            frange(0.,6.,0.3),
    'DeltaEtaLLMet':            frange(0.,6.,0.3),
    'DeltaRLLMet':            frange(0.,6.,0.3),
    'nElectrons': frange(0,10,1),
    'nMuons': frange(0,10,1),
    'mTL3Met': frange(0.,300.,10)
}
def calculateSFAndError( numerator_data, denominator_toScale,additional_fix):
    num_dataErr=(ROOT.Double(0))
    num_data = numerator_data.IntegralAndError(0,-1,num_dataErr)
    
    den_toScaleErr=(ROOT.Double(0))
    den_toScale = denominator_toScale.IntegralAndError(0,-1,den_toScaleErr)

    add_fixError=ROOT.Double(0)
    add_fix = additional_fix.IntegralAndError(0,-1,add_fixError)

    alpha = (num_data-add_fix)/(den_toScale)
    alphaErr =np.sqrt( (num_dataErr/den_toScale)**2. +(add_fixError/den_toScale)**2. +(den_toScaleErr*(num_data-add_fix)/(den_toScale)**2.)**2. )
    
    return [alpha,alphaErr/alpha] if den_toScale else [1.,0.]

def drawTTCR(sampleNames, name, binning=None, binningName="", xTitle=None, yTitle=None):
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()
    
    folder= (name.split("/"))[0]
    
    style.divideByBinWidth = False
    
    style.minimumOne=True
    
    #dataHist = aux.stdHist(dataMuonEG, name, binning)
    dataHist = aux.stdHist(dataLL, name, binning)
    aux.drawOpt(dataHist, "data")
    
    zgHist = aux.stdHist(zgamma, name, binning)
    ttgHist = aux.stdHist(ttgamma, name, binning)
    ttg080Hist = aux.stdHist(ttgamma, name.replace("EM","080EM"), binning)
    ttg80Hist = aux.stdHist(ttgamma, name.replace("EM","80EM"), binning)
    zzHist = aux.stdHist(zz, name, binning)
    wwgHist = aux.stdHist(wwgamma, name, binning)
    wzgHist = aux.stdHist(wzgamma, name, binning)
    dyHist = aux.stdHist(DYjetsNLO, name, binning)
    wjetsHist = aux.stdHist(wjets, name, binning)
    ttHist = aux.stdHist(tt, name, binning)
    tt080Hist = aux.stdHist(tt, name.replace("EM","080EM"), binning)
    tt80Hist = aux.stdHist(tt, name.replace("EM","80EM"), binning)
    singletopHist=aux.stdHist(singletop, name, binning)
    wzHist=aux.stdHist(wz, name, binning)
    wwHist=aux.stdHist(ww, name, binning)
    zz4lHist=aux.stdHist(zz4l, name, binning)
    wgHist=aux.stdHist(wgamma, name, binning)
    
    
    #Scaling
    zg_AvgTopPtWeightHisto = zgamma.getHist(folder+"/weight_topPt")
    zg_AvgNIsrWeightHisto = zgamma.getHist(folder+"/weight_nISR")
    zg_AvgEWKinoWeightHisto = zgamma.getHist(folder+"/weight_EWKinoPairPt")
    zg_AvgleptonWeightHisto = zgamma.getHist(folder+"/weight_leptonPairPt")
    ttg_AvgTopPtWeightHisto = ttgamma.getHist(folder+"/weight_topPt")
    ttg_AvgNIsrWeightHisto = ttgamma.getHist(folder+"/weight_nISR")
    ttg_AvgEWKinoWeightHisto = ttgamma.getHist(folder+"/weight_EWKinoPairPt")
    ttg_AvgleptonWeightHisto = ttgamma.getHist(folder+"/weight_leptonPairPt")
    ttg080_AvgTopPtWeightHisto = ttgamma.getHist(folder.replace("EM","080EM")+"/weight_topPt")
    ttg080_AvgNIsrWeightHisto = ttgamma.getHist(folder.replace("EM","080EM")+"/weight_nISR")
    ttg080_AvgEWKinoWeightHisto = ttgamma.getHist(folder.replace("EM","080EM")+"/weight_EWKinoPairPt")
    ttg080_AvgleptonWeightHisto = ttgamma.getHist(folder.replace("EM","080EM")+"/weight_leptonPairPt")
    ttg80_AvgTopPtWeightHisto = ttgamma.getHist(folder.replace("EM","80EM")+"/weight_topPt")
    ttg80_AvgNIsrWeightHisto = ttgamma.getHist(folder.replace("EM","80EM")+"/weight_nISR")
    ttg80_AvgEWKinoWeightHisto = ttgamma.getHist(folder.replace("EM","80EM")+"/weight_EWKinoPairPt")
    ttg80_AvgleptonWeightHisto = ttgamma.getHist(folder.replace("EM","80EM")+"/weight_leptonPairPt")
    zz_AvgTopPtWeightHisto = zz.getHist(folder+"/weight_topPt")
    zz_AvgNIsrWeightHisto = zz.getHist(folder+"/weight_nISR")
    zz_AvgEWKinoWeightHisto = zz.getHist(folder+"/weight_EWKinoPairPt")
    zz_AvgleptonWeightHisto = zz.getHist(folder+"/weight_leptonPairPt")
    wwg_AvgTopPtWeightHisto = wwgamma.getHist(folder+"/weight_topPt")
    wwg_AvgNIsrWeightHisto = wwgamma.getHist(folder+"/weight_nISR")
    wwg_AvgEWKinoWeightHisto = wwgamma.getHist(folder+"/weight_EWKinoPairPt")
    wwg_AvgleptonWeightHisto = wwgamma.getHist(folder+"/weight_leptonPairPt")
    wzg_AvgTopPtWeightHisto = wzgamma.getHist(folder+"/weight_topPt")
    wzg_AvgNIsrWeightHisto = wzgamma.getHist(folder+"/weight_nISR")
    wzg_AvgEWKinoWeightHisto = wzgamma.getHist(folder+"/weight_EWKinoPairPt")
    wzg_AvgleptonWeightHisto = wzgamma.getHist(folder+"/weight_leptonPairPt")
    dy_AvgTopPtWeightHisto = DYjetsNLO.getHist(folder+"/weight_topPt")
    dy_AvgNIsrWeightHisto = DYjetsNLO.getHist(folder+"/weight_nISR")
    dy_AvgEWKinoWeightHisto = DYjetsNLO.getHist(folder+"/weight_EWKinoPairPt")
    dy_AvgleptonWeightHisto = DYjetsNLO.getHist(folder+"/weight_leptonPairPt")
    wjets_AvgTopPtWeightHisto = wjets.getHist(folder+"/weight_topPt")
    wjets_AvgNIsrWeightHisto = wjets.getHist(folder+"/weight_nISR")
    wjets_AvgEWKinoWeightHisto = wjets.getHist(folder+"/weight_EWKinoPairPt")
    wjets_AvgleptonWeightHisto = wjets.getHist(folder+"/weight_leptonPairPt")
    tt_AvgTopPtWeightHisto = tt.getHist(folder+"/weight_topPt")
    tt_AvgNIsrWeightHisto = tt.getHist(folder+"/weight_nISR")
    tt_AvgEWKinoWeightHisto = tt.getHist(folder+"/weight_EWKinoPairPt")
    tt_AvgleptonWeightHisto = tt.getHist(folder+"/weight_leptonPairPt")
    tt080_AvgTopPtWeightHisto = tt.getHist(folder.replace("EM","080EM")+"/weight_topPt")
    tt080_AvgNIsrWeightHisto = tt.getHist(folder.replace("EM","080EM")+"/weight_nISR")
    tt080_AvgEWKinoWeightHisto = tt.getHist(folder.replace("EM","080EM")+"/weight_EWKinoPairPt")
    tt080_AvgleptonWeightHisto = tt.getHist(folder.replace("EM","080EM")+"/weight_leptonPairPt")
    tt80_AvgTopPtWeightHisto = tt.getHist(folder.replace("EM","80EM")+"/weight_topPt")
    tt80_AvgNIsrWeightHisto = tt.getHist(folder.replace("EM","80EM")+"/weight_nISR")
    tt80_AvgEWKinoWeightHisto = tt.getHist(folder.replace("EM","80EM")+"/weight_EWKinoPairPt")
    tt80_AvgleptonWeightHisto = tt.getHist(folder.replace("EM","80EM")+"/weight_leptonPairPt")
    singletop_AvgTopPtWeightHisto = singletop.getHist(folder+"/weight_topPt")
    singletop_AvgNIsrWeightHisto = singletop.getHist(folder+"/weight_nISR")
    singletop_AvgEWKinoWeightHisto = singletop.getHist(folder+"/weight_EWKinoPairPt")
    singletop_AvgleptonWeightHisto = singletop.getHist(folder+"/weight_leptonPairPt")
    wz_AvgTopPtWeightHisto = wz.getHist(folder+"/weight_topPt")
    wz_AvgNIsrWeightHisto = wz.getHist(folder+"/weight_nISR")
    wz_AvgEWKinoWeightHisto = wz.getHist(folder+"/weight_EWKinoPairPt")
    wz_AvgleptonWeightHisto = wz.getHist(folder+"/weight_leptonPairPt")
    ww_AvgTopPtWeightHisto = ww.getHist(folder+"/weight_topPt")
    ww_AvgNIsrWeightHisto = ww.getHist(folder+"/weight_nISR")
    ww_AvgEWKinoWeightHisto = ww.getHist(folder+"/weight_EWKinoPairPt")
    ww_AvgleptonWeightHisto = ww.getHist(folder+"/weight_leptonPairPt")
    zz4l_AvgTopPtWeightHisto = zz4l.getHist(folder+"/weight_topPt")
    zz4l_AvgNIsrWeightHisto = zz4l.getHist(folder+"/weight_nISR")
    zz4l_AvgEWKinoWeightHisto = zz4l.getHist(folder+"/weight_EWKinoPairPt")
    zz4l_AvgleptonWeightHisto = zz4l.getHist(folder+"/weight_leptonPairPt")
    wg_AvgTopPtWeightHisto = wgamma.getHist(folder+"/weight_topPt")
    wg_AvgNIsrWeightHisto = wgamma.getHist(folder+"/weight_nISR")
    wg_AvgEWKinoWeightHisto = wgamma.getHist(folder+"/weight_EWKinoPairPt")
    wg_AvgleptonWeightHisto = wgamma.getHist(folder+"/weight_leptonPairPt")
    
    histsToScale=[zgHist,ttgHist,ttg080Hist,ttg80Hist,zzHist,wwgHist,wzgHist,dyHist,wjetsHist,ttHist,tt080Hist,tt80Hist,singletopHist,wzHist,wwHist,zz4lHist,wgHist]
    topWeightHists=[zg_AvgTopPtWeightHisto,ttg_AvgTopPtWeightHisto,ttg080_AvgTopPtWeightHisto,ttg80_AvgTopPtWeightHisto,zz_AvgTopPtWeightHisto,wwg_AvgTopPtWeightHisto,wzg_AvgTopPtWeightHisto,dy_AvgTopPtWeightHisto,wjets_AvgTopPtWeightHisto,tt_AvgTopPtWeightHisto,tt080_AvgTopPtWeightHisto,tt80_AvgTopPtWeightHisto,singletop_AvgTopPtWeightHisto,wz_AvgTopPtWeightHisto,ww_AvgTopPtWeightHisto,zz4l_AvgTopPtWeightHisto,wg_AvgTopPtWeightHisto]
    nISRWeightHists=[zg_AvgNIsrWeightHisto,ttg_AvgNIsrWeightHisto,ttg080_AvgNIsrWeightHisto,ttg80_AvgNIsrWeightHisto,zz_AvgNIsrWeightHisto,wwg_AvgNIsrWeightHisto,wzg_AvgNIsrWeightHisto,dy_AvgNIsrWeightHisto,wjets_AvgNIsrWeightHisto,tt_AvgNIsrWeightHisto,tt080_AvgNIsrWeightHisto,tt80_AvgNIsrWeightHisto,singletop_AvgNIsrWeightHisto,wz_AvgNIsrWeightHisto,ww_AvgNIsrWeightHisto,zz4l_AvgNIsrWeightHisto,wg_AvgNIsrWeightHisto]
    EWKinoWeightHists=[zg_AvgEWKinoWeightHisto,ttg_AvgEWKinoWeightHisto,ttg080_AvgEWKinoWeightHisto,tt80_AvgEWKinoWeightHisto,zz_AvgEWKinoWeightHisto,wwg_AvgEWKinoWeightHisto,wzg_AvgEWKinoWeightHisto,dy_AvgEWKinoWeightHisto,wjets_AvgEWKinoWeightHisto,tt_AvgEWKinoWeightHisto,tt080_AvgEWKinoWeightHisto,tt80_AvgEWKinoWeightHisto,singletop_AvgEWKinoWeightHisto,wz_AvgEWKinoWeightHisto,ww_AvgEWKinoWeightHisto,zz4l_AvgEWKinoWeightHisto,wg_AvgEWKinoWeightHisto]
    leptonWeightHists=[zg_AvgleptonWeightHisto,ttg_AvgleptonWeightHisto,ttg080_AvgleptonWeightHisto,ttg80_AvgleptonWeightHisto,zz_AvgleptonWeightHisto,wwg_AvgleptonWeightHisto,wzg_AvgleptonWeightHisto,dy_AvgleptonWeightHisto,wjets_AvgleptonWeightHisto,tt_AvgleptonWeightHisto,tt080_AvgleptonWeightHisto,tt80_AvgleptonWeightHisto,singletop_AvgleptonWeightHisto,wz_AvgleptonWeightHisto,ww_AvgleptonWeightHisto,zz4l_AvgleptonWeightHisto,wg_AvgleptonWeightHisto]
      
    for i in range(len(histsToScale)):
        if topWeightHists[i].Integral>0.:
            histsToScale[i].Scale(topWeightHists[i].GetMean())
        if nISRWeightHists[i].Integral>0.:
            histsToScale[i].Scale(nISRWeightHists[i].GetMean())
        if EWKinoWeightHists[i].Integral>0.:
            histsToScale[i].Scale(EWKinoWeightHists[i].GetMean())
        if leptonWeightHists[i].Integral>0.:
            histsToScale[i].Scale(leptonWeightHists[i].GetMean())


    zgHist.SetLineColor(ROOT.kGreen-3)
    ttgHist.SetLineColor(ROOT.kRed+1)
    #ttg080Hist.SetLineColor(ROOT.kRed+1)
    #ttg80Hist.SetLineColor(ROOT.kRed+1)
    zzHist.SetLineColor(ROOT.kYellow)
    wwgHist.SetLineColor(ROOT.kCyan+3)
    wzgHist.SetLineColor(ROOT.kBlue+3)
    dyHist.SetLineColor(ROOT.kGreen+3)
    wjetsHist.SetLineColor(ROOT.kBlue-9)
    ttHist.SetLineColor(ROOT.kOrange+8)
    singletopHist.SetLineColor(ROOT.kOrange+4)
    wzHist.SetLineColor(ROOT.kAzure-5)
    wwHist.SetLineColor(ROOT.kCyan-2)
    zz4lHist.SetLineColor(ROOT.kOrange-2)
    wgHist.SetLineColor(ROOT.kRed+3)


    #scaleHist=aux.addHists(ttHist,ttgHist)
    scaleHist=aux.addHists(tt080Hist,tt80Hist,ttg080Hist,ttg80Hist)
    fixHist=aux.addHists(zgHist,zzHist,wwgHist,wzgHist,dyHist,wjetsHist,singletopHist,wzHist,wwHist,zz4lHist,wgHist)
    
    
    sf,sferr = calculateSFAndError(dataHist, scaleHist,fixHist)

    #ttHist.Scale(sf)
    #ttgHist.Scale(sf)
    #scaleHist.Scale(sf)


    mcSystUncert = 0.0 # SF, lumi, trigger
    zgSyst = aux.getSysHisto(zgHist, mcSystUncert)
    #ttgSyst = aux.getSysHisto(ttgHist, sferr)
    zzSyst = aux.getSysHisto(zzHist, mcSystUncert)
    wwgSyst = aux.getSysHisto(wwgHist, mcSystUncert)
    wzgSyst = aux.getSysHisto(wzgHist, mcSystUncert)
    dySyst = aux.getSysHisto(dyHist, mcSystUncert)
    wjetsSyst = aux.getSysHisto(wjetsHist, mcSystUncert)
    #ttSyst = aux.getSysHisto(ttHist, sferr)
    singletopSyst=aux.getSysHisto(singletopHist,mcSystUncert)
    wzSyst=aux.getSysHisto(wzHist,mcSystUncert)
    wwSyst=aux.getSysHisto(wwHist,mcSystUncert)
    zz4lSyst=aux.getSysHisto(zz4lHist,mcSystUncert)
    wgSyst=aux.getSysHisto(wgHist,mcSystUncert)

    cutBin=ttHist.FindBin(80)
    
    #ttg080Syst=aux.getSysHisto(ttg080Hist,0.04)
    #ttg80Syst=aux.getSysHisto(ttg80Hist,0.4)
    #tt080Syst=aux.getSysHisto(tt080Hist,0.04)
    #tt80Syst=aux.getSysHisto(tt80Hist,0.4)
    
    
    
    scaleSyst080=aux.getSysHisto(aux.addHists(ttg080Hist,tt080Hist),0.04)
    scaleSyst80=aux.getSysHisto(aux.addHists(ttg80Hist,tt80Hist),0.4)
    #scaleSyst=aux.addHists(tt080Syst,tt80Syst,ttg080Syst,ttg80Syst)
    scaleSyst=aux.addHists(scaleSyst080,scaleSyst80)
    #scaleSyst=aux.getSysHisto(aux.addHists(ttg80Hist,tt80Hist),0.4)
    
    #if "pt_g1" in name:
        #scaleSyst=aux.getSysHistoCut(scaleHist,0.04,0.4,cutBin)
    #else:
        ##scaleSyst=aux.getSysHisto(scaleHist,sferr)
        #scaleSyst=aux.getSysHisto(scaleHist,0.04)
            

    dataSyst= aux.getSysHisto(dataHist,mcSystUncert)

    totStat = aux.addHists(zgHist, ttgHist, zzHist, wwgHist, wzgHist,dyHist,wjetsHist,ttHist,singletopHist,wzHist,wwHist,zz4lHist,wgHist)
    totSyst = aux.addHists(zgSyst, zzSyst, wwgSyst, wzgSyst,dySyst,wjetsSyst,singletopSyst,wzSyst,wwSyst,zz4lSyst,wgSyst,scaleSyst)
    
    totUnc = aux.addHistUncert(totStat, totSyst)
    aux.drawOpt(totUnc, "totUnc")
    aux.drawOpt(totSyst, "sysUnc")
    
    
    signal2 = aux.stdHist(t5bbbbzg_1500_400, name, binning)
    signal1 = aux.stdHist(tching_400, name, binning)
    #signal3 = aux.stdHist(t6ttZg_600_300, name, binning)
    
    signal2_AvgTopPtWeightHisto = t5bbbbzg_1500_400.getHist(folder+"/weight_topPt")
    signal2_AvgNIsrWeightHisto = t5bbbbzg_1500_400.getHist(folder+"/weight_nISR")
    signal2_AvgEWKinoWeightHisto = t5bbbbzg_1500_400.getHist(folder+"/weight_EWKinoPairPt")
    signal2_AvgleptonWeightHisto = t5bbbbzg_1500_400.getHist(folder+"/weight_leptonPairPt")
    signal1_AvgTopPtWeightHisto = tching_600.getHist(folder+"/weight_topPt")
    signal1_AvgNIsrWeightHisto = tching_600.getHist(folder+"/weight_nISR")
    signal1_AvgEWKinoWeightHisto = tching_600.getHist(folder+"/weight_EWKinoPairPt")
    signal1_AvgleptonWeightHisto = tching_600.getHist(folder+"/weight_leptonPairPt")
    #signal3_AvgTopPtWeightHisto = t6ttZg_600_300.getHist(folder+"/weight_topPt")
    #signal3_AvgNIsrWeightHisto = t6ttZg_600_300.getHist(folder+"/weight_nISR")
    #signal3_AvgEWKinoWeightHisto = t6ttZg_600_300.getHist(folder+"/weight_EWKinoPairPt")
    #signal3_AvgleptonWeightHisto = t6ttZg_600_300.getHist(folder+"/weight_leptonPairPt")

    signal2.Scale((1./signal2_AvgTopPtWeightHisto.GetMean()))
    signal2.Scale((1./signal2_AvgNIsrWeightHisto.GetMean()))
    signal2.Scale((1./signal2_AvgEWKinoWeightHisto.GetMean()))
    signal2.Scale((1./signal2_AvgleptonWeightHisto.GetMean()))
    signal1.Scale((1./signal1_AvgTopPtWeightHisto.GetMean()))
    signal1.Scale((1./signal1_AvgNIsrWeightHisto.GetMean()))
    signal1.Scale((1./signal1_AvgEWKinoWeightHisto.GetMean()))
    signal1.Scale((1./signal1_AvgleptonWeightHisto.GetMean()))
    #signal3.Scale((1./signal3_AvgTopPtWeightHisto.GetMean()))
    #signal3.Scale((1./signal3_AvgNIsrWeightHisto.GetMean()))
    #signal3.Scale((1./signal3_AvgEWKinoWeightHisto.GetMean()))
    #signal3.Scale((1./signal3_AvgleptonWeightHisto.GetMean()))
    
    for h in signal1, signal2:
    #for h in signal1, signal2,signal3:
        aux.drawOpt(h, "signal")
    signal1.SetLineColor(ROOT.kBlue+3)
    signal2.SetLineColor(ROOT.kBlue+3)
    signal2.SetLineStyle(2)    
    #signal3.SetLineColor(ROOT.kRed+3)
    #signal3.SetLineStyle(2)    

    m.add(signal2, "T5bbbbZg")
    m.add(signal1, "TChiNg")
    #m.add(signal3, "T6ttZg")
    
    m.addStack(zgHist, "Z#gamma")
    m.addStack(ttgHist, "t#bar{t}#gamma")
    m.addStack(zzHist, "ZZ")
    m.addStack(wwgHist, "WW#gamma")
    m.addStack(wzgHist, "WZ#gamma")
    m.addStack(dyHist, "Drell-Yan/Z")
    m.addStack(wjetsHist, "W+jets")
    m.addStack(ttHist, "t#bar{t}")
    m.addStack(singletopHist, "single t")
    m.addStack(wwHist, "WW")
    m.addStack(wzHist, "WZ")
    m.addStack(zz4lHist, "ZZ(4l)")
    m.addStack(wgHist, "W#gamma")
    
    
    
    m.sortStackByIntegral()


    m.add(dataHist,"Data")
    m.add(totUnc, "Total uncertainty")
    m.add(totSyst, "syst. uncertainty")

    legInfo = "t#bar{t}(+#gamma) Control Region"
    m.leg.SetHeader(legInfo)
    #m.leg.SetY1(.56)
    #m.leg.SetX1(.56)
    #m.leg.SetX1(.46)
    #m.leg.SetX2(.99)
    #m.leg.SetX2(.89)


    m.Draw()

    r = ratio.Ratio("#scale[.9]{#lower[.24]{#splitline{Data/Pred.}{Bkg. frac.}}}", dataHist, totStat, sysHisto=totSyst)
    rMax = 1.5
    #r.draw(0., rMax, m.getStack(), True)
    r.draw(0., rMax, m.getStack())
    
    aux.Label(sim= False, status="Work in Progress")
    #aux.save(name, normal=False, changeMinMax=False)
    directory= "plots_CR_tt/"
    if not os.path.exists(directory):
        os.makedirs(directory)
    aux.save(name.replace("/","_"),folder=directory)
        
    return sf,sferr

def main():
    variables=["eta1","eta2","pt1","pt2","n_jets","n_vtx","phi1","phi2","m_ll","ht","n_photons","pt_g1","met","nElectrons","nMuons","deltaR1_g1","deltaR2_g1","deltaRLL"]
    #variables=["pt_g1"]
    bkgs=[DYjetsNLO,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma,singletop,zz4l,wz,zz,ww]
    #groups=["ControlRegionTT"]
    groups=["CRTT"]
    binnings_=binnings.copy()
    binnings_["pt_g1"]=frange(25,80,5)+frange(80,140,10)+frange(140,200,20)
    toSave={}
    toSave["EM"]={}
    for group in groups:
        for variable in variables:
            toSave["EM"][variable] = drawTTCR("EM",group+"EM/"+variable,binning=binnings_[variable],xTitle=labels[variable][0])
        if not os.path.exists("plots_CR_tt/factors"):
            os.makedirs("plots_CR_tt/factors")
        pkl.dump( toSave, open( "plots_CR_tt/factors/"+group+".pkl", "wb" ) )
        
        
if __name__=="__main__":
    main()        
        
        

