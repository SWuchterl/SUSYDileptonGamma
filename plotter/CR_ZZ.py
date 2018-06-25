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
    'pt1':              frange(20,150,10),
    #'pt2':              frange(20,80,5)+frange(80,160,20),
    'pt2':              frange(20.,100,10.),
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
    'ht':               frange(0., 1000.,50),
    #'met':               [0,25,50,75,100,150,190,230,500],
    #'met':               frange(0,49,10)+frange(50,100,25),
    #'met':               frange(0., 80.,10),
    #'met':               frange(0., 50.,10)+frange(60,120,20),
    'met':               frange(0., 120.,10),
    #'met':               frange(0., 100.,20),
    #'m_ll':             frange(50.,100.,10)+frange(100., 300.,20.),
    'm_ll':             frange(84.,98.,1),
    #'m_ll2':             frange(50.,100.,10)+frange(100., 300.,20.),
    'm_ll2':             frange(60., 120., 5.),
    'm_llll':             frange(60., 120., 5.),
    'm_llg':             frange(0,100,10)+frange(100,200,10)+frange(200,500,50),
    'pt_llg':             frange(0,100,10)+frange(100,200,10)+frange(200,500,50),
    'n_jets':           frange(0.,7.,1),
    'n_photons':         frange(0.,4.,1),
    'n_vtx':            frange(0.,40.,1),
    'pt_g1':            frange(20,100,10)+frange(100,150,25)+frange(150,250,50),
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
    'mTL3Met': frange(0.,300.,10),
    'Fakes': frange(0,5,1)
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

def drawCR(sampleNames, name,datasetToUse, binning=None, binningName="", xTitle=None, yTitle=None):
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()
    
    #folder= (name.split("/"))[0]
    folder= (name.split("/"))[0]+"/"+(name.split("/"))[1]+"/"+(name.split("/"))[2]
    
    style.divideByBinWidth = False
    
    style.minimumOne=True
    
    dataHist = aux.stdHist(datasetToUse, name, binning)
    aux.drawOpt(dataHist, "data")
    
    zgHist = aux.stdHist(zgamma, name, binning)
    ttgHist = aux.stdHist(ttgamma, name, binning)
    zzHist = aux.stdHist(zz, name, binning)
    wwgHist = aux.stdHist(wwgamma, name, binning)
    wzgHist = aux.stdHist(wzgamma, name, binning)
    dyHist = aux.stdHist(DYjetsNLO, name, binning)
    wjetsHist = aux.stdHist(wjets, name, binning)
    ttHist = aux.stdHist(tt, name, binning)
    singletopHist=aux.stdHist(singletop, name, binning)
    wzHist=aux.stdHist(wz, name, binning)
    wwHist=aux.stdHist(ww, name, binning)
    zz4lHist=aux.stdHist(zz4l, name, binning)
    wgHist=aux.stdHist(wgamma, name, binning)
    
    
    #Scaling
    zg_AvgTopPtWeightHisto = zgamma.getHist(folder+"/weight_topPt")
    zg_AvgNIsrWeightHisto = zgamma.getHist(folder+"/weight_nISR")
    zg_AvgEWKinoWeightHisto = zgamma.getHist(folder+"/weight_EWKinoPairPt")
    #zg_AvgleptonWeightHisto = zgamma.getHist(folder+"/weight_leptonPairPt")
    ttg_AvgTopPtWeightHisto = ttgamma.getHist(folder+"/weight_topPt")
    ttg_AvgNIsrWeightHisto = ttgamma.getHist(folder+"/weight_nISR")
    ttg_AvgEWKinoWeightHisto = ttgamma.getHist(folder+"/weight_EWKinoPairPt")
    #ttg_AvgleptonWeightHisto = ttgamma.getHist(folder+"/weight_leptonPairPt")
    zz_AvgTopPtWeightHisto = zz.getHist(folder+"/weight_topPt")
    zz_AvgNIsrWeightHisto = zz.getHist(folder+"/weight_nISR")
    zz_AvgEWKinoWeightHisto = zz.getHist(folder+"/weight_EWKinoPairPt")
    #zz_AvgleptonWeightHisto = zz.getHist(folder+"/weight_leptonPairPt")
    wwg_AvgTopPtWeightHisto = wwgamma.getHist(folder+"/weight_topPt")
    wwg_AvgNIsrWeightHisto = wwgamma.getHist(folder+"/weight_nISR")
    wwg_AvgEWKinoWeightHisto = wwgamma.getHist(folder+"/weight_EWKinoPairPt")
    #wwg_AvgleptonWeightHisto = wwgamma.getHist(folder+"/weight_leptonPairPt")
    wzg_AvgTopPtWeightHisto = wzgamma.getHist(folder+"/weight_topPt")
    wzg_AvgNIsrWeightHisto = wzgamma.getHist(folder+"/weight_nISR")
    wzg_AvgEWKinoWeightHisto = wzgamma.getHist(folder+"/weight_EWKinoPairPt")
    #wzg_AvgleptonWeightHisto = wzgamma.getHist(folder+"/weight_leptonPairPt")
    dy_AvgTopPtWeightHisto = DYjetsNLO.getHist(folder+"/weight_topPt")
    dy_AvgNIsrWeightHisto = DYjetsNLO.getHist(folder+"/weight_nISR")
    dy_AvgEWKinoWeightHisto = DYjetsNLO.getHist(folder+"/weight_EWKinoPairPt")
    #dy_AvgleptonWeightHisto = DYjetsNLO.getHist(folder+"/weight_leptonPairPt")
    wjets_AvgTopPtWeightHisto = wjets.getHist(folder+"/weight_topPt")
    wjets_AvgNIsrWeightHisto = wjets.getHist(folder+"/weight_nISR")
    wjets_AvgEWKinoWeightHisto = wjets.getHist(folder+"/weight_EWKinoPairPt")
    #wjets_AvgleptonWeightHisto = wjets.getHist(folder+"/weight_leptonPairPt")
    tt_AvgTopPtWeightHisto = tt.getHist(folder+"/weight_topPt")
    tt_AvgNIsrWeightHisto = tt.getHist(folder+"/weight_nISR")
    tt_AvgEWKinoWeightHisto = tt.getHist(folder+"/weight_EWKinoPairPt")
    #tt_AvgleptonWeightHisto = tt.getHist(folder+"/weight_leptonPairPt")
    singletop_AvgTopPtWeightHisto = singletop.getHist(folder+"/weight_topPt")
    singletop_AvgNIsrWeightHisto = singletop.getHist(folder+"/weight_nISR")
    singletop_AvgEWKinoWeightHisto = singletop.getHist(folder+"/weight_EWKinoPairPt")
    #singletop_AvgleptonWeightHisto = singletop.getHist(folder+"/weight_leptonPairPt")
    wz_AvgTopPtWeightHisto = wz.getHist(folder+"/weight_topPt")
    wz_AvgNIsrWeightHisto = wz.getHist(folder+"/weight_nISR")
    wz_AvgEWKinoWeightHisto = wz.getHist(folder+"/weight_EWKinoPairPt")
    #wz_AvgleptonWeightHisto = wz.getHist(folder+"/weight_leptonPairPt")
    ww_AvgTopPtWeightHisto = ww.getHist(folder+"/weight_topPt")
    ww_AvgNIsrWeightHisto = ww.getHist(folder+"/weight_nISR")
    ww_AvgEWKinoWeightHisto = ww.getHist(folder+"/weight_EWKinoPairPt")
    #ww_AvgleptonWeightHisto = ww.getHist(folder+"/weight_leptonPairPt")
    zz4l_AvgTopPtWeightHisto = zz4l.getHist(folder+"/weight_topPt")
    zz4l_AvgNIsrWeightHisto = zz4l.getHist(folder+"/weight_nISR")
    zz4l_AvgEWKinoWeightHisto = zz4l.getHist(folder+"/weight_EWKinoPairPt")
    #zz4l_AvgleptonWeightHisto = zz4l.getHist(folder+"/weight_leptonPairPt")
    wg_AvgTopPtWeightHisto = wgamma.getHist(folder+"/weight_topPt")
    wg_AvgNIsrWeightHisto = wgamma.getHist(folder+"/weight_nISR")
    wg_AvgEWKinoWeightHisto = wgamma.getHist(folder+"/weight_EWKinoPairPt")
    #wg_AvgleptonWeightHisto = wgamma.getHist(folder+"/weight_leptonPairPt")
    
    histsToScale=[zgHist,ttgHist,zzHist,wwgHist,wzgHist,dyHist,wjetsHist,ttHist,singletopHist,wzHist,wwHist,zz4lHist,wgHist]
    topWeightHists=[zg_AvgTopPtWeightHisto,ttg_AvgTopPtWeightHisto,zz_AvgTopPtWeightHisto,wwg_AvgTopPtWeightHisto,wzg_AvgTopPtWeightHisto,dy_AvgTopPtWeightHisto,wjets_AvgTopPtWeightHisto,tt_AvgTopPtWeightHisto,singletop_AvgTopPtWeightHisto,wz_AvgTopPtWeightHisto,ww_AvgTopPtWeightHisto,zz4l_AvgTopPtWeightHisto,wg_AvgTopPtWeightHisto]
    nISRWeightHists=[zg_AvgNIsrWeightHisto,ttg_AvgNIsrWeightHisto,zz_AvgNIsrWeightHisto,wwg_AvgNIsrWeightHisto,wzg_AvgNIsrWeightHisto,dy_AvgNIsrWeightHisto,wjets_AvgNIsrWeightHisto,tt_AvgNIsrWeightHisto,singletop_AvgNIsrWeightHisto,wz_AvgNIsrWeightHisto,ww_AvgNIsrWeightHisto,zz4l_AvgNIsrWeightHisto,wg_AvgNIsrWeightHisto]
    EWKinoWeightHists=[zg_AvgEWKinoWeightHisto,ttg_AvgEWKinoWeightHisto,zz_AvgEWKinoWeightHisto,wwg_AvgEWKinoWeightHisto,wzg_AvgEWKinoWeightHisto,dy_AvgEWKinoWeightHisto,wjets_AvgEWKinoWeightHisto,tt_AvgEWKinoWeightHisto,singletop_AvgEWKinoWeightHisto,wz_AvgEWKinoWeightHisto,ww_AvgEWKinoWeightHisto,zz4l_AvgEWKinoWeightHisto,wg_AvgEWKinoWeightHisto]
    #leptonWeightHists=[zg_AvgleptonWeightHisto,ttg_AvgleptonWeightHisto,zz_AvgleptonWeightHisto,wwg_AvgleptonWeightHisto,wzg_AvgleptonWeightHisto,dy_AvgleptonWeightHisto,wjets_AvgleptonWeightHisto,tt_AvgleptonWeightHisto,singletop_AvgleptonWeightHisto,wz_AvgleptonWeightHisto,ww_AvgleptonWeightHisto,zz4l_AvgleptonWeightHisto,wg_AvgleptonWeightHisto]
    
    for i in range(len(histsToScale)):
        if topWeightHists[i].Integral>0.:
            histsToScale[i].Scale(topWeightHists[i].GetMean())
        if nISRWeightHists[i].Integral>0.:
            histsToScale[i].Scale(nISRWeightHists[i].GetMean())
        if EWKinoWeightHists[i].Integral>0.:
            histsToScale[i].Scale(EWKinoWeightHists[i].GetMean())
        #if leptonWeightHists[i].Integral>0.:
            #histsToScale[i].Scale(leptonWeightHists[i].GetMean())


    zgHist.SetLineColor(ROOT.kGreen-3)
    ttgHist.SetLineColor(ROOT.kRed+1)
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


    scaleHist=aux.addHists(zz4lHist)
    fixHist=aux.addHists(ttgHist,wzHist,wwgHist,wzgHist,ttHist,wjetsHist,singletopHist,dyHist,wwHist,zzHist,wgHist,zgHist)
    
    
    sf,sferr = calculateSFAndError(dataHist, scaleHist,fixHist)

    zz4lHist.Scale(sf)
    scaleHist.Scale(sf)


    mcSystUncert = 0.0 # SF, lumi, trigger
    zgSyst = aux.getSysHisto(zgHist, mcSystUncert)
    ttgSyst = aux.getSysHisto(ttgHist, mcSystUncert)
    zzSyst = aux.getSysHisto(zzHist, mcSystUncert)
    wwgSyst = aux.getSysHisto(wwgHist, mcSystUncert)
    wzgSyst = aux.getSysHisto(wzgHist, mcSystUncert)
    dySyst = aux.getSysHisto(dyHist, mcSystUncert)
    wjetsSyst = aux.getSysHisto(wjetsHist, mcSystUncert)
    ttSyst = aux.getSysHisto(ttHist, mcSystUncert)
    singletopSyst=aux.getSysHisto(singletopHist,mcSystUncert)
    wzSyst=aux.getSysHisto(wzHist,mcSystUncert)
    wwSyst=aux.getSysHisto(wwHist,mcSystUncert)
    #zz4lSyst=aux.getSysHisto(zz4lHist,mcSystUncert)
    wgSyst=aux.getSysHisto(wgHist,mcSystUncert)

    cutBin=ttHist.FindBin(80)
    #
    #
    #if "pt_g1" in name:
        #scaleSyst=aux.getSysHistoCut(scaleHist,0.04,0.4,cutBin)
    #else:
    scaleSyst=aux.getSysHisto(scaleHist,sferr)
            

    dataSyst= aux.getSysHisto(dataHist,mcSystUncert)

    totStat = aux.addHists(zgHist, ttgHist, zzHist, wwgHist, wzgHist,dyHist,wjetsHist,ttHist,singletopHist,wzHist,wwHist,zz4lHist,wgHist)
    totSyst = aux.addHists(ttgSyst, zzSyst, wwgSyst, wzgSyst,ttSyst,wjetsSyst,singletopSyst,wwSyst,wzSyst,dySyst,zgSyst,scaleSyst)
    
    totUnc = aux.addHistUncert(totStat, totSyst)
    aux.drawOpt(totUnc, "totUnc")
    aux.drawOpt(totSyst, "sysUnc")
    
    signal2 = aux.stdHist(t5bbbbzg_1500_400, name, binning)
    signal1 = aux.stdHist(tching_400, name, binning)
    signal1 = aux.stdHist(tching_400, name, binning)
    
    signal2_AvgTopPtWeightHisto = t5bbbbzg_1500_400.getHist(folder+"/weight_topPt")
    signal2_AvgNIsrWeightHisto = t5bbbbzg_1500_400.getHist(folder+"/weight_nISR")
    signal2_AvgEWKinoWeightHisto = t5bbbbzg_1500_400.getHist(folder+"/weight_EWKinoPairPt")
    #signal2_AvgleptonWeightHisto = t5bbbbzg_1500_400.getHist(folder+"/weight_leptonPairPt")
    signal1_AvgTopPtWeightHisto = tching_600.getHist(folder+"/weight_topPt")
    signal1_AvgNIsrWeightHisto = tching_600.getHist(folder+"/weight_nISR")
    signal1_AvgEWKinoWeightHisto = tching_600.getHist(folder+"/weight_EWKinoPairPt")
    #signal1_AvgleptonWeightHisto = tching_600.getHist(folder+"/weight_leptonPairPt")

    signal2.Scale((1./signal2_AvgTopPtWeightHisto.GetMean()))
    signal2.Scale((1./signal2_AvgNIsrWeightHisto.GetMean()))
    signal2.Scale((1./signal2_AvgEWKinoWeightHisto.GetMean()))
    #signal2.Scale((1./signal2_AvgleptonWeightHisto.GetMean()))
    signal1.Scale((1./signal1_AvgTopPtWeightHisto.GetMean()))
    signal1.Scale((1./signal1_AvgNIsrWeightHisto.GetMean()))
    signal1.Scale((1./signal1_AvgEWKinoWeightHisto.GetMean()))
    #signal1.Scale((1./signal1_AvgleptonWeightHisto.GetMean()))
    
    for h in signal1, signal2:
        aux.drawOpt(h, "signal")
    signal1.SetLineColor(ROOT.kBlue+3)
    signal2.SetLineColor(ROOT.kBlue+3)
    signal2.SetLineStyle(2)    

    m.add(signal2, "T5bbbbZg")
    m.add(signal1, "TChiNg")
    
    
    
    m.addStack(zgHist, "Z#gamma")
    m.addStack(ttgHist, "t#bar{t}#gamma")
    m.addStack(zzHist, "ZZ(#rightarrowll#nu#nu)")
    m.addStack(wwgHist, "WW#gamma")
    m.addStack(wzgHist, "WZ#gamma")
    m.addStack(dyHist, "Drell-Yan/Z")
    m.addStack(wjetsHist, "W+jets")
    m.addStack(ttHist, "t#bar{t}")
    m.addStack(singletopHist, "single t")
    m.addStack(wwHist, "WW")
    m.addStack(wzHist, "WZ")
    m.addStack(zz4lHist, "ZZ(#rightarrow4l)")
    m.addStack(wgHist, "W#gamma")
    
    
    
    m.sortStackByIntegral()


    m.add(dataHist,"Data")
    m.add(totUnc, "Total uncertainty")
    m.add(totSyst, "syst. uncertainty")

    legInfo = "ZZ Control Region"
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
    directory= "plots_CR_zz/"
    if not os.path.exists(directory):
        os.makedirs(directory)
    aux.save(name.replace("/","_"),folder=directory)
        
    return sf,sferr

def main():
    #variables=["eta1","eta2","pt1","pt2","n_jets","n_vtx","phi1","phi2","m_ll","ht","n_photons","met","nElectrons","nMuons","deltaRLL"]
    #variables=["eta1","eta2","pt1","pt2","n_jets","n_vtx","phi1","phi2","m_ll","m_ll2","m_llll","ht","n_photons","met","nElectrons","nMuons"]
    variables=["eta1","eta2","pt1","pt2","n_jets","phi1","phi2","m_ll","m_ll2","ht","n_photons","met"]
    #variables=["pt_g1"]
    bkgs=[DYjetsNLO,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma,singletop,zz4l,wz,zz,ww]
    #groups=["ControlRegionZZ"]
    groups=["CRZZ"]
    binnings_=binnings.copy()
    #binnings_["pt_g1"]=frange(25,80,5)+frange(80,140,10)+frange(140,200,20)
    #binnings_["eta1"]=frange(0., 2.61, 0.4);
    #binnings_["eta2"]=frange(0., 2.61, 0.4);
    #binnings_["eta3"]=frange(0., 2.61, 0.4);
    #binnings_["eta4"]=frange(0., 2.61, 0.4);
    #binnings_["phi1"]=frange(0., 3.51, 0.25);
    #binnings_["phi2"]=frange(0., 3.51, 0.25);
    #binnings_["phi3"]=frange(0., 3.51, 0.25);
    #binnings_["phi4"]=frange(0., 3.51, 0.25);
    #binnings_["m_ll"]=frange(80., 110., 1.);
    #binnings_["m_ll2"]=frange(60., 120., 5.);
    #binnings_["m_llll"]=frange(100., 250., 5.);
    #binnings_["pt1"]=frange(25,160,10);
    #binnings_["pt2"]=frange(25,160,10);
    #binnings_["pt3"]=frange(25,160,10);
    #binnings_["pt4"]=frange(25,160,10);
    toSave={}
    toSave["EE"]={}
    toSave["MM"]={}
    toSave["LL"]={}
    for group in groups:
        for variable in variables:
            #toSave["EE"][variable] = drawCR("EE",group+"EE/"+variable,dataDoubleEG,binning=binnings_[variable],xTitle=labels[variable][0])
            #toSave["MM"][variable] = drawCR("MM",group+"MM/"+variable,dataDoubleMuon,binning=binnings_[variable],xTitle=labels[variable][0])
            toSave["EE"][variable] = drawCR("EE",group+"/EE/nom/"+variable,dataLL,binning=binnings_[variable],xTitle=labels[variable][0])
            toSave["MM"][variable] = drawCR("MM",group+"/MM/nom/"+variable,dataLL,binning=binnings_[variable],xTitle=labels[variable][0])
            toSave["LL"][variable] = drawCR("LL",group+"/LL/nom/"+variable,dataLL,binning=binnings_[variable],xTitle=labels[variable][0])
        if not os.path.exists("plots_CR_zz/factors"):
            os.makedirs("plots_CR_zz/factors")
        pkl.dump( toSave, open( "plots_CR_zz/factors/"+group+".pkl", "wb" ) )
        
        
        
if __name__=="__main__":
    main()                
        

