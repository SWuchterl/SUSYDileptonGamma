from dataMC import labels,frange,frangeN

import ROOT
from ROOT import *
from array import array
from include import *
import numpy as np
import pickle as pkl
import CR_tt
binnings = {
    'pt1':              frange(0,100,10)+frange(100,200,25)+range(200,350,50),
    'pt2':              frange(0,100,10)+frange(100,200,25),
    'pt3':              range(0, 200, 10),
    'pt4':              range(0, 200, 10),
    'eta1':             frange(0., 2.405, 0.1),
    'eta2':             frange(0., 2.4, 0.1),
    'eta3':             frange(0., 2.4, 0.1),
    'eta4':             frange(0., 2.4, 0.1),
    'phi1':             frange(0., 3.50, 0.1),
    'phi2':             frange(0., 3.50, 0.1),
    'phi3':             frange(0., 3.50, 0.1),
    'phi4':             frange(0., 3.50, 0.1),
    'ht':               frange(0., 1000.,50),
    #'met':               [0,25,50,75,100,150,190,230,500],
    'met':               [75,100,150,190,230,500],
    'm_ll':             frange(50.,100.,10)+frange(100., 300.,20.),
    'm_ll2':             frange(50.,100.,10)+frange(100., 300.,20.),
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
    'mTL3Met': frange(0.,300.,10)
}





def calculateChiSquared(dataHist,mcHists):
    obs=[]
    expected=[]
    #print dataHist,mcHists
    #for bin in range(dataHist.GetNbinsX()+1):
    for bin in range(dataHist.GetNbinsX()):
        bin=bin+1
        bw = dataHist.GetBinWidth(bin) if style.divideByBinWidth else 1.
        #bw2 = mcHists[0].GetBinWidth(bin) if style.divideByBinWidth else 1.
        
        #print dataHist.GetBinLowEdge(bin)
        
        obs.append(dataHist.GetBinContent(bin)*bw)
        expectedTemp=0.
        for histo in mcHists:
            expectedTemp+=(histo.GetBinContent(bin)*bw)
            #print histo,histo.GetBinContent(bin)*bw
        expected.append(expectedTemp)
    degOFr=dataHist.GetNbinsX() - 1.
    chiq2=0.
    #print "obs",obs
    #print "exp", expected
    for i in range(len(obs)):
        chiq2+=((obs[i]-expected[i])**2./(expected[i]))
    #chiq2=chiq2/(degOFr)
    #print chiq2
    return [chiq2,degOFr]


def divideDatasetIntegrals( numerator, denominator, name ):
    numMerged = sum(numerator)
    h_num = numMerged.getHist( name )
    num = h_num.Integral(0,-1)
    denMerged = sum(denominator)
    h_den = denMerged.getHist( name )
    den = h_den.Integral(0,-1)
    return num/den if den else 1.
    
def calculateBeta( data, ScaledTT,toScaleTTG,fix):
    num_dataErr=(ROOT.Double(0))
    num_data = data.IntegralAndError(0,-1,num_dataErr)
    
    den_ScaledTTErr=(ROOT.Double(0))
    den_ScaledTT = ScaledTT.IntegralAndError(0,-1,den_ScaledTTErr)
    
    den_toScaleTTGErr=(ROOT.Double(0))
    den_toScaleTTG = toScaleTTG.IntegralAndError(0,-1,den_toScaleTTGErr)

    add_fixError=ROOT.Double(0)
    add_fix = fix.IntegralAndError(0,-1,add_fixError)

    beta = (num_data-add_fix-den_ScaledTT)/(den_toScaleTTG)
    betaErr =np.sqrt( (num_dataErr/den_toScaleTTG)**2. +(add_fixError/den_toScaleTTG)**2. +(den_ScaledTTErr/den_toScaleTTG)**2. +(den_toScaleTTGErr*(num_data-add_fix-den_ScaledTT)/(den_toScaleTTG)**2.)**2. )
    
    if beta<0.:
        beta=0.0001
    
    return [beta,betaErr/beta] if den_toScaleTTG else [1.,0.]
    
def calculateSFAndError( numerator_data, denominator_toScale,additional_fix, name ):
    numMerged_data = sum(numerator_data)
    h_num_data = numMerged_data.getHist( name )
    num_dataErr=(ROOT.Double(0))
    num_data = h_num_data.IntegralAndError(0,-1,num_dataErr)
    
    denMerged_toScale = sum(denominator_toScale)
    h_den_toScale = denMerged_toScale.getHist( name )
    den_toScaleErr=(ROOT.Double(0))
    den_toScale = h_den_toScale.IntegralAndError(0,-1,den_toScaleErr)
    
    addMerged_fix = sum(additional_fix)
    h_add_fix = addMerged_fix.getHist( name )
    add_fixError=ROOT.Double(0)
    add_fix = h_add_fix.IntegralAndError(0,-1,add_fixError)
    
    alpha = (num_data-add_fix)/(den_toScale)
    
    alphaErr =np.sqrt( (num_dataErr/den_toScale)**2. +(add_fixError/den_toScale)**2. +(den_toScaleErr*(num_data-add_fix)/(den_toScale)**2.)**2. )
    
    return [alpha,alphaErr/alpha] if den_toScale else [1.,0.]


def drawSameHistogram( sampleNames, name, bkg=[], additional=[], binning=None, binningName="", scaleToData=True, xTitle=None, yTitle=None,toScale=[],SF=1.):
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()
    
    style.divideByBinWidth = False
    style.minimumOne=True
        
    scale = SF
    scaleErr=0.
    
    folder="CRTT/EM/nom"
    #print name,folder
    
    if __name__=="__main__":
        dataHist = aux.stdHist(dataLL, name, binning)
    else:
        dataHist = aux.stdHist(dataLL,"CRTT/EM/nom/"+(name.split("/"))[-1],binning)
    aux.drawOpt(dataHist, "data")
    
    zgHist = aux.stdHist(zgamma, name, binning)
    ttgHist = aux.stdHist(ttgamma, name, binning)
    ttg080Hist = aux.stdHist(ttgamma, name.replace("CRTT","CRTT080"), binning)
    ttg80Hist = aux.stdHist(ttgamma, name.replace("CRTT","CRTT80"), binning)
    zzHist = aux.stdHist(zz, name, binning)
    wwgHist = aux.stdHist(wwgamma, name, binning)
    wzgHist = aux.stdHist(wzgamma, name, binning)
    dyHist = aux.stdHist(DYjetsNLO, name, binning)
    wjetsHist = aux.stdHist(wjets, name, binning)
    ttHist = aux.stdHist(tt, name, binning)
    tt080Hist = aux.stdHist(tt, name.replace("CRTT","CRTT080"), binning)
    tt80Hist = aux.stdHist(tt, name.replace("CRTT","CRTT80"), binning)
    if __name__=="__main__":
        singletopHist=aux.stdHist(singletop, name, binning)
    else:
        singletopHist = aux.stdHist(singletop,"CRTT/EM/nom/"+(name.split("/"))[-1],binning)
    #singletopHist=aux.stdHist(singletop, name, binning)
    
    wzHist=aux.stdHist(wz, name, binning)
    wwHist=aux.stdHist(ww, name, binning)
    zz4lHist=aux.stdHist(zz4l, name, binning)
    wgHist=aux.stdHist(wgamma, name, binning)
    
    #print folder+"/weight_topPt"
    #Scaling
    zg_AvgTopPtWeightHisto = zgamma.getHist(folder+"/weight_topPt")
    zg_AvgNIsrWeightHisto = zgamma.getHist(folder+"/weight_nISR")
    zg_AvgEWKinoWeightHisto = zgamma.getHist(folder+"/weight_EWKinoPairPt")
    zg_PDFWeightHisto = zgamma.getHist(folder+"/weight_PDF")
    ttg_AvgTopPtWeightHisto = ttgamma.getHist(folder+"/weight_topPt")
    ttg_AvgNIsrWeightHisto = ttgamma.getHist(folder+"/weight_nISR")
    ttg_AvgEWKinoWeightHisto = ttgamma.getHist(folder+"/weight_EWKinoPairPt")
    ttg_PDFWeightHisto = ttgamma.getHist(folder+"/weight_PDF")
    ttg080_AvgTopPtWeightHisto = ttgamma.getHist(folder.replace("CRTT","CRTT080")+"/weight_topPt")
    ttg080_AvgNIsrWeightHisto = ttgamma.getHist(folder.replace("CRTT","CRTT080")+"/weight_nISR")
    ttg080_AvgEWKinoWeightHisto = ttgamma.getHist(folder.replace("CRTT","CRTT080")+"/weight_EWKinoPairPt")
    ttg080_PDFWeightHisto = ttgamma.getHist(folder.replace("CRTT","CRTT080")+"/weight_PDF")
    ttg80_AvgTopPtWeightHisto = ttgamma.getHist(folder.replace("CRTT","CRTT80")+"/weight_topPt")
    ttg80_AvgNIsrWeightHisto = ttgamma.getHist(folder.replace("CRTT","CRTT80")+"/weight_nISR")
    ttg80_AvgEWKinoWeightHisto = ttgamma.getHist(folder.replace("CRTT","CRTT80")+"/weight_EWKinoPairPt")
    ttg80_PDFWeightHisto = ttgamma.getHist(folder.replace("CRTT","CRTT80")+"/weight_PDF")
    zz_AvgTopPtWeightHisto = zz.getHist(folder+"/weight_topPt")
    zz_AvgNIsrWeightHisto = zz.getHist(folder+"/weight_nISR")
    zz_AvgEWKinoWeightHisto = zz.getHist(folder+"/weight_EWKinoPairPt")
    zz_PDFWeightHisto = zz.getHist(folder+"/weight_PDF")
    wwg_AvgTopPtWeightHisto = wwgamma.getHist(folder+"/weight_topPt")
    wwg_AvgNIsrWeightHisto = wwgamma.getHist(folder+"/weight_nISR")
    wwg_AvgEWKinoWeightHisto = wwgamma.getHist(folder+"/weight_EWKinoPairPt")
    wwg_PDFWeightHisto = wwgamma.getHist(folder+"/weight_PDF")
    wzg_AvgTopPtWeightHisto = wzgamma.getHist(folder+"/weight_topPt")
    wzg_AvgNIsrWeightHisto = wzgamma.getHist(folder+"/weight_nISR")
    wzg_AvgEWKinoWeightHisto = wzgamma.getHist(folder+"/weight_EWKinoPairPt")
    wzg_PDFWeightHisto = wzgamma.getHist(folder+"/weight_PDF")
    dy_AvgTopPtWeightHisto = DYjetsNLO.getHist(folder+"/weight_topPt")
    dy_AvgNIsrWeightHisto = DYjetsNLO.getHist(folder+"/weight_nISR")
    dy_AvgEWKinoWeightHisto = DYjetsNLO.getHist(folder+"/weight_EWKinoPairPt")
    dy_PDFWeightHisto = DYjetsNLO.getHist(folder+"/weight_PDF")
    wjets_AvgTopPtWeightHisto = wjets.getHist(folder+"/weight_topPt")
    wjets_AvgNIsrWeightHisto = wjets.getHist(folder+"/weight_nISR")
    wjets_AvgEWKinoWeightHisto = wjets.getHist(folder+"/weight_EWKinoPairPt")
    wjets_PDFWeightHisto = wjets.getHist(folder+"/weight_PDF")
    tt_AvgTopPtWeightHisto = tt.getHist(folder+"/weight_topPt")
    tt_AvgNIsrWeightHisto = tt.getHist(folder+"/weight_nISR")
    tt_AvgEWKinoWeightHisto = tt.getHist(folder+"/weight_EWKinoPairPt")
    tt_PDFWeightHisto = tt.getHist(folder+"/weight_PDF")
    tt080_AvgTopPtWeightHisto = tt.getHist(folder.replace("CRTT","CRTT080")+"/weight_topPt")
    tt080_AvgNIsrWeightHisto = tt.getHist(folder.replace("CRTT","CRTT080")+"/weight_nISR")
    tt080_AvgEWKinoWeightHisto = tt.getHist(folder.replace("CRTT","CRTT080")+"/weight_EWKinoPairPt")
    tt080_PDFWeightHisto = tt.getHist(folder.replace("CRTT","CRTT080")+"/weight_PDF")
    tt80_AvgTopPtWeightHisto = tt.getHist(folder.replace("CRTT","CRTT80")+"/weight_topPt")
    tt80_AvgNIsrWeightHisto = tt.getHist(folder.replace("CRTT","CRTT80")+"/weight_nISR")
    tt80_AvgEWKinoWeightHisto = tt.getHist(folder.replace("CRTT","CRTT80")+"/weight_EWKinoPairPt")
    tt80_PDFWeightHisto = tt.getHist(folder.replace("CRTT","CRTT80")+"/weight_PDF")
    
    #singletop_AvgTopPtWeightHisto = singletop.getHist(folder+"/weight_topPt")
    #singletop_AvgNIsrWeightHisto = singletop.getHist(folder+"/weight_nISR")
    #singletop_AvgEWKinoWeightHisto = singletop.getHist(folder+"/weight_EWKinoPairPt")
    #singletop_PDFWeightHisto = singletop.getHist(folder+"/weight_PDF")
    if __name__=="__main__":
        singletop_AvgTopPtWeightHisto = singletop.getHist(folder+"/weight_topPt")
        singletop_AvgNIsrWeightHisto = singletop.getHist(folder+"/weight_nISR")
        singletop_AvgEWKinoWeightHisto = singletop.getHist(folder+"/weight_EWKinoPairPt")
        singletop_PDFWeightHisto = singletop.getHist(folder+"/weight_PDF")
    else:
        singletop_AvgTopPtWeightHisto = singletop.getHist("CRTT/EM/nom/"+"/weight_topPt")
        singletop_AvgNIsrWeightHisto = singletop.getHist("CRTT/EM/nom/"+"/weight_nISR")
        singletop_AvgEWKinoWeightHisto = singletop.getHist("CRTT/EM/nom/"+"/weight_EWKinoPairPt")
        singletop_PDFWeightHisto = singletop.getHist("CRTT/EM/nom/"+"/weight_PDF")      
          
    wz_AvgTopPtWeightHisto = wz.getHist(folder+"/weight_topPt")
    wz_AvgNIsrWeightHisto = wz.getHist(folder+"/weight_nISR")
    wz_AvgEWKinoWeightHisto = wz.getHist(folder+"/weight_EWKinoPairPt")
    wz_PDFWeightHisto = wz.getHist(folder+"/weight_PDF")
    ww_AvgTopPtWeightHisto = ww.getHist(folder+"/weight_topPt")
    ww_AvgNIsrWeightHisto = ww.getHist(folder+"/weight_nISR")
    ww_AvgEWKinoWeightHisto = ww.getHist(folder+"/weight_EWKinoPairPt")
    ww_PDFWeightHisto = ww.getHist(folder+"/weight_PDF")
    zz4l_AvgTopPtWeightHisto = zz4l.getHist(folder+"/weight_topPt")
    zz4l_AvgNIsrWeightHisto = zz4l.getHist(folder+"/weight_nISR")
    zz4l_AvgEWKinoWeightHisto = zz4l.getHist(folder+"/weight_EWKinoPairPt")
    zz4l_PDFWeightHisto = zz4l.getHist(folder+"/weight_PDF")
    wg_AvgTopPtWeightHisto = wgamma.getHist(folder+"/weight_topPt")
    wg_AvgNIsrWeightHisto = wgamma.getHist(folder+"/weight_nISR")
    wg_AvgEWKinoWeightHisto = wgamma.getHist(folder+"/weight_EWKinoPairPt")
    wg_PDFWeightHisto = wgamma.getHist(folder+"/weight_PDF")
    
    histsToScale=[zgHist,ttgHist,ttg080Hist,ttg80Hist,zzHist,wwgHist,wzgHist,dyHist,wjetsHist,ttHist,tt080Hist,tt80Hist,singletopHist,wzHist,wwHist,zz4lHist,wgHist]
    topWeightHists=[zg_AvgTopPtWeightHisto,ttg_AvgTopPtWeightHisto,ttg080_AvgTopPtWeightHisto,ttg80_AvgTopPtWeightHisto,zz_AvgTopPtWeightHisto,wwg_AvgTopPtWeightHisto,wzg_AvgTopPtWeightHisto,dy_AvgTopPtWeightHisto,wjets_AvgTopPtWeightHisto,tt_AvgTopPtWeightHisto,tt080_AvgTopPtWeightHisto,tt80_AvgTopPtWeightHisto,singletop_AvgTopPtWeightHisto,wz_AvgTopPtWeightHisto,ww_AvgTopPtWeightHisto,zz4l_AvgTopPtWeightHisto,wg_AvgTopPtWeightHisto]
    nISRWeightHists=[zg_AvgNIsrWeightHisto,ttg_AvgNIsrWeightHisto,ttg080_AvgNIsrWeightHisto,ttg80_AvgNIsrWeightHisto,zz_AvgNIsrWeightHisto,wwg_AvgNIsrWeightHisto,wzg_AvgNIsrWeightHisto,dy_AvgNIsrWeightHisto,wjets_AvgNIsrWeightHisto,tt_AvgNIsrWeightHisto,tt080_AvgNIsrWeightHisto,tt80_AvgNIsrWeightHisto,singletop_AvgNIsrWeightHisto,wz_AvgNIsrWeightHisto,ww_AvgNIsrWeightHisto,zz4l_AvgNIsrWeightHisto,wg_AvgNIsrWeightHisto]
    EWKinoWeightHists=[zg_AvgEWKinoWeightHisto,ttg_AvgEWKinoWeightHisto,ttg080_AvgEWKinoWeightHisto,tt80_AvgEWKinoWeightHisto,zz_AvgEWKinoWeightHisto,wwg_AvgEWKinoWeightHisto,wzg_AvgEWKinoWeightHisto,dy_AvgEWKinoWeightHisto,wjets_AvgEWKinoWeightHisto,tt_AvgEWKinoWeightHisto,tt080_AvgEWKinoWeightHisto,tt80_AvgEWKinoWeightHisto,singletop_AvgEWKinoWeightHisto,wz_AvgEWKinoWeightHisto,ww_AvgEWKinoWeightHisto,zz4l_AvgEWKinoWeightHisto,wg_AvgEWKinoWeightHisto]
    PDFWeightHists=[zg_PDFWeightHisto,ttg_PDFWeightHisto,ttg080_PDFWeightHisto,ttg80_PDFWeightHisto,zz_PDFWeightHisto,wwg_PDFWeightHisto,wzg_PDFWeightHisto,dy_PDFWeightHisto,wjets_PDFWeightHisto,tt_PDFWeightHisto,tt080_PDFWeightHisto,tt80_PDFWeightHisto,singletop_PDFWeightHisto,wz_PDFWeightHisto,ww_PDFWeightHisto,zz4l_PDFWeightHisto,wg_PDFWeightHisto]
      
    for i in range(len(histsToScale)):
        if topWeightHists[i].Integral()>0.:
            histsToScale[i].Scale(1./topWeightHists[i].GetMean())
        if nISRWeightHists[i].Integral()>0.:
            histsToScale[i].Scale(1./nISRWeightHists[i].GetMean())
        if EWKinoWeightHists[i].Integral()>0.:
            histsToScale[i].Scale(1./EWKinoWeightHists[i].GetMean())
        if PDFWeightHists[i].Integral()>0.:
            histsToScale[i].Scale(1./PDFWeightHists[i].GetMean())


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
    
    
    #now scale tt with alpha=SF
    #ttHist.Scale(SF)
    ttgHist.Scale(SF)
    
    fixHist=aux.addHists(zgHist,zzHist,wwgHist,wzgHist,dyHist,wjetsHist,singletopHist,wzHist,wwHist,zz4lHist,wgHist)
    
    #now calculate beta (alpha for scaling tt, beta for ttg, norm fixed)
    #SF2,SF2Err=calculateBeta(dataHist,ttHist,ttgHist,fixHist)
    SF2,SF2Err=calculateBeta(dataHist,ttgHist,ttHist,fixHist)
    
    #now scale ttg with beta
    #ttgHist.Scale(SF2)
    ttHist.Scale(SF2)
    #ttgHist.Scale(SF)
    #and calculate chi2
    
    totalData=dataHist
    #totalMC=aux.addHists(ttHist,ttgHist,fixHist)
    totalMC=[ttHist,ttgHist,zgHist,zzHist,wwgHist,wzgHist,dyHist,wjetsHist,singletopHist,wzHist,wwHist,zz4lHist,wgHist]
    
    #print totalData,totalMC
    #SF2=SF
    return(SF,SF2,calculateChiSquared(totalData,totalMC))
        

      
            
def drawTT(binningToUse_=CR_tt.binnings.copy(),addName=""):
    #variables=["eta1","eta2","pt1","pt2","n_jets","n_vtx","phi1","phi2","m_ll","ht","n_photons","pt_g1","met","nElectrons","nMuons","deltaR1_g1","deltaR2_g1","deltaRLL"]
    #variables=["eta1","eta2","pt1","pt2","met","pt_g1","m_ll","phi1","phi2","n_jets","ht"]
    variables=["pt_g1","pt1","eta1","eta2","pt2","met","m_ll","phi1","phi2","n_jets","ht"]
    #variables=["eta1","eta2","pt1","pt2","pt_g1","m_ll","phi1","phi2","n_jets","ht"]
    #variables=["eta1","eta2","pt1","pt2","pt_g1","phi1","phi2","n_jets","ht"]
    #variables=["pt1"]
    bkgs=[DYjetsNLO,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma,singletop,zz4l,wz,ww]
    groups=["CRTT"]
    binnings_=CR_tt.binnings.copy()
    
    #scaleFactors=frange(0.6,1.5,0.01)
    scaleFactors=frange(0.01,1.9,0.01)
    #scaleFactorsTTG=frange(0.6,1.5,0.01)
    
    
    additionalFolder ="nom/"
    
    
    saveValues={}
    for group in groups:
        saveValues[group]={}
        for variable in variables:
            saveValues[group][variable]={}
            chiValues=[]
            SFtt=[]
            SFttg=[]
            for SF in scaleFactors:
                tempEM=(drawSameHistogram("EM",group+"/EM/"+additionalFolder+variable, bkgs, additional=[dataLL],binning=binningToUse_[variable],xTitle=labels[variable][0],toScale=[tt,ttgamma],SF=SF))
                #chiValues.append(tempEM[0])
                SFtt.append(tempEM[0])
                SFttg.append(tempEM[1])
                chiValues.append(tempEM[2][0])
            gr = TGraph();
            gr2 = TGraph();

            #print SFtt
            #print SFttg
            #print chiValues

            for i in range(len(chiValues)):
                #gr.SetPoint(i,scaleFactors[i],chiValues[i])
                gr.SetPoint(i,SFtt[i],chiValues[i])
                gr2.SetPoint(i,SFttg[i],chiValues[i])
            
            c = TCanvas("canvas","",800,800)

            chiValues=np.array(chiValues)
            scaleFactorsTT=np.array(SFtt)
            scaleFactorsTTG=np.array(SFttg)

            findMinIndex=chiValues.argmin()
            plusOneLeft=chiValues[findMinIndex]+1
            plusOneRight=chiValues[findMinIndex]+1
            print chiValues[findMinIndex],scaleFactorsTT[findMinIndex]
            print chiValues[findMinIndex],scaleFactorsTTG[findMinIndex]
            
            up=0.
            down=0.
            
            gr.Fit("pol4")
            gr.GetFunction("pol4").SetLineColor(kRed)
            gr.GetFunction("pol4").SetLineWidth(1)
            
            fit=gr.GetFunction("pol4")
            
            gr2.Fit("pol4")
            gr2.GetFunction("pol4").SetLineColor(kRed)
            gr2.GetFunction("pol4").SetLineWidth(1)
            
            fit2=gr2.GetFunction("pol4")
            
            chiResult=fit.GetMinimum()
            sfResult=fit.GetMinimumX()
            sfResultUp=fit.GetX(chiResult+1,0.,sfResult)
            sfResultDn=fit.GetX(chiResult+1,sfResult,100.)
            
            chiResult2=fit2.GetMinimum()
            sfResult2=fit2.GetMinimumX()
            sfResult2Up=fit2.GetX(chiResult2+1,0.,sfResult2)
            sfResult2Dn=fit2.GetX(chiResult2+1,sfResult2,100.)
            
            up=abs(sfResult-sfResultUp)
            down=abs(sfResult-sfResultDn)
            up2=abs(sfResult2-sfResult2Up)
            down2=abs(sfResult2-sfResult2Dn)
            
            gr.SetTitle("; Scale Factor #alpha; #chi^{2}")
            gr.SetMarkerStyle(20)
            gr.SetMarkerSize(0.5)
            gr.SetName("grTT")
            
            gr2.SetTitle("; Scale Factor #beta; #chi^{2}")
            gr2.SetMarkerStyle(20)
            gr2.SetMarkerSize(0.5)
            gr2.SetName("grTTG")


            gr.Draw("ACP")
            l = ROOT.TLatex( 0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
            l.SetNDC()
            l.Draw()
            l2 = ROOT.TLatex( 0.47, .65, "#scale[0.66]{#font[52]{#alpha = %.3f^{#plus%.3f}_{#minus%.3f}}}"%(sfResult,up,down))
            l2.SetNDC()
            l2.Draw()
            l3 = ROOT.TLatex( 0.47, .6, "#scale[0.66]{#font[52]{#chi}^{2} = %.2f / %.0f}"%(chiResult,tempEM[1]))
            l3.SetNDC()
            l3.Draw()
            l4 = ROOT.TLatex( 0.47, .7, "#scale[0.66]{#font[52]{t#bar{t}(+#gamma) Control Region}}")
            l4.SetNDC()
            l4.Draw()
            lum = ROOT.TLatex( .62, .95, "%.1f fb^{-1} (%s TeV)"%(aux.intLumi/1000., aux.Label.cmsEnergy) )
            lum.SetNDC()
            lum.Draw()
            
            leg=TLegend(0.5,0.8,0.7,0.9)
            leg.AddEntry("grTT","measured points","p")
            leg.AddEntry(gr.GetFunction("pol4"),"polynomial fit")
            leg.Draw()
            
            c.Update()
            c.SaveAs('tt_free/chiOverlap/TT'+'_'+variable+addName+'.pdf')
            #c.SaveAs('tt_free/chi/TT'+'_'+variable+addName+'.root')
            
            
            c = TCanvas("canvas","",800,800)
            gr2.Draw("ACP")
            l = ROOT.TLatex( 0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
            l.SetNDC()
            l.Draw()
            l2 = ROOT.TLatex( 0.47, .65, "#scale[0.66]{#font[52]{#beta = %.3f^{#plus%.3f}_{#minus%.3f}}}"%(sfResult2,up2,down2))
            l2.SetNDC()
            l2.Draw()
            l3 = ROOT.TLatex( 0.47, .6, "#scale[0.66]{#font[52]{#chi}^{2} = %.2f / %.0f}"%(chiResult2,tempEM[2][1]))
            l3.SetNDC()
            l3.Draw()
            l4 = ROOT.TLatex( 0.47, .7, "#scale[0.66]{#font[52]{t#bar{t}(+#gamma) Control Region}}")
            l4.SetNDC()
            l4.Draw()
            lum = ROOT.TLatex( .62, .95, "%.1f fb^{-1} (%s TeV)"%(aux.intLumi/1000., aux.Label.cmsEnergy) )
            lum.SetNDC()
            lum.Draw()
            
            leg=TLegend(0.5,0.8,0.7,0.9)
            leg.AddEntry("grTTG","measured points","p")
            leg.AddEntry(gr2.GetFunction("pol4"),"polynomial fit")
            leg.Draw()
            
            c.Update()
            c.SaveAs('tt_free/chiOverlap/TTG'+'_'+variable+addName+'.pdf')
            
            bestValue=scaleFactorsTT[findMinIndex]
            erUp=up
            erDown=down
            bestValue2=scaleFactorsTTG[findMinIndex]
            erUp2=up2
            erDown2=down2
            saveValues[group][variable]["valueTT"]=bestValue
            saveValues[group][variable]["erUpTT"]=erUp
            saveValues[group][variable]["erDownTT"]=erDown
            saveValues[group][variable]["valueTTG"]=bestValue2
            saveValues[group][variable]["erUpTTG"]=erUp2
            saveValues[group][variable]["erDownTTG"]=erDown2
    #pkl.dump( saveValues, open( "tt_free/TT_chi"+addName+".pkl", "wb" ) )
    pkl.dump( saveValues, open( "tt_free/TTOverlap_chi"+addName+".pkl", "wb" ) )
    

def main():
    import style
    style.defaultStyle()
        
    nominalBinningTT=CR_tt.binnings.copy()
    nominalBinningTT["eta1"] = frange(0., 2.4, 0.1)
    nominalBinningTT["eta2"] = frange(0., 2.4, 0.1)
    nominalBinningTT["pt1"]  = frange(20.,150,10.)
    nominalBinningTT["pt2"]  = frange(20.,150,10.)
    #nominalBinningTT["pt_g1"]  = frange(20.,200,20.)
    nominalBinningTT["pt_g1"]  = frange(20.,200,10.)
    nominalBinningTT["phi1"]  = frange(0., 3.1, 0.1)
    nominalBinningTT["phi2"]  = frange(0., 3.1, 0.1)
    nominalBinningTT["n_jets"]  = frange(0., 7, 7)
    nominalBinningTT["ht"]  = frange(0., 1000.,50)
    nominalBinningTT["met"]  = frange(0., 450.,50)
    #nominalBinningTT["met"]  = [0,25,50,75,100,150,250,450]
    nominalBinningTT["m_ll"]  = frange(40.,300.,20)
    
    drawTT(binningToUse_=nominalBinningTT,addName="")


if __name__=="__main__":
    main()
