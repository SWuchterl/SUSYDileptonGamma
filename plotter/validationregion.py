#from dataMC import binnings,labels,frange
from dataMC import labels,frange

import ROOT
from ROOT import *
from array import array
from include import *
import numpy as np
import pickle as pkl

binnings = {
    #'pt1':              frange(20,100,10)+frange(100,155,25),
    #'pt1':              frange(20,100,10)+frange(100,155,10),
    'pt1':              frange(25,151,25),
    'pt2':              frange(0,151,25),
    #'pt2':              frange(0, 300, 10),
    'eta1':             frange(0., 2.7, 0.52),
    'eta2':             frange(0., 2.7, 0.52),
    'phi1':             frange(0., 3.51, 0.7),
    'phi2':             frange(0., 3.51, 0.7),
    'ht':               frange(0., 2000.,250),
    #'met':               [0,25,50,75,100,150,190,230,500],
    'met':               frange(100,155,10),
    'm_ll':             frange(80.,111.,10),
    'm_llg':             frange(0,100,10)+frange(100,200,10)+frange(200,500,50),
    'pt_llg':             frange(0,100,10)+frange(100,200,10)+frange(200,500,50),
    'n_jets':           frange(0.,7.,1),
    'n_photons':         frange(0.,4.,1),
    'n_vtx':            frange(0.,40.,1),
    #'pt_g1':            frange(20,100,10)+frange(100,210,50),
    'pt_g1':            frange(25,126,25),
    'eta_g1':           frange(0., 2.60, 0.1),
    'phi_g1':           frange(0., 3.50, 0.1),
    'sigmaIetaIeta_g1': frange(0., 0.04,0.01),
    'sigmaIphiIphi_g1': frange(0., 0.2,0.01),
    'deltaR1_g1':       frange(0., 6.,0.1),
    'deltaR2_g1':       frange(0., 6.,0.1),
    'r9_g1':            frange(0., 1.5,0.05),
    'hOverE_g1':        frange(0., 0.1,0.01),
    'deltaEtaLL':       frange(0.,6.,0.1),
    'deltaPhiLL':       frange(0.,6.,0.1),
    'deltaEtaLLG':      frange(0.,6.,0.1),
    'deltaPhiLLG':      frange(0.,6.,0.01),
    'deltaRLL':         frange(0.,6.,0.12),
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



def divideDatasetIntegrals( numerator, denominator, name ):
    numMerged = sum(numerator)
    h_num = numMerged.getHist( name )
    num = h_num.Integral(0,-1)
    denMerged = sum(denominator)
    h_den = denMerged.getHist( name )
    den = h_den.Integral(0,-1)
    return num/den if den else 1.


def drawSameHistogram( sampleNames, name, bkg=[], additional=[], binning=None, binningName="", scaleToData=False, xTitle=None, yTitle=None, sfZZ=[], sfDY=[],sfTT=[],sfWZ=[] ):
#def drawSameHistogram( sampleNames, name, bkg=[], additional=[], binning=None, binningName="", scaleToData=True, xTitle=None, yTitle=None ):
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()
    
    #style.divideByBinWidth = True
    style.divideByBinWidth = False
    
    #style.minimumOne=False
    style.minimumOne=True
    
    #yTitle=None
    
    
    
    scale = 1.
    if scaleToData: scale = divideDatasetIntegrals( [ i for i in additional if "Data" in i.label ], bkg, name )
    

    #print name

    folder= (name.split("/"))[0]
    #avgTopPtWeight = d.getHist(folder+"/weight_topPt").GetMean()
    #print avgTopPtWeight

    histsForSystError=[]

    for d in bkg[-1::-1]:
        h = d.getHist( name )
        
        avgTopPtWeightHisto = d.getHist(folder+"/weight_topPt")
        avgNIsrWeightHisto = d.getHist(folder+"/weight_nISR")
        avgEWKinoWeightHisto = d.getHist(folder+"/weight_EWKinoPairPt")
        avgleptonWeightHisto = d.getHist(folder+"/weight_leptonPairPt")
        
        if avgTopPtWeightHisto:
            avgTopPtWeight= avgTopPtWeightHisto.GetMean()
        else:
            avgTopPtWeight=1.
        if avgNIsrWeightHisto :
            avgNIsrWeight= avgNIsrWeightHisto.GetMean()
        else:
            avgNIsrWeight=1.
        if avgEWKinoWeightHisto:
            avgEWKinoWeight = avgEWKinoWeightHisto.GetMean()
        else:
            avgEWKinoWeight=1.
        if avgleptonWeightHisto: 
            avgleptonWeight = avgleptonWeightHisto.GetMean()
        else:
            avgleptonWeight=1.
        
        if not h: continue
        if not h.Integral(): continue
        h.Scale(scale)
        
        
        if d in [tt,ttgamma]:
            h.Scale(sfTT[0])
        if d in [DYjetsNLO,zgamma]:
            h.Scale(sfDY[0])
        if d in [zz,zz4l]:
            h.Scale(sfZZ[0])
        if d in [wz]:
            h.Scale(sfWZ[0])
        
        h.Scale(1./avgTopPtWeight)
        h.Scale(1./avgNIsrWeight)
        h.Scale(1./avgEWKinoWeight)
        h.Scale(1./avgleptonWeight)
        
        
        if binning: 
        #if (binning.any()): 
            h = aux.rebin( h, binning )

        aux.appendFlowBin( h )
        h.SetYTitle( aux.getYAxisTitle( h ) )
        if yTitle:
            h.SetYTitle( yTitle )
        else:
            h.SetYTitle( aux.getYAxisTitle( h ) )
        if xTitle:
            h.SetXTitle( xTitle )
        
        m.addStack( h, d.label )

        if d in [tt,ttgamma]:
            mcSystUncert = sfTT[1]
        elif d in [DYjetsNLO,zgamma]:
            mcSystUncert = sfDY[1]
        elif d in [zz,zz4l]:
            mcSystUncert = sfZZ[1]
        elif d in [wz]:
            mcSystUncert = sfWZ[1]
        else:
            mcSystUncert = 0.
        #print d.label,mcSystUncert
        SystHist = aux.getSysHisto(h, mcSystUncert)
        histsForSystError.append(SystHist)
        
        
    totalSysHist=aux.addHists(histsForSystError[0])
    for temp in histsForSystError[1:]:
        totalSysHist=aux.addHists(totalSysHist,temp)
    #totUnc = aux.addHistUncert(totStat, dirSyst)
    totUnc = aux.addHistUncert(totalSysHist)
    #aux.drawOpt(totUnc, "totUnc")
    aux.drawOpt(totUnc, "sysUnc")
    m.add(totUnc, "syst. uncertainty")

    dataHist = None
    for d in additional:
        h = d.getHist( name )
        if not h: continue
        if not h.Integral(): continue
        #if (binning.any()): 
        if (binning): 
            h = aux.rebin( h, binning )
        aux.appendFlowBin( h )

        if h.GetLineColor() == ROOT.kBlack: # data
            h.drawOption_ = "ep"
            h.SetMarkerStyle(20)
            #h.SetMarkerSize(0.5)
            h.SetMarkerSize(0.7)
            #h.SetLineWidth(0.7)
            # disable errors for data, so that ErrorOption is working
            # TODO: kPoisson also for rebinned and scaled histograms
            #if not(binning.any()): h.Sumw2(False)
            if not(binning): h.Sumw2(False)
            h.SetBinErrorOption( ROOT.TH1.kPoisson )
            dataHist = h
        else:
            h.drawOption_ = "hist e"
            h.SetLineWidth(3)

        m.add( h, d.label )

    m.sortStackByIntegral()
    
    
    
    
    if m.Draw():

        # ratio
        hsm = m.hists[0].GetStack().Last()
        if dataHist:
            #r = ratio.Ratio( "Data/MC", dataHist, hsm )
            r = ratio.Ratio( "Data/MC", dataHist, hsm,sysHisto=totUnc  )
            r.draw(0.5,2.)

        info = ""
        #info = "ee" if "EE" in sampleNames else "#mu#mu" if "MM" in sampleNames else ""
        l = aux.Label(info="#scale[0.7]{%s}"%info, sim=not((dataDoubleMuon in additional)or(dataDoubleEG in additional)or(dataHt in additional)or(dataDoubleSF in additional)or(dataMuonEG in additional)))
        #l = aux.Label(info="#scale[0.7]{%s}"%info, sim=False)

        if binningName: binningName = "_"+binningName
        name = name.replace("/","__")
        saveName = "sameHistograms_{}_{}{}".format(sampleNames, name, binningName )
        aux.save("DataMC_"+saveName,folder="plots_VR/" )





def main():
    bkgs=[DYjetsNLO,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma,singletop,zz4l,wz,ww]
    variables=["eta1","eta2","pt1","pt2","n_jets","n_vtx","phi1","phi2","m_ll","ht","n_photons","pt_g1","met"]
    groups=["ValidationRegion"]
    pklZZ = pkl.load( open( "plots_CR/factors/ControlRegionZZ.pkl", "rb" ) )
    pklDY = pkl.load( open( "plots_CR/factors/ControlRegionDY.pkl", "rb" ) )
    pklTT = pkl.load( open( "plots_CR/factors/ControlRegionTT.pkl", "rb" ) )
    pklWZ = pkl.load( open( "plots_CR/factors/ControlRegionWZ.pkl", "rb" ) )
    ZZsf=pklZZ["LL"]["m_ll"]
    DYsf=pklDY["LL"]["pt1"]
    TTsf=pklTT["EM"]["pt1"]
    WZsf=pklWZ["LL"]["pt1"]
    #print ZZsf,DYsf,TTsf
    for group in groups:
        for variable in variables:
            #drawSameHistogram("EE",group+"EE/"+variable, bkgs, additional=[dataDoubleEG],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("MM",group+"MM/"+variable, bkgs, additional=[dataDoubleMuon],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("EM",group+"EM/"+variable, bkgs, additional=[dataMuonEG],binning=binnings[variable],xTitle=labels[variable][0])
            drawSameHistogram("LL",group+"/"+variable, bkgs, additional=[dataDoubleSF],binning=binnings[variable],xTitle=labels[variable][0],sfZZ=ZZsf,sfDY=DYsf,sfTT=TTsf,sfWZ=WZsf)






if __name__=="__main__":
    main()
