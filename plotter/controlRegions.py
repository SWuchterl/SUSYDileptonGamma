from dataMC import labels,frange

import ROOT
from ROOT import *
from array import array
from include import *
import numpy as np
import pickle as pkl

binnings = {
    'pt1':              frange(0,100,10)+frange(100,200,25)+range(200,350,50),
    'pt2':              frange(0,100,10)+frange(100,200,25),
    'pt3':              range(0, 200, 10),
    'pt4':              range(0, 200, 10),
    'eta1':             frange(0., 2.6, 0.1),
    'eta2':             frange(0., 2.6, 0.1),
    'eta3':             frange(0., 2.6, 0.1),
    'eta4':             frange(0., 2.6, 0.1),
    'phi1':             frange(0., 3.50, 0.1),
    'phi2':             frange(0., 3.50, 0.1),
    'phi3':             frange(0., 3.50, 0.1),
    'phi4':             frange(0., 3.50, 0.1),
    'ht':               frange(0., 1000.,50),
    'met':               [0,25,50,75,100,150,190,230,500],
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



def divideDatasetIntegrals( numerator, denominator, name ):
    numMerged = sum(numerator)
    h_num = numMerged.getHist( name )
    num = h_num.Integral(0,-1)
    denMerged = sum(denominator)
    h_den = denMerged.getHist( name )
    den = h_den.Integral(0,-1)
    return num/den if den else 1.
    
def calculateSFAndError( numerator_data, denominator_toScale,additional_fix, name ):
    folder= (name.split("/"))[0]
    
    numMerged_data = sum(numerator_data)
    h_num_data = numMerged_data.getHist( name )
    num_dataErr=(ROOT.Double(0))
    num_data = h_num_data.IntegralAndError(0,-1,num_dataErr)
    
    #denMerged_toScale = sum(denominator_toScale)
    #h_den_toScale = denMerged_toScale.getHist( name )
    #den_toScaleErr=(ROOT.Double(0))
    #den_toScale = h_den_toScale.IntegralAndError(0,-1,den_toScaleErr)
    
    h_den_toScale = denominator_toScale[0].getHist(name).Clone()
    avgTopPtWeightHisto = denominator_toScale[0].getHist(folder+"/weight_topPt")
    avgNIsrWeightHisto = denominator_toScale[0].getHist(folder+"/weight_nISR")
    avgEWKinoWeightHisto = denominator_toScale[0].getHist(folder+"/weight_EWKinoPairPt")
    avgleptonWeightHisto = denominator_toScale[0].getHist(folder+"/weight_leptonPairPt")
        
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
    
    h_den_toScale.Scale(1./avgTopPtWeight)
    h_den_toScale.Scale(1./avgNIsrWeight)
    h_den_toScale.Scale(1./avgEWKinoWeight)
    h_den_toScale.Scale(1./avgleptonWeight)
    for d in denominator_toScale[1:]:
        temp_h = d.getHist(name).Clone()
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
        
        temp_h.Scale(1./avgTopPtWeight)
        temp_h.Scale(1./avgNIsrWeight)
        temp_h.Scale(1./avgEWKinoWeight)
        temp_h.Scale(1./avgleptonWeight)
        
        h_den_toScale.Add(temp_h)

    den_toScaleErr=(ROOT.Double(0))
    den_toScale = h_den_toScale.IntegralAndError(0,-1,den_toScaleErr)
    
    #addMerged_fix = sum(additional_fix)
    #h_add_fix = addMerged_fix.getHist( name )
    #add_fixError=ROOT.Double(0)
    #add_fix = h_add_fix.IntegralAndError(0,-1,add_fixError)
    h_add_fix = additional_fix[0].getHist(name).Clone()
    avgTopPtWeightHisto = additional_fix[0].getHist(folder+"/weight_topPt")
    avgNIsrWeightHisto = additional_fix[0].getHist(folder+"/weight_nISR")
    avgEWKinoWeightHisto = additional_fix[0].getHist(folder+"/weight_EWKinoPairPt")
    avgleptonWeightHisto = additional_fix[0].getHist(folder+"/weight_leptonPairPt")
        
    if avgTopPtWeightHisto.Integral()>0:
        avgTopPtWeight= avgTopPtWeightHisto.GetMean()
    else:
        avgTopPtWeight=1.
    if avgNIsrWeightHisto.Integral()>0 :
        avgNIsrWeight= avgNIsrWeightHisto.GetMean()
    else:
        avgNIsrWeight=1.
    if avgEWKinoWeightHisto.Integral()>0:
        avgEWKinoWeight = avgEWKinoWeightHisto.GetMean()
    else:
        avgEWKinoWeight=1.
    if avgleptonWeightHisto.Integral()>0: 
        avgleptonWeight = avgleptonWeightHisto.GetMean()
    else:
        avgleptonWeight=1.
    
    h_add_fix.Scale(1./avgTopPtWeight)
    h_add_fix.Scale(1./avgNIsrWeight)
    h_add_fix.Scale(1./avgEWKinoWeight)
    h_add_fix.Scale(1./avgleptonWeight)
    for d in additional_fix[1:]:
        temp_h = d.getHist(name).Clone()
        avgTopPtWeightHisto = d.getHist(folder+"/weight_topPt")
        avgNIsrWeightHisto = d.getHist(folder+"/weight_nISR")
        avgEWKinoWeightHisto = d.getHist(folder+"/weight_EWKinoPairPt")
        avgleptonWeightHisto = d.getHist(folder+"/weight_leptonPairPt")
    
        #print d
        #print avgTopPtWeightHisto.Integral(),folder
        
        if avgTopPtWeightHisto.Integral()>0:
            avgTopPtWeight= avgTopPtWeightHisto.GetMean()
        else:
            avgTopPtWeight=1.
        if avgNIsrWeightHisto.Integral()>0 :
            avgNIsrWeight= avgNIsrWeightHisto.GetMean()
        else:
            avgNIsrWeight=1.
        if avgEWKinoWeightHisto.Integral()>0:
            avgEWKinoWeight = avgEWKinoWeightHisto.GetMean()
        else:
            avgEWKinoWeight=1.
        if avgleptonWeightHisto.Integral()>0: 
            avgleptonWeight = avgleptonWeightHisto.GetMean()
        else:
            avgleptonWeight=1.
        
        temp_h.Scale(1./avgTopPtWeight)
        temp_h.Scale(1./avgNIsrWeight)
        temp_h.Scale(1./avgEWKinoWeight)
        temp_h.Scale(1./avgleptonWeight)
        
        h_add_fix.Add(temp_h)

    add_fixError=(ROOT.Double(0))
    add_fix = h_add_fix.IntegralAndError(0,-1,den_toScaleErr)
    
    
    
        
    alpha = (num_data-add_fix)/(den_toScale)
    
    alphaErr =np.sqrt( (num_dataErr/den_toScale)**2. +(add_fixError/den_toScale)**2. +(den_toScaleErr*(num_data-add_fix)/(den_toScale)**2.)**2. )
    
    #return num/den if den else 1.
    #return [alpha,alphaErr] if den_toScale else [1.,0.]
    return [alpha,alphaErr/alpha] if den_toScale else [1.,0.]


def drawSameHistogram( sampleNames, name, bkg=[], additional=[], binning=None, binningName="", scaleToData=True, xTitle=None, yTitle=None,toScale=[] ):
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()
    
    #style.divideByBinWidth = True
    style.divideByBinWidth = False
    
    #style.minimumOne=False
    style.minimumOne=True
    
    #yTitle=None
    
    scale = 1.
    scaleErr=0.
    addBKG=bkg[:]
    for y in toScale:
        addBKG.remove(y)
    if scaleToData: scale,scaleErr = calculateSFAndError( [ i for i in additional if "Data" in i.label ], toScale, addBKG, name )




    folder= (name.split("/"))[0]

    histsForSystError=[]

    #sumToScale=sum(toScale)
    #print "sum",sumToScale.names


    for d in bkg[-1::-1]:
        h = d.getHist( name )
        
        #print "name",d.names
        
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
        if d in toScale:
            h.Scale(scale)
        
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

        if d in toScale:
            mcSystUncert = scaleErr
        else:
            mcSystUncert = 0.

        
        #if d=sumToScale:
        
        SystHist = aux.getSysHisto(h, mcSystUncert)
        histsForSystError.append(SystHist)
        
        
    forSystH=toScale[0].getHist(name)
    avgTopPtWeightHisto = toScale[0].getHist(folder+"/weight_topPt")
    avgNIsrWeightHisto = toScale[0].getHist(folder+"/weight_nISR")
    avgEWKinoWeightHisto = toScale[0].getHist(folder+"/weight_EWKinoPairPt")
    avgleptonWeightHisto = toScale[0].getHist(folder+"/weight_leptonPairPt")
    
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
    
    forSystH.Scale(1./avgTopPtWeight)
    forSystH.Scale(1./avgNIsrWeight)
    forSystH.Scale(1./avgEWKinoWeight)
    forSystH.Scale(1./avgleptonWeight)
    if binning: 
        forSystH = aux.rebin( forSystH, binning )
    aux.appendFlowBin( forSystH )
            
    for d in toScale[1:]:
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
        
        h.Scale(1./avgTopPtWeight)
        h.Scale(1./avgNIsrWeight)
        h.Scale(1./avgEWKinoWeight)
        h.Scale(1./avgleptonWeight)
        
        if binning: 
            h = aux.rebin( h, binning )
        aux.appendFlowBin( h )
        
        forSystH.Add(h)
    
    #totalSysHist=aux.getSysHisto(forSystH,0.4)
    totalSysHist=aux.addHists(aux.getSysHisto(forSystH,0.4))
    

    #totalSysHist=aux.addHists(histsForSystError[0])
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
            r = ratio.Ratio( "Data/MC", dataHist, hsm,sysHisto=totUnc )
            r.draw(0.5,1.5)

        info = ""
        #info = "ee" if "EE" in sampleNames else "#mu#mu" if "MM" in sampleNames else ""
        l = aux.Label(info="#scale[0.7]{%s}"%info, sim=not((dataDoubleMuon in additional)or(dataDoubleEG in additional)or(dataHt in additional)or(dataDoubleSF in additional)or(dataMuonEG in additional)))
        #l = aux.Label(info="#scale[0.7]{%s}"%info, sim=False)

        if binningName: binningName = "_"+binningName
        name = name.replace("/","__")
        saveName = "sameHistograms_{}_{}{}".format(sampleNames, name, binningName )
        aux.save("DataMC_"+saveName,folder="plots_CR/" )

    return [scale,scaleErr]




def drawDY():
    bkgs=[DYjetsNLO,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma,singletop,zz4l,wz,ww]
    variables=["eta1","eta2","pt1","pt2","n_jets","n_vtx","phi1","phi2","m_ll","ht","n_photons","pt_g1","met","nElectrons","nMuons","deltaR1_g1","deltaR2_g1","deltaRLL"]
    groups=["ControlRegionDY"]
    binnings_=binnings.copy()
    binnings_["met"]=frange(0,105,10)
    toSave={}
    toSave["EE"]={}
    toSave["LL"]={}
    toSave["MM"]={}
    for group in groups:
        for variable in variables:
            #drawSameHistogram("EE",group+"EE/"+variable, bkgs, additional=[dataDoubleEG],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("MM",group+"MM/"+variable, bkgs, additional=[dataDoubleMuon],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("EM",group+"EM/"+variable, bkgs, additional=[dataMuonEG],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("LL",group+"/"+variable, bkgs, additional=[dataDoubleSF],binning=binnings[variable],xTitle=labels[variable][0])
            toSave["EE"][variable] = drawSameHistogram("EE",group+"EE/"+variable, bkgs, additional=[dataDoubleEG],binning=binnings_[variable],xTitle=labels[variable][0],toScale=[DYjetsNLO,zgamma])
            toSave["MM"][variable] = drawSameHistogram("MM",group+"MM/"+variable, bkgs, additional=[dataDoubleMuon],binning=binnings_[variable],xTitle=labels[variable][0],toScale=[DYjetsNLO,zgamma])
            #drawSameHistogram("EM",group+"EM/"+variable, bkgs, additional=[dataMuonEG],binning=binnings[variable],xTitle=labels[variable][0],toScale=[DYjetsNLO,zgamma])
            toSave["LL"][variable] = drawSameHistogram("LL",group+"/"+variable, bkgs, additional=[dataDoubleSF],binning=binnings_[variable],xTitle=labels[variable][0],toScale=[DYjetsNLO,zgamma])
        pkl.dump( toSave, open( "plots_CR/factors/"+group+".pkl", "wb" ) )
            
def drawZZ():
    variables=["eta1","eta2","eta3","eta4","pt1","pt3","pt4","pt2","n_jets","n_vtx","phi1","phi2","phi3","phi4","m_ll","m_ll2","ht","n_photons","met","nElectrons","nMuons"]
    bkgs=[DYjetsNLO,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma,singletop,zz4l,wz,ww]
    groups=["ControlRegionZZ"]
    toSave={}
    toSave["LL"]={}
    toSave["EE"]={}
    toSave["MM"]={}
    binnings_=binnings.copy()
    binnings_["eta1"]=frange(0., 2.61, 0.4);
    binnings_["eta2"]=frange(0., 2.61, 0.4);
    binnings_["eta3"]=frange(0., 2.61, 0.4);
    binnings_["eta4"]=frange(0., 2.61, 0.4);
    binnings_["phi1"]=frange(0., 3.51, 0.25);
    binnings_["phi2"]=frange(0., 3.51, 0.25);
    binnings_["phi3"]=frange(0., 3.51, 0.25);
    binnings_["phi4"]=frange(0., 3.51, 0.25);
    binnings_["m_ll"]=frange(80., 110., 1.);
    binnings_["m_ll2"]=frange(60., 120., 5.);
    binnings_["pt1"]=frange(25,160,10);
    binnings_["pt2"]=frange(25,160,10);
    binnings_["pt3"]=frange(25,160,10);
    binnings_["pt4"]=frange(25,160,10);
    for group in groups:
        for variable in variables:
            #drawSameHistogram("EE",group+"EE/"+variable, bkgs, additional=[dataDoubleEG],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("MM",group+"MM/"+variable, bkgs, additional=[dataDoubleMuon],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("EM",group+"EM/"+variable, bkgs, additional=[dataMuonEG],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("LL",group+"/"+variable, bkgs, additional=[dataLL],binning=binnings_[variable],xTitle=labels[variable][0])
            toSave["LL"][variable] = drawSameHistogram("LL",group+"/"+variable, bkgs, additional=[dataLL],binning=binnings_[variable],xTitle=labels[variable][0],toScale=[zz4l])
            toSave["EE"][variable] = drawSameHistogram("EE",group+"EE/"+variable, bkgs, additional=[dataDoubleEG],binning=binnings_[variable],xTitle=labels[variable][0],toScale=[zz4l])
            toSave["MM"][variable] = drawSameHistogram("MM",group+"MM/"+variable, bkgs, additional=[dataDoubleMuon],binning=binnings_[variable],xTitle=labels[variable][0],toScale=[zz4l])
        pkl.dump( toSave, open( "plots_CR/factors/"+group+".pkl", "wb" ) )
def drawWZ():
    variables=["eta1","eta2","eta3","pt1","pt2","pt3","n_jets","n_vtx","phi1","phi2","phi3","m_ll","ht","n_photons","met","nElectrons","nMuons","mTL3Met","met"]
    bkgs=[DYjetsNLO,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma,singletop,zz4l,wz,ww]
    groups=["ControlRegionWZ"]
    toSave={}
    toSave["LL"]={}
    toSave["EE"]={}
    toSave["MM"]={}
    binnings_=binnings.copy()
    binnings_["eta1"]=frange(0., 2.6, 0.2);
    binnings_["eta2"]=frange(0., 2.6, 0.2);
    binnings_["eta3"]=frange(0., 2.6, 0.2);
    binnings_["eta4"]=frange(0., 2.6, 0.2);
    binnings_["phi1"]=frange(0., 3.50, 0.1);
    binnings_["phi2"]=frange(0., 3.50, 0.1);
    binnings_["phi3"]=frange(0., 3.50, 0.1);
    binnings_["phi4"]=frange(0., 3.50, 0.1);
    binnings_["m_ll"]=frange(80., 110., 1.);
    binnings_["m_ll2"]=frange(60., 120., 5.);
    binnings_["pt1"]=frange(0,160,5);
    for group in groups:
        for variable in variables:
            #drawSameHistogram("EE",group+"EE/"+variable, bkgs, additional=[dataDoubleEG],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("MM",group+"MM/"+variable, bkgs, additional=[dataDoubleMuon],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("EM",group+"EM/"+variable, bkgs, additional=[dataMuonEG],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("LL",group+"/"+variable, bkgs, additional=[dataLL],binning=binnings_[variable],xTitle=labels[variable][0])
            toSave["LL"][variable] = drawSameHistogram("LL",group+"/"+variable, bkgs, additional=[dataLL],binning=binnings_[variable],xTitle=labels[variable][0],toScale=[wz])
            toSave["EE"][variable] = drawSameHistogram("EE",group+"EE/"+variable, bkgs, additional=[dataLL],binning=binnings_[variable],xTitle=labels[variable][0],toScale=[wz])
            toSave["MM"][variable] = drawSameHistogram("MM",group+"MM/"+variable, bkgs, additional=[dataLL],binning=binnings_[variable],xTitle=labels[variable][0],toScale=[wz])
        pkl.dump( toSave, open( "plots_CR/factors/"+group+".pkl", "wb" ) )
            
            
def drawTT():
    #variables=["eta1","eta2","pt1","pt2","n_jets","n_vtx","phi1","phi2","m_ll","ht","n_photons","pt_g1","met","nElectrons","nMuons","deltaR1_g1","deltaR2_g1","deltaRLL"]
    variables=["pt_g1"]
    #bkgs=[DYjetsNLO,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma,singletop,zz4l,wz,zz,ww,tt+ttgamma]
    bkgs=[DYjetsNLO,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma,singletop,zz4l,wz,zz,ww]
    groups=["ControlRegionTT"]
    binnings_=binnings.copy()
    binnings_["pt_g1"]=frange(25,80,5)+frange(80,140,10)+frange(140,200,20)
    toSave={}
    toSave["EM"]={}
    for group in groups:
        for variable in variables:
            #drawSameHistogram("EE",group+"EE/"+variable, bkgs, additional=[dataDoubleEG],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("MM",group+"MM/"+variable, bkgs, additional=[dataDoubleMuon],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("EM",group+"EM/"+variable, bkgs, additional=[dataMuonEG],binning=binnings[variable],xTitle=labels[variable][0])
            toSave["EM"][variable] = drawSameHistogram("EM",group+"EM/"+variable, bkgs, additional=[dataMuonEG],binning=binnings_[variable],xTitle=labels[variable][0],toScale=[tt,ttgamma])
            #drawSameHistogram("LL",group+"/"+variable, bkgs, additional=[dataDoubleSF],binning=binnings[variable],xTitle=labels[variable][0])
        pkl.dump( toSave, open( "plots_CR/factors/"+group+".pkl", "wb" ) )


def main():
    #drawDY()
    #drawZZ()
    drawTT()
    #drawWZ()


if __name__=="__main__":
    main()
