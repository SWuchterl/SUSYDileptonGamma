from dataMC import labels,frange

import ROOT
from ROOT import *
from array import array
from include import *
import numpy as np
import pickle as pkl
import CR_tt,CR_DY,CR_WZ,CR_ZZ
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


def drawSameHistogram( sampleNames, name, bkg=[], additional=[], binning=None, binningName="", scaleToData=True, xTitle=None, yTitle=None,toScale=[],SF=1. ):
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()
    
    style.divideByBinWidth = False
    style.minimumOne=True
        
    scale = SF
    scaleErr=0.
    
    
    addBKG=bkg[:]
    for y in toScale:
        addBKG.remove(y)


    folder= (name.split("/"))[0]

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
        
        if d in toScale:
            h.Scale(scale)
        
        h.Scale(1./avgTopPtWeight)
        h.Scale(1./avgNIsrWeight)
        h.Scale(1./avgEWKinoWeight)
        h.Scale(1./avgleptonWeight)
        
        
        if binning: 
            h = aux.rebin( h, binning )
        aux.appendFlowBin( h )
        
        m.addStack( h, d.label )

    dataHist = None
    for d in additional:
        h = d.getHist( name )
        if not h: continue
        if not h.Integral(): continue
        if (binning): 
            h = aux.rebin( h, binning )
        aux.appendFlowBin( h )

        if h.GetLineColor() == ROOT.kBlack: # data
            if not(binning): h.Sumw2(False)
            h.SetBinErrorOption( ROOT.TH1.kPoisson )
            dataHist = h
        else:
            h.drawOption_ = "hist e"
            h.SetLineWidth(3)

        m.add( h, d.label )

        m.sortStackByIntegral()
        
    totalMC=m.histsToStack
    totalData=m.hists[0]
    return(calculateChiSquared(totalData,totalMC))
        





def drawDY():
    bkgs=[DYjetsNLO,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma,singletop,zz4l,wz,ww]
    #variables=["eta1","eta2","pt1","pt2","n_jets","n_vtx","phi1","phi2","m_ll","ht","n_photons","pt_g1","met","nElectrons","nMuons","deltaR1_g1","deltaR2_g1","deltaRLL"]
    #variables=["eta1","eta2","pt1","pt2","met","pt_g1","phi1","phi2","n_jets","ht"]
    variables=["pt1"]
    #groups=["ControlRegionDY"]
    groups=["CRDY"]
    #binnings_=binnings.copy()
    binnings_=CR_DY.binnings.copy()
    #binnings_["met"]=frange(0,105,10)
    scaleFactors=frange(1.,1.151,0.0005)
    #scaleFactors=frange(1.,1.151,0.01)
    #chiValuesEE=[]
    #chiValuesLL=[]
    #chiValuesMM=[]
    saveValuesEE={}
    saveValuesMM={}
    saveValuesLL={}
    for group in groups:
        saveValuesEE[group]={}
        saveValuesMM[group]={}
        saveValuesLL[group]={}
        for variable in variables:
            saveValuesEE[group][variable]={}
            saveValuesMM[group][variable]={}
            saveValuesLL[group][variable]={}
            chiValuesEE=[]
            chiValuesMM=[]
            chiValuesLL=[]
            for SF in scaleFactors:
                tempEE=(drawSameHistogram("EE",group+"EE/"+variable, bkgs, additional=[dataDoubleEG],binning=binnings_[variable],xTitle=labels[variable][0],toScale=[DYjetsNLO,zgamma],SF=SF))
                tempMM=(drawSameHistogram("MM",group+"MM/"+variable, bkgs, additional=[dataDoubleMuon],binning=binnings_[variable],xTitle=labels[variable][0],toScale=[DYjetsNLO,zgamma],SF=SF))
                tempLL=(drawSameHistogram("LL",group+"/"+variable, bkgs, additional=[dataDoubleSF],binning=binnings_[variable],xTitle=labels[variable][0],toScale=[DYjetsNLO,zgamma],SF=SF))
                chiValuesEE.append(tempEE[0])
                chiValuesMM.append(tempMM[0])
                chiValuesLL.append(tempLL[0])
            grEE = TGraph();
            grMM = TGraph();
            grLL = TGraph();

            for i in range(len(chiValuesEE)):
                grEE.SetPoint(i,scaleFactors[i],chiValuesEE[i])
            for i in range(len(chiValuesMM)):
                grMM.SetPoint(i,scaleFactors[i],chiValuesMM[i])
            for i in range(len(chiValuesLL)):
                grLL.SetPoint(i,scaleFactors[i],chiValuesLL[i])
            
            c = TCanvas("canvas","",800,800)

            chiValuesEE=np.array(chiValuesEE)
            chiValuesMM=np.array(chiValuesMM)
            chiValuesLL=np.array(chiValuesLL)
            scaleFactors=np.array(scaleFactors)

            findMinIndexEE=chiValuesEE.argmin()
            plusOneLeftEE=chiValuesEE[findMinIndexEE]+1
            plusOneRightEE=chiValuesEE[findMinIndexEE]+1
            
            findMinIndexMM=chiValuesMM.argmin()
            plusOneLeftMM=chiValuesMM[findMinIndexMM]+1
            plusOneRightMM=chiValuesMM[findMinIndexMM]+1
            
            findMinIndexLL=chiValuesLL.argmin()
            plusOneLeftLL=chiValuesLL[findMinIndexLL]+1
            plusOneRightLL=chiValuesLL[findMinIndexLL]+1
            
            upEE=0.
            downEE=0.
            upMM=0.
            downMM=0.
            upLL=0.
            downLL=0.

            #for i in range(len(chiValuesEE)):
                #grEE.SetPoint(i,scaleFactors[i],chiValuesEE[i]-chiValuesEE[findMinIndexEE])
            #for i in range(len(chiValuesMM)):
                #grMM.SetPoint(i,scaleFactors[i],chiValuesMM[i]-chiValuesMM[findMinIndexMM])
            #for i in range(len(chiValuesLL)):
                #grLL.SetPoint(i,scaleFactors[i],chiValuesLL[i]-chiValuesLL[findMinIndexLL])
            
            for x in scaleFactors[0:findMinIndexEE]:
                y=grEE.Eval(x)
                if(abs(plusOneLeftEE-y)<0.05):
                    print "EE-",abs(x-scaleFactors[findMinIndexEE]),y
                    downEE=abs(x-scaleFactors[findMinIndexEE])
            for x in scaleFactors[findMinIndexEE:-1]:
                y=grEE.Eval(x)
                if(abs(plusOneRightEE-y)<0.05):
                    print "EE+",abs(x-scaleFactors[findMinIndexEE]),y
                    upEE=abs(x-scaleFactors[findMinIndexEE])
                    
            for x in scaleFactors[0:findMinIndexMM]:
                y=grMM.Eval(x)
                if(abs(plusOneLeftMM-y)<0.05):
                    print "MM-",abs(x-scaleFactors[findMinIndexMM]),y
                    downMM=abs(x-scaleFactors[findMinIndexMM])
            for x in scaleFactors[findMinIndexMM:-1]:
                y=grMM.Eval(x)
                if(abs(plusOneRightMM-y)<0.05):
                    print "MM+",abs(x-scaleFactors[findMinIndexMM]),y
                    upMM=abs(x-scaleFactors[findMinIndexMM])
                    
            for x in scaleFactors[0:findMinIndexLL]:
                y=grLL.Eval(x)
                if(abs(plusOneLeftLL-y)<0.05):
                    print "LL-",abs(x-scaleFactors[findMinIndexLL]),y
                    downLL=abs(x-scaleFactors[findMinIndexLL])
            for x in scaleFactors[findMinIndexLL:-1]:
                y=grLL.Eval(x)
                if(abs(plusOneRightLL-y)<0.05):
                    print "LL+",abs(x-scaleFactors[findMinIndexLL]),y
                    upLL=abs(x-scaleFactors[findMinIndexLL])
            
            #grEE.SetTitle("; Scale Factor #alpha; #chi^{2}")
            grEE.SetTitle("; Scale Factor #alpha; #chi^{2}-min")
            grEE.SetMarkerStyle(20)
            grEE.SetMarkerSize(0.5)
            grMM.SetTitle("; Scale Factor #alpha; #chi^{2}-min")
            grMM.SetMarkerStyle(20)
            grMM.SetMarkerSize(0.5)
            grLL.SetTitle("; Scale Factor #alpha; #chi^{2}-min")
            grLL.SetMarkerStyle(20)
            grLL.SetMarkerSize(0.5)
            
            grEE.Draw("ACP")
            l = ROOT.TLatex( 0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
            l.SetNDC()
            l.Draw()
            l2 = ROOT.TLatex( 0.47, .65, "#scale[0.66]{#font[52]{#alpha = %.3f^{#plus%.3f}_{#minus%.3f}}}"%(scaleFactors[findMinIndexEE],upEE,downEE))
            l2.SetNDC()
            l2.Draw()
            l3 = ROOT.TLatex( 0.47, .6, "#scale[0.66]{#font[52]{#chi}^{2} = %.2f / %.0f}"%(chiValuesEE[findMinIndexEE],tempEE[1]))
            l3.SetNDC()
            l3.Draw()
            l4 = ROOT.TLatex( 0.47, .8, "#scale[0.66]{#font[52]{DY Control Region ee}}")
            l4.SetNDC()
            l4.Draw()
            lum = ROOT.TLatex( .62, .95, "%.1f fb^{-1} (%s TeV)"%(aux.intLumi/1000., aux.Label.cmsEnergy) )
            lum.SetNDC()
            lum.Draw()
            c.Update()
            c.SaveAs('plots_CR/chi/DY_EE'+'_'+variable+'.pdf')
            c.SaveAs('plots_CR/chi/DY_EE'+'_'+variable+'.root')
            

            grMM.Draw("ACP")
            l = ROOT.TLatex( 0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
            l.SetNDC()
            l.Draw()
            l2 = ROOT.TLatex( 0.47, .65, "#scale[0.66]{#font[52]{#alpha = %.3f^{#plus%.3f}_{#minus%.3f}}}"%(scaleFactors[findMinIndexMM],upMM,downMM))
            l2.SetNDC()
            l2.Draw()
            l3 = ROOT.TLatex( 0.47, .6, "#scale[0.66]{#font[52]{#chi}^{2} = %.2f / %.0f}"%(chiValuesMM[findMinIndexMM],tempMM[1]))
            l3.SetNDC()
            l3.Draw()
            l4 = ROOT.TLatex( 0.47, .8, "#scale[0.66]{#font[52]{DY Control Region #mu#mu}}")
            l4.SetNDC()
            l4.Draw()
            lum = ROOT.TLatex( .62, .95, "%.1f fb^{-1} (%s TeV)"%(aux.intLumi/1000., aux.Label.cmsEnergy) )
            lum.SetNDC()
            lum.Draw()
            c.Update()
            c.SaveAs('plots_CR/chi/DY_MM'+'_'+variable+'.pdf')
            c.SaveAs('plots_CR/chi/DY_MM'+'_'+variable+'.root')
            

            grLL.Draw("ACP")
            l = ROOT.TLatex( 0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
            l.SetNDC()
            l.Draw()
            l2 = ROOT.TLatex( 0.47, .65, "#scale[0.66]{#font[52]{#alpha = %.3f^{#plus%.3f}_{#minus%.3f}}}"%(scaleFactors[findMinIndexLL],upLL,downLL))
            l2.SetNDC()
            l2.Draw()
            l3 = ROOT.TLatex( 0.47, .6, "#scale[0.66]{#font[52]{#chi}^{2} = %.2f / %.0f}"%(chiValuesLL[findMinIndexLL],tempLL[1]))
            l3.SetNDC()
            l3.Draw()
            l4 = ROOT.TLatex( 0.47, .8, "#scale[0.66]{#font[52]{DY Control Region ee+#mu#mu}}")
            l4.SetNDC()
            l4.Draw()
            lum = ROOT.TLatex( .62, .95, "%.1f fb^{-1} (%s TeV)"%(aux.intLumi/1000., aux.Label.cmsEnergy) )
            lum.SetNDC()
            lum.Draw()
            c.Update()
            c.SaveAs('plots_CR/chi/DY_LL'+'_'+variable+'.pdf')
            c.SaveAs('plots_CR/chi/DY_LL'+'_'+variable+'.root')
            
            bestValueLL=scaleFactors[findMinIndexLL]
            bestValueEE=scaleFactors[findMinIndexEE]
            bestValueMM=scaleFactors[findMinIndexMM]
            erUpLL=upLL
            erUpEE=upEE
            erUpMM=upMM
            erDownLL=downLL
            erDownEE=downEE
            erDownMM=downMM
            saveValuesEE[group][variable]["value"]=bestValueEE
            saveValuesEE[group][variable]["erUp"]=erUpEE
            saveValuesEE[group][variable]["erDown"]=erDownEE
            saveValuesMM[group][variable]["value"]=bestValueMM
            saveValuesMM[group][variable]["erUp"]=erUpMM
            saveValuesMM[group][variable]["erDown"]=erDownMM
            saveValuesLL[group][variable]["value"]=bestValueLL
            saveValuesLL[group][variable]["erUp"]=erUpLL
            saveValuesLL[group][variable]["erDown"]=erDownLL
    pkl.dump( saveValuesEE, open( "plots_CR/chi/DY_chi_EE.pkl", "wb" ) )
    pkl.dump( saveValuesLL, open( "plots_CR/chi/DY_chi_LL.pkl", "wb" ) )
    pkl.dump( saveValuesMM, open( "plots_CR/chi/DY_chi_MM.pkl", "wb" ) )
            
def drawZZ():
    #variables=["eta1","eta2","eta3","eta4","pt1","pt3","pt4","pt2","n_jets","n_vtx","phi1","phi2","phi3","phi4","m_ll","m_ll2","ht","n_photons","met","nElectrons","nMuons"]
    #variables=["eta1"]
    #variables=["eta1","eta2","pt1","pt2","met","m_ll","phi1","phi2","n_jets","ht"]
    variables=["pt1"]
    bkgs=[DYjetsNLO,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma,singletop,zz4l,wz,ww]
    #groups=["ControlRegionZZ"]
    groups=["CRZZ"]
    #scaleFactors=frange(1.2,1.601,0.001)
    scaleFactors=frange(1.2,1.61,0.001)
    chiValues=[]
    #binnings_=binnings.copy()
    binnings_=CR_ZZ.binnings.copy()
    #binnings_["eta1"]=frange(0., 2.61, 0.4);
    #frange(0., 2.405, 0.1)
    #binnings_["eta1"]=frange(0., 2.405, 0.4);
    #binnings_["eta2"]=frange(0., 2.61, 0.4);
    #binnings_["eta3"]=frange(0., 2.61, 0.4);
    #binnings_["eta4"]=frange(0., 2.61, 0.4);
    #binnings_["phi1"]=frange(0., 3.51, 0.25);
    #binnings_["phi2"]=frange(0., 3.51, 0.25);
    #binnings_["phi3"]=frange(0., 3.51, 0.25);
    #binnings_["phi4"]=frange(0., 3.51, 0.25);
    #binnings_["m_ll"]=frange(80., 110., 1.);
    #binnings_["m_ll2"]=frange(60., 120., 5.);
    #binnings_["pt1"]=frange(25,160,10);
    #binnings_["pt2"]=frange(25,160,10);
    #binnings_["pt3"]=frange(25,160,10);
    #binnings_["pt4"]=frange(25,160,10);
    saveValues={}
    for group in groups:
        saveValues[group]={}
        for variable in variables:
            saveValues[group][variable]={}
            chiValues=[]
            for SF in scaleFactors:
                tempLL=(drawSameHistogram("LL",group+"/"+variable, bkgs, additional=[dataLL],binning=binnings_[variable],xTitle=labels[variable][0],toScale=[zz4l],SF=SF))
                chiValues.append(tempLL[0])
            gr = TGraph();

            for i in range(len(chiValues)):
                gr.SetPoint(i,scaleFactors[i],chiValues[i])
            
            c = TCanvas("canvas","",800,800)

            chiValues=np.array(chiValues)
            scaleFactors=np.array(scaleFactors)

            findMinIndex=chiValues.argmin()
            plusOneLeft=chiValues[findMinIndex]+1
            plusOneRight=chiValues[findMinIndex]+1
            print chiValues[findMinIndex],scaleFactors[findMinIndex]
            
            up=0.
            down=0.
            
            for x in scaleFactors[0:findMinIndex]:
                y=gr.Eval(x)
                if(abs(plusOneLeft-y)<0.01):
                    print "-",abs(x-scaleFactors[findMinIndex]),y
                    down=abs(x-scaleFactors[findMinIndex])
            for x in scaleFactors[findMinIndex:-1]:
                y=gr.Eval(x)
                if(abs(plusOneRight-y)<0.01):
                    print "+",abs(x-scaleFactors[findMinIndex]),y
                    up=abs(x-scaleFactors[findMinIndex])
            
            gr.SetTitle("; Scale Factor #alpha; #chi^{2}")
            gr.SetMarkerStyle(20)
            gr.SetMarkerSize(0.5)
            
            gr.Draw("ACP")
            l = ROOT.TLatex( 0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
            l.SetNDC()
            l.Draw()
            l2 = ROOT.TLatex( 0.47, .65, "#scale[0.66]{#font[52]{#alpha = %.3f^{#plus%.3f}_{#minus%.3f}}}"%(scaleFactors[findMinIndex],up,down))
            l2.SetNDC()
            l2.Draw()
            l3 = ROOT.TLatex( 0.47, .6, "#scale[0.66]{#font[52]{#chi}^{2} = %.2f / %.0f}"%(chiValues[findMinIndex],tempLL[1]))
            l3.SetNDC()
            l3.Draw()
            l4 = ROOT.TLatex( 0.47, .8, "#scale[0.66]{#font[52]{ZZ Control Region}}")
            l4.SetNDC()
            l4.Draw()
            lum = ROOT.TLatex( .62, .95, "%.1f fb^{-1} (%s TeV)"%(aux.intLumi/1000., aux.Label.cmsEnergy) )
            lum.SetNDC()
            lum.Draw()
            c.Update()
            c.SaveAs('plots_CR/chi/ZZ'+'_'+variable+'.pdf')
            c.SaveAs('plots_CR/chi/ZZ'+'_'+variable+'.root')
            
            bestValue=scaleFactors[findMinIndex]
            erUp=up
            erDown=down
            saveValues[group][variable]["value"]=bestValue
            saveValues[group][variable]["erUp"]=erUp
            saveValues[group][variable]["erDown"]=erDown
    pkl.dump( saveValues, open( "plots_CR/chi/ZZ_chi.pkl", "wb" ) )
            
            
def drawWZ():
    #variables=["eta1","eta2","eta3","pt1","pt2","pt3","n_jets","n_vtx","phi1","phi2","phi3","m_ll","ht","n_photons","met","nElectrons","nMuons","mTL3Met","met"]
    #variables=["eta1","eta2","pt1","pt2","met","m_ll","phi1","phi2","n_jets","ht"]
    variables=["pt1"]
    bkgs=[DYjetsNLO,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma,singletop,zz4l,wz,ww]
    #groups=["ControlRegionWZ"]
    groups=["CRWZ"]
    #scaleFactors=frange(1.,1.325,0.001)
    scaleFactors=frange(1.1,1.305,0.001)
    chiValues=[]
    #binnings_=binnings.copy()
    binnings_=CR_WZ.binnings.copy()
    #binnings_["eta1"]=frange(0., 2.6, 0.2);
    #binnings_["eta1"]=frange(0., 2.401, 0.4);
    #binnings_["eta2"]=frange(0., 2.6, 0.2);
    #binnings_["eta3"]=frange(0., 2.6, 0.2);
    #binnings_["eta4"]=frange(0., 2.6, 0.2);
    #binnings_["phi1"]=frange(0., 3.50, 0.1);
    #binnings_["phi2"]=frange(0., 3.50, 0.1);
    #binnings_["phi3"]=frange(0., 3.50, 0.1);
    #binnings_["phi4"]=frange(0., 3.50, 0.1);
    #binnings_["m_ll"]=frange(80., 110., 1.);
    #binnings_["m_ll2"]=frange(60., 120., 5.);
    #binnings_["pt1"]=frange(0,160,5);
    saveValues={}
    for group in groups:
        saveValues[group]={}
        for variable in variables:
            saveValues[group][variable]={}
            chiValues=[]
            for SF in scaleFactors:
                tempLL=(drawSameHistogram("LL",group+"/"+variable, bkgs, additional=[dataLL],binning=binnings_[variable],xTitle=labels[variable][0],toScale=[wz],SF=SF))
                chiValues.append(tempLL[0])
            gr = TGraph();

            for i in range(len(chiValues)):
                gr.SetPoint(i,scaleFactors[i],chiValues[i])
            
            c = TCanvas("canvas","",800,800)

            chiValues=np.array(chiValues)
            scaleFactors=np.array(scaleFactors)

            findMinIndex=chiValues.argmin()
            plusOneLeft=chiValues[findMinIndex]+1
            plusOneRight=chiValues[findMinIndex]+1
            print chiValues[findMinIndex],scaleFactors[findMinIndex]
            
            up=0.
            down=0.
            
            for x in scaleFactors[0:findMinIndex]:
                y=gr.Eval(x)
                if(abs(plusOneLeft-y)<0.025):
                    print "-",abs(x-scaleFactors[findMinIndex]),y
                    down=abs(x-scaleFactors[findMinIndex])
            for x in scaleFactors[findMinIndex:-1]:
                y=gr.Eval(x)
                if(abs(plusOneRight-y)<0.025):
                    print "+",abs(x-scaleFactors[findMinIndex]),y
                    up=abs(x-scaleFactors[findMinIndex])
            
            gr.SetTitle("; Scale Factor #alpha; #chi^{2}")
            gr.SetMarkerStyle(20)
            gr.SetMarkerSize(0.5)
            
            gr.Draw("ACP")
            l = ROOT.TLatex( 0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
            l.SetNDC()
            l.Draw()
            l2 = ROOT.TLatex( 0.47, .65, "#scale[0.66]{#font[52]{#alpha = %.3f^{#plus%.3f}_{#minus%.3f}}}"%(scaleFactors[findMinIndex],up,down))
            l2.SetNDC()
            l2.Draw()
            l3 = ROOT.TLatex( 0.47, .6, "#scale[0.66]{#font[52]{#chi}^{2} = %.2f / %.0f}"%(chiValues[findMinIndex],tempLL[1]))
            l3.SetNDC()
            l3.Draw()
            l4 = ROOT.TLatex( 0.47, .8, "#scale[0.66]{#font[52]{WZ Control Region}}")
            l4.SetNDC()
            l4.Draw()
            lum = ROOT.TLatex( .62, .95, "%.1f fb^{-1} (%s TeV)"%(aux.intLumi/1000., aux.Label.cmsEnergy) )
            lum.SetNDC()
            lum.Draw()
            c.Update()
            c.SaveAs('plots_CR/chi/WZ'+'_'+variable+'.pdf')
            c.SaveAs('plots_CR/chi/WZ'+'_'+variable+'.root')
            
            bestValue=scaleFactors[findMinIndex]
            erUp=up
            erDown=down
            saveValues[group][variable]["value"]=bestValue
            saveValues[group][variable]["erUp"]=erUp
            saveValues[group][variable]["erDown"]=erDown
    pkl.dump( saveValues, open( "plots_CR/chi/WZ_chi.pkl", "wb" ) )
            
            
def drawTT():
    #variables=["eta1","eta2","pt1","pt2","n_jets","n_vtx","phi1","phi2","m_ll","ht","n_photons","pt_g1","met","nElectrons","nMuons","deltaR1_g1","deltaR2_g1","deltaRLL"]
    #variables=["eta1","eta2","pt1","pt2","met","pt_g1","m_ll","phi1","phi2","n_jets","ht"]
    variables=["pt1"]
    #variables=["met"]
    bkgs=[DYjetsNLO,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma,singletop,zz4l,wz,ww]
    #groups=["ControlRegionTT"]
    groups=["CRTT"]
    #binnings_=binnings.copy()
    binnings_=CR_tt.binnings.copy()
    #scaleFactors=frange(0.6,1.505,0.001)
    scaleFactors=frange(0.8,1.155,0.001)
    #scaleFactors=frange(0.6,1.505,0.01)
    #chiValues=[]
    
    saveValues={}
    for group in groups:
        saveValues[group]={}
        for variable in variables:
            saveValues[group][variable]={}
            chiValues=[]
            for SF in scaleFactors:
                #tempEM=(drawSameHistogram("EM",group+"EM/"+variable, bkgs, additional=[dataMuonEG],binning=binnings_[variable],xTitle=labels[variable][0],toScale=[tt,ttgamma],SF=SF))
                tempEM=(drawSameHistogram("EM",group+"EM/"+variable, bkgs, additional=[dataLL],binning=binnings_[variable],xTitle=labels[variable][0],toScale=[tt,ttgamma],SF=SF))
                chiValues.append(tempEM[0])
            gr = TGraph();

            for i in range(len(chiValues)):
                gr.SetPoint(i,scaleFactors[i],chiValues[i])
            
            c = TCanvas("canvas","",800,800)

            chiValues=np.array(chiValues)
            scaleFactors=np.array(scaleFactors)

            findMinIndex=chiValues.argmin()
            plusOneLeft=chiValues[findMinIndex]+1
            plusOneRight=chiValues[findMinIndex]+1
            print chiValues[findMinIndex],scaleFactors[findMinIndex]
            
            up=0.
            down=0.
            
            for x in scaleFactors[0:findMinIndex]:
                y=gr.Eval(x)
                if(abs(plusOneLeft-y)<0.03):
                    print "-",abs(x-scaleFactors[findMinIndex]),y
                    down=abs(x-scaleFactors[findMinIndex])
            for x in scaleFactors[findMinIndex:-1]:
                y=gr.Eval(x)
                if(abs(plusOneRight-y)<0.03):
                    print "+",abs(x-scaleFactors[findMinIndex]),y
                    up=abs(x-scaleFactors[findMinIndex])
            
            gr.SetTitle("; Scale Factor #alpha; #chi^{2}")
            gr.SetMarkerStyle(20)
            gr.SetMarkerSize(0.5)


            gr.Draw("ACP")
            l = ROOT.TLatex( 0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
            l.SetNDC()
            l.Draw()
            l2 = ROOT.TLatex( 0.47, .65, "#scale[0.66]{#font[52]{#alpha = %.3f^{#plus%.3f}_{#minus%.3f}}}"%(scaleFactors[findMinIndex],up,down))
            l2.SetNDC()
            l2.Draw()
            l3 = ROOT.TLatex( 0.47, .6, "#scale[0.66]{#font[52]{#chi}^{2} = %.2f / %.0f}"%(chiValues[findMinIndex],tempEM[1]))
            l3.SetNDC()
            l3.Draw()
            l4 = ROOT.TLatex( 0.47, .8, "#scale[0.66]{#font[52]{t#bar{t}(+#gamma) Control Region}}")
            l4.SetNDC()
            l4.Draw()
            lum = ROOT.TLatex( .62, .95, "%.1f fb^{-1} (%s TeV)"%(aux.intLumi/1000., aux.Label.cmsEnergy) )
            lum.SetNDC()
            lum.Draw()
            c.Update()
            c.SaveAs('plots_CR/chi/TT'+'_'+variable+'.pdf')
            c.SaveAs('plots_CR/chi/TT'+'_'+variable+'.root')
            
            bestValue=scaleFactors[findMinIndex]
            erUp=up
            erDown=down
            saveValues[group][variable]["value"]=bestValue
            saveValues[group][variable]["erUp"]=erUp
            saveValues[group][variable]["erDown"]=erDown
    pkl.dump( saveValues, open( "plots_CR/chi/TT_chi.pkl", "wb" ) )
    

def main():
    import style
    style.defaultStyle()
    drawDY()
    drawZZ()
    drawTT()
    drawWZ()


if __name__=="__main__":
    main()
