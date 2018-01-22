import ROOT
from ROOT import *
from array import array
from include import *
import numpy as np

list_of_variables = ['pt1',
                     'pt2',
                     'eta1',
                     'eta2',
                     'met',
                     'phi1',
                     'phi2',
                     'ht',
                     'gen_ht',
                     'm_ll',
                     'm_ll_e',
                     'm_ll_m',
                     'n_jets',
                     'n_vtx',
                     'pt_g1',
                     'eta_g1',
                     'phi_g1',
                     'sigmaIetaIeta_g1',
                     'DeltaEtaLL',
                     'DeltaPhiLL',
                     'DeltaEtaLLG',
                     'DeltaPhiLLG',
                     'DeltaRLL',
                     'DeltaRLLG',
                     'st',
                     'stmet',
                     'zpt'
                     ]


binnings = {
    #'pt1':              np.arange(0., 500., 10.),
    #'pt1':              np.concatenate((np.arange(0,200,10),np.arange(200,1750,50)),axis=0),
    #'pt1':              np.concatenate((np.arange(0,200,10),np.arange(200,350,50)),axis=0),
    'pt1':              np.concatenate((np.arange(0,100,10),np.arange(100,200,25),np.arange(200,350,50)),axis=0),
    'pt2':              np.arange(0., 300, 10),
    'eta1':             np.arange(0., 2.6, 0.1),
    #'eta1':             np.arange(0., 2.6, 0.05),
    'eta2':             np.arange(0., 2.6, 0.1),
    #'phi1':             np.arange(0., 3.50, 0.10),
    'phi1':             np.arange(0., 3.50, 0.35),
    'phi2':             np.arange(0., 3.50, 0.10),
    'ht':               np.arange(0., 5000.,100.),
    #'met':              np.arange(0., 1000.,100.),
    #'met':              np.arange(0., 1000.,500.),
    #'met':              np.concatenate((np.arange(0,400,50),np.arange(400,3000,1300)),axis=0),
    'met':              np.arange(0., 600.,25.),
    #'met':              np.concatenate((np.arange(0,200,50),np.arange(200,300,100),np.arange(300,3000,1400)),axis=0),
    'gen_ht':           np.arange(0., 5000.,100.),
    'm_ll':             np.arange(50., 130.,2.),
    #'m_ll':             np.arange(0., 650.,2.),
    #'m_ll':             np.concatenate((np.arange(0,100,10),np.arange(100,200,10),np.arange(200,500,50)),axis=0),
    'm_{ll#gamma}':             np.concatenate((np.arange(0,100,10),np.arange(100,200,10),np.arange(200,500,50)),axis=0),
    'm_llg':             np.concatenate((np.arange(0,100,10),np.arange(100,200,10),np.arange(200,500,50)),axis=0),
    'pt_llg':             np.concatenate((np.arange(0,100,10),np.arange(100,200,10),np.arange(200,500,50)),axis=0),
    #'m_ll':             np.concatenate((np.arange(0,200,10),np.arange(200,950,50)),axis=0),
    #'m_ll':             np.arange(80,110,1),
    'm_ll_e':           np.arange(0., 650.,10.),
    'm_ll_m':           np.arange(0., 650.,10.),
    'n_jets':           np.arange(0.,10.,1.),
    'n_photons':           np.arange(0.,10.,1.),
    'n_vtx':            np.arange(0.,40.,1.),
    #'pt_g1':            np.array((range(0,100,10)+np.array([100., 150., 200.,300,400.]))),
    #'pt_g1':            np.concatenate((np.arange(0,100,10),np.arange(100,950,50)),axis=0),
    'pt_g1':            np.concatenate((np.arange(0,100,10),np.arange(100,350,50)),axis=0),
    'eta_g1':           np.arange(0., 2.60, 0.10),
    'phi_g1':           np.arange(0., 3.50, 0.10),
    'sigmaIetaIeta_g1': np.arange(0., 0.04,0.001),
    'sigmaIphiIphi_g1': np.arange(0., 0.2,0.01),
    'deltaR1_g1':       np.arange(0., 6.,0.25),
    'deltaR2_g1':       np.arange(0., 6.,0.25),
    'r9_g1':            np.arange(0., 1.5,0.05),
    'hOverE_g1':        np.arange(0., 0.1,0.001),
    'deltaEtaLL':       np.arange(0.,6.,0.1),
    'deltaPhiLL':       np.arange(0.,6.,0.1),
    'deltaEtaLLG':      np.arange(0.,6.,0.1),
    'deltaPhiLLG':      np.arange(0.,6.,0.1),
    'deltaRLL':         np.arange(0.,6.,0.12),
    'deltaRLLG':        np.arange(0.,6.,0.1),
    'st':               np.arange(0.,5000.,100.),
    #'stmet':           np.arange(0.,5000.,100.),
    'stmet':            np.concatenate((np.arange(0,1000,100),np.arange(1000,4000,1000)),axis=0),
    'zpt':              np.arange(0.,2000.,100.),
    'mtll':             np.arange(0.,3000.,100.),
    'mtllg':            np.arange(0.,3000.,100.),
    'mtgmet':            np.concatenate((np.arange(0,200,50),np.arange(200,500,100),np.arange(500,4000,500)),axis=0),
    'mtl1met':            np.concatenate((np.arange(0,200,50),np.arange(200,500,100),np.arange(500,4000,500)),axis=0),
    'mtl2met':            np.concatenate((np.arange(0,200,50),np.arange(200,500,100),np.arange(500,4000,500)),axis=0),
    'mtllmet':            np.concatenate((np.arange(0,200,50),np.arange(200,500,100),np.arange(500,4000,500)),axis=0),
    'mtllgmet':            np.concatenate((np.arange(0,200,50),np.arange(200,500,100),np.arange(500,4000,500)),axis=0),
    'mt2':            np.arange(0.,300.,5.),
    'mzg_exo':            np.concatenate((np.arange(0,100,50),np.arange(100,200,50),np.arange(200,500,50)),axis=0),
    'gammaMotherID':            np.arange(0.,200.,1.)
}
labels = {
    'pt1': ["p_{T}^{leading}[GeV]","Events / 30 GeV"],
    'pt2': ["p_{T}^{trailing}[GeV]"," Events / 30 GeV"],
    'eta1': ["|#eta_{leading}|","Events / 0.1"],
    'eta2': ["|#eta_{trailing}|","Events / 0.1"],
    'phi1': ["|#Phi_{leading}|","Events / 0.35"],
    'phi2': ["|#Phi_{trailing}|","Events / 0.35"],
    'ht': ["H_{T} [GeV]", "Events / 10 GeV"],
    'met': ["E_{T}^{miss} [GeV]", "Events / 10 GeV "],
    'gen_ht': ["H_{T}^{gen} [GeV]", "Events / 10 GeV"],
    #'m_ll': ["m_{ll} from all leptons [GeV]" ,"Events / 65 GeV"],
    'm_ll': ["m_{ll} [GeV]" ,"Events / Bin"],
    'n_jets': ["N Jets", "Events / 1"],
    'n_photons': ["N Photons", "Events / 1"],
    'n_vtx': ["N Vertices", "Events / 1"],
    'pt_g1': ["p_{T}^{#gamma}","Events / 10 GeV"],
    'eta_g1': ["|#Eta|_{#gamma}","Events / 0.1"],
    'phi_g1': ["|#phi|_{#gamma}","Events / 0.1"],
    'sigmaIetaIeta_g1': ["#sigma_{I#etaI#eta}","Events / 0.001"],
    'sigmaIphiIphi_g1': ["#sigma_{I#phiI#phi}","Events / 0.001"],
    'r9_g1': ["r9","Events / 0.001"],
    'hOverE_g1': ["H/E","Events / 0.001"],
    'deltaR1_g1': ["#DeltaR1","Events / 0.001"],
    'deltaR2_g1': ["#DeltaR2","Events / 0.001"],
    'deltaEtaLL':       ["#Delta#Eta_{ll}","Events / 0.1"],
    'deltaPhiLL':       ["#Delta#Phi_{ll}","Events / 0.1"],
    'deltaEtaLLG':      ["#Delta#Eta_{ll,#gamma}","Events / 0.1"],
    'deltaPhiLLG':      ["#Delta#Phi_{ll,#gamma}","Events / 0.1"],
    'deltaRLL':         ["#DeltaR_{ll}","Events / 0.1"],
    'deltaRLLG':        ["#DeltaR_{ll,#gamma}","Events / 0.1"],    
    'st':        ["#itS_T [GeV]","Events / 100 GeV"],    
    'stmet':        [" S_T+E_{T}^{miss} [GeV]","Events / 100 GeV"],    
    'zpt':        ["#Z_{#p_{T}} [GeV]","Events / 50. GeV"],    
    'mtll':        ["m_{T}^{ll} [GeV]","Events / 100. GeV"],    
    'mtllg':        ["m_{T}^{ll#gamma} [GeV]","Events / 100. GeV"],    
    'mtgmet':        ["m_{T}^{#gamma,met} [GeV]","Events / 100. GeV"],    
    'mtl1met':        ["m_{T}^{l1,met} [GeV]","Events / 100. GeV"],    
    'mtl2met':        ["m_{T}^{l2,met} [GeV]","Events / 100. GeV"],    
    'mtllmet':        ["m_{T}^{ll,met} [GeV]","Events / 100. GeV"],    
    'mt2':        ["m_{T}^{2} [GeV]","Events / 100. GeV"],    
    'mtllgmet':        ["m_{T}^{ll#gamma,met} [GeV]","Events / 100. GeV"]  ,  
    'm_{ll#gamma}':        ["m_{ll#gamma} [GeV]","Events / 100. GeV"],    
    'm_llg':        ["m_{ll#gamma} [GeV]","Events / 100. GeV"],    
    'pt_llg':        ["pt_llg",""],
    'mzg_exo':        ["mzg_exo",""],
    'gammaMotherID':        ["gammaMotherID",""]    
}



def drawSameHistogram( sampleNames, name, bkg=[], additional=[], binning=None, binningName="", scaleToData=False, xTitle=None, yTitle=None ):
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()
    
    style.divideByBinWidth = True
    
    
    yTitle=None
    
    scale = 1.
    if scaleToData: scale = divideDatasetIntegrals( [ i for i in additional if "Data" in i.label ], bkg, name )


    for d in bkg[-1::-1]:
        h = d.getHist( name )
        if not h: continue
        if not h.Integral(): continue
        h.Scale(scale)
        if (binning.any()): 
            h = aux.rebin( h, binning )

        aux.appendFlowBin( h )
        #h.SetYTitle( aux.getYAxisTitle( h ) )
        if yTitle:
            h.SetYTitle( yTitle )
        else:
            h.SetYTitle( aux.getYAxisTitle( h ) )
        if xTitle:
            h.SetXTitle( xTitle )
        
        m.addStack( h, d.label )

    dataHist = None
    for d in additional:
        h = d.getHist( name )
        if not h: continue
        if not h.Integral(): continue
        if (binning.any()): 
            h = aux.rebin( h, binning )
        aux.appendFlowBin( h )

        if h.GetLineColor() == ROOT.kBlack: # data
            h.drawOption_ = "ep"
            h.SetMarkerStyle(20)
            h.SetMarkerSize(0.5)
            # disable errors for data, so that ErrorOption is working
            # TODO: kPoisson also for rebinned and scaled histograms
            if not(binning.any()): h.Sumw2(False)
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
            r = ratio.Ratio( "Data/MC", dataHist, hsm )
            r.draw(0.5,1.5)

        #info = ""
        info = "ee" if "EE" in sampleNames else "#mu#mu" if "MM" in sampleNames else ""
        l = aux.Label(info="#scale[0.7]{%s}"%info, sim=((dataDoubleMuon not in additional)or(dataDoubleEG not in additional)or(dataHt not in additional)or(dataSF not in additional)))

        if binningName: binningName = "_"+binningName
        name = name.replace("/","__")
        saveName = "sameHistograms_{}_{}{}".format(sampleNames, name, binningName )
        #aux.save("DataMC_"+saveName )
        #aux.save("DataMC_"+saveName,folder="plots_noVeto/" )
        aux.save("DataMC_"+saveName,folder="plots_EGRegression/" )
        #aux.save("DataMC_"+saveName,folder="plots_LO/" )
        #aux.save("DataMC_"+saveName,folder="plots_mllWeights/" )
        #aux.save("DataMC_"+saveName,folder="plots_NLO/" )

        
def main():
    #bkgs=[DYjets,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma]
    bkgs=[DYjetsLO,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma]
    #bkgs=[DYjetsNLO,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma]
    #bkgs=[DYjetsNLO,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma]
    #bkgs=[DYjets,zgamma]
    #bkgs=[DYjets]
    variables=["eta1","pt1","n_jets","n_vtx","phi1","m_ll","ht","n_photons","pt_g1","m_llg","pt_llg","mzg_exo","gammaMotherID"]
    #variables=["eta1","pt1","n_jets","n_vtx","phi1","m_ll","ht","n_photons","pt_g1","sigmaIetaIeta_g1","sigmaIphiIphi_g1","r9_g1","hOverE_g1","deltaR1_g1","deltaR2_g1","eta_g1","phi_g1"]
    #groups=["sel","onZ","dilep","onZG","mllG110"]
    #groups=["sel","onZ","dilep","onZG","exo"]
    groups=["EGRegression"]
    for group in groups:
        for variable in variables:
            #drawSameHistogram("EE",group+"EE/"+variable, bkgs, additional=[dataDoubleEG],binning=binnings[variable],yTitle=labels[variable][1],xTitle=labels[variable][0])
            drawSameHistogram("EE",group+"EE/"+variable, bkgs, additional=[dataDoubleEG],binning=binnings[variable],xTitle=labels[variable][0])
            drawSameHistogram("MM",group+"MM/"+variable, bkgs, additional=[dataDoubleMuon],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("MM",group+"MM/"+variable, bkgs, additional=[dataDoubleMuon],binning=binnings[variable],yTitle=labels[variable][1],xTitle=labels[variable][0])
            #drawSameHistogram("EE",group+"EE/"+variable, bkgs, additional=[zgamma],binning=binnings[variable],yTitle=labels[variable][1],xTitle=labels[variable][0])
            #drawSameHistogram("MM",group+"MM/"+variable, bkgs, additional=[zgamma],binning=binnings[variable],yTitle=labels[variable][1],xTitle=labels[variable][0])
            #drawSameHistogram("EE",group+"EE/"+variable, bkgs, additional=[zgamma],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("MM",group+"MM/"+variable, bkgs, additional=[zgamma],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("EE",group+"EE/"+variable, bkgs, additional=[DYjets,zgamma],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("MM",group+"MM/"+variable, bkgs, additional=[DYjets,zgamma],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("EE",group+"EE/"+variable, bkgs ,binning=binnings[variable],yTitle=labels[variable][1],xTitle=labels[variable][0])
            #drawSameHistogram("MM",group+"MM/"+variable, bkgs,binning=binnings[variable],yTitle=labels[variable][1],xTitle=labels[variable][0])
def main2():
    bkgs=[DYjets,zgamma,wwgamma,wzgamma,ttgamma,zz,tt,wjets]
    variables=["met","stmet","st","mtll","mtllg","mtl1met","mtl2met","mtllmet","mtllgmet","mtgmet","m_ll","n_vtx","mt2"]
    groups=["onZ"]
    for group in groups:
        for variable in variables:
            drawSameHistogram("EE+signal",group+"EE/"+variable, bkgs, additional=[t5bbbbzg_1800_1700,t5bbbbzg_1800_400,t5bbbbzg_1800_600,tching_1200,tching_400],binning=binnings[variable],yTitle=labels[variable][1],xTitle=labels[variable][0])
            drawSameHistogram("MM+signal",group+"MM/"+variable, bkgs, additional=[t5bbbbzg_1800_1700,t5bbbbzg_1800_400,t5bbbbzg_1800_600,tching_1200,tching_400],binning=binnings[variable],yTitle=labels[variable][1],xTitle=labels[variable][0])

if __name__=="__main__":
    main()
    #main2()
