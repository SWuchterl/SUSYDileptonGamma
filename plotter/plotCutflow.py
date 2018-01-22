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
    'cutflow_fine':              np.concatenate((np.arange(0,200,10),np.arange(200,350,50)),axis=0)
}
labels = {
    #'pt1': ["p_{T}^{leading}[GeV]","Events / 30 GeV"],
    'cutflow_fine': ["p_{T}^{leading}[GeV]","Events / 30 GeV"]
}



def drawSameHistogram( sampleNames, name, bkg=[], additional=[], binning=None, binningName="", scaleToData=False, xTitle=None, yTitle=None ):
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()
    
    #style.divideByBinWidth = True
    
    
    yTitle=None
    
    scale = 1.
    if scaleToData: scale = divideDatasetIntegrals( [ i for i in additional if "Data" in i.label ], bkg, name )


    for d in bkg[-1::-1]:
        h = d.getHist( name )
        if not h: continue
        if not h.Integral(): continue
        h.Scale(scale)
        #if (binning.any()): 
        if (binning): 
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
        #if (binning.any()): 
        if (binning): 
            h = aux.rebin( h, binning )
        aux.appendFlowBin( h )

        if h.GetLineColor() == ROOT.kBlack: # data
            h.drawOption_ = "ep"
            h.SetMarkerStyle(20)
            h.SetMarkerSize(0.5)
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
            r = ratio.Ratio( "Data/MC", dataHist, hsm )
            r.draw(0.5,1.5)

        #info = ""
        info = "ee" if "EE" in sampleNames else "#mu#mu" if "MM" in sampleNames else ""
        l = aux.Label(info="#scale[0.7]{%s}"%info, sim=((dataDoubleMuon not in additional)or(dataDoubleEG not in additional)or(dataHt not in additional)or(dataSF not in additional)))

        if binningName: binningName = "_"+binningName
        name = name.replace("/","__")
        saveName = "sameHistograms_{}_{}{}".format(sampleNames, name, binningName )
        #aux.save("DataMC_"+saveName )
        aux.save("DataMC_"+saveName,folder="plots_CutFlow/" )

        
def main():
    #bkgs=[DYjets,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma]
    #bkgs=[DYjetsLO,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma]
    bkgs=[DYjetsNLO,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma]
    #variables=["eta1","pt1","n_jets","n_vtx","phi1","m_ll","ht","n_photons","pt_g1","m_llg","pt_llg","mzg_exo","gammaMotherID"]
    variables=["cutflow_fine"]
    groups=["cutFlow_Fine_onZ"]
    for group in groups:
        for variable in variables:
            #drawSameHistogram("EE",group+"EE/"+variable, bkgs, additional=[dataDoubleEG],binning=binnings[variable],yTitle=labels[variable][1],xTitle=labels[variable][0])
            #drawSameHistogram("EE",group+"EE/"+variable, bkgs, additional=[dataDoubleEG],binning=binnings[variable],xTitle=labels[variable][0])
            drawSameHistogram("EE",group+"EE/"+variable, bkgs, additional=[dataDoubleEG],xTitle=labels[variable][0])
            drawSameHistogram("MM",group+"MM/"+variable, bkgs, additional=[dataDoubleMuon],xTitle=labels[variable][0])
            #drawSameHistogram("MM",group+"MM/"+variable, bkgs, additional=[dataDoubleMuon],binning=binnings[variable],yTitle=labels[variable][1],xTitle=labels[variable][0])


if __name__=="__main__":
    main()
