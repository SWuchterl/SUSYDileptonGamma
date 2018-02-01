import ROOT
from ROOT import *
from array import array
from include import *
import numpy as np


import datasetsNoVeto

#from dataMC import binnings,labels
from dataMC import binnings,labels

def drawHistos(sampleNames, name, binning=None, binningName="", scaleToData=False, xTitle=None, yTitle=None ):
    style.divideByBinWidth = True
    #aux.drawOpt(dirHist, "data")
    
    nBins =  binning
    
    #print name
    
    zgHist = aux.stdHist(zgamma, name, nBins)
    dyLOHist = aux.stdHist(DYjetsLO, name, nBins)
    dyNLOHist = aux.stdHist(DYjetsNLO, name, nBins)

    zgHist.SetLineColor(ROOT.kGreen-3)
    dyLOHist.SetLineColor(ROOT.kBlue-7)
    dyNLOHist.SetLineColor(ROOT.kRed-7)

    aux.drawOpt(zgHist,"signal")
    aux.drawOpt(dyLOHist,"signal")
    aux.drawOpt(dyNLOHist,"signal")

    #aux.drawOpt(dyHist, "data")
    zgStat = aux.addHists(zgHist)
    dyLOStat = aux.addHists(dyLOHist)
    dyNLOStat = aux.addHists(dyNLOHist)
    
    aux.drawOpt(zgStat,"statUnc")
    aux.drawOpt(dyLOStat,"statUnc")    
    aux.drawOpt(dyNLOStat,"statUnc")    
    zgStat.SetLineColor(ROOT.kGreen-3)
    dyLOStat.SetLineColor(ROOT.kBlue-7)
    dyNLOStat.SetLineColor(ROOT.kRed-7)

    if yTitle:
        dyNLOStat.SetYTitle( yTitle )
    else:
        dyNLOStat.SetYTitle( aux.getYAxisTitle( dyNLOStat ) )
    if xTitle:
        dyNLOStat.SetXTitle( xTitle )
    
    c = ROOT.TCanvas()
    m = multiplot.Multiplot()
    m.add(dyLOHist,"DY LO")
    m.add(dyNLOHist,"DY NLO")
    m.add(zgHist,"ZG")
    m.add(zgStat,"") 
    m.add(dyLOStat,"") 
    m.add(dyNLOStat,"") 
    
    m.maximum = 2.6*m.getMaximum()
    m.minimum = m.getMinimum()
    
    

    
    
    m.Draw()
    
    info = "ee" if "EE" in sampleNames else "#mu#mu" if "MM" in sampleNames else ""
    l = aux.Label(info="#scale[0.7]{%s}"%info, sim=True)

    if binningName: binningName = "_"+binningName
    name = name.replace("/","__")
    saveName = "sameHistograms_{}_{}{}".format(sampleNames, name, binningName )
    #aux.save("DataMC_"+saveName )
    aux.save("DataMC_"+saveName,folder="plots_compareZGDY/" )


def main():
    #variables=["eta1","pt1","n_jets","n_vtx","phi1","m_ll","ht","n_photons","pt_g1","m_llg","pt_llg","mzg_exo","gammaMotherID"]
    #variables=["eta1","pt1","n_jets","phi1","m_ll","ht","n_photons","pt_g1","m_llg","gammaMotherID","sigmaIetaIeta_g1","sigmaIphiIphi_g1","r9_g1","hOverE_g1","deltaR1_g1","deltaR2_g1","eta_g1","phi_g1"]
    #variables=["pt_g1","genPhotonPT","genPhotonPT_Veto","PhotonPT_Veto","genPhotonPT_NoVeto","PhotonPT_NoVeto","VetoCompare"]
    variables=["pt_g1","genPhotonPT","genPhotonPT_Veto","PhotonPT_Veto","genPhotonPT_NoVeto","PhotonPT_NoVeto"]
    #variables=["eta1"]
    #variables=["eta1","pt1","n_jets","n_vtx","phi1","m_ll","ht","n_photons","pt_g1","sigmaIetaIeta_g1","sigmaIphiIphi_g1","r9_g1","hOverE_g1","deltaR1_g1","deltaR2_g1","eta_g1","phi_g1"]
    #groups=["sel","onZ","dilep","onZG","mllG110"]
    #groups=["sel","onZ","dilep","onZG","exo"]
    #groups=["sel","onZ","onZG","exo"]
    groups=["sel"]
    for group in groups:
        #print group
        for variable in variables:
            #print variable
            drawHistos("EE",group+"EE/"+variable,binning=binnings[variable],xTitle=labels[variable][0])
            drawHistos("MM",group+"MM/"+variable,binning=binnings[variable],xTitle=labels[variable][0])
            #drawHistos("EE",group+"EE/"+variable,binning=binnings[variable],xTitle=labels[variable[0])

main()
