import ROOT
from ROOT import *
from array import array
from include import *
import numpy as np


import datasetsNoVeto

from dataMC import binnings,labels

def drawPhotonPT(sampleNames,name, binning=None, binningName="", scaleToData=False, xTitle=None, yTitle=None ):
    style.divideByBinWidth = True
    
    nBins =  binning
    
    dyLOVetoHist = aux.stdHist(DYjetsNLO,name+"PhotonPT_Veto",nBins)
    dyLONoVetoHist = aux.stdHist(DYjetsNLO,name+"PhotonPT_NoVeto",nBins)
    dyNLOVetoHist = aux.stdHist(DYjetsNLO,name+"PhotonPT_Veto",nBins)
    dyNLONoVetoHist = aux.stdHist(DYjetsNLO,name+"PhotonPT_NoVeto",nBins)
    zgHist = aux.stdHist(zgamma,name+"PhotonPT_NoVeto",nBins)
    
    dyLOVetoHist.SetLineColor(ROOT.kBlue-3) 
    dyLONoVetoHist.SetLineColor(ROOT.kBlue+3) 
    dyNLOVetoHist.SetLineColor(ROOT.kRed-3) 
    dyNLONoVetoHist.SetLineColor(ROOT.kRed+3) 
    zgHist.SetLineColor(ROOT.kGreen-3)
    
    aux.drawOpt(dyLOVetoHist,"signal")
    aux.drawOpt(dyLONoVetoHist,"signal")
    aux.drawOpt(dyNLOVetoHist,"signal")
    aux.drawOpt(dyNLONoVetoHist,"signal")
    aux.drawOpt(zgHist,"signal")
    
    dyLOVetoStat = aux.addHists(dyLOVetoHist)
    dyLONoVetoStat = aux.addHists(dyLONoVetoHist)
    dyNLOVetoStat = aux.addHists(dyNLOVetoHist)
    dyNLONoVetoStat = aux.addHists(dyNLONoVetoHist)
    zgStat = aux.addHists(zgHist)
    
    aux.drawOpt(dyLOVetoStat,"statUnc")
    aux.drawOpt(dyLONoVetoStat,"statUnc")
    aux.drawOpt(dyNLOVetoStat,"statUnc")
    aux.drawOpt(dyNLONoVetoStat,"statUnc")
    aux.drawOpt(zgStat,"statUnc")
    
    dyLOVetoStat.SetLineColor(ROOT.kBlue-3) 
    dyLONoVetoStat.SetLineColor(ROOT.kBlue+3) 
    dyNLOVetoStat.SetLineColor(ROOT.kRed-3) 
    dyNLONoVetoStat.SetLineColor(ROOT.kRed+3) 
    zgStat.SetLineColor(ROOT.kGreen-3)

    if yTitle:
        dyLOVetoStat.SetYTitle( yTitle )
    else:
        dyLOVetoStat.SetYTitle( aux.getYAxisTitle( dyLOVetoStat ) )
    if xTitle:
        dyLOVetoStat.SetXTitle( xTitle )
    
    c = ROOT.TCanvas()
    m = multiplot.Multiplot()
    #m.add(dyLOVetoHist,"DY LO Veto")
    #m.add(dyLONoVetoHist,"DY LO NoVeto")
    m.add(dyNLOVetoHist,"DY NLO Veto")
    m.add(dyNLONoVetoHist,"DY NLO NoVeto")
    m.add(zgHist,"ZG")
    m.add(zgStat,"") 
    #m.add(dyLOVetoStat,"") 
    #m.add(dyLONoVetoStat,"") 
    m.add(dyNLOVetoStat,"") 
    m.add(dyNLONoVetoStat,"") 
    
    m.maximum = 2.6*m.getMaximum()
    m.minimum = m.getMinimum()
    
    

    
    
    m.Draw()
    
    info = "ee" if "EE" in sampleNames else "#mu#mu" if "MM" in sampleNames else ""
    l = aux.Label(info="#scale[0.7]{%s}"%info, sim=True)

    if binningName: binningName = "_"+binningName
    name = name.replace("/","__")
    saveName = "photonPtVetoCompare_{}_{}{}".format(sampleNames, name, binningName )
    #aux.save("DataMC_"+saveName )
    aux.save("DataMC_"+saveName,folder="plots_compareZGDY/" )

def drawGenPhotonPT(sampleNames,name, binning=None, binningName="", scaleToData=False, xTitle=None, yTitle=None ):
    style.divideByBinWidth = True
    
    nBins =  binning
    
    dyLOVetoHist = aux.stdHist(DYjetsNLO,name+"genPhotonPT_Veto",nBins)
    dyLONoVetoHist = aux.stdHist(DYjetsNLO,name+"genPhotonPT_NoVeto",nBins)
    dyNLOVetoHist = aux.stdHist(DYjetsNLO,name+"genPhotonPT_Veto",nBins)
    dyNLONoVetoHist = aux.stdHist(DYjetsNLO,name+"genPhotonPT_NoVeto",nBins)
    zgHist = aux.stdHist(zgamma,name+"genPhotonPT_NoVeto",nBins)
    
    dyLOVetoHist.SetLineColor(ROOT.kBlue-3) 
    dyLONoVetoHist.SetLineColor(ROOT.kBlue+3) 
    dyNLOVetoHist.SetLineColor(ROOT.kRed-3) 
    dyNLONoVetoHist.SetLineColor(ROOT.kRed+3) 
    zgHist.SetLineColor(ROOT.kGreen-3)
    
    aux.drawOpt(dyLOVetoHist,"signal")
    aux.drawOpt(dyLONoVetoHist,"signal")
    aux.drawOpt(dyNLOVetoHist,"signal")
    aux.drawOpt(dyNLONoVetoHist,"signal")
    aux.drawOpt(zgHist,"signal")
    
    dyLOVetoStat = aux.addHists(dyLOVetoHist)
    dyLONoVetoStat = aux.addHists(dyLONoVetoHist)
    dyNLOVetoStat = aux.addHists(dyNLOVetoHist)
    dyNLONoVetoStat = aux.addHists(dyNLONoVetoHist)
    zgStat = aux.addHists(zgHist)
    
    aux.drawOpt(dyLOVetoStat,"statUnc")
    aux.drawOpt(dyLONoVetoStat,"statUnc")
    aux.drawOpt(dyNLOVetoStat,"statUnc")
    aux.drawOpt(dyNLONoVetoStat,"statUnc")
    aux.drawOpt(zgStat,"statUnc")
    
    dyLOVetoStat.SetLineColor(ROOT.kBlue-3) 
    dyLONoVetoStat.SetLineColor(ROOT.kBlue+3) 
    dyNLOVetoStat.SetLineColor(ROOT.kRed-3) 
    dyNLONoVetoStat.SetLineColor(ROOT.kRed+3) 
    zgStat.SetLineColor(ROOT.kGreen-3)

    if yTitle:
        dyLOVetoStat.SetYTitle( yTitle )
    else:
        dyLOVetoStat.SetYTitle( aux.getYAxisTitle( dyLOVetoStat ) )
    if xTitle:
        dyLOVetoStat.SetXTitle( xTitle )
    
    c = ROOT.TCanvas()
    m = multiplot.Multiplot()
    #m.add(dyLOVetoHist,"DY LO Veto")
    #m.add(dyLONoVetoHist,"DY LO NoVeto")
    m.add(dyNLOVetoHist,"DY NLO Veto")
    m.add(dyNLONoVetoHist,"DY NLO NoVeto")
    m.add(zgHist,"ZG")
    m.add(zgStat,"") 
    #m.add(dyLOVetoStat,"") 
    #m.add(dyLONoVetoStat,"") 
    m.add(dyNLOVetoStat,"") 
    m.add(dyNLONoVetoStat,"") 
    
    m.maximum = 2.6*m.getMaximum()
    m.minimum = m.getMinimum()
    
    

    
    
    m.Draw()
    
    info = "ee" if "EE" in sampleNames else "#mu#mu" if "MM" in sampleNames else ""
    l = aux.Label(info="#scale[0.7]{%s}"%info, sim=True)

    if binningName: binningName = "_"+binningName
    name = name.replace("/","__")
    saveName = "genPhotonPtVetoCompare_{}_{}{}".format(sampleNames, name, binningName )
    #aux.save("DataMC_"+saveName )
    aux.save("DataMC_"+saveName,folder="plots_compareZGDY/" )


def drawPhotonPTStack(sampleNames,name, binning=None, binningName="", scaleToData=False, xTitle=None, yTitle=None ):
    style.divideByBinWidth = True
    
    nBins =  binning
    
    dyLOVetoHist = aux.stdHist(DYjetsNLO,name+"PhotonPT_Veto",nBins)
    dyLONoVetoHist = aux.stdHist(DYjetsNLO,name+"PhotonPT_NoVeto",nBins)
    dyNLOVetoHist = aux.stdHist(DYjetsNLO,name+"PhotonPT_Veto",nBins)
    dyNLONoVetoHist = aux.stdHist(DYjetsNLO,name+"PhotonPT_NoVeto",nBins)
    dyNLOPTHist = aux.stdHist(DYjetsNLO,name+"pt_g1",nBins)
    zgHist = aux.stdHist(zgamma,name+"pt_g1",nBins)
    
    dyLOVetoHist.SetLineColor(ROOT.kBlue-3) 
    dyLONoVetoHist.SetLineColor(ROOT.kBlue+3) 
    dyNLOVetoHist.SetLineColor(ROOT.kRed-3) 
    dyNLONoVetoHist.SetLineColor(ROOT.kRed+3) 
    dyNLOPTHist.SetLineColor(ROOT.kYellow)
    zgHist.SetLineColor(ROOT.kGreen-3)
    
    #aux.drawOpt(dyLOVetoHist,"signal")
    #aux.drawOpt(dyLONoVetoHist,"signal")
    aux.drawOpt(dyNLOVetoHist,"signal")
    #aux.drawOpt(dyNLONoVetoHist,"signal")
    aux.drawOpt(dyNLOPTHist,"signal")
    #aux.drawOpt(zgHist,"signal")
    
    dyLOVetoStat = aux.addHists(dyLOVetoHist)
    dyLONoVetoStat = aux.addHists(dyLONoVetoHist)
    dyNLOVetoStat = aux.addHists(dyNLOVetoHist)
    dyNLONoVetoStat = aux.addHists(dyNLONoVetoHist)
    zgStat = aux.addHists(zgHist)
    dyNLOPTStat = aux.addHists(dyNLOPTHist)
    
    aux.drawOpt(dyLOVetoStat,"statUnc")
    aux.drawOpt(dyLONoVetoStat,"statUnc")
    aux.drawOpt(dyNLOVetoStat,"statUnc")
    aux.drawOpt(dyNLONoVetoStat,"statUnc")
    aux.drawOpt(zgStat,"statUnc")
    aux.drawOpt(dyNLOPTStat,"statUnc")
    
    dyLOVetoStat.SetLineColor(ROOT.kBlue-3) 
    dyLONoVetoStat.SetLineColor(ROOT.kBlue+3) 
    dyNLOVetoStat.SetLineColor(ROOT.kRed-3) 
    dyNLONoVetoStat.SetLineColor(ROOT.kRed+3) 
    zgStat.SetLineColor(ROOT.kGreen-3)
    dyNLOPTStat.SetLineColor(ROOT.kYellow)

    if yTitle:
        dyLOVetoStat.SetYTitle( yTitle )
    else:
        dyLOVetoStat.SetYTitle( aux.getYAxisTitle( dyLOVetoStat ) )
    if xTitle:
        dyLOVetoStat.SetXTitle( xTitle )
    
    c = ROOT.TCanvas()
    m = multiplot.Multiplot()
    #m.add(dyLOVetoHist,"DY LO Veto")
    #m.add(dyLONoVetoHist,"DY LO NoVeto")
    m.addStack(dyNLONoVetoHist,"DY NLO NoVeto")
    m.add(dyNLOVetoHist,"DY NLO Veto")
    m.add(dyNLOPTHist,"DY NLO")
    m.addStack(zgHist,"ZG")
    #m.add(zgStat,"") 
    #m.add(dyLOVetoStat,"") 
    #m.add(dyLONoVetoStat,"") 
    m.add(dyNLOVetoStat,"") 
    #m.add(dyNLONoVetoStat,"") 
    m.add(dyNLOPTStat,"") 
    
    
    
    m.maximum = 2.6*m.getMaximum()
    m.minimum = m.getMinimum()
    
    

    
    
    m.Draw()
    
    info = "ee" if "EE" in sampleNames else "#mu#mu" if "MM" in sampleNames else ""
    l = aux.Label(info="#scale[0.7]{%s}"%info, sim=True)

    if binningName: binningName = "_"+binningName
    name = name.replace("/","__")
    saveName = "genPhotonPtStack_{}_{}{}".format(sampleNames, name, binningName )
    #aux.save("DataMC_"+saveName )
    aux.save("DataMC_"+saveName,folder="plots_compareZGDY/" )


def main():
    drawPhotonPT("EE","sel"+"EE/",binning=binnings["PhotonPT_Veto"])
    drawPhotonPT("MM","sel"+"MM/",binning=binnings["PhotonPT_Veto"])
    drawGenPhotonPT("EE","sel"+"EE/",binning=binnings["genPhotonPT_Veto"])
    drawGenPhotonPT("MM","sel"+"MM/",binning=binnings["genPhotonPT_Veto"])
    drawPhotonPTStack("EE","sel"+"EE/",binning=binnings["PhotonPT_Veto"])
    drawPhotonPTStack("MM","sel"+"MM/",binning=binnings["PhotonPT_Veto"])

main()
