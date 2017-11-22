import ROOT
from ROOT import *
import numpy as np



complete = TFile("egammaEffi.txt_EGM2D.root","READ")
partial = TFile("scaleFactors.root","READ")

IDScaleFactorHisto = partial.Get("GsfElectronToMVATightTightIP2DSIP3D4")
IsoScaleFactorHisto = partial.Get("MVAVLooseElectronToMini")
ConvMissHitScaleFactorHisto = partial.Get("MVATightElectronToConvVetoIHit0")

comp = complete.Get("EGamma_SF2D")

ptPoints=np.arange(10.,200.,10.)
etaPoints=np.arange(0.,2.5,0.1)

#print ptPoints
#print etaPoints
for pt in ptPoints:
    for eta in etaPoints:
        compSF=comp.GetBinContent(comp.GetXaxis().FindBin(eta),comp.GetYaxis().FindBin(pt))
        compSFErr=comp.GetBinError(comp.GetXaxis().FindBin(eta),comp.GetYaxis().FindBin(pt))
        #partIDSF=IDScaleFactorHisto.GetBinContent(IDScaleFactorHisto.GetYaxis().FindBin(eta),IDScaleFactorHisto.GetXaxis().FindBin(pt))
        #partIsoSF=IsoScaleFactorHisto.GetBinContent(IsoScaleFactorHisto.GetYaxis().FindBin(eta),IsoScaleFactorHisto.GetXaxis().FindBin(pt))
        #partConvSF=ConvMissHitScaleFactorHisto.GetBinContent(ConvMissHitScaleFactorHisto.GetYaxis().FindBin(eta),ConvMissHitScaleFactorHisto.GetXaxis().FindBin(pt))
        partIDSF=IDScaleFactorHisto.GetBinContent(IDScaleFactorHisto.GetXaxis().FindBin(pt),IDScaleFactorHisto.GetYaxis().FindBin(eta))
        partIsoSF=IsoScaleFactorHisto.GetBinContent(IsoScaleFactorHisto.GetXaxis().FindBin(pt),IsoScaleFactorHisto.GetYaxis().FindBin(eta))
        partConvSF=ConvMissHitScaleFactorHisto.GetBinContent(ConvMissHitScaleFactorHisto.GetXaxis().FindBin(pt),ConvMissHitScaleFactorHisto.GetYaxis().FindBin(eta))
        partSF=partIDSF*partIsoSF*partConvSF
        diff=abs(compSF-partSF)
        print pt,eta,compSF,partSF
        #print pt,eta,round(abs(1.-compSF/partSF),2),round(compSFErr/compSF,2)
        #print pt,eta,abs(1.-compSF/partSF)/(compSFErr/compSF)
        
