import ROOT
from ROOT import *

muonTrackScaleFactorEtaHisto = TH1F("muonTrackScaleFactorEtaHisto","muonTrackScaleFactorEta",12,0.,2.4);

muonTrackScaleFactorEtaHisto.SetBinContent(1,0.9970);
muonTrackScaleFactorEtaHisto.SetBinContent(2,0.9977);
muonTrackScaleFactorEtaHisto.SetBinContent(3,0.9981);
muonTrackScaleFactorEtaHisto.SetBinContent(4,0.9978);
muonTrackScaleFactorEtaHisto.SetBinContent(5,0.9980);
muonTrackScaleFactorEtaHisto.SetBinContent(6,0.9972);
muonTrackScaleFactorEtaHisto.SetBinContent(7,0.9962);
muonTrackScaleFactorEtaHisto.SetBinContent(8,0.9955);
muonTrackScaleFactorEtaHisto.SetBinContent(9,0.9958);
muonTrackScaleFactorEtaHisto.SetBinContent(10,0.9939);
muonTrackScaleFactorEtaHisto.SetBinContent(11,0.9929);
muonTrackScaleFactorEtaHisto.SetBinContent(12,0.9873);

f = ROOT.TFile("test.root","RECREATE")
muonTrackScaleFactorEtaHisto.Write()
f.Close()
