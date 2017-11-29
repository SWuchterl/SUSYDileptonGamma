import ROOT
from ROOT import *
from ROOT.TColor import *
from array import array
from include import *
import numpy as np

path = "../minimal/output/"

DY = ROOT.TFile(path+"DYJetsToLL_M-50_ext_hists.root")
DY.ls()
hDY = DY.Get("selMM/")
hDY.SetDirectory(0)
#ZG = ROOT.TFile(path+"ZGTo2LG_ext_hists.root")
ZG = ROOT.TFile(path+"DYJetsToLL_M-50-amcatnloFXFX_ext_hists.root")
hZG = ZG.Get("selMM/")
hZG.SetDirectory(0)
hDY=aux.rebin2d(hDY,range(0,300,2),range(0,400,2))
hZG=aux.rebin2d(hZG,range(0,300,2),range(0,400,2))
nDY = hDY.Integral()
nZG = hZG.Integral()
#print nDY,nZG


c = TCanvas("","",600,600)
c.SetRightMargin(0.18)
hDY.SetStats(0)
hDY.Draw("COLZ")
c.Print("DY.pdf")

c = TCanvas("","",600,600)
c.SetRightMargin(0.18)
hZG.SetStats(0)
hZG.Draw("COLZ")
c.Print("ZG.pdf")
#normalize
hDY.Scale(1./nZG*100.)
hZG.Scale(1./nZG*100.)
c = TCanvas("","",600,600)
c.SetRightMargin(0.18)
hDY.GetZaxis().SetRangeUser(0.,0.3)
hDY.SetStats(0)
hDY.Draw("COLZ")
c.Print("DY_norm.pdf")

c = TCanvas("","",600,600)
c.SetRightMargin(0.18)
hZG.GetZaxis().SetRangeUser(0.,0.3)
hZG.SetStats(0)
hZG.Draw("COLZ")
c.Print("ZG_norm.pdf")
