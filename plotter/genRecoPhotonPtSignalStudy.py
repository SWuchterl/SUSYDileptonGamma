from dataMC import labels, frange

import ROOT
from ROOT import *
from array import array
from include import *
import numpy as np
import pickle as pkl
import os

binnings = {

    'pt_g1':            frange(25, 126, 25),

}


def drawVR(massG, massN, binning=None, binningName="", xTitle=None, yTitle=None, sfZZ=[], sfDY=[], sfTT=[], sfWZ=[]):
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()

    style.divideByBinWidth = False

    style.minimumOne = True

    mG = str(massG)
    mN = str(massN)

    recoHist = aux.getFromFile(
        "/net/data_cms1b/user/swuchterl/tempTrees/AN/SMS-T5bbbbZg_" + mG + "_" + mN + "_myTuple.root", "reco")
    genHist = aux.getFromFile(
        "/net/data_cms1b/user/swuchterl/tempTrees/AN/SMS-T5bbbbZg_" + mG + "_" + mN + "_myTuple.root", "gen")

    recoHist = aux.rebin(recoHist, binning)
    genHist = aux.rebin(genHist, binning)

    aux.drawOpt(recoHist, "signal")
    aux.drawOpt(genHist, "signal")

    recoHist.SetLineColor(ROOT.kRed)

    m.add(recoHist, "reco")
    m.add(genHist, "gen")

    m.Draw()

    r = ratio.Ratio(
        "#scale[.9]{#lower[.24]{reco/gen}}", recoHist, genHist)
    rMax = 2.
    r.draw(0., rMax)

    directory = "plots/"
    if not os.path.exists(directory):
        os.makedirs(directory)

    aux.save("RecoGenPhotonPtStudy" + mG + "_" + mN, folder=directory)

    # aux.save(name.replace("/","_"),folder=directory)


binnings_ = binnings.copy()

binnings____ = frange(20, 1500, 10)
# binnings____ = frange(20, 40, 1)


drawVR(1200, 10, binning=binnings____)
drawVR(1200, 500, binning=binnings____)
