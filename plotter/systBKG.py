from dataMC import labels, frange

import ROOT
from ROOT import *
from array import array
from include import *
import numpy as np
import pickle as pkl
import os

from CR_tt import calculateSFAndError, drawTTCR
from CR_DY import drawDYCR
from CR_WZ import drawCRWZ
from CR_ZZ import drawCRZZ

binnings = {
    # 'pt1':              frange(20,100,10)+frange(100,200,25)+range(200,350,50),
    # 'pt1':              frange(25,160,5),
    # 'pt1':              frange(20,150,10),
    # 'pt1':              frange(20,150,10),
    'pt1':              frange(20., 200., 10.),
    # 'pt2':              frange(20,100,10)+frange(100,200,25),
    'pt2':              frange(20, 150, 10),
    'pt3':              range(0, 200, 10),
    'pt4':              range(0, 200, 10),
    'eta1':             frange(0., 2.45, 0.1),
    # 'eta1':             frange(0., 2.4, 0.05),
    'eta2':             frange(0., 2.45, 0.1),
    # 'eta2':             frange(0., 2.4, 0.05),
    'eta3':             frange(0., 2.45, 0.1),
    'eta4':             frange(0., 2.45, 0.1),
    'phi1':             frange(0., 3.15, 0.1),
    'phi2':             frange(0., 3.15, 0.1),
    'phi3':             frange(0., 3.15, 0.1),
    'phi4':             frange(0., 3.15, 0.1),
    'ht':               frange(0., 1000., 50),
    # 'met':               [0,25,50,75,100,150,190,230,500],
    'met':               [0, 25, 50, 75, 100, 150, 250, 450],
    'm_ll':             frange(50., 100., 10) + frange(100., 300., 20.),
    'm_ll2':             frange(50., 100., 10) + frange(100., 300., 20.),
    'm_llg':             frange(0, 100, 10) + frange(100, 200, 10) + frange(200, 500, 50),
    'pt_llg':             frange(0, 100, 10) + frange(100, 200, 10) + frange(200, 500, 50),
    'n_jets':           frange(0., 7., 1),
    'n_photons':         frange(0., 4., 1),
    'n_vtx':            frange(0., 40., 1),
    # 'pt_g1':            frange(20,100,10)+frange(100,150,25)+frange(150,250,50),
    "pt_g1":            frange(25, 80, 5) + frange(80, 140, 10) + frange(140, 200, 20),
    'eta_g1':           frange(0., 2.60, 0.1),
    'phi_g1':           frange(0., 3.50, 0.1),
    'sigmaIetaIeta_g1': frange(0., 0.04, 0.01),
    'sigmaIphiIphi_g1': frange(0., 0.2, 0.01),
    'deltaR1_g1':       frange(0., 1., 0.05),
    'deltaR2_g1':       frange(0., 1., 0.05),
    'r9_g1':            frange(0., 1.5, 0.05),
    'hOverE_g1':        frange(0., 0.1, 0.01),
    'deltaEtaLL':       frange(0., 6., 0.1),
    'deltaPhiLL':       frange(0., 6., 0.1),
    'deltaEtaLLG':      frange(0., 6., 0.1),
    'deltaPhiLLG':      frange(0., 6., 0.01),
    'deltaRLL':         frange(0., 1., 0.1),
    'deltaRLLG':        frange(0., 6., 0.1),
    'st':               frange(0., 1000., 100),
    'stmet':            frange(0, 1000, 100) + frange(1000, 4000, 500),
    'zpt':              frange(0., 2000., 100),
    'mtll':             frange(0, 200, 25) + frange(200, 500, 50) + frange(500, 1250, 250),
    'mtllg':            frange(0, 200, 25) + frange(200, 500, 50) + frange(500, 1750, 250),
    'mtgmet':            frange(0, 200, 25) + frange(200, 500, 50) + frange(500, 1750, 250),
    'mtl1met':            frange(0, 200, 25) + frange(200, 500, 50) + frange(500, 1750, 250),
    'mtl2met':            frange(0, 200, 25) + frange(200, 500, 50) + frange(500, 1750, 250),
    'mtllmet':            frange(0, 200, 25) + frange(200, 500, 50) + frange(500, 1750, 250),
    'mtllgmet':            frange(0, 200, 25) + frange(200, 500, 50) + frange(500, 1750, 250),
    'mt2':            frange(0., 425., 25.),
    'mzg_exo':            frange(0, 100, 50) + frange(100, 200, 50) + frange(200, 500, 50),
    'gammaMotherID':            frange(0., 200., 1.),
    'genPhotonPT':            np.concatenate((np.arange(0, 100, 10), np.arange(100, 350, 50)), axis=0),
    'genPhotonPT_Veto':            np.concatenate((np.arange(0, 100, 10), np.arange(100, 350, 50)), axis=0),
    'PhotonPT_Veto':            np.concatenate((np.arange(0, 100, 10), np.arange(100, 350, 50)), axis=0),
    'genPhotonPT_NoVeto':            np.concatenate((np.arange(0, 100, 10), np.arange(100, 350, 50)), axis=0),
    'PhotonPT_NoVeto':            np.concatenate((np.arange(0, 100, 10), np.arange(100, 350, 50)), axis=0),
    'VetoCompare':            np.arange(0., 200., 1.),
    'DeltaPhiLLMet':            frange(0., 6., 0.3),
    'DeltaEtaLLMet':            frange(0., 6., 0.3),
    'DeltaRLLMet':            frange(0., 6., 0.3),
    'nElectrons': frange(0, 10, 1),
    'nMuons': frange(0, 10, 1),
    'mTL3Met': frange(0., 300., 10),
    'Fakes': frange(0, 5, 1)
}


pklZZ = pkl.load(open("plots_CR_zz/factors/CRZZ.pkl", "rb"))
pklDY = pkl.load(open("plots_CR_dy/factors/CRDY.pkl", "rb"))
pklTT = pkl.load(open("plots_CR_tt/factors/CRTT.pkl", "rb"))
pklWZ = pkl.load(open("plots_CR_wz/factors/CRWZ.pkl", "rb"))
#sfZZ,sfZZErr = pklZZ["LL"]["m_ll"]
#sfDY,sfDYErr = pklDY["LL"]["eta1"]
#sfTT = pklTT["EM"]["eta1"][0]
#sfWZ,sfWZErr = pklWZ["LL"]["eta1"]

toSave = {}
# toSave["wz"]={}
# toSave["zz"]={}
# toSave["dy"]={}
toSave_yield = {}
binnings_ = binnings.copy()


def doSyst(groups, variables, region, flavorComb):
    print "doSyst", groups, region, flavorComb
    toSave[region] = {}
    toSave_yield[region] = {}
    for group in groups:
        for variable in variables:
            # if(region in ["tt","ttg","dy","zg","zz","wz"]):
            if(region in ["tt", "ttg"]):
                toSave[region]["nom"] = drawTTCR(flavorComb, group + "/" + flavorComb + "/nom/" + variable, binning=binnings_[
                                                 variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["JESu"] = drawTTCR(flavorComb, group + "/" + flavorComb + "/JESu/" + variable, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["JESd"] = drawTTCR(flavorComb, group + "/" + flavorComb + "/JESd/" + variable, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["JERu"] = drawTTCR(flavorComb, group + "/" + flavorComb + "/JERu/" + variable, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["JERd"] = drawTTCR(flavorComb, group + "/" + flavorComb + "/JERd/" + variable, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["lepSFu"] = drawTTCR(flavorComb, group + "/" + flavorComb + "/lepSFu/" + variable, binning=binnings_[
                                                    variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["lepSFd"] = drawTTCR(flavorComb, group + "/" + flavorComb + "/lepSFd/" + variable, binning=binnings_[
                                                    variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["photonSFu"] = drawTTCR(flavorComb, group + "/" + flavorComb + "/photonSFu/" + variable, binning=binnings_[
                                                       variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["photonSFd"] = drawTTCR(flavorComb, group + "/" + flavorComb + "/photonSFd/" + variable, binning=binnings_[
                                                       variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["PUu"] = drawTTCR(flavorComb, group + "/" + flavorComb + "/PUu/" + variable, binning=binnings_[
                                                 variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["PUd"] = drawTTCR(flavorComb, group + "/" + flavorComb + "/PUd/" + variable, binning=binnings_[
                                                 variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["ISRu"] = drawTTCR(flavorComb, group + "/" + flavorComb + "/ISRu/" + variable, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISRUp", "topPt", "ewk"])
                toSave[region]["ISRd"] = drawTTCR(flavorComb, group + "/" + flavorComb + "/ISRd/" + variable, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISRDn", "topPt", "ewk"])
                toSave[region]["EWKu"] = drawTTCR(flavorComb, group + "/" + flavorComb + "/EWKu/" + variable, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewkUp"])
                toSave[region]["EWKd"] = drawTTCR(flavorComb, group + "/" + flavorComb + "/EWKd/" + variable, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewkDn"])
                toSave[region]["scale"] = {}
                for iPDF in range(9):
                    toSave[region]["scale"][str(iPDF)] = drawTTCR(flavorComb, group + "/" + flavorComb + "/" + str(iPDF) + "/" + variable, binning=binnings_[
                        variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk", "pdf" + str(iPDF)])
                toSave[region]["pdf"] = {}
                for iPDF in range(9, 110):
                    toSave[region]["pdf"][str(iPDF)] = drawTTCR(flavorComb, group + "/" + flavorComb + "/" + str(iPDF) + "/" + variable, binning=binnings_[
                        variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk", "pdf" + str(iPDF)])
                    # toSave[region]["pdf"][str(iPDF)]=drawTTCR(flavorComb,group+"/"+flavorComb+"/nom/"+variable,binning=binnings_[variable],xTitle=labels[variable][0],weightsToUse=["nISR","topPt","ewk"])

            if(region in ["zg", "dy"]):
                toSave[region]["nom"] = drawDYCR(flavorComb, group + "/" + flavorComb + "/nom/" + variable, dataLL, binning=binnings_[
                                                 variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["JESu"] = drawDYCR(flavorComb, group + "/" + flavorComb + "/JESu/" + variable, dataLL, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["JESd"] = drawDYCR(flavorComb, group + "/" + flavorComb + "/JESd/" + variable, dataLL, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["JERu"] = drawDYCR(flavorComb, group + "/" + flavorComb + "/JERu/" + variable, dataLL, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["JERd"] = drawDYCR(flavorComb, group + "/" + flavorComb + "/JERd/" + variable, dataLL, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["lepSFu"] = drawDYCR(flavorComb, group + "/" + flavorComb + "/lepSFu/" + variable, dataLL, binning=binnings_[
                                                    variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["lepSFd"] = drawDYCR(flavorComb, group + "/" + flavorComb + "/lepSFd/" + variable, dataLL, binning=binnings_[
                                                    variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["photonSFu"] = drawDYCR(flavorComb, group + "/" + flavorComb + "/photonSFu/" + variable, dataLL, binning=binnings_[
                                                       variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["photonSFd"] = drawDYCR(flavorComb, group + "/" + flavorComb + "/photonSFd/" + variable, dataLL, binning=binnings_[
                                                       variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["PUu"] = drawDYCR(flavorComb, group + "/" + flavorComb + "/PUu/" + variable, dataLL, binning=binnings_[
                                                 variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["PUd"] = drawDYCR(flavorComb, group + "/" + flavorComb + "/PUd/" + variable, dataLL, binning=binnings_[
                                                 variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["ISRu"] = drawDYCR(flavorComb, group + "/" + flavorComb + "/ISRu/" + variable, dataLL, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISRUp", "topPt", "ewk"])
                toSave[region]["ISRd"] = drawDYCR(flavorComb, group + "/" + flavorComb + "/ISRd/" + variable, dataLL, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISRDn", "topPt", "ewk"])
                toSave[region]["EWKu"] = drawDYCR(flavorComb, group + "/" + flavorComb + "/EWKu/" + variable, dataLL, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewkUp"])
                toSave[region]["EWKd"] = drawDYCR(flavorComb, group + "/" + flavorComb + "/EWKd/" + variable, dataLL, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewkDn"])
                toSave[region]["scale"] = {}
                for iPDF in range(9):
                    toSave[region]["scale"][str(iPDF)] = drawDYCR(flavorComb, group + "/" + flavorComb + "/" + str(iPDF) + "/" + variable, dataLL, binning=binnings_[
                        variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk", "pdf" + str(iPDF)])
                toSave[region]["pdf"] = {}
                for iPDF in range(9, 110):
                    toSave[region]["pdf"][str(iPDF)] = drawDYCR(flavorComb, group + "/" + flavorComb + "/" + str(iPDF) + "/" + variable, dataLL, binning=binnings_[
                        variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk", "pdf" + str(iPDF)])
                    # toSave[region]["pdf"][str(iPDF)]=drawDYCR(flavorComb,group+"/"+flavorComb+"/nom/"+variable,dataLL,binning=binnings_[variable],xTitle=labels[variable][0],weightsToUse=["nISR","topPt","ewk"])
            if(region in ["zz", "zz4l"]):
                toSave[region]["nom"] = drawCRZZ(flavorComb, group + "/" + flavorComb + "/nom/" + variable, dataLL, binning=binnings_[
                                                 variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["JESu"] = drawCRZZ(flavorComb, group + "/" + flavorComb + "/JESu/" + variable, dataLL, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["JESd"] = drawCRZZ(flavorComb, group + "/" + flavorComb + "/JESd/" + variable, dataLL, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["JERu"] = drawCRZZ(flavorComb, group + "/" + flavorComb + "/JERu/" + variable, dataLL, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["JERd"] = drawCRZZ(flavorComb, group + "/" + flavorComb + "/JERd/" + variable, dataLL, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["lepSFu"] = drawCRZZ(flavorComb, group + "/" + flavorComb + "/lepSFu/" + variable, dataLL, binning=binnings_[
                                                    variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["lepSFd"] = drawCRZZ(flavorComb, group + "/" + flavorComb + "/lepSFd/" + variable, dataLL, binning=binnings_[
                                                    variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["photonSFu"] = drawCRZZ(flavorComb, group + "/" + flavorComb + "/photonSFu/" + variable, dataLL, binning=binnings_[
                                                       variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["photonSFd"] = drawCRZZ(flavorComb, group + "/" + flavorComb + "/photonSFd/" + variable, dataLL, binning=binnings_[
                                                       variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["PUu"] = drawCRZZ(flavorComb, group + "/" + flavorComb + "/PUu/" + variable, dataLL, binning=binnings_[
                                                 variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["PUd"] = drawCRZZ(flavorComb, group + "/" + flavorComb + "/PUd/" + variable, dataLL, binning=binnings_[
                                                 variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["ISRu"] = drawCRZZ(flavorComb, group + "/" + flavorComb + "/ISRu/" + variable, dataLL, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISRUp", "topPt", "ewk"])
                toSave[region]["ISRd"] = drawCRZZ(flavorComb, group + "/" + flavorComb + "/ISRd/" + variable, dataLL, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISRDn", "topPt", "ewk"])
                toSave[region]["EWKu"] = drawCRZZ(flavorComb, group + "/" + flavorComb + "/EWKu/" + variable, dataLL, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewkUp"])
                toSave[region]["EWKd"] = drawCRZZ(flavorComb, group + "/" + flavorComb + "/EWKd/" + variable, dataLL, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewkDn"])
                toSave[region]["scale"] = {}
                for iPDF in range(9):
                    toSave[region]["scale"][str(iPDF)] = drawCRZZ(flavorComb, group + "/" + flavorComb + "/" + str(iPDF) + "/" + variable, dataLL, binning=binnings_[
                        variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk", "pdf" + str(iPDF)])
                toSave[region]["pdf"] = {}
                for iPDF in range(9, 110):
                    toSave[region]["pdf"][str(iPDF)] = drawCRZZ(flavorComb, group + "/" + flavorComb + "/" + str(iPDF) + "/" + variable, dataLL, binning=binnings_[
                        variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk", "pdf" + str(iPDF)])
                    # toSave[region]["pdf"][str(iPDF)]=drawCRZZ(flavorComb,group+"/"+flavorComb+"/nom/"+variable,dataLL,binning=binnings_[variable],xTitle=labels[variable][0],weightsToUse=["nISR","topPt","ewk"])
            if(region in ["wz"]):
                toSave[region]["nom"] = drawCRWZ(flavorComb, group + "/" + flavorComb + "/nom/" + variable, dataLL, binning=binnings_[
                                                 variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["JESu"] = drawCRWZ(flavorComb, group + "/" + flavorComb + "/JESu/" + variable, dataLL, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["JESd"] = drawCRWZ(flavorComb, group + "/" + flavorComb + "/JESd/" + variable, dataLL, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["JERu"] = drawCRWZ(flavorComb, group + "/" + flavorComb + "/JERu/" + variable, dataLL, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["JERd"] = drawCRWZ(flavorComb, group + "/" + flavorComb + "/JERd/" + variable, dataLL, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["lepSFu"] = drawCRWZ(flavorComb, group + "/" + flavorComb + "/lepSFu/" + variable, dataLL, binning=binnings_[
                                                    variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["lepSFd"] = drawCRWZ(flavorComb, group + "/" + flavorComb + "/lepSFd/" + variable, dataLL, binning=binnings_[
                                                    variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["photonSFu"] = drawCRWZ(flavorComb, group + "/" + flavorComb + "/photonSFu/" + variable, dataLL, binning=binnings_[
                                                       variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["photonSFd"] = drawCRWZ(flavorComb, group + "/" + flavorComb + "/photonSFd/" + variable, dataLL, binning=binnings_[
                                                       variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["PUu"] = drawCRWZ(flavorComb, group + "/" + flavorComb + "/PUu/" + variable, dataLL, binning=binnings_[
                                                 variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["PUd"] = drawCRWZ(flavorComb, group + "/" + flavorComb + "/PUd/" + variable, dataLL, binning=binnings_[
                                                 variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk"])
                toSave[region]["ISRu"] = drawCRWZ(flavorComb, group + "/" + flavorComb + "/ISRu/" + variable, dataLL, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISRUp", "topPt", "ewk"])
                toSave[region]["ISRd"] = drawCRWZ(flavorComb, group + "/" + flavorComb + "/ISRd/" + variable, dataLL, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISRDn", "topPt", "ewk"])
                toSave[region]["EWKu"] = drawCRWZ(flavorComb, group + "/" + flavorComb + "/EWKu/" + variable, dataLL, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewkUp"])
                toSave[region]["EWKd"] = drawCRWZ(flavorComb, group + "/" + flavorComb + "/EWKd/" + variable, dataLL, binning=binnings_[
                                                  variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewkDn"])
                toSave[region]["scale"] = {}
                for iPDF in range(9):
                    toSave[region]["scale"][str(iPDF)] = drawCRWZ(flavorComb, group + "/" + flavorComb + "/" + str(iPDF) + "/" + variable, dataLL, binning=binnings_[
                        variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk", "pdf" + str(iPDF)])
                toSave[region]["pdf"] = {}
                for iPDF in range(9, 110):
                    toSave[region]["pdf"][str(iPDF)] = drawCRWZ(flavorComb, group + "/" + flavorComb + "/" + str(iPDF) + "/" + variable, dataLL, binning=binnings_[
                        variable], xTitle=labels[variable][0], weightsToUse=["nISR", "topPt", "ewk", "pdf" + str(iPDF)])
                    # toSave[region]["pdf"][str(iPDF)]=drawCRWZ(flavorComb,group+"/"+flavorComb+"/nom/"+variable,dataLL,binning=binnings_[variable],xTitle=labels[variable][0],weightsToUse=["nISR","topPt","ewk"])

        allMC = zgamma + ttgamma + zz + wwgamma + wzgamma + \
            DYjetsNLO + wjets + tt + wz + ww + zz4l + wjets
        allMC.label = "MC mix"
        if(region in ["tt", "ttg", "dy", "zg", "zz", "wz"]):
            toSave_yield[region]["nom"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/nom", region, toSave[region]["nom"][0], weightsToUse=["nISR", "topPt", "ewk"])
            toSave_yield[region]["JESu"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/JESu", region, toSave[region]["JESu"][0], weightsToUse=["nISR", "topPt", "ewk"])
            toSave_yield[region]["JESd"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/JESd", region, toSave[region]["JESd"][0], weightsToUse=["nISR", "topPt", "ewk"])
            toSave_yield[region]["JERu"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/JERu", region, toSave[region]["JERu"][0], weightsToUse=["nISR", "topPt", "ewk"])
            toSave_yield[region]["JERd"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/JERd", region, toSave[region]["JERd"][0], weightsToUse=["nISR", "topPt", "ewk"])
            toSave_yield[region]["lepSFu"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/lepSFu", region, toSave[region]["lepSFu"][0], weightsToUse=["nISR", "topPt", "ewk"])
            toSave_yield[region]["lepSFd"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/lepSFd", region, toSave[region]["lepSFd"][0], weightsToUse=["nISR", "topPt", "ewk"])
            toSave_yield[region]["photonSFu"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/photonSFu", region, toSave[region]["photonSFu"][0], weightsToUse=["nISR", "topPt", "ewk"])
            toSave_yield[region]["photonSFd"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/photonSFd", region, toSave[region]["photonSFd"][0], weightsToUse=["nISR", "topPt", "ewk"])
            toSave_yield[region]["PUu"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/PUu", region, toSave[region]["PUu"][0], weightsToUse=["nISR", "topPt", "ewk"])
            toSave_yield[region]["PUd"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/PUd", region, toSave[region]["PUd"][0], weightsToUse=["nISR", "topPt", "ewk"])
            toSave_yield[region]["ISRu"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/ISRu", region, toSave[region]["ISRu"][0], weightsToUse=["nISRUp", "topPt", "ewk"])
            toSave_yield[region]["ISRd"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/ISRd", region, toSave[region]["ISRd"][0], weightsToUse=["nISRDn", "topPt", "ewk"])
            toSave_yield[region]["EWKu"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/EWKu", region, toSave[region]["EWKu"][0], weightsToUse=["nISR", "topPt", "ewkUp"])
            toSave_yield[region]["EWKd"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/EWKd", region, toSave[region]["EWKd"][0], weightsToUse=["nISR", "topPt", "ewkDn"])
            toSave_yield[region]["scale"] = {}
            for iPDF in range(9):
                toSave_yield[region]["scale"][str(iPDF)] = getYield("final_MC", allMC, "xx_0_0/sig/LL/" + str(
                    iPDF), region, toSave[region]["scale"][str(iPDF)][0], weightsToUse=["nISR", "topPt", "ewk", "pdf" + str(iPDF)])
            toSave_yield[region]["pdf"] = {}
            for iPDF in range(9, 110):
                toSave_yield[region]["pdf"][str(iPDF)] = getYield("final_MC", allMC, "xx_0_0/sig/LL/" + str(
                    iPDF), region, toSave[region]["pdf"][str(iPDF)][0], weightsToUse=["nISR", "topPt", "ewk", "pdf" + str(iPDF)])
                #toSave_yield[region]["pdf"][str(iPDF)]=getYield("final_MC", allMC, "xx_0_0/sig/LL/nom",region,toSave[region]["nom"][0],weightsToUse=["nISR","topPt","ewk"])
        else:
            toSave_yield[region]["nom"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/nom", region, 9999., weightsToUse=["nISR", "topPt", "ewk"])
            toSave_yield[region]["JESu"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/JESu", region, 9999, weightsToUse=["nISR", "topPt", "ewk"])
            toSave_yield[region]["JESd"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/JESd", region, 9999, weightsToUse=["nISR", "topPt", "ewk"])
            toSave_yield[region]["JERu"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/JERu", region, 9999, weightsToUse=["nISR", "topPt", "ewk"])
            toSave_yield[region]["JERd"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/JERd", region, 9999, weightsToUse=["nISR", "topPt", "ewk"])
            toSave_yield[region]["lepSFu"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/lepSFu", region, 9999, weightsToUse=["nISR", "topPt", "ewk"])
            toSave_yield[region]["lepSFd"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/lepSFd", region, 9999, weightsToUse=["nISR", "topPt", "ewk"])
            toSave_yield[region]["photonSFu"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/photonSFu", region, 9999, weightsToUse=["nISR", "topPt", "ewk"])
            toSave_yield[region]["photonSFd"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/photonSFd", region, 9999, weightsToUse=["nISR", "topPt", "ewk"])
            toSave_yield[region]["PUu"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/PUu", region, 9999, weightsToUse=["nISR", "topPt", "ewk"])
            toSave_yield[region]["PUd"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/PUd", region, 9999, weightsToUse=["nISR", "topPt", "ewk"])
            toSave_yield[region]["ISRu"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/ISRu", region, 9999, weightsToUse=["nISRUp", "topPt", "ewk"])
            toSave_yield[region]["ISRd"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/ISRd", region, 9999, weightsToUse=["nISRDn", "topPt", "ewk"])
            toSave_yield[region]["EWKu"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/EWKu", region, 9999, weightsToUse=["nISR", "topPt", "ewkUp"])
            toSave_yield[region]["EWKd"] = getYield(
                "final_MC", allMC, "xx_0_0/sig/LL/EWKd", region, 9999, weightsToUse=["nISR", "topPt", "ewkDn"])
            toSave_yield[region]["scale"] = {}
            for iPDF in range(9):
                toSave_yield[region]["scale"][str(iPDF)] = getYield("final_MC", allMC, "xx_0_0/sig/LL/" + str(
                    iPDF), region, 9999, weightsToUse=["nISR", "topPt", "ewk", "pdf" + str(iPDF)])
            toSave_yield[region]["pdf"] = {}
            for iPDF in range(9, 110):
                toSave_yield[region]["pdf"][str(iPDF)] = getYield("final_MC", allMC, "xx_0_0/sig/LL/" + str(
                    iPDF), region, 9999, weightsToUse=["nISR", "topPt", "ewk", "pdf" + str(iPDF)])
                #toSave_yield[region]["pdf"][str(iPDF)]=getYield("final_MC", allMC, "xx_0_0/sig/LL/nom",region,9999.,weightsToUse=["nISR","topPt","ewk"])
        # for key in toSave_yield[region]:
            #print key, toSave_yield[region][key]


def calculateEnvelope(dicti, region, unc):
    kappa = dicti
    tt1 = kappa[region]["nom"]["binMC_1"][0]
    tt2 = kappa[region]["nom"]["binMC_2"][0]
    # tt3=kappa[region]["nom"]["binMC_3"][0]
    tt1JESu = kappa[region][unc + "u"]["binMC_1"][0]
    tt2JESu = kappa[region][unc + "u"]["binMC_2"][0]
    # tt3JESu=kappa[region][unc+"u"]["binMC_3"][0]
    tt1JESd = kappa[region][unc + "d"]["binMC_1"][0]
    tt2JESd = kappa[region][unc + "d"]["binMC_2"][0]
    # tt3JESd=kappa[region][unc+"d"]["binMC_3"][0]

    if tt1 == 0. and tt2 != 0.:
        up2 = abs(tt2JESu - tt2) / abs(tt2)
        down2 = abs(tt2JESd - tt2) / abs(tt2)
        # return "tt1 0 ",max(abs(tt1JESu-tt1),abs(tt1JESd-tt1)),max(up2,down2)
        # return 0.,max(up2,down2)
        return 0., abs((tt2JESu - tt2JESd) / 2. / tt2)
    if tt2 == 0. and tt1 != 0.:
        up1 = abs(tt1JESu - tt1) / abs(tt1)
        down1 = abs(tt1JESd - tt1) / abs(tt1)
        # return "tt2 0 ",max(up1,down1),max(abs(tt2JESu-tt2),abs(tt2JESd-tt2))
        # return max(up1,down1),0.
        return abs((tt1JESu - tt1JESd) / 2. / tt1), 0.
    if tt1 == 0. and tt2 == 0.:
        # return "tt1/2 0 ",max(abs(tt1JESu-tt1),abs(tt1JESd-tt1)),max(abs(tt2JESu-tt2),abs(tt2JESd-tt2))
        return 0., 0.
    up1 = abs(tt1JESu - tt1) / abs(tt1)
    down1 = abs(tt1JESd - tt1) / abs(tt1)
    up2 = abs(tt2JESu - tt2) / abs(tt2)
    down2 = abs(tt2JESd - tt2) / abs(tt2)
    # return np.round(max(up1,down1),3),np.round(max(up2,down2),3)
    return abs((tt1JESu - tt1JESd) / 2. / tt1), abs((tt2JESu - tt2JESd) / 2. / tt2)
    # return abs((tt1JESu-tt1JESd)/2./tt1),abs((tt2JESu-tt2JESd)/2./tt2),abs((tt3JESu-tt3JESd)/2./tt3)


def calculateEnvelopeScale(dicti, region):
    kappa = dicti
    tt1 = kappa[region]["nom"]["binMC_1"][0]
    tt2 = kappa[region]["nom"]["binMC_2"][0]
    # tt3=kappa[region]["nom"]["binMC_3"][0]
    tt1Shift = [kappa[region]["scale"][str(c)]["binMC_1"][0] for c in range(9)]
    tt2Shift = [kappa[region]["scale"][str(c)]["binMC_2"][0] for c in range(9)]
    #tt3Shift=[kappa[region]["scale"][str(c)]["binMC_3"][0] for c in range(9)]

    dev1 = [abs(tt1Shift[i] - tt1) for i in range(9)]
    dev2 = [abs(tt2Shift[i] - tt2) for i in range(9)]
    #dev3=[abs(tt3Shift[i]-tt3) for i in range(9)]

    #print tt1
    #print tt1Shift

    #print dev1
    #print dev2

    max1 = max(dev1)
    max2 = max(dev2)
    # max3=max(dev3)

    if tt1 == 0. and tt2 != 0.:
        # return "tt1 0 ",max1,max2/tt2
        # return 0,max2/abs(tt2)
        return 0, abs((max(tt2Shift) - min(tt2Shift))) / 2. / abs(tt2)
    if tt2 == 0. and tt1 != 0.:
        # return "tt2 0 ",max1/tt1,max2
        # return max1/abs(tt1),0
        return abs((max(tt1Shift) - min(tt1Shift))) / 2. / abs(tt1), 0
    if tt1 == 0. and tt2 == 0.:
        # return "tt1/2 0 ",max1,max2
        return 0, 0
    # return max1/abs(tt1),max2/abs(tt2)
    return abs((max(tt1Shift) - min(tt1Shift))) / 2. / abs(tt1), abs((max(tt2Shift) - min(tt2Shift))) / 2. / abs(tt2)
    # return abs((max(tt1Shift)-min(tt1Shift)))/2./abs(tt1),abs((max(tt2Shift)-min(tt2Shift)))/2./abs(tt2),abs((max(tt3Shift)-min(tt3Shift)))/2./abs(tt3)


def calculateEnvelopePDF(dicti, region):
    kappa = dicti
    tt1 = kappa[region]["nom"]["binMC_1"][0]
    tt2 = kappa[region]["nom"]["binMC_2"][0]
    # tt3=kappa[region]["nom"]["binMC_3"][0]
    tt1Shift = [kappa[region]["pdf"]
                [str(c)]["binMC_1"][0] for c in range(9, 110)]
    tt2Shift = [kappa[region]["pdf"]
                [str(c)]["binMC_2"][0] for c in range(9, 110)]
    #tt3Shift=[kappa[region]["pdf"][str(c)]["binMC_3"][0] for c in range(9,110)]

    a1 = np.array(tt1Shift)
    a2 = np.array(tt2Shift)
    # a3=np.array(tt2Shift)

    mean1 = np.mean(a1)
    mean2 = np.mean(a2)
    # mean3=np.mean(a3)
    std1 = np.std(a1)
    std2 = np.std(a2)
    # std3=np.std(a3)

    # if tt1==0. and tt2!=0.:
    # return "tt1 0 ",max1,max2/tt2
    # if tt2==0. and tt1!=0.:
    # return "tt2 0 ",max1/tt1,max2
    # if tt1==0. and tt2==0.:
    # return "tt1/2 0 ",max1,max2
    if mean1 == 0 and mean2 != 0:
        return 0., abs(std2 / mean2)
    if mean1 != 0 and mean2 == 0:
        return abs(std1 / mean1), 0.
    if mean1 == 0 and mean2 == 0:
        # return 0., 0.
        return 0., 0.
    # return abs(std1/mean1),abs(std2/mean2),abs(std3/mean3)
    return abs(std1 / mean1), abs(std2 / mean2)


def main():
    bkgs = [DYjetsNLO, zgamma, tt, ttgamma, wwgamma,
            wzgamma, zz, wjets, wgamma, singletop, zz4l, wz, ww]

    doSyst(["CRTT"], ["eta1"], "tt", "EM")
    # doSyst(["nothing"],["eta1"],"ww","LL")
    # doSyst(["nothing"],["eta1"],"wzg","LL")
    # doSyst(["nothing"],["eta1"],"wwg","LL")
    # doSyst(["nothing"],["eta1"],"singletop","LL")
    # doSyst(["nothing"],["eta1"],"wjets","LL")
    # doSyst(["nothing"],["eta1"],"wgamma","LL")
    doSyst(["nothing"], ["eta1"], "other", "LL")
    doSyst(["CRTT"], ["eta1"], "ttg", "EM")
    doSyst(["CRWZ"], ["eta1"], "wz", "LL")
    doSyst(["CRDY"], ["eta1"], "dy", "LL")
    doSyst(["CRDY"], ["eta1"], "zg", "LL")
    doSyst(["CRZZ"], ["eta1"], "zz", "LL")
    doSyst(["CRZZ"], ["eta1"], "zz4l", "LL")
    if not os.path.exists("systBKG"):
        os.makedirs("systBKG")
    pkl.dump(toSave, open("systBKG/SF.pkl", "wb"))
    print "dumped systBKG/SF.pkl"
    pkl.dump(toSave_yield, open("systBKG/yields.pkl", "wb"))
    print "dumped systBKG/yields.pkl"
    kappa = pkl.load(open("systBKG/yields.pkl", "rb"))
    print "loaded systBKG/yields.pkl"
    # for key2 in kappa:
    #print key2,kappa
    #kappa=pkl.load(open( "systBKG/yields.pkl", "rb" ) )
    # regions=["tt","ttg","dy","zg","wz","singletop","wwg","wzg","ww","wjets","wgamma"]
    # regions=["tt","ttg","dy","zg","wz","singletop","wwg","wzg","ww","wjets"]
    # regions=["tt","ttg","dy","zg","wz","singletop","wwg","wzg","ww","wjets","zz","zz4l","wgamma"]
    regions = ["tt", "ttg", "dy", "zg", "wz", "zz", "zz4l", "other"]
    # regions=["tt","ttg","dy","zg","wz","zz","zz4l"]
    # regions=["tt","ttg","zg","wz","zz","zz4l"]
    # regions=["tt"]
    uncerts = ["JES", "JER", "PU", "lepSF", "photonSF", "ISR", "EWK"]
    # uncerts=["JES","JER","PU","lepSF","photonSF"]

    toSaveUncert = {}

    for region in regions:
        toSaveUncert[region] = {}
        for unc in uncerts:
            toSaveUncert[region][unc] = {}
            toSaveUncert[region][unc]["binMC_1"], toSaveUncert[region][unc]["binMC_2"] = calculateEnvelope(
                kappa, region, unc)
            # toSaveUncert[region][unc]["binMC_1"],toSaveUncert[region][unc]["binMC_2"],toSaveUncert[region][unc]["binMC_3"]=calculateEnvelope(kappa,region,unc)
            print region, unc, calculateEnvelope(kappa, region, unc)
        # toSaveUncert[region]["scale"]=calculateEnvelopeScale(kappa,region)
        toSaveUncert[region]["scale"] = {}
        toSaveUncert[region]["scale"]["binMC_1"], toSaveUncert[region]["scale"]["binMC_2"] = calculateEnvelopeScale(
            kappa, region)
        # toSaveUncert[region]["scale"]["binMC_1"],toSaveUncert[region]["scale"]["binMC_2"],toSaveUncert[region]["scale"]["binMC_3"]=calculateEnvelopeScale(kappa,region)
        print region, "scale", calculateEnvelopeScale(kappa, region)
        toSaveUncert[region]["pdf"] = {}
        toSaveUncert[region]["pdf"]["binMC_1"], toSaveUncert[region]["pdf"]["binMC_2"] = calculateEnvelopePDF(
            kappa, region)
        # toSaveUncert[region]["pdf"]["binMC_1"],toSaveUncert[region]["pdf"]["binMC_2"],toSaveUncert[region]["pdf"]["binMC_3"]=calculateEnvelopePDF(kappa,region)
        print region, "pdf", calculateEnvelopePDF(kappa, region)

    toSaveUncert["tt+ttg"] = {}
    toSaveUncert["zz+zz4l"] = {}
    toSaveUncert["dy+zg"] = {}
    uncerts2 = uncerts
    uncerts2.append("pdf")
    uncerts2.append("scale")
    for unc in uncerts2:
        toSaveUncert["tt+ttg"][unc] = {}
        toSaveUncert["tt+ttg"][unc]["binMC_1"] = (toSaveUncert["tt"][unc]["binMC_1"] * kappa["tt"]["nom"]["binMC_1"][0] + toSaveUncert["ttg"]
                                                  [unc]["binMC_1"] * kappa["ttg"]["nom"]["binMC_1"][0]) / (kappa["tt"]["nom"]["binMC_1"][0] + kappa["ttg"]["nom"]["binMC_1"][0])
        toSaveUncert["tt+ttg"][unc]["binMC_2"] = (toSaveUncert["tt"][unc]["binMC_2"] * kappa["tt"]["nom"]["binMC_2"][0] + toSaveUncert["ttg"]
                                                  [unc]["binMC_2"] * kappa["ttg"]["nom"]["binMC_2"][0]) / (kappa["tt"]["nom"]["binMC_2"][0] + kappa["ttg"]["nom"]["binMC_2"][0])
        print "tt+ttg", unc, toSaveUncert["tt+ttg"][unc]["binMC_1"], toSaveUncert["tt+ttg"][unc]["binMC_2"]
        #print "tt(before)",toSaveUncert["tt"][unc]["binMC_1"],toSaveUncert["tt"][unc]["binMC_2"]
        #print "tt(yield)",kappa["tt"]["nom"]["binMC_1"][0],kappa["tt"]["nom"]["binMC_2"][0]
        #print "ttg(before)",toSaveUncert["ttg"][unc]["binMC_1"],toSaveUncert["ttg"][unc]["binMC_2"]
        #print "ttg(yield)",kappa["ttg"]["nom"]["binMC_1"][0],kappa["ttg"]["nom"]["binMC_2"][0]
    toSaveUncert["dy+zg"] = {}
    for unc in uncerts:
        toSaveUncert["zz+zz4l"][unc] = {}
        toSaveUncert["zz+zz4l"][unc]["binMC_1"] = (toSaveUncert["zz"][unc]["binMC_1"] * kappa["zz"]["nom"]["binMC_1"][0] + toSaveUncert["zz4l"]
                                                   [unc]["binMC_1"] * kappa["zz4l"]["nom"]["binMC_1"][0]) / (kappa["zz"]["nom"]["binMC_1"][0] + kappa["zz4l"]["nom"]["binMC_1"][0])
        toSaveUncert["zz+zz4l"][unc]["binMC_2"] = (toSaveUncert["zz"][unc]["binMC_2"] * kappa["zz"]["nom"]["binMC_2"][0] + toSaveUncert["zz4l"]
                                                   [unc]["binMC_2"] * kappa["zz4l"]["nom"]["binMC_2"][0]) / (kappa["zz"]["nom"]["binMC_2"][0] + kappa["zz4l"]["nom"]["binMC_2"][0])
        print "zz+zz4l", unc, toSaveUncert["zz+zz4l"][unc]["binMC_1"], toSaveUncert["zz+zz4l"][unc]["binMC_2"]
    for unc in uncerts:
        toSaveUncert["dy+zg"][unc] = {}
        toSaveUncert["dy+zg"][unc]["binMC_1"] = (toSaveUncert["dy"][unc]["binMC_1"] * kappa["dy"]["nom"]["binMC_1"][0] + toSaveUncert["zg"]
                                                 [unc]["binMC_1"] * kappa["zg"]["nom"]["binMC_1"][0]) / (kappa["dy"]["nom"]["binMC_1"][0] + kappa["zg"]["nom"]["binMC_1"][0])
        toSaveUncert["dy+zg"][unc]["binMC_2"] = (toSaveUncert["dy"][unc]["binMC_2"] * kappa["dy"]["nom"]["binMC_2"][0] + toSaveUncert["zg"]
                                                 [unc]["binMC_2"] * kappa["zg"]["nom"]["binMC_2"][0]) / (kappa["dy"]["nom"]["binMC_2"][0] + kappa["zg"]["nom"]["binMC_2"][0])
        print "dy+zg", unc, toSaveUncert["dy+zg"][unc]["binMC_1"], toSaveUncert["dy+zg"][unc]["binMC_2"]
    pkl.dump(toSaveUncert, open("systBKG/unc.pkl", "wb"))
    #hm=pkl.load(open( "systBKG/SF.pkl", "rb" ) )
    #print hm["zg"]["nom"]
    #print hm["zg"]["scale"]


def getYield(name, dirSet, dirDir, BKGName, sfDict, weightsToUse=["nISR", "topPt", "ewk"]):
    style.divideByBinWidth = False
    #nBins = [0,25,50,75,100,150,250,450]
    nBins = [150, 200, 350]  # yes!
    #nBins = [100,150,200,350]

    dirHist = aux.stdHist(dirSet, dirDir + "/met", nBins)

    style.additionalPoissonUncertainty = False

    #zgHist = aux.stdHistWithoutNGen(zgamma, dirDir+"/met", nBins)
    #ttgHist = aux.stdHistWithoutNGen(ttgamma, dirDir+"/met", nBins)
    #zzHist = aux.stdHistWithoutNGen(zz, dirDir+"/met", nBins)
    #wwgHist = aux.stdHistWithoutNGen(wwgamma, dirDir+"/met", nBins)
    #wzgHist = aux.stdHistWithoutNGen(wzgamma, dirDir+"/met", nBins)
    #dyHist = aux.stdHistWithoutNGen(DYjetsNLO, dirDir+"/met", nBins)
    #wjetsHist = aux.stdHistWithoutNGen(wjets, dirDir+"/met", nBins)
    #ttHist = aux.stdHistWithoutNGen(tt, dirDir+"/met", nBins)
    ##singletopHist=aux.stdHistWithoutNGen(singletop, dirDir+"/met", nBins)
    #singletopHist=aux.stdHistWithoutNGen(singletop, "xx_0_0/sig/LL/nom"+"/met", nBins)
    #wzHist=aux.stdHistWithoutNGen(wz, dirDir+"/met", nBins)
    #wwHist=aux.stdHistWithoutNGen(ww, dirDir+"/met", nBins)
    #zz4lHist=aux.stdHistWithoutNGen(zz4l, dirDir+"/met", nBins)
    # wgHist=aux.stdHistWithoutNGen(wgamma,dirDir+"/met",nBins)
    zgHist = aux.stdHistWithoutNGenWithWeights(
        zgamma, dirDir + "/met", weightsToUse, nBins)
    ttgHist = aux.stdHistWithoutNGenWithWeights(
        ttgamma, dirDir + "/met", weightsToUse, nBins)
    zzHist = aux.stdHistWithoutNGenWithWeights(
        zz, dirDir + "/met", weightsToUse, nBins)
    wwgHist = aux.stdHistWithoutNGenWithWeights(
        wwgamma, dirDir + "/met", weightsToUse, nBins)
    wzgHist = aux.stdHistWithoutNGenWithWeights(
        wzgamma, dirDir + "/met", weightsToUse, nBins)
    dyHist = aux.stdHistWithoutNGenWithWeights(
        DYjetsNLO, dirDir + "/met", weightsToUse, nBins)
    wjetsHist = aux.stdHistWithoutNGenWithWeights(
        wjets, dirDir + "/met", weightsToUse, nBins)
    ttHist = aux.stdHistWithoutNGenWithWeights(
        tt, dirDir + "/met", weightsToUse, nBins)
    #singletopHist=aux.stdHistWithoutNGenWithWeights(singletop, dirDir+"/met",weightsToUse, nBins)
    singletopHist = aux.stdHistWithoutNGenWithWeights(
        singletop, "xx_0_0/sig/LL/nom" + "/met", weightsToUse, nBins)
    wzHist = aux.stdHistWithoutNGenWithWeights(
        wz, dirDir + "/met", weightsToUse, nBins)
    wwHist = aux.stdHistWithoutNGenWithWeights(
        ww, dirDir + "/met", weightsToUse, nBins)
    zz4lHist = aux.stdHistWithoutNGenWithWeights(
        zz4l, dirDir + "/met", weightsToUse, nBins)
    wgHist = aux.stdHistWithoutNGenWithWeights(
        wgamma, dirDir + "/met", weightsToUse, nBins)

    otherHist = aux.addHists(wwgHist, wzgHist, wjetsHist,
                             singletopHist, wwHist, wgHist)

    # Scaling
    #zg_AvgTopPtWeightHisto = zgamma.getHistWithoutNGen(dirDir+"/weight_topPt")
    #zg_AvgNIsrWeightHisto = zgamma.getHistWithoutNGen(dirDir+"/weight_nISR")
    #zg_AvgEWKinoWeightHisto = zgamma.getHistWithoutNGen(dirDir+"/weight_EWKinoPairPt")
    #zg_PDFWeightHisto = zgamma.getHistWithoutNGen(dirDir+"/weight_PDF")
    #ttg_AvgTopPtWeightHisto = ttgamma.getHistWithoutNGen(dirDir+"/weight_topPt")
    #ttg_AvgNIsrWeightHisto = ttgamma.getHistWithoutNGen(dirDir+"/weight_nISR")
    #ttg_AvgEWKinoWeightHisto = ttgamma.getHistWithoutNGen(dirDir+"/weight_EWKinoPairPt")
    #ttg_PDFWeightHisto = ttgamma.getHistWithoutNGen(dirDir+"/weight_PDF")
    #zz_AvgTopPtWeightHisto = zz.getHistWithoutNGen(dirDir+"/weight_topPt")
    #zz_AvgNIsrWeightHisto = zz.getHistWithoutNGen(dirDir+"/weight_nISR")
    #zz_AvgEWKinoWeightHisto = zz.getHistWithoutNGen(dirDir+"/weight_EWKinoPairPt")
    #zz_PDFWeightHisto = zz.getHistWithoutNGen(dirDir+"/weight_PDF")
    #wwg_AvgTopPtWeightHisto = wwgamma.getHistWithoutNGen(dirDir+"/weight_topPt")
    #wwg_AvgNIsrWeightHisto = wwgamma.getHistWithoutNGen(dirDir+"/weight_nISR")
    #wwg_AvgEWKinoWeightHisto = wwgamma.getHistWithoutNGen(dirDir+"/weight_EWKinoPairPt")
    #wwg_PDFWeightHisto = wwgamma.getHistWithoutNGen(dirDir+"/weight_PDF")
    #wzg_AvgTopPtWeightHisto = wzgamma.getHistWithoutNGen(dirDir+"/weight_topPt")
    #wzg_AvgNIsrWeightHisto = wzgamma.getHistWithoutNGen(dirDir+"/weight_nISR")
    #wzg_AvgEWKinoWeightHisto = wzgamma.getHistWithoutNGen(dirDir+"/weight_EWKinoPairPt")
    #wzg_PDFWeightHisto = wzgamma.getHistWithoutNGen(dirDir+"/weight_PDF")
    #dy_AvgTopPtWeightHisto = DYjetsNLO.getHistWithoutNGen(dirDir+"/weight_topPt")
    #dy_AvgNIsrWeightHisto = DYjetsNLO.getHistWithoutNGen(dirDir+"/weight_nISR")
    #dy_AvgEWKinoWeightHisto = DYjetsNLO.getHistWithoutNGen(dirDir+"/weight_EWKinoPairPt")
    #dy_PDFWeightHisto = DYjetsNLO.getHistWithoutNGen(dirDir+"/weight_PDF")
    #wjets_AvgTopPtWeightHisto = wjets.getHistWithoutNGen(dirDir+"/weight_topPt")
    #wjets_AvgNIsrWeightHisto = wjets.getHistWithoutNGen(dirDir+"/weight_nISR")
    #wjets_AvgEWKinoWeightHisto = wjets.getHistWithoutNGen(dirDir+"/weight_EWKinoPairPt")
    #wjets_PDFWeightHisto = wjets.getHistWithoutNGen(dirDir+"/weight_PDF")
    #tt_AvgTopPtWeightHisto = tt.getHistWithoutNGen(dirDir+"/weight_topPt")
    #tt_AvgNIsrWeightHisto = tt.getHistWithoutNGen(dirDir+"/weight_nISR")
    #tt_AvgEWKinoWeightHisto = tt.getHistWithoutNGen(dirDir+"/weight_EWKinoPairPt")
    #tt_PDFWeightHisto = tt.getHistWithoutNGen(dirDir+"/weight_PDF")
    #singletop_AvgTopPtWeightHisto = singletop.getHistWithoutNGen(dirDir+"/weight_topPt")
    #singletop_AvgNIsrWeightHisto = singletop.getHistWithoutNGen(dirDir+"/weight_nISR")
    #singletop_AvgEWKinoWeightHisto = singletop.getHistWithoutNGen(dirDir+"/weight_EWKinoPairPt")
    #singletop_PDFWeightHisto = singletop.getHistWithoutNGen(dirDir+"/weight_PDF")
    #singletop_AvgTopPtWeightHisto = singletop.getHistWithoutNGen("xx_0_0/sig/LL/nom"+"/weight_topPt")
    #singletop_AvgNIsrWeightHisto = singletop.getHistWithoutNGen("xx_0_0/sig/LL/nom"+"/weight_nISR")
    #singletop_AvgEWKinoWeightHisto = singletop.getHistWithoutNGen("xx_0_0/sig/LL/nom"+"/weight_EWKinoPairPt")
    #singletop_PDFWeightHisto = singletop.getHistWithoutNGen("xx_0_0/sig/LL/nom"+"/weight_PDF")
    #wz_AvgTopPtWeightHisto = wz.getHistWithoutNGen(dirDir+"/weight_topPt")
    #wz_AvgNIsrWeightHisto = wz.getHistWithoutNGen(dirDir+"/weight_nISR")
    #wz_AvgEWKinoWeightHisto = wz.getHistWithoutNGen(dirDir+"/weight_EWKinoPairPt")
    #wz_PDFWeightHisto = wz.getHistWithoutNGen(dirDir+"/weight_PDF")
    #ww_AvgTopPtWeightHisto = ww.getHistWithoutNGen(dirDir+"/weight_topPt")
    #ww_AvgNIsrWeightHisto = ww.getHistWithoutNGen(dirDir+"/weight_nISR")
    #ww_AvgEWKinoWeightHisto = ww.getHistWithoutNGen(dirDir+"/weight_EWKinoPairPt")
    #ww_PDFWeightHisto = ww.getHistWithoutNGen(dirDir+"/weight_PDF")
    #zz4l_AvgTopPtWeightHisto = zz4l.getHistWithoutNGen(dirDir+"/weight_topPt")
    #zz4l_AvgNIsrWeightHisto = zz4l.getHistWithoutNGen(dirDir+"/weight_nISR")
    #zz4l_AvgEWKinoWeightHisto = zz4l.getHistWithoutNGen(dirDir+"/weight_EWKinoPairPt")
    #zz4l_PDFWeightHisto = zz4l.getHistWithoutNGen(dirDir+"/weight_PDF")
    #wg_AvgTopPtWeightHisto = wgamma.getHistWithoutNGen(dirDir+"/weight_topPt")
    #wg_AvgNIsrWeightHisto = wgamma.getHistWithoutNGen(dirDir+"/weight_nISR")
    #wg_AvgEWKinoWeightHisto = wgamma.getHistWithoutNGen(dirDir+"/weight_EWKinoPairPt")
    #wg_PDFWeightHisto = wgamma.getHistWithoutNGen(dirDir+"/weight_PDF")

    histsToScale = [zgHist, ttgHist, zzHist, wwgHist, wzgHist, dyHist,
                    wjetsHist, ttHist, singletopHist, wzHist, wwHist, zz4lHist, wgHist]
    # topWeightHists=[zg_AvgTopPtWeightHisto,ttg_AvgTopPtWeightHisto,zz_AvgTopPtWeightHisto,wwg_AvgTopPtWeightHisto,wzg_AvgTopPtWeightHisto,dy_AvgTopPtWeightHisto,wjets_AvgTopPtWeightHisto,tt_AvgTopPtWeightHisto,singletop_AvgTopPtWeightHisto,wz_AvgTopPtWeightHisto,ww_AvgTopPtWeightHisto,zz4l_AvgTopPtWeightHisto,wg_AvgTopPtWeightHisto]
    # nISRWeightHists=[zg_AvgNIsrWeightHisto,ttg_AvgNIsrWeightHisto,zz_AvgNIsrWeightHisto,wwg_AvgNIsrWeightHisto,wzg_AvgNIsrWeightHisto,dy_AvgNIsrWeightHisto,wjets_AvgNIsrWeightHisto,tt_AvgNIsrWeightHisto,singletop_AvgNIsrWeightHisto,wz_AvgNIsrWeightHisto,ww_AvgNIsrWeightHisto,zz4l_AvgNIsrWeightHisto,wg_AvgNIsrWeightHisto]
    # EWKinoWeightHists=[zg_AvgEWKinoWeightHisto,ttg_AvgEWKinoWeightHisto,zz_AvgEWKinoWeightHisto,wwg_AvgEWKinoWeightHisto,wzg_AvgEWKinoWeightHisto,dy_AvgEWKinoWeightHisto,wjets_AvgEWKinoWeightHisto,tt_AvgEWKinoWeightHisto,singletop_AvgEWKinoWeightHisto,wz_AvgEWKinoWeightHisto,ww_AvgEWKinoWeightHisto,zz4l_AvgEWKinoWeightHisto,wg_AvgEWKinoWeightHisto]
    # PDFWeightHists=[zg_PDFWeightHisto,ttg_PDFWeightHisto,zz_PDFWeightHisto,wwg_PDFWeightHisto,wzg_PDFWeightHisto,dy_PDFWeightHisto,wjets_PDFWeightHisto,tt_PDFWeightHisto,singletop_PDFWeightHisto,wz_PDFWeightHisto,ww_PDFWeightHisto,zz4l_PDFWeightHisto,wg_PDFWeightHisto]

    # for i in range(len(histsToScale)):
    # if topWeightHists[i].Integral()>0.:
    # histsToScale[i].Scale(1./topWeightHists[i].GetMean())
    # if nISRWeightHists[i].Integral()>0.:
    # histsToScale[i].Scale(1./nISRWeightHists[i].GetMean())
    # if EWKinoWeightHists[i].Integral()>0.:
    # histsToScale[i].Scale(1./EWKinoWeightHists[i].GetMean())
    # if PDFWeightHists[i].Integral()>0.:
    # histsToScale[i].Scale(1./PDFWeightHists[i].GetMean())

    #pklZZ = pkl.load( open( "plots_CR_zz/factors/CRZZ.pkl", "rb" ) )
    #pklDY = pkl.load( open( "plots_CR_dy/factors/CRDY.pkl", "rb" ) )
    #pklTT = pkl.load( open( "plots_CR_tt/factors/CRTT.pkl", "rb" ) )
    #pklWZ = pkl.load( open( "plots_CR_wz/factors/CRWZ.pkl", "rb" ) )
    sfZZ, sfZZErr = pklZZ["LL"]["m_ll"]
    sfDY, sfDYErr = pklDY["LL"]["eta1"]
    sfTT = pklTT["EM"]["eta1"][0]
    sfWZ, sfWZErr = pklWZ["LL"]["eta1"]

    if(BKGName == "tt"):
        sfTT = sfDict
    if(BKGName == "zz"):
        sfZZ = sfDict
    if(BKGName == "wz"):
        sfWz = sfDict
    if(BKGName == "dy"):
        sfDY = sfDict

    zzHist.Scale(sfZZ)
    zz4lHist.Scale(sfZZ)
    dyHist.Scale(sfDY)
    zgHist.Scale(sfDY)
    wzHist.Scale(sfWZ)
    ttHist.Scale(sfTT)
    ttgHist.Scale(sfTT)

    #print "INT", ttgHist.Integral(),ttHist.Integral()

    dirHist = aux.addHists(*histsToScale)

    binDict = {}

    # for bin in range(dirHist.GetNbinsX()-1, dirHist.GetNbinsX()+1):
    for bin in range(dirHist.GetNbinsX() - 2, dirHist.GetNbinsX() + 1):
        binName = "bin{}_{}".format(name.split("_")[1], bin)
        bw = dirHist.GetBinWidth(bin) if style.divideByBinWidth else 1.
        #print (ttHist.GetBinLowEdge(bin))
        #print (ttHist.GetBinContent(bin)*bw)
        #print (ttgHist.GetBinLowEdge(bin))
        #print (ttgHist.GetBinContent(bin)*bw)
        if(BKGName == "ttg"):
            binDict[binName] = ttgHist.GetBinContent(bin) * bw,
        if(BKGName == "zg"):
            binDict[binName] = zgHist.GetBinContent(bin) * bw,
        if(BKGName == "zz"):
            binDict[binName] = zzHist.GetBinContent(bin) * bw,
        if(BKGName == "wwg"):
            binDict[binName] = wwgHist.GetBinContent(bin) * bw,
        if(BKGName == "wzg"):
            binDict[binName] = wzgHist.GetBinContent(bin) * bw,
        if(BKGName == "dy"):
            binDict[binName] = dyHist.GetBinContent(bin) * bw,
        if(BKGName == "wjets"):
            binDict[binName] = wjetsHist.GetBinContent(bin) * bw,
        if(BKGName == "tt"):
            binDict[binName] = ttHist.GetBinContent(bin) * bw,
            # binDict[binName]=ttHist.GetBinContent(bin)*bw+ttgHist.GetBinContent(bin)*bw,
        if(BKGName == "singletop"):
            binDict[binName] = singletopHist.GetBinContent(bin) * bw,
        if(BKGName == "wz"):
            binDict[binName] = wzHist.GetBinContent(bin) * bw,
        if(BKGName == "ww"):
            binDict[binName] = wwHist.GetBinContent(bin) * bw,
        if(BKGName == "wgamma"):
            binDict[binName] = wgHist.GetBinContent(bin) * bw,
        if(BKGName == "zz4l"):
            binDict[binName] = zz4lHist.GetBinContent(bin) * bw,
        if(BKGName == "other"):
            binDict[binName] = otherHist.GetBinContent(bin) * bw,
    return binDict


if __name__ == "__main__":
    main()
