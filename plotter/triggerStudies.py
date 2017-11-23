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
    'pt1':              np.arange(0., 300., 10.),
    'pt2':              np.arange(0., 300, 10),
    'eta1':             np.arange(0., 2.6, 0.1),
    'eta2':             np.arange(0., 2.60, 0.1),
    'phi1':             np.arange(0., 3.50, 0.10),
    'phi2':             np.arange(0., 3.50, 0.10),
    'ht':               np.arange(0., 1000.,10.),
    'met':              np.arange(0., 1000.,10.),
    'gen_ht':           np.arange(0., 1000.,10.),
    'm_ll':             np.arange(0., 650.,10.),
    'm_ll_e':           np.arange(0., 650.,10.),
    'm_ll_m':           np.arange(0., 650.,10.),
    'n_jets':           np.arange(0.,10.,1.),
    'n_photons':           np.arange(0.,10.,1.),
    'n_vtx':            np.arange(0.,40.,4.),
    #'pt_g1':            np.array((range(0,100,10)+np.array([100., 150., 200.,300,400.]))),
    'pt_g1':            np.concatenate((np.arange(0,100,10),np.arange(100,450,50)),axis=0),
    'eta_g1':           np.arange(0., 2.60, 0.10),
    'phi_g1':           np.arange(0., 3.50, 0.10),
    'sigmaIetaIeta_g1': np.arange(0., 0.04,0.001),
    'DeltaEtaLL':       np.arange(0.,6.,0.1),
    'DeltaPhiLL':       np.arange(0.,6.,0.1),
    'DeltaEtaLLG':      np.arange(0.,6.,0.1),
    'DeltaPhiLLG':      np.arange(0.,6.,0.1),
    'DeltaRLL':         np.arange(0.,6.,0.1),
    'DeltaRLLG':        np.arange(0.,6.,0.1),
    #'st':        np.arange(0.,2000.,100.),
    #'stmet':        np.arange(0.,2000.,100.),
    #'zpt':        np.arange(0.,500.,50.),
}
labels = {
    'pt1': ["p_{T}^{leading}[GeV]","Events / 30 GeV"],
    'pt2': ["p_{T}^{trailing}[GeV]"," Events / 30 GeV"],
    'eta1': ["|#Eta_{leading}|","Events / 0.1"],
    'eta2': ["|#Eta_{trailing}|","Events / 0.1"],
    'phi1': ["|#Phi_{leading}|","Events / 0.1"],
    'phi2': ["|#Phi_{trailing}|","Events / 0.1"],
    'ht': ["H_{T} [GeV]", "Events / 10 GeV"],
    'met': ["E_{T}^{miss} [GeV]", "Events / 10 GeV "],
    'gen_ht': ["H_{T}^{gen} [GeV]", "Events / 10 GeV"],
    'm_ll': ["m_{ll} from all leptons [GeV]" ,"Events / 65 GeV"],
    'n_jets': ["N Jets", "Events / 1"],
    'n_photons': ["N Photons", "Events / 1"],
    'n_vtx': ["N Vertices", "Events / 1"],
    'pt_g1': ["p_{T} of leading photon","Events / 10 GeV"],
    'eta_g1': ["|#Eta| of leading photon","Events / 0.1"],
    'phi_g1': ["|#phi| of leading photon","Events / 0.1"],
    'sigmaIetaIeta_g1': ["#sigma_{I#etaI#eta} of leading photon","Events / 0.001"],
    'DeltaEtaLL':       ["#Delta#Eta_{ll}","Events / 0.1"],
    'DeltaPhiLL':       ["#Delta#Phi_{ll}","Events / 0.1"],
    'DeltaEtaLLG':      ["#Delta#Eta_{ll,#gamma}","Events / 0.1"],
    'DeltaPhiLLG':      ["#Delta#Phi_{ll,#gamma}","Events / 0.1"],
    'DeltaRLL':         ["#DeltaR_{ll}","Events / 0.1"],
    'DeltaRLLG':        ["#DeltaR_{ll,#gamma}","Events / 0.1"],    
    #'st':        ["#itS_T [GeV]","Events / 100 GeV"],    
    #'stmet':        ["#it S_T+E_{T}^{miss}" [GeV],"Events / 100 GeV"],    
    #'zpt':        ["Z_{p_T} [GeV]","Events / 50. GeV"],    
}




def main():
    bkgs=[DYjets,zgamma,wwgamma,wzgamma,ttgamma,zz,tt,wjets]
    #variables=["eta1","pt1","n_jets","n_vtx","phi1","m_ll","ht","n_photons"]
    #groups=["dilep","sel","onZ"]
    for group in groups:
        for variable in variables:


main()
