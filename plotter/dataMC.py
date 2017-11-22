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
                     'sigmaIetaIeta_g1'
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
    'sigmaIetaIeta_g1': np.arange(0., 0.04,0.001)
}
labels = {
    'pt1': ["p_{T}^{leading}[GeV]","Events / 30 GeV"],
    'pt2': ["p_{T}^{trailing}[GeV]"," Events / 30 GeV"],
    'eta1': ["#Eta_{leading}","Events / 0.1"],
    'eta2': ["#Eta_{trailing}","Events / 0.1"],
    'phi1': ["#Phi_{leading}","Events / 0.1"],
    'phi2': ["#Phi_{trailing}","Events / 0.1"],
    'ht': ["H_{T} [GeV]", "Events / 10 GeV"],
    'met': ["E_{T}^{miss} [GeV]", "Events / 10 GeV "],
    'gen_ht': ["H_{T}^{gen} [GeV]", "Events / 10 GeV"],
    'm_ll': ["m_{ll} from all leptons [GeV]" ,"Events / 65 GeV"],
    'n_jets': ["N Jets", "Events / 1"],
    'n_photons': ["N Photons", "Events / 1"],
    'n_vtx': ["N Vertices", "Events / 1"],
    'pt_g1': ["p_{T} of leading photon","Events / 10 GeV"],
    'eta_g1': ["#Eta of leading photon","Events / 0.1"],
    'phi_g1': ["#phi of leading photon","Events / 0.1"],
    'sigmaIetaIeta_g1': ["#sigma_{I#etaI#eta} of leading photon","Events / 0.001"]
}



def drawSameHistogram( sampleNames, name, bkg=[], additional=[], binning=None, binningName="", scaleToData=False ):
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()

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
        h.SetYTitle( aux.getYAxisTitle( h ) )
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
            # disable errors for data, so that ErrorOption is working
            # TODO: kPoisson also for rebinned and scaled histograms
            if not(binning.any()): h.Sumw2(False)
            h.SetBinErrorOption( ROOT.TH1.kPoisson )
            dataHist = h
        else:
            h.drawOption_ = "hist"
            h.SetLineWidth(3)

        m.add( h, d.label )

    m.sortStackByIntegral()
    if m.Draw():

        # ratio
        hsm = m.hists[0].GetStack().Last()
        if dataHist:
            r = ratio.Ratio( "Data/MC", dataHist, hsm )
            r.draw(0.5,1.5)

        info = ""
        l = aux.Label(info="#scale[0.7]{%s}"%info, sim=((dataDoubleMuon not in additional)or(dataDoubleEG not in additional)))

        if binningName: binningName = "_"+binningName
        name = name.replace("/","__")
        saveName = "sameHistograms_{}_{}{}".format(sampleNames, name, binningName )
        aux.save( saveName )

        
def main():
    bkgs=[DYjets,zgamma,wwgamma,wzgamma,ttgamma,zz,tt]
    #variables=["eta1","pt1","n_jets","n_vtx","phi1","eta2","phi2","m_ll","ht","n_photons","pt2"]
    variables=["eta1","pt1","n_jets","n_vtx","phi1","m_ll","ht","n_photons"]
    groups=["dilep","sel"]
    for group in groups:
        for variable in variables:
            drawSameHistogram("EE",group+"EE/"+variable, bkgs, additional=[dataDoubleEG],binning=binnings[variable])
            drawSameHistogram("MM",group+"MM/"+variable, bkgs, additional=[dataDoubleMuon],binning=binnings[variable])
            #drawSameHistogram("SF",group+"/"+variable, bkgs, additional=[dataDoubleSF],binning=binnings[variable])
            #drawSameHistogram("HT",group+"/"+variable, bkgs, additional=[dataHt],binning=binnings[variable])

main()