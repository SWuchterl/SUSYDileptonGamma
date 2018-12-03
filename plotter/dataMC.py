import ROOT
from ROOT import *
from array import array
from include import *
import numpy as np


def frange(start, end, step):
    a = []
    tmp = start
    # while(tmp < end):
    while(tmp <= end):
        a.append(tmp)
        tmp += step
    return a


def frangeN(start, end, N):
    a = []
    tmp = start
    step = abs(end - start) / N
    while(tmp <= end):
        a.append(tmp)
        tmp += step
    return a


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
    # 'pt1':              np.arange(0., 500., 10.),
    # 'pt1':              np.concatenate((np.arange(0,200,10),np.arange(200,1750,50)),axis=0),
    # 'pt1':              np.concatenate((np.arange(0,200,10),np.arange(200,350,50)),axis=0),
    # 'pt1':              np.concatenate((np.arange(0,100,10),np.arange(100,200,25),np.arange(200,350,50)),axis=0),
    'pt1':              range(0, 100, 10) + range(100, 200, 25) + range(200, 350, 50),
    # 'pt1':              range(0,100,25)+range(100,200,25)+range(200,350,50),
    # 'pt1':              range(0,400,10),
    # 'pt2':              np.arange(0., 300, 10),
    'pt2':              range(0, 300, 10),
    # 'eta1':             np.arange(0., 2.6, 0.1),
    'eta1':             frange(0., 2.6, 0.1),
    # 'eta1':             np.arange(0., 2.6, 0.05),
    # 'eta2':             np.arange(0., 2.6, 0.1),
    'eta2':             frange(0., 2.6, 0.1),
    # 'phi1':             np.arange(0., 3.50, 0.10),
    # 'phi1':             np.arange(0., 3.50, 0.35),
    'phi1':             frange(0., 3.50, 0.1),
    # 'phi2':             np.arange(0., 3.50, 0.10),
    'phi2':             frange(0., 3.50, 0.1),
    # 'ht':               np.arange(0., 5000.,500.),
    'ht':               frange(0., 5000., 100),
    # 'met':              np.arange(0., 1000.,100.),
    # 'met':              np.arange(0., 1000.,500.),
    # 'met':              np.concatenate((np.arange(0,400,50),np.arange(400,3000,1300)),axis=0),
    # 'met':              np.arange(0., 301.,25.),
    # 'met':              np.concatenate((np.arange(0,200,25),np.arange(200,300,100),np.arange(300,3000,1400)),axis=0),
    # 'met':              np.concatenate((np.arange(0,200,25),np.arange(200,300,100),np.arange(300,700,200)),axis=0),
    # 'met':               range(0,200,25)+[200,300,500],
    'met':               [0, 25, 50, 75, 100, 150, 190, 230, 500],
    # 'met':              np.concatenate((np.linspace(0,100,5),np.arange(200,301,100)),axis=0),
    # 'met':              np.concatenate((np.linspace(0,100,5),np.arange(100,301,200)),axis=0),
    # 'gen_ht':           np.arange(0., 5000.,100.),
    # 'm_ll':             np.arange(50., 130.,2.),
    # 'm_ll':             frange(0., 350.,2.),
    # 'm_ll':             frange(0., 300.,10.),
    # 'm_ll':             frange(50.,100.,10)+frange(100., 300.,20.),
    'm_ll':             frange(50., 130., 2),
    'm_ll2':             frange(50., 100., 10) + frange(100., 300., 20.),
    # 'm_ll':             np.arange(80., 102.,2.),
    # 'm_ll':             np.arange(50., 300.,2.),
    # 'm_ll':             np.arange(50., 300.,10.),
    # 'm_ll':             np.arange(0., 650.,2.),
    # 'm_ll':             np.concatenate((np.arange(0,100,10),np.arange(100,200,10),np.arange(200,500,50)),axis=0),
    # 'm_{ll#gamma}':             np.concatenate((np.arange(0,100,10),np.arange(100,200,10),np.arange(200,500,50)),axis=0),
    # 'm_llg':             np.concatenate((np.arange(0,100,10),np.arange(100,200,10),np.arange(200,500,50)),axis=0),
    'm_llg':             frange(0, 100, 50) + frange(100, 200, 50) + frange(200, 500, 50),
    # 'pt_llg':             np.concatenate((np.arange(0,100,10),np.arange(100,200,10),np.arange(200,500,50)),axis=0),
    'pt_llg':             frange(0, 100, 10) + frange(100, 200, 10) + frange(200, 500, 50),
    # 'm_ll':             np.concatenate((np.arange(0,200,10),np.arange(200,950,50)),axis=0),
    # 'm_ll':             np.arange(80,110,1),
    # 'm_ll_e':           np.arange(0., 650.,10.),
    # 'm_ll_m':           np.arange(0., 650.,10.),
    'n_jets':           frange(0., 10., 1),
    'n_bjets':           frange(0., 10., 1),
    'n_photons':         frange(0., 10., 1),
    'n_vtx':            frange(0., 40., 1),
    # 'pt_g1':            np.array((range(0,100,10)+np.array([100., 150., 200.,300,400.]))),
    # 'pt_g1':            np.concatenate((np.arange(0,100,10),np.arange(100,950,50)),axis=0),
    # 'pt_g1':            np.concatenate((np.arange(0,100,10),np.arange(100,350,50)),axis=0),
    'pt_g1':            frange(20, 100, 10) + frange(100, 150, 25) + frange(150, 250, 50),
    'eta_g1':           frange(0., 2.60, 0.1),
    'phi_g1':           frange(0., 3.50, 0.1),
    'sigmaIetaIeta_g1': frange(0., 0.04, 0.001),
    'sigmaIphiIphi_g1': frange(0., 0.2, 0.01),
    'deltaR1_g1':       frange(0., 1., 0.05),
    'deltaR2_g1':       frange(0., 1., 0.05),
    'r9_g1':            frange(0., 1.5, 0.05),
    'hOverE_g1':        frange(0., 0.1, 0.01),
    'deltaEtaLL':       frange(0., 6., 0.1),
    'deltaPhiLL':       frange(0., 3.2, 0.2),
    'deltaEtaLL_neg':       frange(-6., 6., 0.1),
    'deltaPhiLL_neg':       frange(-6., 6., 0.1),
    'deltaEtaLLG':      frange(0., 6., 0.1),
    'deltaPhiLLG':      frange(0., 3.2, 0.2),
    # 'deltaRLL':         frange(0.,6.,0.12),
    'deltaRLL':         frange(0., 1., 0.1),
    'deltaRLLG':        frange(0., 6., 0.1),
    # 'st':               np.arange(0.,1000.,100.),
    'st':               frange(0., 250., 250) + frange(250., 751., 500.),
    # 'stmet':           np.arange(0.,5000.,100.),
    # 'stmet':            np.concatenate((np.arange(0,1000,100),np.arange(1000,4000,1000)),axis=0),
    'stmet':            frange(0, 250, 50) + frange(250, 1500, 250),
    # 'zpt':              np.arange(0.,2000.,100.),
    # 'zpt':              frange(0.,100.,50)+frange(100,200,100)+frange(200,500,300),
    'zpt':              frange(0., 450., 50),
    # 'mtll':             np.arange(0.,3000.,100.),
    # 'mtll':             np.concatenate((np.arange(0,200,25),np.arange(200,500,50),np.arange(500,1250,250)),axis=0),
    'mtll':             frange(0, 300, 100) + frange(300, 701, 400),
    # 'mtllg':            np.arange(0.,3000.,100.),
    # 'mtllg':            np.concatenate((np.arange(0,200,25),np.arange(200,500,50),np.arange(500,1750,250)),axis=0),
    'mtllg':            frange(0, 300, 100) + frange(300, 500, 200) + frange(500, 1001, 500),
    # 'mtgmet':            np.concatenate((np.arange(0,200,25),np.arange(200,500,50),np.arange(500,1750,250)),axis=0),
    'mtgmet':            frange(0, 200, 25) + frange(200, 500, 50) + frange(500, 1750, 250),
    # 'mtl1met':            np.concatenate((np.arange(0,200,25),np.arange(200,500,50),np.arange(500,1750,250)),axis=0),
    'mtl1met':            frange(0, 200, 100) + frange(200, 500, 100) + frange(500, 1001, 500),
    # 'mtl2met':            np.concatenate((np.arange(0,200,25),np.arange(200,500,50),np.arange(500,1750,250)),axis=0),
    'mtl2met':            frange(0, 200, 100) + frange(200, 500, 100) + frange(500, 1001, 500),
    # 'mtllmet':            np.concatenate((np.arange(0,200,25),np.arange(200,500,50),np.arange(500,1750,250)),axis=0),
    'mtllmet':            frange(0, 200, 25) + frange(200, 500, 50) + frange(500, 1750, 250),
    # 'mtllgmet':            np.concatenate((np.arange(0,200,25),np.arange(200,500,50),np.arange(500,1750,250)),axis=0),
    'mtllgmet':            frange(0, 500, 250) + frange(500, 1001, 500),
    # 'mt2':            np.arange(0.,425.,25.),
    'mt2':            frange(0., 425., 25.),
    # 'mzg_exo':            np.concatenate((np.arange(0,100,50),np.arange(100,200,50),np.arange(200,500,50)),axis=0),
    'mzg_exo':            frange(0, 100, 50) + frange(100, 200, 50) + frange(200, 500, 50),
    'gammaMotherID':            frange(0., 200., 1.),
    'genPhotonPT':            np.concatenate((np.arange(0, 100, 10), np.arange(100, 350, 50)), axis=0),
    'genPhotonPT_Veto':            np.concatenate((np.arange(0, 100, 10), np.arange(100, 350, 50)), axis=0),
    'PhotonPT_Veto':            np.concatenate((np.arange(0, 100, 10), np.arange(100, 350, 50)), axis=0),
    'genPhotonPT_NoVeto':            np.concatenate((np.arange(0, 100, 10), np.arange(100, 350, 50)), axis=0),
    'PhotonPT_NoVeto':            np.concatenate((np.arange(0, 100, 10), np.arange(100, 350, 50)), axis=0),
    'VetoCompare':            np.arange(0., 200., 1.),
    # 'DeltaPhiLLMet':            np.arange(0.,6.,0.3),
    'DeltaPhiLLMet':            frange(0., 3.2, 0.2),
    'DeltaPhiGMet':            frange(0., 6., 0.3),
    # 'DeltaEtaLLMet':            np.arange(0.,6.,0.3),
    'DeltaEtaLLMet':            frange(0., 6., 0.3),
    # 'DeltaRLLMet':            np.arange(0.,6.,0.3)
    'DeltaRLLMet':            frange(0., 6., 0.3),
    'nElectrons': frange(0, 10, 1),
    'nMuons': frange(0, 10, 1),
    'mTL3Met': frange(0., 200., 10),
    'Fakes': frange(0, 5, 1)
}
labels = {
    'pt1': ["p_{T}^{leading lepton} (GeV)", "Events / 30 GeV"],
    'pt2': ["p_{T}^{trailing lepton} (GeV)", " Events / 30 GeV"],
    'pt3': ["p_{T}^{2trailing lepton} (GeV)", " Events / 30 GeV"],
    'pt4': ["p_{T}^{3trailing lepton} (GeV)", " Events / 30 GeV"],
    'jetPt1': ["p_{T}^{leading jet} (GeV)", "Events / 30 GeV"],
    'jetPt2': ["p_{T}^{trailing jet} (GeV)", " Events / 30 GeV"],
    'jetPt3': ["p_{T}^{2trailing jet} (GeV)", " Events / 30 GeV"],
    'jetPt4': ["p_{T}^{3trailing jet} (GeV)", " Events / 30 GeV"],
    'eta1': ["|#eta_{leading lepton}|", "Events / 0.1"],
    'eta2': ["|#eta_{trailing lepton}|", "Events / 0.1"],
    'eta3': ["|#eta_{2trailing lepton}|", "Events / 0.1"],
    'eta4': ["|#eta_{3trailing lepton}|", "Events / 0.1"],
    'jetEta1': ["|#eta_{leading jet}|", "Events / 0.1"],
    'jetEta2': ["|#eta_{trailing jet}|", "Events / 0.1"],
    'jetEta3': ["|#eta_{2trailing jet}|", "Events / 0.1"],
    'jetEta4': ["|#eta_{3trailing jet}|", "Events / 0.1"],
    'phi1': ["|#phi_{leading lepton}|", "Events / 0.35"],
    'phi2': ["|#phi_{trailing lepton}|", "Events / 0.35"],
    'phi3': ["|#phi_{2trailing}|", "Events / 0.35"],
    'phi4': ["|#phi_{3trailing}|", "Events / 0.35"],
    'jetPhi1': ["|#Phi_{leading  jet}|", "Events / 0.35"],
    'jetPhi2': ["|#Phi_{trailing  jet}|", "Events / 0.35"],
    'jetPhi3': ["|#Phi_{2trailing  jet}|", "Events / 0.35"],
    'jetPhi4': ["|#Phi_{3trailing  jet}|", "Events / 0.35"],
    'ht': ["H_{T} (GeV)", "Events / 10 GeV"],
    'met': ["p_{T}^{miss} (GeV)", "Events / 10 GeV "],
    'gen_ht': ["H_{T}^{gen} (GeV)", "Events / 10 GeV"],
    # 'm_ll': ["m_{ll} from all leptons (GeV)" ,"Events / 65 GeV"],
    'm_ll': ["m_{ll} (GeV)", "Events / Bin"],
    'm_llll': ["m_{llll} (GeV)", "Events / Bin"],
    'm_ll2': ["m_{l_{3}l_{4}} (GeV)", "Events / Bin"],
    'n_jets': ["N Jets", "Events / 1"],
    'n_bjets': ["N BJets", "Events / 1"],
    'n_photons': ["N Photons", "Events / 1"],
    'n_vtx': ["N Vertices", "Events / 1"],
    'pt_g1': ["p_{T}^{#gamma} (GeV)", "Events / 10 GeV"],
    'eta_g1': ["|#eta|_{#gamma}", "Events / 0.1"],
    'phi_g1': ["|#phi|_{#gamma}", "Events / 0.1"],
    'sigmaIetaIeta_g1': ["#sigma_{I#etaI#eta}", "Events / 0.001"],
    'sigmaIphiIphi_g1': ["#sigma_{I#phiI#phi}", "Events / 0.001"],
    'r9_g1': ["r9", "Events / 0.001"],
    'hOverE_g1': ["H/E", "Events / 0.001"],
    'deltaR1_g1': ["#DeltaR1", "Events / 0.001"],
    'deltaR2_g1': ["#DeltaR2", "Events / 0.001"],
    'deltaEtaLL':       ["#Delta#Eta_{ll}", "Events / 0.1"],
    'deltaPhiLL':       ["#Delta#Phi_{ll}", "Events / 0.1"],
    'deltaEtaLL_neg':       ["#Delta#Eta_{ll}", "Events / 0.1"],
    'deltaPhiLL_neg':       ["#Delta#Phi_{ll}", "Events / 0.1"],
    'deltaEtaLLG':      ["#Delta#Eta_{ll,#gamma}", "Events / 0.1"],
    'deltaPhiLLG':      ["#Delta#Phi_{ll,#gamma}", "Events / 0.1"],
    'deltaRLL':         ["#DeltaR_{ll}", "Events / 0.1"],
    'deltaRLLG':        ["#DeltaR_{ll,#gamma}", "Events / 0.1"],
    'st':        ["S_{T} (GeV)", "Events / 100 GeV"],
    'stmet':        [" S_{T}+p_{T}^{miss} (GeV)", "Events / 100 GeV"],
    'zpt':        ["p_{T}^{Z} (GeV)", "Events / 50. GeV"],
    'mtll':        ["m_{T}^{ll} (GeV)", "Events / 100. GeV"],
    'mtllg':        ["m_{T}^{ll#gamma} (GeV)", "Events / 100. GeV"],
    'mtgmet':        ["m_{T}^{#gamma,met} (GeV)", "Events / 100. GeV"],
    'mtl1met':        ["m_{T}^{l1,met} (GeV)", "Events / 100. GeV"],
    'mtl2met':        ["m_{T}^{l2,met} (GeV)", "Events / 100. GeV"],
    'mtllmet':        ["m_{T}^{ll,met} (GeV)", "Events / 100. GeV"],
    'mt2':        ["M_{T2} (GeV)", "Events / 100. GeV"],
    'mtllgmet':        ["m_{T}^{ll#gamma,met} (GeV)", "Events / 100. GeV"],
    'm_{ll#gamma}':        ["m_{ll#gamma} (GeV)", "Events / 100. GeV"],
    'm_llg':        ["m_{ll#gamma} (GeV)", "Events / 100. GeV"],
    'pt_llg':        ["pt_llg", ""],
    'mzg_exo':        ["mzg_exo", ""],
    'gammaMotherID':        ["gammaMotherID", ""],
    'genPhotonPT':        ["gammaMotherID", ""],
    'genPhotonPT_Veto':        ["gammaMotherID", ""],
    'PhotonPT_Veto':        ["gammaMotherID", ""],
    'genPhotonPT_NoVeto':        ["gammaMotherID", ""],
    'PhotonPT_NoVeto':        ["gammaMotherID", ""],
    'VetoCompare':        ["gammaMotherID", ""],
    'DeltaPhiLLMet':        ["DeltaPhiLLMet", ""],
    'DeltaPhiGMet':        ["DeltaPhiGMet", ""],
    'DeltaEtaLLMet':        ["DeltaEtaLLMet", ""],
    'DeltaRLLMet':        ["DeltaRLLMet", ""],
    'nElectrons':           ["N_{Electrons}", ""],
    'nMuons':           ["N_{Muons}", ""],
    'mTL3Met':          ["m_{T}^{l3,met}", ""],
    'Fakes': ["", ""],
    'DeltaPhiL1G': ["#Delta#phi_{l1,#gamma}", ""],
    'DeltaPhiL2G': ["#Delta#phi_{l2,#gamma}", ""],
    'DeltaPhiLLG': ["#Delta#phi_{ll,#gamma}", ""],
    'm_l1g': ["#m_{l1,#gamma}", ""],
    'm_l2g': ["#m_{l2,#gamma}", ""]
}


def divideDatasetIntegrals(numerator, denominator, name):
    numMerged = sum(numerator)
    h_num = numMerged.getHist(name)
    num = h_num.Integral(0, -1)
    denMerged = sum(denominator)
    h_den = denMerged.getHist(name)
    den = h_den.Integral(0, -1)
    return num / den if den else 1.


def drawSameHistogram(sampleNames, name, bkg=[], additional=[], binning=None, binningName="", scaleToData=False, xTitle=None, yTitle=None):
    # def drawSameHistogram( sampleNames, name, bkg=[], additional=[], binning=None, binningName="", scaleToData=True, xTitle=None, yTitle=None ):
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()

    #style.divideByBinWidth = True
    style.divideByBinWidth = False

    style.minimumOne = False
    # style.minimumOne=True

    # yTitle=None

    scale = 1.
    if scaleToData:
        scale = divideDatasetIntegrals(
            [i for i in additional if "Data" in i.label], bkg, name)

    #print name

    folder = (name.split("/"))[0] + "/" + (name.split("/"))[1]
    #avgTopPtWeight = d.getHist(folder+"/weight_topPt").GetMean()
    #print avgTopPtWeight
    #print name,folder

    for d in bkg[-1::-1]:
        #h = d.getHist( name )
        h = d.getHistWithWeights(name, arWeightNames=["nISR", "topPt", "ewk"])

        if not h:
            continue
        if not h.Integral():
            continue
        h.Scale(scale)

        if binning:
            # if (binning.any()):
            h = aux.rebin(h, binning)

        aux.appendFlowBin(h)
        h.SetYTitle(aux.getYAxisTitle(h))
        if yTitle:
            h.SetYTitle(yTitle)
        else:
            h.SetYTitle(aux.getYAxisTitle(h))
        if xTitle:
            h.SetXTitle(xTitle)

        m.addStack(h, d.label)

    dataHist = None
    for d in additional:
        #h = d.getHist( name )
        h = d.getHistWithWeights(name, arWeightNames=["nISR", "topPt", "ewk"])
        if not h:
            continue
        if not h.Integral():
            continue
        # if (binning.any()):
        if (binning):
            h = aux.rebin(h, binning)
        aux.appendFlowBin(h)

        if h.GetLineColor() == ROOT.kBlack:  # data
            h.drawOption_ = "ep"
            h.SetMarkerStyle(20)
            # h.SetMarkerSize(0.5)
            h.SetMarkerSize(0.7)
            # h.SetLineWidth(0.7)
            # disable errors for data, so that ErrorOption is working
            # TODO: kPoisson also for rebinned and scaled histograms
            #if not(binning.any()): h.Sumw2(False)
            if not(binning):
                h.Sumw2(False)
            h.SetBinErrorOption(ROOT.TH1.kPoisson)
            dataHist = h
        else:
            h.drawOption_ = "hist e"
            h.SetLineWidth(3)

        m.add(h, d.label)

    m.sortStackByIntegral()

    if m.Draw():

        # ratio
        hsm = m.hists[0].GetStack().Last()
        if dataHist:
            r = ratio.Ratio("Data/MC", dataHist, hsm)
            r.draw(0.5, 1.5)

        info = ""
        # info = "ee" if "EE" in sampleNames else "#mu#mu" if "MM" in sampleNames else ""
        l = aux.Label(info="#scale[0.7]{%s}" % info, sim=not((dataDoubleMuon in additional)or(
            dataDoubleEG in additional)or(dataHt in additional)or(dataDoubleSF in additional)or(dataMuonEG in additional)))
        # l = aux.Label(info="#scale[0.7]{%s}"%info, sim=False)

        if binningName:
            binningName = "_" + binningName
        name = name.replace("/", "__")
        saveName = "sameHistograms_{}_{}{}".format(
            sampleNames, name, binningName)
        aux.save("DataMC_" + saveName)
        #aux.save("DataMC_"+saveName,folder="plots_noVeto/" )
        #aux.save("DataMC_"+saveName,folder="plots_uncorrected/" )
        #aux.save("DataMC_"+saveName,folder="plots_EGRegression/" )
        #aux.save("DataMC_"+saveName,folder="plots_LO/" )
        #aux.save("DataMC_"+saveName,folder="plots_mllWeights/" )
        #aux.save("DataMC_"+saveName,folder="plots_NLO/" )


def main():
    # bkgs=[DYjets,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma]
    # bkgs=[DYjetsLO,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma]
    # bkgs=[DYjetsNLO,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma]
    bkgs = [DYjetsNLO, zgamma, tt, ttgamma, wwgamma,
            wzgamma, zz, wjets, wgamma, singletop, zz4l, ww, wz]
    # bkgs=[DYjetsNLO,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma]
    # bkgs=[DYjetsLO,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma]
    # bkgs=[DYjets,zgamma]
    # bkgs=[DYjets]
    # variables=["eta1","pt1","n_jets","n_vtx","phi1","m_ll","ht","n_photons","pt_g1","met","nElectrons","nMuons"]
    # variables = ["eta1", "pt1", "n_jets", "n_vtx", "phi1", "m_ll", "ht", "n_photons",
    #              "pt_g1", "nElectrons", "nMuons", "deltaR1_g1", "deltaR2_g1", "deltaRLL"]
    # variables=["eta1","pt1","n_jets","n_vtx","phi1","m_ll","ht","n_photons","pt_g1"]
    # variables=["pt1","n_photons","met"]
    variables = ["eta_g1"]
    # variables=["eta1","pt1","n_jets","n_vtx","phi1","m_ll","ht","n_photons","pt_g1","sigmaIetaIeta_g1","sigmaIphiIphi_g1","r9_g1","hOverE_g1","deltaR1_g1","deltaR2_g1","eta_g1","phi_g1"]
    # groups=["sel","onZ","dilep","onZG","mllG110"]
    # groups=["sel","onZ","dilep","onZG","exo"]
    # groups = ["sel", "dilep", "onZ"]
    groups = ["onZ"]
    # groups=["sel"]
    # groups=["sel","onZ"]
    # groups=["sel","dilep"]
    # groups=["EGRegression"]
    for group in groups:
        for variable in variables:
            #drawSameHistogram("EE",group+"EE/"+variable, bkgs, additional=[dataDoubleEG],binning=binnings[variable],yTitle=labels[variable][1],xTitle=labels[variable][0])
            #drawSameHistogram("EE",group+"EE/"+variable, bkgs, additional=[dataDoubleEG],binning=binnings[variable],xTitle=labels[variable][0])
            drawSameHistogram("MM", group + "/MM/" + variable, bkgs, additional=[
                              dataDoubleMuon], binning=binnings[variable], xTitle=labels[variable][0])
            #drawSameHistogram("EM",group+"EM/"+variable, bkgs, additional=[dataMuonEG],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("LL",group+"/"+variable, bkgs, additional=[dataDoubleSF],binning=binnings[variable],xTitle=labels[variable][0])
            drawSameHistogram("EE", group + "/EE/" + variable, bkgs, additional=[
                              dataLL], binning=binnings[variable], xTitle=labels[variable][0])
            #drawSameHistogram("MM",group+"MM/"+variable, bkgs, additional=[dataLL],binning=binnings[variable],xTitle=labels[variable][0])
            drawSameHistogram("EM", group + "/EM/" + variable, bkgs, additional=[
                              dataLL], binning=binnings[variable], xTitle=labels[variable][0])
            drawSameHistogram("LL", group + "/LL/" + variable, bkgs, additional=[
                              dataLL], binning=binnings[variable], xTitle=labels[variable][0])
            #drawSameHistogram("MM",group+"MM/"+variable, bkgs, additional=[dataDoubleMuon],binning=binnings[variable],yTitle=labels[variable][1],xTitle=labels[variable][0])
            #drawSameHistogram("EE",group+"EE/"+variable, bkgs, additional=[zgamma],binning=binnings[variable],yTitle=labels[variable][1],xTitle=labels[variable][0])
            #drawSameHistogram("MM",group+"MM/"+variable, bkgs, additional=[zgamma],binning=binnings[variable],yTitle=labels[variable][1],xTitle=labels[variable][0])
            #drawSameHistogram("EE",group+"EE/"+variable, bkgs, additional=[zgamma],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("MM",group+"MM/"+variable, bkgs, additional=[zgamma],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("EE",group+"EE/"+variable, bkgs, additional=[DYjets,zgamma],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("MM",group+"MM/"+variable, bkgs, additional=[DYjets,zgamma],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("EE",group+"EE/"+variable, bkgs ,binning=binnings[variable],yTitle=labels[variable][1],xTitle=labels[variable][0])
            #drawSameHistogram("MM",group+"MM/"+variable, bkgs,binning=binnings[variable],yTitle=labels[variable][1],xTitle=labels[variable][0])


def main3():
    # bkgs=[DYjetsNLO,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma]
    bkgs = [DYjetsNLO, zgamma, tt, ttgamma, wwgamma,
            wzgamma, zz, wjets, wgamma, singletop, zz4l, wz, ww]
    variables = ["eta1", "pt1", "n_jets", "n_vtx", "phi1", "m_ll", "ht",
                 "n_photons", "pt_g1", "met", "deltaRLL", "deltaR1_g1", "deltaR2_g1"]
    groups = ["onZMet0100"]
    for group in groups:
        for variable in variables:
            drawSameHistogram("EE", group + "/EE/" + variable, bkgs, additional=[
                              dataDoubleEG], binning=binnings[variable], xTitle=labels[variable][0])
            drawSameHistogram("MM", group + "/MM/" + variable, bkgs, additional=[
                              dataDoubleMuon], binning=binnings[variable], xTitle=labels[variable][0])
            drawSameHistogram("LL", group + "/LL/" + variable, bkgs, additional=[
                              dataDoubleSF], binning=binnings[variable], xTitle=labels[variable][0])


def main2():
    # bkgs=[DYjetsNLO,zgamma,wwgamma,wzgamma,ttgamma,zz,tt,wjets]
    # bkgs=[DYjetsNLO,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma]
    # bkgs=[DYjetsNLO,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma,singletop,zz4l,wz,ww]
    bkgs = [DYjetsNLO, zgamma, tt, ttgamma, wwgamma,
            wzgamma, wjets, wgamma, singletop, zz4l, wz, ww]
    # bkgs=[zgamma]
    # variables=["mt2","Fakes","met","stmet","st","mtll","mtllg","mtl1met","mtl2met","mtllmet","mtllgmet","mtgmet","m_llg","m_ll","n_vtx","DeltaPhiLLMet","DeltaEtaLLMet","DeltaRLLMet","pt1","pt_g1","ht","eta1","phi1","n_jets","n_photons","deltaPhiLLG","deltaEtaLLG","deltaRLLG","deltaPhiLL","deltaEtaLL","deltaRLL","zpt","nElectrons","nMuons","deltaR1_g1","deltaR2_g1","DeltaPhiGMet","n_bjets"]
    # variables=["met","pt1","pt_g1","ht","eta1","phi1","n_jets","n_photons"]
    # variables = ["Fakes"]
    # variables = ["eta_g1"]
    variables = ["mt2"]
    # variables=["met"]
    # variables=["hOverE_g1"]
    # groups=["sel","onZ","onZMet","onZMet0100","onZMet100","onZMet200","onZMet100200","onZMet100300","onZMet200300"]
    # groups=["sel","onZ","onZMet0100","onZMet100","onZMet150","onZMet100150"]
    groups = ["onZ"]
    # groups = ["onZMet150"]
    # groups=["onZMet0100"]
    # groups=["sel"]
    # groups=["onZ"]
    # groups=["onZMet100200",]
    binnings_ = binnings.copy()
    binnings_['pt_g1'] = frange(0., 100., 50.) + frange(100., 251., 150.)
    binnings_['met'] = [0, 25, 50, 75, 100, 150, 250, 450]
    binnings_['m_llg'] = [0, 100, 250, 500]
    binnings_['deltaPhiLL'] = [0, 1, 2, 3.5]
    binnings_["mt2"] = frange(80., 300., 10.)
    # binnings_["mt2"] = [0., 100, 400]
    # binnings_['m_llg']=[0,150,250,500]
    for group in groups:
        for variable in variables:
            #drawSameHistogram("EE+signal",group+"EE/"+variable, bkgs, additional=[t5bbbbzg_1800_1700,t5bbbbzg_1800_400,t5bbbbzg_1800_600,tching_1200,tching_400,t5ttttzg_1800_400,t5ttttzg_1800_600],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("MM+signal",group+"MM/"+variable, bkgs, additional=[t5bbbbzg_1800_1700,t5bbbbzg_1800_400,t5bbbbzg_1800_600,tching_1200,tching_400,t5ttttzg_1800_400,t5ttttzg_1800_600],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("EE+signal",group+"EE/"+variable, bkgs, additional=[t5bbbbzg_1500_1400,t5bbbbzg_1500_400,t5bbbbzg_1500_600,tching_600,tching_400],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("MM+signal",group+"MM/"+variable, bkgs, additional=[t5bbbbzg_1500_1400,t5bbbbzg_1500_400,t5bbbbzg_1500_600,tching_600,tching_400],binning=binnings[variable],xTitle=labels[variable][0])
            drawSameHistogram("LL+signal", group + "/LL/" + variable, bkgs, additional=[
                              t5bbbbzg_1500_1400, t5bbbbzg_1500_400, t5bbbbzg_1500_600, tching_600, tching_400, gmsb_240_230, gmsb_290_205], binning=binnings_[variable], xTitle=labels[variable][0])
            #drawSameHistogram("EM+signal",group+"EM/"+variable, bkgs, additional=[t5bbbbzg_1500_1400,t5bbbbzg_1500_400,t5bbbbzg_1500_600,tching_600,tching_400],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("EE+signal",group+"EE/"+variable, bkgs, additional=[t5bbbbzg_1800_1700,t5bbbbzg_1800_400,t5bbbbzg_1800_600,tching_1200,tching_400,ggm_m1550_m2750,ggm_m11200_m21000,ggm_m1800_m2600],binning=binnings[variable],xTitle=labels[variable][0])
            #drawSameHistogram("MM+signal",group+"MM/"+variable, bkgs, additional=[t5bbbbzg_1800_1700,t5bbbbzg_1800_400,t5bbbbzg_1800_600,tching_1200,tching_400,ggm_m1550_m2750,ggm_m11200_m21000,ggm_m1800_m2600],binning=binnings[variable],xTitle=labels[variable][0])


if __name__ == "__main__":
    # main()
    main2()
    # main3()
