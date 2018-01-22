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
    #'pt2':              np.arange(0., 300, 10),
    'pt2':              np.arange(0., 100, 1),
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



def efficiency(dataset, name, savename="", binning=None, binningName=""):
    c = ROOT.TCanvas()
    c.SetLogy(0)
    eff = dataset.getHist(name)
    if eff.UsesWeights(): eff.SetStatisticOption(ROOT.TEfficiency.kFNormal)

    h_pas = eff.GetPassedHistogram()
    h_tot = eff.GetTotalHistogram()

    if name.endswith("_ps"):
        eff2 = dataset.getHist(name.replace("_ps", ""))
        h_tot = eff2.GetTotalHistogram()

    if (binning.any()):
        h_pas = aux.rebin(h_pas, binning, False)
        h_tot = aux.rebin(h_tot, binning, False)

    if name.endswith("_ps"):
        ratio = h_pas.Clone(aux.randomName())
        ratio.Divide(h_tot)
        gr = ROOT.TGraphAsymmErrors(ratio)
        gr.SetTitle(";{};prescaled #varepsilon".format(eff.CreateGraph().GetHistogram().GetXaxis().GetTitle()))
        gr.SetLineColor(1)
        gr.Draw("ap")
    else:
        x = h_pas.Clone(aux.randomName())
        y = h_tot.Clone(aux.randomName())
        eff = ROOT.TEfficiency(h_pas, h_tot)
        h_pas, h_tot = x, y
        eff.Draw()
        ROOT.gPad.Update()
        gr = eff.GetPaintedGraph()

    gr.GetYaxis().SetRangeUser(0., 1.1)
    if name.endswith("emht__ht600__p90"):
        gr.GetYaxis().SetRangeUser(0., 0.1)

    if "eff_pt__p90ht600__ht600" in name or "eff_pt_ee__p90ht600__ht600" in name:
        cutValue = 100
    elif "emht__ht600" in name or "emht__p90ht600" in name:
        cutValue = 700
    elif "emht__ht800" in name:
        cutValue = 900
    elif "eff_pt__ele27__ht600" in name:
        cutValue = 30
    elif "pt1" in name:
        cutValue = 25
    elif "pt2" in name:
        cutValue = 20
    else:
        cutValue = 0

    if cutValue or True:
        bin = h_pas.FindFixBin(cutValue)
        passed = int(h_pas.Integral(bin, -1))
        total = int(h_tot.Integral(bin, -1))
        if not total: return
        e = 1.*passed/total
        if passed<=total:
            conf = ROOT.TEfficiency().GetConfidenceLevel()
            e_up = ROOT.TEfficiency.ClopperPearson(total, passed, conf, True)
            e_dn = ROOT.TEfficiency.ClopperPearson(total, passed, conf, False)
        else:
            passed, epassed = aux.integralAndError(h_pas,bin,-1)
            total, etotal = aux.integralAndError(h_pas,bin,-1)
            ee = e * math.sqrt((epassed/passed)**2 + (etotal/total)**2)
            e_up = e + ee
            e_dn = e - ee
        eLabel = ROOT.TLatex(0.7, .15, "#varepsilon = {:.1f}^{{#plus{:.1f}}}_{{#minus{:.1f}}}%".format(100*e, 100*(e_up-e),100*(e-e_dn)))
        eLabel.SetNDC()
        eLabel.Draw()

        # graphical representation
        l = ROOT.TLine()
        l.SetLineWidth(2)
        l.SetLineColor(ROOT.kRed)
        xmin = gr.GetHistogram().GetXaxis().GetXmin()
        xmax = gr.GetHistogram().GetXaxis().GetXmax()
        l.DrawLine(max(cutValue,xmin), e, xmax, e)
        l.DrawLine(max(cutValue,xmin), e_up, xmax, e_up)
        l.DrawLine(max(cutValue,xmin), e_dn, xmax, e_dn)

        if cutValue > xmin:
            # cut line
            l.SetLineStyle(2)
            ymin = gr.GetYaxis().GetXmin()
            ymax = gr.GetYaxis().GetXmax()
            l.DrawLine(cutValue, ymin, cutValue, ymax)

    saveLumi = aux.intLumi
    if "B-F" in dataset.label:
        aux.intLumi = 20.101e3 # brilcalc lumi -b "STABLE BEAMS" --normtag=/afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt --end 278808 # Up to run F
    l = aux.Label(sim="Data" not in dataset.label)
    aux.intLumi = saveLumi
    name = "_"+name.split("/")[-1]
    if binningName: binningName = "_"+binningName
    aux.save("efficiency_"+savename+name+binningName, log=False)

    if False:
        h_tot.SetLineColor(2)
        h_tot.Draw("hist")
        h_pas.Draw("same e*")
        aux.save("efficiency_"+savename+name+"_raw")



#def efficiencies(dataset, savename=""):
    #names = ["triggerStudies/"+x for x in aux.getObjectNames(dataset.files[0], "triggerStudies", [ROOT.TEfficiency]) ]
    #for name in names:
        #efficiency(dataset, name, savename)
        #if "_pt_" in name and "ele" not in name:
            #efficiency(dataset, name, savename, range(0,80,8) + range(80,108,4) + range(108,300,12) + [300,400,500, 1000], "1")
        #if "_emht__" in name:
            #efficiency(dataset, name, savename, range(500,1001,40) + range(1000,1500,2000), "1")
        #if "_met__" in name:
            #efficiency(dataset, name, savename, range(0, 100, 5)+range(100,201, 10), "1")
        #if "_nIso__" in name:
            #efficiency(dataset, name, savename, [0, .2, .4, .8, 1, 1.2, 1.4, 1.6, 2, 2.2, 2.4, 2.6, 2.8, 3, 4, 5, 6, 7, 8, 9, 10, 15], "1")
        #if "_r9__" in name:
            #efficiency(dataset, name, savename, [.4,.5,.6,.7,.8]+aux.frange(0.8, 1.1, 0.01), "1")
        #if "_sie__" in name:
            #efficiency(dataset, name, savename, [.0045,.007,.0075]+aux.frange(0.0075,0.011,5e-5), "1")

#def main():
    #efficiencies(dataHt, "jetHt")





def main():
    bkgs=[DYjets,zgamma,wwgamma,wzgamma,ttgamma,zz,tt,wjets]
    #variables=["eta1","pt1","n_jets","n_vtx","phi1","m_ll","ht","n_photons"]
    variables=["pt1","pt2"]
    #groups=["trigDilep","trigSel","trigOnZ"]
    groups=["trigDilep"]
    #groups=["trigDilep_ptcuts","trigSel_ptcuts","trigOnZ_ptcuts"]
    for group in groups:
        for variable in variables:
            efficiency(dataHt,group+"EE/"+variable,"dataHT_"+group+"_EE",binning=binnings[variable])
            efficiency(dataHt,group+"MM/"+variable,"dataHT_"+group+"_MM",binning=binnings[variable])
            efficiency(dataHt,group+"MM_ptcuts/"+variable,"dataHT_"+group+"_ptcuts"+"_MM",binning=binnings[variable])
            efficiency(dataHt,group+"EE_ptcuts/"+variable,"dataHT_"+group+"_ptcuts_"+"_EE",binning=binnings[variable])
            efficiency(dataHt,group+"MM_pt1cut/"+variable,"dataHT_"+group+"_pt1cut"+"_MM",binning=binnings[variable])
            efficiency(dataHt,group+"EE_pt1cut/"+variable,"dataHT_"+group+"_pt1cut_"+"_EE",binning=binnings[variable])
            efficiency(dataHt,group+"MM_pt2cut/"+variable,"dataHT_"+group+"_pt2cut"+"_MM",binning=binnings[variable])
            efficiency(dataHt,group+"EE_pt2cut/"+variable,"dataHT_"+group+"_pt2cut_"+"_EE",binning=binnings[variable])
  

main()
