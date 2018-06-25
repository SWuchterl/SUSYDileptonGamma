import ROOT
from ROOT import *
from array import array
from include import *
import numpy as np

def frange(start, end, step):
    a=[]
    tmp = start
    while(tmp < end):
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
    #'pt1':              np.arange(0., 300., 10.),
    'pt1':              frange(0., 310., 10.),
    #'pt2':              np.arange(0., 300, 10),
    #'pt2':              np.arange(0., 100, 1),
    'pt2':              frange(0., 102, 2.),
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
    'pt1': ["p_{T}^{leading #mathcal{l}}(GeV)","Events / 30 GeV"],
    'pt2': ["p_{T}^{trailing #mathcal{l}}(GeV)"," Events / 30 GeV"],
    'eta1': ["|#Eta_{leading}|","Events / 0.1"],
    'eta2': ["|#Eta_{trailing}|","Events / 0.1"],
    'phi1': ["|#Phi_{leading}|","Events / 0.1"],
    'phi2': ["|#Phi_{trailing}|","Events / 0.1"],
    'ht': ["H_{T} (GeV)", "Events / 10 GeV"],
    'met': ["E_{T}^{miss} (GeV)", "Events / 10 GeV "],
    'gen_ht': ["H_{T}^{gen} (GeV)", "Events / 10 GeV"],
    'm_ll': ["m_{ll} from all leptons (GeV)" ,"Events / 65 GeV"],
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
    #'st':        ["#itS_T (GeV)","Events / 100 GeV"],    
    #'stmet':        ["#it S_T+E_{T}^{miss}" (GeV),"Events / 100 GeV"],    
    #'zpt':        ["Z_{p_T} (GeV)","Events / 50. GeV"],    
}



def efficiency(dataset, name, savename="", binning=None, binningName="",additional=[]):

    c = ROOT.TCanvas()
    c.SetLogy(0)
    print name
    eff = dataset.getHist(name)
    if eff.UsesWeights(): eff.SetStatisticOption(ROOT.TEfficiency.kFNormal)
    if additional:
        effAdd = additional[0].getHist(name)
        if effAdd.UsesWeights(): effAdd.SetStatisticOption(ROOT.TEfficiency.kFNormal)

    h_pas = eff.GetPassedHistogram()
    h_tot = eff.GetTotalHistogram()
    if additional:
        h_pasAdd = effAdd.GetPassedHistogram()
        h_totAdd = effAdd.GetTotalHistogram()

    if name.endswith("_ps"):
        eff2 = dataset.getHist(name.replace("_ps", ""))
        h_tot = eff2.GetTotalHistogram()
        if additional:
            eff2Add = additional[0].getHist(name.replace("_ps", ""))
            h_totAdd = eff2Add.GetTotalHistogram()
    
    #if (binning.any()):
    if (binning):
        h_pas = aux.rebin(h_pas, binning, False)
        h_tot = aux.rebin(h_tot, binning, False)

        if additional:
            h_pasAdd = aux.rebin(h_pasAdd, binning, False)
            h_totAdd = aux.rebin(h_totAdd, binning, False)



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
        
        titleForNow = gr.GetXaxis().GetTitle()
        if "leading" in titleForNow:
            titleNow = titleForNow.replace("leading","leading lepton")
        if "trailing" in titleForNow:
            titleNow = titleForNow.replace("trailing","trailing lepton")
        gr.GetXaxis().SetTitle(titleNow)
        eff.SetTitle(";"+titleNow+";"+gr.GetYaxis().GetTitle()+"#varepsilon")
        
        
        if additional:
            xAdd = h_pasAdd.Clone(aux.randomName())
            yAdd = h_totAdd.Clone(aux.randomName())
            effAdd = ROOT.TEfficiency(h_pasAdd, h_totAdd)
            h_pasAdd, h_totAdd = xAdd, yAdd
            effAdd.Draw("same")
            effAdd.SetLineColor(kRed)
            ROOT.gPad.Update()
            grAdd = effAdd.GetPaintedGraph()
            grAdd.GetXaxis().SetTitle(titleNow)
            grAdd.SetTitle(";"+gr.GetXaxis().GetTitle()+";"+titleNow)

    #gr.GetYaxis().SetRangeUser(0., 1.1)
    gr.GetYaxis().SetRangeUser(0.5, 1.1)
    

    if additional:
        #grAdd.GetYaxis().SetRangeUser(0., 1.1)
        grAdd.GetYaxis().SetRangeUser(0.5, 1.1)
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
        #print h_tot.Integral()
        print cutValue
        bin = h_pas.FindFixBin(cutValue)
        print bin
        passed = int(h_pas.Integral(bin, -1))
        total = int(h_tot.Integral(bin, -1))
        print total,h_tot.Integral(bin,-1)
        if not total: return
        e = 1.*passed/total
        if additional:
            #print "huhu"
            binAdd = h_pasAdd.FindFixBin(cutValue)
            passedAdd = int(h_pasAdd.Integral(binAdd, -1))
            totalAdd = int(h_totAdd.Integral(binAdd, -1))
            #print totalAdd
            if not totalAdd: return
            eAdd = 1.*passedAdd/totalAdd
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
        e_upSyst = e+0.03
        e_dnSyst = e-0.03
        
        e_syst=0.03
        if e_upSyst>1.:
            e_upSyst=1.
        #eLabel = ROOT.TLatex(0.57, .17, "#varepsilon_{{Data}} = {:.2f}^{{#plus{:.2f}}}_{{#minus{:.2f}}}%".format(100*e, 100*(e_up-e),100*(e-e_dn)))
        eLabel = ROOT.TLatex(0.47, .17, "#varepsilon_{{Data}} = {:.1f}^{{#plus{:.1f}}}_{{#minus{:.1f}}}%".format(100*e, 100*(e_up-e),100*(e-e_dn)))
        #eLabel = ROOT.TLatex(0.47, .17, "#varepsilon_{{Data}} = {:.1f}^{{#plus{:.1f}}}_{{#minus{:.1f}}}(stat.)^{{#plus{:.1f}}}_{{#minus{:.1f}}}(syst.) %".format(100*e, 100*(e_up-e),100*(e-e_dn),100*(e_upSyst-e),100*(e-e_dnSyst)))
        eLabel.SetNDC()
        eLabel.Draw()
        if additional:
            if passedAdd<=totalAdd:
                confAdd = ROOT.TEfficiency().GetConfidenceLevel()
                e_upAdd = ROOT.TEfficiency.ClopperPearson(totalAdd, passedAdd, confAdd, True)
                e_dnAdd = ROOT.TEfficiency.ClopperPearson(totalAdd, passedAdd, confAdd, False)
            else:
                passedAdd, epassedAdd = aux.integralAndError(h_pasAdd,binAdd,-1)
                totalAdd, etotalAdd = aux.integralAndError(h_pasAdd,binAdd,-1)
                eeAdd = eAdd * math.sqrt((epassedAdd/passedAdd)**2 + (etotalAdd/totalAdd)**2)
                e_upAdd = eAdd + eeAdd
                e_dnAdd = eAdd - eeAdd
            #eLabelAdd = ROOT.TLatex(0.57, .25, "#varepsilon_{{MC}} = {:.2f}^{{#plus{:.2f}}}_{{#minus{:.2f}}}%".format(100*eAdd, 100*(e_upAdd-eAdd),100*(eAdd-e_dnAdd)))
            eLabelAdd = ROOT.TLatex(0.47, .25, "#varepsilon_{{MC}} = {:.2f}^{{#plus{:.2f}}}_{{#minus{:.2f}}}%".format(100*eAdd, 100*(e_upAdd-eAdd),100*(eAdd-e_dnAdd)))
            eLabelAdd.SetTextColor(kRed)
            eLabelAdd.SetNDC()
            eLabelAdd.Draw()

        # graphical representation
        l = ROOT.TLine()
        l.SetLineWidth(1)
        l.SetLineColor(ROOT.kGray+2)
        xmin = gr.GetHistogram().GetXaxis().GetXmin()
        xmax = gr.GetHistogram().GetXaxis().GetXmax()
        l.DrawLine(max(cutValue,xmin), e, xmax, e)
        l.DrawLine(max(cutValue,xmin), e_up, xmax, e_up)
        l.DrawLine(max(cutValue,xmin), e_dn, xmax, e_dn)

        

        nBins=len(binning)
        binningMin=binning[0]
        binningMax=binning[-1]
        
        #print nBins,binningMax,binningMin

        #TH1 h; // the histogram (you should set the number of bins, the title etc)
        tempH = TH1F("","",nBins-1,binningMin,binningMax)
        nPoints = grAdd.GetN()
        #for(int i=0; i < nPoints; ++i) {
           #double x,y;
           #graph.GetPoint(i, x, y);
           #h->Fill(x,y); // ?
        #}
        for i in range(nPoints):
           x=ROOT.Double(0)
           y=ROOT.Double(0)
           grAdd.GetPoint(i, x, y)
           #print x
           #print x,y
           tempH.Fill(x,y) 


        #print tempH.GetBinContent(tempH.FindBin(50))

        #effAddSyst = aux.getSysHisto(effAdd, 0.03)
        effAddSyst = aux.getSysHisto(tempH, 0.03)
        #effAddSyst = aux.getSysHisto(grAdd, 0.03)
        aux.drawOpt(effAddSyst, "sysUnc")
        
        #x= array.array("f",[xmin, xmax]) 
        #y= array.array("f", [e,e]) 
        #ex= array.array("f", [0.,0.])
        #ey= array.array("f", [e*0.03,e*0.03])
        #ge = ROOT.TGraphErrors(2, x, y, ex, ey)
        #ge.SetFillColor(ROOT.kRed-3)
        #ge.SetFillStyle(1001)
        #ge.SetFillStyle(3004)
        #ge.SetLineColor(ROOT.kWhite)
        #ge.Draw("SAME 3")
        effAddSyst.Draw("same e2")

        #DRAW EVERYTHING AGAIN
        #if name.endswith("_ps"):
            #gr.Draw("ap")
        #else:
        #eff.Draw()
        #
        #if additional:
            #effAdd.Draw("same")
            #
        #eLabel.Draw()
        #if additional:
            #eLabelAdd.Draw()
#
        #l.DrawLine(max(cutValue,xmin), e, xmax, e)
        #l.DrawLine(max(cutValue,xmin), e_up, xmax, e_up)
        #l.DrawLine(max(cutValue,xmin), e_dn, xmax, e_dn)

        leg = ROOT.TLegend(.56,.49,.94,.615)
        leg.AddEntry( eff, "Data (ee)" if "EE" in savename else "Data (#mu#mu)" if "MM" in savename else "Data (e#mu)", "epl" )
        if additional:
            leg.AddEntry( effAdd, "MC (ee)" if "EE" in savename else "MC (#mu#mu)" if "MM" in savename else "MC (e#mu)", "epl" )
        linie = ROOT.TLine()
        linie.SetLineWidth(1)
        linie.SetLineColor(ROOT.kGray+2)
        leg.AddEntry(linie,"mean_{Data} #pm 1 #sigma stat.","l")
        systLinie= ROOT.TLine()
        systLinie.SetLineWidth(5)
        systLinie.SetLineColor(ROOT.kRed-3)
        #systLinie.SetFillStyle(3004)
        #leg.AddEntry(ge,"mean #pm 1 #sigma syst.","")
        #leg.AddEntry(ge,"syst. unc.","f")
        leg.AddEntry(effAddSyst,"syst. unc.","f")
        leg.Draw()
        
        leg2 = ROOT.TLegend(.56,.501,.94,.52)
        leg2.SetFillStyle(0)
        linie.SetLineWidth(1)
        leg2.AddEntry(linie," ","l")
        leg2.AddEntry(linie," ","l")
        #leg2.Draw("same")


        if cutValue > xmin:
            # cut line
            l.SetLineStyle(2)
            ymin = gr.GetYaxis().GetXmin()
            ymax = gr.GetYaxis().GetXmax()
            #l.DrawLine(cutValue, ymin, cutValue, ymax)
            l.DrawLine(cutValue, 0.5, cutValue, ymax)




    saveLumi = aux.intLumi
    if "B-F" in dataset.label:
        aux.intLumi = 20.101e3 # brilcalc lumi -b "STABLE BEAMS" --normtag=/afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt --end 278808 # Up to run F
    l = aux.Label(sim="Data" not in dataset.label)
    aux.intLumi = saveLumi
    #textPW="Work in Progress"
    #pw = ROOT.TLatex( 0.2, .887, "#scale[0.76]{#font[52]{%s}}"%textPW )
    #pw.Draw()
    name = "_"+name.split("/")[-1]
    if binningName: binningName = "_"+binningName
    #aux.save("efficiency_"+savename+name+binningName, log=False)
    aux.save("efficiency_"+savename+name+binningName, log=False,folder="triggerStudies/")

    if False:
        h_tot.SetLineColor(2)
        h_tot.Draw("hist")
        h_pas.Draw("same e*")
        #aux.save("efficiency_"+savename+name+"_raw")
        aux.save("efficiency_"+savename+name+"_raw",folder="triggerStudies/")




def efficiencyAllMC(dataset, name, savename="", binning=None, binningName="",additional=[]):
    c = ROOT.TCanvas()
    c.SetLogy(0)
    eff = dataset.getHist(name)
    if eff.UsesWeights(): eff.SetStatisticOption(ROOT.TEfficiency.kFNormal)
    
    #zgamma+ttgamma+zz+wwgamma+wzgamma+DYjetsNLO+wjets+tt+zz4l+wz+ww+singletop+wgamma
    #if additional:
            #effAdd = additional[0].getHist(name)
            #if effAdd.UsesWeights(): effAdd.SetStatisticOption(ROOT.TEfficiency.kFNormal)
            
            
    effAdd_zgamma = zgamma.getHist(name)
    effAdd_ttgamma = ttgamma.getHist(name)
    effAdd_zz = zz.getHist(name)
    effAdd_wwgamma = wwgamma.getHist(name)
    effAdd_wzgamma = wzgamma.getHist(name)
    effAdd_DYjetsNLO = DYjetsNLO.getHist(name)
    effAdd_wjets = wjets.getHist(name)
    effAdd_tt = tt.getHist(name)
    effAdd_zz4l = zz4l.getHist(name)
    effAdd_wz = wz.getHist(name)
    effAdd_ww = ww.getHist(name)
    effAdd_singletop = singletop.getHist(name)
    effAdd_wgamma = wgamma.getHist(name)

    arAdd=[effAdd_zgamma,effAdd_ttgamma,effAdd_zz,effAdd_wwgamma,effAdd_wzgamma,effAdd_DYjetsNLO,effAdd_wjets,effAdd_tt,effAdd_zz4l,effAdd_wz,effAdd_ww,effAdd_singletop,effAdd_wgamma]

    for i in range(len(arAdd)):
        if arAdd[i].UsesWeights(): arAdd[i].SetStatisticOption(ROOT.TEfficiency.kFNormal)

    h_pas = eff.GetPassedHistogram()
    h_tot = eff.GetTotalHistogram()
    #if additional:
        #for iAdd in effAddArray:
            #effAddArrayPas.append(iAdd.GetPassedHistogram())
            #effAddArrayTot.append(iAdd.GetTotalHistogram())

    zgammaPas = effAdd_zgamma.GetPassedHistogram()
    zgammaTot = effAdd_zgamma.GetTotalHistogram()
    ttgammaPas = effAdd_ttgamma.GetPassedHistogram()
    ttgammaTot = effAdd_ttgamma.GetTotalHistogram()
    zzPas = effAdd_zz.GetPassedHistogram()
    zzTot = effAdd_zz.GetTotalHistogram()
    wwgammaPas = effAdd_wwgamma.GetPassedHistogram()
    wwgammaTot = effAdd_wwgamma.GetTotalHistogram()
    wzgammaPas = effAdd_wzgamma.GetPassedHistogram()
    wzgammaTot = effAdd_wzgamma.GetTotalHistogram()
    DYjetsNLOPas = effAdd_DYjetsNLO.GetPassedHistogram()
    DYjetsNLOTot = effAdd_DYjetsNLO.GetTotalHistogram()
    wjetsPas = effAdd_wjets.GetPassedHistogram()
    wjetsTot = effAdd_wjets.GetTotalHistogram()
    ttPas = effAdd_tt.GetPassedHistogram()
    ttTot = effAdd_tt.GetTotalHistogram()
    zz4lPas = effAdd_zz4l.GetPassedHistogram()
    zz4lTot = effAdd_zz4l.GetTotalHistogram()
    wzPas = effAdd_wz.GetPassedHistogram()
    wzTot = effAdd_wz.GetTotalHistogram()
    wwPas = effAdd_ww.GetPassedHistogram()
    wwTot = effAdd_ww.GetTotalHistogram()
    singletopPas = effAdd_singletop.GetPassedHistogram()
    singletopTot = effAdd_singletop.GetTotalHistogram()
    wgammaPas = effAdd_wgamma.GetPassedHistogram()
    wgammaTot = effAdd_wgamma.GetTotalHistogram()

    arAddPas=[zgammaPas,ttgammaPas,zzPas,wwgammaPas,wzgammaPas,DYjetsNLOPas,wjetsPas,ttPas,zz4lPas,wzPas,wwPas,singletopPas,wgammaPas]
    arAddTot=[zgammaTot,ttgammaTot,zzTot,wwgammaTot,wzgammaTot,DYjetsNLOTot,wjetsTot,ttTot,zz4lTot,wzTot,wwTot,singletopTot,wgammaTot]

    #if name.endswith("_ps"):
        #eff2 = dataset.getHist(name.replace("_ps", ""))
        #h_tot = eff2.GetTotalHistogram()
        #if additional:
            #for iAdd in
            #eff2Add = additional[0].getHist(name.replace("_ps", ""))
            #h_totAdd = eff2Add.GetTotalHistogram()

    if (binning):
        h_pas = aux.rebin(h_pas, binning, False)
        h_tot = aux.rebin(h_tot, binning, False)
        #if additional:
            #for iAdd in effAddArrayPas:
                #iAdd = aux.rebin(iAdd, binning, False)
            #for iAdd in effAddArrayTot:
                #iAdd = aux.rebin(iAdd, binning, False)
        for i in range(len(arAddPas)):
            arAddPas[i]=aux.rebin(arAddPas[i],binning,False)
            arAddTot[i]=aux.rebin(arAddTot[i],binning,False)
        
        
        zgammaPas = aux.rebin(zgammaPas,binning,False)
        zgammaTot = aux.rebin(zgammaTot,binning,False)
        ttgammaPas = aux.rebin(ttgammaPas,binning,False)
        ttgammaTot = aux.rebin(ttgammaTot,binning,False)
        zzPas = aux.rebin(zzPas,binning,False)
        zzTot = aux.rebin(zzTot,binning,False)
        wwgammaPas = aux.rebin(wwgammaPas,binning,False)
        wwgammaTot = aux.rebin(wwgammaTot,binning,False)
        wzgammaPas = aux.rebin(wzgammaPas,binning,False)
        wzgammaTot = aux.rebin(wzgammaTot,binning,False)
        DYjetsNLOPas = aux.rebin(DYjetsNLOPas,binning,False)
        DYjetsNLOTot = aux.rebin(DYjetsNLOTot,binning,False)
        wjetsPas = aux.rebin(wjetsPas,binning,False)
        wjetsTot = aux.rebin(wjetsTot,binning,False)
        ttPas = aux.rebin(ttPas,binning,False)
        ttTot = aux.rebin(ttTot,binning,False)
        zz4lPas = aux.rebin(zz4lPas,binning,False)
        zz4lTot = aux.rebin(zz4lTot,binning,False)
        wzPas = aux.rebin(wzPas,binning,False)
        wzTot = aux.rebin(wzTot,binning,False)
        wwPas = aux.rebin(wwPas,binning,False)
        wwTot = aux.rebin(wwTot,binning,False)
        singletopPas = aux.rebin(singletopPas,binning,False)
        singletopTot = aux.rebin(singletopTot,binning,False)
        wgammaPas = aux.rebin(wgammaPas,binning,False)
        wgammaTot = aux.rebin(wgammaTot,binning,False)



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
        
        titleForNow = gr.GetXaxis().GetTitle()
        if "leading" in titleForNow:
            titleNow = titleForNow.replace("leading","leading lepton")
        if "trailing" in titleForNow:
            titleNow = titleForNow.replace("trailing","trailing lepton")
        gr.GetXaxis().SetTitle(titleNow)
        eff.SetTitle(";"+titleNow+";"+gr.GetYaxis().GetTitle()+"#varepsilon")
        
        #if additional:
            #for i in range(len(effAddArrayPas)):
                #print i
                ###xAdd = effAddArrayPas[i].Clone(aux.randomName())
                ###yAdd = effAddArrayTot[i].Clone(aux.randomName())
                #xAddArray.append(effAddArrayPas[i].Clone(aux.randomName()))
                #yAddArray.append(effAddArrayTot[i].Clone(aux.randomName()))
                #effAddArrayTemp.append(ROOT.TEfficiency(effAddArrayPas[i], effAddArrayTot[i]))
                #effAddArrayPas[i], effAddArrayTot[i] = xAddArray[i], yAddArray[i]
                #effAddArray[i].Draw("same")
                #effAddArray[i].SetLineColor(additional[i].color)
                #ROOT.gPad.Update()
                #grAddArray.append(effAddArray[i].GetPaintedGraph())
                #grAddArray[i].GetXaxis().SetTitle(titleNow)
                #grAddArray[i].SetTitle(";"+gr.GetXaxis().GetTitle()+";"+titleNow)
                
                
        #zgamma+ttgamma+zz+wwgamma+wzgamma+DYjetsNLO+wjets+tt+zz4l+wz+ww+singletop+wgamma
        x_zgamma = zgammaPas.Clone(aux.randomName())
        y_zgamma = zgammaTot.Clone(aux.randomName())
        x_ttgamma = ttgammaPas.Clone(aux.randomName())
        y_ttgamma = ttgammaTot.Clone(aux.randomName())
        x_zz = zzPas.Clone(aux.randomName())
        y_zz = zzTot.Clone(aux.randomName())
        x_wwgamma = wwgammaPas.Clone(aux.randomName())
        y_wwgamma = wwgammaTot.Clone(aux.randomName())
        x_wzgamma = wzgammaPas.Clone(aux.randomName())
        y_wzgamma = wzgammaTot.Clone(aux.randomName())
        x_DYjetsNLO = DYjetsNLOPas.Clone(aux.randomName())
        y_DYjetsNLO = DYjetsNLOTot.Clone(aux.randomName())
        x_wjets = wjetsPas.Clone(aux.randomName())
        y_wjets = wjetsTot.Clone(aux.randomName())
        x_tt = ttPas.Clone(aux.randomName())
        y_tt = ttTot.Clone(aux.randomName())
        x_zz4l = zz4lPas.Clone(aux.randomName())
        y_zz4l = zz4lTot.Clone(aux.randomName())
        x_wz = wzPas.Clone(aux.randomName())
        y_wz = wzTot.Clone(aux.randomName())
        x_ww = wwPas.Clone(aux.randomName())
        y_ww = wwTot.Clone(aux.randomName())
        x_singletop = singletopPas.Clone(aux.randomName())
        y_singletop = singletopTot.Clone(aux.randomName())
        x_wgamma = wgammaPas.Clone(aux.randomName())
        y_wgamma = wgammaTot.Clone(aux.randomName())

        eff_zgamma = ROOT.TEfficiency(zgammaPas,zgammaTot)
        eff_ttgamma = ROOT.TEfficiency(ttgammaPas,ttgammaTot)
        eff_zz = ROOT.TEfficiency(zzPas,zzTot)
        eff_wwgamma = ROOT.TEfficiency(wwgammaPas,wwgammaTot)
        eff_wzgamma = ROOT.TEfficiency(wzgammaPas,wzgammaTot)
        eff_DYjetsNLO = ROOT.TEfficiency(DYjetsNLOPas,DYjetsNLOTot)
        eff_wjets = ROOT.TEfficiency(wjetsPas,wjetsTot)
        eff_tt = ROOT.TEfficiency(ttPas,ttTot)
        eff_zz4l = ROOT.TEfficiency(zz4lPas,zz4lTot)
        eff_wz = ROOT.TEfficiency(wzPas,wzTot)
        eff_ww = ROOT.TEfficiency(wwPas,wwTot)
        eff_singletop = ROOT.TEfficiency(singletopPas,singletopTot)
        eff_wgamma = ROOT.TEfficiency(wgammaPas,wgammaTot)

        xAr = [x_zgamma,x_ttgamma,x_zz,x_wwgamma,x_wzgamma,x_DYjetsNLO,x_wjets,x_tt,x_zz4l,x_wz,x_ww,x_singletop,x_wgamma]
        yAr = [y_zgamma,y_ttgamma,y_zz,y_wwgamma,y_wzgamma,y_DYjetsNLO,y_wjets,y_tt,y_zz4l,y_wz,y_ww,y_singletop,y_wgamma]
        effAr = [eff_zgamma,eff_ttgamma,eff_zz,eff_wwgamma,eff_wzgamma,eff_DYjetsNLO,eff_wjets,eff_tt,eff_zz4l,eff_wz,eff_ww,eff_singletop,eff_wgamma]

        #print len(xAr),len(yAr),len(effAr), len(additional)

        for i in range(len(xAr)):
            arAddPas[i],arAddTot[i] = xAr[i],yAr[i]
            effAr[i].SetLineColor(additional[i].color)
            effAr[i].Draw("Same")
            #print effAr[i],xAr[i],yAr[i]
            #print effAr[i].GetEfficiency(10)
            ROOT.gPad.Update()
            #c.Modified()
            #c.Update()

    gr.GetYaxis().SetRangeUser(0.5, 1.1)


    #if additional:
        #for i in range(len(effAddArrayPas)):
            #grAddArray[i].GetYaxis().SetRangeUser(0.5, 1.1)
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
        #if additional:
            #binAdd = h_pasAdd.FindFixBin(cutValue)
            #passedAdd = int(h_pasAdd.Integral(binAdd, -1))
            #totalAdd = int(h_totAdd.Integral(binAdd, -1))
            #if not totalAdd: return
            #eAdd = 1.*passedAdd/totalAdd
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
        e_upSyst = e+0.03
        e_dnSyst = e-0.03
        e_syst=0.03
        if e_upSyst>1.:
            e_upSyst=1.
        eLabel = ROOT.TLatex(0.47, .17, "#varepsilon_{{Data}} = {:.1f}^{{#plus{:.1f}}}_{{#minus{:.1f}}}%".format(100*e, 100*(e_up-e),100*(e-e_dn)))
        eLabel.SetNDC()
        eLabel.Draw()
        #if additional:
            #if passedAdd<=totalAdd:
                #confAdd = ROOT.TEfficiency().GetConfidenceLevel()
                #e_upAdd = ROOT.TEfficiency.ClopperPearson(totalAdd, passedAdd, confAdd, True)
                #e_dnAdd = ROOT.TEfficiency.ClopperPearson(totalAdd, passedAdd, confAdd, False)
            #else:
                #passedAdd, epassedAdd = aux.integralAndError(h_pasAdd,binAdd,-1)
                #totalAdd, etotalAdd = aux.integralAndError(h_pasAdd,binAdd,-1)
                #eeAdd = eAdd * math.sqrt((epassedAdd/passedAdd)**2 + (etotalAdd/totalAdd)**2)
                #e_upAdd = eAdd + eeAdd
                #e_dnAdd = eAdd - eeAdd
            #eLabelAdd = ROOT.TLatex(0.47, .25, "#varepsilon_{{MC}} = {:.2f}^{{#plus{:.2f}}}_{{#minus{:.2f}}}%".format(100*eAdd, 100*(e_upAdd-eAdd),100*(eAdd-e_dnAdd)))
            #eLabelAdd.SetTextColor(kRed)
            #eLabelAdd.SetNDC()
            #eLabelAdd.Draw()

        # graphical representation
        l = ROOT.TLine()
        l.SetLineWidth(1)
        l.SetLineColor(ROOT.kGray+2)
        xmin = gr.GetHistogram().GetXaxis().GetXmin()
        xmax = gr.GetHistogram().GetXaxis().GetXmax()
        #l.DrawLine(max(cutValue,xmin), e, xmax, e)
        #l.DrawLine(max(cutValue,xmin), e_up, xmax, e_up)
        #l.DrawLine(max(cutValue,xmin), e_dn, xmax, e_dn)



        nBins=len(binning)
        binningMin=binning[0]
        binningMax=binning[-1]
        

        tempH = TH1F("","",nBins-1,binningMin,binningMax)
        #nPoints = grAdd.GetN()
        #for i in range(nPoints):
           #x=ROOT.Double(0)
           #y=ROOT.Double(0)
           #grAdd.GetPoint(i, x, y)
           #tempH.Fill(x,y) 



        #effAddSyst = aux.getSysHisto(tempH, 0.03)
        #aux.drawOpt(effAddSyst, "sysUnc")
        #effAddSyst.Draw("same e2")


        leg = ROOT.TLegend(.56,.49,.94,.615)
        leg.AddEntry( eff, "Data (ee)" if "EE" in savename else "Data (#mu#mu)" if "MM" in savename else "Data (e#mu)", "epl" )
        #if additional:
            #leg.AddEntry( effAdd, "MC (ee)" if "EE" in savename else "MC (#mu#mu)" if "MM" in savename else "MC (e#mu)", "epl" )
        linie = ROOT.TLine()
        linie.SetLineWidth(1)
        linie.SetLineColor(ROOT.kGray+2)
        #leg.AddEntry(linie,"mean_{Data} #pm 1 #sigma stat.","l")
        systLinie= ROOT.TLine()
        systLinie.SetLineWidth(5)
        systLinie.SetLineColor(ROOT.kRed-3)
        #leg.AddEntry(effAddSyst,"syst. unc.","f")
        #leg.Draw()
        
        leg2 = ROOT.TLegend(.56,.501,.94,.52)
        leg2.SetFillStyle(0)
        linie.SetLineWidth(1)
        leg2.AddEntry(linie," ","l")
        leg2.AddEntry(linie," ","l")


        if cutValue > xmin:
            # cut line
            l.SetLineStyle(2)
            ymin = gr.GetYaxis().GetXmin()
            ymax = gr.GetYaxis().GetXmax()
            l.DrawLine(cutValue, 0.5, cutValue, ymax)




    saveLumi = aux.intLumi
    if "B-F" in dataset.label:
        aux.intLumi = 20.101e3 # brilcalc lumi -b "STABLE BEAMS" --normtag=/afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt --end 278808 # Up to run F
    l = aux.Label(sim="Data" not in dataset.label)
    aux.intLumi = saveLumi
    #textPW="Work in Progress"
    #pw = ROOT.TLatex( 0.2, .887, "#scale[0.76]{#font[52]{%s}}"%textPW )
    #pw.Draw()
    name = "_"+name.split("/")[-1]
    if binningName: binningName = "_"+binningName
    #aux.save("efficiency_"+savename+name+binningName, log=False)
    aux.save("efficiency_"+savename+name+binningName, log=False,folder="triggerStudies/")

    if False:
        h_tot.SetLineColor(2)
        h_tot.Draw("hist")
        h_pas.Draw("same e*")
        #aux.save("efficiency_"+savename+name+"_raw")
        aux.save("efficiency_"+savename+name+"_raw",folder="triggerStudies/")



def main():

    #bkgs=[DYjets,zgamma,wwgamma,wzgamma,ttgamma,zz,tt,wjets]
    allMC = zgamma+ttgamma+zz+wwgamma+wzgamma+DYjetsNLO+wjets+tt+zz4l+wz+ww+singletop+wgamma
    #variables=["eta1","pt1","n_jets","n_vtx","phi1","m_ll","ht","n_photons"]
    variables=["pt1","pt2"]
    #groups=["trigDilep","trigSel","trigOnZ"]
    groups=["trigDilep"]
    #groups=["trigDilep_ptcuts","trigSel_ptcuts","trigOnZ_ptcuts"]
    for group in groups:
        for variable in variables:
            efficiency(dataHt,group+"/EE/"+variable,"dataHT_"+group+"_EE",binning=binnings[variable],additional=[allMC])
            efficiency(dataHt,group+"/MM/"+variable,"dataHT_"+group+"_MM",binning=binnings[variable],additional=[allMC])
            efficiency(dataHt,group+"/EM/"+variable,"dataHT_"+group+"_EM",binning=binnings[variable],additional=[allMC])
            
            efficiency(dataHt,group+"_ptcuts/EE/"+variable,"dataHT_"+group+"_ptcuts_EE",binning=binnings[variable],additional=[allMC])
            efficiency(dataHt,group+"_ptcuts/MM/"+variable,"dataHT_"+group+"_ptcuts_MM",binning=binnings[variable],additional=[allMC])
            efficiency(dataHt,group+"_ptcuts/EM/"+variable,"dataHT_"+group+"_ptcuts_EM",binning=binnings[variable],additional=[allMC])
            
            newBinnings=binnings.copy()
            newBinnings["pt1"]=frange(0., 302, 5.)
            newBinnings["pt2"]=frange(0., 102, 5.)
            
            #efficiencyAllMC(dataHt,group+"MM_ptcuts/"+variable,"test_"+group+"_ptcuts_MM",binning=newBinnings[variable],additional=[zgamma,ttgamma,zz,wwgamma,wzgamma,DYjetsNLO,wjets,tt,zz4l,wz,ww,singletop,wgamma])
            efficiencyAllMC(dataHt,group+"/MM/"+variable,"test_"+group+"_MM",binning=newBinnings[variable],additional=[zgamma,ttgamma,zz,wwgamma,wzgamma,DYjetsNLO,wjets,tt,zz4l,wz,ww,singletop,wgamma])
            efficiencyAllMC(dataHt,group+"/EM/"+variable,"test_"+group+"_EM",binning=newBinnings[variable],additional=[zgamma,ttgamma,zz,wwgamma,wzgamma,DYjetsNLO,wjets,tt,zz4l,wz,ww,singletop,wgamma])
            efficiencyAllMC(dataHt,group+"/EE/"+variable,"test_"+group+"_EE",binning=newBinnings[variable],additional=[zgamma,ttgamma,zz,wwgamma,wzgamma,DYjetsNLO,wjets,tt,zz4l,wz,ww,singletop,wgamma])


  
if __name__=="__main__":
    main()
