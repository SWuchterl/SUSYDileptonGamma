from dataMC import labels, frange, frangeN

import ROOT
from ROOT import *
from array import array
from include import *
import numpy as np
import pickle as pkl
import CR_tt
import CR_DY
import CR_WZ
import CR_ZZ
binnings = {
    'pt1':              frange(0, 100, 10) + frange(100, 200, 25) + range(200, 350, 50),
    'pt2':              frange(0, 100, 10) + frange(100, 200, 25),
    'pt3':              range(0, 200, 10),
    'pt4':              range(0, 200, 10),
    'eta1':             frange(0., 2.405, 0.1),
    'eta2':             frange(0., 2.4, 0.1),
    'eta3':             frange(0., 2.4, 0.1),
    'eta4':             frange(0., 2.4, 0.1),
    'phi1':             frange(0., 3.50, 0.1),
    'phi2':             frange(0., 3.50, 0.1),
    'phi3':             frange(0., 3.50, 0.1),
    'phi4':             frange(0., 3.50, 0.1),
    'ht':               frange(0., 1000., 50),
    # 'met':               [0,25,50,75,100,150,190,230,500],
    'met':               [75, 100, 150, 190, 230, 500],
    'm_ll':             frange(50., 100., 10) + frange(100., 300., 20.),
    'm_ll2':             frange(50., 100., 10) + frange(100., 300., 20.),
    'm_llg':             frange(0, 100, 10) + frange(100, 200, 10) + frange(200, 500, 50),
    'pt_llg':             frange(0, 100, 10) + frange(100, 200, 10) + frange(200, 500, 50),
    'n_jets':           frange(0., 7., 1),
    'n_photons':         frange(0., 4., 1),
    'n_vtx':            frange(0., 40., 1),
    'pt_g1':            frange(20, 100, 10) + frange(100, 150, 25) + frange(150, 250, 50),
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
    'mTL3Met': frange(0., 300., 10)
}


def calculateChiSquared(dataHist, mcHists):
    obs = []
    expected = []
    #print dataHist,mcHists
    # for bin in range(dataHist.GetNbinsX()+1):
    for bin in range(dataHist.GetNbinsX()):
        bin = bin + 1
        bw = dataHist.GetBinWidth(bin) if style.divideByBinWidth else 1.
        #bw2 = mcHists[0].GetBinWidth(bin) if style.divideByBinWidth else 1.

        #print dataHist.GetBinLowEdge(bin)

        obs.append(dataHist.GetBinContent(bin) * bw)
        expectedTemp = 0.
        for histo in mcHists:
            expectedTemp += (histo.GetBinContent(bin) * bw)
            #print histo,histo.GetBinContent(bin)*bw
        expected.append(expectedTemp)
    degOFr = dataHist.GetNbinsX() - 1.
    chiq2 = 0.
    #print "obs",obs
    #print "exp", expected
    for i in range(len(obs)):
        chiq2 += ((obs[i] - expected[i])**2. / (expected[i]))
    # chiq2=chiq2/(degOFr)
    #print chiq2
    return [chiq2, degOFr]


def divideDatasetIntegrals(numerator, denominator, name):
    numMerged = sum(numerator)
    h_num = numMerged.getHist(name)
    num = h_num.Integral(0, -1)
    denMerged = sum(denominator)
    h_den = denMerged.getHist(name)
    den = h_den.Integral(0, -1)
    return num / den if den else 1.


def calculateSFAndError(numerator_data, denominator_toScale, additional_fix, name):
    numMerged_data = sum(numerator_data)
    h_num_data = numMerged_data.getHist(name)
    num_dataErr = (ROOT.Double(0))
    num_data = h_num_data.IntegralAndError(0, -1, num_dataErr)

    denMerged_toScale = sum(denominator_toScale)
    h_den_toScale = denMerged_toScale.getHist(name)
    den_toScaleErr = (ROOT.Double(0))
    den_toScale = h_den_toScale.IntegralAndError(0, -1, den_toScaleErr)

    addMerged_fix = sum(additional_fix)
    h_add_fix = addMerged_fix.getHist(name)
    add_fixError = ROOT.Double(0)
    add_fix = h_add_fix.IntegralAndError(0, -1, add_fixError)

    alpha = (num_data - add_fix) / (den_toScale)

    alphaErr = np.sqrt((num_dataErr / den_toScale)**2. + (add_fixError / den_toScale)
                       ** 2. + (den_toScaleErr * (num_data - add_fix) / (den_toScale)**2.)**2.)

    return [alpha, alphaErr / alpha] if den_toScale else [1., 0.]


def drawSameHistogram(sampleNames, name, bkg=[], additional=[], binning=None, binningName="", scaleToData=True, xTitle=None, yTitle=None, toScale=[], SF=1.):
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()

    style.divideByBinWidth = False
    style.minimumOne = True

    scale = SF
    scaleErr = 0.

    addBKG = bkg[:]
    for y in toScale:
        addBKG.remove(y)

    #folder= (name.split("/"))[0]
    folder = (name.split("/"))[0] + "/" + \
        (name.split("/"))[1] + "/" + (name.split("/"))[2]
    #print "folder",folder

    #print name,folder

    for d in bkg[-1::-1]:
        h = d.getHistWithWeights(name, ["nISR", "topPt", "ewk"])

        if not h:
            continue
        if not h.Integral():
            continue

        if d in toScale:
            h.Scale(scale)

        if binning:
            h = aux.rebin(h, binning)
        aux.appendFlowBin(h)

        m.addStack(h, d.label)

    dataHist = None
    for d in additional:
        h = d.getHist(name)
        if not h:
            continue
        if not h.Integral():
            continue
        if (binning):
            h = aux.rebin(h, binning)
        aux.appendFlowBin(h)

        if h.GetLineColor() == ROOT.kBlack:  # data
            if not(binning):
                h.Sumw2(False)
            h.SetBinErrorOption(ROOT.TH1.kPoisson)
            dataHist = h
        else:
            h.drawOption_ = "hist e"
            h.SetLineWidth(3)

        m.add(h, d.label)

        m.sortStackByIntegral()

    totalMC = m.histsToStack
    totalData = m.hists[0]
    return(calculateChiSquared(totalData, totalMC))


def drawDY(binningToUse_=CR_DY.binnings.copy(), addName=""):
    bkgs = [DYjetsNLO, zgamma, tt, ttgamma, wwgamma,
            wzgamma, zz, wjets, wgamma, singletop, zz4l, wz, ww]
    # variables=["eta1","eta2","pt1","pt2","n_jets","n_vtx","phi1","phi2","m_ll","ht","n_photons","pt_g1","met","nElectrons","nMuons","deltaR1_g1","deltaR2_g1","deltaRLL"]
    # variables=["eta1","eta2","pt1","pt2","met","pt_g1","phi1","phi2","n_jets","ht"]
    # variables=["eta1","eta2","pt1","pt2","pt_g1","m_ll","phi1","phi2","n_jets","ht"]
    variables = ["eta1", "eta2", "pt1", "pt2",
                 "pt_g1", "phi1", "phi2", "n_jets", "ht", "met"]
    # variables=["pt1"]
    groups = ["CRDY"]
    binnings_ = CR_DY.binnings.copy()

    additionalFolder = "nom/"
    if "JERu" in addName:
        additionalFolder = "JESu/"
    if "JESd" in addName:
        additionalFolder = "JESd/"
    if "JERu" in addName:
        additionalFolder = "JERu/"
    if "JERd" in addName:
        additionalFolder = "JERd/"
    if "lepSFu" in addName:
        additionalFolder = "lepSFu/"
    if "lepSFd" in addName:
        additionalFolder = "lepSFd/"
    if "photonSFu" in addName:
        additionalFolder = "photonSFu/"
    if "photonSFd" in addName:
        additionalFolder = "photonSFd/"

    scaleFactors = frange(1.05, 1.09, 0.002)
    # saveValuesEE={}
    # saveValuesMM={}
    saveValuesLL = {}
    for group in groups:
        # saveValuesEE[group]={}
        # saveValuesMM[group]={}
        saveValuesLL[group] = {}
        for variable in variables:
            # saveValuesEE[group][variable]={}
            # saveValuesMM[group][variable]={}
            saveValuesLL[group][variable] = {}
            # chiValuesEE=[]
            # chiValuesMM=[]
            chiValuesLL = []
            for SF in scaleFactors:
                #tempEE=(drawSameHistogram("EE",group+"/EE/"+additionalFolder+variable, bkgs, additional=[dataDoubleEG],binning=binningToUse_[variable],xTitle=labels[variable][0],toScale=[DYjetsNLO,zgamma],SF=SF))
                #tempMM=(drawSameHistogram("MM",group+"/MM/"+additionalFolder+variable, bkgs, additional=[dataDoubleMuon],binning=binningToUse_[variable],xTitle=labels[variable][0],toScale=[DYjetsNLO,zgamma],SF=SF))
                tempLL = (drawSameHistogram("LL", group + "/LL/" + additionalFolder + variable, bkgs, additional=[
                          dataDoubleSF], binning=binningToUse_[variable], xTitle=labels[variable][0], toScale=[DYjetsNLO, zgamma], SF=SF))
                # chiValuesEE.append(tempEE[0])
                # chiValuesMM.append(tempMM[0])
                chiValuesLL.append(tempLL[0])
            #grEE = TGraph();
            #grMM = TGraph();
            grLL = TGraph()

            # for i in range(len(chiValuesEE)):
            # grEE.SetPoint(i,scaleFactors[i],chiValuesEE[i])
            # for i in range(len(chiValuesMM)):
            # grMM.SetPoint(i,scaleFactors[i],chiValuesMM[i])
            for i in range(len(chiValuesLL)):
                grLL.SetPoint(i, scaleFactors[i], chiValuesLL[i])

            c = TCanvas("canvas", "", 800, 800)

            # chiValuesEE=np.array(chiValuesEE)
            # chiValuesMM=np.array(chiValuesMM)
            chiValuesLL = np.array(chiValuesLL)
            scaleFactors = np.array(scaleFactors)

            # findMinIndexEE=chiValuesEE.argmin()
            # plusOneLeftEE=chiValuesEE[findMinIndexEE]+1
            # plusOneRightEE=chiValuesEE[findMinIndexEE]+1

            # findMinIndexMM=chiValuesMM.argmin()
            # plusOneLeftMM=chiValuesMM[findMinIndexMM]+1
            # plusOneRightMM=chiValuesMM[findMinIndexMM]+1

            findMinIndexLL = chiValuesLL.argmin()
            plusOneLeftLL = chiValuesLL[findMinIndexLL] + 1
            plusOneRightLL = chiValuesLL[findMinIndexLL] + 1

            # upEE=0.
            # downEE=0.
            # upMM=0.
            # downMM=0.
            upLL = 0.
            downLL = 0.

            # grEE.Fit("pol4")
            # grEE.GetFunction("pol4").SetLineColor(kRed)
            # grEE.GetFunction("pol4").SetLineWidth(1)

            # fitEE=grEE.GetFunction("pol4")

            # chiResultEE=fitEE.GetMinimum()
            # sfResultEE=fitEE.GetMinimumX()
            # sfResultUpEE=fitEE.GetX(chiResultEE+1,0.,sfResultEE)
            # sfResultDnEE=fitEE.GetX(chiResultEE+1,sfResultEE,100.)

            # upEE=abs(sfResultEE-sfResultUpEE)
            # downEE=abs(sfResultEE-sfResultDnEE)

            # grMM.Fit("pol4")
            # grMM.GetFunction("pol4").SetLineColor(kRed)
            # grMM.GetFunction("pol4").SetLineWidth(1)

            # fitMM=grMM.GetFunction("pol4")

            # chiResultMM=fitMM.GetMinimum()
            # sfResultMM=fitMM.GetMinimumX()
            # sfResultUpMM=fitMM.GetX(chiResultMM+1,0.,sfResultMM)
            # sfResultDnMM=fitMM.GetX(chiResultMM+1,sfResultMM,100.)

            # upMM=abs(sfResultMM-sfResultUpMM)
            # downMM=abs(sfResultMM-sfResultDnMM)

            grLL.Fit("pol4")
            grLL.GetFunction("pol4").SetLineColor(kRed)
            grLL.GetFunction("pol4").SetLineWidth(1)

            fitLL = grLL.GetFunction("pol4")

            chiResultLL = fitLL.GetMinimum()
            sfResultLL = fitLL.GetMinimumX()
            sfResultUpLL = fitLL.GetX(chiResultLL + 1, 0., sfResultLL)
            sfResultDnLL = fitLL.GetX(chiResultLL + 1, sfResultLL, 100.)

            upLL = abs(sfResultLL - sfResultUpLL)
            downLL = abs(sfResultLL - sfResultDnLL)

            # grEE.SetTitle("; Scale Factor #alpha; #chi^{2}")
            # grEE.SetMarkerStyle(20)
            # grEE.SetMarkerSize(0.5)
            # grEE.SetName("grDYEE")
            # grMM.SetTitle("; Scale Factor #alpha; #chi^{2}")
            # grMM.SetMarkerStyle(20)
            # grMM.SetMarkerSize(0.5)
            # grMM.SetName("grDYMM")

            grLL.SetTitle("; Scale Factor #alpha; #chi^{2}")
            grLL.SetMarkerStyle(20)
            grLL.SetMarkerSize(0.5)
            grLL.SetName("grDYLL")

            # grEE.Draw("ACP")
            # l = ROOT.TLatex( 0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
            # l.SetNDC()
            # l.Draw()
            # l2 = ROOT.TLatex( 0.47, .65, "#scale[0.66]{#font[52]{#alpha = %.3f^{#plus%.3f}_{#minus%.3f}}}"%(sfResultEE,upEE,downEE))
            # l2.SetNDC()
            # l2.Draw()
            # l3 = ROOT.TLatex( 0.47, .6, "#scale[0.66]{#font[52]{#chi}^{2} = %.2f / %.0f}"%(chiResultEE,tempEE[1]))
            # l3.SetNDC()
            # l3.Draw()
            # l4 = ROOT.TLatex( 0.47, .8, "#scale[0.66]{#font[52]{DY Control Region ee}}")
            # l4.SetNDC()
            # l4.Draw()
            #lum = ROOT.TLatex( .62, .95, "%.1f fb^{-1} (%s TeV)"%(aux.intLumi/1000., aux.Label.cmsEnergy) )
            # lum.SetNDC()
            # lum.Draw()

            # leg=TLegend(0.5,0.8,0.7,0.9)
            #leg.AddEntry("grDYEE","measured points","p")
            #leg.AddEntry(grEE.GetFunction("pol4"),"polynomial fit")
            # leg.Draw()

            # c.Update()
            # c.SaveAs('plots_CR/chi/DY_EE'+'_'+variable+addName+'.pdf')
            # c.SaveAs('plots_CR/chi/DY_EE'+'_'+variable+addName+'.root')

            # grMM.Draw("ACP")
            # l = ROOT.TLatex( 0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
            # l.SetNDC()
            # l.Draw()
            # l2 = ROOT.TLatex( 0.47, .65, "#scale[0.66]{#font[52]{#alpha = %.3f^{#plus%.3f}_{#minus%.3f}}}"%(sfResultMM,upMM,downMM))
            # l2.SetNDC()
            # l2.Draw()
            # l3 = ROOT.TLatex( 0.47, .6, "#scale[0.66]{#font[52]{#chi}^{2} = %.2f / %.0f}"%(chiResultMM,tempMM[1]))
            # l3.SetNDC()
            # l3.Draw()
            # l4 = ROOT.TLatex( 0.47, .7, "#scale[0.66]{#font[52]{DY Control Region #mu#mu}}")
            # l4.SetNDC()
            # l4.Draw()
            #lum = ROOT.TLatex( .62, .95, "%.1f fb^{-1} (%s TeV)"%(aux.intLumi/1000., aux.Label.cmsEnergy) )
            # lum.SetNDC()
            # lum.Draw()

            # leg=TLegend(0.5,0.8,0.7,0.9)
            #leg.AddEntry("grDYMM","measured points","p")
            #leg.AddEntry(grMM.GetFunction("pol4"),"polynomial fit")
            # leg.Draw()

            # c.Update()
            # c.SaveAs('plots_CR/chi/DY_MM'+'_'+variable+addName+'.pdf')
            # c.SaveAs('plots_CR/chi/DY_MM'+'_'+variable+addName+'.root')

            grLL.Draw("ACP")
            l = ROOT.TLatex(
                0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
            l.SetNDC()
            l.Draw()
            l2 = ROOT.TLatex(0.47, .65, "#scale[0.66]{#font[52]{#alpha = %.3f^{#plus%.3f}_{#minus%.3f}}}" % (
                sfResultLL, upLL, downLL))
            l2.SetNDC()
            l2.Draw()
            l3 = ROOT.TLatex(
                0.47, .6, "#scale[0.66]{#font[52]{#chi}^{2}/ndf = %.2f / %.0f}" % (chiResultLL, tempLL[1]))
            l3.SetNDC()
            l3.Draw()
            l4 = ROOT.TLatex(
                0.47, .7, "#scale[0.66]{#font[52]{DY/Z(#gamma) Control Region ee+#mu#mu}}")
            l4.SetNDC()
            l4.Draw()
            lum = ROOT.TLatex(.62, .95, "%.1f fb^{-1} (%s TeV)" %
                              (aux.intLumi / 1000., aux.Label.cmsEnergy))
            lum.SetNDC()
            lum.Draw()

            leg = TLegend(0.5, 0.8, 0.7, 0.9)
            leg.AddEntry("grDYLL", "measured points", "p")
            leg.AddEntry(grLL.GetFunction("pol4"), "polynomial fit")
            leg.Draw()

            c.Update()
            c.SaveAs('plots_CR/chi/DY_LL' + '_' + variable + addName + '.pdf')
            c.SaveAs('plots_CR/chi/DY_LL' + '_' + variable + addName + '.root')

            bestValueLL = sfResultLL
            # bestValueEE=sfResultEE
            # bestValueMM=sfResultMM
            erUpLL = upLL
            # erUpEE=upEE
            # erUpMM=upMM
            erDownLL = downLL
            # erDownEE=downEE
            # erDownMM=downMM
            # saveValuesEE[group][variable]["value"]=bestValueEE
            # saveValuesEE[group][variable]["erUp"]=erUpEE
            # saveValuesEE[group][variable]["erDown"]=erDownEE
            # saveValuesMM[group][variable]["value"]=bestValueMM
            # saveValuesMM[group][variable]["erUp"]=erUpMM
            # saveValuesMM[group][variable]["erDown"]=erDownMM
            saveValuesLL[group][variable]["value"] = bestValueLL
            saveValuesLL[group][variable]["erUp"] = erUpLL
            saveValuesLL[group][variable]["erDown"] = erDownLL
    #pkl.dump( saveValuesEE, open( "plots_CR/chi/DY_chi_EE"+addName+".pkl", "wb" ) )
    pkl.dump(saveValuesLL, open(
        "plots_CR/chi/DY_chi_LL" + addName + ".pkl", "wb"))
    #pkl.dump( saveValuesMM, open( "plots_CR/chi/DY_chi_MM"+addName+".pkl", "wb" ) )


def drawZZ(binningToUse_=CR_ZZ.binnings.copy(), addName=""):
    # variables=["eta1","eta2","eta3","eta4","pt1","pt3","pt4","pt2","n_jets","n_vtx","phi1","phi2","phi3","phi4","m_ll","m_ll2","ht","n_photons","met","nElectrons","nMuons"]
    # variables=["eta1"]
    # variables=["eta1","eta2","pt1","pt2","met","m_ll","phi1","phi2","n_jets","ht"]
    # variables=["eta1","eta2","pt1","pt2","pt_g1","m_ll","phi1","phi2","n_jets","ht"]
    variables = ["eta1", "eta2", "pt1", "pt2", "phi1",
                 "phi2", "n_jets", "ht", "m_ll", "m_ll2", "met"]
    # variables=["pt1"]
    bkgs = [DYjetsNLO, zgamma, tt, ttgamma, wwgamma,
            wzgamma, zz, wjets, wgamma, singletop, zz4l, wz, ww]
    groups = ["CRZZ"]
    scaleFactors = frange(0.9, 1.4, 0.05)
    chiValues = []
    binnings_ = CR_ZZ.binnings.copy()

    additionalFolder = "nom/"
    if "JERu" in addName:
        additionalFolder = "JESu/"
    if "JESd" in addName:
        additionalFolder = "JESd/"
    if "JERu" in addName:
        additionalFolder = "JERu/"
    if "JERd" in addName:
        additionalFolder = "JERd/"
    if "lepSFu" in addName:
        additionalFolder = "lepSFu/"
    if "lepSFd" in addName:
        additionalFolder = "lepSFd/"
    if "photonSFu" in addName:
        additionalFolder = "photonSFu/"
    if "photonSFd" in addName:
        additionalFolder = "photonSFd/"

    saveValues = {}
    for group in groups:
        saveValues[group] = {}
        for variable in variables:
            saveValues[group][variable] = {}
            chiValues = []
            for SF in scaleFactors:
                tempLL = (drawSameHistogram("LL", group + "/LL/" + additionalFolder + variable, bkgs, additional=[
                          dataLL], binning=binningToUse_[variable], xTitle=labels[variable][0], toScale=[zz4l], SF=SF))
                chiValues.append(tempLL[0])
            gr = TGraph()

            for i in range(len(chiValues)):
                gr.SetPoint(i, scaleFactors[i], chiValues[i])

            c = TCanvas("canvas", "", 800, 800)

            chiValues = np.array(chiValues)
            scaleFactors = np.array(scaleFactors)

            up = 0.
            down = 0.

            gr.Fit("pol4")
            gr.GetFunction("pol4").SetLineColor(kRed)
            gr.GetFunction("pol4").SetLineWidth(1)

            fit = gr.GetFunction("pol4")

            chiResult = fit.GetMinimum()
            sfResult = fit.GetMinimumX()
            sfResultUp = fit.GetX(chiResult + 1, 0., sfResult)
            sfResultDn = fit.GetX(chiResult + 1, sfResult, 100.)

            up = abs(sfResult - sfResultUp)
            down = abs(sfResult - sfResultDn)

            gr.SetTitle("; Scale Factor #alpha; #chi^{2}")
            gr.SetName("gr")
            gr.SetMarkerStyle(20)
            gr.SetMarkerSize(0.7)

            gr.Draw("ACP")
            l = ROOT.TLatex(
                0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
            l.SetNDC()
            l.Draw()
            l2 = ROOT.TLatex(0.47, .65, "#scale[0.66]{#font[52]{#alpha = %.3f^{#plus%.3f}_{#minus%.3f}}}" % (
                sfResult, up, down))
            l2.SetNDC()
            l2.Draw()
            l3 = ROOT.TLatex(
                0.47, .6, "#scale[0.66]{#font[52]{#chi}^{2}/ndf = %.2f / %.0f}" % (chiResult, tempLL[1]))
            l3.SetNDC()
            l3.Draw()
            l4 = ROOT.TLatex(
                0.47, .7, "#scale[0.66]{#font[52]{ZZ Control Region}}")
            l4.SetNDC()
            l4.Draw()
            lum = ROOT.TLatex(.62, .95, "%.1f fb^{-1} (%s TeV)" %
                              (aux.intLumi / 1000., aux.Label.cmsEnergy))
            lum.SetNDC()
            lum.Draw()

            leg = TLegend(0.5, 0.8, 0.7, 0.9)
            leg.AddEntry("gr", "measured points", "p")
            leg.AddEntry(gr.GetFunction("pol4"), "polynomial fit")
            leg.Draw()

            c.Update()
            c.SaveAs('plots_CR/chi/ZZ' + '_' + variable + addName + '.pdf')
            c.SaveAs('plots_CR/chi/ZZ' + '_' + variable + addName + '.root')

            bestValue = sfResult
            erUp = up
            erDown = down
            saveValues[group][variable]["value"] = bestValue
            saveValues[group][variable]["erUp"] = erUp
            saveValues[group][variable]["erDown"] = erDown
    pkl.dump(saveValues, open("plots_CR/chi/ZZ_chi" + addName + ".pkl", "wb"))


def drawWZ(binningToUse_=CR_WZ.binnings.copy(), addName=""):
    # variables=["eta1","eta2","eta3","pt1","pt2","pt3","n_jets","n_vtx","phi1","phi2","phi3","m_ll","ht","n_photons","met","nElectrons","nMuons","mTL3Met","met"]

    # variables=["eta1","eta2","pt1","pt2","met","m_ll","phi1","phi2","n_jets","ht"]
    # variables=["eta1","eta2","pt1","pt2","pt_g1","m_ll","phi1","phi2","n_jets","ht"]

    variables = ["eta1", "eta2", "pt1", "pt2",
                 "phi1", "phi2", "n_jets", "ht", "m_ll", "met"]
    # variables=["pt1"]
    bkgs = [DYjetsNLO, zgamma, tt, ttgamma, wwgamma,
            wzgamma, zz, wjets, wgamma, singletop, zz4l, wz, ww]
    groups = ["CRWZ"]
    scaleFactors = frange(0.9, 1.4, 0.05)
    chiValues = []
    binnings_ = CR_WZ.binnings.copy()

    additionalFolder = "nom/"
    if "JERu" in addName:
        additionalFolder = "JESu/"
    if "JESd" in addName:
        additionalFolder = "JESd/"
    if "JERu" in addName:
        additionalFolder = "JERu/"
    if "JERd" in addName:
        additionalFolder = "JERd/"
    if "lepSFu" in addName:
        additionalFolder = "lepSFu/"
    if "lepSFd" in addName:
        additionalFolder = "lepSFd/"
    if "photonSFu" in addName:
        additionalFolder = "photonSFu/"
    if "photonSFd" in addName:
        additionalFolder = "photonSFd/"

    saveValues = {}
    for group in groups:
        saveValues[group] = {}
        for variable in variables:
            saveValues[group][variable] = {}
            chiValues = []
            for SF in scaleFactors:
                tempLL = (drawSameHistogram("LL", group + "/LL/" + additionalFolder + variable, bkgs, additional=[
                          dataLL], binning=binningToUse_[variable], xTitle=labels[variable][0], toScale=[wz], SF=SF))
                chiValues.append(tempLL[0])
            gr = TGraph()

            for i in range(len(chiValues)):
                gr.SetPoint(i, scaleFactors[i], chiValues[i])

            c = TCanvas("canvas", "", 800, 800)

            chiValues = np.array(chiValues)
            scaleFactors = np.array(scaleFactors)

            up = 0.
            down = 0.

            gr.Fit("pol4")
            gr.GetFunction("pol4").SetLineColor(kRed)
            gr.GetFunction("pol4").SetLineWidth(1)

            fit = gr.GetFunction("pol4")

            chiResult = fit.GetMinimum()
            sfResult = fit.GetMinimumX()
            sfResultUp = fit.GetX(chiResult + 1, 0., sfResult)
            sfResultDn = fit.GetX(chiResult + 1, sfResult, 100.)

            up = abs(sfResult - sfResultUp)
            down = abs(sfResult - sfResultDn)

            gr.SetTitle("; Scale Factor #alpha; #chi^{2}")
            gr.SetMarkerStyle(20)
            gr.SetMarkerSize(0.7)
            gr.SetName("grWZ")

            gr.Draw("ACP")
            l = ROOT.TLatex(
                0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
            l.SetNDC()
            l.Draw()
            l2 = ROOT.TLatex(0.47, .65, "#scale[0.66]{#font[52]{#alpha = %.3f^{#plus%.3f}_{#minus%.3f}}}" % (
                sfResult, up, down))
            l2.SetNDC()
            l2.Draw()
            l3 = ROOT.TLatex(
                0.47, .6, "#scale[0.66]{#font[52]{#chi}^{2}/ndf = %.2f / %.0f}" % (chiResult, tempLL[1]))
            l3.SetNDC()
            l3.Draw()
            l4 = ROOT.TLatex(
                0.47, .7, "#scale[0.66]{#font[52]{WZ Control Region}}")
            l4.SetNDC()
            l4.Draw()
            lum = ROOT.TLatex(.62, .95, "%.1f fb^{-1} (%s TeV)" %
                              (aux.intLumi / 1000., aux.Label.cmsEnergy))
            lum.SetNDC()
            lum.Draw()

            leg = TLegend(0.5, 0.8, 0.7, 0.9)
            leg.AddEntry("grWZ", "measured points", "p")
            leg.AddEntry(gr.GetFunction("pol4"), "polynomial fit")
            leg.Draw()

            c.Update()
            c.SaveAs('plots_CR/chi/WZ' + '_' + variable + addName + '.pdf')
            c.SaveAs('plots_CR/chi/WZ' + '_' + variable + addName + '.root')

            bestValue = sfResult
            erUp = up
            erDown = down
            saveValues[group][variable]["value"] = bestValue
            saveValues[group][variable]["erUp"] = erUp
            saveValues[group][variable]["erDown"] = erDown
    pkl.dump(saveValues, open("plots_CR/chi/WZ_chi" + addName + ".pkl", "wb"))


def drawTT(binningToUse_=CR_tt.binnings.copy(), addName=""):
    # variables=["eta1","eta2","pt1","pt2","n_jets","n_vtx","phi1","phi2","m_ll","ht","n_photons","pt_g1","met","nElectrons","nMuons","deltaR1_g1","deltaR2_g1","deltaRLL"]
    variables = ["eta1", "eta2", "pt1", "pt2", "met",
                 "pt_g1", "m_ll", "phi1", "phi2", "n_jets", "ht"]
    # variables=["eta1","eta2","pt1","pt2","pt_g1","m_ll","phi1","phi2","n_jets","ht"]
    # variables=["eta1","eta2","pt1","pt2","pt_g1","phi1","phi2","n_jets","ht"]
    # variables=["pt1"]
    bkgs = [DYjetsNLO, zgamma, tt, ttgamma, wwgamma,
            wzgamma, zz, wjets, wgamma, singletop, zz4l, wz, ww]
    groups = ["CRTT"]
    binnings_ = CR_tt.binnings.copy()

    # scaleFactors = frange(0.7, 0.9, 0.01)
    scaleFactors = frange(0.7, 0.95, 0.01)

    #print binningToUse_["pt1"]

    additionalFolder = "nom/"
    if "JERu" in addName:
        additionalFolder = "JESu/"
    if "JESd" in addName:
        additionalFolder = "JESd/"
    if "JERu" in addName:
        additionalFolder = "JERu/"
    if "JERd" in addName:
        additionalFolder = "JERd/"
    if "lepSFu" in addName:
        additionalFolder = "lepSFu/"
    if "lepSFd" in addName:
        additionalFolder = "lepSFd/"
    if "photonSFu" in addName:
        additionalFolder = "photonSFu/"
    if "photonSFd" in addName:
        additionalFolder = "photonSFd/"

    #print additionalF

    saveValues = {}
    for group in groups:
        saveValues[group] = {}
        for variable in variables:
            # saveValues[group][variable+addName]={}
            saveValues[group][variable] = {}
            chiValues = []
            # chiValuesSmaller=[]
            # chiValuesBigger=[]
            for SF in scaleFactors:
                #tempEM=(drawSameHistogram("EM",group+"EM/"+variable, bkgs, additional=[dataMuonEG],binning=binnings_[variable],xTitle=labels[variable][0],toScale=[tt,ttgamma],SF=SF))
                #tempEM=(drawSameHistogram("EM",group+"EM/"+variable, bkgs, additional=[dataLL],binning=binnings_[variable],xTitle=labels[variable][0],toScale=[tt,ttgamma],SF=SF))
                tempEM = (drawSameHistogram("EM", group + "/EM/" + additionalFolder + variable, bkgs, additional=[
                          dataLL], binning=binningToUse_[variable], xTitle=labels[variable][0], toScale=[tt, ttgamma], SF=SF))
                #tempEMsmaller=(drawSameHistogram("EM",group+"EM/"+variable, bkgs, additional=[dataLL],binning=smallerBinning[variable],xTitle=labels[variable][0],toScale=[tt,ttgamma],SF=SF))
                #tempEMbigger=(drawSameHistogram("EM",group+"EM/"+variable, bkgs, additional=[dataLL],binning=biggerBinning[variable],xTitle=labels[variable][0],toScale=[tt,ttgamma],SF=SF))
                chiValues.append(tempEM[0])
            gr = TGraph()

            for i in range(len(chiValues)):
                gr.SetPoint(i, scaleFactors[i], chiValues[i])

            c = TCanvas("canvas", "", 800, 800)

            chiValues = np.array(chiValues)
            scaleFactors = np.array(scaleFactors)

            findMinIndex = chiValues.argmin()
            plusOneLeft = chiValues[findMinIndex] + 1
            plusOneRight = chiValues[findMinIndex] + 1
            print chiValues[findMinIndex], scaleFactors[findMinIndex]

            up = 0.
            down = 0.

            gr.Fit("pol4")
            gr.GetFunction("pol4").SetLineColor(kRed)
            gr.GetFunction("pol4").SetLineWidth(1)

            fit = gr.GetFunction("pol4")

            chiResult = fit.GetMinimum()
            sfResult = fit.GetMinimumX()
            sfResultUp = fit.GetX(chiResult + 1, 0., sfResult)
            sfResultDn = fit.GetX(chiResult + 1, sfResult, 100.)

            up = abs(sfResult - sfResultUp)
            down = abs(sfResult - sfResultDn)

            gr.SetTitle("; Scale Factor #alpha; #chi^{2}")
            gr.SetMarkerStyle(20)
            gr.SetMarkerSize(0.5)
            gr.SetName("grTT")

            gr.Draw("ACP")
            l = ROOT.TLatex(
                0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
            l.SetNDC()
            l.Draw()
            l2 = ROOT.TLatex(0.47, .65, "#scale[0.66]{#font[52]{#alpha = %.3f^{#plus%.3f}_{#minus%.3f}}}" % (
                sfResult, up, down))
            l2.SetNDC()
            l2.Draw()
            l3 = ROOT.TLatex(
                0.47, .6, "#scale[0.66]{#font[52]{#chi}^{2}/ndf = %.2f / %.0f}" % (chiResult, tempEM[1]))
            l3.SetNDC()
            l3.Draw()
            l4 = ROOT.TLatex(
                0.47, .7, "#scale[0.66]{#font[52]{t#bar{t}(#gamma) Control Region}}")
            l4.SetNDC()
            l4.Draw()
            lum = ROOT.TLatex(.62, .95, "%.1f fb^{-1} (%s TeV)" %
                              (aux.intLumi / 1000., aux.Label.cmsEnergy))
            lum.SetNDC()
            lum.Draw()

            leg = TLegend(0.5, 0.8, 0.7, 0.9)
            leg.AddEntry("grTT", "measured points", "p")
            leg.AddEntry(gr.GetFunction("pol4"), "polynomial fit")
            leg.Draw()

            c.Update()
            c.SaveAs('plots_CR/chi/TT' + '_' + variable + addName + '.pdf')
            c.SaveAs('plots_CR/chi/TT' + '_' + variable + addName + '.root')

            bestValue = scaleFactors[findMinIndex]
            erUp = up
            erDown = down
            # saveValues[group][variable+addName]["value"]=bestValue
            # saveValues[group][variable+addName]["erUp"]=erUp
            # saveValues[group][variable+addName]["erDown"]=erDown
            saveValues[group][variable]["value"] = bestValue
            saveValues[group][variable]["erUp"] = erUp
            saveValues[group][variable]["erDown"] = erDown
    pkl.dump(saveValues, open("plots_CR/chi/TT_chi" + addName + ".pkl", "wb"))


def main():
    import style
    style.defaultStyle()

    smallerBinningTT = CR_tt.binnings.copy()
    smallerBinningTT["eta1"] = frange(0., 2.4, 0.05)
    smallerBinningTT["eta2"] = frange(0., 2.4, 0.05)
    smallerBinningTT["pt1"] = frange(25., 150, 5.)
    smallerBinningTT["pt2"] = frange(20., 150, 5.)
    #smallerBinningTT["pt_g1"]  = frange(20,200,10)
    smallerBinningTT["pt_g1"] = frange(20, 200, 5)
    smallerBinningTT["phi1"] = frange(0., 3.1, 0.05)
    smallerBinningTT["phi2"] = frange(0., 3.1, 0.05)
    smallerBinningTT["n_jets"] = frange(0., 7, 7)
    smallerBinningTT["ht"] = frange(0., 1000., 25)
    smallerBinningTT["met"] = frange(0., 450., 25)
    smallerBinningTT["m_ll"] = frange(50., 300., 10.)

    nominalBinningTT = CR_tt.binnings.copy()
    nominalBinningTT["eta1"] = frange(0., 2.4, 0.1)
    nominalBinningTT["eta2"] = frange(0., 2.4, 0.1)
    nominalBinningTT["pt1"] = frange(20., 150, 10.)
    nominalBinningTT["pt2"] = frange(20., 150, 10.)
    #nominalBinningTT["pt_g1"]  = frange(20.,200,20.)
    nominalBinningTT["pt_g1"] = frange(20., 200, 10.)
    nominalBinningTT["phi1"] = frange(0., 3.1, 0.1)
    nominalBinningTT["phi2"] = frange(0., 3.1, 0.1)
    nominalBinningTT["n_jets"] = frange(0., 7, 7)
    nominalBinningTT["ht"] = frange(0., 1000., 50)
    nominalBinningTT["met"] = frange(0., 450., 50)
    #nominalBinningTT["met"]  = [0,25,50,75,100,150,250,450]
    nominalBinningTT["m_ll"] = frange(40., 300., 20)

    biggerBinningTT = CR_tt.binnings.copy()
    biggerBinningTT["eta1"] = frange(0., 2.4, 0.2)
    biggerBinningTT["eta2"] = frange(0., 2.4, 0.2)
    biggerBinningTT["pt1"] = frange(25., 150, 25)
    biggerBinningTT["pt2"] = frange(25., 150, 25)
    #biggerBinningTT["pt_g1"]  = frange(20.,220,40)
    biggerBinningTT["pt_g1"] = frange(20., 220, 20)
    biggerBinningTT["phi1"] = frange(0., 3.2, 0.2)
    biggerBinningTT["phi2"] = frange(0., 3.2, 0.2)
    biggerBinningTT["n_jets"] = frange(0., 8, 4)
    biggerBinningTT["ht"] = frange(0., 1000., 100)
    biggerBinningTT["met"] = frange(0., 500., 100)
    biggerBinningTT["m_ll"] = frange(50., 300., 50.)

    # drawTT(binningToUse_=nominalBinningTT, addName="")
    # drawTT(binningToUse_=smallerBinningTT, addName="_small")
    # drawTT(binningToUse_=biggerBinningTT, addName="_big")

    # drawTT(binningToUse_=nominalBinningTT,addName="_JESu")
    # drawTT(binningToUse_=nominalBinningTT,addName="_JESd")
    # drawTT(binningToUse_=nominalBinningTT,addName="_JERu")
    # drawTT(binningToUse_=nominalBinningTT,addName="_JERd")

    # drawTT(binningToUse_=nominalBinningTT,addName="_lepSFu")
    # drawTT(binningToUse_=nominalBinningTT,addName="_lepSFd")
    # drawTT(binningToUse_=nominalBinningTT,addName="_photonSFu")
    # drawTT(binningToUse_=nominalBinningTT,addName="_photonSFd")
    # drawTT(binningToUse_=nominalBinningTT,addName="_ISRu")
    # drawTT(binningToUse_=nominalBinningTT,addName="_ISRd")
    # drawTT(binningToUse_=nominalBinningTT,addName="_EWKu")
    # drawTT(binningToUse_=nominalBinningTT,addName="_EWKd")
    # drawTT(binningToUse_=nominalBinningTT,addName="_PUu")
    # drawTT(binningToUse_=nominalBinningTT,addName="_PUd")

    # drawTT(binningToUse_=nominalBinningTT,addName="_pdf")

    # drawDY()

    smallerBinningDY = CR_DY.binnings.copy()
    smallerBinningDY["eta1"] = frange(0., 2.4, 0.05)
    smallerBinningDY["eta2"] = frange(0., 2.4, 0.05)
    smallerBinningDY["pt1"] = frange(25., 150, 5.)
    smallerBinningDY["pt2"] = frange(20., 150, 5.)
    smallerBinningDY["pt_g1"] = frange(20, 200, 10)
    smallerBinningDY["phi1"] = frange(0., 3.1, 0.05)
    smallerBinningDY["phi2"] = frange(0., 3.1, 0.05)
    smallerBinningDY["n_jets"] = frange(0., 7, 7)
    smallerBinningDY["ht"] = frange(0., 1000., 25)
    smallerBinningDY["met"] = frange(0, 100, 5)
    smallerBinningDY["m_ll"] = frange(80., 100., 1)

    nominalBinningDY = CR_DY.binnings.copy()
    nominalBinningDY["eta1"] = frange(0., 2.4, 0.1)
    nominalBinningDY["eta2"] = frange(0., 2.4, 0.1)
    nominalBinningDY["pt1"] = frange(20., 150, 10.)
    nominalBinningDY["pt2"] = frange(20., 150, 10.)
    nominalBinningDY["pt_g1"] = frange(20., 200, 20.)
    nominalBinningDY["phi1"] = frange(0., 3.1, 0.1)
    nominalBinningDY["phi2"] = frange(0., 3.1, 0.1)
    nominalBinningDY["n_jets"] = frange(0., 7, 7)
    nominalBinningDY["ht"] = frange(0., 1000., 50)
    nominalBinningDY["met"] = frange(0, 100, 10)
    #nominalBinningDY["met"]  = [0,25,50,75,100,150,250,450]
    nominalBinningDY["m_ll"] = frange(80., 100., 2)

    biggerBinningDY = CR_DY.binnings.copy()
    biggerBinningDY["eta1"] = frange(0., 2.4, 0.2)
    biggerBinningDY["eta2"] = frange(0., 2.4, 0.2)
    biggerBinningDY["pt1"] = frange(25., 150, 25)
    biggerBinningDY["pt2"] = frange(25., 150, 25)
    biggerBinningDY["pt_g1"] = frange(20., 220, 40)
    biggerBinningDY["phi1"] = frange(0., 3.2, 0.2)
    biggerBinningDY["phi2"] = frange(0., 3.2, 0.2)
    biggerBinningDY["n_jets"] = frange(0., 8, 4)
    biggerBinningDY["ht"] = frange(0., 1000., 100)
    biggerBinningDY["met"] = frange(0, 100, 20)
    biggerBinningTT["m_ll"] = frange(80., 100., 4)

    # drawDY(binningToUse_=nominalBinningDY, addName="")
    # drawDY(binningToUse_=smallerBinningDY, addName="_small")
    # drawDY(binningToUse_=biggerBinningDY, addName="_big")

    # drawDY(binningToUse_=nominalBinningDY,addName="_JESu")
    # drawDY(binningToUse_=nominalBinningDY,addName="_JESd")
    # drawDY(binningToUse_=nominalBinningDY,addName="_JERu")
    # drawDY(binningToUse_=nominalBinningDY,addName="_JERd")

    # drawDY(binningToUse_=nominalBinningDY,addName="_lepSFu")
    # drawDY(binningToUse_=nominalBinningDY,addName="_lepSFd")
    # drawDY(binningToUse_=nominalBinningDY,addName="_photonSFu")
    # drawDY(binningToUse_=nominalBinningDY,addName="_photonSFd")

    # drawZZ()

    smallerBinningZZ = CR_ZZ.binnings.copy()
    smallerBinningZZ["eta1"] = frange(0., 2.4, 0.05)
    smallerBinningZZ["eta2"] = frange(0., 2.4, 0.05)
    smallerBinningZZ["pt1"] = frange(25., 150, 5.)
    smallerBinningZZ["pt2"] = frange(20., 100, 5.)
    smallerBinningZZ["phi1"] = frange(0., 3.1, 0.05)
    smallerBinningZZ["phi2"] = frange(0., 3.1, 0.05)
    smallerBinningZZ["n_jets"] = frange(0., 7, 7)
    smallerBinningZZ["ht"] = frange(0., 1000., 25)
    smallerBinningZZ["met"] = frange(0., 80., 5)
    smallerBinningZZ["m_ll"] = frange(84., 98., 1)
    smallerBinningZZ["m_ll2"] = frange(60., 120., 2.)

    nominalBinningZZ = CR_ZZ.binnings.copy()
    nominalBinningZZ["eta1"] = frange(0., 2.4, 0.1)
    nominalBinningZZ["eta2"] = frange(0., 2.4, 0.1)
    nominalBinningZZ["pt1"] = frange(20., 150, 10.)
    nominalBinningZZ["pt2"] = frange(20., 100, 10.)
    nominalBinningZZ["phi1"] = frange(0., 3.1, 0.1)
    nominalBinningZZ["phi2"] = frange(0., 3.1, 0.1)
    nominalBinningZZ["n_jets"] = frange(0., 7, 7)
    nominalBinningZZ["ht"] = frange(0., 1000., 50)
    nominalBinningZZ["met"] = frange(0., 80., 10.)
    #nominalBinningZZ["met"]  = [0,25,50,75,100,150,250,450]
    nominalBinningZZ["m_ll"] = frange(84., 98., 2.)
    nominalBinningZZ["m_ll2"] = frange(60., 120., 5.)

    biggerBinningZZ = CR_ZZ.binnings.copy()
    biggerBinningZZ["eta1"] = frange(0., 2.4, 0.2)
    biggerBinningZZ["eta2"] = frange(0., 2.4, 0.2)
    biggerBinningZZ["pt1"] = frange(25., 150, 25)
    biggerBinningZZ["pt2"] = frange(25., 100, 25)
    biggerBinningZZ["phi1"] = frange(0., 3.2, 0.2)
    biggerBinningZZ["phi2"] = frange(0., 3.2, 0.2)
    biggerBinningZZ["n_jets"] = frange(0., 8, 4)
    biggerBinningZZ["ht"] = frange(0., 1000., 100)
    biggerBinningZZ["met"] = frange(0., 80., 20)
    biggerBinningZZ["m_ll"] = frange(82., 98., 4.)
    biggerBinningZZ["m_ll2"] = frange(60., 120., 10.)

    drawZZ(binningToUse_=nominalBinningZZ, addName="")
    drawZZ(binningToUse_=smallerBinningZZ, addName="_small")
    drawZZ(binningToUse_=biggerBinningZZ, addName="_big")

    # drawZZ(binningToUse_=nominalBinningZZ,addName="_JESu")
    # drawZZ(binningToUse_=nominalBinningZZ,addName="_JESd")
    # drawZZ(binningToUse_=nominalBinningZZ,addName="_JERu")
    # drawZZ(binningToUse_=nominalBinningZZ,addName="_JERd")

    # drawZZ(binningToUse_=nominalBinningZZ,addName="_lepSFu")
    # drawZZ(binningToUse_=nominalBinningZZ,addName="_lepSFd")
    # drawZZ(binningToUse_=nominalBinningZZ,addName="_photonSFu")
    # drawZZ(binningToUse_=nominalBinningZZ,addName="_photonSFd")

    # drawTT(binningToUse_=nominalBinningTT,addName="")
    # drawTT(binningToUse_=smallerBinningTT,addName="_small")
    # drawTT(binningToUse_=biggerBinningTT,addName="_big")
    # drawWZ()

    smallerBinningWZ = CR_WZ.binnings.copy()
    smallerBinningWZ["eta1"] = frange(0., 2.4, 0.05)
    smallerBinningWZ["eta2"] = frange(0., 2.4, 0.05)
    smallerBinningWZ["pt1"] = frange(25., 150, 5.)
    smallerBinningWZ["pt2"] = frange(20., 150, 5.)
    smallerBinningWZ["pt_g1"] = frange(20, 200, 10)
    smallerBinningWZ["phi1"] = frange(0., 3.1, 0.05)
    smallerBinningWZ["phi2"] = frange(0., 3.1, 0.05)
    smallerBinningWZ["n_jets"] = frange(0., 7, 7)
    smallerBinningWZ["ht"] = frange(0., 1000., 25)
    smallerBinningWZ["met"] = frange(50., 450., 25)
    smallerBinningWZ["m_ll"] = frange(76., 106., 1)

    nominalBinningWZ = CR_WZ.binnings.copy()
    nominalBinningWZ["eta1"] = frange(0., 2.4, 0.1)
    nominalBinningWZ["eta2"] = frange(0., 2.4, 0.1)
    nominalBinningWZ["pt1"] = frange(20., 150, 10.)
    nominalBinningWZ["pt2"] = frange(20., 150, 10.)
    nominalBinningWZ["phi1"] = frange(0., 3.1, 0.1)
    nominalBinningWZ["phi2"] = frange(0., 3.1, 0.1)
    nominalBinningWZ["n_jets"] = frange(0., 7, 7)
    nominalBinningWZ["ht"] = frange(0., 1000., 50)
    nominalBinningWZ["met"] = frange(50., 450., 50)
    #nominalBinningWZ["met"]  = [0,25,50,75,100,150,250,450]
    nominalBinningWZ["m_ll"] = frange(76., 106., 2.)

    biggerBinningWZ = CR_WZ.binnings.copy()
    biggerBinningWZ["eta1"] = frange(0., 2.4, 0.2)
    biggerBinningWZ["eta2"] = frange(0., 2.4, 0.2)
    biggerBinningWZ["pt1"] = frange(25., 150, 25)
    biggerBinningWZ["pt2"] = frange(25., 150, 25)
    biggerBinningWZ["pt_g1"] = frange(20., 220, 40)
    biggerBinningWZ["phi1"] = frange(0., 3.2, 0.2)
    biggerBinningWZ["phi2"] = frange(0., 3.2, 0.2)
    biggerBinningWZ["n_jets"] = frange(0., 8, 4)
    biggerBinningWZ["ht"] = frange(0., 1000., 100)
    biggerBinningWZ["met"] = frange(0., 500., 100)
    biggerBinningWZ["m_ll"] = frange(75., 110., 5.)

    drawWZ(binningToUse_=nominalBinningWZ, addName="")
    drawWZ(binningToUse_=smallerBinningWZ, addName="_small")
    drawWZ(binningToUse_=biggerBinningWZ, addName="_big")

    # drawWZ(binningToUse_=nominalBinningWZ,addName="_JESu")
    # drawWZ(binningToUse_=nominalBinningWZ,addName="_JESd")
    # drawWZ(binningToUse_=nominalBinningWZ,addName="_JERu")
    # drawWZ(binningToUse_=nominalBinningWZ,addName="_JERd")

    # drawWZ(binningToUse_=nominalBinningWZ,addName="_lepSFu")
    # drawWZ(binningToUse_=nominalBinningWZ,addName="_lepSFd")
    # drawWZ(binningToUse_=nominalBinningWZ,addName="_photonSFu")
    # drawWZ(binningToUse_=nominalBinningWZ,addName="_photonSFd")


if __name__ == "__main__":
    main()
