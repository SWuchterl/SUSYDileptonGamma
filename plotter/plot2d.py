import ROOT
from ROOT import *
from ROOT.TColor import *
from array import array
from include import *
import numpy as np

import style


from dataMC import frange


def drawSameHistogram(sampleNames, name, bkg=[], additional=[], binning=None, binningName="", scaleToData=False, xTitle=None, yTitle=None):
    style.defaultStyle()
    style.style2d()
    gStyle.SetCanvasDefH(200)
    gStyle.SetCanvasDefW(600)
    can = ROOT.TCanvas("", "", 600, 300)
    # can = ROOT.TCanvas()

    scale = 1.

    folder = (name.split("/"))[0] + "/" + (name.split("/"))[1]

    print name
    arWeightsToUse = ["nISR", "topPt", "ewk"]
    zgHist = aux.stdHist2dWithWeights(zgamma, name, arWeightsToUse, binning)
    ttgHist = aux.stdHist2dWithWeights(ttgamma, name, arWeightsToUse, binning)
    zzHist = aux.stdHist2dWithWeights(zz, name, arWeightsToUse, binning)
    wwgHist = aux.stdHist2dWithWeights(wwgamma, name, arWeightsToUse, binning)
    wzgHist = aux.stdHist2dWithWeights(wzgamma, name, arWeightsToUse, binning)
    dyHist = aux.stdHist2dWithWeights(DYjetsNLO, name, arWeightsToUse, binning)
    wjetsHist = aux.stdHist2dWithWeights(wjets, name, arWeightsToUse, binning)
    ttHist = aux.stdHist2dWithWeights(tt, name, arWeightsToUse, binning)
    singletopHist = aux.stdHist2dWithWeights(
        singletop, name, arWeightsToUse, binning)

    wzHist = aux.stdHist2dWithWeights(wz, name, arWeightsToUse, binning)
    wwHist = aux.stdHist2dWithWeights(ww, name, arWeightsToUse, binning)
    zz4lHist = aux.stdHist2dWithWeights(zz4l, name, arWeightsToUse, binning)
    wgHist = aux.stdHist2dWithWeights(wgamma, name, arWeightsToUse, binning)

    #print folder+"/weight_topPt"

    # histsToScale=[zgHist,ttgHist,ttg080Hist,ttg80Hist,zzHist,wwgHist,wzgHist,dyHist,wjetsHist,ttHist,tt080Hist,tt80Hist,singletopHist,wzHist,wwHist,zz4lHist,wgHist]
    #
    # histsToScale = [zgHist, ttgHist, ttg080Hist, ttg80Hist, wwgHist, wzgHist, dyHist,
    #                 wjetsHist, ttHist, tt080Hist, tt80Hist, singletopHist, wzHist, wwHist, zz4lHist, wgHist]
    histsToScale = [zgHist, ttgHist, wwgHist, wzgHist, dyHist,
                    wjetsHist, ttHist,  singletopHist, wzHist, wwHist, zz4lHist, wgHist]

    signal3 = aux.stdHist2d(t5bbbbzg_1500_400, name, binning)
    # signal1 = aux.stdHist2d(tching_400, name, binning)
    signal1 = aux.stdHist2d(tching_600, name, binning)
    signal2 = aux.stdHist2d(gmsb_240_230, name, binning)
    # signal4 = aux.stdHist2d(gmsb_290_205, name, binning)
    signal4 = aux.stdHist2d(gmsb_415_355, name, binning)
    signal5 = aux.stdHist2d(t5bbbbzg_1500_1400, name, binning)

    bkg = aux.addHists(zgHist, ttgHist, ttHist, wwgHist, wzgHist, dyHist,
                       wjetsHist, singletopHist, wzHist, wwHist, zz4lHist, wgHist)

    bkg.SetDirectory(0)

    # print bkg.Integral()
    # binx1 = bkg.GetXaxis().FindBin(150)
    # binx2 = bkg.GetXaxis().FindBin(500)
    # biny1 = bkg.GetYaxis().FindBin(100)
    # biny2 = bkg.GetYaxis().FindBin(500)
    # print bkg.Integral(binx1, binx2, biny1, biny2)

    # print bkg

    # h=signal1.Clone()

    # for i in range(signal1.GetNbinsX()):
    # for j in range(signal1.GetNbinsY()):
    # s=signal1.GetBinContent(i,j)
    # b=bkg.GetBinContent(i,j)

    # if (s+b)==0:
    # value=0
    # else:
    ##print s,b
    # value=s/np.sqrt(s+b)
    # h.SetBinContent(i,j,value)

    #print bkg.GetNcells(),bkg.GetNbinsX(),bkg.GetNbinsY()
    #print signal1.GetNcells(),signal1.GetNbinsX (),signal1.GetNbinsY()

    #bkg.SetYTitle( aux.getYAxisTitle( bkg ) )
    if yTitle:
        bkg.SetYTitle(yTitle)
    if xTitle:
        bkg.SetXTitle(xTitle)
    #signal1.SetYTitle( aux.getYAxisTitle( signal1 ) )
    if yTitle:
        signal1.SetYTitle(yTitle)
        signal2.SetYTitle(yTitle)
        signal3.SetYTitle(yTitle)
        signal4.SetYTitle(yTitle)
        signal5.SetYTitle(yTitle)
    if xTitle:
        signal1.SetXTitle(xTitle)
        signal2.SetXTitle(xTitle)
        signal3.SetXTitle(xTitle)
        signal4.SetXTitle(xTitle)
        signal5.SetXTitle(xTitle)

    # h.SetZTitle("Events/Bin")
    bkg.SetZTitle("Events/Bin")
    signal1.SetZTitle("Events/Bin")
    signal2.SetZTitle("Events/Bin")
    signal3.SetZTitle("Events/Bin")
    signal4.SetZTitle("Events/Bin")
    signal5.SetZTitle("Events/Bin")
    # bkg.SetZTitle("#frac{s}#sqrt{{s+b}}")
    # bkg.SetZTitle("s/#sqrt{s+b}")
    # signal1.SetZTitle("s/#sqrt{s+b}")

    # h.GetZaxis().SetRangeUser(0.,0.3)
    #bkg.Draw("same colz")

    #dataHist = None
    # for d in additional:
    #h = d.getHist2D( name )
    # if not h: continue
    # if not h.Integral(): continue
    # if (binning):
    #h = aux.rebin2d( h, binning[0],binning[1] )
    #aux.appendFlowBin2d( h )

    sigS1 = signal1.Clone()
    sigS2 = signal1.Clone()
    sigS3 = signal1.Clone()
    sigS4 = signal1.Clone()
    sigS5 = signal1.Clone()

    for i in range(signal1.GetNbinsX()):
        for j in range(signal1.GetNbinsY()):
            s1 = signal1.GetBinContent(i, j)
            s2 = signal2.GetBinContent(i, j)
            s3 = signal3.GetBinContent(i, j)
            s4 = signal4.GetBinContent(i, j)
            s5 = signal5.GetBinContent(i, j)
            b = bkg.GetBinContent(i, j)

            if b < 0:
                b = 0

            # value1 = 0. if (s1 + b) == 0 else s1 / np.sqrt(s1 + b)
            # value2 = 0. if (s2 + b) == 0 else s2 / np.sqrt(s2 + b)
            # value3 = 0. if (s3 + b) == 0 else s3 / np.sqrt(s3 + b)
            # value4 = 0. if (s4 + b) == 0 else s4 / np.sqrt(s4 + b)
            # value5 = 0. if (s5 + b) == 0 else s5 / np.sqrt(s5 + b)
            # value1 = 0. if (b) == 0 else ROOT.TMath.Sqrt(2 *
            #                                              (s1 + b) * ROOT.TMath.Log(1 + s1 / b) - 2 * s1)
            # value2 = 0. if (b) == 0 else ROOT.TMath.Sqrt(2 *
            #                                              (s2 + b) * ROOT.TMath.Log(1 + s2 / b) - 2 * s2)
            # value3 = 0. if (b) == 0 else ROOT.TMath.Sqrt(2 *
            #                                              (s3 + b) * ROOT.TMath.Log(1 + s3 / b) - 2 * s3)
            # value4 = 0. if (b) == 0 else ROOT.TMath.Sqrt(2 *
            #                                              (s4 + b) * ROOT.TMath.Log(1 + s4 / b) - 2 * s4)
            # value5 = 0. if (b) == 0 else ROOT.TMath.Sqrt(2 *
            #                                              (s5 + b) * ROOT.TMath.Log(1 + s5 / b) - 2 * s5)
            value1 = s1
            value2 = s2
            value3 = s3
            value4 = s4
            value5 = s5

            sigS1.SetBinContent(i, j, value1)
            sigS2.SetBinContent(i, j, value2)
            sigS3.SetBinContent(i, j, value3)
            sigS4.SetBinContent(i, j, value4)
            sigS5.SetBinContent(i, j, value5)

    sigS1.SetZTitle("Events / bin")
    sigS2.SetZTitle("Events / bin")
    sigS3.SetZTitle("Events / bin")
    sigS4.SetZTitle("Events / bin")
    sigS5.SetZTitle("Events / bin")
    # sigS1.SetZTitle("s/#sqrt{s+b}")
    # sigS2.SetZTitle("s/#sqrt{s+b}")
    # sigS3.SetZTitle("s/#sqrt{s+b}")
    # sigS4.SetZTitle("s/#sqrt{s+b}")
    # sigS5.SetZTitle("s/#sqrt{s+b}")

    bkg.GetZaxis().SetRangeUser(0., bkg.GetMaximum() * 1.1)

    signal1.GetZaxis().SetRangeUser(0., signal1.GetMaximum() * 1.1)
    signal2.GetZaxis().SetRangeUser(0., signal2.GetMaximum() * 1.1)
    signal3.GetZaxis().SetRangeUser(0., signal3.GetMaximum() * 1.1)
    signal4.GetZaxis().SetRangeUser(0., signal4.GetMaximum() * 1.1)
    signal5.GetZaxis().SetRangeUser(0., signal5.GetMaximum() * 1.1)

    sigS1.GetZaxis().SetRangeUser(0., sigS1.GetMaximum() * 1.1)
    sigS2.GetZaxis().SetRangeUser(0., sigS2.GetMaximum() * 1.1)
    sigS3.GetZaxis().SetRangeUser(0., sigS3.GetMaximum() * 1.1)
    sigS4.GetZaxis().SetRangeUser(0., sigS4.GetMaximum() * 1.1)
    sigS5.GetZaxis().SetRangeUser(0., sigS5.GetMaximum() * 1.1)

    #print sigS1.GetMaximum()
    #print sigS2.GetMaximum()
    #print sigS3.GetMaximum()
    #print sigS4.GetMaximum()
    #print sigS5.GetMaximum()

    info = ""
    l = aux.Label(info="#scale[0.7]{%s}" % info, sim=not((dataDoubleMuon in additional)or(
        dataDoubleEG in additional)or(dataHt in additional)or(dataDoubleSF in additional)or(dataMuonEG in additional)))

    if binningName:
        binningName = "_" + binningName
    name = name.replace("/", "__")
    saveName = "sameHistograms_{}_{}{}".format(sampleNames, name, binningName)
    #aux.save("DataMC_"+saveName+"_bkg_",folder="plots_2d/" ,log=False)

    style.style2d()
    s = style.style2d()
    s.SetPadLeftMargin(0.18)
    l = ROOT.TLatex(
        0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
    l2 = ROOT.TLatex(0.31, .88, "#scale[0.76]{#font[52]{Simulation}}")
    l3 = ROOT.TLatex(0.51, .88, "#scale[0.46]{#font[12]{Background}}")
    cr = ROOT.TLatex(0.22, 0.5, "#scale[0.4]{DY/Z(#gamma) Control Region}")
    vr = ROOT.TLatex(0.45, 0.15, "#scale[0.4]{Validation Region}")
    sr = ROOT.TLatex(0.51, 0.8, "#scale[0.4]{Signal Region}")
    l.SetNDC()
    l2.SetNDC()
    l3.SetNDC()
    cr.SetNDC()
    cr.SetTextAngle(90)
    vr.SetNDC()
    sr.SetNDC()
    lineMet100 = TLine(100, 0, 100, 1000)
    lineMet150 = TLine(150, 100, 150, 1000)
    # lineMet200 = TLine(200, 0, 200, 500)
    lineMt2100 = TLine(150, 100, 1000, 100)
    lineMet100.SetLineStyle(9)
    lineMet150.SetLineStyle(9)
    # lineMet200.SetLineStyle(10)
    lineMt2100.SetLineStyle(9)

    lum = ROOT.TLatex(.62, .95, "%.1f fb^{-1} (%s TeV)" %
                      (aux.intLumi / 1000., aux.Label.cmsEnergy))
    lum.SetNDC()

    # can = ROOT.TCanvas()
    can = ROOT.TCanvas("", "", 600, 500)

    can.SetLogz()
    # bkg.SetMaximum(100)
    # bkg.SetMinimum(0.01)
    bkg.SetMinimum(0.0001)
    bkg.Draw("colz")
    lum.Draw()
    l.Draw()
    l2.Draw()
    l3.Draw()
    cr.Draw()
    vr.Draw()
    sr.Draw()
    lineMet100.Draw()
    lineMet150.Draw()
    lineMt2100.Draw()
    aux.save("DataMC_" + saveName + "_bkg_", folder="plots_2d/", log=False)
    # can.SaveAs('plots_2d/root'+'_'+saveName+'.root')

    # can = ROOT.TCanvas()
    can = ROOT.TCanvas("", "", 600, 500)

    can.SetLogz()
    signal1.Draw("same colz")
    # aux.save("DataMC_" + saveName + "_tching400_",
    lum.Draw()
    l.Draw()
    l2.Draw()
    l3 = ROOT.TLatex(
        0.51, .88, "#scale[0.46]{#font[12]{TChiZG m(NLSP)=600 GeV}}")
    l3.SetNDC()
    l3.Draw()
    cr.Draw()
    vr.Draw()
    sr.Draw()
    lineMet100.Draw()
    lineMet150.Draw()
    lineMt2100.Draw()
    aux.save("DataMC_" + saveName + "_tching600_",
             folder="plots_2d/", log=False)
    can = ROOT.TCanvas()
    can.SetLogz()
    signal2.Draw("same colz")
    lum.Draw()
    l.Draw()
    l2.Draw()
    l3 = ROOT.TLatex(0.51, .88, "#scale[0.46]{#font[12]{GMSB 240 230}}")
    l3.SetNDC()
    l3.Draw()
    cr.Draw()
    vr.Draw()
    sr.Draw()
    lineMet100.Draw()
    lineMet150.Draw()
    lineMt2100.Draw()
    aux.save("DataMC_" + saveName + "_gmsb_240_230_",
             folder="plots_2d/", log=False)
    # can = ROOT.TCanvas()
    can = ROOT.TCanvas("", "", 600, 500)
    can.SetLogz()
    signal3.Draw("same colz")
    lum.Draw()
    l.Draw()
    l2.Draw()
    l3 = ROOT.TLatex(0.51, .88, "#scale[0.46]{#font[12]{T5bbbbZG 1500 400}}")
    l3.SetNDC()
    l3.Draw()
    cr.Draw()
    vr.Draw()
    sr.Draw()
    lineMet100.Draw()
    lineMet150.Draw()
    lineMt2100.Draw()
    aux.save("DataMC_" + saveName + "_t5bbbbzg_1500_400_",
             folder="plots_2d/", log=False)
    # can = ROOT.TCanvas()
    can = ROOT.TCanvas("", "", 600, 500)
    can.SetLogz()
    signal4.Draw("same colz")
    lum.Draw()
    l.Draw()
    l2.Draw()
    l3 = ROOT.TLatex(0.51, .88, "#scale[0.46]{#font[12]{GMSB 290 205}}")
    l3.SetNDC()
    l3.Draw()
    cr.Draw()
    vr.Draw()
    sr.Draw()
    lineMet100.Draw()
    lineMet150.Draw()
    lineMt2100.Draw()
    aux.save("DataMC_" + saveName + "_gmsb_290_205_",
             folder="plots_2d/", log=False)
    # can = ROOT.TCanvas()
    can = ROOT.TCanvas("", "", 600, 500)
    can.SetLogz()
    signal5.Draw("same colz")
    lum.Draw()
    l.Draw()
    l2.Draw()
    l3 = ROOT.TLatex(0.51, .88, "#scale[0.46]{#font[12]{T5bbbbZG 1500 1400}}")
    l3.SetNDC()
    l3.Draw()
    cr.Draw()
    vr.Draw()
    sr.Draw()
    lineMet100.Draw()
    lineMet150.Draw()
    lineMt2100.Draw()
    aux.save("DataMC_" + saveName + "_t5bbbbzg_1500_1400_",
             folder="plots_2d/", log=False)

    # can = ROOT.TCanvas()
    can = ROOT.TCanvas("", "", 600, 500)
    can.SetLogz()
    sigS1.Draw("same colz")
    lum.Draw()
    l.Draw()
    l2.Draw()
    l3 = ROOT.TLatex(0.51, .88, "#scale[0.46]{#font[12]{TChiZG 600}}")
    l3.SetNDC()
    l3.Draw()
    cr.Draw()
    vr.Draw()
    sr.Draw()
    lineMet100.Draw()
    lineMet150.Draw()
    lineMt2100.Draw()
    # aux.save("DataMC_" + saveName + "_SIG_tching400_",
    aux.save("DataMC_" + saveName + "_SIG_tching600_",
             folder="plots_2d/", log=False)
    # can = ROOT.TCanvas()
    can = ROOT.TCanvas("", "", 600, 500)
    can.SetLogz()
    sigS2.Draw("same colz")
    lum.Draw()
    l.Draw()
    l2.Draw()
    l3 = ROOT.TLatex(0.51, .88, "#scale[0.46]{#font[12]{GMSB 240 230}}")
    l3.SetNDC()
    l3.Draw()
    cr.Draw()
    vr.Draw()
    sr.Draw()
    lineMet100.Draw()
    lineMet150.Draw()
    lineMt2100.Draw()
    aux.save("DataMC_" + saveName + "_SIG_gmsb_240_230_",
             folder="plots_2d/", log=False)
    # can = ROOT.TCanvas()
    can = ROOT.TCanvas("", "", 600, 500)
    can.SetLogz()
    sigS3.Draw("same colz")
    lum.Draw()
    l.Draw()
    l2.Draw()
    l3 = ROOT.TLatex(
        0.51, .88, "#scale[0.46]{#font[12]{T5bbbbzg m(#tilde{g})=1500 m(#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}})=400}}")
    l3.SetNDC()
    l3.Draw()
    cr.Draw()
    vr.Draw()
    sr.Draw()
    lineMet100.Draw()
    lineMet150.Draw()
    lineMt2100.Draw()
    aux.save("DataMC_" + saveName + "_SIG_t5bbbbzg_1500_400_",
             folder="plots_2d/", log=False)
    # can = ROOT.TCanvas()
    can = ROOT.TCanvas("", "", 600, 500)
    can.SetLogz()
    sigS4.Draw("same colz")
    lum.Draw()
    l.Draw()
    l2.Draw()
    # l3 = ROOT.TLatex(0.61, .88, "#scale[0.46]{#font[12]{GMSB 290 205}}")
    l3 = ROOT.TLatex(
        0.51, .88, "#scale[0.46]{#font[12]{GMSB m(#tilde{W})=415 m(#tilde{B})=355}}")
    l3.SetNDC()
    l3.Draw()
    cr.Draw()
    vr.Draw()
    sr.Draw()
    lineMet100.Draw()
    lineMet150.Draw()
    lineMt2100.Draw()
    # aux.save("DataMC_" + saveName + "_SIG_gmsb_290_205_",
    #          folder="plots_2d/", log=False)
    aux.save("DataMC_" + saveName + "_SIG_gmsb_415_355_",
             folder="plots_2d/", log=False)
    # can = ROOT.TCanvas()
    can = ROOT.TCanvas("", "", 600, 500)
    can.SetLogz()
    sigS5.Draw("same colz")
    lum.Draw()
    l.Draw()
    l2.Draw()
    l3 = ROOT.TLatex(0.51, .88, "#scale[0.46]{#font[12]{T5bbbbZG 1500 1400}}")
    l3.SetNDC()
    l3.Draw()
    cr.Draw()
    vr.Draw()
    sr.Draw()
    lineMet100.Draw()
    lineMet150.Draw()
    lineMt2100.Draw()
    aux.save("DataMC_" + saveName + "_SIG_t5bbbbzg_1500_1400_",
             folder="plots_2d/", log=False)


def main2():
    # bkgs=[DYjetsNLO,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma,singletop,zz4l,wz,ww]
    bkgs = [ttgamma]
    # bkgs=[zgamma]
    # variables=["Fakes","met","stmet","st","mtll","mtllg","mtl1met","mtl2met","mtllmet","mtllgmet","mtgmet","m_llg","m_ll","n_vtx","DeltaPhiLLMet","DeltaEtaLLMet","DeltaRLLMet","pt1","pt_g1","ht","eta1","phi1","n_jets","n_photons","deltaPhiLLG","deltaEtaLLG","deltaRLLG","deltaPhiLL","deltaEtaLL","deltaRLL","zpt","nElectrons","nMuons","deltaR1_g1","deltaR2_g1","DeltaPhiGMet","n_bjets"]
    # variables=["MetZpt","MetMllg","MetDeltaPhiLL"]
    # variables=["MetDeltaPhiLL"]
    variables = ["MetMt2"]
    groups = ["onZ"]
    binnings_ = {

        'MetZpt': [frange(0, 500, 50), frange(0, 500, 50)],

        # 'MetMt2': [frange(0, 300, 5), frange(0, 300, 5)],
        'MetMt2': [frange(0, 1000, 10), frange(0, 1000, 10)],

        'MetMllg': [[150, 200, 300], [100, 200, 500]],
        'MetDeltaPhiLL': [[150, 200, 300], [0, 1., 2., 3.5]],
    }
    labels_ = {
        'MetZpt': ["p_{T}^{miss} (GeV)", "p_{T}^{Z} (GeV) "],
        'MetMllg': ["p_{T}^{miss} (GeV)", "m_{ll#gamma}(GeV) "],
        'MetDeltaPhiLL': ["p_{T}^{miss} (GeV)", "#Delta#Phi_{ll}"],
        'MetMt2': ["p_{T}^{miss} (GeV)", "M_{T2} (GeV)"]
    }
    for group in groups:
        for variable in variables:
            #drawSameHistogram("LL+signal",group+"/LL/"+variable, bkgs, additional=[t5bbbbzg_1500_1400,t5bbbbzg_1500_400,t5bbbbzg_1500_600,tching_600,tching_400,gmsb_240_230,gmsb_290_205],binning=binnings_[variable],xTitle=labels_[variable][0],yTitle=labels_[variable])
            drawSameHistogram("LL+signal", group + "/LL/" + variable, bkgs, additional=[],
                              binning=binnings_[variable], xTitle=labels_[variable][0], yTitle=labels_[variable][1])
            #drawSameHistogram("LL+signal",group+"/LL/"+variable, bkgs, additional=[],xTitle=labels_[variable][0],yTitle=labels_[variable][1])


if __name__ == "__main__":
    main2()
