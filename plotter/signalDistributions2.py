import ROOT
from ROOT import *
from array import array
from include import *
import numpy as np
import pickle as pkl


def getDatacardUncertFromHist(h, b):
    c = h.GetBinContent(b)
    return 1. + max(0, h.GetBinError(b) / c if c else 0)


def finalDistributionSignalHist(name, dirSet, dirDir):
    #style.divideByBinWidth = True
    style.divideByBinWidth = False

    #nBins = range(0,200,10)+[200, 250, 300, 350, 450, 600, 800]
    #nBins = range(0,200,25)+[200, 300,500]
    #nBins = [0,100,150]+[150,190, 230,500]
    #nBins = [0,100,150,190, 230,500]
    #nBinsOld = [0,25,50,75,100,150,190,230,500]
    #nBinsDataOld = [0,25,50,75,100]
    #nBins = [0,25,50,75,100,150,190,230,500]
    #nBinsData = [0,25,50,75,100]
    # according to Knuts Tool ->250
    #nBins = [0,25,50,75,100,150,220,350]
    #nBins = [0,25,50,75,100,150,200,350]
    #nBins = [0,30,50,80,100,150,200,350]
    nBins = [150, 200, 350]  # yes!
    #nBins = [100,150,200,350]
    #nBins = [150,220,350]
    #nBins = [0,25,50,75,100,150,250,500,700]
    # nBinsData = [0,25,50,75,100]????
    #nBinsData = [0,25,50,75,100,150]
    # nBinsData = [0, 30, 50, 80, 100, 150]
    nBinsData = [0, 30, 50, 80, 100, 150, 200, 350]  # unblind

    # direct stuff

    print dirDir

    dataHist = aux.stdHistWithoutNGen(dataLL, dirDir + "/met_RAW", nBins)

    dirHist = aux.stdHist(dirSet, dirDir + "/met", nBins)

    #dataHist = aux.stdHist(dataDoubleSF,dirDir+"/met",nBins)
    #dataHist = aux.stdHist(dataLL,dirDir+"/met",nBins)
    #dataHist = aux.rebin(dataHist,nBinsData)
    # dataHist.SetYTitle(aux.getYAxisTitle(dataHist))

    style.additionalPoissonUncertainty = False
    aux.drawOpt(dirHist, "data")
    aux.drawOpt(dataHist, "data")

    #aux.drawOpt(dataHist, "data")
    #gjetHist, gjetSyst, info = gjetPrediction(dirHist, preSet, zg+wg+ttg+wjets+ttjets_nlo+znunu, "met", nBins, weight, name+"_divByBinWidth" if style.divideByBinWidth else name)
    # gjetHist.SetLineColor(rwth.myLightBlue)
    # gjetHist.GetXaxis().SetTitle("#it{p}_{T}^{miss} (GeV)")

    #eHist = aux.stdHist(preSetElectron, preDirElectron+"/met", nBins)
    # eHist.GetXaxis().SetTitle("#it{p}_{T}^{miss} (GeV)")
    #eHist.Scale( 0.0267 if dirSet == data else 0.0154 )
    # eHist.SetLineColor(rwth.myYellow)
    #eSyst = aux.getSysHisto(eHist, 0.3)

    zgHist = aux.stdHistWithoutNGenWithWeights(
        zgamma, dirDir + "/met", ["nISR", "topPt", "ewk"], nBins)
    ttgHist = aux.stdHistWithoutNGenWithWeights(
        ttgamma, dirDir + "/met", ["nISR", "topPt", "ewk"], nBins)
    ttg080Hist = aux.stdHistWithoutNGenWithWeights(ttgamma, dirDir.replace(
        "sig", "sig080") + "/met", ["nISR", "topPt", "ewk"], nBins)
    ttg80Hist = aux.stdHistWithoutNGenWithWeights(ttgamma, dirDir.replace(
        "sig", "sig80") + "/met", ["nISR", "topPt", "ewk"], nBins)
    zzHist = aux.stdHistWithoutNGenWithWeights(
        zz, dirDir + "/met", ["nISR", "topPt", "ewk"], nBins)
    wwgHist = aux.stdHistWithoutNGenWithWeights(
        wwgamma, dirDir + "/met", ["nISR", "topPt", "ewk"], nBins)
    wzgHist = aux.stdHistWithoutNGenWithWeights(
        wzgamma, dirDir + "/met", ["nISR", "topPt", "ewk"], nBins)
    dyHist = aux.stdHistWithoutNGenWithWeights(
        DYjetsNLO, dirDir + "/met", ["nISR", "topPt", "ewk"], nBins)
    wjetsHist = aux.stdHistWithoutNGenWithWeights(
        wjets, dirDir + "/met", ["nISR", "topPt", "ewk"], nBins)
    ttHist = aux.stdHistWithoutNGenWithWeights(
        tt, dirDir + "/met", ["nISR", "topPt", "ewk"], nBins)
    tt080Hist = aux.stdHistWithoutNGenWithWeights(tt, dirDir.replace(
        "sig", "sig080") + "/met", ["nISR", "topPt", "ewk"], nBins)
    tt80Hist = aux.stdHistWithoutNGenWithWeights(tt, dirDir.replace(
        "sig", "sig80") + "/met", ["nISR", "topPt", "ewk"], nBins)
    singletopHist = aux.stdHistWithoutNGenWithWeights(
        singletop, dirDir + "/met", ["nISR", "topPt", "ewk"], nBins)
    wzHist = aux.stdHistWithoutNGenWithWeights(
        wz, dirDir + "/met", ["nISR", "topPt", "ewk"], nBins)
    wwHist = aux.stdHistWithoutNGenWithWeights(
        ww, dirDir + "/met", ["nISR", "topPt", "ewk"], nBins)
    zz4lHist = aux.stdHistWithoutNGenWithWeights(
        zz4l, dirDir + "/met", ["nISR", "topPt", "ewk"], nBins)
    wgHist = aux.stdHistWithoutNGenWithWeights(
        wgamma, dirDir + "/met", ["nISR", "topPt", "ewk"], nBins)

    #print ttg080Hist.Integral(),ttg80Hist.Integral(),ttHist.Integral()

    # addedHistos:
    final_dyHist = aux.addHists(dyHist, zgHist)
    final_ttHist = aux.addHists(ttHist, ttgHist)
    final_tt080Hist = aux.addHists(tt080Hist, ttg080Hist)
    final_tt80Hist = aux.addHists(tt80Hist, ttg80Hist)
    final_zzHist = aux.addHists(zzHist, zz4lHist)
    final_wzHist = aux.addHists(wzHist)
    final_otherHist = aux.addHists(
        wwgHist, wzgHist, wjetsHist, singletopHist, wwHist, wgHist)

    # Scaling

    # histsToScale=[zgHist,ttgHist,zzHist,wwgHist,wzgHist,dyHist,wjetsHist,ttHist,singletopHist,wzHist,wwHist,zz4lHist,wgHist]
    histsToScale = [final_dyHist, final_ttHist,
                    final_zzHist, final_wzHist, final_otherHist]

    pklZZ = pkl.load(open("plots_CR_zz/factors/CRZZ.pkl", "rb"))
    pklDY = pkl.load(open("plots_CR_dy/factors/CRDY.pkl", "rb"))
    pklTT = pkl.load(open("plots_CR_tt/factors/CRTT.pkl", "rb"))
    pklWZ = pkl.load(open("plots_CR_wz/factors/CRWZ.pkl", "rb"))
    sfZZ, sfZZErr = pklZZ["LL"]["m_ll"]
    sfDY, sfDYErr = pklDY["LL"]["eta1"]
    sfTT = pklTT["EM"]["eta1"][0]
    sfTTErr = pklTT["EM"]["eta1"][1]
    sfTT080Err = 0.04
    sfTT80Err = 0.4
    sfWZ, sfWZErr = pklWZ["LL"]["eta1"]

    zzHist.Scale(sfZZ)
    zz4lHist.Scale(sfZZ)
    dyHist.Scale(sfDY)
    zgHist.Scale(sfDY)
    wzHist.Scale(sfWZ)
    ttHist.Scale(sfTT)
    tt80Hist.Scale(sfTT)
    tt080Hist.Scale(sfTT)
    ttgHist.Scale(sfTT)
    ttg80Hist.Scale(sfTT)
    ttg080Hist.Scale(sfTT)

    final_tt080Hist.Scale(sfTT)
    final_tt80Hist.Scale(sfTT)
    final_ttHist.Scale(sfTT)
    final_zzHist.Scale(sfZZ)
    final_wzHist.Scale(sfWZ)
    final_dyHist.Scale(sfDY)

    #print "INT", ttgHist.Integral(),ttHist.Integral()

    dirHist = aux.addHists(*histsToScale)

    zgHist.SetLineColor(ROOT.kGreen - 3)
    ttgHist.SetLineColor(ROOT.kRed + 1)
    ttg80Hist.SetLineColor(ROOT.kRed + 1)
    ttg080Hist.SetLineColor(ROOT.kRed + 1)
    zzHist.SetLineColor(ROOT.kYellow)
    wwgHist.SetLineColor(ROOT.kCyan + 3)
    wzgHist.SetLineColor(ROOT.kBlue + 3)
    dyHist.SetLineColor(ROOT.kGreen + 3)
    wjetsHist.SetLineColor(ROOT.kBlue - 9)
    ttHist.SetLineColor(ROOT.kOrange + 8)
    tt080Hist.SetLineColor(ROOT.kOrange + 8)
    tt80Hist.SetLineColor(ROOT.kOrange + 8)
    singletopHist.SetLineColor(ROOT.kOrange + 4)
    # wzHist.SetLineColor(ROOT.kAzure-5)
    wzHist.SetLineColor(ROOT.kAzure - 6)
    wwHist.SetLineColor(ROOT.kCyan - 3)
    zz4lHist.SetLineColor(ROOT.kOrange - 2)
    wgHist.SetLineColor(ROOT.kRed + 3)

    # final_otherHist.SetLineColor(ROOT.kCyan+3)
    # final_ttHist.SetLineColor(ROOT.kOrange+8)
    # final_ttHist.SetLineColor(ROOT.kRed+1)
    final_ttHist.SetLineColor(ROOT.kOrange + 7)
    final_zzHist.SetLineColor(ROOT.kOrange - 2)
    # final_wzHist.SetLineColor(ROOT.kAzure-5)
    final_wzHist.SetLineColor(ROOT.kAzure - 6)
    # final_wzHist.SetLineColor(ROOT.kBlue+3)
    # final_dyHist.SetLineColor(ROOT.kGreen-3)
    final_dyHist.SetLineColor(ROOT.kGreen + 3)
    # final_dyHist.SetLineColor(ROOT.kGreen+2)
    final_otherHist.SetLineColor(ROOT.kGray + 2)

    # mcSystUncert = 0.04 # SF, lumi, trigger
    # mcSystUncert = 0.2 # SF, lumi, trigger
    mcSystUncert = 0.026  # lumi

    # ttgSFUnc =

    import pickle
    systPkl = pickle.load(open("systBKG/unc.pkl", "rb"))
    import uncertHelper as hlp
    # zgScaleSyst=hlp.getScaleUncertHisto(zgHist,[0.,systPkl["zg"]["scale"]["binMC_1"],systPkl["zg"]["scale"]["binMC_2"]])
    # zgScaleSyst=hlp.getScaleUncertHisto(zgHist,[0.,systPkl["zg"]["scale"]["binMC_1"],systPkl["zg"]["scale"]["binMC_2"]],zgHist.Integral()/sum(zgamma.ngens))
    # ttgScaleSyst=hlp.getScaleUncertHisto(ttgHist,[0.,systPkl["ttg"]["scale"]["binMC_1"],systPkl["ttg"]["scale"]["binMC_2"]])
    # ttg080ScaleSyst=hlp.getScaleUncertHisto(ttg080Hist,[0.,systPkl["ttg"]["scale"]["binMC_1"],systPkl["ttg"]["scale"]["binMC_2"]])
    # ttg80ScaleSyst=hlp.getScaleUncertHisto(ttg80Hist,[0.,systPkl["ttg"]["scale"]["binMC_1"],systPkl["ttg"]["scale"]["binMC_2"]])
    # zzScaleSyst=hlp.getScaleUncertHisto(zzHist,[0.,systPkl["zz"]["scale"]["binMC_1"],systPkl["zz"]["scale"]["binMC_2"]])
    # wwgScaleSyst=hlp.getScaleUncertHisto(wwgHist,[0.,systPkl["wwg"]["scale"]["binMC_1"],systPkl["wwg"]["scale"]["binMC_2"],0.1],wwgHist.Integral()/sum(wwgamma.ngens))
    # wzgScaleSyst=hlp.getScaleUncertHisto(wzgHist,[0.,systPkl["wzg"]["scale"]["binMC_1"],systPkl["wzg"]["scale"]["binMC_2"],0.1],wzgHist.Integral()/sum(wzgamma.ngens))
    # dyScaleSyst=hlp.getScaleUncertHisto(dyHist,[0.,systPkl["dy"]["scale"]["binMC_1"],systPkl["dy"]["scale"]["binMC_2"]],dyHist.Integral()/DYjetsNLO.ngens[0])
    # wjetsScaleSyst=hlp.getScaleUncertHisto(wjetsHist,[0.,systPkl["wjets"]["scale"]["binMC_1"],systPkl["wjets"]["scale"]["binMC_2"],0.1],wjetsHist.Integral()/sum(wjets.ngens))
    # ttScaleSyst=hlp.getScaleUncertHisto(ttHist,[0.,systPkl["tt"]["scale"]["binMC_1"],systPkl["tt"]["scale"]["binMC_2"]])
    # ttScaleSyst=hlp.getScaleUncertHisto(tt80Hist,[0.,systPkl["tt"]["scale"]["binMC_1"],systPkl["tt"]["scale"]["binMC_2"]])
    # ttScaleSyst=hlp.getScaleUncertHisto(tt080Hist,[0.,systPkl["tt"]["scale"]["binMC_1"],systPkl["tt"]["scale"]["binMC_2"]])
    # singletopScaleSyst=hlp.getScaleUncertHisto(singletopHist,[0.,systPkl["singletop"]["scale"]["binMC_1"],systPkl["singletop"]["scale"]["binMC_2"],0.1],singletopHist.Integral()/sum(singletop.ngens))
    wzScaleSyst = hlp.getScaleUncertHisto(
        wzHist, [0., systPkl["wz"]["scale"]["binMC_1"], systPkl["wz"]["scale"]["binMC_2"]])
    # wwScaleSyst=hlp.getScaleUncertHisto(wwHist,[0.,systPkl["ww"]["scale"]["binMC_1"],systPkl["ww"]["scale"]["binMC_2"],0.1],wwHist.Integral()/sum(ww.ngens))
    # zz4lScaleSyst=hlp.getScaleUncertHisto(zz4lHist,[0.,systPkl["zz4l"]["scale"]["binMC_1"],systPkl["zz4l"]["scale"]["binMC_2"],])
    # wgScaleSyst=hlp.getScaleUncertHisto(wgHist,[0.,systPkl["wgamma"]["scale"]["binMC_1"],systPkl["wgamma"]["scale"]["binMC_2"],0.1],wgHist.Integral()/sum(wgamma.ngens))
    otherScaleSyst = hlp.getScaleUncertHisto(final_otherHist, [
                                             0., systPkl["other"]["scale"]["binMC_1"], systPkl["other"]["scale"]["binMC_2"]], final_otherHist.Integral() / sum(otherBKG.ngens))
    final_ttScaleSyst = hlp.getScaleUncertHisto(final_ttHist, [
                                                0., systPkl["tt+ttg"]["scale"]["binMC_1"], systPkl["tt+ttg"]["scale"]["binMC_2"]])
    final_tt080ScaleSyst = hlp.getScaleUncertHisto(final_tt080Hist, [
                                                   0., systPkl["tt+ttg"]["scale"]["binMC_1"], systPkl["tt+ttg"]["scale"]["binMC_2"]])
    final_tt80ScaleSyst = hlp.getScaleUncertHisto(final_tt80Hist, [
                                                  0., systPkl["tt+ttg"]["scale"]["binMC_1"], systPkl["tt+ttg"]["scale"]["binMC_2"]])
    final_dyScaleSyst = hlp.getScaleUncertHisto(final_dyHist, [
                                                0., systPkl["dy+zg"]["scale"]["binMC_1"], systPkl["dy+zg"]["scale"]["binMC_2"]], final_dyHist.Integral() / (sum(DYjetsNLO.ngens) + sum(zgamma.ngens)))
    final_zzScaleSyst = hlp.getScaleUncertHisto(final_zzHist, [
                                                0., systPkl["zz+zz4l"]["scale"]["binMC_1"], systPkl["zz+zz4l"]["scale"]["binMC_2"]])

    # zgPDFSyst=hlp.getPDFUncertHisto(zgHist,[0.,systPkl["zg"]["pdf"]["binMC_1"],systPkl["zg"]["pdf"]["binMC_2"]],zgHist.Integral()/sum(zgamma.ngens))
    # ttgPDFSyst=hlp.getPDFUncertHisto(ttgHist,[0.,systPkl["ttg"]["pdf"]["binMC_1"],systPkl["ttg"]["pdf"]["binMC_2"]])
    # ttg080PDFSyst=hlp.getPDFUncertHisto(ttg080Hist,[0.,systPkl["ttg"]["pdf"]["binMC_1"],systPkl["ttg"]["pdf"]["binMC_2"]])
    # ttg80PDFSyst=hlp.getPDFUncertHisto(ttg80Hist,[0.,systPkl["ttg"]["pdf"]["binMC_1"],systPkl["ttg"]["pdf"]["binMC_2"]])
    # zzPDFSyst=hlp.getPDFUncertHisto(zzHist,[0.,systPkl["zz"]["pdf"]["binMC_1"],systPkl["zz"]["pdf"]["binMC_2"]])
    # wwgPDFSyst=hlp.getPDFUncertHisto(wwgHist,[0.,systPkl["wwg"]["pdf"]["binMC_1"],systPkl["wwg"]["pdf"]["binMC_2"],0.1],wwgHist.Integral()/sum(wwgamma.ngens))
    # wzgPDFSyst=hlp.getPDFUncertHisto(wzgHist,[0.,systPkl["wzg"]["pdf"]["binMC_1"],systPkl["wzg"]["pdf"]["binMC_2"],0.1],wzgHist.Integral()/sum(wzgamma.ngens))
    # dyPDFSyst=hlp.getPDFUncertHisto(dyHist,[0.,systPkl["dy"]["pdf"]["binMC_1"],systPkl["dy"]["pdf"]["binMC_2"]],dyHist.Integral()/DYjetsNLO.ngens[0])
    # wjetsPDFSyst=hlp.getPDFUncertHisto(wjetsHist,[0.,systPkl["wjets"]["pdf"]["binMC_1"],systPkl["wjets"]["pdf"]["binMC_2"],0.1],wjetsHist.Integral()/sum(wjets.ngens))
    # ttPDFSyst=hlp.getPDFUncertHisto(ttHist,[0.,systPkl["tt"]["pdf"]["binMC_1"],systPkl["tt"]["pdf"]["binMC_2"]])
    # ttPDFSyst=hlp.getPDFUncertHisto(tt80Hist,[0.,systPkl["tt"]["pdf"]["binMC_1"],systPkl["tt"]["pdf"]["binMC_2"]])
    # ttPDFSyst=hlp.getPDFUncertHisto(tt080Hist,[0.,systPkl["tt"]["pdf"]["binMC_1"],systPkl["tt"]["pdf"]["binMC_2"]])
    # singletopPDFSyst=hlp.getPDFUncertHisto(singletopHist,[0.,systPkl["singletop"]["pdf"]["binMC_1"],systPkl["singletop"]["pdf"]["binMC_2"],0.1],singletopHist.Integral()/sum(singletop.ngens))
    wzPDFSyst = hlp.getPDFUncertHisto(
        wzHist, [0., systPkl["wz"]["pdf"]["binMC_1"], systPkl["wz"]["pdf"]["binMC_2"]])
    # wwPDFSyst=hlp.getPDFUncertHisto(wwHist,[0.,systPkl["ww"]["pdf"]["binMC_1"],systPkl["ww"]["pdf"]["binMC_2"],0.1],wwHist.Integral()/sum(ww.ngens))
    # zz4lPDFSyst=hlp.getPDFUncertHisto(zz4lHist,[0.,systPkl["zz4l"]["pdf"]["binMC_1"],systPkl["zz4l"]["pdf"]["binMC_2"]])
    # wgPDFSyst=hlp.getPDFUncertHisto(wgHist,[0.,systPkl["wgamma"]["pdf"]["binMC_1"],systPkl["wgamma"]["pdf"]["binMC_2"],0.1],wgHist.Integral()/sum(wgamma.ngens))
    otherPDFSyst = hlp.getPDFUncertHisto(final_otherHist, [
                                         0., systPkl["other"]["pdf"]["binMC_1"], systPkl["other"]["pdf"]["binMC_2"]], final_otherHist.Integral() / sum(otherBKG.ngens))
    final_ttPDFSyst = hlp.getScaleUncertHisto(final_ttHist, [
                                              0., systPkl["tt+ttg"]["pdf"]["binMC_1"], systPkl["tt+ttg"]["pdf"]["binMC_2"]])
    final_tt080PDFSyst = hlp.getScaleUncertHisto(final_tt080Hist, [
                                                 0., systPkl["tt+ttg"]["pdf"]["binMC_1"], systPkl["tt+ttg"]["pdf"]["binMC_2"]])
    final_tt80PDFSyst = hlp.getScaleUncertHisto(final_tt80Hist, [
                                                0., systPkl["tt+ttg"]["pdf"]["binMC_1"], systPkl["tt+ttg"]["pdf"]["binMC_2"]])
    final_dyPDFSyst = hlp.getScaleUncertHisto(final_dyHist, [
                                              0., systPkl["dy+zg"]["pdf"]["binMC_1"], systPkl["dy+zg"]["pdf"]["binMC_2"]], final_dyHist.Integral() / (sum(DYjetsNLO.ngens) + sum(zgamma.ngens)))
    final_zzPDFSyst = hlp.getScaleUncertHisto(final_zzHist, [
                                              0., systPkl["zz+zz4l"]["pdf"]["binMC_1"], systPkl["zz+zz4l"]["pdf"]["binMC_2"]])

    # zgPUSyst=hlp.getScaleUncertHisto(zgHist,[0.,systPkl["zg"]["PU"]["binMC_1"],systPkl["zg"]["PU"]["binMC_2"]],zgHist.Integral()/sum(zgamma.ngens))
    # ttgPUSyst=hlp.getScaleUncertHisto(ttgHist,[0.,systPkl["ttg"]["PU"]["binMC_1"],systPkl["ttg"]["PU"]["binMC_2"]])
    # ttg080PUSyst=hlp.getScaleUncertHisto(ttg080Hist,[0.,systPkl["ttg"]["PU"]["binMC_1"],systPkl["ttg"]["PU"]["binMC_2"]])
    # ttg80PUSyst=hlp.getScaleUncertHisto(ttg80Hist,[0.,systPkl["ttg"]["PU"]["binMC_1"],systPkl["ttg"]["PU"]["binMC_2"]])
    # zzPUSyst=hlp.getScaleUncertHisto(zzHist,[0.,systPkl["zz"]["PU"]["binMC_1"],systPkl["zz"]["PU"]["binMC_2"]])
    # wwgPUSyst=hlp.getScaleUncertHisto(wwgHist,[0.,systPkl["wwg"]["PU"]["binMC_1"],systPkl["wwg"]["PU"]["binMC_2"],0.1],wwgHist.Integral()/sum(wwgamma.ngens))
    # wzgPUSyst=hlp.getScaleUncertHisto(wzgHist,[0.,systPkl["wzg"]["PU"]["binMC_1"],systPkl["wzg"]["PU"]["binMC_2"],0.1],wzgHist.Integral()/sum(wzgamma.ngens))
    # dyPUSyst=hlp.getScaleUncertHisto(dyHist,[0.,systPkl["dy"]["PU"]["binMC_1"],systPkl["dy"]["PU"]["binMC_2"]],dyHist.Integral()/DYjetsNLO.ngens[0])
    # wjetsPUSyst=hlp.getScaleUncertHisto(wjetsHist,[0.,systPkl["wjets"]["PU"]["binMC_1"],systPkl["wjets"]["PU"]["binMC_2"],0.1],wjetsHist.Integral()/sum(wjets.ngens))
    # ttPUSyst=hlp.getScaleUncertHisto(ttHist,[0.,systPkl["tt"]["PU"]["binMC_1"],systPkl["tt"]["PU"]["binMC_2"]])
    # tt80PUSyst=hlp.getScaleUncertHisto(tt80Hist,[0.,systPkl["tt"]["PU"]["binMC_1"],systPkl["tt"]["PU"]["binMC_2"]])
    # tt080PUSyst=hlp.getScaleUncertHisto(tt080Hist,[0.,systPkl["tt"]["PU"]["binMC_1"],systPkl["tt"]["PU"]["binMC_2"]])
    # singletopPUSyst=hlp.getScaleUncertHisto(zgHist,[0.,systPkl["singletop"]["PU"]["binMC_1"],systPkl["singletop"]["PU"]["binMC_2"],0.1],singletopHist.Integral()/sum(singletop.ngens))
    wzPUSyst = hlp.getScaleUncertHisto(
        wzHist, [0., systPkl["wz"]["PU"]["binMC_1"], systPkl["wz"]["PU"]["binMC_2"]])
    # wwPUSyst=hlp.getScaleUncertHisto(wwHist,[0.,systPkl["ww"]["PU"]["binMC_1"],systPkl["ww"]["PU"]["binMC_2"],0.1],wwHist.Integral()/sum(ww.ngens))
    # zz4lPUSyst=hlp.getScaleUncertHisto(zz4lHist,[0.,systPkl["zz4l"]["PU"]["binMC_1"],systPkl["zz4l"]["PU"]["binMC_2"]])
    # wgPUSyst=hlp.getScaleUncertHisto(wgHist,[0.,systPkl["wgamma"]["PU"]["binMC_1"],systPkl["wgamma"]["PU"]["binMC_2"],0.1],wgHist.Integral()/sum(wgamma.ngens))
    otherPUSyst = hlp.getScaleUncertHisto(final_otherHist, [
                                          0., systPkl["other"]["PU"]["binMC_1"], systPkl["other"]["PU"]["binMC_2"]], final_otherHist.Integral() / sum(otherBKG.ngens))
    final_ttPUSyst = hlp.getScaleUncertHisto(final_ttHist, [
                                             0., systPkl["tt+ttg"]["PU"]["binMC_1"], systPkl["tt+ttg"]["PU"]["binMC_2"]])
    final_tt080PUSyst = hlp.getScaleUncertHisto(final_tt080Hist, [
                                                0., systPkl["tt+ttg"]["PU"]["binMC_1"], systPkl["tt+ttg"]["PU"]["binMC_2"]])
    final_tt80PUSyst = hlp.getScaleUncertHisto(final_tt80Hist, [
                                               0., systPkl["tt+ttg"]["PU"]["binMC_1"], systPkl["tt+ttg"]["PU"]["binMC_2"]])
    final_dyPUSyst = hlp.getScaleUncertHisto(final_dyHist, [
                                             0., systPkl["dy+zg"]["PU"]["binMC_1"], systPkl["dy+zg"]["PU"]["binMC_2"]], final_dyHist.Integral() / (sum(DYjetsNLO.ngens) + sum(zgamma.ngens)))
    final_zzPUSyst = hlp.getScaleUncertHisto(final_zzHist, [
                                             0., systPkl["zz+zz4l"]["PU"]["binMC_1"], systPkl["zz+zz4l"]["PU"]["binMC_2"]])

    # zgJESSyst=hlp.getScaleUncertHisto(zgHist,[0.,systPkl["zg"]["JES"]["binMC_1"],systPkl["zg"]["JES"]["binMC_2"]],zgHist.Integral()/sum(zgamma.ngens))
    # ttgJESSyst=hlp.getScaleUncertHisto(ttgHist,[0.,systPkl["ttg"]["JES"]["binMC_1"],systPkl["ttg"]["JES"]["binMC_2"]])
    # ttg80JESSyst=hlp.getScaleUncertHisto(ttg80Hist,[0.,systPkl["ttg"]["JES"]["binMC_1"],systPkl["ttg"]["JES"]["binMC_2"]])
    # ttg080JESSyst=hlp.getScaleUncertHisto(ttg080Hist,[0.,systPkl["ttg"]["JES"]["binMC_1"],systPkl["ttg"]["JES"]["binMC_2"]])
    # zzJESSyst=hlp.getScaleUncertHisto(zzHist,[0.,systPkl["zz"]["JES"]["binMC_1"],systPkl["zz"]["JES"]["binMC_2"]])
    # wwgJESSyst=hlp.getScaleUncertHisto(wwgHist,[0.,systPkl["wwg"]["JES"]["binMC_1"],systPkl["wwg"]["JES"]["binMC_2"],0.1],wwgHist.Integral()/sum(wwgamma.ngens))
    # wzgJESSyst=hlp.getScaleUncertHisto(wzgHist,[0.,systPkl["wzg"]["JES"]["binMC_1"],systPkl["wzg"]["JES"]["binMC_2"],0.1],wzgHist.Integral()/sum(wzgamma.ngens))
    # dyJESSyst=hlp.getScaleUncertHisto(dyHist,[0.,systPkl["dy"]["JES"]["binMC_1"],systPkl["dy"]["JES"]["binMC_2"]],dyHist.Integral()/DYjetsNLO.ngens[0])
    # wjetsJESSyst=hlp.getScaleUncertHisto(wjetsHist,[0.,systPkl["wjets"]["JES"]["binMC_1"],systPkl["wjets"]["JES"]["binMC_2"],0.1],wjetsHist.Integral()/sum(wjets.ngens))
    # ttJESSyst=hlp.getScaleUncertHisto(ttHist,[0.,systPkl["tt"]["JES"]["binMC_1"],systPkl["tt"]["JES"]["binMC_2"]])
    # tt80JESSyst=hlp.getScaleUncertHisto(tt080Hist,[0.,systPkl["tt"]["JES"]["binMC_1"],systPkl["tt"]["JES"]["binMC_2"]])
    # tt080JESSyst=hlp.getScaleUncertHisto(tt80Hist,[0.,systPkl["tt"]["JES"]["binMC_1"],systPkl["tt"]["JES"]["binMC_2"]])
    # singletopJESSyst=hlp.getScaleUncertHisto(singletopHist,[0.,systPkl["singletop"]["JES"]["binMC_1"],systPkl["singletop"]["JES"]["binMC_2"],0.1],singletopHist.Integral()/sum(singletop.ngens))
    wzJESSyst = hlp.getScaleUncertHisto(
        wzHist, [0., systPkl["wz"]["JES"]["binMC_1"], systPkl["wz"]["JES"]["binMC_2"]])
    # wwJESSyst=hlp.getScaleUncertHisto(wwHist,[0.,systPkl["ww"]["JES"]["binMC_1"],systPkl["ww"]["JES"]["binMC_2"],0.1],wwHist.Integral()/sum(ww.ngens))
    # zz4lJESSyst=hlp.getScaleUncertHisto(zz4lHist,[0.,systPkl["zz4l"]["JES"]["binMC_1"],systPkl["zz4l"]["JES"]["binMC_2"]])
    # wgJESSyst=hlp.getScaleUncertHisto(wgHist,[0.,systPkl["wgamma"]["JES"]["binMC_1"],systPkl["wgamma"]["JES"]["binMC_2"],0.1],wgHist.Integral()/sum(wgamma.ngens))
    otherJESSyst = hlp.getScaleUncertHisto(final_otherHist, [
                                           0., systPkl["other"]["JES"]["binMC_1"], systPkl["other"]["JES"]["binMC_2"]], final_otherHist.Integral() / sum(otherBKG.ngens))
    final_ttJESSyst = hlp.getScaleUncertHisto(final_ttHist, [
                                              0., systPkl["tt+ttg"]["JES"]["binMC_1"], systPkl["tt+ttg"]["JES"]["binMC_2"]])
    final_tt080JESSyst = hlp.getScaleUncertHisto(final_tt080Hist, [
                                                 0., systPkl["tt+ttg"]["JES"]["binMC_1"], systPkl["tt+ttg"]["JES"]["binMC_2"]])
    final_tt80JESSyst = hlp.getScaleUncertHisto(final_tt80Hist, [
                                                0., systPkl["tt+ttg"]["JES"]["binMC_1"], systPkl["tt+ttg"]["JES"]["binMC_2"]])
    final_dyJESSyst = hlp.getScaleUncertHisto(final_dyHist, [
                                              0., systPkl["dy+zg"]["JES"]["binMC_1"], systPkl["dy+zg"]["JES"]["binMC_2"]], final_dyHist.Integral() / (sum(DYjetsNLO.ngens) + sum(zgamma.ngens)))
    final_zzJESSyst = hlp.getScaleUncertHisto(final_zzHist, [
                                              0., systPkl["zz+zz4l"]["JES"]["binMC_1"], systPkl["zz+zz4l"]["JES"]["binMC_2"]])

    # zgJERSyst=hlp.getScaleUncertHisto(zgHist,[0.,systPkl["zg"]["JER"]["binMC_1"],systPkl["zg"]["JER"]["binMC_2"]],zgHist.Integral()/sum(zgamma.ngens))
    # ttgJERSyst=hlp.getScaleUncertHisto(ttgHist,[0.,systPkl["ttg"]["JER"]["binMC_1"],systPkl["ttg"]["JER"]["binMC_2"]])
    # zzJERSyst=hlp.getScaleUncertHisto(zzHist,[0.,systPkl["zz"]["JER"]["binMC_1"],systPkl["zz"]["JER"]["binMC_2"]])
    # wwgJERSyst=hlp.getScaleUncertHisto(wwgHist,[0.,systPkl["wwg"]["JER"]["binMC_1"],systPkl["wwg"]["JER"]["binMC_2"],0.1],wwgHist.Integral()/sum(wwgamma.ngens))
    # wzgJERSyst=hlp.getScaleUncertHisto(wzgHist,[0.,systPkl["wzg"]["JER"]["binMC_1"],systPkl["wzg"]["JER"]["binMC_2"],0.1],wzgHist.Integral()/sum(wzgamma.ngens))
    # dyJERSyst=hlp.getScaleUncertHisto(dyHist,[0.,systPkl["dy"]["JER"]["binMC_1"],systPkl["dy"]["JER"]["binMC_2"]],dyHist.Integral()/DYjetsNLO.ngens[0])
    # wjetsJERSyst=hlp.getScaleUncertHisto(wjetsHist,[0.,systPkl["wjets"]["JER"]["binMC_1"],systPkl["wjets"]["JER"]["binMC_2"],0.1],wjetsHist.Integral()/sum(wjets.ngens))
    # ttJERSyst=hlp.getScaleUncertHisto(ttHist,[0.,systPkl["tt"]["JER"]["binMC_1"],systPkl["tt"]["JER"]["binMC_2"]])
    # singletopJERSyst=hlp.getScaleUncertHisto(singletopHist,[0.,systPkl["singletop"]["JER"]["binMC_1"],systPkl["singletop"]["JER"]["binMC_2"],0.1],singletopHist.Integral()/sum(singletop.ngens))
    wzJERSyst = hlp.getScaleUncertHisto(
        wzHist, [0., systPkl["wz"]["JER"]["binMC_1"], systPkl["wz"]["JER"]["binMC_2"]])
    # wwJERSyst=hlp.getScaleUncertHisto(wwHist,[0.,systPkl["ww"]["JER"]["binMC_1"],systPkl["ww"]["JER"]["binMC_2"],0.1],wwHist.Integral()/sum(ww.ngens))
    # zz4lJERSyst=hlp.getScaleUncertHisto(zz4lHist,[0.,systPkl["zz4l"]["JER"]["binMC_1"],systPkl["zz4l"]["JER"]["binMC_2"]])
    # wgJERSyst=hlp.getScaleUncertHisto(wgHist,[0.,systPkl["wgamma"]["JER"]["binMC_1"],systPkl["wgamma"]["JER"]["binMC_2"],0.1],wgHist.Integral()/sum(wgamma.ngens))
    otherJERSyst = hlp.getScaleUncertHisto(final_otherHist, [
                                           0., systPkl["other"]["JER"]["binMC_1"], systPkl["other"]["JER"]["binMC_2"]], final_otherHist.Integral() / sum(otherBKG.ngens))
    final_ttJERSyst = hlp.getScaleUncertHisto(final_ttHist, [
                                              0., systPkl["tt+ttg"]["JER"]["binMC_1"], systPkl["tt+ttg"]["JER"]["binMC_2"]])
    final_tt080JERSyst = hlp.getScaleUncertHisto(final_tt080Hist, [
                                                 0., systPkl["tt+ttg"]["JER"]["binMC_1"], systPkl["tt+ttg"]["JER"]["binMC_2"]])
    final_tt80JERSyst = hlp.getScaleUncertHisto(final_tt80Hist, [
                                                0., systPkl["tt+ttg"]["JER"]["binMC_1"], systPkl["tt+ttg"]["JER"]["binMC_2"]])
    final_dyJERSyst = hlp.getScaleUncertHisto(final_dyHist, [
                                              0., systPkl["dy+zg"]["JER"]["binMC_1"], systPkl["dy+zg"]["JER"]["binMC_2"]], final_dyHist.Integral() / (sum(DYjetsNLO.ngens) + sum(zgamma.ngens)))
    final_zzJERSyst = hlp.getScaleUncertHisto(final_zzHist, [
                                              0., systPkl["zz+zz4l"]["JER"]["binMC_1"], systPkl["zz+zz4l"]["JER"]["binMC_2"]])

    # zglepSFSyst=hlp.getScaleUncertHisto(zgHist,[0.,systPkl["zg"]["lepSF"]["binMC_1"],systPkl["zg"]["lepSF"]["binMC_2"]],zgHist.Integral()/sum(zgamma.ngens))
    # ttglepSFSyst=hlp.getScaleUncertHisto(ttgHist,[0.,systPkl["ttg"]["lepSF"]["binMC_1"],systPkl["ttg"]["lepSF"]["binMC_2"]])
    # zzlepSFSyst=hlp.getScaleUncertHisto(wwgHist,[0.,systPkl["zz"]["lepSF"]["binMC_1"],systPkl["zz"]["lepSF"]["binMC_2"]])
    # wwglepSFSyst=hlp.getScaleUncertHisto(wwHist,[0.,systPkl["wwg"]["lepSF"]["binMC_1"],systPkl["wwg"]["lepSF"]["binMC_2"],0.1],wwgHist.Integral()/sum(wwgamma.ngens))
    # wzglepSFSyst=hlp.getScaleUncertHisto(wzgHist,[0.,systPkl["wzg"]["lepSF"]["binMC_1"],systPkl["wzg"]["lepSF"]["binMC_2"],0.1],wzgHist.Integral()/sum(wzgamma.ngens))
    # dylepSFSyst=hlp.getScaleUncertHisto(dyHist,[0.,systPkl["dy"]["lepSF"]["binMC_1"],systPkl["dy"]["lepSF"]["binMC_2"]],dyHist.Integral()/DYjetsNLO.ngens[0])
    # wjetslepSFSyst=hlp.getScaleUncertHisto(wjetsHist,[0.,systPkl["wjets"]["lepSF"]["binMC_1"],systPkl["wjets"]["lepSF"]["binMC_2"],0.1],wjetsHist.Integral()/sum(wjets.ngens))
    # ttlepSFSyst=hlp.getScaleUncertHisto(ttHist,[0.,systPkl["tt"]["lepSF"]["binMC_1"],systPkl["tt"]["lepSF"]["binMC_2"]])
    # singletoplepSFSyst=hlp.getScaleUncertHisto(singletopHist,[0.,systPkl["singletop"]["lepSF"]["binMC_1"],systPkl["singletop"]["lepSF"]["binMC_2"],0.1],singletopHist.Integral()/sum(singletop.ngens))
    wzlepSFSyst = hlp.getScaleUncertHisto(
        wzHist, [0., systPkl["wz"]["lepSF"]["binMC_1"], systPkl["wz"]["lepSF"]["binMC_2"]])
    # wwlepSFSyst=hlp.getScaleUncertHisto(wwHist,[0.,systPkl["ww"]["lepSF"]["binMC_1"],systPkl["ww"]["lepSF"]["binMC_2"],0.1],wwHist.Integral()/sum(ww.ngens))
    # zz4llepSFSyst=hlp.getScaleUncertHisto(zz4lHist,[0.,systPkl["zz4l"]["lepSF"]["binMC_1"],systPkl["zz4l"]["lepSF"]["binMC_2"]])
    # wglepSFSyst=hlp.getScaleUncertHisto(wgHist,[0.,systPkl["wgamma"]["lepSF"]["binMC_1"],systPkl["wgamma"]["lepSF"]["binMC_2"],0.1],wgHist.Integral()/sum(wgamma.ngens))
    otherlepSFSyst = hlp.getScaleUncertHisto(final_otherHist, [
                                             0., systPkl["other"]["lepSF"]["binMC_1"], systPkl["other"]["lepSF"]["binMC_2"]], final_otherHist.Integral() / sum(otherBKG.ngens))
    final_ttlepSFSyst = hlp.getScaleUncertHisto(final_ttHist, [
                                                0., systPkl["tt+ttg"]["lepSF"]["binMC_1"], systPkl["tt+ttg"]["lepSF"]["binMC_2"]])
    final_tt080lepSFSyst = hlp.getScaleUncertHisto(final_tt080Hist, [
                                                   0., systPkl["tt+ttg"]["lepSF"]["binMC_1"], systPkl["tt+ttg"]["lepSF"]["binMC_2"]])
    final_tt80lepSFSyst = hlp.getScaleUncertHisto(final_tt80Hist, [
                                                  0., systPkl["tt+ttg"]["lepSF"]["binMC_1"], systPkl["tt+ttg"]["lepSF"]["binMC_2"]])
    final_dylepSFSyst = hlp.getScaleUncertHisto(final_dyHist, [
                                                0., systPkl["dy+zg"]["lepSF"]["binMC_1"], systPkl["dy+zg"]["lepSF"]["binMC_2"]], final_dyHist.Integral() / (sum(DYjetsNLO.ngens) + sum(zgamma.ngens)))
    final_zzlepSFSyst = hlp.getScaleUncertHisto(final_zzHist, [
                                                0., systPkl["zz+zz4l"]["lepSF"]["binMC_1"], systPkl["zz+zz4l"]["lepSF"]["binMC_2"]])

    #zgphotonSFSyst=hlp.getScaleUncertHisto(zgHist,[0.,      systPkl["zg"]["photonSF"]["binMC_1"],systPkl["zg"]["photonSF"]["binMC_2"]],zgHist.Integral()/sum(zgamma.ngens))
    #ttgphotonSFSyst=hlp.getScaleUncertHisto(ttgHist,[0.,    systPkl["ttg"]["photonSF"]["binMC_1"],systPkl["ttg"]["photonSF"]["binMC_2"]])
    #zzphotonSFSyst=hlp.getScaleUncertHisto(zzHist,[0.,      systPkl["zz"]["photonSF"]["binMC_1"],systPkl["zz"]["photonSF"]["binMC_2"]])
    #wwgphotonSFSyst=hlp.getScaleUncertHisto(wwgHist,[0.,    systPkl["wwg"]["photonSF"]["binMC_1"],systPkl["wwg"]["photonSF"]["binMC_2"],0.1],wwgHist.Integral()/sum(wwgamma.ngens))
    #wzgphotonSFSyst=hlp.getScaleUncertHisto(wzgHist,[0.,    systPkl["wzg"]["photonSF"]["binMC_1"],systPkl["wzg"]["photonSF"]["binMC_2"],0.1],wzgHist.Integral()/sum(wzgamma.ngens))
    #dyphotonSFSyst=hlp.getScaleUncertHisto(dyHist,[0.,      systPkl["dy"]["photonSF"]["binMC_1"],systPkl["dy"]["photonSF"]["binMC_2"]],dyHist.Integral()/DYjetsNLO.ngens[0])
    # wjetsphotonSFSyst=hlp.getScaleUncertHisto(wjetsHist,[0.,systPkl["wjets"]["photonSF"]["binMC_1"],systPkl["wjets"]["photonSF"]["binMC_2"],0.1],wjetsHist.Integral()/sum(wjets.ngens))
    #ttphotonSFSyst=hlp.getScaleUncertHisto(ttHist,[0.,      systPkl["tt"]["photonSF"]["binMC_1"],systPkl["tt"]["photonSF"]["binMC_2"]])
    # singletopphotonSFSyst=hlp.getScaleUncertHisto(singletopHist,[0.,systPkl["singletop"]["photonSF"]["binMC_1"],systPkl["singletop"]["photonSF"]["binMC_2"],0.1],singletopHist.Integral()/sum(singletop.ngens))
    wzphotonSFSyst = hlp.getScaleUncertHisto(
        wzHist, [0.,      systPkl["wz"]["photonSF"]["binMC_1"], systPkl["wz"]["photonSF"]["binMC_2"]])
    #wwphotonSFSyst=hlp.getScaleUncertHisto(wwHist,[0.,      systPkl["ww"]["photonSF"]["binMC_1"],systPkl["ww"]["photonSF"]["binMC_2"],0.1],wwHist.Integral()/sum(ww.ngens))
    #zz4lphotonSFSyst=hlp.getScaleUncertHisto(zz4lHist,[0.,  systPkl["zz4l"]["photonSF"]["binMC_1"],systPkl["zz4l"]["photonSF"]["binMC_2"]])
    #wgphotonSFSyst=hlp.getScaleUncertHisto(wgHist,[0.,      systPkl["wgamma"]["photonSF"]["binMC_1"],systPkl["wgamma"]["photonSF"]["binMC_2"],0.1],wgHist.Integral()/sum(wgamma.ngens))
    otherphotonSFSyst = hlp.getScaleUncertHisto(final_otherHist, [
                                                0.,      systPkl["other"]["photonSF"]["binMC_1"], systPkl["other"]["photonSF"]["binMC_2"]], final_otherHist.Integral() / sum(otherBKG.ngens))
    final_ttphotonSFSyst = hlp.getScaleUncertHisto(final_ttHist, [
                                                   0., systPkl["tt+ttg"]["photonSF"]["binMC_1"], systPkl["tt+ttg"]["photonSF"]["binMC_2"]])
    final_tt080photonSFSyst = hlp.getScaleUncertHisto(final_tt080Hist, [
                                                      0., systPkl["tt+ttg"]["photonSF"]["binMC_1"], systPkl["tt+ttg"]["photonSF"]["binMC_2"]])
    final_tt80photonSFSyst = hlp.getScaleUncertHisto(final_tt80Hist, [
                                                     0., systPkl["tt+ttg"]["photonSF"]["binMC_1"], systPkl["tt+ttg"]["photonSF"]["binMC_2"]])
    final_dyphotonSFSyst = hlp.getScaleUncertHisto(final_dyHist, [
                                                   0., systPkl["dy+zg"]["photonSF"]["binMC_1"], systPkl["dy+zg"]["photonSF"]["binMC_2"]], final_dyHist.Integral() / (sum(DYjetsNLO.ngens) + sum(zgamma.ngens)))
    final_zzphotonSFSyst = hlp.getScaleUncertHisto(final_zzHist, [
                                                   0., systPkl["zz+zz4l"]["photonSF"]["binMC_1"], systPkl["zz+zz4l"]["photonSF"]["binMC_2"]])

    # zgEWKSyst=hlp.getScaleUncertHisto(zgHist,[0.,systPkl["zg"]["EWK"]["binMC_1"],systPkl["zg"]["EWK"]["binMC_2"]],zgHist.Integral()/sum(zgamma.ngens))
    # ttgEWKSyst=hlp.getScaleUncertHisto(ttgHist,[0.,systPkl["ttg"]["EWK"]["binMC_1"],systPkl["ttg"]["EWK"]["binMC_2"]])
    # zzEWKSyst=hlp.getScaleUncertHisto(zzHist,[0.,systPkl["zz"]["EWK"]["binMC_1"],systPkl["zz"]["EWK"]["binMC_2"]])
    # wwgEWKSyst=hlp.getScaleUncertHisto(wwgHist,[0.,systPkl["wwg"]["EWK"]["binMC_1"],systPkl["wwg"]["EWK"]["binMC_2"],0.1],wwgHist.Integral()/sum(wwgamma.ngens))
    # wzgEWKSyst=hlp.getScaleUncertHisto(wzgHist,[0.,systPkl["wzg"]["EWK"]["binMC_1"],systPkl["wzg"]["EWK"]["binMC_2"],0.1],wzgHist.Integral()/sum(wzgamma.ngens))
    # dyEWKSyst=hlp.getScaleUncertHisto(dyHist,[0.,systPkl["dy"]["EWK"]["binMC_1"],systPkl["dy"]["EWK"]["binMC_2"]],dyHist.Integral()/DYjetsNLO.ngens[0])
    # wjetsEWKSyst=hlp.getScaleUncertHisto(wjetsHist,[0.,systPkl["wjets"]["EWK"]["binMC_1"],systPkl["wjets"]["EWK"]["binMC_2"],0.1],wjetsHist.Integral()/sum(wjets.ngens))
    # ttEWKSyst=hlp.getScaleUncertHisto(ttHist,[0.,systPkl["tt"]["EWK"]["binMC_1"],systPkl["tt"]["EWK"]["binMC_2"]])
    # singletopEWKSyst=hlp.getScaleUncertHisto(singletopHist,[0.,systPkl["singletop"]["EWK"]["binMC_1"],systPkl["singletop"]["EWK"]["binMC_2"],0.1],singletopHist.Integral()/sum(singletop.ngens))
    wzEWKSyst = hlp.getScaleUncertHisto(
        wzHist, [0., systPkl["wz"]["EWK"]["binMC_1"], systPkl["wz"]["EWK"]["binMC_2"]])
    # wwEWKSyst=hlp.getScaleUncertHisto(wwHist,[0.,systPkl["ww"]["EWK"]["binMC_1"],systPkl["ww"]["EWK"]["binMC_2"],0.1],wwHist.Integral()/sum(ww.ngens))
    # zz4lEWKSyst=hlp.getScaleUncertHisto(zz4lHist,[0.,systPkl["zz4l"]["EWK"]["binMC_1"],systPkl["zz4l"]["EWK"]["binMC_2"]])
    # wgEWKSyst=hlp.getScaleUncertHisto(wgHist,[0.,systPkl["wgamma"]["EWK"]["binMC_1"],systPkl["wgamma"]["EWK"]["binMC_2"],0.1],wgHist.Integral()/sum(wgamma.ngens))
    otherEWKSyst = hlp.getScaleUncertHisto(final_otherHist, [
                                           0., systPkl["other"]["EWK"]["binMC_1"], systPkl["other"]["EWK"]["binMC_2"]], final_otherHist.Integral() / sum(otherBKG.ngens))
    final_ttEWKSyst = hlp.getScaleUncertHisto(final_ttHist, [
                                              0., systPkl["tt+ttg"]["EWK"]["binMC_1"], systPkl["tt+ttg"]["EWK"]["binMC_2"]])
    final_tt080EWKSyst = hlp.getScaleUncertHisto(final_tt080Hist, [
                                                 0., systPkl["tt+ttg"]["EWK"]["binMC_1"], systPkl["tt+ttg"]["EWK"]["binMC_2"]])
    final_tt80EWKSyst = hlp.getScaleUncertHisto(final_tt80Hist, [
                                                0., systPkl["tt+ttg"]["EWK"]["binMC_1"], systPkl["tt+ttg"]["EWK"]["binMC_2"]])
    final_dyEWKSyst = hlp.getScaleUncertHisto(final_dyHist, [
                                              0., systPkl["dy+zg"]["EWK"]["binMC_1"], systPkl["dy+zg"]["EWK"]["binMC_2"]], final_dyHist.Integral() / (sum(DYjetsNLO.ngens) + sum(zgamma.ngens)))
    final_zzEWKSyst = hlp.getScaleUncertHisto(final_zzHist, [
                                              0., systPkl["zz+zz4l"]["EWK"]["binMC_1"], systPkl["zz+zz4l"]["EWK"]["binMC_2"]])

    # zgISRSyst=hlp.getScaleUncertHisto(zgHist,[0.,systPkl["zg"]["ISR"]["binMC_1"],systPkl["zg"]["ISR"]["binMC_2"]],zgHist.Integral()/sum(zgamma.ngens))
    # ttgISRSyst=hlp.getScaleUncertHisto(ttgHist,[0.,systPkl["ttg"]["ISR"]["binMC_1"],systPkl["ttg"]["ISR"]["binMC_2"]])
    # zzISRSyst=hlp.getScaleUncertHisto(zzHist,[0.,systPkl["zz"]["ISR"]["binMC_1"],systPkl["zz"]["ISR"]["binMC_2"]])
    # wwgISRSyst=hlp.getScaleUncertHisto(wwgHist,[0.,systPkl["wwg"]["ISR"]["binMC_1"],systPkl["wwg"]["ISR"]["binMC_2"],0.1],wwgHist.Integral()/sum(wwgamma.ngens))
    # wzgISRSyst=hlp.getScaleUncertHisto(wzgHist,[0.,systPkl["wzg"]["ISR"]["binMC_1"],systPkl["wzg"]["ISR"]["binMC_2"],0.1],wzgHist.Integral()/sum(wzgamma.ngens))
    # dyISRSyst=hlp.getScaleUncertHisto(dyHist,[0.,systPkl["dy"]["ISR"]["binMC_1"],systPkl["dy"]["ISR"]["binMC_2"]],dyHist.Integral()/DYjetsNLO.ngens[0])
    # wjetsISRSyst=hlp.getScaleUncertHisto(wjetsHist,[0.,systPkl["wjets"]["ISR"]["binMC_1"],systPkl["wjets"]["ISR"]["binMC_2"],0.1],wjetsHist.Integral()/sum(wjets.ngens))
    # ttISRSyst=hlp.getScaleUncertHisto(ttHist,[0.,systPkl["tt"]["ISR"]["binMC_1"],systPkl["tt"]["ISR"]["binMC_2"]])
    # singletopISRSyst=hlp.getScaleUncertHisto(singletopHist,[0.,systPkl["singletop"]["ISR"]["binMC_1"],systPkl["singletop"]["ISR"]["binMC_2"],0.1],singletopHist.Integral()/sum(singletop.ngens))
    wzISRSyst = hlp.getScaleUncertHisto(
        wzHist, [0., systPkl["wz"]["ISR"]["binMC_1"], systPkl["wz"]["ISR"]["binMC_2"]])
    # wwISRSyst=hlp.getScaleUncertHisto(wwHist,[0.,systPkl["ww"]["ISR"]["binMC_1"],systPkl["ww"]["ISR"]["binMC_2"],0.1],wwHist.Integral()/sum(ww.ngens))
    # zz4lISRSyst=hlp.getScaleUncertHisto(zz4lHist,[0.,systPkl["zz4l"]["ISR"]["binMC_1"],systPkl["zz4l"]["ISR"]["binMC_2"]])
    # wgISRSyst=hlp.getScaleUncertHisto(wgHist,[0.,systPkl["wgamma"]["ISR"]["binMC_1"],systPkl["wgamma"]["ISR"]["binMC_2"],0.1],wgHist.Integral()/sum(wgamma.ngens))
    otherISRSyst = hlp.getScaleUncertHisto(final_otherHist, [
                                           0., systPkl["other"]["ISR"]["binMC_1"], systPkl["other"]["ISR"]["binMC_2"]], final_otherHist.Integral() / sum(otherBKG.ngens))
    final_ttISRSyst = hlp.getScaleUncertHisto(final_ttHist, [
                                              0., systPkl["tt+ttg"]["ISR"]["binMC_1"], systPkl["tt+ttg"]["ISR"]["binMC_2"]])
    final_tt080ISRSyst = hlp.getScaleUncertHisto(final_tt080Hist, [
                                                 0., systPkl["tt+ttg"]["ISR"]["binMC_1"], systPkl["tt+ttg"]["ISR"]["binMC_2"]])
    final_tt80ISRSyst = hlp.getScaleUncertHisto(final_tt80Hist, [
                                                0., systPkl["tt+ttg"]["ISR"]["binMC_1"], systPkl["tt+ttg"]["ISR"]["binMC_2"]])
    final_dyISRSyst = hlp.getScaleUncertHisto(final_dyHist, [
                                              0., systPkl["dy+zg"]["ISR"]["binMC_1"], systPkl["dy+zg"]["ISR"]["binMC_2"]], final_dyHist.Integral() / (sum(DYjetsNLO.ngens) + sum(zgamma.ngens)))
    final_zzISRSyst = hlp.getScaleUncertHisto(final_zzHist, [
                                              0., systPkl["zz+zz4l"]["ISR"]["binMC_1"], systPkl["zz+zz4l"]["ISR"]["binMC_2"]])

    #zgSyst = aux.getSysHistoWithMeanWeight(zgHist, mcSystUncert,zgHist.Integral()/sum(zgamma.ngens))
    #ttgSyst = aux.getSysHisto(ttgHist, mcSystUncert)
    #zzSyst = aux.getSysHisto(zzHist, mcSystUncert)
    #wwgSyst = aux.getSysHistoWithMeanWeight(wwgHist, mcSystUncert,wwgHist.Integral()/sum(wwgamma.ngens))
    #wzgSyst = aux.getSysHistoWithMeanWeight(wzgHist, mcSystUncert,wzgHist.Integral()/sum(wzgamma.ngens))
    #dySyst = aux.getSysHisto(dyHist, mcSystUncert)
    #dySyst = aux.getSysHistoWithMeanWeight(dyHist, mcSystUncert,dyHist.Integral()/DYjetsNLO.ngens[0])
    #wjetsSyst = aux.getSysHistoWithMeanWeight(wjetsHist, mcSystUncert,wjetsHist.Integral()/sum(wjets.ngens))
    #ttSyst = aux.getSysHisto(ttHist, mcSystUncert)
    # singletopSyst=aux.getSysHistoWithMeanWeight(singletopHist,mcSystUncert,singletopHist.Integral()/sum(singletop.ngens))
    wzSyst = aux.getSysHisto(wzHist, mcSystUncert)
    # wwSyst=aux.getSysHistoWithMeanWeight(wwHist,mcSystUncert,wwHist.Integral()/sum(ww.ngens))
    # zz4lSyst=aux.getSysHisto(zz4lHist,mcSystUncert)
    # wgSyst=aux.getSysHistoWithMeanWeight(wgHist,mcSystUncert,wgHist.Integral()/sum(wgamma.ngens))
    otherSyst = aux.getSysHistoWithMeanWeight(
        final_otherHist, mcSystUncert, final_otherHist.Integral() / sum(otherBKG.ngens))
    final_ttSyst = aux.getSysHisto(final_ttHist, mcSystUncert)
    final_tt080Syst = aux.getSysHisto(final_tt080Hist, mcSystUncert)
    final_tt80Syst = aux.getSysHisto(final_tt80Hist, mcSystUncert)
    final_dySyst = aux.getSysHistoWithMeanWeight(
        final_dyHist, mcSystUncert, final_dyHist.Integral() / (sum(DYjetsNLO.ngens) + sum(zgamma.ngens)))
    final_zzSyst = aux.getSysHisto(final_zzHist, mcSystUncert)

    #zgTriggerSyst = aux.getSysHistoWithMeanWeight(zgHist, 0.03,zgHist.Integral()/sum(zgamma.ngens))
    #ttgTriggerSyst = aux.getSysHisto(ttgHist, 0.03)
    #zzTriggerSyst = aux.getSysHisto(zzHist, 0.03)
    #wwgTriggerSyst = aux.getSysHistoWithMeanWeight(wwgHist, 0.03,wwgHist.Integral()/sum(wwgamma.ngens))
    #wzgTriggerSyst = aux.getSysHistoWithMeanWeight(wzgHist, 0.03,wzgHist.Integral()/sum(wzgamma.ngens))
    #dyTriggerSyst = aux.getSysHistoWithMeanWeight(dyHist, 0.03,dyHist.Integral()/DYjetsNLO.ngens[0])
    #wjetsTriggerSyst = aux.getSysHistoWithMeanWeight(wjetsHist, 0.03,wjetsHist.Integral()/sum(wjets.ngens))
    #ttTriggerSyst = aux.getSysHisto(ttHist, 0.03)
    # singletopTriggerSyst=aux.getSysHistoWithMeanWeight(singletopHist,0.03,singletopHist.Integral()/sum(singletop.ngens))
    wzTriggerSyst = aux.getSysHisto(wzHist, 0.03)
    # wwTriggerSyst=aux.getSysHistoWithMeanWeight(wwHist,0.03,wwHist.Integral()/sum(ww.ngens))
    # zz4lTriggerSyst=aux.getSysHisto(zz4lHist,0.03)
    # wgTriggerSyst=aux.getSysHistoWithMeanWeight(wgHist,0.03,final_otherHist.Integral()/sum(otherBKG.ngens))
    otherTriggerSyst = aux.getSysHistoWithMeanWeight(
        final_otherHist, 0.03, final_otherHist.Integral() / sum(otherBKG.ngens))
    final_ttTriggerSyst = aux.getSysHisto(final_ttHist, 0.03)
    final_tt080TriggerSyst = aux.getSysHisto(final_tt080Hist, 0.03)
    final_tt80TriggerSyst = aux.getSysHisto(final_tt80Hist, 0.03)
    final_dyTriggerSyst = aux.getSysHistoWithMeanWeight(
        final_dyHist, 0.03, final_dyHist.Integral() / (sum(DYjetsNLO.ngens) + sum(zgamma.ngens)))
    final_zzTriggerSyst = aux.getSysHisto(final_zzHist, 0.03)

    #zgXSECSyst = aux.getSysHistoWithMeanWeight(zgHist, 0.,zgHist.Integral()/sum(zgamma.ngens))
    #ttgXSECSyst = aux.getSysHisto(ttgHist, 0.)
    #zzXSECSyst = aux.getSysHisto(zzHist, 0.)
    #wwgXSECSyst = aux.getSysHistoWithMeanWeight(wwgHist, 0.5,wwgHist.Integral()/sum(wwgamma.ngens))
    #wzgXSECSyst = aux.getSysHistoWithMeanWeight(wzgHist, 0.5,wzgHist.Integral()/sum(wzgamma.ngens))
    #dyXSECSyst = aux.getSysHistoWithMeanWeight(dyHist, 0.,dyHist.Integral()/DYjetsNLO.ngens[0])
    #wjetsXSECSyst = aux.getSysHistoWithMeanWeight(wjetsHist, 0.5,wjetsHist.Integral()/sum(wjets.ngens))
    #ttXSECSyst = aux.getSysHisto(ttHist, 0.)
    # singletopXSECSyst=aux.getSysHistoWithMeanWeight(singletopHist,0.5,singletopHist.Integral()/sum(singletop.ngens))
    # wzXSECSyst=aux.getSysHisto(wzHist,0.)
    # wwXSECSyst=aux.getSysHistoWithMeanWeight(wwHist,0.5,wwHist.Integral()/sum(ww.ngens))
    # zz4lXSECSyst=aux.getSysHisto(zz4lHist,0.)
    # wgXSECSyst=aux.getSysHistoWithMeanWeight(wgHist,0.5,final_otherHist.Integral()/sum(otherBKG.ngens))
    otherXSECSyst = aux.getSysHistoWithMeanWeight(
        final_otherHist, 0.5, final_otherHist.Integral() / sum(otherBKG.ngens))

    #zgSFSyst = aux.getSysHistoWithMeanWeight(zgHist, sfDYErr,zgHist.Integral()/sum(zgamma.ngens))
    #ttgSFSyst = aux.getSysHisto(ttgHist, 0.04)
    #ttg80SFSyst = aux.getSysHisto(ttg80Hist, 0.4)
    #ttg80SFSyst = aux.getSysHisto(ttg80Hist, 0.04)
    #ttg080SFSyst = aux.getSysHisto(ttg080Hist, 0.4)
    #ttg080SFSyst = aux.getSysHisto(ttg080Hist, 0.04)
    #zzSFSyst = aux.getSysHisto(zzHist, sfZZErr)
    #wwgSFSyst = aux.getSysHistoWithMeanWeight(wwgHist, 0.,wwgHist.Integral()/sum(wwgamma.ngens))
    #wzgSFSyst = aux.getSysHistoWithMeanWeight(wzgHist, 0.,wzgHist.Integral()/sum(wzgamma.ngens))
    #dySFSyst = aux.getSysHistoWithMeanWeight(dyHist, sfDYErr,dyHist.Integral()/DYjetsNLO.ngens[0])
    #wjetsSFSyst = aux.getSysHistoWithMeanWeight(wjetsHist, 0.,wjetsHist.Integral()/sum(wjets.ngens))
    #ttSFSyst = aux.getSysHisto(ttHist,0.04)
    #tt80SFSyst = aux.getSysHisto(tt80Hist,0.4)
    #tt080SFSyst = aux.getSysHisto(tt080Hist,0.4)
    #tt80SFSyst = aux.getSysHisto(tt80Hist,0.04)
    #tt080SFSyst = aux.getSysHisto(tt080Hist,0.04)
    # singletopSFSyst=aux.getSysHistoWithMeanWeight(singletopHist,0.,singletopHist.Integral()/sum(singletop.ngens))
    # wzSFSyst=aux.getSysHisto(wzHist,sfWZErr)
    wzSFSyst = aux.getSysHisto(wzHist, sfWZErr / sfWZ)
    # wwSFSyst=aux.getSysHistoWithMeanWeight(wwHist,0.,wwHist.Integral()/sum(ww.ngens))
    # zz4lSFSyst=aux.getSysHisto(zz4lHist,sfZZErr)
    # wgSFSyst=aux.getSysHistoWithMeanWeight(wgHist,0.,wgHist.Integral()/sum(wgamma.ngens))
    # otherSFSyst=aux.getSysHistoWithMeanWeight(final_otherHist,0.,final_otherHist.Integral()/sum(otherBKG.ngens))
    #final_tt080SFSyst= aux.getSysHisto(final_tt080Hist,0.04)
    final_tt080SFSyst = aux.getSysHisto(final_tt080Hist, sfTTErr)
    #final_tt80SFSyst= aux.getSysHisto(final_tt80Hist,0.4)
    final_tt80SFSyst = aux.getSysHisto(final_tt80Hist, 0.2)
    #final_ttSFSyst= aux.getSysHisto(final_tt80Hist,0.4)
    #final_zzSFSyst= aux.getSysHisto(final_zzHist,sfZZErr)
    final_zzSFSyst = aux.getSysHisto(final_zzHist, sfZZErr / sfZZ)
    #final_dySFSyst= aux.getSysHistoWithMeanWeight(final_dyHist,sfDYErr,final_dyHist.Integral()/(sum(DYjetsNLO.ngens)+sum(zgamma.ngens)))
    final_dySFSyst = aux.getSysHistoWithMeanWeight(
        final_dyHist, sfDYErr / sfDY, final_dyHist.Integral() / (sum(DYjetsNLO.ngens) + sum(zgamma.ngens)))

    # ttgTotSFSyst=aux.addHists(ttg080SFSyst,ttg80SFSyst)
    # ttTotSFSyst=aux.addHists(tt080SFSyst,tt80SFSyst)
    # final_ttTotSFSyst=aux.addHists(final_tt080SFSyst,final_tt80SFSyst)
    #print sfTT,sfTTErr,sfTTErr/sfTT
    final_ttTotSFSyst = aux.getSysHisto(final_ttHist, sfTTErr / sfTT)
    # ttgTotSFSyst=ttgSFSyst
    # ttTotSFSyst=ttSFSyst

    # zgTotSyst=aux.addUncertaintiesQuadratic([zgSyst,zgTriggerSyst,zgphotonSFSyst,zglepSFSyst,zgJESSyst,zgJERSyst,zgPUSyst,zgScaleSyst,zgSFSyst,zgPDFSyst,zgEWKSyst,zgISRSyst,zgXSECSyst])
    # ttgTotSyst=aux.addUncertaintiesQuadratic([ttgSyst,ttgTriggerSyst,ttgphotonSFSyst,ttglepSFSyst,ttgJESSyst,ttgJERSyst,ttgPUSyst,ttgScaleSyst,ttgSFSyst])
    # ttgTotSyst=aux.addUncertaintiesQuadratic([ttgSyst,ttgTriggerSyst,ttgphotonSFSyst,ttglepSFSyst,ttgJESSyst,ttgJERSyst,ttgPUSyst,ttgScaleSyst,ttgTotSFSyst,ttgPDFSyst,ttgEWKSyst,ttgISRSyst,ttgXSECSyst])
    # zzTotSyst=aux.addUncertaintiesQuadratic([zzSyst,zzTriggerSyst,zzphotonSFSyst,zzlepSFSyst,zzJESSyst,zzJERSyst,zzPUSyst,zzScaleSyst,zzSFSyst,zzPDFSyst,zzEWKSyst,zzISRSyst,zzXSECSyst])
    # wwgTotSyst=aux.addUncertaintiesQuadratic([wwgSyst,wwgTriggerSyst,wwgphotonSFSyst,wwglepSFSyst,wwgJESSyst,wwgJERSyst,wwgPUSyst,wwgScaleSyst,wwgSFSyst,wwgPDFSyst,wwgEWKSyst,wwgISRSyst,wwgXSECSyst])
    # wzgTotSyst=aux.addUncertaintiesQuadratic([wzgSyst,wzgTriggerSyst,wzgphotonSFSyst,wzglepSFSyst,wzgJESSyst,wzgJERSyst,wzgPUSyst,wzgScaleSyst,wzgSFSyst,wzgPDFSyst,wzgEWKSyst,wzgISRSyst,wzgXSECSyst])
    #print "---------"
    # dyTotSyst=aux.addUncertaintiesQuadratic([dySyst,dyTriggerSyst,dyphotonSFSyst,dylepSFSyst,dyJESSyst,dyJERSyst,dyPUSyst,dyScaleSyst,dySFSyst,dyPDFSyst,dyEWKSyst,dyISRSyst,dyXSECSyst])
    #print "---------"
    # wjetsTotSyst=aux.addUncertaintiesQuadratic([wjetsSyst,wjetsTriggerSyst,wjetsphotonSFSyst,wjetslepSFSyst,wjetsJESSyst,wjetsJERSyst,wjetsPUSyst,wjetsScaleSyst,wjetsSFSyst,wjetsPDFSyst,wjetsEWKSyst,wjetsISRSyst,wjetsXSECSyst])
    # ttTotSyst=aux.addUncertaintiesQuadratic([ttSyst,ttTriggerSyst,ttphotonSFSyst,ttlepSFSyst,ttJESSyst,ttJERSyst,ttPUSyst,ttScaleSyst,ttSFSyst])
    # ttTotSyst=aux.addUncertaintiesQuadratic([ttSyst,ttTriggerSyst,ttphotonSFSyst,ttlepSFSyst,ttJESSyst,ttJERSyst,ttPUSyst,ttScaleSyst,ttTotSFSyst,ttPDFSyst,ttEWKSyst,ttISRSyst,ttXSECSyst])
    # singletopTotSyst=aux.addUncertaintiesQuadratic([singletopSyst,singletopTriggerSyst,singletopphotonSFSyst,singletoplepSFSyst,singletopJESSyst,singletopJERSyst,singletopPUSyst,singletopScaleSyst,singletopSFSyst,singletopPDFSyst,singletopEWKSyst,singletopISRSyst,singletopXSECSyst])
    # wzTotSyst=aux.addUncertaintiesQuadratic([wzSyst,wzTriggerSyst,wzphotonSFSyst,wzlepSFSyst,wzJESSyst,wzJERSyst,wzPUSyst,wzScaleSyst,wzSFSyst,wzPDFSyst,wzEWKSyst,wzISRSyst,wzXSECSyst])
    # wwTotSyst=aux.addUncertaintiesQuadratic([wwSyst,wwTriggerSyst,wwphotonSFSyst,wwlepSFSyst,wwJESSyst,wwJERSyst,wwPUSyst,wwScaleSyst,wwSFSyst,wwPDFSyst,wwEWKSyst,wwISRSyst,wwXSECSyst])
    # zz4lTotSyst=aux.addUncertaintiesQuadratic([zz4lSyst,zz4lTriggerSyst,zz4lphotonSFSyst,zz4llepSFSyst,zz4lJESSyst,zz4lJERSyst,zz4lPUSyst,zz4lScaleSyst,zz4lSFSyst,zz4lPDFSyst,zz4lEWKSyst,zz4lISRSyst,zz4lXSECSyst])
    # wgTotSyst=aux.addUncertaintiesQuadratic([wgSyst,wgTriggerSyst,wgphotonSFSyst,wglepSFSyst,wgJESSyst,wgJERSyst,wgPUSyst,wgScaleSyst,wgSFSyst,wgPDFSyst,wgEWKSyst,wgISRSyst,wgXSECSyst])
    # wgTotSyst=aux.addUncertaintiesQuadratic([otherSyst,otherTriggerSyst,otherphotonSFSyst,otherlepSFSyst,otherJESSyst,otherJERSyst,otherPUSyst,otherScaleSyst,otherSFSyst,otherPDFSyst,otherEWKSyst,otherISRSyst,otherXSECSyst])
    wzTotSyst = aux.addUncertaintiesQuadratic([wzSyst, wzTriggerSyst, wzphotonSFSyst, wzlepSFSyst,
                                               wzJESSyst, wzJERSyst, wzPUSyst, wzScaleSyst, wzSFSyst, wzPDFSyst, wzEWKSyst, wzISRSyst])
    zzTotSyst = aux.addUncertaintiesQuadratic([final_zzSyst, final_zzTriggerSyst, final_zzphotonSFSyst, final_zzlepSFSyst, final_zzJESSyst,
                                               final_zzJERSyst, final_zzPUSyst, final_zzScaleSyst, final_zzSFSyst, final_zzPDFSyst, final_zzEWKSyst, final_zzISRSyst])
    dyTotSyst = aux.addUncertaintiesQuadratic([final_dySyst, final_dyTriggerSyst, final_dyphotonSFSyst, final_dylepSFSyst, final_dyJESSyst,
                                               final_dyJERSyst, final_dyPUSyst, final_dyScaleSyst, final_dySFSyst, final_dyPDFSyst, final_dyEWKSyst, final_dyISRSyst])
    ttTotSyst = aux.addUncertaintiesQuadratic([final_ttSyst, final_ttTriggerSyst, final_ttphotonSFSyst, final_ttlepSFSyst, final_ttJESSyst,
                                               final_ttJERSyst, final_ttPUSyst, final_ttScaleSyst, final_ttTotSFSyst, final_ttPDFSyst, final_ttEWKSyst, final_ttISRSyst])
    otherTotSyst = aux.addUncertaintiesQuadratic([otherSyst, otherTriggerSyst, otherphotonSFSyst, otherlepSFSyst,
                                                  otherJESSyst, otherJERSyst, otherPUSyst, otherScaleSyst, otherPDFSyst, otherEWKSyst, otherISRSyst, otherXSECSyst])

    # for bin in range(dyTotSyst.GetNbinsX()-2, dyTotSyst.GetNbinsX()+1):
    #c = ttgTotSyst.GetBinContent(bin)
    # return 1.+max(0,h.GetBinError(b)/c if c else 0)
    #print "-1",dySyst.GetBinError(bin)
    #print "-2",dyTriggerSyst.GetBinError(bin)
    #print "-3",dyphotonSFSyst.GetBinError(bin)
    #print "-4",dylepSFSyst.GetBinError(bin)
    #print "-5",dyJESSyst.GetBinError(bin)
    #print "-6",dyJERSyst.GetBinError(bin)
    #print "-7",dyPUSyst.GetBinError(bin)
    #print "-8",dySFSyst.GetBinError(bin)
    #print "-9",dyPDFSyst.GetBinError(bin)
    #print "-10",dyEWKSyst.GetBinError(bin)
    #print "-11",dyISRSyst.GetBinError(bin)
    #print "-12",dyXSECSyst.GetBinError(bin)
    #print "-13",dyTotSyst.GetBinError(bin)
    #print ttgTotSyst.GetBinError(bin)/c,c

    #dirSyst= aux.getSysHisto(dirHist,mcSystUncert)
    #dirSyst= aux.getSysHisto(aux.addHists(zgSyst,ttgSyst,zzSyst,wwgSyst,wzgSyst,dySyst,wjetsSyst,ttSyst,singletopSyst,wzSyst,wwSyst,zz4lSyst,wgSyst),mcSystUncert)
    #dirSyst= aux.getSysHisto(aux.addHists(zgSyst,ttgSyst,zzSyst,dySyst,ttSyst,wzSyst,zz4lSyst,otherSyst),mcSystUncert)
    # dirSyst= aux.getSysHisto(aux.addHists(,,zzSyst,dySyst,ttSyst,wzSyst,,otherSyst),mcSystUncert)

    #zgSyst = aux.addUncertaintiesQuadratic([zgSyst,zgPdfUnc,zgScaleUnc,zgJesUnc,zgPuUnc])
    #wgSyst = aux.addUncertaintiesQuadratic([wgSyst,wgPdfUnc,wgScaleUnc,wgJesUnc,wgPuUnc])
    #tgSyst = aux.addUncertaintiesQuadratic([tgSyst,tgPdfUnc,tgScaleUnc,tgJesUnc,tgPuUnc])

    #totStat = aux.addHists(gjetHist, eHist, zgHist, wgHist, tgHist)
    #totSyst = aux.addHists(gjetSyst, eSyst, zgSyst, wgSyst, tgSyst)
    #totStat = aux.addHists(zgHist, ttgHist, zzHist, wwgHist, wzgHist,dyHist,wjetsHist,ttHist)
    #totSyst = aux.addHists(zgSyst, ttgSyst, zzSyst, wwgSyst, wzgSyst,dySyst,wjetsSyst,ttSyst)
    #totStat = aux.addHists(zgHist, ttgHist, zzHist, wwgHist, wzgHist,dyHist,wjetsHist,ttHist,singletopHist,wzHist,wwHist,zz4lHist,wgHist)
    #totSyst = aux.addHists(zgSyst, ttgSyst, zzSyst, wwgSyst, wzgSyst,dySyst,wjetsSyst,ttSyst,singletopSyst,wzSyst,wwSyst,zz4lSyst,wjetsSyst,wgSyst)
    #totSyst = aux.addHists(zgTotSyst, ttgTotSyst, zzTotSyst, wwgTotSyst, wzgTotSyst,dyTotSyst,wjetsTotSyst,ttTotSyst,singletopTotSyst,wzTotSyst,wwTotSyst,zz4lTotSyst,wjetsTotSyst,wgTotSyst)
    #totSyst = aux.addHists(zgTotSyst, ttgTotSyst, zzTotSyst, wwgTotSyst, wzgTotSyst,dyTotSyst,wjetsTotSyst,ttTotSyst,singletopTotSyst,wzTotSyst,wwTotSyst,zz4lTotSyst,wgTotSyst)
    totStat = aux.addHists(final_zzHist, final_dyHist,
                           final_ttHist, wzHist, final_otherHist)
    totSyst = aux.addHists(dyTotSyst, zzTotSyst,
                           ttTotSyst, wzTotSyst, otherTotSyst)
    #totSyst = dyXSECSyst

    #signal2 = aux.stdHistWithoutNGenWithWeights(t5bbbbzg_1500_400,"Zg_0_0/sig/LL/nom"+"/met",["nISR","topPt","ewk"],nBins)
    #signal2 = aux.stdHistWithoutNGenWithWeights(t5bbbbzg_1500_1400,"Zg_0_0/sig/LL/nom"+"/met",["topPt"],nBins)
    #signal2 = aux.stdHist(tching_400, dirDir+"/met", nBins)
    #signal2 = aux.stdHist(t5bbbbzg_1500_400, dirDir+"/met", nBins)
    #signal2 = aux.stdHistWithoutNGen(gmsb_290_205, dirDir+"/met", nBins)
    #signal2 = aux.stdHistWithoutNGen(gmsb_290_205, "GMSB_0_0/sig/LL/nom"+"/met", nBins)
    signal2 = aux.stdHistWithoutNGenWithWeights(
        gmsb_290_205, "GMSB_0_0/sig/LL/nom" + "/met", ["nISR", "topPt", "ewk"], nBins)
    #signal1 = aux.stdHist(tching_400, dirDir+"/met", nBins)
    #signal1 = aux.stdHist(tching_600, dirDir+"/met", nBins)
    #signal1 = aux.stdHistWithoutNGen(tching_600, dirDir+"/met", nBins)
    #signal1 = aux.stdHistWithoutNGen(tching_600, "Ng_0_0/sig/LL/nom"+"/met", nBins)
    signal1 = aux.stdHistWithoutNGenWithWeights(
        tching_600, "Ng_0_0/sig/LL/nom" + "/met", ["nISR", "topPt", "ewk"], nBins)
    #signal1 = aux.stdHistWithWeights(tching_600,"onZMet150/LL"+"/met",["nISR","topPt","ewk"],nBins)
    # signal1.Scale(192644)
    #print signal1.Integral()

    signal1GenMetHisto = aux.stdHistWithoutNGenWithWeights(
        tching_600, "Ng_0_0/sig/LL/genmet" + "/met", ["nISR", "topPt", "ewk"], nBins)
    signal1lepSFuHisto = aux.stdHistWithoutNGenWithWeights(
        tching_600, "Ng_0_0/sig/LL/lepSFu" + "/met", ["nISR", "topPt", "ewk"], nBins)
    signal1lepSFdHisto = aux.stdHistWithoutNGenWithWeights(
        tching_600, "Ng_0_0/sig/LL/lepSFd" + "/met", ["nISR", "topPt", "ewk"], nBins)
    signal1photonSFuHisto = aux.stdHistWithoutNGenWithWeights(
        tching_600, "Ng_0_0/sig/LL/photonSFu" + "/met", ["nISR", "topPt", "ewk"], nBins)
    signal1photonSFdHisto = aux.stdHistWithoutNGenWithWeights(
        tching_600, "Ng_0_0/sig/LL/photonSFd" + "/met", ["nISR", "topPt", "ewk"], nBins)
    signal1JESuSFHisto = aux.stdHistWithoutNGenWithWeights(
        tching_600, "Ng_0_0/sig/LL/JESu" + "/met", ["nISR", "topPt", "ewk"], nBins)
    signal1JESdSFHisto = aux.stdHistWithoutNGenWithWeights(
        tching_600, "Ng_0_0/sig/LL/JESd" + "/met", ["nISR", "topPt", "ewk"], nBins)
    signal1JERuSFHisto = aux.stdHistWithoutNGenWithWeights(
        tching_600, "Ng_0_0/sig/LL/JERu" + "/met", ["nISR", "topPt", "ewk"], nBins)
    signal1JERdSFHisto = aux.stdHistWithoutNGenWithWeights(
        tching_600, "Ng_0_0/sig/LL/JERd" + "/met", ["nISR", "topPt", "ewk"], nBins)
    signal1NoPUuHisto = aux.stdHistWithoutNGenWithWeights(
        tching_600, "Ng_0_0/sig/LL/NoPUu" + "/met", ["nISR", "topPt", "ewk"], nBins)
    signal1NoPUdHisto = aux.stdHistWithoutNGenWithWeights(
        tching_600, "Ng_0_0/sig/LL/NoPUd" + "/met", ["nISR", "topPt", "ewk"], nBins)
    signal1NoPUHisto = aux.stdHistWithoutNGenWithWeights(
        tching_600, "Ng_0_0/sig/LL/NoPU" + "/met", ["nISR", "topPt", "ewk"], nBins)

    signal1EWKuHisto = aux.stdHistWithoutNGenWithWeights(
        tching_600, "Ng_0_0/sig/LL/EWKu" + "/met", ["nISR", "topPt", "ewkUp"], nBins)
    signal1EWKdHisto = aux.stdHistWithoutNGenWithWeights(
        tching_600, "Ng_0_0/sig/LL/EWKd" + "/met", ["nISR", "topPt", "ewkDn"], nBins)

    signal1ISRuHisto = aux.stdHistWithoutNGenWithWeights(
        tching_600, "Ng_0_0/sig/LL/ISRu" + "/met", ["nISRUp", "topPt", "ewk"], nBins)
    signal1ISRdHisto = aux.stdHistWithoutNGenWithWeights(
        tching_600, "Ng_0_0/sig/LL/ISRd" + "/met", ["nISRDn", "topPt", "ewk"], nBins)

    signal1NoPUuHisto.Scale(signal1.GetEntries() /
                            signal1NoPUuHisto.GetEntries())
    signal1NoPUdHisto.Scale(signal1.GetEntries() /
                            signal1NoPUdHisto.GetEntries())

    signal1scaleHists = dict([(i, aux.stdHistWithoutNGenWithWeights(tching_600, "Ng_0_0/sig/LL/" + str(
        i) + "/met", ["nISR", "topPt", "ewk", "pdf" + str(i)], nBins)) for i in range(0, 9)])
    #signal1PDFHists = dict([(i,aux.stdHistWithoutNGenWithWeights(tching_600,"Ng_0_0/sig/LL/"+str(i)+"/met",["nISR","topPt","ewk","pdf"+str(i)],nBins)) for i in range(9,110)])
    #print [h.Integral() for h in signal1scaleHists.values()]

    #signal1Dummy= aux.stdHistWithoutNGenWithWeights(tching_600,"Ng_0_0/sig/LL/nom"+"/met",["nISR","topPt","ewk"],nBins)
    #print tching_600.weights
    #print signal1Dummy.Integral()

    # for bin in range)

    #print signal1.GetEntries()

    for h in signal1, signal2:
        aux.drawOpt(h, "signal")
    #    h.Add(totStat)
    signal1.SetLineColor(ROOT.kBlue + 3)
    # signal2.SetLineColor(ROOT.kRed-3)
    signal2.SetLineColor(ROOT.kBlue + 3)
    signal2.SetLineStyle(2)

    #signal1_pre = aux.createHistoFromDatasetTree(t5wg_1600_100, "met*{}".format(info["shift"]), weight, nBins, "tr_jControl/simpleTree")
    # signal1_pre.Scale(info["scale"])

    #totSyst = gjetSyst

    #totUnc = aux.addHistUncert(totStat, totSyst)
    #totUnc = aux.addHistUncert(totStat, dirSyst)
    totUnc = aux.addHistUncert(totStat, totSyst)
    #totUnc = aux.addHistUncert(totStat)

    for bin in range(totUnc.GetNbinsX() - 1, totUnc.GetNbinsX() + 1):
        print "totUnc", totUnc.GetBinError(bin)

    #totUnc = aux.addHistUncert(dirSyst)
    aux.drawOpt(totUnc, "totUnc")
    #totUnc = dirSyst
    #aux.drawOpt(totUnc, "totUnc")

    c = ROOT.TCanvas()
    m = multiplot.Multiplot()
    # if dirSet == data:

    # if dirSet == dataDoubleSF:
    #m.add(dirHist, "Data")
    m.add(dataHist, "Data")
    # else:
    #m.add(dirHist, "Direct simulation")


#   m.add(signal1_pre, "contamination")
    #m.add(signal1, "T5bbbbZg")
    #m.add(signal2, "TChiNg")
    #m.add(signal2, "T5bbbbZg")
    m.add(signal2, "GMSB")
    # m.add(signal1, "TChiNg")
    m.add(signal1, "TChiZG")

    # m.addStack(zgHist, "Z#gamma")
    # m.addStack(ttgHist, "t#bar{t}#gamma")
    # m.addStack(zzHist, "ZZ(#rightarrowll#nu#nu)")
    # m.addStack(wwgHist, "WW#gamma")
    # m.addStack(wzgHist, "WZ#gamma")
    #m.addStack(dyHist, "Drell-Yan/Z")
    #m.addStack(wjetsHist, "W+jets")
    # m.addStack(ttHist, "t#bar{t}")
    #m.addStack(singletopHist, "single t")
    #m.addStack(wwHist, "WW")
    #m.addStack(wzHist, "WZ")
    # m.addStack(zz4lHist, "ZZ(#rightarrow4l)")
    # m.addStack(wgHist, "W#gamma")

    # different plot style:add bkgs
    # m.addStack(aux.addHists(zgHist,dyHist), "DY/Z#gamma")
    # m.addStack(aux.addHists(ttgHist,ttHist), "t#bar{t}(+#gamma)")
    #m.addStack(aux.addHists(zzHist,zz4lHist), "ZZ")
    #m.addStack(aux.addHists(wwgHist,wzgHist,wjetsHist,singletopHist,wwHist,wgHist), "other")
    m.addStack(final_dyHist, "DY/Z#gamma")
    m.addStack(final_ttHist, "t#bar{t}(+#gamma)")
    m.addStack(final_zzHist, "ZZ")
    m.addStack(final_otherHist, "other")
    m.addStack(wzHist, "WZ")

    #m.histsToStack = sorted( m.histsToStack, key=lambda x: x.Integral(2,-1) )
    #m.histsToStack = sorted( m.histsToStack, key=lambda x: x.Integral(1,-1) )
    #m.histsToStack = sorted( m.histsToStack, key=lambda x: x.Integral(5,-1) )
    m.histsToStack = sorted(m.histsToStack, key=lambda x: x.Integral(0, -1))

    # m.add(dataHist,"Data")

    aux.drawOpt(totSyst, "sysUnc")

    m.add(totUnc, "Total uncertainty")
    m.add(totSyst, "syst. uncertainty")
    #m.maximum = 2.6*m.getMaximum()
    m.maximum = 100.6 * m.getMaximum()
    m.minimum = m.getMinimum()
    #m.minimum = 10
    #if "final_lowEMHT" in name: m.minimum = 4e-2
    #if "final_highEMHT" in name: m.minimum = 2e-3
    # legInfo = "#it{H}_{T}^{#gamma} < 2TeV" if "lowEMHT" in name else "2TeV < #it{H}_{T}^{#gamma}"
    #if "ee" in name: legInfo += ", EE"
    # legInfo += ", |#Delta#phi|>0.3"
    #legInfo = "DiMu" if "MM" in name else "DiEle"
    legInfo = "ee+#mu#mu"
    m.leg.SetHeader(legInfo)
    # m.leg.SetY1(.56)
    # m.leg.SetX1(.56)
    # m.leg.SetX1(.46)
    # m.leg.SetX2(.99)
    # m.leg.SetX2(.89)

    # for pie chart calc
    pie = {}
    for bin in range(dirHist.GetNbinsX() - 1, dirHist.GetNbinsX() + 1):
        binName = "bin{}_{}".format(name.split("_")[1], bin)
        pie[binName] = {}
        pie[binName]["tt"] = final_ttHist.GetBinContent(bin)
        pie[binName]["dy"] = final_dyHist.GetBinContent(
            bin) if final_dyHist.GetBinContent(bin) > 0. else 0.
        pie[binName]["zz"] = final_zzHist.GetBinContent(bin)
        pie[binName]["wz"] = final_wzHist.GetBinContent(bin)
        pie[binName]["other"] = final_otherHist.GetBinContent(
            bin) if final_otherHist.GetBinContent(bin) > 0. else 0.

    # import matplotlib.pyplot as plt

    # Data to plot
    labels = 'tt', 'dy', 'zz', 'wz', 'other'
    sizes1 = [pie["binMC_1"][key] for key in pie["binMC_1"]]
    sizes2 = [pie["binMC_2"][key] for key in pie["binMC_2"]]
    # colors = ['red', 'green', 'yellow', 'blue','grey']
    # explode = (0.1, 0.1, 0.1, 0.1,0.1)

    # Plot
    # plt.rcParams['font.size'] = 20.0
    # plt.pie(sizes1, explode=explode, labels=labels, colors=colors,autopct='%1.1f%%', shadow=True, startangle=140)

    # plt.axis('equal')
    # plt.show()
    # plt.savefig("pie1.pdf")
    # plt.clf()
    # plt.pie(sizes2, explode=explode, labels=labels, colors=colors,autopct='%1.1f%%', shadow=True, startangle=140)

    # plt.show()
    # plt.axis('equal')
    # plt.savefig("pie2.pdf")

# BUILD PIE CHARTS
    # final_ttHist.SetLineColor(ROOT.kOrange+7)
    # final_zzHist.SetLineColor(ROOT.kOrange-2)
    # final_wzHist.SetLineColor(ROOT.kAzure-6)
    # final_dyHist.SetLineColor(ROOT.kGreen+3)
    # final_otherHist.SetLineColor(ROOT.kGray+2)

    c2 = TCanvas()
    pie1 = ROOT.TPie("pie1", "", 5)
    pie2 = ROOT.TPie("pie2", "", 5)
    # pie1=ROOT.TPie("pie1",r"100GeV<p_{T}^{miss}<200GeV",5)
    pie1.SetEntryVal(0, sizes1[0])
    pie1.SetEntryFillColor(0, ROOT.kOrange + 7)
    pie1.SetEntryLabel(0, r"t#bar{t}(+#gamma)")
    pie1.SetEntryVal(1, sizes1[1])
    pie1.SetEntryFillColor(1, ROOT.kGreen + 3)
    pie1.SetEntryLabel(1, r"DY/Z(+#gamma)")
    pie1.SetEntryVal(2, sizes1[2])
    pie1.SetEntryFillColor(2, ROOT.kOrange - 2)
    pie1.SetEntryLabel(2, r"ZZ")
    pie1.SetEntryVal(3, sizes1[3])
    pie1.SetEntryFillColor(3, ROOT.kAzure - 6)
    pie1.SetEntryLabel(3, r"WZ")
    pie1.SetEntryVal(4, sizes1[4])
    pie1.SetEntryFillColor(4, ROOT.kGray + 2)
    pie1.SetEntryLabel(4, r"other")

    pie2.SetEntryVal(0, sizes2[0])
    pie2.SetEntryFillColor(0, ROOT.kOrange + 7)
    pie2.SetEntryLabel(0, r"t#bar{t}(+#gamma)")
    pie2.SetEntryVal(1, sizes2[1])
    pie2.SetEntryFillColor(1, ROOT.kGreen + 3)
    pie2.SetEntryLabel(1, r"DY/Z(+#gamma)")
    pie2.SetEntryVal(2, sizes2[2])
    pie2.SetEntryFillColor(2, ROOT.kOrange - 2)
    pie2.SetEntryLabel(2, r"ZZ")
    pie2.SetEntryVal(3, sizes2[3])
    pie2.SetEntryFillColor(3, ROOT.kAzure - 6)
    pie2.SetEntryLabel(3, r"WZ")
    pie2.SetEntryVal(4, sizes2[4])
    pie2.SetEntryFillColor(4, ROOT.kGray + 2)
    pie2.SetEntryLabel(4, r"other")

    # pie1.SetLabelFormat("#scale[0.76]{#splitline{%txt}{(%perc)}}")
    pie1.SetLabelFormat("#scale[0.76]{%perc}")
    pie1.SetRadius(.3)
    pie1.SetX(0.4)
    pie1.SetY(0.35)
    pie2.SetLabelFormat("#scale[0.76]{%perc}")
    pie2.SetRadius(.3)
    pie2.SetX(0.4)
    pie2.SetY(0.35)

    pie1.Draw()

    l = ROOT.TLatex(
        0.05, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
    l2 = ROOT.TLatex(0.07, .88, "#scale[0.76]{#font[52]{Simulation}}")
    l.SetNDC()
    l2.SetNDC()
    l.Draw()
    l2.Draw()
    pieleg = pie1.MakeLegend()
    pieleg.SetHeader(r"150 GeV < p_{T}^{miss} < 200 GeV")
    c2.SaveAs("pie1.pdf")

    pie2.Draw()

    l.SetNDC()
    l2.SetNDC()
    l.Draw()
    l2.Draw()
    pieleg2 = pie1.MakeLegend()
    pieleg2.SetHeader(r"200 GeV < p_{T}^{miss}")
    c2.SaveAs("pie2.pdf")

    c = ROOT.TCanvas()
    m.Draw()

    # draw other labels

    l = ROOT.TLine()
    l.SetLineStyle(2)
    # l.SetLineColor(ROOT.kGray+2)
    # l.SetLineColor(ROOT.kRed-2)
    l.SetLineColor(ROOT.kBlack)
    text = ROOT.TLatex()
    text.SetTextSize(0.8 * text.GetTextSize())
    l.DrawLine(100, 0, 100, totUnc.GetBinContent(totUnc.FindBin(100)))
    # text.SetTextAngle(90)
    text.SetTextAngle(0)
    #text.DrawLatexNDC(.23,.315, "Normalization")
    #text.DrawLatexNDC(.22,.2, "CR")
    #text.DrawLatexNDC(.22,.4, "CR")
    # text.SetTextAngle(0)
    #text.DrawLatexNDC(.311,.315, "Validation")
    #text.DrawLatexNDC(.361,.2, "VR")
    #text.DrawLatexNDC(.681,.2, "SR")
    #text.DrawLatexNDC(.361,.4, "VR")
    #text.DrawLatexNDC(.681,.4, "SR")
    # if "final" in name:
    #l.DrawLine(150, 0, 150, totUnc.GetBinContent(totUnc.FindBin(150)))
    #l.DrawLine(100, 0, 100, totUnc.GetBinContent(totUnc.FindBin(100)))

    # r = ratio.Ratio("#scale[.9]{#lower[.24]{#splitline{Data/Pred.}{Bkg. frac.}}}", dirHist, totStat)
    r = ratio.Ratio("#scale[.9]{#lower[.24]{Bkg. frac.}}",
                    dataHist, totStat, sysHisto=totSyst)
    # r = ratio.Ratio("#scale[.9]{#lower[.24]{Bkg. frac.}}", dirHist, totStat)
    hsm = m.hists[0].GetStack().Last()
    # r = ratio.Ratio("#scale[.9]{#lower[.24]{Bkg. frac.}}", dataHist, hsm,sysHisto=totSyst)
    rMax = 2.5
    # rMax = 1.1
    #if name == "final_lowEMHT": rMax = 1.6
    #if name == "final_highEMHT": rMax = 3.6
    r.draw(0., rMax, m.getStack(), True)
    #r.draw(0., rMax, m.getStack())

    #aux.Label(sim= not dirSet==data, status="" if "allMC" not in name else "Private Work")
    #aux.Label(sim= not dirSet==dataDoubleSF, status="" if "allMC" not in name else "Private Work")
    #aux.Label(sim= not dirSet==dataDoubleSF, status="" if "MC" not in name else "Private Work")
    #aux.Label(sim= not dirSet==dataDoubleSF, status="" if "MC" not in name else "Work in Progress")
    aux.Label(sim=False, status="" if "MC" not in name else "Work in Progress")
    #aux.save(name, normal=False, changeMinMax=False)
    aux.save(name, normal=False, changeMinMax=True)

    #if name == "final_lowEMHT": dc = limitTools.MyDatacard()
    #if name == "final_MC_MM": dc = limitTools.MyDatacard()
    if name == "final_MC":
        dc = limitTools.MyDatacard()
    # elif name == "final_highEMHT": dc = limitTools.MyDatacard("testDatacard.txt")
    else:
        return
    # for bin in range(dirHist.GetNbinsX()-2, dirHist.GetNbinsX()+1):
    # for bin in range(dirHist.GetNbinsX()-1, dirHist.GetNbinsX()+1):
    # for bin in range(dirHist.GetNbinsX()-2, dirHist.GetNbinsX()+1):
    for bin in range(dirHist.GetNbinsX() - 1, dirHist.GetNbinsX() + 1):  # 2bins
        # for bin in range(dirHist.GetNbinsX()-2, dirHist.GetNbinsX()+1):#3bins
        # for bin in range(dirHist.GetNbinsX()-3, dirHist.GetNbinsX()+1):
        binName = "bin{}_{}".format(name.split("_")[1], bin)
        #print dirHist.GetBinWidth(bin)
        bw = dirHist.GetBinWidth(bin) if style.divideByBinWidth else 1.
        print "lowEdge", (signal1.GetBinLowEdge(bin))
        print "width", (signal1.GetBinWidth(bin))
        #print "signal",(signal1.GetBinContent(bin)*bw)
        #print [h.GetBinContent(bin) for h in signal1scaleHists.values()]
        #print 1.+(max([h.GetBinContent(bin) for h in signal1scaleHists.values()])-min([h.GetBinContent(bin) for h in signal1scaleHists.values()]))/2./signal1.GetBinContent(bin) if signal1.GetBinContent(bin) else 1
        #print "content",(dirHist.GetBinContent(bin))
        #print ttgHist.ProfileX().GetBinEntries(bin)
        #print "ttg",ttgHist.GetBinError(bin)
        #print "zg",zgHist.GetBinError(bin)
        #print "zz",zzHist.GetBinError(bin)
        #print "wwg",wwgHist.GetBinError(bin)
        #print "wzg",wzgHist.GetBinError(bin)
        #print "dy",dyHist.GetBinError(bin)
        #print "wjets",wjetsHist.GetBinError(bin)
        #print "tt",ttHist.GetBinError(bin)
        #print "final tt SFerr",final_ttTotSFSyst.GetBinError(bin)
        #print "final total Err",totUnc.GetBinError(bin)
        #print "singletop",singletopHist.GetBinError(bin)
        #print "wz",wzHist.GetBinError(bin)
        #print "ww",wwHist.GetBinError(bin)
        #print "wg",wgHist.GetBinError(bin)
        #print "zz4l",zz4lHist.GetBinError(bin)
        print "other", final_otherHist.GetBinContent(
            bin), final_otherHist.GetBinError(bin)

        print "stat. unc", (getDatacardUncertFromHist(
            dirHist, bin) - 1.) * 100., "%"
        #print (signal1.GetBinContent(bin-1)*bw)
        #print (signal1.GetBinContent(bin+1)*bw)
        # dc.addBin(binName, int(round(dirHist.GetBinContent(bin) * bw)),
        dc.addBin(binName, int(round(dataHist.GetBinContent(bin) * bw)),
                  {
            # [zgSyst,zgTriggerSyst,zgphotonSFSyst,zglepSFSyst,zgJESSyst,zgJERSyst,zgPUSyst,zgScaleSyst,zgSFSyst]
            # "signal": (signal1.GetBinContent(bin)-totStat.GetBinContent(bin))*bw,
            # "signal": (signal1.GetBinContent(bin)-totStat.GetBinContent(bin)-signal1_pre.GetBinContent(bin))*bw,
            "signal": (signal1.GetBinContent(bin)) * bw,
            # "ttg": ttgHist.GetBinContent(bin)*bw,
            # "zg": 0. if zgHist.GetBinContent(bin)*bw<0.000001 else zgHist.GetBinContent(bin)*bw,
            # "zz": zzHist.GetBinContent(bin)*bw,
            # "wwg": wwgHist.GetBinContent(bin)*bw,
            # "wzg": wzgHist.GetBinContent(bin)*bw,
            # "dy": dyHist.GetBinContent(bin)*bw,
            # "wjets": wjetsHist.GetBinContent(bin)*bw,
            # "tt": ttHist.GetBinContent(bin)*bw,
            # "singletop": singletopHist.GetBinContent(bin)*bw,
            # "wz": wzHist.GetBinContent(bin)*bw,
            # "ww": wwHist.GetBinContent(bin)*bw,
            # "wg": wgHist.GetBinContent(bin)*bw,
            # "zz4l": zz4lHist.GetBinContent(bin)*bw,
            "tt": final_ttHist.GetBinContent(bin) * bw,
            "dy": final_dyHist.GetBinContent(bin) * bw,
            "zz": final_zzHist.GetBinContent(bin) * bw,
            "wz": wzHist.GetBinContent(bin) * bw,
            "other": final_otherHist.GetBinContent(bin) * bw,
            #                "cont": signal1_pre.GetBinContent(bin)*bw
        }, {
            # "ttgStat_"+binName: {"tt": getDatacardUncertFromHist(aux.addHists(),bin)},
            # "ttgStat_"+binName: {"ttg": getDatacardUncertFromHist(ttgHist,bin)},
            # "zgStat_"+binName: {"zg": getDatacardUncertFromHist(zgHist,bin)},
            # "zzStat_"+binName: {"zz": getDatacardUncertFromHist(zzHist,bin)},
            # "wwgStat_"+binName: {"wwg": getDatacardUncertFromHist(wwgHist,bin)},
            # "wzgStat_"+binName: {"wzg": getDatacardUncertFromHist(wzgHist,bin)},
            # "dyStat_"+binName: {"dy": getDatacardUncertFromHist(dyHist,bin)},
            # "wjetsStat_"+binName: {"wjets": getDatacardUncertFromHist(wjetsHist,bin)},
            # "ttStat_"+binName: {"tt": getDatacardUncertFromHist(ttHist,bin)},
            # "singletopStat_"+binName: {"singletop": getDatacardUncertFromHist(singletopHist,bin)},
            # "wzStat_"+binName: {"wz": getDatacardUncertFromHist(wzHist,bin)},
            # "wwStat_"+binName: {"ww": getDatacardUncertFromHist(wwHist,bin)},
            # "wgStat_"+binName: {"wg": getDatacardUncertFromHist(wgHist,bin)},
            # "zz4lStat_"+binName: {"zz4l": getDatacardUncertFromHist(zz4lHist,bin)},
            "otherStat" + binName: {"other": getDatacardUncertFromHist(final_otherHist, bin)},
            "ttStat" + binName: {"tt": getDatacardUncertFromHist(final_ttHist, bin)},
            "dyStat" + binName: {"dy": getDatacardUncertFromHist(final_dyHist, bin)},
            "zzStat" + binName: {"zz": getDatacardUncertFromHist(final_zzHist, bin)},
            "wzStat" + binName: {"wz": getDatacardUncertFromHist(final_wzHist, bin)},
            "signalStat_" + binName: {"signal": getDatacardUncertFromHist(signal1, bin)},
            # "pdf": {
            # "wg": getDatacardUncertFromHist(wgPdfUnc,bin),
            # "zg": getDatacardUncertFromHist(zgPdfUnc,bin),
            # "tg": getDatacardUncertFromHist(tgPdfUnc,bin)},
            # "scale": {
            # "wg": getDatacardUncertFromHist(wgScaleUnc,bin),
            # "zg": getDatacardUncertFromHist(zgScaleUnc,bin),
            # "tg": getDatacardUncertFromHist(tgScaleUnc,bin)},
            "scale": {
                "signal": 1. + (max([h.GetBinContent(bin) for h in signal1scaleHists.values()]) - min([h.GetBinContent(bin) for h in signal1scaleHists.values()])) / 2. / signal1.GetBinContent(bin) if signal1.GetBinContent(bin) else 1,
                # "ttg": systPkl["ttg"]["scale"][binName]+1.,
                # "zg": systPkl["zg"]["scale"][binName]+1.,
                # "zz": systPkl["zz"]["scale"][binName]+1.,
                # "wwg": systPkl["wwg"]["scale"][binName]+1.,
                # "wzg": systPkl["wzg"]["scale"][binName]+1.,
                # "dy": systPkl["dy"]["scale"][binName]+1.,
                # "wjets": systPkl["wjets"]["scale"][binName]+1.,
                # "tt": systPkl["tt"]["scale"][binName]+1.,
                # "singletop": systPkl["singletop"]["scale"][binName]+1.,
                # "wz": systPkl["wz"]["scale"][binName]+1.,
                # "ww": systPkl["ww"]["scale"][binName]+1.,
                # "wg": systPkl["wgamma"]["scale"][binName]+1.,
                # "zz4l": systPkl["zz4l"]["scale"][binName]+1.},
                "other": systPkl["other"]["scale"][binName] + 1.,
                "wz": systPkl["wz"]["scale"][binName] + 1.,
                # "tt": ((systPkl["tt+ttgg"]["scale"][binName]*ttgHist.GetBinContent(bin)*bw+systPkl["tt"]["scale"][binName]*ttHist.GetBinContent(bin)*bw)/(ttHist.GetBinContent(bin)*bw+ttgHist.GetBinContent(bin)))+1.,
                "tt": systPkl["tt+ttg"]["scale"][binName] + 1.,
                "dy": systPkl["dy+zg"]["scale"][binName] + 1.,
                "zz": systPkl["zz+zz4l"]["scale"][binName] + 1.},
            "pdf": {
                "signal": 1.,
                # "ttg": systPkl["ttg"]["pdf"][binName]+1.,
                # "zg": systPkl["zg"]["pdf"][binName]+1.,
                # "zz": systPkl["zz"]["pdf"][binName]+1.,
                # "wwg": systPkl["wwg"]["pdf"][binName]+1.,
                # "wzg": systPkl["wzg"]["pdf"][binName]+1.,
                # "dy": systPkl["dy"]["pdf"][binName]+1.,
                # "wjets": systPkl["wjets"]["pdf"][binName]+1.,
                # "tt": systPkl["tt"]["pdf"][binName]+1.,
                # "singletop": systPkl["singletop"]["pdf"][binName]+1.,
                # "wz": systPkl["wz"]["pdf"][binName]+1.,
                # "ww": systPkl["ww"]["pdf"][binName]+1.,
                # "wg": systPkl["wgamma"]["pdf"][binName]+1.,
                # "other": systPkl["other"]["pdf"][binName]+1.,
                # "zz4l": systPkl["zz4l"]["pdf"][binName]+1.},
                "other": systPkl["other"]["pdf"][binName] + 1.,
                "wz": systPkl["wz"]["pdf"][binName] + 1.,
                "tt": systPkl["tt+ttg"]["pdf"][binName] + 1.,
                "dy": systPkl["dy+zg"]["pdf"][binName] + 1.,
                "zz": systPkl["zz+zz4l"]["pdf"][binName] + 1., },
            "JES": {
                "signal": 1. + abs(signal1JESuSFHisto.GetBinContent(bin) - signal1JESdSFHisto.GetBinContent(bin)) / 2. / signal1.GetBinContent(bin),
                # "ttg": systPkl["ttg"]["JES"][binName]+1.,
                # "zg": systPkl["zg"]["JES"][binName]+1.,
                # "zz": systPkl["zz"]["JES"][binName]+1.,
                # "wwg": systPkl["wwg"]["JES"][binName]+1.,
                # "wzg": systPkl["wzg"]["JES"][binName]+1.,
                # "dy": systPkl["dy"]["JES"][binName]+1.,
                # "wjets": systPkl["wjets"]["JES"][binName]+1.,
                # "tt": systPkl["tt"]["JES"][binName]+1.,
                # "singletop": systPkl["singletop"]["JES"][binName]+1.,
                # "wz": systPkl["wz"]["JES"][binName]+1.,
                # "ww": systPkl["ww"]["JES"][binName]+1.,
                # "wg": systPkl["wgamma"]["JES"][binName]+1.,
                # "other": systPkl["other"]["JES"][binName]+1.,
                # "zz4l": systPkl["zz4l"]["JES"][binName]+1.},
                "other": systPkl["other"]["JES"][binName] + 1.,
                "wz": systPkl["wz"]["JES"][binName] + 1.,
                "tt": systPkl["tt+ttg"]["JES"][binName] + 1.,
                "dy": systPkl["dy+zg"]["JES"][binName] + 1.,
                "zz": systPkl["zz+zz4l"]["JES"][binName] + 1., },
            "JER": {
                "signal": 1. + abs(signal1JERuSFHisto.GetBinContent(bin) - signal1JERdSFHisto.GetBinContent(bin)) / 2. / signal1.GetBinContent(bin),
                # "ttg": systPkl["ttg"]["JER"][binName]+1.,
                # "zg": systPkl["zg"]["JER"][binName]+1.,
                # "zz": systPkl["zz"]["JER"][binName]+1.,
                # "wwg": systPkl["wwg"]["JER"][binName]+1.,
                # "wzg": systPkl["wzg"]["JER"][binName]+1.,
                # "dy": systPkl["dy"]["JER"][binName]+1.,
                # "wjets": systPkl["wjets"]["JER"][binName]+1.,
                # "tt": systPkl["tt"]["JER"][binName]+1.,
                # "singletop": systPkl["singletop"]["JER"][binName]+1.,
                # "wz": systPkl["wz"]["JER"][binName]+1.,
                # "ww": systPkl["ww"]["JER"][binName]+1.,
                # "wg": systPkl["wgamma"]["JER"][binName]+1.,
                # "other": systPkl["other"]["JER"][binName]+1.,
                # "zz4l": systPkl["zz4l"]["JER"][binName]+1.},
                "other": systPkl["other"]["JER"][binName] + 1.,
                "wz": systPkl["wz"]["JER"][binName] + 1.,
                "tt": systPkl["tt+ttg"]["JER"][binName] + 1.,
                "dy": systPkl["dy+zg"]["JER"][binName] + 1.,
                "zz": systPkl["zz+zz4l"]["JER"][binName] + 1., },
            "PU": {
                "signal": 1. + abs(signal1NoPUuHisto.GetBinContent(bin) - signal1NoPUdHisto.GetBinContent(bin)) / 2 / signal1.GetBinContent(bin) if signal1.GetBinContent(bin) else 1,
                # "ttg": systPkl["ttg"]["PU"][binName]+1.,
                # "zg": systPkl["zg"]["PU"][binName]+1.,
                # "zz": systPkl["zz"]["PU"][binName]+1.,
                # "wwg": systPkl["wwg"]["PU"][binName]+1.,
                # "wzg": systPkl["wzg"]["PU"][binName]+1.,
                # "dy": systPkl["dy"]["PU"][binName]+1.,
                # "wjets": systPkl["wjets"]["PU"][binName]+1.,
                # "tt": systPkl["tt"]["PU"][binName]+1.,
                # "singletop": systPkl["singletop"]["PU"][binName]+1.,
                # "wz": systPkl["wz"]["PU"][binName]+1.,
                # "ww": systPkl["ww"]["PU"][binName]+1.,
                # "wg": systPkl["wgamma"]["PU"][binName]+1.,
                # "other": systPkl["other"]["PU"][binName]+1.,
                # "zz4l": systPkl["zz4l"]["PU"][binName]+1.},
                "other": systPkl["other"]["PU"][binName] + 1.,
                "wz": systPkl["wz"]["PU"][binName] + 1.,
                "tt": systPkl["tt+ttg"]["PU"][binName] + 1.,
                "dy": systPkl["dy+zg"]["PU"][binName] + 1.,
                "zz": systPkl["zz+zz4l"]["PU"][binName] + 1., },
            "lepSF": {
                "signal": 1. + abs(signal1lepSFuHisto.GetBinContent(bin) - signal1lepSFdHisto.GetBinContent(bin)) / 2. / signal1.GetBinContent(bin),
                # "ttg": systPkl["ttg"]["lepSF"][binName]+1.,
                # "zg": systPkl["zg"]["lepSF"][binName]+1.,
                # "zz": systPkl["zz"]["lepSF"][binName]+1.,
                # "wwg": systPkl["wwg"]["lepSF"][binName]+1.,
                # "wzg": systPkl["wzg"]["lepSF"][binName]+1.,
                # "dy": systPkl["dy"]["lepSF"][binName]+1.,
                # "wjets": systPkl["wjets"]["lepSF"][binName]+1.,
                # "tt": systPkl["tt"]["lepSF"][binName]+1.,
                # "singletop": systPkl["singletop"]["lepSF"][binName]+1.,
                # "wz": systPkl["wz"]["lepSF"][binName]+1.,
                # "ww": systPkl["ww"]["lepSF"][binName]+1.,
                # "wg": systPkl["wgamma"]["lepSF"][binName]+1.,
                # "other": systPkl["other"]["lepSF"][binName]+1.,
                # "zz4l": systPkl["zz4l"]["lepSF"][binName]+1.},
                "other": systPkl["other"]["lepSF"][binName] + 1.,
                "wz": systPkl["wz"]["lepSF"][binName] + 1.,
                "tt": systPkl["tt+ttg"]["lepSF"][binName] + 1.,
                "dy": systPkl["dy+zg"]["lepSF"][binName] + 1.,
                "zz": systPkl["zz+zz4l"]["lepSF"][binName] + 1., },
            "photonSF": {
                "signal": 1. + abs(signal1photonSFuHisto.GetBinContent(bin) - signal1photonSFdHisto.GetBinContent(bin)) / 2. / signal1.GetBinContent(bin),
                # "ttg": systPkl["ttg"]["photonSF"][binName]+1.,
                # "zg": systPkl["zg"]["photonSF"][binName]+1.,
                # "zz": systPkl["zz"]["photonSF"][binName]+1.,
                # "wwg": systPkl["wwg"]["photonSF"][binName]+1.,
                # "wzg": systPkl["wzg"]["photonSF"][binName]+1.,
                # "dy": systPkl["dy"]["photonSF"][binName]+1.,
                # "wjets": systPkl["wjets"]["photonSF"][binName]+1.,
                # "tt": systPkl["tt"]["photonSF"][binName]+1.,
                # "singletop": systPkl["singletop"]["photonSF"][binName]+1.,
                # "wz": systPkl["wz"]["photonSF"][binName]+1.,
                # "ww": systPkl["ww"]["photonSF"][binName]+1.,
                # "wg": systPkl["wgamma"]["photonSF"][binName]+1.,
                # "other": systPkl["other"]["photonSF"][binName]+1.,
                # "zz4l": systPkl["zz4l"]["photonSF"][binName]+1.},
                "other": systPkl["other"]["photonSF"][binName] + 1.,
                "wz": systPkl["wz"]["photonSF"][binName] + 1.,
                "tt": systPkl["tt+ttg"]["photonSF"][binName] + 1.,
                "dy": systPkl["dy+zg"]["photonSF"][binName] + 1.,
                "zz": systPkl["zz+zz4l"]["photonSF"][binName] + 1., },
            "EWK": {
                "signal": 1. + abs(signal1EWKuHisto.GetBinContent(bin) - signal1EWKdHisto.GetBinContent(bin)) / 2. / signal1.GetBinContent(bin),
                # "ttg": systPkl["ttg"]["EWK"][binName]+1.,
                # "zg": systPkl["zg"]["EWK"][binName]+1.,
                # "zz": systPkl["zz"]["EWK"][binName]+1.,
                # "wwg": systPkl["wwg"]["EWK"][binName]+1.,
                # "wzg": systPkl["wzg"]["EWK"][binName]+1.,
                # "dy": systPkl["dy"]["EWK"][binName]+1.,
                # "wjets": systPkl["wjets"]["EWK"][binName]+1.,
                # "tt": systPkl["tt"]["EWK"][binName]+1.,
                # "singletop": systPkl["singletop"]["EWK"][binName]+1.,
                # "wz": systPkl["wz"]["EWK"][binName]+1.,
                # "ww": systPkl["ww"]["EWK"][binName]+1.,
                # "wg": systPkl["wgamma"]["EWK"][binName]+1.,
                # "other": systPkl["other"]["EWK"][binName]+1.,
                # "zz4l": systPkl["zz4l"]["EWK"][binName]+1.},
                "other": systPkl["other"]["EWK"][binName] + 1.,
                "wz": systPkl["wz"]["EWK"][binName] + 1.,
                "tt": systPkl["tt+ttg"]["EWK"][binName] + 1.,
                "dy": systPkl["dy+zg"]["EWK"][binName] + 1.,
                "zz": systPkl["zz+zz4l"]["EWK"][binName] + 1., },
            "ISR": {
                "signal": 1.001 + abs(signal1ISRuHisto.GetBinContent(bin) - signal1ISRdHisto.GetBinContent(bin)) / 2. / signal1.GetBinContent(bin),
                # "ttg": systPkl["ttg"]["ISR"][binName]+1.,
                # "zg": systPkl["zg"]["ISR"][binName]+1.,
                # "zz": systPkl["zz"]["ISR"][binName]+1.,
                # "wwg": systPkl["wwg"]["ISR"][binName]+1.,
                # "wzg": systPkl["wzg"]["ISR"][binName]+1.,
                # "dy": systPkl["dy"]["ISR"][binName]+1.,
                # "wjets": systPkl["wjets"]["ISR"][binName]+1.,
                # "tt": systPkl["tt"]["ISR"][binName]+1.,
                # "singletop": systPkl["singletop"]["ISR"][binName]+1.,
                # "wz": systPkl["wz"]["ISR"][binName]+1.,
                # "ww": systPkl["ww"]["ISR"][binName]+1.,
                # "wg": systPkl["wgamma"]["ISR"][binName]+1.,
                # "other": systPkl["other"]["ISR"][binName]+1.,
                # "zz4l": systPkl["zz4l"]["ISR"][binName]+1.},
                "other": systPkl["other"]["ISR"][binName] + 1.,
                "wz": systPkl["wz"]["ISR"][binName] + 1.,
                "tt": systPkl["tt+ttg"]["ISR"][binName] + 1.,
                "dy": systPkl["dy+zg"]["ISR"][binName] + 1.,
                "zz": systPkl["zz+zz4l"]["ISR"][binName] + 1., },
            "lumi": {
                "signal": 1.026,
                # "ttg": 1.026,
                # "zg": 1.026,
                # "zz": 1.026,
                # "wwg": 1.026,
                # "wzg": 1.026,
                # "dy": 1.026,
                # "wjets": 1.026,
                # "tt": 1.026,
                # "singletop": 1.026,
                # "wz": 1.026,
                # "ww": 1.026,
                # "wg": 1.026,
                "other": 1.026,
                "tt": 1.026,
                "dy": 1.026,
                "wz": 1.026,
                # "zz4l": 1.026},
                "zz": 1.026},
            "genmet": {
                "signal": 1. + abs(signal1.GetBinContent(bin) - signal1GenMetHisto.GetBinContent(bin)) / 2 / signal1.GetBinContent(bin) if signal1.GetBinContent(bin) else 1},
            # "ttg": 1.,
            # "zg": 1.,
            # "zz": 1.,
            # "wwg": 1.,
            # "wzg": 1.,
            # "dy": 1.,
            # "wjets": 1.,
            # "tt": 1.,
            # "singletop": 1.,
            # "wz": 1.,
            # "ww": 1.,
            # "wg": 1.,
            # "other": 1.,
            # "zz4l": 1.},
            # "zz": 1.},
            # "pu": {
            # "wg": getDatacardUncertFromHist(wgPuUnc,bin),
            # "zg": getDatacardUncertFromHist(zgPuUnc,bin),
            # "tg": getDatacardUncertFromHist(tgPuUnc,bin)},
            # "jes": {
            # "wg": getDatacardUncertFromHist(wgJesUnc,bin),
            # "zg": getDatacardUncertFromHist(zgJesUnc,bin),
            # "tg": getDatacardUncertFromHist(tgJesUnc,bin)},
            "SFDY": {
                # "zg": (1.+sfDYErr),
                # "zg": getDatacardUncertFromHist(zgSFSyst,bin),
                "dy": (1. + sfDYErr / sfDY)},
            "SFTT": {
                # "tt": getDatacardUncertFromHist(),
                # "ttg": getDatacardUncertFromHist()},
                # "tt": (1.+0.04),
                # "ttg": (1.+0.04)},
                # "tt": getDatacardUncertFromHist(ttTotSFSyst,bin),
                # "ttg": getDatacardUncertFromHist(ttgTotSFSyst,bin)},
                "tt": getDatacardUncertFromHist(final_ttTotSFSyst, bin)},
            "SFZZ": {
                # "zz": (1.+sfZZErr),
                # "zz4l": (1.+sfZZErr)},
                "zz": (1. + sfZZErr / sfZZ)},
            "SFWZ": {
                "wz": (1. + sfWZErr / sfWZ)},
            # "xSecWWG": {
            # "wwg": (1.5)},
            # "xSecWZG": {
            # "wzg": (1.5)},
            # "xSecWJets": {
            # "wjets": (1.5)},
            # "xSecSingletop": {
            # "singletop": (1.5)},
            # "xSecWW": {
            # "ww": (1.5)},
            # "xSecWG": {
            # "wg": (1.5)},
            "xSecOTHER": {
                "other": (1.5)},
            # "dataMC": {
            # "signal": 1.05,
            # "ttg": 1.2,
            # "zg": 1.2,
            # "zz": 1.2,
            # "wwg": 1.2,
            # "wzg": 1.2,
            # "dy": 1.2,
            # "wjets": 1.2,
            # "tt": 1.2,
            # "singletop": 1.2,
            # "wz": 1.2,
            # "ww": 1.2,
            # "wg": 1.2,
            # "zz4l": 1.2},
            "trigger": {
                "signal": 1.03,
                # "ttg": 1.03,
                # "zg": 1.03,
                "zz": 1.03,
                # "wwg": 1.03,
                # "wzg": 1.03,
                "dy": 1.03,
                # "wjets": 1.03,
                "tt": 1.03,
                # "singletop": 1.03,
                "wz": 1.03,
                # "ww": 1.03,
                # "wg": 1.03,
                # "other": 1.03,
                "other": 1.03}
            # "zz4l": 1.03},
            # "isr": {"signal": 1.001},
            # "genMet": {"signal": 1.001},
            # jes, jer splitting
        }
        )
    dc.write("testDatacard.txt")
    dc.write("limitCalculations/testDatacard.txt")


def main():
    #allMC = gjets+qcd+zg+wg+ttg+wjets+ttjets_nlo+znunu
    allMC = zgamma + ttgamma + zz + wwgamma + wzgamma + \
        DYjetsNLO + wjets + tt + singletop + wz + ww + zz4l + wjets
    allMC.label = "MC mix"
    #finalDistributionSignalHist("final_lowEMHT", data, "signal_lowEMHT", dataHt, data, "signal_lowEMHT_eControl")
    #finalDistributionSignalHist("allMC_lowEMHT", allMc, "signal_lowEMHT", allMc, allMc, "signal_lowEMHT_eControl")
    #finalDistributionSignalHist("final_MC_MM", allMC, "onZMM")
    #finalDistributionSignalHist("final_MC_EE", allMC, "onZEE")
    #finalDistributionSignalHist("final_MC", allMC, "onZMet100")
    #finalDistributionSignalHist("final_MC", allMC, "onZMet150")
    #finalDistributionSignalHist("final_MC", allMC, "onZ/LL")
    #finalDistributionSignalHist("final_MC", allMC, "onZMet150/LL")
    finalDistributionSignalHist("final_MC", allMC, "xx_0_0/sig/LL/nom")


main()
