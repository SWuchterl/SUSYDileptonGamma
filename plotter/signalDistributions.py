import ROOT
from ROOT import *
from array import array
from include import *
import numpy as np
import pickle as pkl

def getDatacardUncertFromHist(h,b):
    c = h.GetBinContent(b)
    return 1.+max(0,h.GetBinError(b)/c if c else 0)


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
    #according to Knuts Tool ->250
    nBins = [0,25,50,75,100,150,250,450]
    #nBinsData = [0,25,50,75,100]????
    nBinsData = [0,25,50,75,100,150]
    
    #nBins = [200, 300,500]
    #nBins = np.concatenate((np.arange(0,200,50),np.arange(200,300,100),np.arange(300,3000,1400)),axis=0)
    #nBins = np.concatenate((np.arange(0,200,25),np.arange(200,300,100),np.arange(300,600,200)),axis=0)
    #nBins = np.concatenate((np.arange(0,200,25),np.arange(200,500,300)),axis=0)
    #print nBins

    # direct stuff
    
    dirHist = aux.stdHist(dirSet, dirDir+"/met", nBins)
    
    #dataHist = aux.stdHist(dataDoubleSF,dirDir+"/met",nBins)
    dataHist = aux.stdHist(dataLL,dirDir+"/met",nBins)
    dataHist = aux.rebin(dataHist,nBinsData)
    dataHist.SetYTitle(aux.getYAxisTitle(dataHist))

    style.additionalPoissonUncertainty = False
    aux.drawOpt(dirHist, "data")


    aux.drawOpt(dataHist, "data")
    #gjetHist, gjetSyst, info = gjetPrediction(dirHist, preSet, zg+wg+ttg+wjets+ttjets_nlo+znunu, "met", nBins, weight, name+"_divByBinWidth" if style.divideByBinWidth else name)
    #gjetHist.SetLineColor(rwth.myLightBlue)
    #gjetHist.GetXaxis().SetTitle("#it{p}_{T}^{miss} (GeV)")


    #eHist = aux.stdHist(preSetElectron, preDirElectron+"/met", nBins)
    #eHist.GetXaxis().SetTitle("#it{p}_{T}^{miss} (GeV)")
    #eHist.Scale( 0.0267 if dirSet == data else 0.0154 )
    #eHist.SetLineColor(rwth.myYellow)
    #eSyst = aux.getSysHisto(eHist, 0.3)

#zgamma+ttgamma+zz+wwgamma+wzgamma+DYjetsNLO+wjets+tt+singletop+wz+ww+zz4l+wjets

    zgHist = aux.stdHist(zgamma, dirDir+"/met", nBins)
    ttgHist = aux.stdHist(ttgamma, dirDir+"/met", nBins)
    zzHist = aux.stdHist(zz, dirDir+"/met", nBins)
    wwgHist = aux.stdHist(wwgamma, dirDir+"/met", nBins)
    wzgHist = aux.stdHist(wzgamma, dirDir+"/met", nBins)
    dyHist = aux.stdHist(DYjetsNLO, dirDir+"/met", nBins)
    wjetsHist = aux.stdHist(wjets, dirDir+"/met", nBins)
    ttHist = aux.stdHist(tt, dirDir+"/met", nBins)
    #ttg080Hist = aux.stdHist(ttgamma, name.replace("VR","VR080"), binning)
    #ttg80Hist = aux.stdHist(ttgamma, name.replace("VR","VR80"), binning)
    singletopHist=aux.stdHist(singletop, dirDir+"/met", nBins)
    wzHist=aux.stdHist(wz, dirDir+"/met", nBins)
    wwHist=aux.stdHist(ww, dirDir+"/met", nBins)
    zz4lHist=aux.stdHist(zz4l, dirDir+"/met", nBins)
    wgHist=aux.stdHist(wgamma,dirDir+"/met",nBins)
    
    #Scaling
    zg_AvgTopPtWeightHisto = zgamma.getHist(dirDir+"/weight_topPt")
    zg_AvgNIsrWeightHisto = zgamma.getHist(dirDir+"/weight_nISR")
    zg_AvgEWKinoWeightHisto = zgamma.getHist(dirDir+"/weight_EWKinoPairPt")
    #zg_AvgleptonWeightHisto = zgamma.getHist(dirDir+"/weight_leptonPairPt")
    ttg_AvgTopPtWeightHisto = ttgamma.getHist(dirDir+"/weight_topPt")
    ttg_AvgNIsrWeightHisto = ttgamma.getHist(dirDir+"/weight_nISR")
    ttg_AvgEWKinoWeightHisto = ttgamma.getHist(dirDir+"/weight_EWKinoPairPt")
    #ttg_AvgleptonWeightHisto = ttgamma.getHist(dirDir+"/weight_leptonPairPt")
    zz_AvgTopPtWeightHisto = zz.getHist(dirDir+"/weight_topPt")
    zz_AvgNIsrWeightHisto = zz.getHist(dirDir+"/weight_nISR")
    zz_AvgEWKinoWeightHisto = zz.getHist(dirDir+"/weight_EWKinoPairPt")
    #zz_AvgleptonWeightHisto = zz.getHist(dirDir+"/weight_leptonPairPt")
    wwg_AvgTopPtWeightHisto = wwgamma.getHist(dirDir+"/weight_topPt")
    wwg_AvgNIsrWeightHisto = wwgamma.getHist(dirDir+"/weight_nISR")
    wwg_AvgEWKinoWeightHisto = wwgamma.getHist(dirDir+"/weight_EWKinoPairPt")
    #wwg_AvgleptonWeightHisto = wwgamma.getHist(dirDir+"/weight_leptonPairPt")
    wzg_AvgTopPtWeightHisto = wzgamma.getHist(dirDir+"/weight_topPt")
    wzg_AvgNIsrWeightHisto = wzgamma.getHist(dirDir+"/weight_nISR")
    wzg_AvgEWKinoWeightHisto = wzgamma.getHist(dirDir+"/weight_EWKinoPairPt")
    #wzg_AvgleptonWeightHisto = wzgamma.getHist(dirDir+"/weight_leptonPairPt")
    dy_AvgTopPtWeightHisto = DYjetsNLO.getHist(dirDir+"/weight_topPt")
    dy_AvgNIsrWeightHisto = DYjetsNLO.getHist(dirDir+"/weight_nISR")
    dy_AvgEWKinoWeightHisto = DYjetsNLO.getHist(dirDir+"/weight_EWKinoPairPt")
    #dy_AvgleptonWeightHisto = DYjetsNLO.getHist(dirDir+"/weight_leptonPairPt")
    wjets_AvgTopPtWeightHisto = wjets.getHist(dirDir+"/weight_topPt")
    wjets_AvgNIsrWeightHisto = wjets.getHist(dirDir+"/weight_nISR")
    wjets_AvgEWKinoWeightHisto = wjets.getHist(dirDir+"/weight_EWKinoPairPt")
    #wjets_AvgleptonWeightHisto = wjets.getHist(dirDir+"/weight_leptonPairPt")
    tt_AvgTopPtWeightHisto = tt.getHist(dirDir+"/weight_topPt")
    tt_AvgNIsrWeightHisto = tt.getHist(dirDir+"/weight_nISR")
    tt_AvgEWKinoWeightHisto = tt.getHist(dirDir+"/weight_EWKinoPairPt")
    #tt_AvgleptonWeightHisto = tt.getHist(dirDir+"/weight_leptonPairPt")
    singletop_AvgTopPtWeightHisto = singletop.getHist(dirDir+"/weight_topPt")
    singletop_AvgNIsrWeightHisto = singletop.getHist(dirDir+"/weight_nISR")
    singletop_AvgEWKinoWeightHisto = singletop.getHist(dirDir+"/weight_EWKinoPairPt")
    #singletop_AvgleptonWeightHisto = singletop.getHist(dirDir+"/weight_leptonPairPt")
    wz_AvgTopPtWeightHisto = wz.getHist(dirDir+"/weight_topPt")
    wz_AvgNIsrWeightHisto = wz.getHist(dirDir+"/weight_nISR")
    wz_AvgEWKinoWeightHisto = wz.getHist(dirDir+"/weight_EWKinoPairPt")
    #wz_AvgleptonWeightHisto = wz.getHist(dirDir+"/weight_leptonPairPt")
    ww_AvgTopPtWeightHisto = ww.getHist(dirDir+"/weight_topPt")
    ww_AvgNIsrWeightHisto = ww.getHist(dirDir+"/weight_nISR")
    ww_AvgEWKinoWeightHisto = ww.getHist(dirDir+"/weight_EWKinoPairPt")
    #ww_AvgleptonWeightHisto = ww.getHist(dirDir+"/weight_leptonPairPt")
    zz4l_AvgTopPtWeightHisto = zz4l.getHist(dirDir+"/weight_topPt")
    zz4l_AvgNIsrWeightHisto = zz4l.getHist(dirDir+"/weight_nISR")
    zz4l_AvgEWKinoWeightHisto = zz4l.getHist(dirDir+"/weight_EWKinoPairPt")
    #zz4l_AvgleptonWeightHisto = zz4l.getHist(dirDir+"/weight_leptonPairPt")
    wg_AvgTopPtWeightHisto = wgamma.getHist(dirDir+"/weight_topPt")
    wg_AvgNIsrWeightHisto = wgamma.getHist(dirDir+"/weight_nISR")
    wg_AvgEWKinoWeightHisto = wgamma.getHist(dirDir+"/weight_EWKinoPairPt")
    #wg_AvgleptonWeightHisto = wgamma.getHist(dirDir+"/weight_leptonPairPt")
    
    histsToScale=[zgHist,ttgHist,zzHist,wwgHist,wzgHist,dyHist,wjetsHist,ttHist,singletopHist,wzHist,wwHist,zz4lHist,wgHist]
    topWeightHists=[zg_AvgTopPtWeightHisto,ttg_AvgTopPtWeightHisto,zz_AvgTopPtWeightHisto,wwg_AvgTopPtWeightHisto,wzg_AvgTopPtWeightHisto,dy_AvgTopPtWeightHisto,wjets_AvgTopPtWeightHisto,tt_AvgTopPtWeightHisto,singletop_AvgTopPtWeightHisto,wz_AvgTopPtWeightHisto,ww_AvgTopPtWeightHisto,zz4l_AvgTopPtWeightHisto,wg_AvgTopPtWeightHisto]
    nISRWeightHists=[zg_AvgNIsrWeightHisto,ttg_AvgNIsrWeightHisto,zz_AvgNIsrWeightHisto,wwg_AvgNIsrWeightHisto,wzg_AvgNIsrWeightHisto,dy_AvgNIsrWeightHisto,wjets_AvgNIsrWeightHisto,tt_AvgNIsrWeightHisto,singletop_AvgNIsrWeightHisto,wz_AvgNIsrWeightHisto,ww_AvgNIsrWeightHisto,zz4l_AvgNIsrWeightHisto,wg_AvgNIsrWeightHisto]
    EWKinoWeightHists=[zg_AvgEWKinoWeightHisto,ttg_AvgEWKinoWeightHisto,zz_AvgEWKinoWeightHisto,wwg_AvgEWKinoWeightHisto,wzg_AvgEWKinoWeightHisto,dy_AvgEWKinoWeightHisto,wjets_AvgEWKinoWeightHisto,tt_AvgEWKinoWeightHisto,singletop_AvgEWKinoWeightHisto,wz_AvgEWKinoWeightHisto,ww_AvgEWKinoWeightHisto,zz4l_AvgEWKinoWeightHisto,wg_AvgEWKinoWeightHisto]
    #leptonWeightHists=[zg_AvgleptonWeightHisto,ttg_AvgleptonWeightHisto,zz_AvgleptonWeightHisto,wwg_AvgleptonWeightHisto,wzg_AvgleptonWeightHisto,dy_AvgleptonWeightHisto,wjets_AvgleptonWeightHisto,tt_AvgleptonWeightHisto,singletop_AvgleptonWeightHisto,wz_AvgleptonWeightHisto,ww_AvgleptonWeightHisto,zz4l_AvgleptonWeightHisto,wg_AvgleptonWeightHisto]
      
    #print len(histsToScale),len(topWeightHists),len(nISRWeightHists),len(EWKinoWeightHists),len(leptonWeightHists)
    
    for i in range(len(histsToScale)):
        if topWeightHists[i].Integral>0.:
            histsToScale[i].Scale(topWeightHists[i].GetMean())
        if nISRWeightHists[i].Integral>0.:
            histsToScale[i].Scale(nISRWeightHists[i].GetMean())
        if EWKinoWeightHists[i].Integral>0.:
            histsToScale[i].Scale(EWKinoWeightHists[i].GetMean())
        #if leptonWeightHists[i].Integral>0.:
            #histsToScale[i].Scale(leptonWeightHists[i].GetMean())
            

    pklZZ = pkl.load( open( "plots_CR_zz/factors/CRZZ.pkl", "rb" ) )
    pklDY = pkl.load( open( "plots_CR_dy/factors/CRDY.pkl", "rb" ) )
    pklTT = pkl.load( open( "plots_CR_tt/factors/CRTT.pkl", "rb" ) )
    pklWZ = pkl.load( open( "plots_CR_wz/factors/CRWZ.pkl", "rb" ) )
    sfZZ,sfZZErr = pklZZ["LL"]["m_ll"]
    sfDY,sfDYErr = pklDY["LL"]["eta1"]
    sfTT = pklTT["EM"]["eta1"][0]
    sfTT080Err=0.04
    sfTT80Err=0.4
    sfWZ,sfWZErr = pklWZ["LL"]["eta1"]


    zzHist.Scale(sfZZ)
    zz4lHist.Scale(sfZZ)
    dyHist.Scale(sfDY)
    zgHist.Scale(sfDY)
    wzHist.Scale(sfWZ)
    ttHist.Scale(sfTT)
    ttgHist.Scale(sfTT)
    
    dirHist = aux.addHists(*histsToScale)



    zgHist.SetLineColor(ROOT.kGreen-3)
    ttgHist.SetLineColor(ROOT.kRed+1)
    zzHist.SetLineColor(ROOT.kYellow)
    wwgHist.SetLineColor(ROOT.kCyan+3)
    wzgHist.SetLineColor(ROOT.kBlue+3)
    dyHist.SetLineColor(ROOT.kGreen+3)
    wjetsHist.SetLineColor(ROOT.kBlue-9)
    ttHist.SetLineColor(ROOT.kOrange+8)
    singletopHist.SetLineColor(ROOT.kOrange+4)
    wzHist.SetLineColor(ROOT.kAzure-5)
    wwHist.SetLineColor(ROOT.kCyan-3)
    zz4lHist.SetLineColor(ROOT.kOrange-2)
    wgHist.SetLineColor(ROOT.kRed+3)

    #zgPdfUnc = pdfUncertainty(zg+znunu, dirDir, nBins)
    #wgPdfUnc = pdfUncertainty(wg+wjets, dirDir, nBins)
    #tgPdfUnc = pdfUncertainty(ttjets_nlo+ttg, dirDir, nBins)

    #zgPuUnc = puUncertainty(zg+znunu, dirDir, nBins)
    #wgPuUnc = puUncertainty(wg+wjets, dirDir, nBins)
    #tgPuUnc = puUncertainty(ttjets_nlo+ttg, dirDir, nBins)

    #zgScaleUnc = scaleUncertainty(zg+znunu, dirDir, nBins)
    #wgScaleUnc = scaleUncertainty(wg+wjets, dirDir, nBins)
    #tgScaleUnc = scaleUncertainty(ttjets_nlo+ttg, dirDir, nBins)

    #zgJesUnc = jecUncertainty(zg+znunu, dirDir, nBins)
    #wgJesUnc = jecUncertainty(wg+wjets, dirDir, nBins)
    #tgJesUnc = jecUncertainty(ttjets_nlo+ttg, dirDir, nBins)

    #mcSystUncert = 0.04 # SF, lumi, trigger
    mcSystUncert = 0.2 # SF, lumi, trigger
    
    #ttgSFUnc = 
    
    zgSyst = aux.getSysHisto(zgHist, mcSystUncert)
    ttgSyst = aux.getSysHisto(ttgHist, mcSystUncert)
    zzSyst = aux.getSysHisto(zzHist, mcSystUncert)
    wwgSyst = aux.getSysHisto(wwgHist, mcSystUncert)
    wzgSyst = aux.getSysHisto(wzgHist, mcSystUncert)
    dySyst = aux.getSysHisto(dyHist, mcSystUncert)
    wjetsSyst = aux.getSysHisto(wjetsHist, mcSystUncert)
    ttSyst = aux.getSysHisto(ttHist, mcSystUncert)
    singletopSyst=aux.getSysHisto(singletopHist,mcSystUncert)
    wzSyst=aux.getSysHisto(wzHist,mcSystUncert)
    wwSyst=aux.getSysHisto(wwHist,mcSystUncert)
    zz4lSyst=aux.getSysHisto(zz4lHist,mcSystUncert)
    wgSyst=aux.getSysHisto(wgHist,mcSystUncert)

    #dirSyst= aux.getSysHisto(dirHist,mcSystUncert)
    dirSyst= aux.getSysHisto(aux.addHists(zgSyst,ttgSyst,zzSyst,wwgSyst,wzgSyst,dySyst,wjetsSyst,ttSyst,singletopSyst,wzSyst,wwSyst,zz4lSyst,wgSyst),mcSystUncert)

    #zgSyst = aux.addUncertaintiesQuadratic([zgSyst,zgPdfUnc,zgScaleUnc,zgJesUnc,zgPuUnc])
    #wgSyst = aux.addUncertaintiesQuadratic([wgSyst,wgPdfUnc,wgScaleUnc,wgJesUnc,wgPuUnc])
    #tgSyst = aux.addUncertaintiesQuadratic([tgSyst,tgPdfUnc,tgScaleUnc,tgJesUnc,tgPuUnc])

    #totStat = aux.addHists(gjetHist, eHist, zgHist, wgHist, tgHist)
    #totSyst = aux.addHists(gjetSyst, eSyst, zgSyst, wgSyst, tgSyst)
    #totStat = aux.addHists(zgHist, ttgHist, zzHist, wwgHist, wzgHist,dyHist,wjetsHist,ttHist)
    #totSyst = aux.addHists(zgSyst, ttgSyst, zzSyst, wwgSyst, wzgSyst,dySyst,wjetsSyst,ttSyst)
    totStat = aux.addHists(zgHist, ttgHist, zzHist, wwgHist, wzgHist,dyHist,wjetsHist,ttHist,singletopHist,wzHist,wwHist,zz4lHist,wjetsHist,wgHist)
    totSyst = aux.addHists(zgSyst, ttgSyst, zzSyst, wwgSyst, wzgSyst,dySyst,wjetsSyst,ttSyst,singletopSyst,wzSyst,wwSyst,zz4lSyst,wjetsSyst,wgSyst)

    #signal1 = aux.stdHist(t5bbbbzg_1500_400, dirDir+"/met", nBins)
    #signal2 = aux.stdHist(tching_400, dirDir+"/met", nBins)
    #signal2 = aux.stdHist(t5bbbbzg_1500_400, dirDir+"/met", nBins)
    signal2 = aux.stdHist(gmsb_290_205, dirDir+"/met", nBins)
    #signal1 = aux.stdHist(tching_400, dirDir+"/met", nBins)
    signal1 = aux.stdHist(tching_600, dirDir+"/met", nBins)
    
    signal2_AvgTopPtWeightHisto = t5bbbbzg_1500_400.getHist(dirDir+"/weight_topPt")
    signal2_AvgNIsrWeightHisto = t5bbbbzg_1500_400.getHist(dirDir+"/weight_nISR")
    signal2_AvgEWKinoWeightHisto = t5bbbbzg_1500_400.getHist(dirDir+"/weight_EWKinoPairPt")
    #signal2_AvgleptonWeightHisto = t5bbbbzg_1500_400.getHist(dirDir+"/weight_leptonPairPt")
    signal1_AvgTopPtWeightHisto = tching_600.getHist(dirDir+"/weight_topPt")
    signal1_AvgNIsrWeightHisto = tching_600.getHist(dirDir+"/weight_nISR")
    signal1_AvgEWKinoWeightHisto = tching_600.getHist(dirDir+"/weight_EWKinoPairPt")
    #signal1_AvgleptonWeightHisto = tching_600.getHist(dirDir+"/weight_leptonPairPt")

#
    signal2.Scale((1./signal2_AvgTopPtWeightHisto.GetMean()))
    signal2.Scale((1./signal2_AvgNIsrWeightHisto.GetMean()))
    signal2.Scale((1./signal2_AvgEWKinoWeightHisto.GetMean()))
    #signal2.Scale((1./signal2_AvgleptonWeightHisto.GetMean()))
    signal1.Scale((1./signal1_AvgTopPtWeightHisto.GetMean()))
    signal1.Scale((1./signal1_AvgNIsrWeightHisto.GetMean()))
    signal1.Scale((1./signal1_AvgEWKinoWeightHisto.GetMean()))
    #signal1.Scale((1./signal1_AvgleptonWeightHisto.GetMean()))



    
    #print signal1.GetEntries()
    
    for h in signal1, signal2:
        aux.drawOpt(h, "signal")
    #    h.Add(totStat)
    signal1.SetLineColor(ROOT.kBlue+3)
    #signal2.SetLineColor(ROOT.kRed-3)
    signal2.SetLineColor(ROOT.kBlue+3)
    signal2.SetLineStyle(2)

    #signal1_pre = aux.createHistoFromDatasetTree(t5wg_1600_100, "met*{}".format(info["shift"]), weight, nBins, "tr_jControl/simpleTree")
    #signal1_pre.Scale(info["scale"])

    #totSyst = gjetSyst

    #totUnc = aux.addHistUncert(totStat, totSyst)
    totUnc = aux.addHistUncert(totStat, dirSyst)
    #totUnc = aux.addHistUncert(dirSyst)
    aux.drawOpt(totUnc, "totUnc")
    #totUnc = dirSyst
    #aux.drawOpt(totUnc, "totUnc")
    
    c = ROOT.TCanvas()
    m = multiplot.Multiplot()
    #if dirSet == data:
    
    #if dirSet == dataDoubleSF:
        #m.add(dirHist, "Data")
    #else:
        #m.add(dirHist, "Direct simulation")
        
        
#   m.add(signal1_pre, "contamination")
    #m.add(signal1, "T5bbbbZg")
    #m.add(signal2, "TChiNg")
    #m.add(signal2, "T5bbbbZg")
    m.add(signal2, "GMSB")
    m.add(signal1, "TChiNg")
    #m.addStack(eHist, "e#rightarrow#gamma")
    m.addStack(zgHist, "Z#gamma")
    m.addStack(ttgHist, "t#bar{t}#gamma")
    m.addStack(zzHist, "ZZ(#rightarrowll#nu#nu)")
    m.addStack(wwgHist, "WW#gamma")
    m.addStack(wzgHist, "WZ#gamma")
    m.addStack(dyHist, "Drell-Yan/Z")
    m.addStack(wjetsHist, "W+jets")
    m.addStack(ttHist, "t#bar{t}")
    m.addStack(singletopHist, "single t")
    m.addStack(wwHist, "WW")
    m.addStack(wzHist, "WZ")
    m.addStack(zz4lHist, "ZZ(#rightarrow4l)")
    m.addStack(wgHist, "W#gamma")
    
    
    
    #m.histsToStack = sorted( m.histsToStack, key=lambda x: x.Integral(2,-1) )
    #m.histsToStack = sorted( m.histsToStack, key=lambda x: x.Integral(1,-1) )
    m.histsToStack = sorted( m.histsToStack, key=lambda x: x.Integral(5,-1) )


    m.add(dataHist,"Data")

    m.add(totUnc, "Total uncertainty")
    m.maximum = 2.6*m.getMaximum()
    m.minimum = m.getMinimum()
    #if "final_lowEMHT" in name: m.minimum = 4e-2
    #if "final_highEMHT" in name: m.minimum = 2e-3
    #legInfo = "#it{H}_{T}^{#gamma} < 2TeV" if "lowEMHT" in name else "2TeV < #it{H}_{T}^{#gamma}"
    #if "ee" in name: legInfo += ", EE"
    #legInfo += ", |#Delta#phi|>0.3"
    #legInfo = "DiMu" if "MM" in name else "DiEle"
    legInfo = "ee+#mu#mu"
    m.leg.SetHeader(legInfo)
    #m.leg.SetY1(.56)
    #m.leg.SetX1(.56)
    #m.leg.SetX1(.46)
    #m.leg.SetX2(.99)
    #m.leg.SetX2(.89)


    m.Draw()

    # draw other labels

    l = ROOT.TLine()
    l.SetLineStyle(2)
    #l.SetLineColor(ROOT.kGray+2)
    #l.SetLineColor(ROOT.kRed-2)
    l.SetLineColor(ROOT.kBlack)
    text = ROOT.TLatex()
    text.SetTextSize(0.8*text.GetTextSize())
    l.DrawLine(100, 0, 100, totUnc.GetBinContent(totUnc.FindBin(100)))
    #text.SetTextAngle(90)
    text.SetTextAngle(0)
    #text.DrawLatexNDC(.23,.315, "Normalization")
    #text.DrawLatexNDC(.22,.2, "CR")
    text.DrawLatexNDC(.22,.4, "CR")
    #text.SetTextAngle(0)
    #text.DrawLatexNDC(.311,.315, "Validation")
    #text.DrawLatexNDC(.361,.2, "VR")
    #text.DrawLatexNDC(.681,.2, "SR")
    text.DrawLatexNDC(.361,.4, "VR")
    text.DrawLatexNDC(.681,.4, "SR")
    if "final" in name:
        l.DrawLine(150, 0, 150, totUnc.GetBinContent(totUnc.FindBin(150)))
        #l.DrawLine(100, 0, 100, totUnc.GetBinContent(totUnc.FindBin(100)))



    #r = ratio.Ratio("#scale[.9]{#lower[.24]{#splitline{Data/Pred.}{Bkg. frac.}}}", dirHist, totStat)
    r = ratio.Ratio("#scale[.9]{#lower[.24]{Bkg. frac.}}", dirHist, totStat)
    hsm = m.hists[0].GetStack().Last()
    #r = ratio.Ratio("#scale[.9]{#lower[.24]{Bkg. frac.}}", dataHist, hsm,sysHisto=totSyst)
    #rMax = 1.5
    rMax = 1.1
    #if name == "final_lowEMHT": rMax = 1.6
    #if name == "final_highEMHT": rMax = 3.6
    r.draw(0., rMax, m.getStack(), True)
    #r.draw(0., rMax, m.getStack())

    #aux.Label(sim= not dirSet==data, status="" if "allMC" not in name else "Private Work")
    #aux.Label(sim= not dirSet==dataDoubleSF, status="" if "allMC" not in name else "Private Work")
    #aux.Label(sim= not dirSet==dataDoubleSF, status="" if "MC" not in name else "Private Work")
    #aux.Label(sim= not dirSet==dataDoubleSF, status="" if "MC" not in name else "Work in Progress")
    aux.Label(sim= False, status="" if "MC" not in name else "Work in Progress")
    aux.save(name, normal=False, changeMinMax=False)

    #if name == "final_lowEMHT": dc = limitTools.MyDatacard()
    #if name == "final_MC_MM": dc = limitTools.MyDatacard()
    if name == "final_MC": dc = limitTools.MyDatacard()
    #elif name == "final_highEMHT": dc = limitTools.MyDatacard("testDatacard.txt")
    else: return
    #for bin in range(dirHist.GetNbinsX()-2, dirHist.GetNbinsX()+1):
    #for bin in range(dirHist.GetNbinsX()-1, dirHist.GetNbinsX()+1):
    #for bin in range(dirHist.GetNbinsX()-2, dirHist.GetNbinsX()+1):
    for bin in range(dirHist.GetNbinsX()-1, dirHist.GetNbinsX()+1):
    #for bin in range(dirHist.GetNbinsX()-3, dirHist.GetNbinsX()+1):
        binName = "bin{}_{}".format(name.split("_")[1],bin)
        #print dirHist.GetBinWidth(bin)
        bw = dirHist.GetBinWidth(bin) if style.divideByBinWidth else 1.
        print (signal1.GetBinLowEdge(bin))
        #print (signal1.GetBinWidth(bin))
        print (signal1.GetBinContent(bin)*bw)
        #print "content",(dirHist.GetBinContent(bin))
        #print "stat. unc",(getDatacardUncertFromHist(dirHist,bin)-1.)*100.,"%"
        #print (signal1.GetBinContent(bin-1)*bw)
        #print (signal1.GetBinContent(bin+1)*bw)
        dc.addBin(binName, int(round(dirHist.GetBinContent(bin)*bw)),
            {
                #"signal": (signal1.GetBinContent(bin)-totStat.GetBinContent(bin))*bw,
                "signal": (signal1.GetBinContent(bin))*bw,
                #"signal": (signal1.GetBinContent(bin)-totStat.GetBinContent(bin)-signal1_pre.GetBinContent(bin))*bw,
                "ttg": ttgHist.GetBinContent(bin)*bw,
                "zg": zgHist.GetBinContent(bin)*bw,
                "zz": zzHist.GetBinContent(bin)*bw,
                "wwg": wwgHist.GetBinContent(bin)*bw,
                "wzg": wzgHist.GetBinContent(bin)*bw,
                "dy": dyHist.GetBinContent(bin)*bw,
                "wjets": wjetsHist.GetBinContent(bin)*bw,
                "tt": ttHist.GetBinContent(bin)*bw,
                "singletop": singletopHist.GetBinContent(bin)*bw,
                "wz": wzHist.GetBinContent(bin)*bw,
                "ww": wwHist.GetBinContent(bin)*bw,
                "wg": wgHist.GetBinContent(bin)*bw,
                "zz4l": zz4lHist.GetBinContent(bin)*bw,
#                "cont": signal1_pre.GetBinContent(bin)*bw
            }, {
                "ttgStat_"+binName: {"ttg": getDatacardUncertFromHist(ttgHist,bin)},
                "zgStat_"+binName: {"zg": getDatacardUncertFromHist(zgHist,bin)},
                "zzStat_"+binName: {"zz": getDatacardUncertFromHist(zzHist,bin)},
                "wwgStat_"+binName: {"wwg": getDatacardUncertFromHist(wwgHist,bin)},
                "wzgStat_"+binName: {"wzg": getDatacardUncertFromHist(wzgHist,bin)},
                "dyStat_"+binName: {"dy": getDatacardUncertFromHist(dyHist,bin)},
                "wjetsStat_"+binName: {"wjets": getDatacardUncertFromHist(wjetsHist,bin)},
                "ttStat_"+binName: {"tt": getDatacardUncertFromHist(ttHist,bin)},
                "singletopStat_"+binName: {"singletop": getDatacardUncertFromHist(singletopHist,bin)},
                "wzStat_"+binName: {"wz": getDatacardUncertFromHist(wzHist,bin)},
                "wwStat_"+binName: {"ww": getDatacardUncertFromHist(wwHist,bin)},
                "wgStat_"+binName: {"wg": getDatacardUncertFromHist(wgHist,bin)},
                "zz4lStat_"+binName: {"zz4l": getDatacardUncertFromHist(zz4lHist,bin)},
                "signalStat_"+binName: {"signal": getDatacardUncertFromHist(signal1,bin)},
                #"gqcdSyst": {"gqcd": getDatacardUncertFromHist(gjetSyst,bin)},
                #"eleSyst": {"ele": getDatacardUncertFromHist(eSyst,bin)},
                #"pdf": {
                    #"wg": getDatacardUncertFromHist(wgPdfUnc,bin),
                    #"zg": getDatacardUncertFromHist(zgPdfUnc,bin),
                    #"tg": getDatacardUncertFromHist(tgPdfUnc,bin)},
                #"scale": {
                    #"wg": getDatacardUncertFromHist(wgScaleUnc,bin),
                    #"zg": getDatacardUncertFromHist(zgScaleUnc,bin),
                    #"tg": getDatacardUncertFromHist(tgScaleUnc,bin)},
                "lumi": {
                    "signal": 1.026,
                    "ttg": 1.026,
                    "zg": 1.026,
                    "zz": 1.026,
                    "wwg": 1.026,
                    "wzg": 1.026,
                    "dy": 1.026,
                    "wjets": 1.026,
                    "tt": 1.026,
                    "singletop": 1.026,
                    "wz": 1.026,
                    "ww": 1.026,
                    "wg": 1.026,
                    "zz4l": 1.026},
                #"pu": {
                    #"wg": getDatacardUncertFromHist(wgPuUnc,bin),
                    #"zg": getDatacardUncertFromHist(zgPuUnc,bin),
                    #"tg": getDatacardUncertFromHist(tgPuUnc,bin)},
                #"jes": {
                    #"wg": getDatacardUncertFromHist(wgJesUnc,bin),
                    #"zg": getDatacardUncertFromHist(zgJesUnc,bin),
                    #"tg": getDatacardUncertFromHist(tgJesUnc,bin)},
                "SFDY": {
                    "zg": (1.+sfDYErr),
                    "dy": (1.+sfDYErr)},
                "SFTT": {
                    #"tt": getDatacardUncertFromHist(),
                    #"ttg": getDatacardUncertFromHist()},
                    "tt": (1.+0.04),
                    "ttg": (1.+0.04)},
                "SFZZ": {
                    "zz": (1.+sfZZErr),
                    "zz4l": (1.+sfZZErr)},
                "SFWZ": {
                    "wz": (1.+sfWZErr)},
                "dataMC": {
                    "signal": 1.05,
                    "ttg": 1.2,
                    "zg": 1.2,
                    "zz": 1.2,
                    "wwg": 1.2,
                    "wzg": 1.2,
                    "dy": 1.2,
                    "wjets": 1.2,
                    "tt": 1.2,
                    "singletop": 1.2,
                    "wz": 1.2,
                    "ww": 1.2,
                    "wg": 1.2,
                    "zz4l": 1.2},
                #"trigger": {
                    #"signal": 1.004,
                    #"wg": 1.004,
                    #"zg": 1.004,
                    #"tg": 1.004},
                #"isr": {"signal": 1.001},
                #"genMet": {"signal": 1.001},
                # jes, jer splitting
            }
        )
    dc.write("testDatacard.txt")
    dc.write("limitCalculations/testDatacard.txt")
    
def main():
    #allMC = gjets+qcd+zg+wg+ttg+wjets+ttjets_nlo+znunu
    allMC = zgamma+ttgamma+zz+wwgamma+wzgamma+DYjetsNLO+wjets+tt+singletop+wz+ww+zz4l+wjets
    allMC.label = "MC mix"
    #finalDistributionSignalHist("final_lowEMHT", data, "signal_lowEMHT", dataHt, data, "signal_lowEMHT_eControl")
    #finalDistributionSignalHist("allMC_lowEMHT", allMc, "signal_lowEMHT", allMc, allMc, "signal_lowEMHT_eControl")
    #finalDistributionSignalHist("final_MC_MM", allMC, "onZMM")
    #finalDistributionSignalHist("final_MC_EE", allMC, "onZEE")
    #finalDistributionSignalHist("final_MC", allMC, "onZMet100")
    #finalDistributionSignalHist("final_MC", allMC, "onZMet150")
    finalDistributionSignalHist("final_MC", allMC, "onZ/LL")

main()
