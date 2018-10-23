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
    #nBinsData = [0,25,50,75,100,150]
    
    #nBins2d=[[0,25,75,100,150,200,500],[0.,2.,3.5]]
    nBins2d=[[0,25,75,100,150,200,500],[0.,100.,300]]
    #nBins2d=[[0,25,75,100,150,210,500],[0.,100.,300]]
    #[[150,250,500],[0,2.,3.5]]
    #nBins = [200, 300,500]
    #nBins = np.concatenate((np.arange(0,200,50),np.arange(200,300,100),np.arange(300,3000,1400)),axis=0)
    #nBins = np.concatenate((np.arange(0,200,25),np.arange(200,300,100),np.arange(300,600,200)),axis=0)
    #nBins = np.concatenate((np.arange(0,200,25),np.arange(200,500,300)),axis=0)
    #print nBins

    # direct stuff
    
    #dirHist = aux.stdHist2d(dirSet, dirDir+"/MetDeltaPhiLL", nBins2d)
    
    #dataHist = aux.stdHist(dataDoubleSF,dirDir+"/met",nBins)
    #dataHist = aux.stdHist(dataLL,dirDir+"/met",nBins)
    #dataHist = aux.rebin(dataHist,nBinsData)
    #dataHist.SetYTitle(aux.getYAxisTitle(dataHist))

    style.additionalPoissonUncertainty = False
    #aux.drawOpt(dirHist, "data")


    #aux.drawOpt(dataHist, "data")
    #gjetHist, gjetSyst, info = gjetPrediction(dirHist, preSet, zg+wg+ttg+wjets+ttjets_nlo+znunu, "met", nBins, weight, name+"_divByBinWidth" if style.divideByBinWidth else name)
    #gjetHist.SetLineColor(rwth.myLightBlue)
    #gjetHist.GetXaxis().SetTitle("#it{p}_{T}^{miss} (GeV)")


    #eHist = aux.stdHist(preSetElectron, preDirElectron+"/met", nBins)
    #eHist.GetXaxis().SetTitle("#it{p}_{T}^{miss} (GeV)")
    #eHist.Scale( 0.0267 if dirSet == data else 0.0154 )
    #eHist.SetLineColor(rwth.myYellow)
    #eSyst = aux.getSysHisto(eHist, 0.3)

#zgamma+ttgamma+zz+wwgamma+wzgamma+DYjetsNLO+wjets+tt+singletop+wz+ww+zz4l+wjets

    #zgHist = aux.stdHist2d(zgamma, dirDir+"/MetDeltaPhiLL", nBins2d)
    #ttgHist = aux.stdHist2d(ttgamma, dirDir+"/MetDeltaPhiLL", nBins2d)
    #zzHist = aux.stdHist2d(zz, dirDir+"/MetDeltaPhiLL", nBins2d)
    #wwgHist = aux.stdHist2d(wwgamma, dirDir+"/MetDeltaPhiLL", nBins2d)
    #wzgHist = aux.stdHist2d(wzgamma, dirDir+"/MetDeltaPhiLL", nBins2d)
    #dyHist = aux.stdHist2d(DYjetsNLO, dirDir+"/MetDeltaPhiLL", nBins2d)
    #wjetsHist = aux.stdHist2d(wjets, dirDir+"/MetDeltaPhiLL", nBins2d)
    #ttHist = aux.stdHist2d(tt, dirDir+"/MetDeltaPhiLL", nBins2d)
    #singletopHist=aux.stdHist2d(singletop, dirDir+"/MetDeltaPhiLL", nBins2d)
    #wzHist=aux.stdHist2d(wz, dirDir+"/MetDeltaPhiLL", nBins2d)
    #wwHist=aux.stdHist2d(ww, dirDir+"/MetDeltaPhiLL", nBins2d)
    #zz4lHist=aux.stdHist2d(zz4l, dirDir+"/MetDeltaPhiLL", nBins2d)
    #wgHist=aux.stdHist2d(wgamma,dirDir+"/MetDeltaPhiLL",nBins2d)
    zgHist = aux.stdHist2d(zgamma, dirDir+"/MetMt2", nBins2d)
    ttgHist = aux.stdHist2d(ttgamma, dirDir+"/MetMt2", nBins2d)
    zzHist = aux.stdHist2d(zz, dirDir+"/MetMt2", nBins2d)
    wwgHist = aux.stdHist2d(wwgamma, dirDir+"/MetMt2", nBins2d)
    wzgHist = aux.stdHist2d(wzgamma, dirDir+"/MetMt2", nBins2d)
    dyHist = aux.stdHist2d(DYjetsNLO, dirDir+"/MetMt2", nBins2d)
    wjetsHist = aux.stdHist2d(wjets, dirDir+"/MetMt2", nBins2d)
    ttHist = aux.stdHist2d(tt, dirDir+"/MetMt2", nBins2d)
    singletopHist=aux.stdHist2d(singletop, dirDir+"/MetMt2", nBins2d)
    wzHist=aux.stdHist2d(wz, dirDir+"/MetMt2", nBins2d)
    wwHist=aux.stdHist2d(ww, dirDir+"/MetMt2", nBins2d)
    zz4lHist=aux.stdHist2d(zz4l, dirDir+"/MetMt2", nBins2d)
    wgHist=aux.stdHist2d(wgamma,dirDir+"/MetMt2",nBins2d)
    
    
    #zgHist.SaveAs("plots_signal2d/histos/zg.root")
    #ttgHist.SaveAs("plots_signal2d/histos/ttg.root")
    #zzHist.SaveAs("plots_signal2d/histos/zz.root")
    #wwgHist.SaveAs("plots_signal2d/histos/wwg.root")
    #wzgHist.SaveAs("plots_signal2d/histos/wzg.root")
    #dyHist.SaveAs("plots_signal2d/histos/dy.root")
    #wjetsHist.SaveAs("plots_signal2d/histos/wjets.root")
    #ttHist.SaveAs("plots_signal2d/histos/tt.root")
    #singletopHist.SaveAs("plots_signal2d/histos/singletop.root")
    #wzHist.SaveAs("plots_signal2d/histos/wz.root")
    #wwHist.SaveAs("plots_signal2d/histos/ww.root")
    #zz4lHist.SaveAs("plots_signal2d/histos/zz4l.root")
    #wgHist.SaveAs("plots_signal2d/histos/wg.root")

    #for bin in range(ttgHist.GetNbinsX()-1, ttgHist.GetNbinsX()+1):
    #for bin in range(0, ttgHist.GetNbinsX()+1):
        #print bin,ttgHist.GetXaxis().GetBinLowEdge(bin),ttgHist.
    
    #Scaling
    #zg_AvgTopPtWeightHisto = zgamma.getHist(dirDir+"/weight_topPt")
    #zg_AvgNIsrWeightHisto = zgamma.getHist(dirDir+"/weight_nISR")
    #zg_AvgEWKinoWeightHisto = zgamma.getHist(dirDir+"/weight_EWKinoPairPt")
    #ttg_AvgTopPtWeightHisto = ttgamma.getHist(dirDir+"/weight_topPt")
    #ttg_AvgNIsrWeightHisto = ttgamma.getHist(dirDir+"/weight_nISR")
    #ttg_AvgEWKinoWeightHisto = ttgamma.getHist(dirDir+"/weight_EWKinoPairPt")
    #zz_AvgTopPtWeightHisto = zz.getHist(dirDir+"/weight_topPt")
    #zz_AvgNIsrWeightHisto = zz.getHist(dirDir+"/weight_nISR")
    #zz_AvgEWKinoWeightHisto = zz.getHist(dirDir+"/weight_EWKinoPairPt")
    #wwg_AvgTopPtWeightHisto = wwgamma.getHist(dirDir+"/weight_topPt")
    #wwg_AvgNIsrWeightHisto = wwgamma.getHist(dirDir+"/weight_nISR")
    #wwg_AvgEWKinoWeightHisto = wwgamma.getHist(dirDir+"/weight_EWKinoPairPt")
    #wzg_AvgTopPtWeightHisto = wzgamma.getHist(dirDir+"/weight_topPt")
    #wzg_AvgNIsrWeightHisto = wzgamma.getHist(dirDir+"/weight_nISR")
    #wzg_AvgEWKinoWeightHisto = wzgamma.getHist(dirDir+"/weight_EWKinoPairPt")
    #dy_AvgTopPtWeightHisto = DYjetsNLO.getHist(dirDir+"/weight_topPt")
    #dy_AvgNIsrWeightHisto = DYjetsNLO.getHist(dirDir+"/weight_nISR")
    #dy_AvgEWKinoWeightHisto = DYjetsNLO.getHist(dirDir+"/weight_EWKinoPairPt")
    #wjets_AvgTopPtWeightHisto = wjets.getHist(dirDir+"/weight_topPt")
    #wjets_AvgNIsrWeightHisto = wjets.getHist(dirDir+"/weight_nISR")
    #wjets_AvgEWKinoWeightHisto = wjets.getHist(dirDir+"/weight_EWKinoPairPt")
    #tt_AvgTopPtWeightHisto = tt.getHist(dirDir+"/weight_topPt")
    #tt_AvgNIsrWeightHisto = tt.getHist(dirDir+"/weight_nISR")
    #tt_AvgEWKinoWeightHisto = tt.getHist(dirDir+"/weight_EWKinoPairPt")
    #singletop_AvgTopPtWeightHisto = singletop.getHist(dirDir+"/weight_topPt")
    #singletop_AvgNIsrWeightHisto = singletop.getHist(dirDir+"/weight_nISR")
    #singletop_AvgEWKinoWeightHisto = singletop.getHist(dirDir+"/weight_EWKinoPairPt")
    #wz_AvgTopPtWeightHisto = wz.getHist(dirDir+"/weight_topPt")
    #wz_AvgNIsrWeightHisto = wz.getHist(dirDir+"/weight_nISR")
    #wz_AvgEWKinoWeightHisto = wz.getHist(dirDir+"/weight_EWKinoPairPt")
    #ww_AvgTopPtWeightHisto = ww.getHist(dirDir+"/weight_topPt")
    #ww_AvgNIsrWeightHisto = ww.getHist(dirDir+"/weight_nISR")
    #ww_AvgEWKinoWeightHisto = ww.getHist(dirDir+"/weight_EWKinoPairPt")
    #zz4l_AvgTopPtWeightHisto = zz4l.getHist(dirDir+"/weight_topPt")
    #zz4l_AvgNIsrWeightHisto = zz4l.getHist(dirDir+"/weight_nISR")
    #zz4l_AvgEWKinoWeightHisto = zz4l.getHist(dirDir+"/weight_EWKinoPairPt")
    #wg_AvgTopPtWeightHisto = wgamma.getHist(dirDir+"/weight_topPt")
    #wg_AvgNIsrWeightHisto = wgamma.getHist(dirDir+"/weight_nISR")
    #wg_AvgEWKinoWeightHisto = wgamma.getHist(dirDir+"/weight_EWKinoPairPt")
    
    histsToScale=[zgHist,ttgHist,zzHist,wwgHist,wzgHist,dyHist,wjetsHist,ttHist,singletopHist,wzHist,wwHist,zz4lHist,wgHist]
    #topWeightHists=[zg_AvgTopPtWeightHisto,ttg_AvgTopPtWeightHisto,zz_AvgTopPtWeightHisto,wwg_AvgTopPtWeightHisto,wzg_AvgTopPtWeightHisto,dy_AvgTopPtWeightHisto,wjets_AvgTopPtWeightHisto,tt_AvgTopPtWeightHisto,singletop_AvgTopPtWeightHisto,wz_AvgTopPtWeightHisto,ww_AvgTopPtWeightHisto,zz4l_AvgTopPtWeightHisto,wg_AvgTopPtWeightHisto]
    #nISRWeightHists=[zg_AvgNIsrWeightHisto,ttg_AvgNIsrWeightHisto,zz_AvgNIsrWeightHisto,wwg_AvgNIsrWeightHisto,wzg_AvgNIsrWeightHisto,dy_AvgNIsrWeightHisto,wjets_AvgNIsrWeightHisto,tt_AvgNIsrWeightHisto,singletop_AvgNIsrWeightHisto,wz_AvgNIsrWeightHisto,ww_AvgNIsrWeightHisto,zz4l_AvgNIsrWeightHisto,wg_AvgNIsrWeightHisto]
    #EWKinoWeightHists=[zg_AvgEWKinoWeightHisto,ttg_AvgEWKinoWeightHisto,zz_AvgEWKinoWeightHisto,wwg_AvgEWKinoWeightHisto,wzg_AvgEWKinoWeightHisto,dy_AvgEWKinoWeightHisto,wjets_AvgEWKinoWeightHisto,tt_AvgEWKinoWeightHisto,singletop_AvgEWKinoWeightHisto,wz_AvgEWKinoWeightHisto,ww_AvgEWKinoWeightHisto,zz4l_AvgEWKinoWeightHisto,wg_AvgEWKinoWeightHisto]
    #histsToScale=[ttgHist]
    #topWeightHists=[ttg_AvgTopPtWeightHisto]
    #nISRWeightHists=[ttg_AvgNIsrWeightHisto]
    #EWKinoWeightHists=[ttg_AvgEWKinoWeightHisto]
      
    #print len(histsToScale),len(topWeightHists),len(nISRWeightHists),len(EWKinoWeightHists),len(leptonWeightHists)
    
    #for i in range(len(histsToScale)):
        #if topWeightHists[i].Integral()>0.:
            #histsToScale[i].Scale(1./topWeightHists[i].GetMean())
        #if nISRWeightHists[i].Integral()>0.:
            #histsToScale[i].Scale(1./nISRWeightHists[i].GetMean())
        #if EWKinoWeightHists[i].Integral()>0.:
            #histsToScale[i].Scale(1./EWKinoWeightHists[i].GetMean())
            
    
    
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
    
    #print "INT", ttgHist.Integral(),ttHist.Integral()
    
    #dirHist = aux.addHists(*histsToScale)


    #for bin in range(dirHist.GetNbinsX()-1, dirHist.GetNbinsX()+1):
        #print dirHist.getXaxis().GetBinLowEdge(bin)


    #dirHistLowMet=aux.getProjection(zgHist,"y",5,5)
    #dirHistHighMet=aux.getProjection(zgHist,"y",6,-1)


    #aux.drawOpt(dirHistLowMet, "data")
    #aux.drawOpt(dirHistHighMet, "data")


    #nBins2d=[[0,25,75,100,150,200,500],[0.,100.,300]]

    #zgHistLowMet=aux.getProjection(zgHist,"y",5,5)
    #ttgHistLowMet=aux.getProjection(ttgHist,"y",5,5)
    #zzHistLowMet=aux.getProjection(zzHist,"y",5,5)
    #wwgHistLowMet=aux.getProjection(wwgHist,"y",5,5)
    #wzgHistLowMet=aux.getProjection(wzgHist,"y",5,5)
    #dyHistLowMet=aux.getProjection(dyHist,"y",5,5)
    #wjetsHistLowMet=aux.getProjection(wjetsHist,"y",5,5)
    #ttHistLowMet=aux.getProjection(ttHist,"y",5,5)
    #singletopHistLowMet=aux.getProjection(singletopHist,"y",5,5)
    #wzHistLowMet=aux.getProjection(wzHist,"y",5,5)
    #wwHistLowMet=aux.getProjection(wwHist,"y",5,5)
    #zz4lHistLowMet=aux.getProjection(zz4lHist,"y",5,5)
    #wgHistLowMet=aux.getProjection(wgHist,"y",5,5)
    
    #zgHistHighMet=aux.getProjection(zgHist,"y",6,-1)
    #ttgHistHighMet=aux.getProjection(ttgHist,"y",6,-1)
    #zzHistHighMet=aux.getProjection(zzHist,"y",6,-1)
    #wwgHistHighMet=aux.getProjection(wwgHist,"y",6,-1)
    #wzgHistHighMet=aux.getProjection(wzgHist,"y",6,-1)
    #dyHistHighMet=aux.getProjection(dyHist,"y",6,-1)
    #wjetsHistHighMet=aux.getProjection(wjetsHist,"y",6,-1)
    #ttHistHighMet=aux.getProjection(ttHist,"y",6,-1)
    #singletopHistHighMet=aux.getProjection(singletopHist,"y",6,-1)
    #wzHistHighMet=aux.getProjection(wzHist,"y",6,-1)
    #wwHistHighMet=aux.getProjection(wwHist,"y",6,-1)
    #zz4lHistHighMet=aux.getProjection(zz4lHist,"y",6,-1)
    #wgHistHighMet=aux.getProjection(wgHist,"y",6,-1)
    #zgHistLowMet=aux.getProjection(zgHist,"x",1,1)
    #ttgHistLowMet=aux.getProjection(ttgHist,"x",1,1)
    #zzHistLowMet=aux.getProjection(zzHist,"x",1,1)
    #wwgHistLowMet=aux.getProjection(wwgHist,"x",1,1)
    #wzgHistLowMet=aux.getProjection(wzgHist,"x",1,1)
    #dyHistLowMet=aux.getProjection(dyHist,"x",1,1)
    #wjetsHistLowMet=aux.getProjection(wjetsHist,"x",1,1)
    #ttHistLowMet=aux.getProjection(ttHist,"x",1,1)
    #singletopHistLowMet=aux.getProjection(singletopHist,"x",1,1)
    #wzHistLowMet=aux.getProjection(wzHist,"x",1,1)
    #wwHistLowMet=aux.getProjection(wwHist,"x",1,1)
    #zz4lHistLowMet=aux.getProjection(zz4lHist,"x",1,1)
    #wgHistLowMet=aux.getProjection(wgHist,"x",1,1)
    
    zgHistHighMt2=aux.getProjection(zgHist,"x",2,-1)
    ttgHistHighMt2=aux.getProjection(ttgHist,"x",2,-1)
    zzHistHighMt2=aux.getProjection(zzHist,"x",2,-1)
    wwgHistHighMt2=aux.getProjection(wwgHist,"x",2,-1)
    wzgHistHighMt2=aux.getProjection(wzgHist,"x",2,-1)
    dyHistHighMt2=aux.getProjection(dyHist,"x",2,-1)
    wjetsHistHighMt2=aux.getProjection(wjetsHist,"x",2,-1)
    ttHistHighMt2=aux.getProjection(ttHist,"x",2,-1)
    singletopHistHighMt2=aux.getProjection(singletopHist,"x",2,-1)
    wzHistHighMt2=aux.getProjection(wzHist,"x",2,-1)
    wwHistHighMt2=aux.getProjection(wwHist,"x",2,-1)
    zz4lHistHighMt2=aux.getProjection(zz4lHist,"x",2,-1)
    wgHistHighMt2=aux.getProjection(wgHist,"x",2,-1)


    #zgHistLowMet.SetLineColor(ROOT.kGreen-3)
    #ttgHistLowMet.SetLineColor(ROOT.kRed+1)
    #zzHistLowMet.SetLineColor(ROOT.kYellow)
    #wwgHistLowMet.SetLineColor(ROOT.kCyan+3)
    #wzgHistLowMet.SetLineColor(ROOT.kBlue+3)
    #dyHistLowMet.SetLineColor(ROOT.kGreen+3)
    #wjetsHistLowMet.SetLineColor(ROOT.kBlue-9)
    #ttHistLowMet.SetLineColor(ROOT.kOrange+8)
    #singletopHistLowMet.SetLineColor(ROOT.kOrange+4)
    #wzHistLowMet.SetLineColor(ROOT.kAzure-5)
    #wwHistLowMet.SetLineColor(ROOT.kCyan-3)
    #zz4lHistLowMet.SetLineColor(ROOT.kOrange-2)
    #wgHistLowMet.SetLineColor(ROOT.kRed+3)
    
    #totStatLowMet = aux.addHists(zgHistLowMet, ttgHistLowMet, zzHistLowMet, wwgHistLowMet, wzgHistLowMet,dyHistLowMet,wjetsHistLowMet,ttHistLowMet,singletopHistLowMet,wzHistLowMet,wwHistLowMet,zz4lHistLowMet,wjetsHistLowMet,wgHistLowMet)
    
    zgHistHighMt2.SetLineColor(ROOT.kGreen-3)
    ttgHistHighMt2.SetLineColor(ROOT.kRed+1)
    zzHistHighMt2.SetLineColor(ROOT.kYellow)
    wwgHistHighMt2.SetLineColor(ROOT.kCyan+3)
    wzgHistHighMt2.SetLineColor(ROOT.kBlue+3)
    dyHistHighMt2.SetLineColor(ROOT.kGreen+3)
    wjetsHistHighMt2.SetLineColor(ROOT.kBlue-9)
    ttHistHighMt2.SetLineColor(ROOT.kOrange+8)
    singletopHistHighMt2.SetLineColor(ROOT.kOrange+4)
    wzHistHighMt2.SetLineColor(ROOT.kAzure-5)
    wwHistHighMt2.SetLineColor(ROOT.kCyan-3)
    zz4lHistHighMt2.SetLineColor(ROOT.kOrange-2)
    wgHistHighMt2.SetLineColor(ROOT.kRed+3)

    totStatHighMt2 = aux.addHists(zgHistHighMt2, ttgHistHighMt2, zzHistHighMt2, wwgHistHighMt2, wzgHistHighMt2,dyHistHighMt2,wjetsHistHighMt2,ttHistHighMt2,singletopHistHighMt2,wzHistHighMt2,wwHistHighMt2,zz4lHistHighMt2,wjetsHistHighMt2,wgHistHighMt2)




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
    
    #zgSyst = aux.getSysHisto(zgHist, mcSystUncert)
    #ttgSyst = aux.getSysHisto(ttgHist, mcSystUncert)
    #zzSyst = aux.getSysHisto(zzHist, mcSystUncert)
    #wwgSyst = aux.getSysHisto(wwgHist, mcSystUncert)
    #wzgSyst = aux.getSysHisto(wzgHist, mcSystUncert)
    #dySyst = aux.getSysHisto(dyHist, mcSystUncert)
    #wjetsSyst = aux.getSysHisto(wjetsHist, mcSystUncert)
    #ttSyst = aux.getSysHisto(ttHist, mcSystUncert)
    #singletopSyst=aux.getSysHisto(singletopHist,mcSystUncert)
    #wzSyst=aux.getSysHisto(wzHist,mcSystUncert)
    #wwSyst=aux.getSysHisto(wwHist,mcSystUncert)
    #zz4lSyst=aux.getSysHisto(zz4lHist,mcSystUncert)
    #wgSyst=aux.getSysHisto(wgHist,mcSystUncert)

    #dirSyst= aux.getSysHisto(aux.addHists(zgSyst,ttgSyst,zzSyst,wwgSyst,wzgSyst,dySyst,wjetsSyst,ttSyst,singletopSyst,wzSyst,wwSyst,zz4lSyst,wgSyst),mcSystUncert)

    #zgSyst = aux.addUncertaintiesQuadratic([zgSyst,zgPdfUnc,zgScaleUnc,zgJesUnc,zgPuUnc])
    #wgSyst = aux.addUncertaintiesQuadratic([wgSyst,wgPdfUnc,wgScaleUnc,wgJesUnc,wgPuUnc])
    #tgSyst = aux.addUncertaintiesQuadratic([tgSyst,tgPdfUnc,tgScaleUnc,tgJesUnc,tgPuUnc])


    #totStat = aux.addHists(zgHist, ttgHist, zzHist, wwgHist, wzgHist,dyHist,wjetsHist,ttHist,singletopHist,wzHist,wwHist,zz4lHist,wjetsHist,wgHist)
    #totSyst = aux.addHists(zgSyst, ttgSyst, zzSyst, wwgSyst, wzgSyst,dySyst,wjetsSyst,ttSyst,singletopSyst,wzSyst,wwSyst,zz4lSyst,wjetsSyst,wgSyst)


    signal2 = aux.stdHist2d(gmsb_290_205, dirDir+"/MetMt2", nBins2d)
    signal3 = aux.stdHist2d(t5bbbbzg_1500_1400, dirDir+"/MetMt2", nBins2d)
    signal1 = aux.stdHist2d(tching_600, dirDir+"/MetMt2", nBins2d)
    
    #signal3_AvgTopPtWeightHisto = t5bbbbzg_1500_1400.getHist(dirDir+"/weight_topPt")
    #signal3_AvgNIsrWeightHisto = t5bbbbzg_1500_1400.getHist(dirDir+"/weight_nISR")
    #signal3_AvgEWKinoWeightHisto = t5bbbbzg_1500_1400.getHist(dirDir+"/weight_EWKinoPairPt")
    #signal2_AvgleptonWeightHisto = t5bbbbzg_1500_400.getHist(dirDir+"/weight_leptonPairPt")
    #signal1_AvgTopPtWeightHisto = tching_600.getHist(dirDir+"/weight_topPt")
    #signal1_AvgNIsrWeightHisto = tching_600.getHist(dirDir+"/weight_nISR")
    #signal1_AvgEWKinoWeightHisto = tching_600.getHist(dirDir+"/weight_EWKinoPairPt")
    #signal1_AvgleptonWeightHisto = tching_600.getHist(dirDir+"/weight_leptonPairPt")
    #signal2_AvgTopPtWeightHisto = gmsb_290_205.getHist(dirDir+"/weight_topPt")
    #signal2_AvgNIsrWeightHisto = gmsb_290_205.getHist(dirDir+"/weight_nISR")
    #signal2_AvgEWKinoWeightHisto = gmsb_290_205.getHist(dirDir+"/weight_EWKinoPairPt")
    #signal2_AvgleptonWeightHisto = gmsb_290_205.getHist(dirDir+"/weight_leptonPairPt")

#
    #signal2.Scale((1./signal2_AvgTopPtWeightHisto.GetMean()))
    #signal2.Scale((1./signal2_AvgNIsrWeightHisto.GetMean()))
    #signal2.Scale((1./signal2_AvgEWKinoWeightHisto.GetMean()))
    #signal2.Scale((1./signal2_AvgleptonWeightHisto.GetMean()))
    #signal1.Scale((1./signal1_AvgTopPtWeightHisto.GetMean()))
    #signal1.Scale((1./signal1_AvgNIsrWeightHisto.GetMean()))
    #signal1.Scale((1./signal1_AvgEWKinoWeightHisto.GetMean()))
    #signal1.Scale((1./signal1_AvgleptonWeightHisto.GetMean()))
    #signal3.Scale((1./signal3_AvgTopPtWeightHisto.GetMean()))
    #signal3.Scale((1./signal3_AvgNIsrWeightHisto.GetMean()))
    #signal3.Scale((1./signal3_AvgEWKinoWeightHisto.GetMean()))
    #signal3.Scale((1./signal3_AvgleptonWeightHisto.GetMean()))

    #signal1LowMet=aux.getProjection(signal1,"y",5,5)
    #signal1HighMet=aux.getProjection(signal1,"y",6,-1)
    #signal2LowMet=aux.getProjection(signal2,"y",5,5)
    #signal2HighMet=aux.getProjection(signal2,"y",6,-1)
    #signal3LowMet=aux.getProjection(signal3,"y",5,5)
    #signal3HighMet=aux.getProjection(signal3,"y",6,-1)
    #signal1LowMet=aux.getProjection(signal1,"x",1,1)
    signal1HighMt2=aux.getProjection(signal1,"x",2,-1)
    #signal2LowMet=aux.getProjection(signal2,"x",1,1)
    signal2HighMt2=aux.getProjection(signal2,"x",2,-1)
    #signal3LowMet=aux.getProjection(signal3,"x",1,1)
    signal3HighMt2=aux.getProjection(signal3,"x",2,-1)

    
    #print signal1.GetEntries()
    
    #for h in signal1LowMet, signal2LowMet,signal2HighMet,signal1HighMet,signal3HighMet,signal3LowMet:
    for h in signal1HighMt2, signal2HighMt2,signal3HighMt2:
        aux.drawOpt(h, "signal")
    #    ###h.Add(totStat)
    #signal1LowMet.SetLineColor(ROOT.kBlue+3)
    signal1HighMt2.SetLineColor(ROOT.kBlue+3)
    #signal3LowMet.SetLineColor(ROOT.kRed-3)
    signal3HighMt2.SetLineColor(ROOT.kRed-3)
    #signal2LowMet.SetLineColor(ROOT.kBlue+3)
    #signal2LowMet.SetLineStyle(2)
    #signal3LowMet.SetLineStyle(2)
    signal3HighMt2.SetLineStyle(2)
    signal2HighMt2.SetLineColor(ROOT.kBlue+3)
    signal2HighMt2.SetLineStyle(2)

    #signal1_pre = aux.createHistoFromDatasetTree(t5wg_1600_100, "met*{}".format(info["shift"]), weight, nBins, "tr_jControl/simpleTree")
    #signal1_pre.Scale(info["scale"])

    #totSyst = gjetSyst

    #totUnc = aux.addHistUncert(totStat, totSyst)
    #totUnc = aux.addHistUncert(totStat, dirSyst)
    #totUnc = aux.addHistUncert(dirSyst)
    #aux.drawOpt(totUnc, "totUnc")
    #totUnc = dirSyst
    #aux.drawOpt(totUnc, "totUnc")
    
    c = ROOT.TCanvas()
    #c2 = ROOT.TCanvas()
    #m = multiplot.Multiplot()
    m2 = multiplot.Multiplot()
    #if dirSet == data:
    
    #if dirSet == dataDoubleSF:
        #m.add(dirHist, "Data")
    #else:
        #m.add(dirHist, "Direct simulation")
        
        
#   m.add(signal1_pre, "contamination")
    #m.add(signal1, "T5bbbbZg")
    #m.add(signal2, "TChiNg")
    #m.add(signal2, "T5bbbbZg")
    #m.add(signal2LowMet, "GMSB")
    #m.add(signal1LowMet, "TChiNg")
    #m.add(signal3LowMet, "t5bbbbzg_1500_1400")
    #m.addStack(zgHistLowMet, "Z#gamma")
    #m.addStack(ttgHistLowMet, "t#bar{t}#gamma")
    #m.addStack(zzHistLowMet, "ZZ(#rightarrowll#nu#nu)")
    #m.addStack(wwgHistLowMet, "WW#gamma")
    #m.addStack(wzgHistLowMet, "WZ#gamma")
    #m.addStack(dyHistLowMet, "Drell-Yan/Z")
    #m.addStack(wjetsHistLowMet, "W+jets")
    #m.addStack(ttHistLowMet, "t#bar{t}")
    #m.addStack(singletopHistLowMet, "single t")
    #m.addStack(wwHistLowMet, "WW")
    #m.addStack(wzHistLowMet, "WZ")
    #m.addStack(zz4lHistLowMet, "ZZ(#rightarrow4l)")
    #m.addStack(wgHistLowMet, "W#gamma")
    
    #aux.drawOpt(totStatLowMet, "totUnc")
    
    #m.add(totStatLowMet, "statistical uncertainty")
    
    m2.add(signal2HighMt2, "GMSB")
    m2.add(signal1HighMt2, "TChiNg")
    m2.add(signal3HighMt2, "t5bbbbzg_1500_1400")
    m2.addStack(zgHistHighMt2, "Z#gamma")
    m2.addStack(ttgHistHighMt2, "t#bar{t}#gamma")
    m2.addStack(zzHistHighMt2, "ZZ(#rightarrowll#nu#nu)")
    m2.addStack(wwgHistHighMt2, "WW#gamma")
    m2.addStack(wzgHistHighMt2, "WZ#gamma")
    m2.addStack(dyHistHighMt2, "Drell-Yan/Z")
    m2.addStack(wjetsHistHighMt2, "W+jets")
    m2.addStack(ttHistHighMt2, "t#bar{t}")
    m2.addStack(singletopHistHighMt2, "single t")
    m2.addStack(wwHistHighMt2, "WW")
    m2.addStack(wzHistHighMt2, "WZ")
    m2.addStack(zz4lHistHighMt2, "ZZ(#rightarrow4l)")
    m2.addStack(wgHistHighMt2, "W#gamma")
    
    aux.drawOpt(totStatHighMt2, "totUnc")
    
    m2.add(totStatHighMt2, "statistical uncertainty")    
    
    #m.histsToStack = sorted( m.histsToStack, key=lambda x: x.Integral(2,-1) )
    #m.histsToStack = sorted( m.histsToStack, key=lambda x: x.Integral(1,-1) )
    #m.histsToStack = sorted( m.histsToStack, key=lambda x: x.Integral(5,-1) )
    #m.histsToStack = sorted( m.histsToStack, key=lambda x: x.Integral(0,-1) )
    m2.histsToStack = sorted( m2.histsToStack, key=lambda x: x.Integral(0,-1) )


    #m.add(dataHist,"Data")

    #m.add(totUnc, "Total uncertainty")
    #m.maximum = 2.6*m.getMaximum()
    #m.maximum = 5.*m.getMaximum()
    #m.maximum = 25.
    #m.minimum = m.getMinimum()
    #m2.maximum = 25.
    #m2.maximum = 2.6*m2.getMaximum()
    m2.maximum = 5.*m2.getMaximum()
    m2.minimum = m2.getMinimum()
    #if "final_lowEMHT" in name: m.minimum = 4e-2
    #if "final_highEMHT" in name: m.minimum = 2e-3
    #legInfo = "#it{H}_{T}^{#gamma} < 2TeV" if "lowEMHT" in name else "2TeV < #it{H}_{T}^{#gamma}"
    #if "ee" in name: legInfo += ", EE"
    #legInfo += ", |#Delta#phi|>0.3"
    #legInfo = "DiMu" if "MM" in name else "DiEle"
    legInfo = "ee+#mu#mu"
    #m.leg.SetHeader(legInfo)
    m2.leg.SetHeader(legInfo)
    #m.leg.SetY1(.56)
    #m.leg.SetX1(.56)
    #m.leg.SetX1(.46)
    #m.leg.SetX2(.99)
    #m.leg.SetX2(.89)


    #m.Draw()

    # draw other labels

    #l = ROOT.TLine()
    #l.SetLineStyle(2)
    #l.SetLineColor(ROOT.kGray+2)
    #l.SetLineColor(ROOT.kRed-2)
    #l.SetLineColor(ROOT.kBlack)
    #text = ROOT.TLatex()
    #text.SetTextSize(0.8*text.GetTextSize())
    #l.DrawLine(100, 0, 100, totUnc.GetBinContent(totUnc.FindBin(100)))
    #text.SetTextAngle(90)
    #text.SetTextAngle(0)
    #text.DrawLatexNDC(.23,.315, "Normalization")
    #text.DrawLatexNDC(.22,.2, "CR")
    #text.DrawLatexNDC(.22,.4, "CR")
    #text.SetTextAngle(0)
    #text.DrawLatexNDC(.311,.315, "Validation")
    #text.DrawLatexNDC(.361,.2, "VR")
    #text.DrawLatexNDC(.681,.2, "SR")
    #text.DrawLatexNDC(.361,.4, "VR")
    #text.DrawLatexNDC(.681,.4, "SR")
    #if "final" in name:
        #l.DrawLine(150, 0, 150, totUnc.GetBinContent(totUnc.FindBin(150)))
        #l.DrawLine(100, 0, 100, totUnc.GetBinContent(totUnc.FindBin(100)))



    #r = ratio.Ratio("#scale[.9]{#lower[.24]{#splitline{Data/Pred.}{Bkg. frac.}}}", dirHist, totStat)
    #r = ratio.Ratio("#scale[.9]{#lower[.24]{Bkg. frac.}}", dirHist, totStat)
    #hsm = m.hists[0].GetStack().Last()
    #r = ratio.Ratio("#scale[.9]{#lower[.24]{Bkg. frac.}}", dataHist, hsm,sysHisto=totSyst)
    #rMax = 1.5
    #rMax = 1.1
    #if name == "final_lowEMHT": rMax = 1.6
    #if name == "final_highEMHT": rMax = 3.6
    #r.draw(0., rMax, m.getStack(), True)
    #r.draw(0., rMax, m.getStack())

    #aux.Label(sim= not dirSet==data, status="" if "allMC" not in name else "Private Work")
    #aux.Label(sim= not dirSet==dataDoubleSF, status="" if "allMC" not in name else "Private Work")
    #aux.Label(sim= not dirSet==dataDoubleSF, status="" if "MC" not in name else "Private Work")
    #aux.Label(sim= not dirSet==dataDoubleSF, status="" if "MC" not in name else "Work in Progress")
    #aux.Label(sim= False, status="" if "MC" not in name else "Work in Progress")
    #aux.save(name+"_lowMet", normal=False, changeMinMax=False,folder="plots_signal2d/")
    c = ROOT.TCanvas()
    m2.Draw()
    aux.save(name+"_highMt2", normal=False, changeMinMax=False,folder="plots_signal2d/")


    #if name == "final_MC": dc = limitTools.MyDatacard()
    #else: return
    #for bin in range(dirHist.GetNbinsX()-1, dirHist.GetNbinsX()+1):
        #binName = "bin{}_{}".format(name.split("_")[1],bin)
        #bw = dirHist.GetBinWidth(bin) if style.divideByBinWidth else 1.
        #print (signal1.GetBinLowEdge(bin))
        #print (signal1.GetBinContent(bin)*bw)
        #dc.addBin(binName, int(round(dirHist.GetBinContent(bin)*bw)),
            #{
                #"signal": (signal1.GetBinContent(bin))*bw,
                #"ttg": ttgHist.GetBinContent(bin)*bw,
                #"zg": zgHist.GetBinContent(bin)*bw,
                #"zz": zzHist.GetBinContent(bin)*bw,
                #"wwg": wwgHist.GetBinContent(bin)*bw,
                #"wzg": wzgHist.GetBinContent(bin)*bw,
                #"dy": dyHist.GetBinContent(bin)*bw,
                #"wjets": wjetsHist.GetBinContent(bin)*bw,
                #"tt": ttHist.GetBinContent(bin)*bw,
                #"singletop": singletopHist.GetBinContent(bin)*bw,
                #"wz": wzHist.GetBinContent(bin)*bw,
                #"ww": wwHist.GetBinContent(bin)*bw,
                #"wg": wgHist.GetBinContent(bin)*bw,
                #"zz4l": zz4lHist.GetBinContent(bin)*bw,
            #}, {
                #"ttgStat_"+binName: {"ttg": getDatacardUncertFromHist(ttgHist,bin)},
                #"zgStat_"+binName: {"zg": getDatacardUncertFromHist(zgHist,bin)},
                #"zzStat_"+binName: {"zz": getDatacardUncertFromHist(zzHist,bin)},
                #"wwgStat_"+binName: {"wwg": getDatacardUncertFromHist(wwgHist,bin)},
                #"wzgStat_"+binName: {"wzg": getDatacardUncertFromHist(wzgHist,bin)},
                #"dyStat_"+binName: {"dy": getDatacardUncertFromHist(dyHist,bin)},
                #"wjetsStat_"+binName: {"wjets": getDatacardUncertFromHist(wjetsHist,bin)},
                #"ttStat_"+binName: {"tt": getDatacardUncertFromHist(ttHist,bin)},
                #"singletopStat_"+binName: {"singletop": getDatacardUncertFromHist(singletopHist,bin)},
                #"wzStat_"+binName: {"wz": getDatacardUncertFromHist(wzHist,bin)},
                #"wwStat_"+binName: {"ww": getDatacardUncertFromHist(wwHist,bin)},
                #"wgStat_"+binName: {"wg": getDatacardUncertFromHist(wgHist,bin)},
                #"zz4lStat_"+binName: {"zz4l": getDatacardUncertFromHist(zz4lHist,bin)},
                #"signalStat_"+binName: {"signal": getDatacardUncertFromHist(signal1,bin)},
                #"lumi": {
                    #"signal": 1.026,
                    #"ttg": 1.026,
                    #"zg": 1.026,
                    #"zz": 1.026,
                    #"wwg": 1.026,
                    #"wzg": 1.026,
                    #"dy": 1.026,
                    #"wjets": 1.026,
                    #"tt": 1.026,
                    #"singletop": 1.026,
                    #"wz": 1.026,
                    #"ww": 1.026,
                    #"wg": 1.026,
                    #"zz4l": 1.026},
                #"SFDY": {
                    #"zg": (1.+sfDYErr),
                    #"dy": (1.+sfDYErr)},
                #"SFTT": {
                    #"tt": (1.+0.04),
                    #"ttg": (1.+0.04)},
                #"SFZZ": {
                    #"zz": (1.+sfZZErr),
                    #"zz4l": (1.+sfZZErr)},
                #"SFWZ": {
                    #"wz": (1.+sfWZErr)},
                #"dataMC": {
                    #"signal": 1.05,
                    #"ttg": 1.2,
                    #"zg": 1.2,
                    #"zz": 1.2,
                    #"wwg": 1.2,
                    #"wzg": 1.2,
                    #"dy": 1.2,
                    #"wjets": 1.2,
                    #"tt": 1.2,
                    #"singletop": 1.2,
                    #"wz": 1.2,
                    #"ww": 1.2,
                    #"wg": 1.2,
                    #"zz4l": 1.2},
            #}
        #)
    #dc.write("testDatacard.txt")
    #dc.write("limitCalculations/testDatacard.txt")
    
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
    #finalDistributionSignalHist("final_MC", allMC, "onZ/LL")
    finalDistributionSignalHist("final_MC", allMC, "onZMet150/LL")

main()
