import ROOT
from ROOT import *
from array import array
from include import *
import numpy as np

def getDatacardUncertFromHist(h,b):
    c = h.GetBinContent(b)
    return 1.+max(0,h.GetBinError(b)/c if c else 0)


def finalDistributionSignalHist(name, dirSet, dirDir):
    #style.divideByBinWidth = True
    style.divideByBinWidth = False

    #nBins = range(0,200,10)+[200, 250, 300, 350, 450, 600, 800]
    nBins = np.concatenate((np.arange(0,200,50),np.arange(200,300,100),np.arange(300,3000,1400)),axis=0)

    # direct stuff
    dirHist = aux.stdHist(dirSet, dirDir+"/met", nBins)
    style.additionalPoissonUncertainty = False
    aux.drawOpt(dirHist, "data")


    #gjetHist, gjetSyst, info = gjetPrediction(dirHist, preSet, zg+wg+ttg+wjets+ttjets_nlo+znunu, "met", nBins, weight, name+"_divByBinWidth" if style.divideByBinWidth else name)
    #gjetHist.SetLineColor(rwth.myLightBlue)
    #gjetHist.GetXaxis().SetTitle("#it{p}_{T}^{miss} (GeV)")


    #eHist = aux.stdHist(preSetElectron, preDirElectron+"/met", nBins)
    #eHist.GetXaxis().SetTitle("#it{p}_{T}^{miss} (GeV)")
    #eHist.Scale( 0.0267 if dirSet == data else 0.0154 )
    #eHist.SetLineColor(rwth.myYellow)
    #eSyst = aux.getSysHisto(eHist, 0.3)


    zgHist = aux.stdHist(zgamma, dirDir+"/met", nBins)
    ttgHist = aux.stdHist(ttgamma, dirDir+"/met", nBins)
    zzHist = aux.stdHist(zz, dirDir+"/met", nBins)
    wwgHist = aux.stdHist(wwgamma, dirDir+"/met", nBins)
    wzgHist = aux.stdHist(wzgamma, dirDir+"/met", nBins)
    dyHist = aux.stdHist(DYjets, dirDir+"/met", nBins)
    wjetsHist = aux.stdHist(wjets, dirDir+"/met", nBins)
    ttHist = aux.stdHist(tt, dirDir+"/met", nBins)

    zgHist.SetLineColor(ROOT.kGreen-3)
    ttgHist.SetLineColor(ROOT.kRed+1)
    zzHist.SetLineColor(ROOT.kYellow)
    wwgHist.SetLineColor(ROOT.kCyan-7)
    wzgHist.SetLineColor(ROOT.kMagenta+2)
    dyHist.SetLineColor(ROOT.kBlue-7)
    wjetsHist.SetLineColor(ROOT.kGreen+3)
    ttHist.SetLineColor(ROOT.kOrange+8)

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
    zgSyst = aux.getSysHisto(zgHist, mcSystUncert)
    ttgSyst = aux.getSysHisto(ttgHist, mcSystUncert)
    zzSyst = aux.getSysHisto(zzHist, mcSystUncert)
    wwgSyst = aux.getSysHisto(wwgHist, mcSystUncert)
    wzgSyst = aux.getSysHisto(wzgHist, mcSystUncert)
    dySyst = aux.getSysHisto(dyHist, mcSystUncert)
    wjetsSyst = aux.getSysHisto(wjetsHist, mcSystUncert)
    ttSyst = aux.getSysHisto(ttHist, mcSystUncert)

    #zgSyst = aux.addUncertaintiesQuadratic([zgSyst,zgPdfUnc,zgScaleUnc,zgJesUnc,zgPuUnc])
    #wgSyst = aux.addUncertaintiesQuadratic([wgSyst,wgPdfUnc,wgScaleUnc,wgJesUnc,wgPuUnc])
    #tgSyst = aux.addUncertaintiesQuadratic([tgSyst,tgPdfUnc,tgScaleUnc,tgJesUnc,tgPuUnc])

    #totStat = aux.addHists(gjetHist, eHist, zgHist, wgHist, tgHist)
    #totSyst = aux.addHists(gjetSyst, eSyst, zgSyst, wgSyst, tgSyst)
    totStat = aux.addHists(zgHist, ttgHist, zzHist, wwgHist, wzgHist,dyHist,wjetsHist,ttHist)
    totSyst = aux.addHists(zgSyst, ttgHist, zzHist, wwgHist, wzgHist,dyHist,wjetsHist,ttHist)

    signal1 = aux.stdHist(t5bbbbzg_1800_400, dirDir+"/met", nBins)
    signal2 = aux.stdHist(tching_400, dirDir+"/met", nBins)
    for h in signal1, signal2:
        aux.drawOpt(h, "signal")
    #    h.Add(totStat)
    signal1.SetLineColor(ROOT.kBlue+3)
    signal2.SetLineColor(ROOT.kRed-3)
    signal2.SetLineStyle(2)

    #signal1_pre = aux.createHistoFromDatasetTree(t5wg_1600_100, "met*{}".format(info["shift"]), weight, nBins, "tr_jControl/simpleTree")
    #signal1_pre.Scale(info["scale"])

    #totSyst = gjetSyst

    totUnc = aux.addHistUncert(totStat, totSyst)
    aux.drawOpt(totUnc, "totUnc")

    c = ROOT.TCanvas()
    m = multiplot.Multiplot()
    #if dirSet == data:
    if dirSet == dataDoubleSF:
        m.add(dirHist, "Data")
    else:
        m.add(dirHist, "Direct simulation")
#   m.add(signal1_pre, "contamination")
    m.add(signal1, "T5Zg 1800 400")
    m.add(signal2, "TChiNg 400")
    #m.addStack(eHist, "e#rightarrow#gamma")
    m.addStack(zgHist, "#gammaZ")
    m.addStack(ttgHist, "#gammat#bar{t}")
    m.addStack(zzHist, "ZZ")
    m.addStack(wwgHist, "WW#gamma")
    m.addStack(wzgHist, "WZ#gamma")
    m.addStack(dyHist, "DY")
    m.addStack(wjetsHist, "W+jets")
    m.addStack(ttHist, "t#bar{t}")


    m.add(totUnc, "Total uncertainty")
    m.maximum = 2.6*m.getMaximum()
    m.minimum = m.getMinimum()
    #if "final_lowEMHT" in name: m.minimum = 4e-2
    #if "final_highEMHT" in name: m.minimum = 2e-3
    #legInfo = "#it{H}_{T}^{#gamma} < 2TeV" if "lowEMHT" in name else "2TeV < #it{H}_{T}^{#gamma}"
    #if "ee" in name: legInfo += ", EE"
    #legInfo += ", |#Delta#phi|>0.3"
    legInfo = "DiMu" if "MM" in name else "DiEle"
    m.leg.SetHeader(legInfo)
    m.leg.SetY1(.56)
    m.leg.SetX1(.56)
    m.leg.SetX2(.99)


    m.Draw()

    # draw other labels

    l = ROOT.TLine()
    l.SetLineStyle(2)
    l.SetLineColor(ROOT.kGray+2)
    text = ROOT.TLatex()
    text.SetTextSize(0.8*text.GetTextSize())
    #l.DrawLine(100, 0, 100, totUnc.GetBinContent(totUnc.FindBin(100)))
    text.SetTextAngle(90)
    #text.DrawLatexNDC(.23,.315, "Normalization")
    text.SetTextAngle(0)
    #text.DrawLatexNDC(.311,.315, "Validation")
    if "final" in name:
        l.DrawLine(300, 0, 300, totUnc.GetBinContent(totUnc.FindBin(300)))



    r = ratio.Ratio("#scale[.9]{#lower[.24]{#splitline{Data/Pred.}{Bkg. frac.}}}", dirHist, totStat, totSyst)
    rMax = 1.5
    #if name == "final_lowEMHT": rMax = 1.6
    #if name == "final_highEMHT": rMax = 3.6
    #r.draw(0., rMax, m.getStack(), True)

    #aux.Label(sim= not dirSet==data, status="" if "allMC" not in name else "Private Work")
    aux.Label(sim= not dirSet==dataDoubleSF, status="" if "allMC" not in name else "Private Work")
    aux.save(name, normal=False, changeMinMax=False)

    #if name == "final_lowEMHT": dc = limitTools.MyDatacard()
    if name == "final_MC_MM": dc = limitTools.MyDatacard()
    #elif name == "final_highEMHT": dc = limitTools.MyDatacard("testDatacard.txt")
    else: return
    for bin in range(dirHist.GetNbinsX()-2, dirHist.GetNbinsX()+1):
        binName = "bin{}_{}".format(name.split("_")[1],bin)
        bw = dirHist.GetBinWidth(bin) if style.divideByBinWidth else 1
        dc.addBin(binName, int(round(dirHist.GetBinContent(bin)*bw)),
            {
                "signal": (signal1.GetBinContent(bin)-totStat.GetBinContent(bin))*bw,
                #"signal": (signal1.GetBinContent(bin)-totStat.GetBinContent(bin)-signal1_pre.GetBinContent(bin))*bw,
                "ttg": ttgHist.GetBinContent(bin)*bw,
                "zg": zgHist.GetBinContent(bin)*bw,
                "zz": zzHist.GetBinContent(bin)*bw,
                "wwg": wwgHist.GetBinContent(bin)*bw,
                "wzg": wzgHist.GetBinContent(bin)*bw,
                "dy": dyHist.GetBinContent(bin)*bw,
                "wjets": wjetsHist.GetBinContent(bin)*bw,
                "tt": ttHist.GetBinContent(bin)*bw,
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
                    "tt": 1.026},
                #"pu": {
                    #"wg": getDatacardUncertFromHist(wgPuUnc,bin),
                    #"zg": getDatacardUncertFromHist(zgPuUnc,bin),
                    #"tg": getDatacardUncertFromHist(tgPuUnc,bin)},
                #"jes": {
                    #"wg": getDatacardUncertFromHist(wgJesUnc,bin),
                    #"zg": getDatacardUncertFromHist(zgJesUnc,bin),
                    #"tg": getDatacardUncertFromHist(tgJesUnc,bin)},
                "dataMC": {
                    "signal": 1.05,
                    #"ttg": 1.025,
                    "ttg": 1.2,
                    "zg": 1.2,
                    "zz": 1.2,
                    "wwg": 1.2,
                    "wzg": 1.2,
                    "dy": 1.2,
                    "wjets": 1.2,
                    "tt": 1.2},
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
    
def main():
    #allMC = gjets+qcd+zg+wg+ttg+wjets+ttjets_nlo+znunu
    allMC = zgamma+ttgamma+zz+wwgamma+wzgamma+DYjets+wjets+tt
    allMC.label = "MC mix"
    #finalDistributionSignalHist("final_lowEMHT", data, "signal_lowEMHT", dataHt, data, "signal_lowEMHT_eControl")
    #finalDistributionSignalHist("allMC_lowEMHT", allMc, "signal_lowEMHT", allMc, allMc, "signal_lowEMHT_eControl")
    finalDistributionSignalHist("final_MC_MM", allMC, "onZMM")
    finalDistributionSignalHist("final_MC_EE", allMC, "onZEE")

main()
