#!/usr/bin/env python2

from include import *
import DatacardParser
import multiprocessing
import glob
import numpy as np
import pickle as pkl

contTT_tching = pkl.load(open("contamination/TT_tching.pkl", "rb"))
contTT_gmsb = pkl.load(open("contamination/TT_gmsb.pkl", "rb"))
contTT_t5zg = pkl.load(open("contamination/TT_t5zg.pkl", "rb"))
contTT_m1m2 = pkl.load(open("contamination/TT_m1m2.pkl", "rb"))
contTT_m1m3 = pkl.load(open("contamination/TT_m1m3.pkl", "rb"))
contTT_t6 = pkl.load(open("contamination/TT_t6.pkl", "rb"))
contDY_tching = pkl.load(open("contamination/DY_tching.pkl", "rb"))
contDY_gmsb = pkl.load(open("contamination/DY_gmsb.pkl", "rb"))
contDY_t5zg = pkl.load(open("contamination/DY_t5zg.pkl", "rb"))
contDY_m1m2 = pkl.load(open("contamination/DY_m1m2.pkl", "rb"))
contDY_m1m3 = pkl.load(open("contamination/DY_m1m3.pkl", "rb"))
contDY_t6 = pkl.load(open("contamination/DY_t6.pkl", "rb"))
contWZ_tching = pkl.load(open("contamination/WZ_tching.pkl", "rb"))
contWZ_gmsb = pkl.load(open("contamination/WZ_gmsb.pkl", "rb"))
contWZ_t5zg = pkl.load(open("contamination/WZ_t5zg.pkl", "rb"))
contWZ_m1m2 = pkl.load(open("contamination/WZ_m1m2.pkl", "rb"))
contWZ_m1m3 = pkl.load(open("contamination/WZ_m1m3.pkl", "rb"))
contWZ_t6 = pkl.load(open("contamination/WZ_t6.pkl", "rb"))
contZZ_tching = pkl.load(open("contamination/ZZ_tching.pkl", "rb"))
contZZ_gmsb = pkl.load(open("contamination/ZZ_gmsb.pkl", "rb"))
contZZ_t5zg = pkl.load(open("contamination/ZZ_t5zg.pkl", "rb"))
contZZ_m1m2 = pkl.load(open("contamination/ZZ_m1m2.pkl", "rb"))
contZZ_m1m3 = pkl.load(open("contamination/ZZ_m1m3.pkl", "rb"))
contZZ_t6 = pkl.load(open("contamination/ZZ_t6.pkl", "rb"))


"""
###############################################################################
# check consistency
###############################################################################
def checkHistogramConsistency2d(h1, h2, relTolerance=1e-6):
    n1 = h1.GetEntries()
    n2 = h2.GetEntries()
    if n1 != n2: print "Number of entries not the same", n1, n2
    nx1 = h1.GetNbinsX()
    nx2 = h2.GetNbinsX()
    ny1 = h1.GetNbinsY()
    ny2 = h2.GetNbinsY()
    if nx1 != nx2: print "Number of x bins is not the same:", nx1, nx2
    if ny1 != ny2: print "Number of y bins is not the same:", ny1, ny2
    for x,y in aux.loopH(h1):
        c1 = h1.GetBinContent(x,y)
        c2 = h2.GetBinContent(x,y)
        if c1+c2 and abs(c1-c2)/(c1+c2)/2 > relTolerance:
                         print "Not same bin content in bin {}:{}:".format(
                             x,y), c1, c2
def checkHistogramConsistency1d(h1, h2, relTolerance=1e-6):
    n1 = h1.GetEntries()
    n2 = h2.GetEntries()
    if n1 != n2: print "Number of entries not the same", n1, n2
    nx1 = h1.GetNbinsX()
    nx2 = h2.GetNbinsX()
    if nx1 != nx2: print "Number of x bins is not the same:", nx1, nx2
    for x in aux.loopH(h1):
        c1 = h1.GetBinContent(x)
        c2 = h2.GetBinContent(x)
        if c1+c2 and abs(c1-c2)/(c1+c2)/2 > relTolerance:
                         print "Not same bin content in bin {}:".format(
                             x), c1, c2
def checkConsistency(datacardFile, signalScan, treeFile):
    hScan = aux.getFromFile(signalScan, "Wg_1600_100/signal_lowEMHT/met")
    hScan.Scale(aux.intLumi*aux.getXsecSMSglu(1600))
    # hPlot = aux.getFromFile(treeFile, "signal_lowEMHT/met")
    hPlot = t5wg_1600_100.getHist("signal_lowEMHT/met")
    print "Plot vs scan differences"
    sameHists = checkHistogramConsistency1d(hScan, hPlot, 1e-3)
    xBins = [350, 450, 600]
    hScanRebinned = aux.rebinX(hScan, xBins)
    dc = limitTools.MyDatacard(datacardFile)
    print "Datacard versus scan:"
    print dc.exp["binlowEMHT_24"]["signal"], hScanRebinned.GetBinContent(1)
    print dc.exp["binlowEMHT_25"]["signal"], hScanRebinned.GetBinContent(2)
    print dc.exp["binlowEMHT_26"]["signal"], hScanRebinned.GetBinContent(3)
    return
###############################################################################
# end check consistency
###############################################################################
"""


def checkUpToDateInputSignal(inputSignal):
    t1 = os.path.getmtime(inputSignal)
    toCheck = []
    # if "T5" in inputSignal: toCheck = ["SMS-T5Wg_signalScan.root", "SMS-T5Wg_mGo2150To2500_signalScan.root"]
    if "T5" in inputSignal:
        toCheck = ["SMS-T5Wg_signalScan.root",
                   "SMS-T5Wg_mGo2150To2500_signalScan.root"]
    if "T6" in inputSignal:
        toCheck = ["SMS-T6Wg_signalScan.root",
                   "SMS-T6Wg_mSq1850To2150_signalScan.root"]
    tCheck = max(
        [0] + [os.path.getmtime(os.path.join(os.path.dirname(inputSignal), x)) for x in toCheck])
    if tCheck:
        if tCheck > t1:
            print "please rerun hadd"
    elif "TChi" in inputSignal:
        return
    else:
        print "please add files"


def proceedWithWeakScan(outputDir, scanName, xsecFile):

    c = ROOT.TCanvas()
    axisHist = ROOT.TH1F("", "", 10, 300, 1300)
    axisHist.SetStats(0)

    axisHist.Draw("axis")
    scanRes = {}
    # print outputDir
    # print glob.glob("{}/*.txt.limit".format(outputDir))
    for fname in glob.glob("{}/*.txt.limit".format(outputDir)):
        # print "huhu"
        m = re.match(".*_(\d+)_(\d+).txt.limit", fname)
        with open(fname) as f:
            scanRes[int(m.group(1))] = limitTools.infoFromOut(f.read())
        # print limitTools.infoFromOut(f.read())

    defaultGr = ROOT.TGraph(len(scanRes))
    graphs = dict((x, defaultGr.Clone(x))
                  for x in ["obs", "exp", "exp1up", "exp1dn", "exp2up", "exp2dn"])
    obsGr = ROOT.TGraph()
    expGr = ROOT.TGraph()
    expGr2Up = ROOT.TGraph()
    expGr2Dn = ROOT.TGraph()
    expGrUp = ROOT.TGraph()
    expGrDn = ROOT.TGraph()
    xsecGr = ROOT.TGraph()
    xsecGrUp = ROOT.TGraph()
    xsecGrDn = ROOT.TGraph()
    exp1sigma = ROOT.TGraphAsymmErrors()
    exp2sigma = ROOT.TGraphAsymmErrors()

    x = []
    y = []
    ymin1 = []
    ymin2 = []
    ymax1 = []
    ymax2 = []

    for i, m in enumerate(sorted(scanRes)):
        xsec, xsec_unc = aux.getXsecInfoSMS(m, xsecFile)
        xsecGr.SetPoint(i, m, xsec)
        xsecGrUp.SetPoint(i, m, xsec * (1 + xsec_unc))
        xsecGrDn.SetPoint(i, m, xsec * (1 - xsec_unc))
        obsGr.SetPoint(i, m, xsec * scanRes[m]["obs"])
        expR = scanRes[m]["exp"]
        # exp1sigma.SetPoint(i, m, xsec * expR)
        expGr.SetPoint(i, m, xsec * expR)
        expGr2Up.SetPoint(i, m, xsec * (scanRes[m]["exp2up"]))
        expGr2Dn.SetPoint(i, m, xsec * (scanRes[m]["exp2dn"]))
        expGrUp.SetPoint(i, m, xsec * (scanRes[m]["exp1up"]))
        expGrDn.SetPoint(i, m, xsec * (scanRes[m]["exp1dn"]))
        # x.append(m)
        # y.append(xsec * expR)
        # ymin1.append(xsec * (scanRes[m]["exp1dn"]))
        # ymin2.append(xsec * (scanRes[m]["exp2dn"]))
        # ymax1.append(xsec * (scanRes[m]["exp1up"]))
        # ymax2.append(xsec * (scanRes[m]["exp2up"]))
        # exp2sigma.SetPoint(i, m, xsec * expR)
        # exp1sigma.SetPointEYhigh(i, xsec * (scanRes[m]["exp1up"] - expR))
        # exp2sigma.SetPointEYhigh(i, xsec * (scanRes[m]["exp2up"] - expR))
        # exp1sigma.SetPointEYlow(i, xsec * (expR - scanRes[m]["exp1dn"]))
        # exp2sigma.SetPointEYlow(i, xsec * (expR - scanRes[m]["exp2dn"]))
        for name in graphs:
            graphs[name].SetPoint(i, m, xsec * scanRes[m][name])
    writeDict(graphs, outputDir + "/Graphs2d.root")

    # beautify
    obsGr.SetLineWidth(2)

    gs = ROOT.TGraphSmooth()

    # obsGr2 = ROOT.TGraph(gs.Approx(obsGr))
    obsGr2 = ROOT.TGraph(gs.SmoothSuper(obsGr))
    # obsGr2 = gs.Approx(obsGr)
    # obsGr2 = gs.Approx(obsGr)
    # obsGr2.SetDirectory(0)
    obsGr2.SetLineWidth(2)
    # obsGr2 = gs.SmoothSuper(obsGr)
    expGr2 = ROOT.TGraph(gs.SmoothSuper(expGr))
    expGr2Up2 = ROOT.TGraph(gs.SmoothSuper(expGr2Up))
    expGr2Dn2 = ROOT.TGraph(gs.SmoothSuper(expGr2Dn))
    expGr2Up = ROOT.TGraph(gs.SmoothSuper(expGrUp))
    expGr2Dn = ROOT.TGraph(gs.SmoothSuper(expGrDn))

    for i in range(expGr2.GetN()):
        x.append(expGr2.GetX()[i])
        ymax2.append(expGr2Up2.GetY()[i])
        ymin2.append(expGr2Dn2.GetY()[i])
    for i in range(expGr.GetN()):
        # x.append(expGr.GetX()[i])
        ymax1.append(expGr2Up.GetY()[i])
        ymin1.append(expGr2Dn.GetY()[i])

    n = expGr2Up.GetN()
    grshade2 = ROOT.TGraph(2 * n)
    for i in range(n):
        grshade2.SetPoint(i, x[i], ymax2[i])
        grshade2.SetPoint(n + i, x[n - i - 1], ymin2[n - i - 1])

    grshade = ROOT.TGraph(2 * n)
    for i in range(n):
        grshade.SetPoint(i, x[i], ymax1[i])
        grshade.SetPoint(n + i, x[n - i - 1], ymin1[n - i - 1])

    # xSecGrUp2 = gs.Approx(xSecGrUp)
    # xSecGrDown = gs.Approx(xSecGrDown)
    # xSecGr2 = gs.Approx(xSecGr)
    # xSecGrUp2 = gs.Approx(xsecGrUp)
    # xSecGrDown = gs.Approx(xsecGrDn)
    # xSecGr2 = gs.Approx(xsecGr)
    # gss = ROOT.TGraphSmooth()
    # exp1sigma = gss.Approx(exp1sigma,"linear")
    # exp1sigma = gss.SmoothLowess(exp1sigma)
    # exp2sigma2 = gs.Approx(exp2sigma)
    # exp2sigma = gs.Approx(exp2sigma)
    # exp2sigma = gs.SmoothSuper(exp2sigma)

    for g in xsecGr, xsecGrUp, xsecGrDn:
        g.SetLineColor(ROOT.kBlue)
    xsecGrUp.SetLineStyle(2)
    xsecGrDn.SetLineStyle(2)

    expGr2.SetLineStyle(2)
    expGr2.SetLineWidth(2)
    expGr2.SetLineColor(2)

    exp2sigma.SetFillColor(ROOT.kOrange)
    exp2sigma.SetLineColor(exp2sigma.GetFillColor())
    exp1sigma.SetFillColor(ROOT.kGreen + 1)
    exp1sigma.SetLineColor(2)
    exp1sigma.SetLineStyle(2)
    exp1sigma.SetLineWidth(2)

    # lsp_ = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{}}#scale[0.85]{_{1}}"
    lsp_ = "#lower[-0.12]{#tilde{G}}#lower[0.2]{#scale[0.85]{}}#scale[0.85]{_{1}}"
    # lsp_0 = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
    lsp_0 = "NLSP"
    lsp_pm = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}"

    # scanParticle = lsp_0 if "NG" in scanName else lsp_
    scanParticle = lsp_0 if "NG" in scanName else lsp_
    # exp2sigma.SetMaximum(2)
    # exp2sigma.SetTitle(
    #     ";m_{{{}}} (GeV);95% CL cross section upper limit (pb)".format(scanParticle))
    # exp2sigma.GetXaxis().SetLimits(300,1300)
    # exp1sigma.GetXaxis().SetLimits(300,1300)
    # xsecGrUp.GetXaxis().SetLimits(300,1300)
    # xsecGrDn.GetXaxis().SetLimits(300,1300)
    axisHist.SetTitle(
        ";m_{{{}}} (GeV);95% CL cross section upper limit (pb)".format(scanParticle))
    # draw
    # exp2sigma.Draw("ap3")
    # exp2sigma.Draw("ap3 same")
    # exp2sigma.Draw("p3 same")
    # exp1sigma.Draw("3 same")
    # exp1sigma.Draw("xc")
    # exp1sigma.Draw("xc same")

    # xsecGr.Draw("xc")
    # xsecGrUp.Draw("xc")
    # xsecGrDn.Draw("xc")
    # xsecGr.Draw("xc same")
    # xsecGrUp.Draw("xc same")
    # xsecGrDn.Draw("xc same")

    defaultGr.GetXaxis().SetLimits(300, 1300)
    # exp2sigma.GetXaxis().SetLimits(300, 1300)
    # exp1sigma.GetXaxis().SetLimits(300, 1300)
    xsecGrUp.GetXaxis().SetLimits(300, 1300)
    xsecGrDn.GetXaxis().SetLimits(300, 1300)
    obsGr.GetXaxis().SetLimits(300, 1300)
    obsGr2.GetXaxis().SetLimits(300, 1300)

    defaultGr.GetXaxis().SetRangeUser(300, 1300)
    # exp2sigma.GetXaxis().SetRangeUser(300, 1300)
    # exp1sigma.GetXaxis().SetRangeUser(300, 1300)
    xsecGrUp.GetXaxis().SetRangeUser(300, 1300)
    xsecGrDn.GetXaxis().SetRangeUser(300, 1300)
    obsGr.GetXaxis().SetRangeUser(300, 1300)
    obsGr2.GetXaxis().SetRangeUser(300, 1300)

    # ROOT.gPad.Update()

    # obsGr.Draw("xcp")  # no observed

    # grshade2.SetFillStyle(3013)
    grshade2.SetFillColor(ROOT.kOrange)
    grshade2.Draw("f")
    grshade.SetFillColor(ROOT.kGreen + 1)
    grshade.Draw("f")

    xsecGr.Draw("xc same")
    xsecGrUp.Draw("xc same")
    xsecGrDn.Draw("xc same")

    obsGr2.Draw("xcp")
    expGr2.Draw("xcp")

    # legend
    exp1sigmaClone = exp1sigma.Clone()
    exp1sigmaClone.SetLineColor(exp1sigma.GetFillColor())
    exp1sigmaClone.SetLineStyle(1)
    leg = ROOT.TLegend(.45, .59, .94, .92)
    leg.SetFillStyle(0)
    leg.AddEntry(obsGr, "Observed limit", "l")
    leg.AddEntry(exp2sigma, "Expected limit #pm 1(2) s.d._{experiment}", "f")
    leg.AddEntry(xsecGr, "Signal cross section #pm s.d._{theory}", "l")
    leg.Draw()

    leg2 = ROOT.TLegend(.45, .66, .94, .85)
    leg2.SetFillStyle(0)
    leg2.AddEntry(None, "", "")
    leg2.AddEntry(exp1sigmaClone, "", "f")
    leg2.AddEntry(None, "", "")
    leg2.Draw()

    leg3 = leg.Clone()
    leg3.Clear()
    leg3.AddEntry(None, "", "")
    leg3.AddEntry(exp1sigma, "", "l")
    leg3.AddEntry(None, "", "")
    leg3.Draw()

    leg4 = ROOT.TLegend(.45, .63, .94, .662)
    leg4.SetFillStyle(0)
    leg4.AddEntry(xsecGrUp, "", "l")
    leg4.AddEntry(xsecGrUp, "", "l")
    leg4.Draw()

    t = ROOT.TLatex()
    if "WG" in scanName:
        info = "pp#rightarrow%s(#tilde{G}#gamma) %s(#tilde{G}W^{#pm})" % (
            lsp_0, lsp_pm)
        t.DrawLatexNDC(.2, .18, info)
    if "NG" in scanName:
        # t.DrawLatexNDC(.18,.3, "pp#rightarrow{1}{1}/{1}{0}".format(lsp_0,lsp_pm))
        # t.DrawLatexNDC(.18,.25, "{}#rightarrow{}+soft".format(lsp_pm,lsp_0))
        # t.DrawLatexNDC(.18,.18, "{0}(#tilde{{G}}#gamma) {0}(#tilde{{G}}H/Z)".format(lsp_0))
        t.DrawLatexNDC(.24, .86, "#scale[0.76]{pp#rightarrow#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}/#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}}")
        t.DrawLatexNDC(.24, .81, "#scale[0.76]{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}#rightarrow#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}+soft}")
        t.DrawLatexNDC(.24, .76, "#scale[0.76]{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}(#tilde{G}#gamma) #lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}(#tilde{G}Z)}")
    aux.Label(sim=True)
    ROOT.gPad.SetLogy()
    axisHist.Draw("same axis")
    axisHist.GetXaxis().SetRangeUser(300, 1300)
    axisHist.GetXaxis().SetLimits(300, 1300)
    c.Update()
    c.Modified()
    ROOT.gPad.Modified()
    # c.SaveAs("testLIMIT.pdf")
    aux.save("{}_limit".format(scanName), log=False)


def getPointFromDir(name):
    m = re.match("(.*)_(.*)_(.*)", name)
    combi, m1, m2 = m.groups()
    m1, m2 = int(m1), int(m2)
    return combi, m1, m2


def writeDict(d, filename):
    f = ROOT.TFile(filename, "recreate")
    for name, ob in d.iteritems():
        if ob:
            ob.Write(name)
    f.Close()


def readDict(filename):
    f = ROOT.TFile(filename)
    tmpDir = f.GetDirectory(path)
    d = {}
    for element in tmpDir.GetListOfKeys():
        obj = element.ReadObj()
        obj = ROOT.gROOT.CloneObject(obj)
        d[element.GetName()] = obj
    return d


def writeSMSLimitConfig(infile, configName):
    text = """
HISTOGRAM {0} obs_hist
EXPECTED {0} exp exp1up exp1dn kRed kOrange
OBSERVED {0} obs obs1up obs1dn kBlack kGray
PRELIMINARY Private Work
LUMI {1:.1f}
ENERGY 13
""".format(infile, aux.intLumi / 1e3)
    # t2 = "OBSERVED ../../master/singlePhoton/PlotsSMS/config/SUS14004/2015-01-09-limits/SMS_T5wg/ROOT/SMS_T5wg_gluino_chi1_Exclusion_witXsecLimit.root Expected_limit Expected_limit_up Expected_limit_dn kBlack kGray"
    # text = "\n".join([l for l in text.split("\n") if "OBSERVED" not in l]+[t2, ""])
    with open(configName, "w+") as f:
        f.write(text)


def getXsecFile(name):
    xsecFile = ""
    if "T5" in name:
        xsecFile = "data/xSec_SMS_Gluino_13TeV.pkl"
    elif "T6" in name:
        xsecFile = "data/xSec_SMS_Squark_13TeV.pkl"
    elif "TChiWG" in name:
        xsecFile = "data/xSec_SMS_N2C1_13TeV.pkl"
    elif "TChiNG" in name:
        xsecFile = "data/xSec_SMS_TChiNG_13TeV.pkl"
    elif "M2" in name:
        xsecFile = "data/xsec_GGM_M1_M2.pkl"
    elif "M3" in name:
        xsecFile = "data/xsec_GGM_M1_M3.pkl"
    elif "GMSB" in name:
        xsecFile = "data/xSec_GMSB.pkl"
    else:
        print "Do not know which cross section belongs to", name
    return xsecFile


def getMultiScanName(inputSignal):
    if not "M1" in inputSignal and not "GMSB" in inputSignal:
        return os.path.basename(inputSignal).split("_")[0][4:]
    elif "GMSB" in inputSignal:
        return os.path.basename(inputSignal).split("_")[0]
    else:
        return os.path.basename(inputSignal).split("_")[0] + os.path.basename(inputSignal).split("_")[2][:2] + os.path.basename(inputSignal).split("_")[3][:2]


def getScanName(inputSignal, combi):
    return getMultiScanName(inputSignal).replace("Wg", combi)

# def getSignalUncertainties(inputSignal, dirname):


def getSignalUncertainties(inputSignal, dirname, datacardtemp):
    # print inputSignal,dirname
    # xBins = [200, 300, 500]
    # xBins = np.concatenate((np.arange(200,300,100),np.arange(300,600,200)),axis=0)
    # xBins = [200, 300]
    # xBins = [150,190,230]
    # xBins = [100,150,190,230]
    # xBins = [150,190,230]
    # xBins = [150,250]
    xBins = [150, 200]
    # xBins = [100,150,200]
    # print xBins
    out = {}
    # print dirname
    # print inputSignal

    if "Ng" in dirname:
        mass = dirname.split("_")[1]
        massstr2 = ""
    else:
        if "GMSB" in dirname:
            mass = dirname.split("_")[1]
            massstr2 = "_" + dirname.split("_")[2].split("/")[0]
        if "Zg" in dirname:
            mass = dirname.split("_")[1]
            massstr2 = "_" + dirname.split("_")[2].split("/")[0]
        if "GGM" in dirname:
            mass = dirname.split("_")[1]
            massstr2 = "_" + dirname.split("_")[2].split("/")[0]
    # print mass + massstr2
    f = ROOT.TFile(inputSignal)
    hNominal = f.Get(dirname + "/LL/nom/met")
    hNominalOld = f.Get(dirname + "/LL/nom/met")
    hNominal = aux.rebinX(hNominal, xBins)
    # aux.appendFlowBin( hNominal, under=False, over=True ) #yes or no? ask knut

    # for weight Histo:
    # weightHisto=f.Get("weightHisto_"+str(mass))
    topPtWeight = aux.getWeightForWeights(
        inputSignal, histoName="weightHisto_" + str(mass) + massstr2, whichWeight="pu_mc_toppt")
    ewkWeight = aux.getWeightForWeights(
        inputSignal, histoName="weightHisto_" + str(mass) + massstr2, whichWeight="pu_mc_ewk")
    ewkUpWeight = aux.getWeightForWeights(
        inputSignal, histoName="weightHisto_" + str(mass) + massstr2, whichWeight="pu_mc_ewkUp")
    ewkDnWeight = aux.getWeightForWeights(
        inputSignal, histoName="weightHisto_" + str(mass) + massstr2, whichWeight="pu_mc_ewkDn")
    nisrWeight = aux.getWeightForWeights(
        inputSignal, histoName="weightHisto_" + str(mass) + massstr2, whichWeight="pu_mc_nisr")
    nisrUpWeight = aux.getWeightForWeights(
        inputSignal, histoName="weightHisto_" + str(mass) + massstr2, whichWeight="pu_mc_nisrUp")
    nisrDnWeight = aux.getWeightForWeights(
        inputSignal, histoName="weightHisto_" + str(mass) + massstr2, whichWeight="pu_mc_nisrDn")
    # print "nbins",hNominal.GetNbinsX()

    # Scale because of Reweightings(TopPt, nISR....)

    hNominal.Scale(topPtWeight * ewkWeight * nisrWeight)

    hGenMet = f.Get(dirname + "/LL/genmet/met")
    hGenMet = aux.rebinX(hGenMet, xBins)
    hGenMet.Scale(topPtWeight * ewkWeight * nisrWeight)
    out["genmet"] = [1. + abs(hNominal.GetBinContent(b) - hGenMet.GetBinContent(b)) / 2 /
                     hNominal.GetBinContent(b) if hNominal.GetBinContent(b) else 1 for b in range(1, 3)]
    # out["genmet"] = [1.+abs(hNominal.GetBinContent(b)-hGenMet.GetBinContent(b))/2/hNominal.GetBinContent(b) if hNominal.GetBinContent(b) else 1 for b in range(1,4)]

    hPuUp = f.Get(dirname + "/LL/NoPUu/met")
    hPuDn = f.Get(dirname + "/LL/NoPUd/met")
    hPuUp.Scale(topPtWeight * ewkWeight * nisrWeight)
    hPuDn.Scale(topPtWeight * ewkWeight * nisrWeight)
    if hPuUp.GetEntries():
        hPuUp.Scale(hNominal.GetEntries() / hPuUp.GetEntries())
    if hPuDn.GetEntries():
        hPuDn.Scale(hNominal.GetEntries() / hPuDn.GetEntries())
    hPuUp = aux.rebinX(hPuUp, xBins)
    hPuDn = aux.rebinX(hPuDn, xBins)
    out["PU"] = [1. + abs(hPuUp.GetBinContent(b) - hPuDn.GetBinContent(b)) / 2 /
                 hNominal.GetBinContent(b) if hNominal.GetBinContent(b) else 1 for b in range(1, 3)]
    # out["PU"] = [1.+abs(hPuUp.GetBinContent(b)-hPuDn.GetBinContent(b))/2/hNominal.GetBinContent(b) if hNominal.GetBinContent(b) else 1 for b in range(1,4)]

    hJesUp = f.Get(dirname + "/LL/JESu/met")
    hJesDn = f.Get(dirname + "/LL/JESd/met")
    hJesUp = aux.rebinX(hJesUp, xBins)
    hJesDn = aux.rebinX(hJesDn, xBins)
    hJesUp.Scale(topPtWeight * ewkWeight * nisrWeight)
    hJesDn.Scale(topPtWeight * ewkWeight * nisrWeight)
    hJerUp = f.Get(dirname + "/LL/JERu/met")
    hJerDn = f.Get(dirname + "/LL/JERd/met")
    hJerUp = aux.rebinX(hJerUp, xBins)
    hJerDn = aux.rebinX(hJerDn, xBins)
    hJerUp.Scale(topPtWeight * ewkWeight * nisrWeight)
    hJerDn.Scale(topPtWeight * ewkWeight * nisrWeight)
    out["JES"] = [1. + (abs(hJesUp.GetBinContent(b) - hJesDn.GetBinContent(b)) / 2 /
                        hNominal.GetBinContent(b)) if hNominal.GetBinContent(b) else 1 for b in range(1, 3)]
    # out["JES"] = [1.+(abs(hJesUp.GetBinContent(b)-hJesDn.GetBinContent(b))/2/hNominal.GetBinContent(b)) if hNominal.GetBinContent(b) else 1 for b in range(1,4)]
    out["JER"] = [1. + (abs(hJerUp.GetBinContent(b) - hJerDn.GetBinContent(b)) / 2 /
                        hNominal.GetBinContent(b)) if hNominal.GetBinContent(b) else 1 for b in range(1, 3)]
    # out["JER"] = [1.+(abs(hJerUp.GetBinContent(b)-hJerDn.GetBinContent(b))/2/hNominal.GetBinContent(b)) if hNominal.GetBinContent(b) else 1 for b in range(1,4)]
    # out["jes"] = [1.+math.sqrt((abs(hJesUp.GetBinContent(b)-hJesDn.GetBinContent(b))/2/hNominal.GetBinContent(b))**2 \
    # + (abs(hJerUp.GetBinContent(b)-hJerDn.GetBinContent(b))/2/hNominal.GetBinContent(b))**2) if hNominal.GetBinContent(b) else 1 for b in range(1,4)]

    # hIsr = f.Get(dirname+"/met_isr")
    hIsrUp = f.Get(dirname + "/LL/ISRu/met")
    hIsrDn = f.Get(dirname + "/LL/ISRu/met")
    # hIsr.Scale(hNominal.Integral(0,-1)/hIsr.Integral(0,-1))
    # hIsrUp.Scale(hNominal.Integral(0,-1)/hIsrUp.Integral(0,-1))
    # hIsrDn.Scale(hNominal.Integral(0,-1)/hIsrDn.Integral(0,-1))
    hIsrUp.Scale(topPtWeight * ewkWeight * nisrUpWeight)
    hIsrDn.Scale(topPtWeight * ewkWeight * nisrDnWeight)
    # hIsr = aux.rebinX(hIsr, xBins)
    hIsrUp = aux.rebinX(hIsrUp, xBins)
    hIsrDn = aux.rebinX(hIsrDn, xBins)
    # out["ISR"] = [1.+abs(hIsrUp.GetBinContent(b)-hIsrDn.GetBinContent(b))/2/hNominal.GetBinContent(b) if hNominal.GetBinContent(b) else 1 for b in range(1,4)]
    out["ISR"] = [1. + abs(hIsrUp.GetBinContent(b) - hIsrDn.GetBinContent(b)) / 2 /
                  hNominal.GetBinContent(b) if hNominal.GetBinContent(b) else 1 for b in range(1, 3)]

    # hIsr = f.Get(dirname+"/met_isr")
    hEwkUp = f.Get(dirname + "/LL/EWKu/met")
    hEwkDn = f.Get(dirname + "/LL/EWKu/met")
    # hIsr.Scale(hNominal.Integral(0,-1)/hIsr.Integral(0,-1))
    # hIsrUp.Scale(hNominal.Integral(0,-1)/hIsrUp.Integral(0,-1))
    # hIsrDn.Scale(hNominal.Integral(0,-1)/hIsrDn.Integral(0,-1))
    hEwkUp.Scale(topPtWeight * ewkUpWeight * nisrWeight)
    hEwkDn.Scale(topPtWeight * ewkDnWeight * nisrWeight)
    # hIsr = aux.rebinX(hIsr, xBins)
    hEwkUp = aux.rebinX(hEwkUp, xBins)
    hEwkDn = aux.rebinX(hEwkDn, xBins)
    # out["EWK"] = [1.+abs(hEwkUp.GetBinContent(b)-hEwkDn.GetBinContent(b))/2/hNominal.GetBinContent(b) if hNominal.GetBinContent(b) else 1 for b in range(1,4)]
    out["EWK"] = [1. + abs(hEwkUp.GetBinContent(b) - hEwkDn.GetBinContent(b)) / 2 /
                  hNominal.GetBinContent(b) if hNominal.GetBinContent(b) else 1 for b in range(1, 3)]

    # scaleHists = dict([(i,f.Get(dirname+"/met_weight{}".format(i))) for i in range(1,9)])
    if not (("GMSB" in dirname) or ("GGM" in dirname)):
        scaleHists = dict(
            [(i, f.Get(dirname + "/LL/{}/met".format(i))) for i in range(0, 9)])
        for hname, h in scaleHists.iteritems():
            # h.Scale(hNominal.Integral(0,-1)/h.Integral(0,-1))
            h.Scale(topPtWeight * ewkWeight * nisrWeight)
            h.Scale(aux.getWeightForWeights(inputSignal, histoName="weightHisto_" +
                                            str(mass) + massstr2, whichWeight="pu_mc_pdf_" + str(hname)))
            scaleHists[hname] = aux.rebin(h, xBins)
        out["scale"] = [1. + (max([h.GetBinContent(b) for h in scaleHists.values()]) - min([h.GetBinContent(b)
                                                                                            for h in scaleHists.values()])) / 2. / hNominal.GetBinContent(b) if hNominal.GetBinContent(b) else 1 for b in range(1, 3)]
        # out["scale"] = [1.+(max([h.GetBinContent(b) for h in scaleHists.values()])-min([h.GetBinContent(b) for h in scaleHists.values()]))/2./hNominal.GetBinContent(b) if hNominal.GetBinContent(b) else 1 for b in range(1,4)]

    hLepSFUp = f.Get(dirname + "/LL/lepSFu/met")
    hLepSFDn = f.Get(dirname + "/LL/lepSFd/met")
    hLepSFUp = aux.rebinX(hLepSFUp, xBins)
    hLepSFDn = aux.rebinX(hLepSFDn, xBins)
    hLepSFUp.Scale(topPtWeight * ewkWeight * nisrWeight)
    hLepSFDn.Scale(topPtWeight * ewkWeight * nisrWeight)
    # 1.+abs(hLepSFUp.GetBinContent(bin)-hLepSFDn.GetBinContent(bin))/2./hNominal.GetBinContent(bin)
    out["lepSF"] = [1. + (abs(hLepSFUp.GetBinContent(b) - hLepSFDn.GetBinContent(b)) / 2 /
                          hNominal.GetBinContent(b)) if hNominal.GetBinContent(b) else 1 for b in range(1, 3)]
    # out["lepSF"]=[1.+(abs(hLepSFUp.GetBinContent(b)-hLepSFDn.GetBinContent(b))/2/hNominal.GetBinContent(b)) if hNominal.GetBinContent(b) else 1 for b in range(1,4)]

    hPhotonSFUp = f.Get(dirname + "/LL/photonSFu/met")
    hPhotonSFDn = f.Get(dirname + "/LL/photonSFd/met")
    hPhotonSFUp = aux.rebinX(hPhotonSFUp, xBins)
    hPhotonSFDn = aux.rebinX(hPhotonSFDn, xBins)
    hPhotonSFUp.Scale(topPtWeight * ewkWeight * nisrWeight)
    hPhotonSFDn.Scale(topPtWeight * ewkWeight * nisrWeight)
    # 1.+abs(signal1SFuHisto.GetBinContent(bin)-signal1lepSFdHisto.GetBinContent(bin))/2./signal1.GetBinContent(bin)
    out["lepSF"] = [1. + (abs(hPhotonSFUp.GetBinContent(b) - hPhotonSFDn.GetBinContent(b)) /
                          2 / hNominal.GetBinContent(b)) if hNominal.GetBinContent(b) else 1 for b in range(1, 3)]
    # out["lepSF"]=[1.+(abs(hPhotonSFUp.GetBinContent(b)-hPhotonSFDn.GetBinContent(b))/2/hNominal.GetBinContent(b)) if hNominal.GetBinContent(b) else 1 for b in range(1,4)]

    out["stat"] = [1. + hNominal.GetBinError(b) / hNominal.GetBinContent(
        b) if hNominal.GetBinContent(b) else 1 for b in range(1, 3)]
    # out["stat"] = [1.+hNominal.GetBinError(b)/hNominal.GetBinContent(b) if hNominal.GetBinContent(b) else 1 for b in range(1,4)]
    # out["stat"] = [1.+hNominal.GetBinError(b)/hNominal.GetBinContent(b) if hNominal.GetBinContent(b) else 1 for b in range(1,3)]
    # out["stat"] = [1.+hNominal.GetBinError(b)/hNominal.GetBinContent(b) if hNominal.GetBinContent(b) else 1 for b in range(1,5)]

    # hCont = f.Get(dirname+"/met_contamination")
    # hCont = aux.rebinX(hCont, xBins)
    # hCRTT = f.Get(dirname.replace("sig","CRTT")+"/EM/nom/met")
    # hCRTT.Scale(topPtWeight*ewkWeight*nisrWeight)
    # hCRTT=aux.rebinX(hCRTT, xBins)
    # hCRDY = f.Get(dirname.replace("sig","CRDY")+"/LL/nom/met")
    # hCRDY.Scale(topPtWeight*ewkWeight*nisrWeight)
    # hCRDY=aux.rebinX(hCRDY, xBins)
    # hCRZZ = f.Get(dirname.replace("sig","CRZZ")+"/LL/nom/met")
    # hCRZZ.Scale(topPtWeight*ewkWeight*nisrWeight)
    # hCRZZ=aux.rebinX(hCRZZ, xBins)
    # hCRWZ = f.Get(dirname.replace("sig","CRWZ")+"/LL/nom/met")
    # hCRWZ.Scale(topPtWeight*ewkWeight*nisrWeight)
    # hCRWZ=aux.rebinX(hCRWZ, xBins)

    obsYields = datacardtemp.returnObs()
    obsBin1 = obsYields["binMC_1"]
    obsBin2 = obsYields["binMC_2"]
    # obsBin3=obsYields["binMC_3"]

    # print dirname

    if("Ng") in dirname:
        # print obsBin1,obsBin2
        # print contTT_tching[int(mass)],contDY_tching[int(mass)],contZZ_tching[int(mass)],contWZ_tching[int(mass)]
        # print contTT_tching[int(mass)]+contDY_tching[int(mass)]+contZZ_tching[int(mass)]+contWZ_tching[int(mass)]
        contBin1 = (contTT_tching[int(mass)] + contDY_tching[int(mass)] +
                    contZZ_tching[int(mass)] + contWZ_tching[int(mass)]) * obsBin1
        contBin2 = (contTT_tching[int(mass)] + contDY_tching[int(mass)] +
                    contZZ_tching[int(mass)] + contWZ_tching[int(mass)]) * obsBin2
        # contBin3=(contTT_tching[int(mass)]+contDY_tching[int(mass)]+contZZ_tching[int(mass)]+contWZ_tching[int(mass)])*obsBin3
    if("GMSB") in dirname:
        mass2 = dirname.split("_")[2].split("/")[0]
        contBin1 = (contTT_gmsb[int(mass)][int(mass2)] + contDY_gmsb[int(mass)][int(mass2)] +
                    contZZ_gmsb[int(mass)][int(mass2)] + contWZ_gmsb[int(mass)][int(mass2)]) * obsBin1
        contBin2 = (contTT_gmsb[int(mass)][int(mass2)] + contDY_gmsb[int(mass)][int(mass2)] +
                    contZZ_gmsb[int(mass)][int(mass2)] + contWZ_gmsb[int(mass)][int(mass2)]) * obsBin2
        # contBin3=(contTT_gmsb[int(mass)][int(mass2)]+contDY_gmsb[int(mass)][int(mass2)]+contZZ_gmsb[int(mass)][int(mass2)]+contWZ_gmsb[int(mass)][int(mass2)])*obsBin3
    if("M1-200") in inputSignal:
        mass2 = dirname.split("_")[2].split("/")[0]
        contBin1 = (contTT_m1m2[int(mass)][int(mass2)] + contDY_m1m2[int(mass)][int(mass2)] +
                    contZZ_m1m2[int(mass)][int(mass2)] + contWZ_m1m2[int(mass)][int(mass2)]) * obsBin1
        contBin2 = (contTT_m1m2[int(mass)][int(mass2)] + contDY_m1m2[int(mass)][int(mass2)] +
                    contZZ_m1m2[int(mass)][int(mass2)] + contWZ_m1m2[int(mass)][int(mass2)]) * obsBin2
        # contBin3=(contTT_gmsb[int(mass)][int(mass2)]+contDY_gmsb[int(mass)][int(mass2)]+contZZ_gmsb[int(mass)][int(mass2)]+contWZ_gmsb[int(mass)][int(mass2)])*obsBin3
    if("M1-50") in inputSignal:
        mass2 = dirname.split("_")[2].split("/")[0]
        contBin1 = (contTT_m1m3[int(mass)][int(mass2)] + contDY_m1m3[int(mass)][int(mass2)] +
                    contZZ_m1m3[int(mass)][int(mass2)] + contWZ_m1m3[int(mass)][int(mass2)]) * obsBin1
        contBin2 = (contTT_m1m3[int(mass)][int(mass2)] + contDY_m1m3[int(mass)][int(mass2)] +
                    contZZ_m1m3[int(mass)][int(mass2)] + contWZ_m1m3[int(mass)][int(mass2)]) * obsBin2
        # contBin3=(contTT_gmsb[int(mass)][int(mass2)]+contDY_gmsb[int(mass)][int(mass2)]+contZZ_gmsb[int(mass)][int(mass2)]+contWZ_gmsb[int(mass)][int(mass2)])*obsBin3
    # if("T5") in dirname:
    if("T5") in inputSignal:
        mass2 = dirname.split("_")[2].split("/")[0]
        contBin1 = (contTT_t5zg[int(mass)][int(mass2)] + contDY_t5zg[int(mass)][int(mass2)] +
                    contZZ_t5zg[int(mass)][int(mass2)] + contWZ_t5zg[int(mass)][int(mass2)]) * obsBin1
        contBin2 = (contTT_t5zg[int(mass)][int(mass2)] + contDY_t5zg[int(mass)][int(mass2)] +
                    contZZ_t5zg[int(mass)][int(mass2)] + contWZ_t5zg[int(mass)][int(mass2)]) * obsBin2
    if("T6") in inputSignal:
        mass2 = dirname.split("_")[2].split("/")[0]
        contBin1 = (contTT_t6[int(mass)][int(mass2)] + contDY_t6[int(mass)][int(mass2)] +
                    contZZ_t6[int(mass)][int(mass2)] + contWZ_t6[int(mass)][int(mass2)]) * obsBin1
        contBin2 = (contTT_t6[int(mass)][int(mass2)] + contDY_t6[int(mass)][int(mass2)] +
                    contZZ_t6[int(mass)][int(mass2)] + contWZ_t6[int(mass)][int(mass2)]) * obsBin2
        # contBin3=(contTT_t5zg[int(mass)][int(mass2)]+contDY_t5zg[int(mass)][int(mass2)]+contZZ_t5zg[int(mass)][int(mass2)]+contWZ_t5zg[int(mass)][int(mass2)])*obsBin3

    # cont = [hCont.GetBinContent(b) for b in range(1,4)]
    # cont = [0. for b in range(1,4)]
    # cont = [0. for b in range(1,3)]
    # print contBin1,contBin2
    cont = [contBin1, contBin2]
    # cont = [contBin1,contBin2,contBin3]
    # cont = [(hCRTT.GetBinContent(b)+hCRDY.GetBinContent(b)+hCRZZ.GetBinContent(b)+hCRWZ.GetBinContent(b)) for b in range(1,3)]
    # cont = [0. for b in range(1,3)]
    # cont = [0. for b in range(1,5)]
    # if "TChiNG" in inputSignal:
    # acc = [hNominal.GetBinContent(b)*2. for b in range(1,4)]
    # else:
    # acc = [hNominal.GetBinContent(b) for b in range(1,4)]
    acc = [hNominal.GetBinContent(b) for b in range(1, 3)]
    # acc = [hNominal.GetBinContent(b) for b in range(1,4)]
    # acc[2] = hNominalOld.Integral(3,100)

    # acc = [hNominal.GetBinContent(b) for b in range(1,3)]
    # acc = [hNominal.GetBinContent(b) for b in range(1,5)]
    # print acc
    # print [hNominal.GetBinLowEdge(b) for b in range(1,3)]
    # print [hNominal.GetEntries(b) for b in range(1,3)]

    print "acc", acc
    print "cont", cont
    # print "out",out
    return acc, cont, out


def writeDataCards(outputDir, inputData, inputSignal, combi="", xsecFile=""):
    f = ROOT.TFile(inputSignal)
    print inputSignal, combi
    dirs = [k.GetName() for k in f.GetListOfKeys()
            if k.GetName().startswith(combi)]
    # print dirs
    # dirs = ["Zg_1800_400","Zg_800_500"] # cross check
    # dirs = ["Zg_1500_400"] # cross check
    # dirs = ["Ng_400_0"] # cross check
    # dirs = ["Ng_600_0"] # cross check
    # binNames = ["binlowEMHT_24", "binlowEMHT_25", "binlowEMHT_26", "binhighEMHT_24", "binhighEMHT_25", "binhighEMHT_26"]
    # binNames = ["met"]
    # binNames = ["binMC_9","binMC_10"]
    # binNames = ["binMC_3","binMC_4","binMC_5"]
    # binNames = ["binMC_6","binMC_7","binMC_8"]
    # binNames = ["binMC_6","binMC_7"]
    # binNames = ["binMC_1","binMC_2"]
    binNames = ["binMC_1", "binMC_2", "binMC_3"]
    # binNames = ["binMC_2","binMC_3","binMC_4","binMC_5"]

    dc = limitTools.MyDatacard(inputData)
    for d in dirs:
        # print d
        combi2, m1, m2 = getPointFromDir(d)
        newD = combi2 + "_" + str(m2) + "_" + str(m1)
        newM2 = m1
        newM1 = m2
        # print newD
        if not "M1" in inputSignal and not "GMSB" in inputSignal:
            xsec = aux.getXsecInfoSMS(m1, xsecFile)[0]
        else:
            # print m1,m2
            xsec = aux.getXsecInfoGGM(m1, m2, xsecFile)[0]
        # acc, cont, syst = getSignalUncertainties(inputSignal, d+"/signal")
        # acc, cont, syst = getSignalUncertainties(inputSignal, d+"/sig")
        acc, cont, syst = getSignalUncertainties(inputSignal, d + "/sig", dc)
        # acc, cont, syst = getSignalUncertainties(inputSignal, d+"/signal_TChiNG_gg")
        # acc, cont, syst = getSignalUncertainties(inputSignal, d+"/signal_TChiNG_gz")
        # acc, cont, syst = getSignalUncertainties(inputSignal, d+"/signal_TChiNG_zz")
        # acc, cont, syst = getSignalUncertainties(inputSignal, d+"/HIGHsignal")
        # acc2, cont2, syst2 = getSignalUncertainties(inputSignal, d+"/signal_highEMHT")
        # acc.extend(acc2)
        # cont.extend(cont2)
        for a, b in syst.iteritems():
            # syst[a] = syst[a]+syst2[a]
            syst[a] = syst[a]

        # print xsec
        # print acc

        acc = [a * aux.intLumi * xsec for a in acc]
        # cont = [a*aux.intLumi*xsec for a in cont]
        cont = [a for a in cont]

        # print "acc", acc,cont
        print "xsec", xsec

        acc = [a - b for a, b in zip(acc, cont)]

        # print "final acc",acc

        systs = {}
        for sName, valueList in syst.iteritems():
            if sName == "stat":
                for binName, unc in zip(binNames, valueList):
                    systs["signalStat_{}".format(binName)] = {binName: unc}
            else:
                systs[sName] = dict(zip(binNames, valueList))

        dc.newSignal(dict(zip(binNames, acc)), systs)
        # if not "GMSB" in inputSignal:
        # dcName = "{}/{}.txt".format(outputDir, d)
        dcName = "{}/{}.txt".format(outputDir, d)
        # else:
        # dcName = "{}/{}.txt".format(outputDir, newD)
        dc.write(dcName)


def callMultiCombine(outputDir):
    files = glob.glob("{}/*.txt".format(outputDir))
    p = multiprocessing.Pool()
    p.map(limitTools.callCombine, files)


def callMultiSignificance(outputDir):
    files = glob.glob("{}/*.txt".format(outputDir))
    p = multiprocessing.Pool()
    p.map(limitTools.callCombineSignificance, files)


def clearWrongCombineOutputs(outputDir):
    files = glob.glob("{}/*.txt.limit".format(outputDir))
    for fname in files:
        with open(fname) as f:
            rInfo = limitTools.infoFromOut(f.read())
        if rInfo["obs"] == 0:
            os.remove(fname)
    files = glob.glob("{}/*.txt.significance".format(outputDir))
    for fname in files:
        with open(fname) as f:
            rInfo = limitTools.significanceFromOut(f.read())
        if not rInfo:
            os.remove(fname)


def build2dGraphsLimit(outputDir):
    files = glob.glob("{}/*.txt.limit".format(outputDir))
    defaultGr = ROOT.TGraph2D(len(files))
    # defaultHist = ROOT.TH2F("", "", 26, 200, 1500, 26, 200, 1500)
    # defaultHist = ROOT.TH2F("", "", 46, 200, 2500, 46, 200, 2500)  # t5zg
    defaultHist = ROOT.TH2F("", "", 50, 0, 2500, 50, 0, 2500)  # t5zg....
    # defaultHist = ROOT.TH2F("", "", 60, 1000, 2500, 108, 100, 3000)
    # defaultHist = ROOT.TH2F("", "", 30, 1000, 2500, 480, 100, 2500)
    # defaultHist = ROOT.TH2F("", "", 30, 1000, 2500, 240, 100, 2500)
    # defaultHist = ROOT.TH2F("","",150,0,1500,150,0,1500)
    # defaultHist = ROOT.TH2F("", "", 40, 200, 1200, 40, 200, 1200)  # gmsb
    # defaultHist = ROOT.TH2F("", "", 60, 200, 1200, 60, 200, 1200)  # gmsb2
    # defaultHist = ROOT.TH2F("", "", 80, 200, 1200, 80, 200, 1200)  # gmsb2
    # defaultHist = ROOT.TH2F("", "", 10, 400, 1400, 15, 0, 1500)  # t6
    # defaultHist = ROOT.TH2F("","",32,205,1005,32,215,1015)
    # defaultHist = ROOT.TH2F("","",36,200,1100,36,200,1100)
    # defaultHist = ROOT.TH2F("","",18,200,1100,18,200,1100)
    graphs = dict((x, defaultGr.Clone(x))
                  for x in ["obs", "exp", "exp1up", "exp1dn", "exp2up", "exp2dn"])
    hists = dict((x, defaultHist.Clone(x))
                 for x in ["obs", "exp", "exp1up", "exp1dn", "exp2up", "exp2dn"])
    for g in graphs.values():
        g.SetDirectory(0)
    for ifile, _file in enumerate(files):
        m = re.match(".*_(\d+)_(\d+).txt.limit", _file)
        m1 = int(m.group(1))
        m2 = int(m.group(2))
        if "GMSB" in outputDir:
            m2 = int(m.group(1))
            m1 = int(m.group(2))
        with open(_file) as f:
            rInfo = limitTools.infoFromOut(f.read())
            # print rInfo
        for name, gr in graphs.iteritems():
            graphs[name].SetPoint(ifile, m1, m2, rInfo[name])
            hists[name].Fill(m1, m2, rInfo[name])
    writeDict(graphs, outputDir + "/saved_graphs2d_limit.root")
    writeDict(hists, outputDir + "/saved_hists2d_limit.root")
    return graphs


def build2dGraphsSignificance(outputDir):
    files = glob.glob("{}/*.txt.significance".format(outputDir))
    defaultGr = ROOT.TGraph2D(len(files))
    defaultGr.SetDirectory(0)
    for ifile, _file in enumerate(files):
        m = re.match(".*_(\d+)_(\d+).txt.significance", _file)
        m1 = int(m.group(1))
        m2 = int(m.group(2))
        with open(_file) as f:
            r = limitTools.significanceFromOut(f.read())
        defaultGr.SetPoint(ifile, m1, m2, r)
    writeDict({"significance": defaultGr}, outputDir +
              "/saved_graphs2d_significance.root")


def latexScanName(scanName):
    if scanName == "T6gg":
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        return "pp #rightarrow #tilde{q}#tilde{q}, #tilde{q} #rightarrow q%s, %s #rightarrow #gamma#tilde{G}" % (lsp_s, lsp_s)
    elif scanName == "T6Wg":
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        return "pp #rightarrow #tilde{q}#tilde{q}, #tilde{q} #rightarrow q%s, %s #rightarrow #gamma/W^{#pm}#tilde{G}" % (lsp_s, lsp_s)
    elif scanName == "T5gg":
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        return "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow qq%s, %s #rightarrow #gamma#tilde{G}" % (lsp_s, lsp_s)
    elif scanName == "T5Wg":
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        return "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow qq%s, %s #rightarrow #gamma/W^{#pm}#tilde{G}" % (lsp_s, lsp_s)
    return scanName


def drawSignalScanHist(h, scanName, saveName):
    smsScan = sms.sms(scanName)
    # print scanName
    style.style2d()
    h.SetTitle("")
    if "T5" in scanName:
        h.GetXaxis().SetTitle("m_{#tilde{g}} (GeV)")
    if "T6" in scanName:
        h.GetXaxis().SetTitle("m_{#tilde{q}} (GeV)")
    if "gg" in scanName:
        h.GetYaxis().SetTitle("m_{#tilde{#chi}^{0}_{1}} (GeV)")
    if "Wg" in scanName:
        h.GetYaxis().SetTitle("m_{#tilde{#chi}^{0/#pm}_{1}} (GeV)")
    hname = h.GetName()
    # print hname
    if hname == "significance":
        h.SetMinimum(-3)
        h.SetMaximum(3)
        style.setPaletteBWR()
        # ROOT.gStyle.SetPalette(ROOT.kTemperatureMap)
        h.GetXaxis().SetRangeUser(smsScan.Xmin, smsScan.Xmax)
        h.GetZaxis().SetNdivisions(6, 0, 0)
        h.GetZaxis().SetTitle("Significance (s.d.)")
    elif hname == "xsec":
        h.Scale(1000)
        h.GetZaxis().SetTitle("signal cross section (fb)")
    else:
        uncerts = {"isr": "ISR uncert", "jes": "Jet energy uncert", "pu": "Pile-up",
                   "scale": "Renorm.+Factor. scale uncert", "genMet": "FastSim p_{T}^{miss}"}
        a = h.GetName().split("_")
        if len(a) == 3:
            var, emhtSel, metSel = h.GetName().split("_")
        else:
            if len(a) == 4:
                var, emhtSel, metSel, dummy1 = h.GetName().split("_")
            if len(a) == 5:
                var, emhtSel, metSel, dummy1, dummy2 = h.GetName().split("_")
        if var in uncerts.keys():
            h.SetMinimum(0)
            h.SetMaximum(0.1)
        emhtSel = emhtSel.replace(
            "binlowEMHT", "700-2000").replace("binhighEMHT", "2000-#infty")
        metSel = metSel.replace(
            "24", "350-450").replace("25", "450-600").replace("26", "600-#infty")
        var = uncerts[var] if var in uncerts else var
        h.GetZaxis().SetTitle("{} ({},{})".format(var, emhtSel, metSel))
    h.GetZaxis().SetTitleOffset(1.42)
    h.GetYaxis().SetTitleOffset(h.GetYaxis().GetTitleOffset() * .95)

    c = ROOT.TCanvas()
    h.Draw("colz")
    aux.drawDiagonal(h, smsScan.Xmin)
    aux.Label2D(info="#scale[.76]{{{}}}".format(
        latexScanName(scanName)), status="Preliminary")
    aux.save("{}_{}".format(scanName, saveName))
    style.defaultStyle()


def drawSignalScanGraph(h, scanName, saveName):
    smsScan = sms.sms(scanName)
    # print scanName
    # h.SaveAs("shit.root")
    s = style.style2d()
    h.SetTitle("")
    if "T5" in scanName:
        h.GetXaxis().SetTitle("m_{#tilde{g}} (GeV)")
    if "T6" in scanName:
        h.GetXaxis().SetTitle("m_{#tilde{q}} (GeV)")
    if "gg" in scanName:
        h.GetYaxis().SetTitle("m_{#tilde{#chi}^{0}_{1}} (GeV)")
    if "Wg" in scanName:
        h.GetYaxis().SetTitle("m_{#tilde{#chi}^{0/#pm}_{1}} (GeV)")
    if "Zg" in scanName:
        h.GetYaxis().SetTitle("m_{#tilde{#chi}^{0}_{1}} (GeV)")
    if "GMSB" in scanName:
        h.GetYaxis().SetTitle("m_{#tilde{W}} (GeV)")
        h.GetXaxis().SetTitle("m_{#tilde{B}} (GeV)")
    hname = h.GetName()
    print hname
    if hname == "significance":
        h.SetMinimum(-3)
        h.SetMaximum(3)
        style.setPaletteBWR()
        # ROOT.gStyle.SetPalette(ROOT.kTemperatureMap)
        h.GetXaxis().SetRangeUser(smsScan.Xmin, smsScan.Xmax)
        h.GetZaxis().SetNdivisions(6, 0, 0)
        h.GetZaxis().SetTitle("Significance (s.d.)")
    elif hname == "xsec":
        h.Scale(1000)
        h.GetZaxis().SetTitle("signal cross section (fb)")
    else:
        uncerts = {"ISR": "ISR uncert.",
                   "JES": "Jet energy scale uncert.",
                   "JER": "Jet energy resolution uncert.",
                   "PU": "Pile-up uncert.",
                   "scale": "Renorm.+Factor. scale uncert.",
                   "genmet": "FastSim p_{T}^{miss} uncert.",
                   "PDF": "PDF uncert.",
                   "lepSF": "Lepton ID and Reco. uncert.",
                   "photonSF": "Photon ID and Reco. uncert.",
                   "signalStat_binMC_1": "stat. uncert.",
                   "signalStat_binMC_2": "stat. uncert."}
        a = h.GetName().split("_")
        if len(a) == 3:
            var, emhtSel, metSel = h.GetName().split("_")
        else:
            if len(a) == 4:
                var, emhtSel, metSel, dummy1 = h.GetName().split("_")
            if len(a) == 5:
                var, emhtSel, metSel, dummy1, dummy2 = h.GetName().split("_")
                var = var + "_" + dummy1 + "_" + dummy2
        # print var, len(a)
        if var in uncerts.keys():
            h.SetMinimum(0)
            # h.SetMaximum(0.1)
            # h.SetMaximum(h.GetZmax() * 1.1)
        emhtSel = emhtSel.replace(
            "binlowEMHT", "700-2000").replace("binhighEMHT", "2000-#infty")
        metSel = metSel.replace(
            "24", "350-450").replace("25", "450-600").replace("26", "600-#infty")
        if metSel == "1":
            metSel = metSel.replace("1", "150 - 200")
        if metSel == "2":
            metSel = metSel.replace("2", "200-#infty")
        var = uncerts[var] if var in uncerts else var
        # h.GetZaxis().SetTitle("{} ({},{})".format(var, emhtSel, metSel))
        # h.GetZaxis().SetTitle(
        # "#scale[.85]{} ({},{})".format(var, emhtSel, metSel))
        # h.GetZaxis().SetTitle("#scale[.85]{}{} ({})}".format(var, metSel))
        h.GetZaxis().SetTitle("{} ({})".format(var, metSel))
    # h.GetZaxis().SetTitleOffset(1.42)
    h.GetZaxis().SetTitleOffset(2.05)
    s.SetPadRightMargin(0.25)
    # h.GetYaxis().SetTitleOffset(h.GetYaxis().GetTitleOffset() * .95)
    h.GetYaxis().SetTitleOffset(h.GetYaxis().GetTitleOffset() * 1.1)

    if "Zg" in scanName:
        h.GetXaxis().SetRangeUser(800, 2500)
        h.GetYaxis().SetRangeUser(0, 2500)
        h.GetYaxis().SetTitle("m_{#tilde{#chi}^{0}_{1}} (GeV)")
        h.GetXaxis().SetTitle("m_{#tilde{g}} (GeV)")
    if "GMSB" in scanName:
        h.GetXaxis().SetRangeUser(200, 1000)
        h.GetYaxis().SetRangeUser(200, 1000)
        h.GetYaxis().SetTitle("m_{#tilde{W}} (GeV)")
        h.GetXaxis().SetTitle("m_{#tilde{B}} (GeV)")

    h.GetYaxis().SetTitleOffset(h.GetYaxis().GetTitleOffset() * .9)

    h.SetMinimum(h.GetZmin() * 0.95)
    h.SetMaximum(h.GetZmax() * 1.05)

    c = ROOT.TCanvas()
    # h.Draw("COLZ")
    # c.SetPhi(360)
    # c.SetTheta(90)
    # h.SetMarkerStyle(21)
    # h.SetMarkerSize(1.5)

    print var, h.GetZmin(), h.GetZmax()

    h.Draw("COLZ")

    c.Update()
    c.Modified()
    # h.Draw("")
    # h.SaveAs("shit2.root")
    # c.SaveAs("shit3.root")
    # c.SaveAs("shit3.pdf")
    # aux.drawDiagonal(h, smsScan.Xmin)
    # aux.Label2D(info="#scale[.76]{{{}}}".format(
    #     latexScanName(scanName)), status="Preliminary")
    aux.Label2D(info="#scale[.76]{{{}}}".format(
        latexScanName(scanName)), status="", sim=True)
    # aux.Label2D(info="#scale[.76]{{{}}}".format(
    #     latexScanName(scanName)), status="Private Work", sim=True)
    # aux.save("{}_{}".format(scanName,saveName))
    # aux.save("{}_{}".format(scanName,saveName),changeMinMax=False)
    aux.save("{}_{}".format(scanName, saveName), endings=[
             ".pdf"], changeMinMax=False)
    # aux.save("{}_{}".format(scanName, saveName), endings=[
    #          ".pdf", ".root"], changeMinMax=False)
    style.defaultStyle()


def drawSignalScanGraph1D(h, scanName, saveName):
    smsScan = sms.sms(scanName)
    # print scanName
    # style.style2d()
    h.SetTitle("")
    hname = h.GetName()
    # print hname
    if hname == "significance":
        h.SetMinimum(-3)
        h.SetMaximum(3)
        style.setPaletteBWR()
        # ROOT.gStyle.SetPalette(ROOT.kTemperatureMap)
        h.GetXaxis().SetRangeUser(smsScan.Xmin, smsScan.Xmax)
        h.GetZaxis().SetNdivisions(6, 0, 0)
        h.GetZaxis().SetTitle("Significance (s.d.)")
    elif hname == "xsec":
        h.Scale(1000)
        h.GetZaxis().SetTitle("signal cross section (fb)")
    else:
        uncerts = {"ISR": "ISR uncert.",
                   "EWK": "ISR uncert.",
                   "JES": "Jet energy scale uncert.",
                   "JER": "Jet energy resolution uncert.",
                   "PU": "Pile-up uncert.",
                   "scale": "Renorm.+Factor. scale uncert.",
                   "genmet": "FastSim p_{T}^{miss} uncert.",
                   "PDF": "PDF uncert.",
                   "lepSF": "Lepton ID and Reco. uncert.",
                   "photonSF": "Photon ID and Reco. uncert.",
                   "signalStat_binMC_1": "stat. uncert.",
                   "signalStat_binMC_2": "stat. uncert."}
        a = h.GetName().split("_")
        if len(a) == 3:
            var, emhtSel, metSel = h.GetName().split("_")
        else:
            if len(a) == 4:
                var, emhtSel, metSel, dummy1 = h.GetName().split("_")
            if len(a) == 5:
                var, emhtSel, metSel, dummy1, dummy2 = h.GetName().split("_")
        if var in uncerts.keys():
            h.SetMinimum(0)
            h.SetMaximum(1)
            # h.SetMaximum(0.1)
        emhtSel = emhtSel.replace(
            "binlowEMHT", "700-2000").replace("binhighEMHT", "2000-#infty")
        metSel = metSel.replace(
            "24", "350-450").replace("25", "450-600").replace("26", "600-#infty")
        if metSel == "1":
            metSel = metSel.replace("1", "150 - 200")
        if metSel == "2":
            metSel = metSel.replace("2", "200-#infty")
        var = uncerts[var] if var in uncerts else var
        # h.GetZaxis().SetTitle("{} ({},{})".format(var, emhtSel, metSel))
        # h.GetZaxis().SetTitle(
        # "#scale[.85]{} ({},{})".format(var, emhtSel, metSel))
        # h.GetZaxis().SetTitle("#scale[.85]{}{} ({})}".format(var, metSel))
        h.GetYaxis().SetTitle("{} ({})".format(var, metSel))
    # h.GetYaxis().SetTitleOffset(1.42)
    # h.GetYaxis().SetTitleOffset(h.GetYaxis().GetTitleOffset() * .95)
    h.GetYaxis().SetTitleOffset(h.GetYaxis().GetTitleOffset() * 1.1)

    if "Zg" in scanName:
        h.GetYaxis().SetTitle("m_{#tilde{#chi}^{0}_{1}} (GeV)")
        h.GetXaxis().SetTitle("m_{#tilde{g}} (GeV)")
    if "T6Zg" in scanName:
        h.GetYaxis().SetTitle("m_{#tilde{#chi}^{0}_{1}} (GeV)")
        h.GetXaxis().SetTitle("m_{#tilde{t}} (GeV)")
    if "GMSB" in scanName:
        h.GetYaxis().SetTitle("m_{#tilde{W}} (GeV)")
        h.GetXaxis().SetTitle("m_{#tilde{B}} (GeV)")
    if "NG" in scanName:
        # h.GetYaxis().SetTitle("m_{#tilde{W}} (GeV)")
        # h.GetXaxis().SetTitle("m_{#tilde{#chi}^{0}_{1}} (GeV)")
        h.GetXaxis().SetTitle("m_{NLSP} (GeV)")

    h.GetYaxis().SetTitleOffset(h.GetYaxis().GetTitleOffset() * .9)
    c = ROOT.TCanvas()

    print var, ROOT.TMath.MinElement(
        h.GetN(), h.GetY()), ROOT.TMath.MaxElement(h.GetN(), h.GetY())
    min = ROOT.TMath.MinElement(h.GetN(), h.GetY())
    max = ROOT.TMath.MaxElement(h.GetN(), h.GetY())
    h.SetMinimum(min - (0.1 * min))
    h.SetMaximum(max + (0.1 * max))

    h.Draw("A*")
    # aux.drawDiagonal(h, smsScan.Xmin)
    # aux.Label2D(info="#scale[.76]{{{}}}".format(
    #     latexScanName(scanName)), status="Preliminary")
    aux.Label(info="#scale[.76]{{{}}}".format(
        latexScanName(scanName)), status="Private Work", sim=True)
    # aux.save("{}_{}".format(scanName,saveName))
    # aux.save("{}_{}".format(scanName,saveName),changeMinMax=False)
    aux.save("{}_{}".format(scanName, saveName), endings=[
             ".pdf"], changeMinMax=False)
    # aux.save("{}_{}".format(scanName, saveName), endings=[
    #          ".pdf", ".root"], changeMinMax=False)
    style.defaultStyle()


def drawLimitInput(outputDir, scanName, xsecFile):
    files = glob.glob("{}/*.txt".format(outputDir))
    dc = limitTools.MyDatacard(files[0])
    defaultGr = ROOT.TGraph2D(len(files))
    graphs = dict((x, defaultGr.Clone(x))
                  for x in ["Acceptance_" + b for b in dc.bins] + ["xsec"])
    # uncerts = "isr", "jes", "pu", "scale", "genMet"
    # graphs.update(dict([(x+"_"+y,defaultGr.Clone(x+"_"+y)) for x in uncerts for y in dc.bins]))
    for ifile, f in enumerate(files):
        m = re.match(".*_(\d+)_(\d+).txt", f)
        m1 = int(m.group(1))
        m2 = int(m.group(2))
        if not "M1" in outputDir and not "GMSB" in outputDir:
            xsec = aux.getXsecInfoSMS(m1, xsecFile)[0]
        else:
            xsec = aux.getXsecInfoGGM(m1, m2, xsecFile)[0]
        # print m1, m2, xsec, xsecFile
        if "GMSB" in outputDir:
            graphs["xsec"].SetPoint(ifile, m2, m1, xsec)
        else:
            graphs["xsec"].SetPoint(ifile, m1, m2, xsec)
        # print xsec
        dc = limitTools.MyDatacard(f)
        systDict = dict([(l[0], l) for l in dc.systs])
        for b in dc.bins:
            if "GMSB" in outputDir:
                graphs["Acceptance_" +
                       b].SetPoint(ifile, m2, m1, dc.exp[b]["signal"] / (xsec * aux.intLumi))
            else:
                if xsec == 0.:
                    graphs["Acceptance_" + b].SetPoint(ifile, m1, m2, 0.)
                else:
                    graphs["Acceptance_" + b].SetPoint(
                        ifile, m1, m2, dc.exp[b]["signal"] / (xsec * aux.intLumi))

            # for unc in uncerts:
                # graphs[unc+"_"+b].SetPoint(ifile, m1, m2, systDict[unc][4][b]["signal"]-1)
    for name, gr in graphs.iteritems():
        h = gr.GetHistogram()
        drawSignalScanHist(h, scanName, name)


def drawSignalUncertainties(outputDir, scanName, xsecFile):
    files = glob.glob("{}/*.txt".format(outputDir))
    dc = limitTools.MyDatacard(files[0])
    # defaultGr = ROOT.TGraph2D(len(files))
    defaultGr = ROOT.TGraph2D()
    # print len(files)
    # graphs = dict((x,defaultGr.Clone(x)) for x in ["statError_"+b for b in dc.bins]+["xsec"])
    # uncs=["signalStat_binMC_1","signalStat_binMC_2","JES","JER","PU","EWK","ISR","genmet","lepSF","photonSF","scale","pdf"]
    if "GMSB" in scanName:
        uncs = ["signalStat_binMC_1", "signalStat_binMC_2", "JES", "JER",
                "PU", "genmet", "lepSF", "photonSF", "scale", "pdf"]
        # uncs = ["signalStat_binMC_1", "signalStat_binMC_2"]
        # uncs = ["signalStat_binMC_1"]
    else:
        uncs = ["signalStat_binMC_1", "signalStat_binMC_2", "JES", "JER",
                "PU", "ISR", "genmet", "lepSF", "photonSF", "scale", "pdf"]
    # binar = ["binMC_1"]
    # uncs=["signalStat_binMC_1","signalStat_binMC_2","JES","JER","PU","genmet","lepSF","photonSF","scale","pdf"]
    # uncs=[uncs[0]]
    # uncs=[uncs[5]]
    # graphs = dict((x,defaultGr.Clone(x)) for x in ["statError_"+b for b in dc.bins])
    # uncerts = "isr", "jes", "pu", "scale", "genMet"
    # graphs.update(dict([(x+"_"+y,defaultGr.Clone(x+"_"+y)) for x in uncerts for y in dc.bins]))
    graphs = dict([(x + "_" + y, defaultGr.Clone(x + "_" + y))
                   for x in uncs for y in dc.bins])
    # graphs = dict([(x + "_" + y, defaultGr.Clone(x + "_" + y))
    #                for x in uncs for y in binar])
    # graphs = dict([(x + "_" + y, ROOT.TGraph(x + "_" + y))
    #                for x in uncs for y in dc.bins])
    # print graphs
    for ifile, f in enumerate(files):
        m = re.match(".*_(\d+)_(\d+).txt", f)
        m1 = int(m.group(1))
        m2 = int(m.group(2))
        # if not "M1" in outputDir and not "GMSB" in outputDir:
        # xsec = aux.getXsecInfoSMS(m1, xsecFile)[0]
        # else:
        # xsec = aux.getXsecInfoGGM(m1,m2,xsecFile)[0]
        # if "GMSB" in outputDir:
        # graphs["xsec"].SetPoint(ifile, m2, m1, xsec)
        # else:
        # graphs["xsec"].SetPoint(ifile, m1, m2, xsec)
        # print m1, m2, ifile
        dc = limitTools.MyDatacard(f)
        systDict = dict([(l[0], l) for l in dc.systs])

        for b in dc.bins:
            # for b in binar:
            # for line in dc.systs:
                # if ("signalStat_"+b in line[0]):
                    # if "GMSB" in outputDir:
                        # graphs["statError_"+b].SetPoint(ifile, m2, m1, line[4][b]["signal"]-1.)
                    # else:
                        # graphs["statError_"+b].SetPoint(ifile, m1, m2, line[4][b]["signal"]-1.)

            for unc in uncs:
                if "GMSB" in outputDir:
                    # print m1, m2, ifile, systDict[unc][4][b]["signal"], unc, b
                    # print m1, m2, ifile, unc, b
                    if (not systDict.get(unc)):
                        continue
                    if systDict[unc][4][b]["signal"] < 1.:
                        # print systDict[unc][4]
                        # print m2,m1
                        # graphs[unc + "_" + b].SetPoint(ifile, m2, m1, 0.)
                        graphs[unc + "_" + b].SetPoint(ifile, m2, m1, 0.001)
                        # graphs[unc + "_" + b].SetPoint(ifile, m1, m2, 0.)
                    else:
                        # print graphs[unc + "_" +
                        #              b], ifile, m2, m1, (systDict[unc][4][b]["signal"] - 1.)
                        graphs[unc + "_" + b].SetPoint(ifile, m2,
                                                       m1, (systDict[unc][4][b]["signal"] - 1.))
                        # graphs[unc + "_" + b].SetPoint(ifile, m1,
                        #                                m2, systDict[unc][4][b]["signal"] - 1.)
                else:
                    # print "ERROR"
                    # graphs[unc+"_"+b].SetPoint(ifile, m1, m2, systDict[unc][4][b]["signal"]-1.)
                    if (not systDict.get(unc)):
                        # continue
                        graphs[unc + "_" + b].SetPoint(ifile, m1, m2, 0.)
                    else:
                        if systDict[unc][4][b]["signal"] < 1.:
                            # print systDict[unc][4]
                            # print m2,m1
                            graphs[unc + "_" + b].SetPoint(ifile, m1, m2, 0.)
                        else:
                            graphs[unc + "_" + b].SetPoint(
                                ifile, m1, m2, systDict[unc][4][b]["signal"] - 1.)
    writeDict(graphs, "testSHIT.root")
    for name, gr in graphs.iteritems():
        # h = gr.GetHistogram()
        drawSignalScanGraph(gr, scanName, name)


def drawSignalUncertainties1D(outputDir, scanName, xsecFile):
    files = glob.glob("{}/*.txt".format(outputDir))
    dc = limitTools.MyDatacard(files[0])
    defaultGr = ROOT.TGraph(len(files))
    # graphs = dict((x,defaultGr.Clone(x)) for x in ["statError_"+b for b in dc.bins]+["xsec"])
    # uncs=["signalStat_binMC_1","signalStat_binMC_2","JES","JER","PU","EWK","ISR","genmet","lepSF","photonSF","scale","pdf"]
    uncs = ["signalStat_binMC_1", "signalStat_binMC_2", "JES", "JER",
            "PU", "ISR", "EWK", "genmet", "lepSF", "photonSF", "scale", "pdf"]
    # uncs=["signalStat_binMC_1","signalStat_binMC_2","JES","JER","PU","genmet","lepSF","photonSF","scale","pdf"]
    # uncs=[uncs[0]]
    # uncs=[uncs[5]]
    # graphs = dict((x,defaultGr.Clone(x)) for x in ["statError_"+b for b in dc.bins])
    # uncerts = "isr", "jes", "pu", "scale", "genMet"
    # graphs.update(dict([(x+"_"+y,defaultGr.Clone(x+"_"+y)) for x in uncerts for y in dc.bins]))
    graphs = dict([(x + "_" + y, defaultGr.Clone(x + "_" + y))
                   for x in uncs for y in dc.bins])
    for ifile, f in enumerate(files):
        m = re.match(".*_(\d+)_(\d+).txt", f)
        m1 = int(m.group(1))
        m2 = int(m.group(2))
        # if not "M1" in outputDir and not "GMSB" in outputDir:
        # xsec = aux.getXsecInfoSMS(m1, xsecFile)[0]
        # else:
        # xsec = aux.getXsecInfoGGM(m1,m2,xsecFile)[0]
        # if "GMSB" in outputDir:
        # graphs["xsec"].SetPoint(ifile, m2, m1, xsec)
        # else:
        # graphs["xsec"].SetPoint(ifile, m1, m2, xsec)
        dc = limitTools.MyDatacard(f)
        systDict = dict([(l[0], l) for l in dc.systs])
        for b in dc.bins:
            # for line in dc.systs:
                # if ("signalStat_"+b in line[0]):
                    # if "GMSB" in outputDir:
                        # graphs["statError_"+b].SetPoint(ifile, m2, m1, line[4][b]["signal"]-1.)
                    # else:
                        # graphs["statError_"+b].SetPoint(ifile, m1, m2, line[4][b]["signal"]-1.)

            for unc in uncs:
                if "GMSB" in outputDir:
                    if (not systDict.get(unc)):
                        continue
                    if systDict[unc][4][b]["signal"] < 1.:
                        # print systDict[unc][4]
                        # print m2,m1
                        graphs[unc + "_" + b].SetPoint(ifile, m2, m1, 0.)
                    else:
                        graphs[unc + "_" + b].SetPoint(ifile, m2,
                                                       m1, systDict[unc][4][b]["signal"] - 1.)
                else:
                    # graphs[unc+"_"+b].SetPoint(ifile, m1, m2, systDict[unc][4][b]["signal"]-1.)
                    if (not systDict.get(unc)):
                        # continue
                        # graphs[unc + "_" + b].SetPoint(ifile, m1, m2, 0.)
                        graphs[unc + "_" + b].SetPoint(ifile, m1, 0.)
                    else:
                        if systDict[unc][4][b]["signal"] < 1.:
                            # print systDict[unc][4]
                            # print m2,m1
                            # graphs[unc + "_" + b].SetPoint(ifile, m1, m2, 0.)
                            graphs[unc + "_" + b].SetPoint(ifile, m1, 0.)
                        else:
                            graphs[unc + "_" + b].SetPoint(
                                # ifile, m1, m2, systDict[unc][4][b]["signal"] - 1.)
                                ifile, m1,  systDict[unc][4][b]["signal"] - 1.)

    for name, gr in graphs.iteritems():
        # h = gr.GetHistogram()
        drawSignalScanGraph1D(gr, scanName, name)


def getXsecLimitHistDelaunay(gr):
    grScaled = scaleObsWithXsec(gr)
    grScaled.SetNpx(500)
    grScaled.SetNpy(500)
    h = grScaled.GetHistogram()
    h = h.Clone(aux.randomName())
    return h


def build2dGraphs(outputDir, xsecFile):
    build2dGraphsLimit(outputDir)
    # build2dGraphsSignificance(outputDir)


def getObsUncertainty(gr2d, xsecFile):
    gr2dup = gr2d.Clone(aux.randomName())
    gr2ddn = gr2d.Clone(aux.randomName())
    gr2dup.SetDirectory(0)
    gr2ddn.SetDirectory(0)
    points = [(gr2d.GetX()[i], gr2d.GetY()[i], gr2d.GetZ()[i])
              for i in range(gr2d.GetN())]
    for ip, (x, y, z) in enumerate(points):
        if "M1" in xsecFile or "GMSB" in xsecFile:
            if "GMSB" in xsecFile:
                xsec, unc = aux.getXsecInfoGGM(y, x, xsecFile)
            else:
                xsec, unc = aux.getXsecInfoGGM(x, y, xsecFile)
        else:
            xsec, unc = aux.getXsecInfoSMS(x, xsecFile)
        gr2dup.SetPoint(ip, x, y, z * (1 - unc / 100))
        gr2ddn.SetPoint(ip, x, y, z * (1 + unc / 100))
    obsUp = limitTools.getContour(gr2dup)
    obsUp.SetName("obs1up")
    obsDn = limitTools.getContour(gr2ddn)
    obsDn.SetName("obs1dn")
    return {"obs1up": obsUp, "obs1dn": obsDn}


def getObsUncertainty2(gr2d, xsecFile):
    # gr2d = limitTools.getFineHisto2(gr2d)
    defaultHist = ROOT.TH2F("", "", 40, 200, 1200, 40, 200, 1200)  # gmsb
    gr2dup = gr2d.Clone(aux.randomName())
    gr2ddn = gr2d.Clone(aux.randomName())
    gr2dup.SetDirectory(0)
    gr2ddn.SetDirectory(0)
    h2dup = defaultHist.Clone(aux.randomName())
    h2ddn = defaultHist.Clone(aux.randomName())
    h2dup.SetDirectory(0)
    h2ddn.SetDirectory(0)
    points = [(gr2d.GetX()[i], gr2d.GetY()[i], gr2d.GetZ()[i])
              for i in range(gr2d.GetN())]
    for ip, (x, y, z) in enumerate(points):
        if "M1" in xsecFile or "GMSB" in xsecFile:
            if "GMSB" in xsecFile:
                xsec, unc = aux.getXsecInfoGGM(y, x, xsecFile)
                realUnc = aux.getXsecuncGMSB(y, x, xsecFile)[1]
                # unc = unc * 100
                # print xsec, unc
                unc = realUnc
                # print unc
                # print unc, realUnc
            else:
                xsec, unc = aux.getXsecInfoGGM(x, y, xsecFile)
        else:
            xsec, unc = aux.getXsecInfoSMS(x, xsecFile)
        gr2dup.SetPoint(ip, x, y, z * (1 - unc / 100))
        gr2ddn.SetPoint(ip, x, y, z * (1 + unc / 100))
        if "GMSB" in xsecFile:
            h2dup.Fill(x, y, z * (1 - unc * 100. / 100))
            h2ddn.Fill(x, y, z * (1 + unc * 100. / 100))
        else:
            h2dup.Fill(x, y, z * (1 - unc / 100))
            h2ddn.Fill(x, y, z * (1 + unc / 100))
    # obsUp = limitTools.getContour2(limitTools.getFineHisto(gr2dup))
    obsUp = limitTools.getContour(limitTools.getFineHisto2(gr2dup))
    obsUp.SetName("obs1up")
    # obsDn = limitTools.getContour2(h2ddn)
    # obsDn = limitTools.getContour2(limitTools.getFineHisto(gr2ddn))
    obsDn = limitTools.getContour(limitTools.getFineHisto2(gr2ddn))
    obsDn.SetName("obs1dn")
    return {"obs1up": obsUp, "obs1dn": obsDn}


def interpolateAlongY(h2):
    for xbin, ybin in aux.loopH(h2):
        if not h2.GetBinContent(xbin, ybin):
            ybinUp, cUp = 0, 0
            for ybinUp in range(ybin, h2.GetNbinsY()):
                cUp = h2.GetBinContent(xbin, ybinUp)
                if cUp:
                    break
            ybinDn, cDn = 0, 0
            for ybinDn in range(ybin, 0, -1):
                cDn = h2.GetBinContent(xbin, ybinDn)
                if cDn:
                    break
            if not cUp or not cDn:
                continue
            h2.SetBinContent(xbin, ybin, ((cUp - cDn) * ybin + (cDn * ybinUp -
                                                                cUp * ybinDn)) / (ybinUp - ybinDn))  # linear interpolation


def interpolateAlongX(h2):
    for xbin, ybin in aux.loopH(h2):
        if not h2.GetBinContent(xbin, ybin):
            xbinUp, cUp = 0, 0
            for xbinUp in range(xbin, h2.GetNbinsX()):
                cUp = h2.GetBinContent(xbinUp, ybin)
                if cUp:
                    break
            xbinDn, cDn = 0, 0
            for xbinDn in range(xbin, 0, -1):
                cDn = h2.GetBinContent(xbinDn, ybin)
                if cDn:
                    break
            if not cUp or not cDn:
                continue
            h2.SetBinContent(xbin, ybin, ((cUp - cDn) * xbin + (cDn * xbinUp -
                                                                cUp * xbinDn)) / (xbinUp - xbinDn))  # linear interpolation


def interpolateHoles(h2):
    for xbin, ybin in aux.loopH(h2):
        c = h2.GetBinContent(xbin, ybin)
        if abs(c) > 1e-5:
            continue
        cNord = h2.GetBinContent(xbin, ybin + 1)
        cSout = h2.GetBinContent(xbin, ybin - 1)
        cWest = h2.GetBinContent(xbin + 1, ybin)
        cEast = h2.GetBinContent(xbin - 1, ybin)
        if cNord and cSout and cWest and cEast:
            h2.SetBinContent(xbin, ybin, sum(
                [cNord, cSout, cWest, cEast]) / 4.)
    return h2


def getXsecLimitHist(gr2d, h, xsecFile):
    h.SetDirectory(0)  # or the next line will overwrite the hist?
    points = [(gr2d.GetX()[i], gr2d.GetY()[i], gr2d.GetZ()[i])
              for i in range(gr2d.GetN())]
    for x, y, z in points:
        if "M1" in xsecFile or "GMSB" in xsecFile:
            if "GMSB" in xsecFile:
                xsec = aux.getXsecInfoGGM(y, x, xsecFile)[0]
            else:
                xsec = aux.getXsecInfoGGM(x, y, xsecFile)[0]
        else:
            xsec = aux.getXsecInfoSMS(x, xsecFile)[0]
        # print x, y, z, xsec, z * xsec
        # h.SetBinContent(h.FindBin(x, y), z * xsec)
    return h
    # return limitTools.getFineHisto2(h)


def getXsecLimitHist2(gr2d, h, xsecFile):
    gr = ROOT.TGraph2D(gr2d)
    h.SetDirectory(0)  # or the next line will overwrite the hist?
    # points = [(gr.GetX()[i], gr.GetY()[i], gr.GetZ()[i])
    #           for i in range(gr.GetN())]
    for i in range(gr.GetN()):
        x = gr.GetX()[i]
        y = gr.GetY()[i]
        z = gr.GetZ()[i]
    # for x, y, z in points:
        if "M1" in xsecFile or "GMSB" in xsecFile:
            if "GMSB" in xsecFile:
                xsec = aux.getXsecInfoGGM(y, x, xsecFile)[0]
            else:
                xsec = aux.getXsecInfoGGM(x, y, xsecFile)[0]
        else:
            xsec = aux.getXsecInfoSMS(x, xsecFile)[0]
        # print x, y, z, xsec, z * xsec
        gr.SetPoint(i, x, y, z * xsec)
        # h.SetBinContent(h.FindBin(x, y), z * xsec)
    grToUse = limitTools.getFineHisto2(gr)
    points = [(gr.GetX()[i], gr.GetY()[i], gr.GetZ()[i])
              for i in range(gr.GetN())]
    for x, y, z in points:
        h.SetBinContent(h.FindBin(x, y), z)
    # h.Smooth()
    return h
    # return limitTools.getFineHisto2(h)


def getHistForModel(model):
    h = ROOT.TH2F()
    if "T5" in model:
        # h = ROOT.TH2F("", "", 35, 775, 2525, 500, 0, 2500)
        # h = ROOT.TH2F("", "", 30, 1000, 2500, 480, 100, 2500)
        h = ROOT.TH2F("", "", 140, 1000, 2600, 480, 100, 2600)  # toUse
        # h = ROOT.TH2F("", "", 280, 0, 2600, 480, 0, 2600)
    elif "T6" in model:
        # h = ROOT.TH2F("", "", 24, 975, 2175, 220, 0, 2200)
        h = ROOT.TH2F("", "", 15, 0, 1500, 30, 0, 1500)
    elif "BR" in model:
        h = ROOT.TH2F("", "", 25, 300, 1300, 50, 0, 100)
    elif "M1M2" in model:
        h = ROOT.TH2F("", "", 26, 200, 1500, 26, 200, 1500)
    elif "GMSB" in model:
        # h = ROOT.TH2F("", "", 36, 200, 1100, 36, 200, 1100)
        h = ROOT.TH2F("", "", 140, 200, 1100, 140, 200, 1100)
    else:
        print "Not specified model", model
    smsScan = sms.sms(model)
    # h.SetTitle("{};{};{};95% CL upper limit on cross section (fb)".format(
    #     model, smsScan.sParticle, smsScan.LSP))
    h.SetTitle("{};{};{};".format(
        model, smsScan.sParticle, smsScan.LSP))
    h.SetMinimum(0)
    return h


# def smoothContour(gr, neighbors=5, sigma=.5):  # use this for T5
def smoothContour(gr, neighbors=5, sigma=.01):
    # def smoothContour(gr, neighbors=10, sigma=.001):
    # def smoothContour(gr, neighbors=1, sigma=0.01):
    # def smoothContour(gr, neighbors=2, sigma=9):
    fgaus = ROOT.TF1("fgaus", "gaus", -10, 10)
    fgaus.SetParameters(1, 0, 1)
    weights = [fgaus.Eval(i * sigma) for i in range(neighbors)]
    out = gr.Clone(aux.randomName())
    out.Set(0)
    n = gr.GetN()
    Xs = [gr.GetX()[i] for i in range(n)]
    Ys = [gr.GetY()[i] for i in range(n)]
    # n = Ys.index(max(Ys)) + 1
    Xs = Xs[0:n]
    Ys = Ys[0:n]
    for i, (x, y) in enumerate(zip(Xs, Ys)):
        pNeigh = min(neighbors, i + 1, n - i)
        newX, ws, newY = 0, 0, 0
        for j in range(pNeigh):
            if j:
                newX += (Xs[i - j] + Xs[i + j]) * weights[j]
                newY += (Ys[i - j] + Ys[i + j]) * weights[j]
                ws += 2 * weights[j]
            else:
                newX += x * weights[0]
                newY += y * weights[0]
                ws += weights[0]
        # print i, x, y, newX, newY, ws
        out.SetPoint(i, newX / ws, newY / ws)
    return out
    # return gr


def build1dGraphs(outputDir, xsecFile, scanName):
    graphs = readDict(outputDir + "/saved_graphs2d_limit.root")
    hists = readDict(outputDir + "/saved_hists2d_limit.root")
    if "GMSB" in scanName:
        # for name, gr in graphs.iteritems():
        #     gr.SetNpx(100)
        #     gr.SetNpy(100)
        # toDraw = dict([(name, limitTools.getContour2(limitTools.getShiftedHisto(gr)))
        #                for name, gr in graphs.iteritems()])
        # toDraw = dict([(name, limitTools.getContourTry(gr))
        # toDraw = dict([(name, limitTools.getContour2(limitTools.getFineHisto(gr)))
        # toDraw = dict([(name, limitTools.getContour2(limitTools.getFineHisto2(gr)))
        toDraw = dict([(name, limitTools.getContour(limitTools.getFineHisto2(gr)))  # yep
                       # toDraw = dict([(name, limitTools.getContour(gr))

                       # toDraw = dict([(name, limitTools.getContour2(gr.GetHistogram()))
                       for name, gr in graphs.iteritems()])
        # toDraw = dict([(name, limitTools.getContourAll(gr))
        #                for name, gr in graphs.iteritems()])
        # toDraw = dict([(name, limitTools.getContour(gr))
        #                for name, gr in graphs.iteritems()])
        # toDraw = dict([(name, limitTools.getContour2(h))
        #                for name, h in hists.iteritems()])
    else:
        toDraw = dict([(name, limitTools.getContour(gr))
                       for name, gr in graphs.iteritems()])
    writeDict(toDraw, "testcontours.root")
    if "GMSB" in scanName:
        toDraw.update(getObsUncertainty2(graphs["obs"], xsecFile))
    else:
        toDraw.update(getObsUncertainty(graphs["obs"], xsecFile))
    # print toDraw
    if not "GMSB" in scanName:
        # if "GMSB" in scanName:
        # a = 1
        toDraw = dict([(n, smoothContour(gr)) for n, gr in toDraw.iteritems()])
    else:
        gs = ROOT.TGraphSmooth()
        # obsGr2 = ROOT.TGraph(gs.Approx(obsGr))
        # toDraw = dict([(n, ROOT.TGraph(gs.Approx(gr)) for n, gr in toDraw.iteritems()])
        # toDraw = dict([(n, dummy(gr)) for n, gr in toDraw.iteritems()])
        toDraw = dict([(n, smoothContour(gr)) for n, gr in toDraw.iteritems()])
    # toDraw["obs_hist"] = getXsecLimitHistDelaunay(graphs["obs"])
    # toDraw["obs_hist"] = getXsecLimitHist(
    toDraw["obs_hist"] = getXsecLimitHist2(
        # graphs["obs"], getHistForModel(scanName), xsecFile)
        graphs["obs"], getHistForModel(scanName), xsecFile)
    interpolateAlongY(toDraw["obs_hist"])
    interpolateAlongX(toDraw["obs_hist"])
    writeDict(toDraw, outputDir + "/saved_graphs1d_limit.root")


def dummy(a):
    gs = ROOT.TGraphSmooth()
    return ROOT.TGraph(gs.Approx(a))
    # return ROOT.TGraph(gs.Approx(a, "constant"))
    # return ROOT.TGraph(gs.SmoothLowess(a))
    # return ROOT.TGraph(gs.SmoothKern(a, "box"))


def drawSignificance(outputDir, scanName):
    graphs = readDict(outputDir + "/saved_graphs2d_significance.root")
    h = graphs["significance"].GetHistogram()
    h.SetName("significance")
    smsScan = sms.sms(scanName)
    h.SetTitle("{};{};{};Significance (s.d.)".format(
        scanName, smsScan.sParticle, smsScan.LSP))
    interpolateHoles(h)
    writeDict({"significance": h}, outputDir +
              "/saved_graphs1d_significance.root")
    drawSignalScanHist(h, scanName, "significance")


def signalScan(combi, version, inputData, inputSignal):
    print "help"
    # checkUpToDateInputSignal(inputSignal)
    scanName = getScanName(inputSignal, combi)
    scanNameOld = getScanName(inputSignal, combi)
    if "BR" in combi:
        scanName = scanName + "_BR"
    print scanName
    # outputDir = "limitCalculations/{}_{}_BRlimits_gg".format(scanName, version)
    # outputDir = "limitCalculations/{}_{}_BRlimits_gz".format(scanName, version)
    # outputDir = "limitCalculations/{}_{}_BRlimits_zz".format(scanName, version)
    # outputDir = "limitCalculations/{}_{}_BRlimits_scaled".format(scanName, version)
    if "BR" in combi:
        # outputDir = "limitCalculations/{}_{}_BRlimits_scaled".format(scanNameOld, version)
        outputDir = "limitCalculations/TChiNG_v11_BRlimits_scaled"

    outputDir = "limitCalculations/{}_{}".format(scanName, version)
    print outputDir
    # outputDir = "limitCalculations/{}_{}_higher".format(scanName, version)
    if not os.path.isdir(outputDir):
        os.mkdir(outputDir)
    xsecFile = getXsecFile(inputSignal)
    print xsecFile
    # writeDataCards(outputDir, inputData, inputSignal, combi, xsecFile)
    # drawLimitInput(outputDir, scanName, xsecFile)
    # drawSignalUncertainties(outputDir, scanName, xsecFile)
    # drawSignalUncertainties1D(outputDir, scanName, xsecFile)
    # callMultiCombine(outputDir)
    # callMultiSignificance(outputDir)
    # clearWrongCombineOutputs(outputDir)
    # if "TChi" in inputSignal:
    #     return proceedWithWeakScan(outputDir, scanName, xsecFile)
    build2dGraphs(outputDir, xsecFile)
    build1dGraphs(outputDir, xsecFile, scanName)
    # drawSignificance(outputDir, scanName)
    writeSMSLimitConfig(outputDir + "/saved_graphs1d_limit.root",
                        "../smsPlotter/config/SUS16047/%s_SUS16047.cfg" % scanName)
    writeSMSLimitConfig("/.automount/home/home__home4/institut_1b/swuchterl/SUSYDileptonGamma/plotter/" +
                        outputDir + "/saved_graphs1d_limit.root", "../smsPlotter/config/SUS16047/%s_SUS16047.cfg" % scanName)
    subprocess.call(["python2", "../smsPlotter/python/makeSMSplots.py",
                     "../smsPlotter/config/SUS16047/%s_SUS16047.cfg" % scanName, "plots/%s_limits_" % scanName])


if __name__ == "__main__":
    # signalScan("Zg", "v11", "limitCalculations/testDatacard.txt", "../minimal/output/SMS-T5bbbbZg_hists.root")
    # signalScan("Zg", "v11", "limitCalculations/testDatacard.txt", "../minimal/output_noVeto/SMS-T5bbbbZg_hists.root")
    # signalScan("Zg", "v12", "limitCalculations/testDatacard.txt", "../minimal/output_noVeto/SMS-T5bbbbZg_hists.root")
    # signalScan("Ng", "v11", "limitCalculations/testDatacard.txt", "../minimal/output/SMS-TChiNG_BF50N50G_hists.root")
    # signalScan("GGM", "v11", "limitCalculations/testDatacard.txt", "../minimal/output/GGM_GravitinoLSP_M1-200to1500_M2-200to1500_hists.root")
    # signalScan("GGM", "v13", "limitCalculations/testDatacard.txt", "../minimal/output/GGM_GravitinoLSP_M1-200to1500_M2-200to1500_hists.root")
    # signalScan("GMSB", "v13", "limitCalculations/testDatacard.txt", "../minimal/output/GMSB_GravitinoLSP_N1decays_hists.root")
    # signalScan("GMSB", "v13_scaleDoubleLumi", "limitCalculations/testDatacard.txt", "../minimal/output/GMSB_GravitinoLSP_N1decays_hists.root")
    # signalScan("Ng_BR", "v11", "limitCalculations/testDatacard.txt", "../minimal/output/SMS-TChiNG_BF50N50G_hists.root")
    # signalScan("Ng", "v11", "limitCalculations/testDatacard.txt", "../minimal/output/SMS-TChiNG_BF50N50G_hists.root")
    # signalScan("Ng", "v12", "limitCalculations/testDatacard.txt", "../minimal/output_noVeto/SMS-TChiNG_BF50N50G_hists.root")
    # signalScan("Ng", "v13_scaleDoubleLumi", "limitCalculations/testDatacard.txt", "../minimal/output/SMS-TChiNG_BF50N50G_hists.root")
    # signalScan("Zg", "v13", "limitCalculations/testDatacard.txt", "../minimal/output/SMS-T5bbbbZg_hists.root")
    # signalScan("Zg", "v13_normal16", "limitCalculations/testDatacard.txt", "../minimal/output/SMS-T5bbbbZg_hists.root")
    # first final try
    # signalScan("Ng", "v14", "limitCalculations/testDatacard.txt", "../myAnalyzer/output_TopPt_noNIsr/SMS-TChiNG_BF50N50G_hists.root")
    # signalScan("Ng", "v14", "limitCalculations/testDatacard.txt", "../myAnalyzer/output/SMS-TChiNG_BF50N50G_hists.root")
    # signalScan("Ng", "v15", "limitCalculations/testDatacard.txt", "../myAnalyzer/output_2/SMS-TChiNG_BF50N50G_hists.root")
    # signalScan("Ng", "v16", "limitCalculations/testDatacard.txt",
    #            "../myAnalyzer/output_2/SMS-TChiNG_BF50N50G_hists.root")
    # signalScan("Ng", "v17", "limitCalculations/testDatacard.txt",
    #            "../myAnalyzer/output_AN/SMS-TChiNG_BF50N50G_hists.root")
    # signalScan("Ng", "v18", "limitCalculations/testDatacard.txt",
    #            "../myAnalyzer/output_AN/SMS-TChiNG_BF50N50G_hists.root")
    # signalScan("GMSB", "v14", "limitCalculations/testDatacard.txt", "../myAnalyzer/output_TopPt_noNIsr/GMSB_GravitinoLSP_N1decays_hists.root")
    # signalScan("GMSB", "v15", "limitCalculations/testDatacard.txt", "../myAnalyzer/output_2/GMSB_GravitinoLSP_N1decays_hists.root")
    # signalScan("GMSB", "v16", "limitCalculations/testDatacard.txt",
    #            "../myAnalyzer/output_2/GMSB_GravitinoLSP_N1decays_hists.root")
    # signalScan("GMSB", "v17", "limitCalculations/testDatacard.txt",
    #            "../myAnalyzer/output_AN/GMSB_GravitinoLSP_N1decays_hists.root")
    # signalScan("GMSB", "v18", "limitCalculations/testDatacard.txt",
    #            "../myAnalyzer/output_AN/GMSB_GravitinoLSP_N1decays_hists.root")
    # signalScan("Zg", "v14", "limitCalculations/testDatacard.txt", "../myAnalyzer/output_TopPt_noNIsr/SMS-T5bbbbZg_hists.root")
    # signalScan("Zg", "v14", "limitCalculations/testDatacard.txt", "../myAnalyzer/output/SMS-T5bbbbZg_hists.root")
    # signalScan("Zg", "v16", "limitCalculations/testDatacard.txt",
    #            "../myAnalyzer/output_2/SMS-T5bbbbZg_hists.root")
    # signalScan("Zg", "v17", "limitCalculations/testDatacard.txt",
    #            "../myAnalyzer/output_AN/SMS-T5bbbbZg_hists.root")
    signalScan("Zg", "v18", "limitCalculations/testDatacard.txt",
               "../myAnalyzer/output_AN/SMS-T5bbbbZg_hists.root")
    # signalScan("Zg", "v17", "limitCalculations/testDatacard.txt",
    #            "../myAnalyzer/output_AN/SMS-T6ttZg_hists.root")
    # signalScan("GGM", "v16", "limitCalculations/testDatacard.txt", "../myAnalyzer/output_2/GGM_GravitinoLSP_M1-200to1500_M2-200to1500_hists.root")

    # signalScan("gg", "v11", "limitCalculations/observation_v7.txt", "../histogramProducer/SMS-T6Wg_signalScan_combined.root")
    # signalScan("Wg", "v11", "limitCalculations/observation_v7.txt", "../histogramProducer/SMS-T6Wg_signalScan_combined.root")
    # signalScan("gg", "v11", "limitCalculations/observation_v7.txt", "../histogramProducer/SMS-T5Wg_signalScan_combined.root")
    # signalScan("Wg", "v11", "limitCalculations/observation_v7.txt", "../histogramProducer/SMS-T5Wg_signalScan_combined.root")

    # signalScan("Wg", "v10", "limitCalculations/observation_v7.txt", "../histogramProducer/SMS-T5Wg_signalScan_combined.root")
    # signalScan("gg", "v10", "limitCalculations/observation_v7.txt", "../histogramProducer/SMS-T5Wg_signalScan_combined.root")
    # signalScan("Wg", "v10", "limitCalculations/observation_v7.txt", "../histogramProducer/SMS-T6Wg_signalScan_combined.root")
    # signalScan("gg", "v10", "limitCalculations/observation_v7.txt", "../histogramProducer/SMS-T6Wg_signalScan_combined.root")
    # signalScan("Wg", "v10", "limitCalculations/observation_v7.txt", "../histogramProducer/SMS-TChiNG_signalScan.root")
    # signalScan("Wg", "v10", "limitCalculations/observation_v7.txt", "../histogramProducer/SMS-TChiWG_signalScan.root")
