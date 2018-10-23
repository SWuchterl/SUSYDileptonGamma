import ROOT
from ROOT import *
from array import array
from include import *
import numpy as np


binnings = {
    'cutflow':   np.arange(0, 10, 10),
    # 'cutflow_fine':   np.arange(0, 5, 5)
    'cutflow_fine':   np.arange(0, 7, 7)
}
labels = {
    'cutflow': ["p_{T}^{leading}[GeV]", "Events / 30 GeV"],
    'cutflow_fine': ["p_{T}^{leading}[GeV]", "Events / 30 GeV"]
}


def drawSameHistogram(sampleNames, name, bkg=[], additional=[], binning=None, binningName="", scaleToData=False, xTitle=None, yTitle=None):
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()

    #style.divideByBinWidth = True

    style.minimumOne = False

    yTitle = None

    scale = 1.
    if scaleToData:
        scale = divideDatasetIntegrals(
            [i for i in additional if "Data" in i.label], bkg, name)

    for d in bkg[-1::-1]:
        h = d.getHist(name)
        if not h:
            continue
        if not h.Integral():
            continue
        h.Scale(scale)
        # if (binning.any()):
        if (binning):
            h = aux.rebin(h, binning)

        aux.appendFlowBin(h)
        #h.SetYTitle( aux.getYAxisTitle( h ) )
        if yTitle:
            h.SetYTitle(yTitle)
        else:
            h.SetYTitle(aux.getYAxisTitle(h))
        # if xTitle:
            #h.SetXTitle( xTitle )
        # h.GetXaxis().SetRange(0,10)
        h.LabelsDeflate("X")

        m.addStack(h, d.label)

    dataHist = None
    for d in additional:
        h = d.getHist(name)
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
            h.SetMarkerSize(0.5)
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

        # h.GetXaxis().SetRange(0,10)
        h.LabelsDeflate("X")

        m.add(h, d.label)

    # m.sortStackByIntegral()

    if m.Draw():

        # ratio
        if bkg != []:
            hsm = m.hists[0].GetStack().Last()
        if dataHist:
            r = ratio.Ratio("Data/MC", dataHist, hsm)
            r.draw(0.5, 1.5)

        #info = ""
        info = "ee" if "EE" in sampleNames else "#mu#mu" if "MM" in sampleNames else ""
        l = aux.Label(info="#scale[0.7]{%s}" % info, sim=((dataDoubleMuon not in additional)or(
            dataDoubleEG not in additional)or(dataHt not in additional)or(dataSF not in additional)))

        m.leg.SetNColumns(2)

        if binningName:
            binningName = "_" + binningName
        name = name.replace("/", "__")
        saveName = "sameHistograms_{}_{}{}".format(
            sampleNames, name, binningName)
        aux.save("DataMC_" + saveName,
                 folder="plots_CutFlow/", changeMinMax=True)


def main():
    # bkgs=[DYjetsLO,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma]
    # bkgs=[DYjetsNLO,zgamma,tt,ttgamma,wwgamma,wzgamma,zz,wjets,wgamma]
    bkgs = [zgamma, DYjetsNLO, tt, ttgamma, wwgamma,
            wzgamma, zz, wjets, wgamma, wjets, zz4l, ww, wz]
    # variables=["cutflow_fine","cutflow"]
    # variables = ["cutflow_fine"]
    variables = ["cutflow"]
    # groups=["cutFlow_Fine_onZ","cutFlow_onZ"]
    # groups = ["cutFlow_Fine_onZ"]
    groups = ["cutFlow_onZ"]
    for group in groups:
        for variable in variables:
            drawSameHistogram("EE+Data", group + "EE/" +
                              variable, bkgs, additional=[dataDoubleEG])
            drawSameHistogram("MM+Data", group + "MM/" +
                              variable, bkgs, additional=[dataDoubleMuon])
            drawSameHistogram("EM+Data", group + "EM/" +
                              variable, bkgs, additional=[dataMuonEG])
            #drawSameHistogram("EE",group+"EE/"+variable, bkgs, additional=[t5bbbbzg_1800_1700,t5bbbbzg_1800_400,t5bbbbzg_1800_600,tching_1200,tching_400],binning=binnings[variable])
            #drawSameHistogram("MM",group+"MM/"+variable, bkgs, additional=[t5bbbbzg_1800_1700,t5bbbbzg_1800_400,t5bbbbzg_1800_600,tching_1200,tching_400],binning=binnings[variable])
            #drawSameHistogram("EE",group+"EE/"+variable, bkgs, additional=[t5bbbbzg_1800_1700,t5bbbbzg_1800_400,t5bbbbzg_1800_600,tching_1200,tching_400])
            #drawSameHistogram("MM",group+"MM/"+variable, bkgs, additional=[t5bbbbzg_1800_1700,t5bbbbzg_1800_400,t5bbbbzg_1800_600,tching_1200,tching_400])
            # drawSameHistogram("EE", group + "EE/" +
            #                   variable, bkgs, additional=[])
            # drawSameHistogram("MM", group + "MM/" +
            #                   variable, bkgs, additional=[])
            # drawSameHistogram("EM", group + "EM/" +
            #                   variable, bkgs, additional=[])
            # drawSameHistogram("LL", group + "LL/" + variable,
            #                   bkgs, additional=[dataLL])
            drawSameHistogram("EEsignal", group + "EE/" + variable, bkg=[], additional=[
                              t5bbbbzg_1500_1400, t5bbbbzg_1500_400, t5bbbbzg_1500_600, tching_1200, tching_400])
            drawSameHistogram("MMsignal", group + "MM/" + variable, bkg=[], additional=[
                              t5bbbbzg_1500_1400, t5bbbbzg_1500_400, t5bbbbzg_1500_600, tching_1200, tching_400])
            drawSameHistogram("EMsignal", group + "EM/" + variable, bkg=[], additional=[
                              t5bbbbzg_1500_1400, t5bbbbzg_1500_400, t5bbbbzg_1500_600, tching_1200, tching_400])
            #drawSameHistogram("EEsignal",group+"EE/"+variable, bkg=[], additional=[t5bbbbzg_1800_1700,t5bbbbzg_1800_400,t5bbbbzg_1800_600,tching_1200,tching_400,ggm_m1550_m2750,ggm_m11200_m21000,ggm_m1800_m2600])
            #drawSameHistogram("MMsignal",group+"MM/"+variable, bkg=[], additional=[t5bbbbzg_1800_1700,t5bbbbzg_1800_400,t5bbbbzg_1800_600,tching_1200,tching_400,ggm_m1550_m2750,ggm_m11200_m21000,ggm_m1800_m2600])


if __name__ == "__main__":
    main()
