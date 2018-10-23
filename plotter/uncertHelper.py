from include import *


def scaleUncertainty(dataset, dirName, nBins, saveName=""):
    hNominal = aux.stdHist(dataset, dirName + "/met", nBins)
    hList = []
    for iWeight in range(9):
        h = aux.stdHist(dataset, dirName +
                        "/met_weight_{}".format(iWeight), nBins)
        hList.append(h.Clone(aux.randomName()))
    hUp, hDn = aux.getEnvelopeHists(hList)
    hSys = aux.getSystFromEnvelopes(hNominal, hUp, hDn)

    if saveName:
        c = ROOT.TCanvas()
        m = multiplot.Multiplot()
        hNominal.SetLineColor(ROOT.kBlack)
        m.add(hNominal)
        scales = [(1, 1), (1, 2), (1, .5), (2, 1), (2, 2),
                  (2, .5), (.5, 1), (.5, 2), (.5, .5)]
        for iWeight, (mur, muf) in enumerate(scales):
            h = hList[iWeight]
            h.drawOption_ = "hist"
            if iWeight / 3 == 0:
                h.SetLineColor(ROOT.kBlack)
            if iWeight / 3 == 1:
                h.SetLineColor(ROOT.kRed)
            if iWeight / 3 == 2:
                h.SetLineColor(ROOT.kBlue)
            if iWeight % 3 == 1:
                h.SetLineStyle(2)
            if iWeight % 3 == 2:
                h.SetLineStyle(3)
            m.add(h, "#mu_{{r}}={} #mu_{{f}}={}".format(mur, muf))
        m.Draw()
        denominator = hNominal.Clone(aux.randomName())
        for b in aux.loopH(denominator):
            denominator.SetBinError(b, 0)
        r = ratio.Ratio("scale uncert.", hNominal, denominator, hSys)
        r.draw(.5, 1.5)
        l = aux.Label(sim=True, info=dataset.label)
        aux.save("scaleUncertainty_" + saveName, normal=False)
    return hSys


def getScaleUncertHisto(nomHisto, binUncs, meanWeight=1.):
    nom = nomHisto.Clone()
    for bin in range(nom.GetNbinsX() + 2):
        #print "binnr",bin,nom.GetBinLowEdge(bin),binUncs
        c = nom.GetBinContent(bin)
        if c > 1e-10:
            if (bin >= len(binUncs)):
                e = 0.
            else:
                e = binUncs[bin] * c
            nom.SetBinError(bin, e)
        elif nom.GetBinContent(bin - 1) or nom.GetBinContent(bin + 1):
            # check if option "width" should be used
            # TODO: check if the weight agrees with the lumi+pu weight
            #meanWeight = nom.Integral(0,-1)/nom.GetEntries()
            poissonZeroError = 1.8410216450098775
            e = meanWeight * poissonZeroError
            e /= nom.GetBinWidth(bin) if style.divideByBinWidth else 1.
            nom.SetBinError(bin, e)
    return nom


def getPDFUncertHisto(nomHisto, binUncs, meanWeight=1.):
    nom = nomHisto.Clone()
    for bin in range(nom.GetNbinsX() + 2):
        #print "binnr",bin,nom.GetBinLowEdge(bin)
        c = nom.GetBinContent(bin)
        if c > 1e-10:
            if (bin >= len(binUncs)):
                e = 0.
            else:
                e = binUncs[bin] * c
            nom.SetBinError(bin, e)
        elif nom.GetBinContent(bin - 1) or nom.GetBinContent(bin + 1):
            # check if option "width" should be used
            # TODO: check if the weight agrees with the lumi+pu weight
            #meanWeight = nom.Integral(0,-1)/nom.GetEntries()
            poissonZeroError = 1.8410216450098775
            e = meanWeight * poissonZeroError
            e /= nom.GetBinWidth(bin) if style.divideByBinWidth else 1.
            nom.SetBinError(bin, e)
    return nom
