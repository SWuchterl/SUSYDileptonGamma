#!/usr/bin/env python2

import re
import subprocess
import ROOT
import optparse
import os
from Datacard import Datacard


def infoFromOut(out):
    infos = {"obs": 0, "exp": 0, "exp1up": 0,
             "exp1dn": 0, "exp2up": 0, "exp2dn": 0}
    for line in out.split("\n"):
        if line.startswith("NLL at global minimum of asimov:"):
            infos["rMinNLL"] = float(
                re.match(".*\(r = (.*)\).*", line).group(1))
        if line.startswith("Observed Limit: r < "):
            infos["obs"] = float(line.split("<")[1])
        if line.startswith("Expected  2.5%: r < "):
            infos["exp2dn"] = float(line.split("<")[1])
        if line.startswith("Expected 16.0%: r < "):
            infos["exp1dn"] = float(line.split("<")[1])
        if line.startswith("Expected 50.0%: r < "):
            infos["exp"] = float(line.split("<")[1])
        if line.startswith("Expected 84.0%: r < "):
            infos["exp1up"] = float(line.split("<")[1])
        if line.startswith("Expected 97.5%: r < "):
            infos["exp2up"] = float(line.split("<")[1])
    # if "rMinNLL" not in infos or infos["rMinNLL"] > 1.99999:
        #infos = { "obs":0, "exp":0, "exp1up":0, "exp1dn":0, "exp2up":0, "exp2dn":0 }
    return infos


def significanceFromOut(out):
    for line in out.split("\n"):
        if line.startswith("Significance:"):
            return float(line.split()[1])
    return 0


def guessSignalPoint(name):
    m = re.match(".*_(\d+)_(\d+).*", name)
    if m:
        return int(m.group(1)), int(m.group(2))
    else:
        print "could not determine signal point for", name
        return 0, 0


def guessScanName(name):
    short = "unknown"
    if "T5Wg" in name:
        short = "T5Wg"
    if "T5gg" in name:
        short = "T5gg"
    return short


def infosFromDatacard(name):
    return infoFromOut(callCombine(name))


def callSaveCombine(name, method, ending, options=[]):
    nameLimit = name + ending
    if not os.path.isfile(nameLimit) or os.path.getmtime(name) > os.path.getmtime(nameLimit):
        with open(nameLimit, "w+") as f:
            bn = os.path.basename(name)
            out = subprocess.check_output(
                ["combine", "-M", method, name, "-n", name] + options, stderr=subprocess.STDOUT)
            f.write(out)
            outputFile = "higgsCombine{}.{}.mH120.root".format(bn, method)
            if os.path.isfile(outputFile):
                os.remove(outputFile)
    with open(nameLimit) as f:
        out = f.read()
    return out


def callCombineHybrid(name):
    callSaveCombine(name, "HybridNew", ".limitHybrid", [
                    "--frequentist", "--testStat", "LHC", "-H", "ProfileLikelihood"])


def callCombine(name):
    nameLimit = name + ".limit"
    if not os.path.isfile(nameLimit) or os.path.getmtime(name) > os.path.getmtime(nameLimit):
        with open(nameLimit, "w+") as f:
            bn = os.path.basename(name)
            #print "calculating limit for", name
            out = subprocess.check_output(
                ["combine", "-M", "Asymptotic", name, "-n", bn], stderr=subprocess.STDOUT)
            f.write(out)
            outputFile = "higgsCombine{}.Asymptotic.mH120.root".format(bn)
            if os.path.isfile(outputFile):
                os.remove(outputFile)
    with open(nameLimit) as f:
        out = f.read()
    return out


def callCombineSignificance(name):
    nameLimit = name + ".significance"
    if not os.path.isfile(nameLimit) or os.path.getmtime(name) > os.path.getmtime(nameLimit):
        with open(nameLimit, "w+") as f:
            out = subprocess.check_output(["combine", "-M", "ProfileLikelihood", "--significance",
                                           "--uncapped", "1", "--rMin", "-5", name], stderr=subprocess.STDOUT)
            f.write(out)
    with open(nameLimit) as f:
        out = f.read()
    return out


def writeDictHere(d, filename):
    f = ROOT.TFile(filename, "recreate")
    for name, ob in d.iteritems():
        if ob:
            ob.Write(name)
    f.Close()


def getContour(gr2d):
    c = ROOT.TCanvas()
    # gr2d.Smooth()
    # gr2d.GetZaxis().SetRangeUser(0.5, 1.05)
    gr2d.SetMaximum(1.05)
    # gr2d.SetMinimum(0.5)
    gr2d.Draw("tri1")
    ROOT.gPad.Update()
    c.Update()
    contours = gr2d.GetContourList(1.)
    if not contours:
        print "Could not find contour"
        contours = [ROOT.TGraph()]
    contoursN = [(c, c.GetN()) for c in contours]

    contoursN = sorted(contoursN, key=lambda x: x[1])
    fil = ROOT.TFile("testCONT.root", "RECREATE")
    i = 0
    for cont, nc in contoursN:
        cont.Write("cont" + str(nc) + "_" + gr2d.GetName() +
                   "_" + str(i) + "_" + str(gr2d))
        i = i + 1
    print "_---"
    print gr2d.GetName()
    print contoursN
    if ("exp1up" in gr2d.GetName()):
        return contoursN[-2][0]
    # if ("exp1dn" in gr2d.GetName()):
    #     return contoursN[-2][0]
    return contoursN[-1][0]


def getContourTry(gr2d):
    c = ROOT.TCanvas()
    # gr2d.GetZaxis().SetRangeUser(0.5, 1.05)

    gr2d.SetNpx(150)
    gr2d.SetNpy(150)

    gr2d.SetMaximum(1.05)
    # gr2d.SetMinimum(0.5)
    gr2d.Draw("tri1")
    ROOT.gPad.Update()
    c.Update()
    contours = gr2d.GetContourList(1.)
    if not contours:
        print "Could not find contour"
        contours = [ROOT.TGraph()]
    contoursN = [(c, c.GetN()) for c in contours]

    contoursN = sorted(contoursN, key=lambda x: x[1])
    fil = ROOT.TFile("testCONT.root", "RECREATE")
    i = 0
    for cont, nc in contoursN:
        cont.Write("cont" + str(nc) + "_" + gr2d.GetName() +
                   "_" + str(i) + "_" + str(gr2d))
        i = i + 1
    print "_---"
    print gr2d.GetName()
    print contoursN
    if ("exp1up" in gr2d.GetName()):
        return contoursN[-2][0]
    return contoursN[-1][0]


def getContourAll(gr2d):
    from signalScan import smoothContour
    c = ROOT.TCanvas()
    gr2d.GetZaxis().SetRangeUser(0.5, 1.05)
    gr2d.SetMaximum(1.05)
    gr2d.SetMinimum(0.5)
    gr2d.Draw("tri1")
    ROOT.gPad.Update()
    contours = gr2d.GetContourList(1.)
    if not contours:
        print "Could not find contour"
        contours = [ROOT.TGraph()]
    contoursN = [(c, c.GetN()) for c in contours]

    contoursN = sorted(contoursN, key=lambda x: x[1])
    fil = ROOT.TFile(gr2d.GetName() + "_testCONT.root", "RECREATE")
    for cont, nc in contoursN:
        smoothContour(cont).Write("contSMOOTH" + str(nc) + gr2d.GetName())
        cont.Write("cont" + str(nc) + gr2d.GetName())
    # return contoursN[-1][0]
    # coll = ROOT.TCollection()
    coll = ROOT.TObjArray()
    # coll = []

    for cont, nc in contoursN:
        # coll.Add(cont)
        coll.Add(smoothContour(cont))
        # coll.append(cont)
    conFinal = ROOT.TGraph()
    conFinal.Merge(coll)
    # return contoursN[0:-1][0]
    return conFinal


def getFineHisto(gr):
    gr.SetMaximum(1.05)
    hsInter = ROOT.TH2F("", "", 140, gr.GetXaxis().GetXmin(), gr.GetXaxis(
    ).GetXmax(), 140, gr.GetYaxis().GetXmin(), gr.GetYaxis().GetXmax())
    for i in range(1, hsInter.GetXaxis().GetNbins() + 1):
        for j in range(1, hsInter.GetYaxis().GetNbins() + 1):
            # hsInter.SetBinContent(i, j, gr.Eval(
            hsInter.SetBinContent(i, j, gr.Interpolate(
                hsInter.GetXaxis().GetBinCenter(i), hsInter.GetYaxis().GetBinCenter(j)))
    # return getContourHS(hsInter),hsInter
    hsInter.SetDirectory(0)
    f = ROOT.TFile("dummyFineHisto_" + str(gr.GetName()) + ".root", "RECREATE")
    hsInter.Write()
    return hsInter


def getFineHisto2(gr):
    xIn = []
    yIn = []
    zIn = []
    for a in range(gr.GetN()):
        xIn.append(gr.GetX()[a])
        yIn.append(gr.GetY()[a])
        zIn.append(gr.GetZ()[a])
    # print xIn
    # print yIn
    # print zIn
    from scipy import interpolate
    import matplotlib.pyplot as plt
    import numpy as np
    # X, Y = np.meshgrid(xIn, yIn)
    # plt.pcolormesh(X, Y, zIn)
    # plt.show()

    # f = interpolate.interp2d(xIn, yIn, zIn, kind='linear')
    f = interpolate.Rbf(xIn, yIn, zIn, function='multiquadric')
    # print np.min(xIn)

    xOut = np.arange(np.min(xIn), np.max(
        xIn), (np.max(xIn) - np.min(xIn)) / 140)
    yOut = np.arange(np.min(yIn), np.max(
        yIn), (np.max(yIn) - np.min(yIn)) / 140)
    zOut = f(xOut, yOut)
    # print xOut
    # print yOut
    # print zOut
    grOut = ROOT.TGraph2D()
    k = 0
    for i in range(len(xOut)):
        for j in range(len(yOut)):
            # grOut.SetPoint(i, xOut[i], yOut[j], zOut[i])
            if(yOut[j] >= (xOut[i])):
                z = f(xOut[i], yOut[j])
                grOut.SetPoint(k, xOut[i], yOut[j], z)
                k = k + 1

    # hsInter.SetDirectory(0)
    # f = ROOT.TFile("dummyFineHisto_" + str(gr.GetName()) + ".root", "RECREATE")
    f = ROOT.TFile("dummyFineGraph_" +
                   str(grOut.GetName()) + ".root", "RECREATE")
    grOut.Write()
    # hsInter.Write()
    # return hsInter
    # hsOut = grOut.GetHistogram()
    # hsOut.SetDirectory(0)
    # hsOut.Write()
    # return getFineHisto(grOut)
    return grOut


def getShiftedHisto(gr):
    # gr2 = gr.Clone()
    gr2 = ROOT.TGraph2D()
    for a in range(gr.GetN()):
        # x = ROOT.Double(0)
        # y = ROOT.Double(0)
        # z = ROOT.Double(0)
        x = gr.GetX()[a]
        y = gr.GetY()[a]
        z = gr.GetZ()[a]
        # gr.GetPoint(a, x, y, z)
        gr2.SetPoint(a, x, y - x, z)
    # hsInter = ROOT.TH2F("", "", 200, gr.GetXaxis().GetXmin(), gr.GetXaxis(
    # ).GetXmax(), 200, gr.GetYaxis().GetXmin(), gr.GetYaxis().GetXmax())
    # for i in range(1, hsInter.GetXaxis().GetNbins() + 1):
    #     for j in range(1, hsInter.GetYaxis().GetNbins() + 1):
    #         hsInter.SetBinContent(i, j - i, gr.Interpolate(
    #             hsInter.GetXaxis().GetBinCenter(i), hsInter.GetYaxis().GetBinCenter(j)))
    hsInter = ROOT.TH2F("", "", 100, gr2.GetXaxis().GetXmin(), gr2.GetXaxis(
    ).GetXmax(), 100, gr2.GetYaxis().GetXmin(), gr2.GetYaxis().GetXmax())
    for i in range(1, hsInter.GetXaxis().GetNbins() + 1):
        for j in range(1, hsInter.GetYaxis().GetNbins() + 1):
            # hsInter.SetBinContent(i, j, gr.Eval(
            hsInter.SetBinContent(i, j - i, gr2.Interpolate(
                hsInter.GetXaxis().GetBinCenter(i), hsInter.GetYaxis().GetBinCenter(j)))
    # return getContourHS(hsInter),hsInter
    hsInter.SetDirectory(0)
    f = ROOT.TFile("dummyFineHistoShift_" +
                   str(gr.GetName()) + ".root", "RECREATE")
    hsInter.Write()
    # emptyH = hsInter.Clone()
    print "sadsad"
    graphClone = ROOT.TGraph2D(hsInter)
    print "asdsada"
    graphClone.Write()
    for i in range(graphClone.GetN()):
        x = graphClone.GetX()[a]
        y = graphClone.GetY()[a]
        z = graphClone.GetZ()[a]
        graphClone.SetPoint(i, x, y + x, z)
    print "asdsada"

    # for nBinx in range(emptyH.GetNbinsX()):
    #     for nBiny in range(emptyH.GetNbinsY()):
    #         nBin = emptyH.GetBin(nBinx, nBiny, 0)
    #         emptyH.SetBinContent(nBin, 0)
    #
    # emptyH.Write()
    # graphClone.Write()
    # retHist = graphClone.GetHistogram()
    # print "asdsada"
    # retHist.Write()
    # return emptyH
    return retHist


def getContour2(h2d):
    #print h2d
    # hist=gr2d.GetHistogram()
    #hist=ROOT.TH2F("", "", 36, 200, 1100, 36, 200, 1100)
    hist = h2d.Clone()
    #print hist
    for i in xrange(hist.GetNbinsX() + 1):
        for j in xrange(hist.GetNbinsY() + 1):
            if hist.GetBinContent(i, j) < 1:
                hist.SetBinContent(i, j, 0.5)
            else:
                hist.SetBinContent(i, j, 5.)
    #print hist
    hist.SetDirectory(0)
    c = ROOT.TCanvas()
    c.cd()
    hist.SetContour(1)
    # hist.SetContourLevel(0, 1.)
    hist.SetContourLevel(0, 5.)
    hist.Draw("CONT Z LIST")
    c.Update()
    conts = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
    f = ROOT.TFile("dummyconts_" + str(hist.GetName()) + ".root", "RECREATE")
    # f = ROOT.TFile(hs2.GetName()"_dummyconts_" + str(hist.GetName()) + ".root", "RECREATE")
    for c2 in conts:
        c2.Write(c2.GetName())
        # return c2
    #print conts.At(0).At(0)
    # conts.SetDirectory(0)
    if"exp1dn" in h2d.GetName():
        return f.Get("TList;1")
    if"exp1up" in h2d.GetName():
        return f.Get("TList;1")
    else:
        return f.Get("TList;1")

    #c  = ROOT.TCanvas()
    # gr2d.GetZaxis().SetRangeUser(0.5,1.05)
    # gr2d.SetMaximum(1.05)
    # gr2d.SetMinimum(0.5)
    # gr2d.Draw("tri1")
    # ROOT.gPad.Update()
    #contours = gr2d.GetContourList(1.)
    # if not contours:
        #print "Could not find contour"
        #contours = [ROOT.TGraph()]
    #contoursN = [(c,c.GetN()) for c in contours]

    #contoursN = sorted( contoursN, key=lambda x: x[1] )
    # fil=ROOT.TFile("testCONT.root","RECREATE")
    # for cont,nc in contoursN:
        # cont.Write("cont"+str(nc)+gr2d.GetName())
    # return contoursN[-1][0]


class MyDatacard(Datacard):
    def __init__(self, dc=""):
        if not dc:
            self.bins = []
            self.obs = {}
            #self.processes = ['signal', 'gqcd', 'ele', 'zg', 'wg', 'ttg']
            self.processes = []
            self.signals = ['signal']
            #self.isSignal = {'wg': False, 'signal': True, 'zg': False, 'gqcd': False, 'ele': False, 'ttg': False}
            self.isSignal = {}
            self.keyline = []
            self.exp = {}
            #self.systs = [(x, False, 'lnN', [], {}) for x in "lumi", "jec", "pdf", "gqcdSyst", "eleSyst", "wgSyst", "zgSyst", "ttgSyst"]
            self.systs = []
            self.shapeMap = {}
            self.hasShape = False
            self.flatParamNuisances = {}
            self.rateParams = {}
            self.rateParamsOrder = []
        elif isinstance(dc, Datacard):
            self.bins = dc.bins
            self.obs = dc.obs
            self.processes = dc.processes
            self.signals = dc.signals
            self.isSignal = dc.isSignal
            self.keyline = dc.keyline
            self.exp = dc.exp
            self.systs = dc.systs
            self.shapeMap = dc.shapeMap
            self.hasShape = dc.hasShape
            self.flatParamNuisances = dc.flatParamNuisances
            self.rateParams = dc.rateParams
            self.rateParamsOrder = dc.rateParamsOrder
        elif dc.endswith(".txt"):
            import DatacardParser
            options, b = DatacardParser.addDatacardParserOptions(
                optparse.OptionParser())
            mydc = MyDatacard(DatacardParser.parseCard(file(dc), options))
            self.__init__(mydc)
        else:
            print "Do not know how to initialize MyDatacard with", dc

    def __str__(self):
        print "Class MyDatacard"
        print "bins               ", self.bins
        print "obs                ", self.obs
        print "processes          ", self.processes
        print "signals            ", self.signals
        print "isSignal           ", self.isSignal
        print "keyline            ", self.keyline
        print "exp                ", self.exp
        for s in self.systs:
            print "systs              ", s
        #print "shapeMap           ", self.shapeMap
        #print "hasShape           ", self.hasShape
        #print "flatParamNuicances ", self.flatParamNuisances
        #print "rateParams         ", self.rateParams
        #print "rateParamsOrder    ", self.rateParamsOrder
        return ""

    def _getProcessNumbers(self):
        counter = 1
        processNumbers = {}
        for a, b in self.isSignal.iteritems():
            if b:
                processNumbers[a] = 0
            else:
                processNumbers[a] = counter
                counter += 1
        return processNumbers

    def write(self, filename=""):
        maxInfoLen = max([len(line[0]) + 7 for line in self.systs])
        out = ""
        out += "\nimax " + str(len(self.bins))
        out += "\njmax *"
        out += "\nkmax *"
        out += "\n\nbin         " + \
            ("{:>15}" * len(self.bins)).format(*self.bins)
        out += "\nobservation " + \
            ("{:>15}" * len(self.bins)).format(*
                                               [str(int(self.obs[x])) for x in self.bins])
        out += "\n\n"

        # create table for syst uncerts
        binNames, processNames, processNumbers = zip(*self.keyline)
        table = []
        table.append(["bin", ""] + list(binNames))
        table.append(["process", ""] + list(processNames))
        processNumbers = self._getProcessNumbers()
        table.append(["process", ""] + [str(processNumbers[x])
                                        for x in processNames])
        table.append(["rate", ""] + [str(round(self.exp[bN][processNames[i]], 3))
                                     for i, bN in enumerate(binNames)])
        for line in self.systs:
            relUncerts = [line[4][bN][processNames[i]]
                          for i, bN in enumerate(binNames)]
            table.append([line[0], line[2]] + ["-" if x == 1 or x ==
                                               0 else str(round(x, 3)) for x in relUncerts])
        # format lengts of strings
        columnWidths = [max([len(i) for i in line]) for line in zip(*table)]
        for irow, row in enumerate(table):
            for icol, col in enumerate(row):
                table[irow][icol] = "{{:>{}}}".format(
                    columnWidths[icol] + 1).format(col)
        # append table to output
        for row in table:
            out += ''.join(row) + "\n"

        if filename:
            with open(filename, "wb") as f:
                f.write(out)
                print "Writing to file:", filename
        else:
            print out

    def addBin(self, name, obs, bkgRates, bkgUncertainties):
        if self.processes:
            if self.processes != bkgRates.keys():
                print "ERROR: Old processes", self.processes, " New processes:", bkgRates.keys()
        else:
            self.processes = bkgRates.keys()
            self.isSignal = dict([(r, r == "signal") for r in self.processes])
        self.bins.append(name)
        self.obs[name] = obs
        # TODD check order, take ordered dict???
        self.keyline.extend([(name, process, process == "signal")
                             for process in bkgRates.keys()])
        self.exp[name] = bkgRates

        systDict = dict([(l[0], l) for l in self.systs])
        for source, newUncerts in bkgUncertainties.iteritems():
            for p in self.processes:
                if p not in newUncerts:
                    newUncerts[p] = 0
            if source not in systDict:
                systDict[source] = (source, False, "lnN", [], dict(
                    [(b, dict([(r, 0) for r in self.processes])) for b in self.bins]))
            systDict[source][4][name] = newUncerts
        for source, line in systDict.iteritems():
            for bin in self.bins:
                if bin not in line[4]:
                    systDict[source][4][bin] = dict(
                        [(r, 0) for r in self.processes])
        self.systs = sorted(systDict.values())

    def newSignal(self, exp, unc):
        #print self.exp
        for bName, newRate in exp.iteritems():
            self.exp[bName]["signal"] = newRate
        systDict = dict([(l[0], l) for l in self.systs])
        for uncName, uncertaintyDict in unc.iteritems():
            for binName, u in uncertaintyDict.iteritems():
                systDict[uncName][4][binName]["signal"] = u
        self.systs = sorted(systDict.values())

    def returnObs(self):
        return self.obs

    def limit(self):
        self.write("/tmp/tmpDataCard.txt")
        return infosFromDatacard("/tmp/tmpDataCard.txt")

    def limitFast(self):
        infos = {"obs": 0, "exp": 0, "exp1up": 0,
                 "exp1dn": 0, "exp2up": 0, "exp2dn": 0}

        for bin in self.bins:
            obs = self.obs[bin]
            bkg = sum(
                [b for a, b in self.exp[bin].iteritems() if a is not "signal"])
            signal = self.exp[bin]["signal"]
            if not signal:
                continue
            err = ROOT.TMath.Sqrt(bkg)
            r = signal / abs(err - abs(obs - bkg))
            infos["obs"] = max(infos["obs"], r)
            r_exp = signal / err
            if r_exp > infos["exp"]:
                infos["exp"] = r_exp
                infos["exp1up"] = (signal + err) / err
                infos["exp2up"] = (signal + 2 * err) / err
                infos["exp1dn"] = (signal - err) / err
                infos["exp2dn"] = (signal - 2 * err) / err
        return infos

    def setExpection(self):
        for binName, expDict in self.exp.iteritems():
            totRate = 0
            for process, rate in expDict.iteritems():
                totRate += 0 if process == "signal" else rate
            self.obs[binName] = round(totRate)


class Limit:
    def __init__(self, datacardname):
        self.datacardname = datacardname
        self.limitfilename = datacardname + ".limit"

        self.obs = -1
        self.exp = -1
        self.expUp = -1
        self.expDn = -1
        self.exp2Up = -1
        self.exp2Dn = -1
        self.rMinNLL = -1
        self.error = None

        self.dc = MyDatacard(datacardname)

    def calculate(self):
        if not os.path.isfile(self.limitfilename) or os.path.getmtime(self.datacardname) > os.path.getmtime(self.limitfilename):
            with open(self.limitfilename, "w+") as f:
                bn = os.path.basename(self.datacardname)
                out = subprocess.check_output(
                    ["combine", "-M", "Asymptotic", self.datacardname, "-n", bn], stderr=subprocess.STDOUT)
                f.write(out)
                outputFile = "higgsCombine{}.Asymptotic.mH120.root".format(bn)
                if os.path.isfile(outputFile):
                    os.remove(outputFile)
        with open(self.limitfilename) as f:
            out = f.read()

    def getInfo(self):
        if not os.path.isfile(self.limitfilename):
            self.calculate()
        with open(self.limitfilename) as f:
            out = f.read()
        lines = out.split("\n")
        for line in lines:
            if line.startswith("NLL at global minimum of asimov:"):
                self.rMinNLL = float(
                    re.match(".*\(r = (.*)\).*", line).group(1))
            if line.startswith("Observed Limit: r < "):
                self.obs = float(line.split("<")[1])
            if line.startswith("Expected  2.5%: r < "):
                self.exp2Dn = float(line.split("<")[1])
            if line.startswith("Expected 16.0%: r < "):
                self.expDn = float(line.split("<")[1])
            if line.startswith("Expected 50.0%: r < "):
                self.exp = float(line.split("<")[1])
            if line.startswith("Expected 84.0%: r < "):
                self.expUp = float(line.split("<")[1])
            if line.startswith("Expected 97.5%: r < "):
                self.exp2Up = float(line.split("<")[1])
            if "ERROR" in line:
                self.error = line


if False:
    inFileName = "limitCalculations/observation_v3.txt"
    dc = MyDatacard(inFileName)
    dc.setExpection()
    dc.write("test.txt")

    # dc.newSignal({
    #    "bin25": (12,1.1),
    #    "bin26": (1, 1.1),
    #    "bin27": (3, 1.2)
    # })

    #print dc
    # dc.write()
    #print dc.limit()
