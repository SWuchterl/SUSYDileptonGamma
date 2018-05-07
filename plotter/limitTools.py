#!/usr/bin/env python2

import re
import subprocess
import ROOT
import optparse
import os
from Datacard import Datacard

def infoFromOut(out):
    infos = { "obs":0, "exp":0, "exp1up":0, "exp1dn":0, "exp2up":0, "exp2dn":0 }
    for line in out.split("\n"):
        if line.startswith("NLL at global minimum of asimov:"): infos["rMinNLL"] = float(re.match(".*\(r = (.*)\).*",line).group(1))
        if line.startswith("Observed Limit: r < "): infos["obs"]   = float(line.split("<")[1])
        if line.startswith("Expected  2.5%: r < "): infos["exp2dn"] = float(line.split("<")[1])
        if line.startswith("Expected 16.0%: r < "): infos["exp1dn"] = float(line.split("<")[1])
        if line.startswith("Expected 50.0%: r < "): infos["exp"]    = float(line.split("<")[1])
        if line.startswith("Expected 84.0%: r < "): infos["exp1up"] = float(line.split("<")[1])
        if line.startswith("Expected 97.5%: r < "): infos["exp2up"] = float(line.split("<")[1])
    #if "rMinNLL" not in infos or infos["rMinNLL"] > 1.99999:
        #infos = { "obs":0, "exp":0, "exp1up":0, "exp1dn":0, "exp2up":0, "exp2dn":0 }
    return infos

def significanceFromOut(out):
    for line in out.split("\n"):
        if line.startswith("Significance:"):
            return float(line.split()[1])
    return 0

def guessSignalPoint(name):
    m = re.match(".*_(\d+)_(\d+).*",name)
    if m:
        return int(m.group(1)), int(m.group(2))
    else:
        print "could not determine signal point for", name
        return 0, 0

def guessScanName(name):
    short = "unknown"
    if "T5Wg" in name: short = "T5Wg"
    if "T5gg" in name: short = "T5gg"
    return short

def infosFromDatacard(name):
    return infoFromOut(callCombine(name))

def callSaveCombine(name, method, ending, options=[]):
    nameLimit = name + ending
    if not os.path.isfile(nameLimit) or os.path.getmtime(name)>os.path.getmtime(nameLimit):
        with open(nameLimit, "w+") as f:
            bn = os.path.basename(name)
            out = subprocess.check_output(["combine", "-M", method, name, "-n", name] + options, stderr=subprocess.STDOUT)
            f.write(out)
            outputFile = "higgsCombine{}.{}.mH120.root".format(bn, method)
            if os.path.isfile(outputFile): os.remove(outputFile)
    with open(nameLimit) as f:
        out = f.read()
    return out

def callCombineHybrid(name):
    callSaveCombine(name, "HybridNew", ".limitHybrid", ["--frequentist", "--testStat", "LHC", "-H", "ProfileLikelihood"])

def callCombine(name):
    nameLimit = name + ".limit"
    if not os.path.isfile(nameLimit) or os.path.getmtime(name)>os.path.getmtime(nameLimit):
        with open(nameLimit, "w+") as f:
            bn = os.path.basename(name)
            #print "calculating limit for", name
            out = subprocess.check_output(["combine", "-M", "Asymptotic", name, "-n", bn], stderr=subprocess.STDOUT)
            f.write(out)
            outputFile = "higgsCombine{}.Asymptotic.mH120.root".format(bn)
            if os.path.isfile(outputFile): os.remove(outputFile)
    with open(nameLimit) as f:
        out = f.read()
    return out

def callCombineSignificance(name):
    nameLimit = name + ".significance"
    if not os.path.isfile(nameLimit) or os.path.getmtime(name)>os.path.getmtime(nameLimit):
        with open(nameLimit, "w+") as f:
            out = subprocess.check_output(["combine", "-M", "ProfileLikelihood", "--significance", "--uncapped", "1", "--rMin", "-5", name], stderr=subprocess.STDOUT)
            f.write(out)
    with open(nameLimit) as f:
        out = f.read()
    return out

def getContour( gr2d ):
    c  = ROOT.TCanvas()
    gr2d.Draw()
    ROOT.gPad.Update()
    contours = gr2d.GetContourList(1.)
    if not contours:
        print "Could not find contour"
        contours = [ROOT.TGraph()]
    contoursN = [(c,c.GetN()) for c in contours]
    contoursN = sorted( contoursN, key=lambda x: x[1] )
    return contoursN[-1][0]

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
            options, b = DatacardParser.addDatacardParserOptions(optparse.OptionParser())
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
            if b: processNumbers[a] = 0
            else:
                processNumbers[a] = counter
                counter += 1
        return processNumbers


    def write(self, filename=""):
        maxInfoLen = max( [len(line[0])+7 for line in self.systs])
        out = ""
        out += "\nimax "+str(len(self.bins))
        out += "\njmax *"
        out += "\nkmax *"
        out += "\n\nbin         " + ("{:>15}"*len(self.bins)).format(*self.bins)
        out += "\nobservation " + ("{:>15}"*len(self.bins)).format(*[str(int(self.obs[x])) for x in self.bins])
        out += "\n\n"

        # create table for syst uncerts
        binNames, processNames, processNumbers = zip(*self.keyline)
        table = []
        table.append(["bin", ""] + list(binNames))
        table.append(["process", ""] + list(processNames))
        processNumbers = self._getProcessNumbers()
        table.append(["process", ""] + [str(processNumbers[x]) for x in processNames])
        table.append(["rate", ""] + [str(round(self.exp[bN][processNames[i]],3)) for i, bN in enumerate(binNames)])
        for line in self.systs:
            relUncerts = [line[4][bN][processNames[i]] for i, bN in enumerate(binNames)]
            table.append([line[0], line[2]] + ["-" if x==1 or x==0 else str(round(x,3)) for x in relUncerts])
        # format lengts of strings
        columnWidths = [max([len(i) for i in line]) for line in zip(*table)]
        for irow, row in enumerate(table):
            for icol, col in enumerate(row):
                table[irow][icol] = "{{:>{}}}".format(columnWidths[icol]+1).format(col)
        # append table to output
        for row in table: out += ''.join(row) + "\n"

        if filename:
            with open(filename, "wb") as f:
                f.write(out)
                print "Writing to file:", filename
        else:
            print out

    def addBin(self, name, obs, bkgRates, bkgUncertainties):
        if self.processes:
            if self.processes != bkgRates.keys(): print "ERROR: Old processes", self.processes, " New processes:", bkgRates.keys()
        else:
            self.processes = bkgRates.keys()
            self.isSignal = dict([(r, r=="signal") for r in self.processes])
        self.bins.append(name)
        self.obs[name] = obs
        self.keyline.extend([(name, process, process=="signal") for process in bkgRates.keys()]) # TODD check order, take ordered dict???
        self.exp[name] = bkgRates

        systDict = dict([(l[0],l) for l in self.systs])
        for source, newUncerts in bkgUncertainties.iteritems():
            for p in self.processes:
                if p not in newUncerts:
                    newUncerts[p] = 0
            if source not in systDict:
                systDict[source] = (source, False, "lnN", [], dict([(b,dict([(r,0) for r in self.processes])) for b in self.bins]))
            systDict[source][4][name] = newUncerts
        for source, line in systDict.iteritems():
            for bin in self.bins:
                if bin not in line[4]:
                    systDict[source][4][bin] = dict([(r,0) for r in self.processes])
        self.systs = sorted(systDict.values())

    def newSignal(self, exp, unc):
        #print self.exp
        for bName, newRate in exp.iteritems():
            self.exp[bName]["signal"] = newRate
        systDict = dict([(l[0],l) for l in self.systs])
        for uncName, uncertaintyDict in unc.iteritems():
            for binName, u in uncertaintyDict.iteritems():
                systDict[uncName][4][binName]["signal"] = u
        self.systs = sorted(systDict.values())


    def limit(self):
        self.write("/tmp/tmpDataCard.txt")
        return infosFromDatacard("/tmp/tmpDataCard.txt")

    def limitFast(self):
        infos = { "obs":0, "exp":0, "exp1up":0, "exp1dn":0, "exp2up":0, "exp2dn":0 }

        for bin in self.bins:
            obs = self.obs[bin]
            bkg = sum([b for a,b in self.exp[bin].iteritems() if a is not "signal"])
            signal = self.exp[bin]["signal"]
            if not signal: continue
            err = ROOT.TMath.Sqrt(bkg)
            r = signal/abs(err - abs(obs - bkg))
            infos["obs"] = max(infos["obs"], r)
            r_exp = signal/err
            if r_exp > infos["exp"]:
                infos["exp"] = r_exp
                infos["exp1up"] = (signal+err)/err
                infos["exp2up"] = (signal+2*err)/err
                infos["exp1dn"] = (signal-err)/err
                infos["exp2dn"] = (signal-2*err)/err
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
        self.limitfilename = datacardname+".limit"

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
        if not os.path.isfile(self.limitfilename) or os.path.getmtime(self.datacardname)>os.path.getmtime(self.limitfilename):
            with open(self.limitfilename, "w+") as f:
                bn = os.path.basename(self.datacardname)
                out = subprocess.check_output(["combine", "-M", "Asymptotic", self.datacardname, "-n", bn], stderr=subprocess.STDOUT)
                f.write(out)
                outputFile = "higgsCombine{}.Asymptotic.mH120.root".format(bn)
                if os.path.isfile(outputFile): os.remove(outputFile)
        with open(self.limitfilename) as f:
            out = f.read()

    def getInfo(self):
        if not os.path.isfile(self.limitfilename): self.calculate()
        with open(self.limitfilename) as f:
            out = f.read()
        lines = out.split("\n")
        for line in lines:
            if line.startswith("NLL at global minimum of asimov:"): self.rMinNLL = float(re.match(".*\(r = (.*)\).*",line).group(1))
            if line.startswith("Observed Limit: r < "): self.obs    = float(line.split("<")[1])
            if line.startswith("Expected  2.5%: r < "): self.exp2Dn = float(line.split("<")[1])
            if line.startswith("Expected 16.0%: r < "): self.expDn  = float(line.split("<")[1])
            if line.startswith("Expected 50.0%: r < "): self.exp    = float(line.split("<")[1])
            if line.startswith("Expected 84.0%: r < "): self.expUp  = float(line.split("<")[1])
            if line.startswith("Expected 97.5%: r < "): self.exp2Up = float(line.split("<")[1])
            if "ERROR" in line: self.error = line



if False:
    inFileName = "limitCalculations/observation_v3.txt"
    dc = MyDatacard(inFileName)
    dc.setExpection()
    dc.write("test.txt")

    #dc.newSignal({
    #    "bin25": (12,1.1),
    #    "bin26": (1, 1.1),
    #    "bin27": (3, 1.2)
    #})

    #print dc
    #dc.write()
    #print dc.limit()
