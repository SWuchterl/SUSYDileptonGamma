# -*- coding: utf-8 -*-
import array
import ConfigParser
import itertools
from math import *
import numpy
import pickle
import re
import ROOT
import style
#import multiplot

binCfg = ConfigParser.SafeConfigParser()
#binCfg.readfp(open('rebin.cfg'))

#ROOT.gErrorIgnoreLevel = kPrint, kInfo, kWarning, kError, kBreak, kSysError, kFatal
#ROOT.gErrorIgnoreLevel = ROOT.kBreak
saveStuff = []

def frange(x, y, jump):
    while x < y:
        yield x
        x += jump

def frange(x, y, jump):
    return [ x+jump*i for i in range(int(floor((y-x)/jump)))]

def getXsecInfoSMS( mother_mass, pklfilename ):
    info = (0,0)
    with open( pklfilename, 'rb') as f:
        data = pickle.load( f )
        if mother_mass in data:
            info = data[mother_mass]
        else:
            print "could not find mass {} in file {}".format(mother_mass, pklfilename)
    return info
def getXsecInfoGGM( mother_mass1, mother_mass2, pklfilename ):
    info = (0,0)
    with open( pklfilename, 'rb') as f:
        data = pickle.load( f )
        if mother_mass1 in data:
            info = data[mother_mass1][mother_mass2]
        else:
            print "could not find mass {} in file {}".format(mother_mass, pklfilename)
    return info

def getXsecSMSglu( mother_mass ):
    return getXsecInfoSMS( mother_mass, "data/xSec_SMS_Gluino_13TeV.pkl" )[0]
def getXsecSMSsq( mother_mass ):
    return getXsecInfoSMS( mother_mass, "data/xSec_SMS_Squark_13TeV.pkl" )[0]
def getXsecSMSstop( mother_mass ):
    return getXsecInfoSMS( mother_mass, "data/xSec_SMS_StopSbottom_13TeV.pkl" )[0]
def getXsecSMSneu( mother_mass ):
    return getXsecInfoSMS( mother_mass, "data/xSec_SMS_C1C1_13TeV.pkl" )[0]

def getXsecFromName( name ):
    m = re.match( "[^_]*_(\d+)_[^\d]*.*", name )
    if not m: m = re.match( ".*mGluino-(\d+)_.*", name )
    firstNumber = int( m.groups()[0] ) if m else 0

    if "T5" in name:
        return getXsecSMSglu( firstNumber )
    if "TChi" in name:
        return getXsecSMSneu( firstNumber )

    if "T2ttgg" in name:
        return getXsecSMSstop( firstNumber )

    print "cout not guess xsec from file name"
    return 0

def getSignalLabel( dataset ):
    m = re.match( "([^_]+)_(\d+)_(\d+)", dataset )
    if m:
        return "{} m(#tilde{{g}})={} m(#tilde{{#chi}}^{{0}}_{{1}})={}".format(*m.groups())
    return dataset

def getDatasetFromKey(key):
    m = re.match( "T([^_]+)_(\d+)_(\d+)", key )
    if m: return "SMS-T{}_mGluino-{}_mNeutralino-{}".format(*m.groups())
    return key


def getProjections( h2, axis="x", scale=True ):
    if axis == "x":
        a1 = h2.GetXaxis()
        a2 = h2.GetYaxis()
    elif axis == "y":
        a1 = h2.GetYaxis()
        a2 = h2.GetXaxis()
    else:
        print "please use 'x' or 'y' as argument"
    hs = []
    label = a2.GetTitle()

    for ybin in range( a2.GetNbins()+2 ):

        ylow = a2.GetBinLowEdge(ybin)
        yhigh = a2.GetBinUpEdge(ybin)
        name = "{} #leq {} < {}".format( ylow, label, yhigh )
        if ybin == 0: name = "{} < {}".format( label, yhigh )
        if ybin == a2.GetNbins()+1: name = "{} #leq {}".format( ylow, label )

        h = h2.ProjectionX( name, ybin, ybin ) if axis == "x" else h2.ProjectionY( name, ybin, ybin )
        h.GetXaxis().SetTitleOffset(h2.GetXaxis().GetTitleOffset())
        h.SetLineColor(getPaletteColor(1.*ybin/(a2.GetNbins()+2)))
        h.SetTitle(";{};{}".format(h.GetXaxis().GetTitle(),h2.GetZaxis().GetTitle()))
        if not h.GetEntries(): continue
        if scale: h.Scale( 1./h.GetEntries() )
        hs.append( h )

    return hs

def drawContributions(stack, mini=0, maxi=1.1, title="Fractions"):
    total = stack.GetStack().Last()
    b = total.GetNbinsX()-1
    newStack = ROOT.THStack()
    hists = [h.Clone(h.GetName()+randomName()) for h in stack.GetStack()]
    for ih in range(len(hists)-1, 0, -1):
        hists[ih].Add(hists[ih-1],-1)
    for h in hists:
        h.Divide(total)
        h.SetFillColor(h.GetFillColor())
        newStack.Add(h)
    newStack.SetTitle(";{};{}".format(h.GetXaxis().GetTitle(),title))
    newStack.SetMinimum(mini)
    newStack.SetMaximum(maxi)
    newStack.Draw("hist")
    saveStuff.append(newStack)

def write2File( obj2Write, name, fname ):
    obj = obj2Write.Clone()
    if isinstance( obj, ROOT.TH1 ):
        for a in obj.GetXaxis(),obj.GetYaxis():
            pass
            #a.UnZoom() #TODO: find out why this affects original histogram as well

    f = ROOT.TFile( fname,"update")
    obj.Write( name, ROOT.TObject.kWriteDelete )
    f.Close()

def writeWeight( obj, name, sampleName ):
    print "Info in <auxiliary::writeWeight>: {} written to file {}".format(name, "weights.root" )
    write2File( obj, name.replace("h_","weight_{}_".format(sampleName)), "weights.root" )

def getFromFile( filename, histoname ):
    f = ROOT.TFile( filename )
    h = f.Get( histoname )
    if not h:
        if ROOT.gErrorIgnoreLevel < ROOT.kBreak:
            print "Object {} not found in file {}".format(histoname, filename)
        return
    h = ROOT.gROOT.CloneObject( h )
    if isinstance( h, ROOT.TH1 ) and not h.GetSumw2N():
        h.Sumw2()
    h.drawOption_ = ""
    return h

def writeNumber(filename, name, number):
    f = ROOT.TFile(filename, "update")
    t = ROOT.TTree(name, "Title")
    v = numpy.zeros(1, dtype=float)
    t.Branch("var", v, "var/D")
    v[0] = number
    t.Fill()
    t.Write(name, ROOT.TObject.kWriteDelete)
    f.Close()

def readNumber(filename, name):
    t = ROOT.TChain(name)
    t.AddFile(filename)
    for e in t:
        return e.var
    print "Could not return {} from {}".format(name, filename)

def dataLikeName(name):
    return "Data" in name or "Pseudodata" in name or "Direct" in name

def getNgen( filename ):
    h = getFromFile( filename, "hCutFlow" )
    return h.GetBinContent(2)

def getObjectNames( filename, path="", objects=[ROOT.TH1] ):
    f = ROOT.TFile( filename )
    tmpDir = f.GetDirectory( path )

    outList = []
    for element in tmpDir.GetListOfKeys():
        if element.GetName() == "rawEff_vs_run": continue
        obj = element.ReadObj()

        if any( [ isinstance( obj, o ) for o in objects ] ):
            outList.append( element.GetName() )

    return outList

def getDirNames(filename):
    f = ROOT.TFile(filename)
    allDirs = [ e.GetName() for e in f.GetListOfKeys() if isinstance(e.ReadObj(),ROOT.TDirectory) ]
    return allDirs

def inerate(axis):
    return range(axis.GetNbins()+2)

def qa(l):
    # quadratic addition
    return math.sqrt(sum([i**2 for i in l]))

def getBinning( axis ):
    binning = []
    for i in range(axis.GetNbins()+1):
        binning.append( axis.GetBinUpEdge(i) )
    return binning

def checkRebinningConsistence( axis, newBinning ):
    oldBinning = getBinning( axis )
    # get rid of unprecise floats:
    oldBinning = [ round(i,5) for i in oldBinning ]
    newBinning = [ round(i,5) for i in newBinning ]
    # ignore new bin edges out of range of old binning
    newBinning = [ i for i in newBinning if i>= oldBinning[0] ]
    for i in newBinning:
        if i not in oldBinning: print "New bin edge is not compatible with old binning", i, "old binning:", oldBinning

def rebin2d( h, binEdgesX=None, binEdgesY=None ):
    if not binEdgesX and not binEdgesY: return h
    # Check consistency with old binning
    if binEdgesX:
        checkRebinningConsistence( h.GetXaxis(), binEdgesX )
    else:
        binEdgesX = [ h.GetXaxis().GetBinLowEdge(bin+1) for bin in range(h.GetNbinsX())]

    if binEdgesY:
        checkRebinningConsistence( h.GetYaxis(), binEdgesY )
    else:
        binEdgesY = [ h.GetYaxis().GetBinLowEdge(bin+1) for bin in range(h.GetNbinsY())]

    # Create
    binEdgesXArr = array.array( 'd', binEdgesX )
    binEdgesYArr = array.array( 'd', binEdgesY )
    hnew = ROOT.TH2F(h.GetName()+randomName(),h.GetTitle(), len(binEdgesX)-1, binEdgesXArr, len(binEdgesY)-1, binEdgesYArr )
    hnew.GetYaxis().SetTitleOffset(h.GetYaxis().GetTitleOffset())

    # GetProperties
    hnew.drawOption_ = h.drawOption_ if hasattr( h, "drawOption_" ) else ""
    hnew.SetTitle("{};{};{}".format(h.GetTitle(),h.GetXaxis().GetTitle(),h.GetYaxis().GetTitle()))

    # Fill
    for xbin in range(h.GetNbinsX()+2):
        x = h.GetXaxis().GetBinCenter(xbin)
        for ybin in range(h.GetNbinsY()+2):
            y = h.GetYaxis().GetBinCenter(ybin)
            newBin = hnew.FindFixBin(x,y)
            hnew.SetBinContent(newBin, hnew.GetBinContent(newBin)+h.GetBinContent(xbin,ybin))
            hnew.SetBinError( newBin, sqrt(hnew.GetBinError(newBin)**2 + h.GetBinError(xbin,ybin)**2 ) )
    return hnew

def rebin3d( h, binEdgesX=None, binEdgesY=None, binEdgesZ=None ):
    if not binEdgesX and not binEdgesY and not binEdgesZ: return h
    # Check consistency with old binning
    if binEdgesX:
        checkRebinningConsistence( h.GetXaxis(), binEdgesX )
    else:
        binEdgesX = [ h.GetXaxis().GetBinLowEdge(bin+1) for bin in range(h.GetNbinsX())]
    if binEdgesY:
        checkRebinningConsistence( h.GetYaxis(), binEdgesY )
    else:
        binEdgesY = [ h.GetYaxis().GetBinLowEdge(bin+1) for bin in range(h.GetNbinsY())]
    if binEdgesZ:
        checkRebinningConsistence( h.GetZaxis(), binEdgesZ )
    else:
        binEdgesZ = [ h.GetZaxis().GetBinLowEdge(bin+1) for bin in range(h.GetNbinsZ())]

    # Create
    binEdgesXArr = array.array( 'd', binEdgesX )
    binEdgesYArr = array.array( 'd', binEdgesY )
    binEdgesZArr = array.array( 'd', binEdgesZ )
    hnew = ROOT.TH3F(h.GetName(),h.GetTitle(), len(binEdgesX)-1, binEdgesXArr, len(binEdgesY)-1, binEdgesYArr, len(binEdgesZ)-1, binEdgesZArr )

    # GetProperties
    hnew.drawOption_ = h.drawOption_ if hasattr( h, "drawOption_" ) else ""
    hnew.SetTitle("{};{};{};{}".format(h.GetTitle(),h.GetXaxis().GetTitle(),h.GetYaxis().GetTitle(),h.GetZaxis().GetTitle()))

    # Fill
    for xbin in range(h.GetNbinsX()+2):
        x = h.GetXaxis().GetBinCenter(xbin)
        for ybin in range(h.GetNbinsY()+2):
            y = h.GetYaxis().GetBinCenter(ybin)
            for zbin in range(h.GetNbinsZ()+2):
                z = h.GetZaxis().GetBinCenter(zbin)
                newBin = hnew.FindFixBin(x,y,z)
                hnew.SetBinContent(newBin, hnew.GetBinContent(newBin) + h.GetBinContent(xbin,ybin,zbin))
                hnew.SetBinError(newBin, sqrt(hnew.GetBinError(newBin)**2 + h.GetBinError(xbin,ybin,zbin)**2 ) )
    return hnew



def rebin( h, binEdges, scale=True ):
    if not binEdges: return h
    #if not binEdges.any(): return h
    checkRebinningConsistence( h.GetXaxis(), binEdges )
    binEdgesArr = array.array( 'd', binEdges )
    #print binEdgesArr
    hnew = h.Rebin( len(binEdges)-1, "new", binEdgesArr )
    hnew.drawOption_ = h.drawOption_ if hasattr( h, "drawOption_" ) else ""
    #if scale: hnew.Scale( 1., "width" )
    if style.divideByBinWidth: hnew.Scale(1.,"width")
    return hnew

def rebinX(h, binEdgesX=None, binEdgesY=None, binEdgesZ=None):
    nDim = h.GetDimension()
    if nDim == 1:
        return rebin(h, binEdgesX, False)
    elif nDim == 2:
        return rebin2d(h, binEdgesX, binEdgesY)
    elif nDim == 3:
        return rebin3d(h, binEdgesX, binEdgesY, binEdgesZ)
    else:
        print "Do not know what to do with dimension", nDim

def absHistWeighted( origHist ):
    origNbins = origHist.GetNbinsX()
    origXmin = origHist.GetBinLowEdge(1)
    origXmax = origHist.GetBinLowEdge(origHist.GetNbinsX()+1)
    if origXmin + origXmax > 0.1:
        if origXmin:
            print "cant handle assymetric histograms"
        # else: print "already symetric?"
        return origHist

    h = origHist.Clone()
    newN = int(ceil(origNbins/2.))
    # TODO: doesn not work, only first bin in filled
    h = rebin( h, [ origXmax *i/newN for i in range(newN+1)] )

    for origBin in range( origNbins+2 ):
        newBin = int(abs(origBin - (origNbins+1.)/2)) + 1

        c1 = origHist.GetBinContent( origBin )
        e1 = origHist.GetBinError( origBin )
        c2 = h.GetBinContent( newBin )
        e2 = h.GetBinError( newBin )

        if e1 and e2:
            h.SetBinContent( newBin, ( c1*e1**-2 + c2*e2**-2 )/(e1**-2 + e2**-2) )
            h.SetBinError( newBin, 1./sqrt( e1**-2 + e2**-2 ) )

        else:
            h.SetBinContent( newBin, origHist.GetBinContent(origBin) )
            h.SetBinError( newBin, origHist.GetBinError(origBin) )

    return h

def randomName():
    """
    Generate a random string. This function is useful to give ROOT objects
    different names to avoid overwriting.
    """
    from random import randint
    from sys import maxint
    return "%x"%(randint(0, maxint))


def TH1F_binning( name, title, binEdges ):
    # Wrapper for the TH1F constructor for variable binning
    binEdgesArray = array.array( "d", binEdges )
    return ROOT.TH1F( name, ";%s;Events"%xVar.title, len(binEdgesArray)-1, binEdgesArray )

def createHistoFromTree( tree, variable, weight, nBins=20, xmin=0, xmax=0 ):
    name = randomName()
    if isinstance( nBins, list ):
        result = TH1F_binning( name, variable, nBins )
    else:
        result = ROOT.TH1F( name, variable, nBins, xmin, xmax )
    tree.Draw("%s>>%s"%(variable, name), weight, "goff")
    if style.divideByBinWidth: result.Scale(1, "width")
    return result

def sumSq( *items ):
    return sqrt( sum( [i**2 for i in items ] ) )

def mergeBins( h, dest, source ):
    h.SetBinContent( dest, h.GetBinContent( dest ) + h.GetBinContent( source ) )
    h.SetBinError( dest, sumSq( h.GetBinError( dest ), h.GetBinError( source ) ) )
    h.SetBinContent( source, 0 )
    h.SetBinError( source, 0 )

def appendFlowBin( h, under=True, over=True ):
    if under:
        mergeBins( h, 1, 0 )
    if over:
        mergeBins( h, h.GetNbinsX(), h.GetNbinsX()+1 )

def appendFlowBin3d(h, mergeX=True, mergeY=True, mergeZ=True):
    if mergeZ:
        for xb in range(h.GetNbinsX()+2):
            for yb in range(h.GetNbinsY()+2):
                mergeBins(h, h.GetBin(xb,yb,1), h.GetBin(xb,yb,0))
                mergeBins(h, h.GetBin(xb,yb,h.GetNbinsZ()), h.GetBin(xb,yb,h.GetNbinsZ()+1))
    if mergeY:
        for xb in range(h.GetNbinsX()+2):
            for zb in range(h.GetNbinsZ()+2):
                mergeBins(h, h.GetBin(xb,1,zb), h.GetBin(xb,0,zb))
                mergeBins(h, h.GetBin(xb,h.GetNbinsY(),zb), h.GetBin(xb,h.GetNbinsY()+1,zb))
    if mergeX:
        for yb in range(h.GetNbinsY()+2):
            for zb in range(h.GetNbinsZ()+2):
                mergeBins(h, h.GetBin(1,yb,zb), h.GetBin(0,yb,zb))
                mergeBins(h, h.GetBin(h.GetNbinsX(),yb,zb), h.GetBin(h.GetNbinsX()+1,yb,zb))


def appendFlowBin2d(h, mergeX=True, mergeY=True):
    if mergeY:
        for b in range(h.GetNbinsX()+2):
            mergeBins(h, h.GetBin(b,1), h.GetBin(b,0))
            mergeBins(h, h.GetBin(b,h.GetNbinsY()), h.GetBin(b,h.GetNbinsY()+1))
    if mergeX:
        for b in range(h.GetNbinsY()+2):
            mergeBins(h, h.GetBin(1,b), h.GetBin(0,b))
            mergeBins(h, h.GetBin(h.GetNbinsX(),b), h.GetBin(h.GetNbinsX()+1,b))

def integralAndError( h, binx1=0, binx2=-1, bins=True ):
    if not bins:
        binx1 = h.FindFixBin(binx1)
        binx2 = h.FindFixBin(binx2)
    e = ROOT.Double()
    c = h.IntegralAndError(binx1,binx2,e)
    return c,e


def getValAndError( val, err, sig=2 ):
    digit = sig - int(floor(log10(err))) - 1
    return ( round(val,digit), round(err,digit) )

def getValAndErrorStr( val, err, sig=2, digits=False, tex=False ):
    # if digits, treat 'sig' as number of siginifcant digits
    out = ""
    if digits:
        if sig == 0:
            out = "{:d} #pm {:d}".format(int(round(val)),int(round(err)))
        elif sig > 0:
            out = "{{:.{0:d}f}} #pm {{:.{0:d}f}}".format(sig).format(val,err)
        else:
            print "What are negative digits even meaning???"
    else:
        out = "{} #pm {}".format( getValAndError( val, err, sig ) )
    if tex:
        out = "${}$".format(out.replace("#","\\").replace(" ",""))
    return out


def getYAxisTitle( histo ):
    # returns e.g.: "Events / 10 GeV"
    if not style.divideByBinWidth:
        return "Events / Bin"
    yTitle = "Events"
    xaxis = histo.GetXaxis()
    binW = xaxis.GetBinWidth(1)
    binWmean = (xaxis.GetXmax()-xaxis.GetXmin())/xaxis.GetNbins()
    unit = "GeV" if "GeV" in xaxis.GetTitle() else None

    if abs(binW-binWmean) < 1e-6: #assume constant bin size
    #if abs(binW-binWmean) > -1000000000.: #assume constant bin size
        if abs(binW-1) < 1e-6:
            return yTitle

        # get two significant digits
        binW = getValAndError( 0, binW )[1]
        if binW.is_integer():
            binW = int(binW)
        if unit:
            return yTitle + " / " + str(binW) + " " + unit
        else:
            return yTitle + " / " + str(binW)
    else: # assume variable bin size
        if unit:
            return yTitle + " / " + unit
        else:
            return yTitle

def setYAxisTitle(h):
    h.SetYTitle(getYAxisTitle(h))

def getROC( hSig, hBkg, highX=True ):
    # highX: signal is at high values of the variable

    nRocBins = hSig.GetNbinsX()

    sigEff = []
    bkgEff = []

    sigDen = hSig.Integral()
    bkgDen = hBkg.Integral()
    if not sigDen or not bkgDen:
        print "Warning, signal or background histogram has no integral"
        return

    for i in range(1, nRocBins+1 ):
        if  highX:
            sigNum = hSig.Integral(i, nRocBins+1)
            bkgNum = hBkg.Integral(i, nRocBins+1)
        else:
            sigNum = hSig.Integral(1, i)
            bkgNum = hBkg.Integral(1, i)

        sigEff.append( sigNum / sigDen )
        bkgEff.append( bkgNum / bkgDen )

    rocGraph = ROOT.TGraph( nRocBins, numpy.array(bkgEff), numpy.array(sigEff) )
    rocGraph.SetTitle(";#varepsilon_{bkg};#varepsilon_{sig}")
    return rocGraph


def automaticRebinner( hlist, minEvents=3 ):
    # Ereates an array of bin edges on a list of histograms, such that each histogram
    # has at least 'maxEvents' events in each bin.
    if minEvents == 0: minEvents = 1e-10
    out = []
    tmp = [0]*len(hlist)
    nBins = hlist[0].GetNbinsX()

    # check overflow bin
    overflowList = [ h.GetBinContent( nBins+1 ) for h in hlist ]
    if min(overflowList) > minEvents:
        out.append( hlist[0].GetBinLowEdge(nBins+1) )

    for bin in range(nBins, -1, -1 ):

        contents = [ h.GetBinContent( bin ) for h in hlist ]
        tmp = [sum(x) for x in zip(tmp, contents)]

        # check upper bondary of last entry
        if not out:
            if contents != [0]*len(hlist):
                out.append( hlist[0].GetBinLowEdge(bin+1) )

        else:
            if min(contents) > minEvents:
                out.append( hlist[0].GetBinLowEdge(bin) )
                tmp = [0]*len(hlist)

    print out[::-1]

def getMinimum(hists):
    # Reset fMinimum, else 'GetMinimum' will return setted minimum
    return min([h.GetMinimum(0) for h in hists if not h.SetMinimum()])


def setMinMaxForLog():
    primitivesOnCanvas = [i for i in ROOT.gPad.GetCanvas().GetListOfPrimitives()]
    histograms = []
    stackedHistograms = []
    for h in primitivesOnCanvas:
        if isinstance(h, ROOT.THStack):
            stackedHistograms.append(h)
            for sh in h.GetStack():
                histograms.append(sh)
        elif isinstance(h, ROOT.TH1):
            histograms.append(h)
    maxC = max([h.GetMaximum() for h in histograms])
    minC = getMinimum(histograms)
    unity = 1./maxBinWidth(histograms[0]) if style.divideByBinWidth else 1.
    minimum = max([unity,minC]) if style.minimumOne else minC
    #minimum /= 9.
    #minimum /= 90.
    for i in histograms:
        i.SetMaximum(2.5*maxC)
        #i.SetMaximum(100.*maxC)
        #i.SetMaximum(100000.*maxC)
        i.SetMinimum(minimum)
    for s in stackedHistograms:
        s.SetMinimum(minimum)
        #s.SetMaximum(2.5*maxC)
        s.SetMaximum(10.*maxC)
    ROOT.gPad.Update()


def modifySaveName(name):
    replacements = {"(": "", ")": "", "&":"AND", ".":"p", "/":"DIV", "<":"", ">":"", "*":"TIMES", "$": "DOLLAR"}
    for a, b in replacements.iteritems():
        name = name.replace(a, b)
    return name

def save( name, folder="plots/", endings=[".pdf"], normal=True, log=True, changeMinMax=True ):
    name = modifySaveName(name)
    if normal:
        for ending in endings:
            ROOT.gPad.GetCanvas().SaveAs( folder+name+ending )
    if log:
        allH2s = [i for i in ROOT.gPad.GetCanvas().GetListOfPrimitives() if isinstance(i,ROOT.TH2)]
        if allH2s:
            ROOT.gPad.GetCanvas().SetLogz()
        else:
            if changeMinMax: setMinMaxForLog()
            ROOT.gPad.GetCanvas().SetLogy()
        for ending in endings:
            ROOT.gPad.GetCanvas().SaveAs( folder + name + "_log" + ending )


def getBinningsFromName( name ):
    out = { "": None }
    # get histogram name
    if binCfg.has_section(name):
        for binningName, binning in binCfg.items( name ):
            binning = [ float(x) for x in binning.split(" ") ]
            out[binningName] = binning
    return out

def drange(start, stop, n):
    out = [start]
    step = 1.*(stop-start)/n
    while out[-1] < stop:
        out.append( out[-1] + step )
    return out

def interpolate2D( h ):
    for xbin in range(1,h.GetNbinsX()+1):
        x = h.GetXaxis().GetBinCenter(xbin)
        for ybin in range(1,h.GetNbinsY()+1):
            y = h.GetYaxis().GetBinCenter(ybin)
            c = h.GetBinContent(xbin,ybin)
            if c: continue
            ctop = h.GetBinContent(xbin,ybin+1)
            cbot = h.GetBinContent(xbin,ybin-1)
            crig = h.GetBinContent(xbin+1,ybin)
            clef = h.GetBinContent(xbin-1,ybin)
            intPoints = []
            if cbot and ctop: intPoints.extend( [cbot, ctop])
            if crig and clef: intPoints.extend( [crig, clef])
            newC = sum(intPoints)/len(intPoints) if len(intPoints) else 0
            h.SetBinContent(xbin,ybin, newC )
    return h

def diagonalFlip( original ):
    # original, flipped are both TH2
    flipped = original.Clone(original.GetName()+"flipped")
    flipped.SetTitle("{};{};{}".format(
            original.GetTitle(),
            original.GetYaxis().GetTitle(),
            original.GetXaxis().GetTitle()
        ))

    for xbin in range(original.GetNbinsX()+2):
        for ybin in range(original.GetNbinsY()+2):
            flipped.SetBinContent( ybin, xbin, original.GetBinContent(xbin,ybin) )
            flipped.SetBinError( ybin, xbin, original.GetBinError(xbin,ybin) )
    return flipped

def drawDiagonal(h2, xmin=None):
    if not xmin: xmin = h2.GetXaxis().GetXmin()
    xmax = h2.GetXaxis().GetXmax()
    ymin = h2.GetYaxis().GetXmin()
    ymax = h2.GetYaxis().GetXmax()
    l = ROOT.TLine()
    l.SetLineColor(ROOT.kGray+2)
    l.SetLineStyle(2)
    l.DrawLine(max(xmin,ymin),max(xmin,ymin),min(xmax,ymax),min(xmax,ymax))

def addHists(*histograms):
    out = histograms[0].Clone()
    for h in histograms[1:]:
        out.Add(h)
    return out


def addHistUncert(*histograms):
    out = histograms[0].Clone()
    for h in histograms[1:]:
        for bin in range(out.GetNbinsX()+2):
            out.SetBinError(bin, sqrt( out.GetBinError(bin)**2 + h.GetBinError(bin)**2 ))
    return out

def maxBinWidth( h ):
    return max([ h.GetBinWidth(bin) for bin in range(h.GetNbinsX()+2) ])

def metricPrefix( n ):
    for unit in ["","K","M","G","T"]:
        if abs(n) < 1000:
            return "{:3.1f}{}".format(n, unit)
        n /= 1000
    return "{:3.1f}".format(n)

def loopH2( h2 ):
    return [(xbin,ybin) for xbin in range(h2.GetNbinsX()+2) for ybin in range(h2.GetNbinsY()+2)]

def loopH3( h3 ):
    return [(xbin,ybin,zbin) for xbin in range(h3.GetNbinsX()+2) for ybin in range(h3.GetNbinsY()+2) for zbin in range(h3.GetNbinsZ()+2)]

def loopH(h):
    if isinstance(h,ROOT.TH3):
        return loopH3(h)
    elif isinstance(h,ROOT.TH2):
        return loopH2(h)
    else:
        return range(h.GetNbinsX()+2)


def printH2(h2,flow=True):
    for xbin, ybin in loopH2(h2):
        if not flow and ( not xbin or not xbin or xbin==h2.GetNbinsX()+1 or ybin==h2.GetNbinsY()+1 ): continue
        print xbin, ybin, h2.GetBinContent(xbin,ybin), "±", h2.GetBinError(xbin,ybin)

def stdHist(dataset, name, binning=None, xCut=True, cut1=0, cut2=1e8):
    h = dataset.getHist(name)
    if not h: return
    if isinstance(h, ROOT.TH2):
        if xCut: h = h.ProjectionY(randomName(), h.GetXaxis().FindFixBin(cut1), h.GetYaxis().FindFixBin(cut2))
        else:    h = h.ProjectionX(randomName(), h.GetYaxis().FindFixBin(cut1), h.GetXaxis().FindFixBin(cut2))
    if binning: h = rebin(h, binning)
    #if (binning.any()): h = rebin(h, binning)
    appendFlowBin(h)
    h.SetYTitle(getYAxisTitle(h))
    return h

def getEnvelopeHists(hlist):
    hUp = hlist[0].Clone(randomName())
    hDn = hlist[0].Clone(randomName())

    for h in hlist:
        for b in loopH(h):
            c = h.GetBinContent(b)
            if c>hUp.GetBinContent(b): hUp.SetBinContent(b,c)
            if c<hDn.GetBinContent(b): hDn.SetBinContent(b,c)
    return hUp, hDn

def integerContent(h, scaledByWidth=False):
    for bin in loopH(h):
        c = h.GetBinContent(bin)
        if scaledByWidth:
            c *= h.GetBinWidth(bin)
        if abs(c-int(round(c)))>1e-5:
            return False
    return True

def drawOpt(h, style):
    if style == "data":
        h.SetLineColor(ROOT.kBlack)
        h.SetMarkerColor(ROOT.kBlack)
        h.SetMarkerStyle(20)
        #h.SetMarkerSize(0.7)
        h.SetMarkerSize(1.)
        h.drawOption_="pz"
        if isinstance(h, ROOT.TH1):
            h.SetBinErrorOption(ROOT.TH1.kPoisson)
            h.drawOption_="e0p0"
            if integerContent(h):
                h.Sumw2(False) # kPoisson uncertainties are drawn
    elif style == "pre":
        h.SetLineColor(ROOT.kBlack)
        h.drawOption_ = "hist"
    elif style == "signal":
        h.SetLineWidth(3)
        h.drawOption_ = "hist"
    elif style == "statUnc":
        h.SetLineWidth(5)
        h.SetMarkerStyle(0)
        h.SetLineColor(ROOT.kGray+2)
        h.drawOption_ = "e2x0"
        #h.SetFillStyle(3254)
        #h.drawOption_ = "e2"
    elif style == "totUnc":
        h.SetFillStyle(3254)
        h.SetMarkerSize(0)
        h.SetFillColor(ROOT.kBlack)
        h.drawOption_ = "e2"
    elif style == "sysUnc":
        h.SetFillStyle(3245)
        h.SetMarkerSize(0)
        h.SetFillColor(ROOT.kRed)
        h.drawOption_ = "e2"
    elif style == "sys":
        c = h.GetLineColor()
        h.SetFillColor(c)
        h.SetMarkerColor(c)
        h.SetFillStyle(3333)
        h.drawOption_ = "e2"
    else:
        print "Do not know what to do with draw option", style

def getPaletteColor(f):
    # f should be a fraction from 0 to 1, such that the whole palette is used
    return ROOT.TColor.GetColorPalette(int(f*ROOT.TColor.GetNumberOfColors()))

def getPoissonUnc(n):
    """http://prd.aps.org/abstract/PRD/v86/i1/e010001 (p.399)
    For the case of Poisson distributed n, the upper and lower limits on
    the mean value ν can be found from the Neyman procedure where the upper
    and lower limits are at confidence levels of 1 − α low and 1 − α up ,
    respectively, and F χ −1 2 is the quantile of the χ distribution
    (inverse of the cumulative distribution). The quantiles F χ −1 2 can be
    obtained from standard tables or from the ROOT routine TMath::ChisquareQuantile.
    For central confidence intervals at confidence level 1 − α, set α lo = α up = α/2.
    """

    # calculate x = 1-alpha ( approx 68% )
    x = ROOT.TMath.Erf(1./sqrt(2))
    alpha = 1-x

    # for central confidence intervals, alpha_lo =alpha_up = alpha/2
    alpha_lo = alpha/2
    alpha_up = alpha/2

    # confidence interval is [ xlo, xup ]
    xlo = 0.5 * ROOT.TMath.ChisquareQuantile( alpha_lo, 2*n )
    xup = 0.5 * ROOT.TMath.ChisquareQuantile( 1-alpha_up, 2*(n+1) )
    return n-xlo, xup-n

def getSysHisto(h, relUncert):
    hsys = h.Clone(randomName())
    for bin in range(hsys.GetNbinsX()+2):
        c = hsys.GetBinContent(bin)
        if c > 1e-10:
            e = relUncert*c
            hsys.SetBinError(bin, e)
        elif hsys.GetBinContent(bin-1) or hsys.GetBinContent(bin+1):
            # check if option "width" should be used
            # TODO: check if the weight agrees with the lumi+pu weight
            meanWeight = hsys.Integral(0,-1)/hsys.GetEntries()
            poissonZeroError = 1.8410216450098775
            e = meanWeight*poissonZeroError
            e /= hsys.GetBinWidth(bin) if style.divideByBinWidth else 1.
            hsys.SetBinError(bin, e)
    return hsys
    
def getSysHistoCut(h, relUncert1,relUncert2,cutBin):
    hsys = h.Clone(randomName())
    for bin in range(hsys.GetNbinsX()+2):
        if bin<cutBin:
            relUncert=relUncert1
        else:
            relUncert=relUncert2
        c = hsys.GetBinContent(bin)
        if c > 1e-10:
            e = relUncert*c
            hsys.SetBinError(bin, e)
        elif hsys.GetBinContent(bin-1) or hsys.GetBinContent(bin+1):
            # check if option "width" should be used
            # TODO: check if the weight agrees with the lumi+pu weight
            meanWeight = hsys.Integral(0,-1)/hsys.GetEntries()
            poissonZeroError = 1.8410216450098775
            e = meanWeight*poissonZeroError
            e /= hsys.GetBinWidth(bin) if style.divideByBinWidth else 1.
            hsys.SetBinError(bin, e)
    return hsys

def iterate(h, axis="x"):
    if axis == "x":
        return range(h.GetNbinsX()+2)
    elif axis == "y":
        return range(h.GetNbinsY()+2)
    elif axis == "z":
        return range(h.GetNbinsZ()+2)
    else:
        print "Please specify correct axis"

def myMatch( regex, string ):
    m = re.match( regex, string )
    return [ m ] if m else []

def getAxis(h, ax="x"):
    ax = ax.lower()
    if ax=="x":
        return h.GetXaxis()
    elif ax=="y":
        return h.GetYaxis()
    elif ax=="z":
        return h.GetZaxis()
    else:
        print "do not know what do do with ", ax

def getProjection(h, ax="x", cutBin1=0, cutBin2=-1):
    ax = ax.lower()
    if ax == "x":
        return h.ProjectionX(randomName(), cutBin1, cutBin2)
    elif ax == "y":
        return h.ProjectionY(randomName(), cutBin1, cutBin2)
    else:
        print "do not know what do do with ", ax

def getSystFromDifference(h1, h2, changeStyle=True):
    out = h1.Clone(randomName())
    if changeStyle: drawOpt(out, "sys")
    for bin in range(out.GetNbinsX()+2):
        c1 = h1.GetBinContent(bin)
        c2 = h2.GetBinContent(bin)
        out.SetBinContent(bin, 0.5*(c1+c2))
        out.SetBinError(bin, 0.5*abs(c1-c2))
    return out

def getSystFromEnvelopes(h0, hUp, hDn, changeStyle=True):
    out = h0.Clone(randomName())
    if changeStyle: drawOpt(out, "sys")
    for bin in range(out.GetNbinsX()+2):
        c = h0.GetBinContent(bin)
        up = hUp.GetBinContent(bin)
        dn = hDn.GetBinContent(bin)
        out.SetBinError(bin, max(abs(up-c), abs(dn-c)))
    return out

def getSystFromVariance(h0, hList):
    import numpy
    out = h0.Clone(randomName())
    for bin in loopH(h0):
        cs = []
        for h in hList: cs.append(h.GetBinContent(bin))
        mean = numpy.mean(cs)
        out.SetBinError(bin, max(numpy.std(cs),abs(mean-out.GetBinContent(bin))))
    return out

def dataCardToLatexTable(filename):
    import DatacardParser
    from optparse import OptionParser
    parser = OptionParser()
    DatacardParser.addDatacardParserOptions(parser)
    (options, args) = parser.parse_args()
    f = open(filename)
    dc = DatacardParser.parseCard(f,options)

    columns = [(" ", "SM", "Data", "Signal")]
    for bin in dc.bins:
        obs = int(dc.obs[bin])
        vals = {}
        for p in dc.processes:
            vals[p] = [dc.exp[bin][p],0]
        for uncert in dc.systs:
            for p, err in uncert[4][bin].iteritems():
                if err>1e-6:
                    vals[p][1] += (err-1)**2
        columns.append( [ \
            bin, \
            "{}\\pm{}".format(vals['bg'][0], vals['bg'][0]*sqrt(vals['bg'][1])), \
            str(obs), \
            "{}\\pm{}".format(vals['sig'][0], vals['sig'][0]*sqrt(vals['sig'][1]))])
    content = '\n'.join([' & '.join(x) for x in zip(*columns)])
    print "\\begin{tabular}{%s}\n"%('c'*len(columns)) + content + "\n\\end{tabular}\n"

############################### tree stuff ####################################

def createHistoFromTree(tree, variable, weight="", nBins=20, firstBin=None, lastBin=None ):
    """
    tree: tree to create histo from
    variable: variable to plot (must be a branch of the tree)
    weight: weights to apply (e.g. "var1*(var2 > 15)" will use weights from var1 and cut on var2 > 15
    nBins, firstBin, lastBin: number of bins, first bin and last bin (same as in TH1F constructor)
    nBins: if nBins is a list, and to a int, a user binned plot will be generated
    returns: histogram
    """
    from ROOT import TH1F
    name = randomName()
    if isinstance( nBins, list ):
        import array
        xBins = array.array('d', nBins )
        result = TH1F(name, variable, len(nBins)-1, xBins)
        result.Sumw2()
        tree.Draw("%s>>%s"%(variable, name), weight, "goff")
    elif firstBin==None and lastBin==None:
        import ROOT
        tree.Draw("%s>>%s(%s,,)"%(variable,name,nBins), weight, "goff")
        result = ROOT.gDirectory.Get( name )
        if isinstance( result, ROOT.TTree ):
            print "Warning, no entries"
            return ROOT.TH1F()
        result.Sumw2() # applying the errors here is perhaps not entirely correct
    else:
        result = TH1F(name, variable, nBins, firstBin, lastBin)
        result.Sumw2()
        tree.Draw("%s>>%s"%(variable, name), weight, "goff")
    if style.divideByBinWidth: result.Scale(1., "width")
    appendFlowBin(result)
    result.SetTitle("")
    if variable.startswith("met"):
        result.SetXTitle("#it{E}_{T}^{miss} (GeV)")
    if variable.startswith("metRaw"):
        result.SetXTitle("uncorrected #it{E}_{T}^{miss} (GeV)")
    setYAxisTitle(result)
    return result

def createHistoFromDatasetTree(dset, variable, weight, nBins, treename="tr/simpleTree"):
    tree = ROOT.TChain(treename)
    for f in dset.files: tree.Add(f)
    h = createHistoFromTree(tree, variable, weight, nBins)
    h.SetLineColor(dset.color)
    return h

def createHistoFromDatasetTreeWeighted(dset, variable, weight, nBins, treename="tr/simpleTree", weightStr="1.44225-0.000563764*emht+1.17322e-07*emht*emht"):
    tree = ROOT.TChain(treename)
    for f in dset.files: tree.Add(f)
    hUnw = createHistoFromTree(tree, variable, weight, nBins)
    hW = createHistoFromTree(tree, variable, weight+"*"+weightStr, nBins)
    hW.SetLineColor(dset.color)
    hW.Scale(hUnw.Integral(0,-1)/hW.Integral(0,-1))
    hSys = hW.Clone(randomName())
    drawOpt(hSys, "sys")
    for bin in range(hW.GetNbinsX()+2):
        hSys.SetBinError(bin,abs(hW.GetBinContent(bin)-hUnw.GetBinContent(bin)))
    return hW, hSys

def beautifyCutString(cut):
    cut = cut.replace("emht","#it{EMH}_{T}").replace("&&",",").replace("pt", "#it{p}_{T}")
    return cut

def blind(h, val=100):
    # Set bin content and error to zero for all bins larger than 'val'
    for b in range(h.FindFixBin(val), h.GetNbinsX()+2):
        h.SetBinContent(b, 0)
        h.SetBinError(b, 0)

def addPoissonUncertainty(hist):
    mediumWeight = hist.Integral(0,-1)/hist.GetEntries()
    for b in loopH(hist):
        c = hist.GetBinContent(b)
        if not c:
            hist.SetBinError(b, mediumWeight * 1.84102164458)
        elif abs((c-mediumWeight)/c) < 1e-2: # one effective entry
            hist.SetBinError(b, mediumWeight * 2.6378596228)

def addUncertaintiesQuadratic(hlist):
    out = hlist[0].Clone(randomName())
    for bin in loopH(out):
        out.SetBinError(bin, sqrt(sum([h.GetBinError(bin)**2 for h in hlist])))
    return out

intLumi = 35.867e3
#intLumi = 35.867e3 + 41.e3

class Label:
    # Create labels
    # Usage:
    # * With Labels(), all default labels will be printed
    # * With Labels(False), the method is only initiated and labels can be modified before calling the 'draw' method

    cmsEnergy = 13 #TeV

    def draw( self ):
        varDict = vars( self )
        for varName, obj in varDict.iteritems():
            if isinstance( obj, ROOT.TLatex ):
                obj.SetNDC()
                obj.Draw()

    #def __init__( self, drawAll=True, sim=False, status="Private Work", info="" ):
    def __init__( self, drawAll=True, sim=False, status="Work in Progress", info="" ):
        saveStuff.append(self)
        #if status == "Private Work":
        if status == "Work in Progress":
            if sim:
                #self.cms = ROOT.TLatex( 0.2, .887, "#scale[0.76]{#font[52]{Private Work Simulation}}" )
                #self.cms = ROOT.TLatex( 0.2, .95, "#scale[0.76]{#font[52]{Private Work Simulation}}" )
                #self.cms = ROOT.TLatex( 0.2, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress Simulation}}" )
                #self.cms = ROOT.TLatex( 0.15, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress Simulation}}" )
                self.cms = ROOT.TLatex( 0.16, .92, "#splitline{#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}}{#scale[0.76]{#font[52]{  Simulation}}}" )
            else:
                #self.pub = ROOT.TLatex( 0.2, .887, "#font[61]{CMS} #scale[0.76]{#font[52]{%s}}"%status )
                self.pub = ROOT.TLatex( 0.15, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{%s}}"%status )
        else:
            if sim:
                self.cms = ROOT.TLatex( 0.2, .887, "#font[61]{CMS} #scale[0.76]{#font[52]{Simulation}}" )
            else:
                self.cms = ROOT.TLatex( 0.2, .887, "#font[61]{CMS}" )
            self.pub = ROOT.TLatex( 0.2, .857, "#scale[0.76]{#font[52]{%s}}"%status )
        self.lum = ROOT.TLatex( .62, .95, "%.1f fb^{-1} (%s TeV)"%(intLumi/1000., self.cmsEnergy) )
        if info: self.info = ROOT.TLatex( .15, .95, info )
        #if info: self.info = ROOT.TLatex( .85, .85, info )

        if drawAll:
            self.draw()

class Label2D(Label):
    def __init__(self, drawAll=True, sim=False, status="", info=""):
        saveStuff.append(self)
        cmsText = ""
        if status != "Private Work": cmsText += "#font[61]{CMS}"
        if sim: cmsText += " #scale[0.76]{#font[52]{Simulation}}"
        if status: cmsText += "#scale[0.76]{#font[52]{%s}}"%status
        self.cms = ROOT.TLatex(.15, .95, cmsText)
        if info: self.info = ROOT.TLatex(0.2, .895, info)
        self.lum = ROOT.TLatex( .47, .95, "%.1f fb^{-1} (%s TeV)"%(intLumi/1000., self.cmsEnergy) )
        if drawAll: self.draw()
