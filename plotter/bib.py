import array
import ConfigParser
import itertools
from math import *
import numpy
import pickle
import re
import ROOT
import style

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


def rebin( h, binEdges, scale=True ):
    #if not binEdges: return h
    checkRebinningConsistence( h.GetXaxis(), binEdges )
    binEdgesArr = array.array( 'd', binEdges )
    hnew = h.Rebin( len(binEdges)-1, "new", binEdgesArr )
    hnew.drawOption_ = h.drawOption_ if hasattr( h, "drawOption_" ) else ""
    #if scale: hnew.Scale( 1., "width" )
    if style.divideByBinWidth: hnew.Scale(1.,"width")
    return hnew
