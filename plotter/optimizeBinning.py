#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
if sys.version_info[:2] == (2,6):
    print "Initialize correct python version first!"
    sys.exit()

import ConfigParser
import ROOT
import math
import argparse
import re
from random import randint
from sys import maxint

# private libs
import ratio
import style
import multiplot
from rwthColors import rwth
import limitTools

import auxiliary as aux
from datasets import *

def writeSimplifiedDatacard( bkgInts, sigInts, name="tmpDataCard.txt" ):
    # assuming no observation, 50% bkg uncert, 15% signal uncert
    nBins = len(bkgInts)
    header = """
imax %d
jmax 1
kmax 2
bin %s
observation %s\n"""%(nBins," ".join(["b%d"%i for i in range(nBins)])," ".join([str(int(i)) for i in bkgInts]) )
    infos = [["bin", "process", "process", "rate", "signal lnN", "bgUnc lnN"]]
    for iBin in range(len(bkgInts)):
        bkgInfos = ["b%d"%iBin, "bg", 1, bkgInts[iBin], "-", 1.5 ]
        sigInfos = ["b%d"%iBin, "sig", 0, sigInts[iBin], 1.15, "-" ]
        bkgInfos = [ str(i) for i in bkgInfos if i is not None ]
        sigInfos = [ str(i) for i in sigInfos if i is not None ]
        infos.append(bkgInfos)
        infos.append(sigInfos)
    for l in zip(*infos):
        header += " ".join(l) + "\n"

    with open(name,"w") as f:
        f.write(header)

bkgSet = zgamma+ ttgamma+ zz+ wwgamma+ wzgamma+DYjetsNLO+wjets+tt+singletop+wz+ww+zz4l+wjets
#sigSet = signal["T5Wg_1550_100"]
sigSet = tching_600
#sigSet = gmsb_240_230
#sigSet = gmsb_290_205
#sigSet = t5bbbbzg_1500_600
#sigSet = t5bbbbzg_1500_1400

#bkgHist = bkgSet.getHist("tr/met")
#sigHist = sigSet.getHist("tr/met")
#bkgHist = bkgSet.getHist("onZMet150/LL/met")
bkgHist = bkgSet.getHistWithoutNGen("xx_0_0/sig/LL/nom/met")
#sigHist = sigSet.getHist("onZMet150/LL/met")
sigHist = sigSet.getHistWithoutNGen("Ng_0_0/sig/LL/nom/met")
#bkgHist = bkgSet.getHist("onZMet150/LL/m_llg")
#sigHist = sigSet.getHist("onZMet150/LL/m_llg")
#bkgHist = bkgSet.getHist("onZMet150/LL/zpt")
#sigHist = sigSet.getHist("onZMet150/LL/zpt")
#bkgHist = bkgSet.getHist("onZMet150/LL/deltaPhiLL")
#sigHist = sigSet.getHist("onZMet150/LL/deltaPhiLL")

def frange(start, end, step):
    a=[]
    tmp = start
    while(tmp < end):
        a.append(tmp)
        tmp += step
    return a

#rebinning=frange(150.,5001.,5.)
#rebinning=frange(0.,1501.,25.)
rebinning=frange(0.,1501.,50.)
#rebinning=frange(0.,4.5,0.1)

bkgHist=aux.rebin(bkgHist,rebinning)
sigHist=aux.rebin(sigHist,rebinning)

binBoarders = [ bkgHist.GetNbinsX()+10 ]
#print binBoarders
oldR = 1e20 # maxint wouly be nicer


#for bin in range( sigHist.GetNbinsX()+2, 0, -1 ):
for bin in range( sigHist.GetNbinsX()+2, 0, -1 ):

    bkgInts = [ bkgHist.Integral(binBoarders[i+1],binBoarders[i]) for i in range(len(binBoarders)-1) ]
    sigInts = [ sigHist.Integral(binBoarders[i+1],binBoarders[i]) for i in range(len(binBoarders)-1) ]

    bkgInts.append( bkgHist.Integral(bin,binBoarders[-1]) )
    sigInts.append( sigHist.Integral(bin,binBoarders[-1]) )



    if min(bkgInts) <1e-6: continue
    if min(sigInts) <1e-6: continue


    writeSimplifiedDatacard( bkgInts, sigInts,name = "tmpDataCard_"+sigSet.names[0]+"_.txt" )
    r = limitTools.infosFromDatacard("tmpDataCard_"+sigSet.names[0]+"_.txt")["exp"]
    print bin,bkgHist.GetBinLowEdge(bin),r
    if (r - oldR)/r>0.05: #change must me larger than 5%
        print "append"
        binBoarders.append(bin)
    oldR = r

print "final r =", r
print binBoarders
metBoarders = [ sigHist.GetBinLowEdge(i) for i in binBoarders ]
metBoarders.reverse()
print metBoarders[0:-1]
