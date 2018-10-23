#!/usr/bin/env python2

from include import *
import DatacardParser
import multiprocessing
import glob
import numpy as np
import pickle as pkl


def getSignalUncertainties(datacard):
    uncDict = {}
    bkgs = ["tt", "dy", "zz", "wz", "other"]

    binNames = ["binMC_1", "binMC_2"]

    dc = limitTools.MyDatacard(datacard)
    systDict = dict([(l[0], l) for l in dc.systs])
    uncs = ["JES", "JER", "PU", "ISR", "genmet", "lepSF", "photonSF", "scale",
            "pdf", "SFDY", "SFTT", "SFZZ", "SFWZ", "xSecOTHER", "lumi", "trigger", "ttStatbinMC_1", "ttStatbinMC_2", "wzStatbinMC_1", "wzStatbinMC_2", "zzStatbinMC_1", "zzStatbinMC_2", "dyStatbinMC_1", "dyStatbinMC_2", "otherStatbinMC_1", "otherStatbinMC_2"]
    systDict = dict([(l[0], l) for l in dc.systs])
    for bkg in bkgs:
        uncDict[bkg] = {}
        for b in dc.bins:
            uncDict[bkg][b] = {}
            for unc in uncs:
                if (not systDict.get(unc)):
                    uncDict[bkg][b][unc] = 0.
                else:
                    uncDict[bkg][b][unc] = systDict[unc][4][b][bkg]
            # uncDict[bkg][b]["yield"]=systDict
            # print  dict([(l[0], l) for l in dc.systs])   # print uncDict
        for binName in binNames:
            uncDict[bkg][binName]["rate"] = dc.exp[binName][bkg]
            # print dc.exp[binName][bkg]
    return uncDict


a = getSignalUncertainties("testDatacard.txt")
# for b in a:
#     for c in a[b]:
#         print b, c, a[b][c]
tt_yield1 = a["tt"]["binMC_1"]["rate"]
tt_yield2 = a["tt"]["binMC_2"]["rate"]
dy_yield1 = a["dy"]["binMC_1"]["rate"]
dy_yield2 = a["dy"]["binMC_2"]["rate"]
zz_yield1 = a["zz"]["binMC_1"]["rate"]
zz_yield2 = a["zz"]["binMC_2"]["rate"]
wz_yield1 = a["wz"]["binMC_1"]["rate"]
wz_yield2 = a["wz"]["binMC_2"]["rate"]
other_yield1 = a["other"]["binMC_1"]["rate"]
other_yield2 = a["other"]["binMC_2"]["rate"]

tt_error1 = 0.
tt_error2 = 0.
dy_error1 = 0.
dy_error2 = 0.
zz_error1 = 0.
zz_error2 = 0.
wz_error1 = 0.
wz_error2 = 0.
other_error1 = 0.
other_error2 = 0.

for unc in ["JES", "JER", "PU", "ISR", "genmet", "lepSF", "photonSF", "scale",
            "pdf", "SFDY", "SFTT", "SFZZ", "SFWZ", "xSecOTHER", "lumi", "trigger"]:
    tt_error1 = tt_error1 + \
        (tt_yield1 * (a["tt"]["binMC_1"][unc] - 1.)
         )**2. if a["tt"]["binMC_1"][unc] > 0 else tt_error1 + 0
    tt_error2 = tt_error2 + \
        (tt_yield2 * (a["tt"]["binMC_2"][unc] - 1.)
         )**2. if a["tt"]["binMC_2"][unc] > 0 else tt_error2 + 0
    dy_error1 = dy_error1 + \
        (dy_yield1 * (a["dy"]["binMC_1"][unc] - 1.)
         )**2. if a["dy"]["binMC_1"][unc] > 0 else dy_error1 + 0
    dy_error2 = dy_error2 + \
        (dy_yield2 * (a["dy"]["binMC_2"][unc] - 1.)
         )**2. if a["dy"]["binMC_2"][unc] > 0 else dy_error2 + 0
    zz_error1 = zz_error1 + \
        (zz_yield1 * (a["zz"]["binMC_1"][unc] - 1.)
         )**2. if a["zz"]["binMC_1"][unc] > 0 else zz_error1 + 0
    zz_error2 = zz_error2 + \
        (zz_yield2 * (a["zz"]["binMC_2"][unc] - 1.)
         )**2. if a["zz"]["binMC_2"][unc] > 0 else zz_error2 + 0
    wz_error1 = wz_error1 + \
        (wz_yield1 * (a["wz"]["binMC_1"][unc] - 1.)
         )**2. if a["wz"]["binMC_1"][unc] > 0 else wz_error1 + 0
    wz_error2 = wz_error2 + \
        (wz_yield2 * (a["wz"]["binMC_2"][unc] - 1.)
         )**2. if a["wz"]["binMC_2"][unc] > 0 else wz_error2 + 0
    other_error1 = other_error1 + \
        (other_yield1 * (a["other"]["binMC_1"][unc] - 1.)
         )**2. if a["other"]["binMC_1"][unc] > 0 else other_error1 + 0
    other_error2 = other_error2 + \
        (other_yield2 * (a["other"]["binMC_2"][unc] - 1.)
         )**2. if a["other"]["binMC_2"][unc] > 0 else other_error2 + 0

tt_error1 = np.sqrt(tt_error1)
tt_error2 = np.sqrt(tt_error2)
dy_error1 = np.sqrt(dy_error1)
dy_error2 = np.sqrt(dy_error2)
zz_error1 = np.sqrt(zz_error1)
zz_error2 = np.sqrt(zz_error2)
wz_error1 = np.sqrt(wz_error1)
wz_error2 = np.sqrt(wz_error2)
other_error1 = np.sqrt(other_error1)
other_error2 = np.sqrt(other_error2)

tt_stat1 = tt_yield1 * (a["tt"]["binMC_1"]["ttStatbinMC_1"] - 1.)
tt_stat2 = tt_yield2 * (a["tt"]["binMC_2"]["ttStatbinMC_2"] - 1.)
dy_stat1 = dy_yield1 * (a["dy"]["binMC_1"]["dyStatbinMC_1"] - 1.)
dy_stat2 = dy_yield2 * (a["dy"]["binMC_2"]["dyStatbinMC_2"] - 1.)
zz_stat1 = zz_yield1 * (a["zz"]["binMC_1"]["zzStatbinMC_1"] - 1.)
zz_stat2 = zz_yield2 * (a["zz"]["binMC_2"]["zzStatbinMC_2"] - 1.)
wz_stat1 = wz_yield1 * (a["wz"]["binMC_1"]["wzStatbinMC_1"] - 1.)
wz_stat2 = wz_yield2 * (a["wz"]["binMC_2"]["wzStatbinMC_2"] - 1.)
other_stat1 = other_yield1 * (a["other"]["binMC_1"]["otherStatbinMC_1"] - 1.)
other_stat2 = other_yield2 * (a["other"]["binMC_2"]["otherStatbinMC_2"] - 1.)


print "tt 1: ", np.round(
    tt_yield1, 3), " +- ", np.round(tt_stat1, 3), " (stat) +- ", np.round(tt_error1, 3), "(syst)"
print "tt 2: ", np.round(
    tt_yield2, 3), " +- ", np.round(tt_stat2, 3), " (stat) +- ", np.round(tt_error2, 3), "(syst)"
print "dy 1: ", np.round(
    dy_yield1, 3), " +- ", np.round(dy_stat1, 3), " (stat) +- ", np.round(dy_error1, 3), "(syst)"
print "dy 2: ", np.round(
    dy_yield2, 3), " +- ", np.round(dy_stat2, 3), " (stat) +- ", np.round(dy_error2, 3), "(syst)"
print "zz 1: ", np.round(
    zz_yield1, 3), " +- ", np.round(zz_stat1, 3), " (stat) +- ", np.round(zz_error1, 3), "(syst)"
print "zz 2: ", np.round(
    zz_yield2, 3), " +- ", np.round(zz_stat2, 3), " (stat) +- ", np.round(zz_error2, 3), "(syst)"
print "wz 1: ", np.round(
    wz_yield1, 3), " +- ", np.round(wz_stat1, 3), " (stat) +- ", np.round(wz_error1, 3), "(syst)"
print "wz 2: ", np.round(
    wz_yield2, 3), " +- ", np.round(wz_stat2, 3), " (stat) +- ", np.round(wz_error2, 3), "(syst)"
print "other 1: ", np.round(
    other_yield1, 3), " +- ", np.round(other_stat1, 3), " (stat) +- ", np.round(other_error1, 3), "(syst)"
print "other 2: ", np.round(
    other_yield2, 3), " +- ", np.round(other_stat2, 3), " (stat) +- ", np.round(other_error2, 3), "(syst)"
print "total 1: ", np.round(
    tt_yield1 + dy_yield1 + zz_yield1 + wz_yield1 + other_yield1, 3), " +- ", np.round(np.sqrt(other_stat1**2. + tt_stat1**2. + dy_stat1**2. + zz_stat1**2. + wz_stat1**2.), 3), " (stat) +- ", np.round(np.sqrt(other_error1**2. + tt_error1**2. + dy_error1**2. + zz_error1**2. + wz_error1**2.), 3), "(syst)"
print "total 2: ", np.round(
    tt_yield2 + dy_yield2 + zz_yield2 + wz_yield2 + other_yield2, 3), " +- ", np.round(np.sqrt(other_stat2**2. + tt_stat2**2. + dy_stat2**2. + zz_stat2**2. + wz_stat2**2.), 3), " (stat) +- ", np.round(np.sqrt(other_error2**2. + tt_error2**2. + dy_error2**2. + zz_error2**2. + wz_error2**2.), 3), "(syst)"
