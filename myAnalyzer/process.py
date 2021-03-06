#!/usr/bin/env python2
import argparse
import os
import multiprocessing
import glob
import subprocess

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# own libray
import run

import createConfig


ds={
    #"met": [
        #"MET_Run2016B-03Feb2017_ver2-v2_myTuple.root",
        #"MET_Run2016C-03Feb2017-v1_myTuple.root",
        #"MET_Run2016D-03Feb2017-v1_myTuple.root",
        #"MET_Run2016E-03Feb2017-v1_myTuple.root",
        #"MET_Run2016F-03Feb2017-v1_myTuple.root",
        #"MET_Run2016G-03Feb2017-v1_myTuple.root",
        #"MET_Run2016H-03Feb2017_ver2-v1_myTuple.root",
        #"MET_Run2016H-03Feb2017_ver3-v1_myTuple.root" 
        #],
    #"ee": [
        #"DoubleEG_Run2016B-03Feb2017_ver2-v2_myTuple.root",
        #"DoubleEG_Run2016C-03Feb2017-v1_myTuple.root",
        #"DoubleEG_Run2016D-03Feb2017-v1_myTuple.root",
        #"DoubleEG_Run2016E-03Feb2017-v1_myTuple.root",
        #"DoubleEG_Run2016F-03Feb2017-v1_myTuple.root",
        #"DoubleEG_Run2016G-03Feb2017-v1_myTuple.root",
        #"DoubleEG_Run2016H-03Feb2017_ver2-v1_myTuple.root",
        #"DoubleEG_Run2016H-03Feb2017_ver3-v1_myTuple.root" 
        #],
    #"mm": [
        #"DoubleMuon_Run2016B-03Feb2017_ver2-v2_myTuple.root",
        #"DoubleMuon_Run2016C-03Feb2017-v1_myTuple.root",
        #"DoubleMuon_Run2016D-03Feb2017-v1_myTuple.root",
        #"DoubleMuon_Run2016E-03Feb2017-v1_myTuple.root",
        #"DoubleMuon_Run2016F-03Feb2017-v1_myTuple.root",
        #"DoubleMuon_Run2016G-03Feb2017-v1_myTuple.root",
        #"DoubleMuon_Run2016H-03Feb2017_ver2-v1_myTuple.root",
        #"DoubleMuon_Run2016H-03Feb2017_ver3-v1_myTuple.root"
        #],
    #"em": [
        #"MuonEG_Run2016B-03Feb2017_ver2-v2_myTuple.root",
        #"MuonEG_Run2016C-03Feb2017-v1_myTuple.root",
        #"MuonEG_Run2016D-03Feb2017-v1_myTuple.root",
        #"MuonEG_Run2016E-03Feb2017-v1_myTuple.root",
        #"MuonEG_Run2016F-03Feb2017-v1_myTuple.root",
        #"MuonEG_Run2016G-03Feb2017-v1_myTuple.root",
        #"MuonEG_Run2016H-03Feb2017_ver2-v1_myTuple.root",
        #"MuonEG_Run2016H-03Feb2017_ver3-v1_myTuple.root"
        #],
    #"ht": [
        #"JetHT_Run2016B-03Feb2017_ver2-v2_myTuple.root",
        #"JetHT_Run2016C-03Feb2017-v1_myTuple.root",
        #"JetHT_Run2016D-03Feb2017-v1_myTuple.root",
        #"JetHT_Run2016E-03Feb2017-v1_myTuple.root",
        #"JetHT_Run2016F-03Feb2017-v1_myTuple.root",
        #"JetHT_Run2016G-03Feb2017-v1_myTuple.root",
        #"JetHT_Run2016H-03Feb2017_ver2-v1_myTuple.root",
        #"JetHT_Run2016H-03Feb2017_ver3-v1_myTuple.root"
        #],
    "dy": [
        "DYJetsToLL_M-50-amcatnloFXFX_ext_myTuple.root"
        ],
    "ttg": [
        "TTGamma_Dilept-amcatnlo_myTuple.root",
        "TTGamma_Hadronic-amcatnlo_myTuple.root",
        "TTGamma_SingleLeptFromT-amcatnlo_myTuple.root",
        "TTGamma_SingleLeptFromTbar-amcatnlo_myTuple.root",
        #####"TTGJets_Total_myTuple.root"
        ],
    "wwg": [
        "WWG-amcatnlo_ext_myTuple.root"
        ],
    "ww": [
        "WWTo2L2Nu_myTuple.root"
        ],
    "wzg": [
        "WZG-amcatnlo_myTuple.root"
        ],
    "wz": [
        "WZTo3LNu_Total_myTuple.root"
        ],
    "zg": [
        ######"ZGTo2LG_ext_myTuple.root",
        "ZGTo2LG_PtG-130_myTuple.root",
        #####"ZGTo2LG_myTuple.root"
        "ZGTo2LG_Total_myTuple.root"
        ],
    "zz": [
        #######"ZZTo2L2Nu_powheg_ext1_myTuple.root",
        #######"ZZTo2L2Nu_powheg_myTuple.root",
        "ZZTo2L2Nu_powheg_Total_myTuple.root",
        ######"ZZTo4L_powheg_ext1_myTuple.root",
        #######"ZZTo4L_powheg_myTuple.root"
        "ZZTo4L_powheg_Total_myTuple.root"
        ],
    "tt":[
        "TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_myTuple.root",
        #############"../mediumIDPOG_TopPt_NIsr/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_myTuple.root",
        ############"TT_myTuple.root"
        ],
    "singletop":[
        "ST_s-channel_4f_leptonDecays-amcatnlo_myTuple.root",
        "ST_t-channel_antitop_4f_inclusiveDecaysV2_myTuple.root",  #just 9
        "ST_t-channel_top_4f_inclusiveDecaysV2_myTuple.root",     #just 9
        "ST_tW_antitop_5f_NoFullyHadronicDecays_ext_myTuple.root", ##NO PDF WEIGHTS
        "ST_tWll_5f_LO-MadGraph_myTuple.root",
        "ST_tW_top_5f_NoFullyHadronicDecays_ext_myTuple.root"##NO PDF WEIGHTS
    ],
    "wjets":[
        #########"WJetsToLNu-amcatnloFXFX_ext_myTuple.root",
        "WJetsToLNu-amcatnloFXFX_Total_myTuple.root",
        ########"WJetsToLNu-amcatnloFXFX_myTuple.root"#,
        ],
    "wgamma":[
        "WGToLNuG-amcatnloFXFX_ext_myTuple.root"
    ],
"signal":[
        #"SMS-T5bbbbZg_1800_1700_myTuple.root",
        #"SMS-T5bbbbZg_1800_400_myTuple.root",
        #"SMS-T5bbbbZg_1800_600_myTuple.root",
        "SMS-T5bbbbZg_1500_1400_myTuple.root",
        "SMS-T5bbbbZg_1500_400_myTuple.root",
        "SMS-T5bbbbZg_1500_600_myTuple.root",
        #"SMS-T5bbbbZg_myTuple.root",
        #"SMS-T6ttZg_myTuple.root",
        #"SMS-T6ttZg_600_200_myTuple.root",
        #"SMS-T6ttZg_400_200_myTuple.root",
        #"SMS-TChiNG_BF50N50G_myTuple.root",
        #"GMSB_GravitinoLSP_N1decays_myTuple.root",
        "SMS-TChiNG_400_myTuple.root",
        "SMS-TChiNG_600_myTuple.root",
        "SMS-TChiNG_1200_myTuple.root",
        "SMS-GMSB_240_230_myTuple.root",
        "SMS-GMSB_290_205_myTuple.root",
        "SMS-GMSB_415_355_myTuple.root",
        #"GGM_GravitinoLSP_M1-200to1500_M2-200to1500_myTuple.root",
        #"GGM_GravitinoLSP_M1-50to1500_M3-1000to2500_myTuple.root"
        ]
}
    
    
#dir = "../scratch/v12/"
#dir = "../scratch/tempTrees/tightMVAID/"
#dir = "../scratch/tempTrees/mediumID/"
#dir = "../scratch/tempTrees/mediumIDPOG/"
dir = "../scratch/tempTrees/AN/"
#dir = "../scratch/tempTrees/AN_rishi/"
#dir = "../scratch/tempTrees/AN/ht/"
#dir = "../scratch/tempTrees/AN/htPure/"
#dir = "../scratch/tempTrees/AN/met/"
#dir = "../scratch/tempTrees/mediumIDPOG/met/"
#dir = "../scratch/tempTrees/mediumIDPOG_noTopPt_noNIsr/ht/"
#dir = "../scratch/tempTrees/mediumIDPOG_noTopPt_noNIsr/htPure/"
#dir = "../scratch/tempTrees/mediumIDPOG_noTopPt_noNIsr/"
#dir = "../scratch/tempTrees/mediumIDPOG_TopPt_NIsr/"
#dir = "../scratch/tempTrees/mediumIDPOG_noISRtop/"
#dir = "../scratch/tempTrees/mediumIDPOG/ht/"

#############################################
# Select datasets to process
#############################################

# compile only
run.run()

parser = argparse.ArgumentParser()
parser.add_argument('datasets', nargs='+', default=["all"], help="all "+' '.join(ds.keys()))
parser.add_argument('--condor', action="store_true")
#parser.add_argument('--ext', action='store_true')
args = parser.parse_args()


if args.datasets == ["all"]:
    toProcess = [x for sublist in ds.values() for x in sublist]
#elif args.datasets == ["2"]:
    #toProcess = [x for sublist in ds.values() for x in sublist if not x.startswith("GJet") and not x.startswith("QCD")]
else:
    toProcess = ds[args.datasets[0]]
    for n in args.datasets[1:]:
        toProcess += ds[n]

print toProcess
print len(toProcess)

if args.condor:
    #extStr = "--ext" if args.ext else ""
    extStr = ""
    for x in toProcess:
        with open("submitCondor.txt","w") as f:
            f.write("""
Universe   = vanilla
Executable = run.sh
Arguments  = {2}{0} {1}
Log        = logs/{0}.log
Output     = logs/{0}.out
Error      = logs/{0}.error
Queue
""".format(x, extStr, dir))
        subprocess.call(["condor_submit", "submitCondor.txt"])

else: # local processing
    files = [dir+x for x in toProcess]
    files.sort(key=os.path.getsize, reverse=True)
    #p = multiprocessing.Pool(8)
    p = multiprocessing.Pool(24)
    p.map(run.runExt, files)
