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
    "ee": [
        "DoubleEG_Run2016B-03Feb2017_ver2-v2_nTuple.root",
        "DoubleEG_Run2016C-03Feb2017-v1_nTuple.root",
        "DoubleEG_Run2016D-03Feb2017-v1_nTuple.root",
        "DoubleEG_Run2016E-03Feb2017-v1_nTuple.root",
        "DoubleEG_Run2016F-03Feb2017-v1_nTuple.root",
        "DoubleEG_Run2016G-03Feb2017-v1_nTuple.root",
        "DoubleEG_Run2016H-03Feb2017_ver2-v1_nTuple.root",
        "DoubleEG_Run2016H-03Feb2017_ver3-v1_nTuple.root" 
        ],
    "mm": [
        "DoubleMuon_Run2016B-03Feb2017_ver2-v2_nTuple.root",
        "DoubleMuon_Run2016C-03Feb2017-v1_nTuple.root",
        "DoubleMuon_Run2016D-03Feb2017-v1_nTuple.root",
        "DoubleMuon_Run2016E-03Feb2017-v1_nTuple.root",
        "DoubleMuon_Run2016F-03Feb2017-v1_nTuple.root",
        "DoubleMuon_Run2016G-03Feb2017-v1_nTuple.root",
        "DoubleMuon_Run2016H-03Feb2017_ver2-v1_nTuple.root",
        "DoubleMuon_Run2016H-03Feb2017_ver3-v1_nTuple.root"
        ],
    "em": [
        "MuonEG_Run2016B-03Feb2017_ver2-v2_nTuple.root",
        "MuonEG_Run2016C-03Feb2017-v1_nTuple.root",
        "MuonEG_Run2016D-03Feb2017-v1_nTuple.root",
        "MuonEG_Run2016E-03Feb2017-v1_nTuple.root",
        "MuonEG_Run2016F-03Feb2017-v1_nTuple.root",
        "MuonEG_Run2016G-03Feb2017-v1_nTuple.root",
        "MuonEG_Run2016H-03Feb2017_ver2-v1_nTuple.root",
        "MuonEG_Run2016H-03Feb2017_ver3-v1_nTuple.root"
        ],
    "ht": [
        "JetHT_Run2016B-03Feb2017_ver2-v2_nTuple.root",
        "JetHT_Run2016C-03Feb2017-v1_nTuple.root",
        "JetHT_Run2016D-03Feb2017-v1_nTuple.root",
        "JetHT_Run2016E-03Feb2017-v1_nTuple.root",
        "JetHT_Run2016F-03Feb2017-v1_nTuple.root",
        "JetHT_Run2016G-03Feb2017-v1_nTuple.root",
        "JetHT_Run2016H-03Feb2017_ver2-v1_nTuple.root",
        "JetHT_Run2016H-03Feb2017_ver3-v1_nTuple.root"
        ],
    "dy": [
        "DYJetsToLL_M-50-amcatnloFXFX_ext_nTuple.root"
        ],
    "ttg": [
        "TTGamma_Dilept-amcatnlo_nTuple.root"
        ],
    "wwg": [
        "WWG-amcatnlo_ext_nTuple.root"
        ],
    "ww": [
        "WWTo2L2Nu_nTuple.root"
        ],
    "wzg": [
        "WZG-amcatnlo_nTuple.root"
        ],
    "wz": [
        "WZTo3LNu_Total_nTuple.root"
        ],
    "zg": [
        ##############"ZGTo2LG_ext_nTuple.root",
        "ZGTo2LG_PtG-130_nTuple.root",
        ########"ZGTo2LG_nTuple.root"
        "ZGTo2LG_Total_nTuple.root"
        ],
    "zz": [
        #########"ZZTo2L2Nu_powheg_ext1_nTuple.root",
        ###########"ZZTo2L2Nu_powheg_nTuple.root",
        "ZZTo2L2Nu_powheg_Total_nTuple.root",
        ###########"ZZTo4L_powheg_ext1_nTuple.root",
        ############"ZZTo4L_powheg_nTuple.root"
        "ZZTo4L_powheg_Total_nTuple.root"
        ],
    "tt":[
        "TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_nTuple.root"
        ],
    "singletop":[
        "ST_s-channel_4f_leptonDecays-amcatnlo_nTuple.root",
        "ST_t-channel_antitop_4f_inclusiveDecaysV2_nTuple.root",
        "ST_t-channel_top_4f_inclusiveDecaysV2_nTuple.root",
        "ST_tW_antitop_5f_NoFullyHadronicDecays_ext_nTuple.root",
        "ST_tWll_5f_LO-MadGraph_nTuple.root",
        "ST_tW_top_5f_NoFullyHadronicDecays_ext_nTuple.root"
    ],
    "wjets":[
        #############"WJetsToLNu-amcatnloFXFX_ext_nTuple.root",
        "WJetsToLNu-amcatnloFXFX_Total_nTuple.root",
        #############"WJetsToLNu-amcatnloFXFX_nTuple.root"#,
        ],
    "wgamma":[
        "WGToLNuG-amcatnloFXFX_ext_nTuple.root"
    ],
"signal":[
        ######"SMS-T5bbbbZg_1800_1700_nTuple.root",
        ######"SMS-T5bbbbZg_1800_400_nTuple.root",
        ######"SMS-T5bbbbZg_1800_600_nTuple.root",
        "SMS-T5bbbbZg_1500_1400_nTuple.root",
        "SMS-T5bbbbZg_1500_400_nTuple.root",
        "SMS-T5bbbbZg_1500_600_nTuple.root",
        "SMS-T5bbbbZg_nTuple.root",
        "SMS-TChiNG_BF50N50G_nTuple.root",
        "GMSB_GravitinoLSP_N1decays_nTuple.root",
        "SMS-TChiNG_400_nTuple.root",
        "SMS-TChiNG_600_nTuple.root",
        "SMS-TChiNG_1200_nTuple.root",
        "GGM_GravitinoLSP_M1-200to1500_M2-200to1500_nTuple.root",
        "GGM_GravitinoLSP_M1-50to1500_M3-1000to2500_nTuple.root"
        ]
}
    
    
#dir = "../scratch/v12/"
dir = "../scratch/v13/"

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
    p = multiprocessing.Pool(8)
    p.map(run.runExt, files)
