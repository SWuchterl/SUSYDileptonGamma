#!/usr/bin/env python2
import suppressor
with suppressor.suppress_stdout_stderr():
    import ROOT
import ROOT
import argparse
import os
#from configs import runConfig as config
import ConfigParser




def run(infile="", selector="myAnalyzer.cc", ext=False):
    # load libraries
    ROOT.gSystem.Load("pluginTreeWriterTreeWriterAuto.so")
    #ROOT.gSystem.Load("MT2Functor_cc.so")
    #ROOT.gSystem.Load("RoccoR_cc.so")
    #ROOT.gSystem.Load("rochcor2016_cc.so")
    #lib = "AutoDict_map_int_pair_int_int____cxx.so"
    #if not os.path.isfile(lib):
        #ROOT.gInterpreter.GenerateDictionary("map<int,pair<int,int> >", "map")
    #ROOT.gSystem.Load("AutoDict_map_int_pair_int_int____cxx.so")

    if infile:
        ch = ROOT.TChain("Tree")
        ch.AddFile(infile)
        extName = infile.replace("_nTuple", "_extNEIN_nTuple")
        ch.Process(selector+"+")
    elif ROOT.TSelector.GetSelector(selector+"++"):
        print "Compiled TSelector"
    else:
        raise Exception ('TSelector could not be compiled!')

def runExt(infile="", selector="myAnalyzer.cc"):
    # wrapper
    run(infile, selector, True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('file', default="", nargs="?")
    parser.add_argument('--ext', action='store_true')
    parser.add_argument('--signal', action='store_true')

    import createConfig

    args = parser.parse_args()
    signalScans = ["nothing"]
    

    if os.path.basename(args.file) in signalScans or args.signal:
        run(args.file, "SignalScan.cc")
    else:
        run(args.file, ext=args.ext)
