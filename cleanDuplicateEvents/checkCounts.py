#!/usr/bin/env python2

import ROOT
ROOT.gErrorIgnoreLevel=ROOT.kError
from ROOT import *

#nTuple trees to clean:
tree_path="../scratch/v01/"
#EM_data=["B","C"]
#MM_data=["B","C"]
#EE_data=["DoubleEG_Run2016B-03Feb2017_ver2-v2_nTuple.root","DoubleEG_Run2016C-03Feb2017-v1_nTuple.root"]
EE_data=["DoubleEG_Run2016C-03Feb2017-v1_nTuple.root"]

def readTree( filename, treename ):
    """
    filename: name of file containing the tree
    treename: name of the tree
    returns: TChain Object
    """
    tree = ROOT.TChain( treename )
    tree.AddFile( filename )
    return tree

def checkDoubleEleTrig(event):
    boolDoubleEleTrig = (
        event.HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v 
        or event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v 
        or event.HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v 
        or event.HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v
        )
    return boolDoubleEleTrig
    
def checkDoubleMuTrig(event):    
    boolDoubleMutrig = (
        event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v 
        or event.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v 
        or event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v 
        or event.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v 
        or event.HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v 
        or event.HLT_Mu27_TkMu8_v 
        or event.HLT_Mu30_TkMu11_v)
    return boolDoubleMutrig


def checkMuEleTrig(event):    
    boolDoubleMutrig = (
        event.HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v  
        or event.HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v  
        or event.HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v  
        or event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v  
        or event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v  
        or event.HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v  
        or event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v  
        or event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v  
        or event.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v  
        or event.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v  
        or event.HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v  
        or event.HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v)
    return boolDoubleMutrig
    


def searchFile(filename, treename):
    print "Now investigating ",filename.split('/')[-1]
    tree = readTree(filename,treename)
    countEE=0
    countEE_woMM=0
    countEE_woEM=0
    countEE_woBoth=0
    countMM=0
    countEM=0    
    length=tree.GetEntries()
    count=0
    for evt in tree:
        if(count % 100000 == 0):
            print count,"of",length, "|",float(count)/float(length)*100.,"%"
        count+=1
        if checkDoubleEleTrig(evt):
            countEE+=1
            if checkDoubleMuTrig(evt):
                countMM+=1
            else:
                countEE_woMM+=1
            if checkMuEleTrig(evt):
                countEM+=1
            else:
                countEE_woEM+=1
            if (not checkDoubleMuTrig(evt) and not checkMuEleTrig(evt)):
                countEE_woBoth+=1
            
    print "Has %d EE events"%(countEE)
    print "Has %d MM+EE events | %.2f pc"%(countMM,float(countMM)/float(countEE)*100.)
    print "Has %d EM+EE events | %.2f pc"%(countEM,float(countEM)/float(countEE)*100.)
    print "Has %d EE w/o MM events | %.2f pc"%(countEE_woMM,float(countEE_woMM)/float(countEE)*100.)
    print "Has %d EE w/o EM events | %.2f pc"%(countEE_woEM,float(countEE_woEM)/float(countEE)*100.)
    print "Has %d EE w/o EM+MM events| %.2f pc "%(countEE_woBoth,float(countEE_woBoth)/float(countEE)*100.)
    
    
if __name__=="__main__":
    #filenames=sys.argv[1:]
    #thefile = open('.txt', 'w')
    for item in EE_data:
        searchFile(tree_path+item,"TreeWriter/eventTree")
