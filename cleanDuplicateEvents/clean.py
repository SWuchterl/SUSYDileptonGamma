#!/usr/bin/env python2

import ROOT
ROOT.gErrorIgnoreLevel=ROOT.kError
from ROOT import *

#nTuple trees to clean:
tree_path="../scratch/v01/"
#EM_data=["B","C"]
#MM_data=["B","C"]
EE_data=["DoubleEG_Run2016B-03Feb2017_ver2-v2_nTuple.root","DoubleEG_Run2016C-03Feb2017-v1_nTuple.root"]

class Event():
    def __init__(self, runNo, evtNo, lumNo):
        self.runNo=runNo
        self.evtNo=evtNo
        self.lumNo=lumNo
    def __eq__(self,other):
        if isinstance(other,self.__class__):
            return self.runNo == other.runNo and self.evtNo == other.evtNo and self.lumNo == other.lumNo
        else:
            return False
    def __ne__(self,other):
        return not self.__eq__(other)

    def __hash__(self):
        return self.evtNo

    def __str__(self):
        return "(r%d, e%d, l%d)"%(self.runNo,self.evtNo,self.lumNo)
        

def readTree( filename, treename ):
    """
    filename: name of file containing the tree
    treename: name of the tree
    returns: TChain Object
    """
    tree = ROOT.TChain( treename )
    tree.AddFile( filename )
    return tree

def filterFiles(filename, filename2, treename,treename2):
    print "Now investigating ",filename.split('/')[-1]," vs ",filename2.split('/')[-1]
    duplicates=set()
    evtList1=set()
    tree = readTree(filename,treename)
    tree2 = readTree(filename2,treename2)
    length1=tree.GetEntries()
    length2=tree2.GetEntries()
    iDoubleEvents=0
    count1=0
    count2=0
    for evt in tree:
        event = Event(evt.runNo, evt.evtNo, evt.lumNo)
        evtList1.add(event)
        if (count1 % 100000==0):
            print count1,"of",length1, "|",float(count1)/float(length1)*100.,"%"
        count1+=1
    for evt2 in tree2:
        if (count2 % 100000==0):
            print count2,"of",length2, "|",float(count2)/float(length2)*100.,"%"
        count2+=1
        Event2=Event(evt2.runNo,evt2.evtNo,evt2.lumNo)
        if(Event2 in evtList1):
            duplicates.add(Event2)
            iDoubleEvents+=1
    if iDoubleEvents>0:
        print "Has %d/%d duplicated events"%(iDoubleEvents,tree.GetEntries())
    else:
        print "Has NO duplicated events"
    return duplicates

if __name__=="__main__":
    #filenames=sys.argv[1:]
    listEE=[]
    listMM=[]
    listEM=[]
    listEE.append(filterFiles(tree_path+EE_data[0],tree_path+EE_data[1],"TreeWriter/eventTree","TreeWriter/eventTree"))
    #print listEE
    thefile = open('test.txt', 'w')
    for item in listEE:
        print>>thefile, item
    #for filename in EE_data:
        #listEE.append(filterFile(filename,"TreeWriter/eventTree"))
    #for filename in MM_data:
        #listMM.append(filterFile(filename,"TreeWriter/eventTree"))
    #for filename in EM_data:
        #listEM.append(filterFile(filename,"TreeWriter/eventTree"))
