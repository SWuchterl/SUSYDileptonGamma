#!/usr/bin/env python2
import ConfigParser
import pickle as pkl

triggerValues=pkl.load(open("../plotter/triggerStudies/triggerFactors.pkl"))

generalSettings = {
    #"eventPercentage" : 1.,
    #"eventPercentage" : 0.1,
    "eventPercentage" : 100.,
    #"eventPercentage" : 30.,
    #"eventPercentage" : 10.,
    "doSignalScan" : 1,
    "doSignalScanTChiNGSplit" : 0,
    "veto" : 6,      # 5=PromptFinalState, 6=isPromptFinalState, 7= None
    "doCutFlow" : 1,
    "doCutFlowFine" : 1,
    "htStudies": 0,
    "htStudiesPure": 0,
    "metStudies": 0,
    "pdfSplit": 0,
    "pdfStart" :0,
    "pdfStep": 10,
    "eeWeight": triggerValues["ee"]["factor"],
    "emWeight": triggerValues["em"]["factor"],
    "mmWeight": triggerValues["mm"]["factor"],
    "eeEff": triggerValues["ee"]["data"],
    "emEff": triggerValues["em"]["data"],
    "mmEff": triggerValues["mm"]["data"]
}
#outputFolder = "standard"
outputFolder = "AN"
#outputFolder = "AN_rishi"
#outputFolder = "standard2"
#outputFolder = "noTopPt_noNIsr"
#outputFolder = "TopPt_noNIsr"
#outputFolder = "TopPt_NIsr"
#outputFolder = "medium"
#outputFolder = "mediumPOG"
#outputFolder = "overlap"
#outputFolder = "overlapNoISRtop"#here
#outputFolder = "all" #overlap removed, no  top weight, all pdf, 2d hists
#outputFolder = "overlapPrompt"
#outputFolder = "overlapNormal"


#outputFolder = "danilo"
#outputFolder = "uncorrected" 
#outputFolder = "EGRegression"
#outputFolder = "noVeto" 
#outputFolder = "signalScan"

selectionsToProcess = {
    "DILEP" : 1,
    "SEL" : 1,
    "ONZ" : 1,
    "EXO" : 0,
    "ONZG" : 0,
    "ONZMET" : 1,
    "ABOVEZG" : 0,
    "EGRegression" : 0,
    "PHOTON" : 0,
    "UNCUT" : 0,
    "CONTROL": 1,
    "VALIDATION": 1
}

outputFolders = {
    "standard" : "",
    "standard2" : "_2",
    "AN" : "_AN",
    "AN_rishi" : "_AN_rishi",
    "noVeto" : "_noVeto",
    "uncorrected" : "_uncorrected",
    "danilo" : "_danilo",
    "EGRegression" : "_EGRegression",
    "signalScan" : "_signalScan",
    "medium": "_medium",
    "mediumPOG": "_mediumPOG",
    "overlap": "_overlap",
    "overlapNoISRtop": "_overlapNoISRtop",
    "overlapPrompt": "_overlapPrompt",
    "overlapNormal": "_overlapNormal",
    "all": "_all",
    "noTopPt_noNIsr": "_noTopPt_noNIsr",
    "TopPt_noNIsr": "_TopPt_noNIsr",
    "TopPt_NIsr": "_TopPt_NIsr"
}




config = ConfigParser.ConfigParser()

config.add_section("generalSettings")
for key in generalSettings:
    config.set("generalSettings",key,str(generalSettings[key]))

config.add_section("selectionsToProcess")
for section in selectionsToProcess:
    config.set("selectionsToProcess",section,str(selectionsToProcess[section]))

config.set("generalSettings","outputFolder",str(outputFolders[outputFolder]))




with open('example.ini', 'wb') as configfile:
    config.write(configfile)

print "Generated example.ini config file."
