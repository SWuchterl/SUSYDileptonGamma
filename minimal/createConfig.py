#!/usr/bin/env python2
import ConfigParser

generalSettings = {
    #"eventPercentage" : 0.1,
    "eventPercentage" : 100.,
    "doSignalScan" : 1,
    "doSignalScanTChiNGSplit" : 1,
    #"veto" : 6,      # 5=PromptFinalState, 6=isPromptFinalState, 7= None
    "veto" : 7,      # 5=PromptFinalState, 6=isPromptFinalState, 7= None
    "doCutFlow" : 1,
    "doCutFlowFine" : 1
}
#outputFolder = "standard" ## old variant  before v12
#outputFolder = "danilo"
outputFolder = "uncorrected" ## will now be without veto
#outputFolder = "EGRegression"
#outputFolder = "noVeto" ## for now right folder
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
    "noVeto" : "_noVeto",
    "uncorrected" : "_uncorrected",
    "danilo" : "_danilo",
    "EGRegression" : "_EGRegression",
    "signalScan" : "_signalScan"
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
