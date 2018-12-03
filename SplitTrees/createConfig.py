#!/usr/bin/env python2
import ConfigParser

generalSettings = {
    "eventPercentage": 100.,
     #"eventPercentage" : 1.,
    # "eventPercentage" : 50.,
    # "eventPercentage" : 10.,
    # "eventPercentage" : 5.,
    # "eventPercentage" : 0.5,
    "veto": 6,      # 5=PromptFinalState, 6=isPromptFinalState, 7= None,
    "ht": 1,
    "met": 0,
    "htpure": 0
}

#outputFolder = "standard"
#outputFolder = "medium"
# outputFolder = "mediumPOG"
#outputFolder = "mediumPOG_noISRtop"
#outputFolder = "mediumIDPOG_noTopPt_noNIsr"
#outputFolder = "mediumIDPOG_TopPt_noNIsr"
#outputFolder = "mediumIDPOG_TopPt_NIsr"
#outputFolder = "AN"
outputFolder = "addtrigger"
# outputFolder = "LesyaCheck"
#outputFolder = "AN_rishi"


# outputPath="/user/swuchterl/tempTrees"
outputPath = "/net/data_cms1b/user/swuchterl/tempTrees"

outputFolders = {
    "standard": outputPath + "/tightMVAID",
    "AN": outputPath + "/AN",
    "AN_rishi": outputPath + "/AN_rishi",
    "medium": outputPath + "/mediumID",
    "mediumPOG": outputPath + "/mediumIDPOG",
    "mediumIDPOG_noTopPt_noNIsr": outputPath + "/mediumIDPOG_noTopPt_noNIsr",
    "mediumIDPOG_TopPt_noNIsr": outputPath + "/mediumIDPOG_TopPt_noNIsr",
    "mediumIDPOG_TopPt_NIsr": outputPath + "/mediumIDPOG_TopPt_NIsr",
    "addtrigger": outputPath + "/addtrigger",
    "LesyaCheck": outputPath + "/LesyaCheck",
}


config = ConfigParser.ConfigParser()

config.add_section("generalSettings")
for key in generalSettings:
    config.set("generalSettings", key, str(generalSettings[key]))


config.set("generalSettings", "outputFolder", str(outputFolders[outputFolder]))


with open('example.ini', 'wb') as configfile:
    config.write(configfile)

print "Generated example.ini config file."
