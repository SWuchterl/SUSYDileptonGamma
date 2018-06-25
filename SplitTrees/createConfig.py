#!/usr/bin/env python2
import ConfigParser

generalSettings = {
    "eventPercentage" : 100.,
    #"eventPercentage" : 50.,
    #"eventPercentage" : 10.,
    #"eventPercentage" : 5.,
    "veto" : 6,      # 5=PromptFinalState, 6=isPromptFinalState, 7= None,
    "ht": 1
}

#outputFolder = "standard" 
#outputFolder = "medium" 
outputFolder = "mediumPOG" 


#outputPath="/user/swuchterl/tempTrees"
outputPath="/net/data_cms1b/user/swuchterl/tempTrees"

outputFolders = {
    "standard" : outputPath+"/tightMVAID",
    "medium": outputPath+"/mediumID",
    "mediumPOG": outputPath+"/mediumIDPOG"
}




config = ConfigParser.ConfigParser()

config.add_section("generalSettings")
for key in generalSettings:
    config.set("generalSettings",key,str(generalSettings[key]))




config.set("generalSettings","outputFolder",str(outputFolders[outputFolder]))




with open('example.ini', 'wb') as configfile:
    config.write(configfile)

print "Generated example.ini config file."
