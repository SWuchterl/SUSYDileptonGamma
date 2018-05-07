#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import argparse
import re
import pickle
import numpy as np

#def getXsecsStrong(filename):
    #xSecs = {}
    #with open( args.filename ) as f:
        #for line in f.readlines():
            #match = re.match("\s*(\d+)\s*GeV(.*)±(.*)%.*", line )
            #if match:
                #m_str, xSec_str, uncert_str = match.groups()
                #xSecs[ int(m_str) ] = ( float(xSec_str), float( uncert_str ) )
    #return xSecs

def getXsecsWeak(filename):
    xSecs = {}
    #with open( args.filename ) as f:
    with open( filename ) as f:
        for line in f.readlines():
            match = re.match("\s*(\d+)\s*([^\s]*)\s*([^\s]*)\s*", line )
            if match:
                m_str, xSec_str, uncert_str = match.groups()
                xSecs[ int(m_str) ] = ( float(xSec_str)/1000, float( uncert_str ) / float( xSec_str ) )
    return xSecs
    
    
if __name__ == "__main__":
    #parser = argparse.ArgumentParser()
    #parser.add_argument("filename")
    #args = parser.parse_args()

    filename1 = "xSec_SMS_C1C1_13TeV.txt"
    filename2 = "xSec_SMS_N2C1_13TeV.txt"
    filename3 = "xSec_SMS_TChiNG_13TeV.txt"

    xSec1 = getXsecsWeak(filename1)
    xSec2 = getXsecsWeak(filename2)
    
    #print xSec1
    #print xSec2
    
    newxSec={}

    for key in xSec1:
        newxSec[key] = (xSec1[key][0]+xSec2[key][0],np.sqrt(xSec1[key][1]**2.+xSec2[key][1]**2.))
        #newxSec[key] = np.sqrt(xSec1[key][1]**2.+xSec2[key][1]**2.)
        
    #print newxSec

    #filename = args.filename
    #if "Gluino" in filename or "Squark" in filename:
        #xSecs = getXsecsStrong(filename)
    #elif "N2C1" in filename or "C1C1" in filename:
        #xSecs = getXsecsWeak(filename)
#
    output = open(filename3.replace("txt", "pkl" ), 'wb')
    pickle.dump(newxSec, output)
    output.close()
