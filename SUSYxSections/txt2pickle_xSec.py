#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import argparse
import re
import pickle

def getXsecsStrong(filename):
    xSecs = {}
    with open( args.filename ) as f:
        for line in f.readlines():
            match = re.match("\s*(\d+)\s*GeV(.*)±(.*)%.*", line )
            if match:
                m_str, xSec_str, uncert_str = match.groups()
                xSecs[ int(m_str) ] = ( float(xSec_str), float( uncert_str ) )
    return xSecs

def getXsecsWeak(filename):
    xSecs = {}
    with open( args.filename ) as f:
        for line in f.readlines():
            match = re.match("\s*(\d+)\s*([^\s]*)\s*([^\s]*)\s*", line )
            if match:
                m_str, xSec_str, uncert_str = match.groups()
                xSecs[ int(m_str) ] = ( float(xSec_str)/1000, float( uncert_str ) / float( xSec_str ) )
    return xSecs



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    args = parser.parse_args()

    filename = args.filename
    if "Gluino" in filename or "Squark" in filename:
        xSecs = getXsecsStrong(filename)
    elif "N2C1" in filename or "C1C1" in filename:
        xSecs = getXsecsWeak(filename)

    output = open( args.filename.replace("txt", "pkl" ), 'wb')
    pickle.dump(xSecs, output)
    output.close()
