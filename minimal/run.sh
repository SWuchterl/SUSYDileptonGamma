#!/bin/bash -e
export SCRAM_ARCH=slc6_amd64_gcc530
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh



thispath=`pwd`
cd ~/CMSSW/treewriter/v22/CMSSW_8_0_26_patch1/src
eval `scramv1 runtime -sh`
cd $thispath

./run.py "$@"
