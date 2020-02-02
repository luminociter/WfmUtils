#!/bin/bash
cd /afs/cern.ch/work/e/egkougko/private
echo "I am here for the asetup: "
pwd
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh --quiet
localSetupROOT 6.18.04-x86_64-centos7-gcc8-opt --quiet
echo "ROOT setup OK"
echo Running on host `hostname`
cd ./LGADUtils/Macros
echo "I am here for to run the macro: "
pwd
root -l -b ./Carbon/CNM_10478_W5S1005_1e14_n/-10C/W5S1005_1e14n-10C-135V.C+
echo "Done!!!"
