cd /afs/cern.ch/work/e/egkougko/private
echo "I am here for the asetup: "
pwd
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh; sleep 1
source ${ATLAS_LOCAL_ROOT_BASE}/packageSetups/atlasLocalROOTSetup.sh --rootVer=6.18.04-x86_64-centos7-gcc8-opt
echo "ROOT setup OK"
cd /afs/cern.ch/work/e/egkougko/private/LGADUtils/macros
echo "I am here for to run the macro: "
pwd
root -l -b Run.C+
echo "Done!!!"