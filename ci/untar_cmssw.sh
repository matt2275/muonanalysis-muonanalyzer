tar xzvf build.tar.gz
cd CMSSW_*/src
scram b projectrename
eval $(scram runtime -sh)
export CMS_PATH=/cvmfs/cms-ib.cern.ch/week0
CMSSW_VERSION_SPLIT=($(echo $CMSSW_VERSION | tr _ \ ))
export CMSSW_MAJOR=${CMSSW_VERSION_SPLIT[1]}
export CMSSW_MINOR=${CMSSW_VERSION_SPLIT[2]}
