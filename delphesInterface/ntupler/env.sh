OLDDIR=`pwd`

#DAnalysis
DANALYSISPATH=/afs/cern.ch/work/c/calderon/public/YR_Delphes/DAnalysis
#DELPHES_PATH=/afs/cern.ch/work/c/calderon/private/delphes
#CMSSW_PATH=/afs/cern.ch/work/c/calderon/private/CMSSW_9_1_0_pre3

#CMSSW and Delphes
cd $DANALYSISPATH/../
cd CMSSW_9_1_0_pre3
eval `scramv1 runtime -sh`
cd ..
cd delphes
source DelphesEnv.sh
export DELPHES_PATH=`pwd`
cd $OLDDIR

#export DELPHES_PATH=$DELPHES_PATH
export PYTHIA8=$CMSSW_RELEASE_BASE/../../../external/pythia8/212-ikhhed3
export LD_LIBRARY_PATH=$PYTHIA8/lib:$LD_LIBRARY_PATH
export DANALYSISPATH=$DANALYSISPATH
export LD_LIBRARY_PATH=$DANALYSISPATH:$LD_LIBRARY_PATH
export PATH=$PATH:$DANALYSISPATH
export LD_LIBRARY_PATH=/afs/cern.ch/work/c/calderon/public/YR_Delphes/DAnalysis:$LD_LIBRARY_PATH
export PATH=/afs/cern.ch/work/c/calderon/public/YR_Delphes/DAnalysis:$PATH
