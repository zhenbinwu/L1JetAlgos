#!/bin/csh -v

set SCRAM = DELSCR
set CMSSW = DELDIR
set EXE   = DELEXE
set OUTPUT = OUTDIR

#============================================================================#
#-----------------------------   Setup the env   ----------------------------#
#============================================================================#
echo "============  Running on" $HOST " ============"
cd ${_CONDOR_SCRATCH_DIR}
source /cvmfs/cms.cern.ch/cmsset_default.csh
setenv SCRAM_ARCH ${SCRAM}
eval `scramv1 project CMSSW ${CMSSW}`
cd ${CMSSW}
eval `scramv1 runtime -csh` # cmsenv is an alias not on the workers
echo "CMSSW: "$CMSSW_BASE

xrdcp "root://cmseos.fnal.gov/${OUTPUT}/CMSSW.tgz" .
tar -xzf CMSSW.tgz
scram b 

if ($? == 0) then
  rm CMSSW.tgz
endif

cd ${_CONDOR_SCRATCH_DIR}
foreach tarfile (`ls *gz FileList/*gz`)
  echo $tarfile
  tar -xzf $tarfile 
end

setenv PYTHONPATH  ${_CONDOR_SCRATCH_DIR}
setenv XDG_CONFIG_HOME ${_CONDOR_SCRATCH_DIR}/.config
setenv XDG_CACHE_HOME ${_CONDOR_SCRATCH_DIR}/.cache
#============================================================================#
#--------------------------   To Run the Process   --------------------------#
#============================================================================#
echo $EXE $argv
echo $PYTHONPATH
ls $PYTHONPATH
python $EXE $argv

if ($? == 0) then
  foreach tarfile (`ls *gz FileList/*gz`)
    tar -tf $tarfile  | xargs rm -r
  end
  foreach outfile (`ls *root`)
    echo "Copying ${outfile} to ${OUTPUT}"
    xrdcp $outfile "root://cmseos.fnal.gov/${OUTPUT}"
    if ($? == 0) then
      rm $outfile
    endif
  end
endif
