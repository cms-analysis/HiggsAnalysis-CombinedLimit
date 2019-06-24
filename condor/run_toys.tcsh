#!/bin/tcsh

set inputRoot2016 = $1
set inputRoot2017 = $2
set signalType = $3
set mass = $4
set year = $5
set dataType = $6
set rVal = $7
set seed = $8
set numToys = $9
set iterations = $10
set doToyS = $11

set base_dir = `pwd`
printf "\n\n base dir is $base_dir\n\n"

source /cvmfs/cms.cern.ch/cmsset_default.csh
setenv SCRAM_ARCH slc6_amd64_gcc530

printf "\n\n ls output\n"
ls -l

printf "\n\n Get the code needed .\n\n"
cmsrel CMSSW_8_1_0
cd CMSSW_8_1_0/src
eval `scramv1 runtime -csh`
git clone https://github.com/StealthStop/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
scram b clean
scram b -j8

printf "\n\n ls output\n"
ls -l

printf "\n\n output of uname -s : "
uname -s
printf "\n\n"

setenv LD_LIBRARY_PATH ${PWD}:${LD_LIBRARY_PATH}
printf "\n\n LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}\n\n"

printf "\n\n ls output\n"
ls -l

printf "\n\n Attempting to run Fit executable.\n\n"
mkdir ${inputRoot2016}
xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/FitInputs/${inputRoot2016}/njets_for_Aron.root     ${inputRoot2016}/.
xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/FitInputs/${inputRoot2016}/ttbar_systematics.root  ${inputRoot2016}/.

mkdir ${inputRoot2017}
xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/FitInputs/${inputRoot2017}/njets_for_Aron.root     ${inputRoot2017}/.
xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/FitInputs/${inputRoot2017}/ttbar_systematics.root  ${inputRoot2017}/.

eval `scramv1 runtime -csh`

combineCards.py Y16=Card2016.txt Y17=Card2017.txt > CardCombo.txt
root -l -q -b 'make_MVA_8bin_ws.C("2016","'${inputRoot2016}'","'${signalType}'","'${mass}'","'${dataType}'")'
root -l -q -b 'make_MVA_8bin_ws.C("2017","'${inputRoot2017}'","'${signalType}'","'${mass}'","'${dataType}'")'
text2workspace.py Card${year}.txt -o ws_${year}_${signalType}_${mass}.root -m ${mass} --keyword-value MODEL=${signalType}

if ($doToyS == 1) then
    printf "\n\n Running sig. toys\n"
    combine -M HybridNew --LHCmode LHC-significance ws_${year}_${signalType}_${mass}.root -m ${mass} --keyword-value MODEL=${signalType} -n ${year} --saveToys --saveHybridResult -T ${numToys} -s ${seed} --fullBToys -i ${iterations}
else
    printf "\n\n Running limit toys\n"
    combine -M HybridNew --LHCmode LHC-limits       ws_${year}_${signalType}_${mass}.root -m ${mass} --keyword-value MODEL=${signalType} -n ${year} --saveToys --saveHybridResult -T ${numToys} -s ${seed} --fullBToys --singlePoint ${rVal} --clsAcc 0 
endif

printf "\n\n ls output\n"
ls -l

mv *.root ${base_dir}
mv log*.txt ${base_dir}

cd ${base_dir}

printf "\n\n ls output\n"
ls -l
