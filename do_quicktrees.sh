#!/bin/bash
version="v34"
today=`date +"0%Y-%m-%d"`
#set -x
inputversion="v16_CiC"
inputfolder="/afs/cern.ch/work/o/obondu/public/forRadion/plotTrees/v16/2014-10-08_selection_noRegression_noMassCut_v16_CiC"
i=-1
controlSampleWeights="scales_2D_pt_data_4GeVbinning.root"
##### PREPARE MGGJJ-FIT TREES
kinfitlabel[0]="noKinFit"
kinfitlabel[1]="withKinFit"
kinfitjet[0]="base"
kinfitjet[1]="kin"
for ikin in `seq 0 1`
do
outfolder="${version}_fitToMggjj_${kinfitlabel[${ikin}]}_CIC_LAST"
mkdir -p ${outfolder}
#for sample in `echo "Radion Graviton MSSM ggh_m125_powheg_8TeV vbf_m125_8TeV wzh_m125_8TeV_wh wzh_m125_8TeV_zh tth_m125_8TeV bbh_m125_8TeV Data DataCS diphojet_sherpa_8TeV ggHH_8TeV qcd_30_8TeV_ff qcd_40_8TeV_ff qcd_30_8TeV_pf qcd_40_8TeV_pf gjet_20_8TeV_pf gjet_40_8TeV_pf DYJetsToLL LNuGG_FSR_8TeV LNuGG_ISR_8TeV ttGG_8TeV tGG_8TeV TTGJets_8TeV ZGToLLG_8TeV"`
for sample in `echo ""`
do
for mass in `echo "400 450 500 550 600 650 700 800 900 1000 1100"`
do
intree=${sample}
outtree=${sample}
itype="1"
removeUndefinedBtagSF=0
applyPhotonIDControlSample=0
suffix=""
extraline=""
# extraline="--applyMjjCut 0 --applyMggCut 0 --applyMtotCut 0"
if [ "${sample}" == "Radion" ]
then
intree="${sample}_m${mass}_8TeV"
outtree="${sample}_m${mass}_8TeV"
itype="-${mass}"
removeUndefinedBtagSF=0
elif [ "${sample}" == "Graviton" ]
then
if [ "${mass}" == "500" ] || [ "${mass}" == "700" ] || [ "${mass}" == "1000" ]
then
intree="${sample}_m${mass}_8TeV"
outtree="${sample}_m${mass}_8TeV"
itype="-${mass}"
removeUndefinedBtagSF=0
else
continue
fi
elif [ "${sample}" == "Data" ]
then
itype="0"
elif [ "${sample}" == "DataCS" ]
then
itype="0"
intree="Data"
applyPhotonIDControlSample=1
suffix="controlSample_"
fi
i=$((${i} + 1))
line[${i}]=""
line[${i}]="${line[${i}]} --inputfile ${inputfolder}/${intree}_noRegression_noMassCut_${suffix}${inputversion}.root"
line[${i}]="${line[${i}]} --inputtree ${intree}"
line[${i}]="${line[${i}]} --outputtree TCVARS"
line[${i}]="${line[${i}]} --outputfile ${outfolder}/${outtree}_m${mass}.root"
line[${i}]="${line[${i}]} --type ${itype}"
line[${i}]="${line[${i}]} --whichJet ${kinfitjet[${ikin}]}"
line[${i}]="${line[${i}]} --fitStrategy mggjj"
line[${i}]="${line[${i}]} --cutLevel 0"
line[${i}]="${line[${i}]} --removeUndefinedBtagSF ${removeUndefinedBtagSF}"
line[${i}]="${line[${i}]} --massCutVersion 4"
line[${i}]="${line[${i}]} --applyPhotonIDControlSample ${applyPhotonIDControlSample}"
line[${i}]="${line[${i}]} --controlSampleWeights ${controlSampleWeights}"
line[${i}]="${line[${i}]} ${extraline}"
log[${i}]="${outfolder}/${outtree}_m${mass}.eo"
# echo -e "i= ${i}\tline= ${line[${i}]}"
done
done
done
##### PREPARE MGG-FIT TREES
kinfitlabel[0]="noKinFit"
kinfitlabel[1]="withKinFit"
kinfitjet[0]="base"
kinfitjet[1]="kin"
for ikin in `seq 0 1`
do
outfolder="${version}_fitToMgg_${kinfitlabel[${ikin}]}_CIC_LAST"
mkdir -p ${outfolder}
#for sample in `echo "Radion Graviton MSSM ggh_m125_powheg_8TeV vbf_m125_8TeV wzh_m125_8TeV_wh wzh_m125_8TeV_zh tth_m125_8TeV bbh_m125_8TeV Data DataCS diphojet_sherpa_8TeV ggHH_8TeV qcd_30_8TeV_ff qcd_40_8TeV_ff qcd_30_8TeV_pf qcd_40_8TeV_pf gjet_20_8TeV_pf gjet_40_8TeV_pf DYJetsToLL LNuGG_FSR_8TeV LNuGG_ISR_8TeV ttGG_8TeV tGG_8TeV TTGJets_8TeV ZGToLLG_8TeV"`
#for sample in `echo "ggHH_8TeV ggh_m125_powheg_8TeV vbf_m125_8TeV wzh_m125_8TeV_wh wzh_m125_8TeV_zh tth_m125_8TeV bbh_m125_8TeV diphojet_sherpa_8TeV qcd_30_8TeV_ff qcd_40_8TeV_ff qcd_30_8TeV_pf qcd_40_8TeV_pf gjet_20_8TeV_pf gjet_40_8TeV_pf DYJetsToLL LNuGG_FSR_8TeV LNuGG_ISR_8TeV ttGG_8TeV tGG_8TeV TTGJets_8TeV ZGToLLG_8TeV"`
for sample in `echo "Data"`
do
for mass in `echo "260 270 300 350 400 450 500"`
do
intree=${sample}
outtree=${sample}
itype="1"
removeUndefinedBtagSF=0
applyPhotonIDControlSample=0
suffix=""
extraline=""
if [ "${sample}" == "Radion" ]
then
if [ "${mass}" == "260" ]
then
continue
else
intree="${sample}_m${mass}_8TeV"
outtree="${sample}_m${mass}_8TeV"
itype="-${mass}"
removeUndefinedBtagSF=0
fi
elif [ "${sample}" == "Graviton" ]
then
if [ "${mass}" == "260" ] || [ "${mass}" == "270" ] || [ "${mass}" == "350" ] || [ "${mass}" == "400" ] || [ "${mass}" == "450" ]
then
continue
else
intree="${sample}_m${mass}_8TeV"
outtree="${sample}_m${mass}_8TeV"
itype="-${mass}"
removeUndefinedBtagSF=0
fi
elif [ "${sample}" == "MSSM" ]
then
if [ "${mass}" == "270" ] || [ "${mass}" == "400" ] || [ "${mass}" == "450" ] || [ "${mass}" == "500" ]
then
continue
else
intree="${sample}_m${mass}_8TeV"
outtree="${sample}_m${mass}_8TeV"
itype="-${mass}"
removeUndefinedBtagSF=0
fi
elif [ "${sample}" == "Data" ]
then
itype="0"
# extraline="--applyMtotCut 0"
elif [ "${sample}" == "DataCS" ]
then
itype="0"
applyPhotonIDControlSample=1
intree="Data"
suffix="controlSample_"
elif [ "${sample}" == "ggHH_8TeV" ]
then
itype="-2"
extraline="--applyMtotCut 0"
fi
i=$((${i} + 1))
line[${i}]=""
line[${i}]="${line[${i}]} --inputfile ${inputfolder}/${intree}_noRegression_noMassCut_${suffix}${inputversion}.root"
line[${i}]="${line[${i}]} --inputtree ${intree}"
line[${i}]="${line[${i}]} --outputtree TCVARS"
line[${i}]="${line[${i}]} --outputfile ${outfolder}/${outtree}_m${mass}.root"
line[${i}]="${line[${i}]} --type ${itype}"
line[${i}]="${line[${i}]} --whichJet ${kinfitjet[${ikin}]}"
line[${i}]="${line[${i}]} --fitStrategy mgg"
line[${i}]="${line[${i}]} --cutLevel 0"
line[${i}]="${line[${i}]} --mass ${mass}"
line[${i}]="${line[${i}]} --removeUndefinedBtagSF ${removeUndefinedBtagSF}"
line[${i}]="${line[${i}]} --massCutVersion 4"
line[${i}]="${line[${i}]} --applyPhotonIDControlSample ${applyPhotonIDControlSample}"
line[${i}]="${line[${i}]} --controlSampleWeights ${controlSampleWeights}"
line[${i}]="${line[${i}]} ${extraline}"
log[${i}]="${outfolder}/${outtree}_m${mass}.eo"
# echo -e "i= ${i}\tline= ${line[${i}]}"
done
done
done # end of loop on kinfit scenarii
#### PRODUCE EVERYTHING
itot=${i}
for iline in `seq 0 ${itot}`
do
echo -e "iline= ${iline} / ${itot}\t\t${line[${iline}]}"
./quickTrees.exe ${line[${iline}]} &> ${log[${iline}]}
done
