#!/bin/bash

version="v34"
today=`date +"0%Y-%m-%d"`
#set -x

inputversion="v16_CiC"
inputfolder="/afs/cern.ch/work/o/obondu/public/forRadion/plotTrees/v16/2014-10-08_selection_noRegression_noMassCut_v16_CiC"

# IMPORTANT NOTES:
# FOR NOW THE DEFAULT IS NO REGRESSION
# As of today (Oct. 13) regression will be produced as well by default for low mass when trees are available
 

# WHICH ANALYSIS TO PROCESS
doNonResonant=1
doResonantLowMass=0
doResonantHighMass=0
# WHICH SAMPLES TO PROCESS (default is also running dataCS and diphoton-sherpa, minimum is data + signal)
doTheStrictMinimum=0


# OTHER GLOBAL SETTINGS, IN MOST USE CASE YOU SHOULD NOT TOUCH THIS
cutLevel=0
massCutVersion=4 # From Summer 14 cut update
controlSampleWeights="scales_2D_pt_data_4GeVbinning.root"
applyFTR14001=1

# INITIALIZATION OF THE PROCESSING LIST
i=-1

if [ ${doResonantHighMass} == 1 ]
then
    ##### PREPARE MGGJJ-FIT TREES
    kinfitlabel[0]="noKinFit"
    kinfitlabel[1]="withKinFit"
    kinfitjet[0]="base"
    kinfitjet[1]="kin"
    for ikin in `seq 0 1`
    do
        outfolder="${version}_fitToMggjj_${kinfitlabel[${ikin}]}"
        mkdir -p ${outfolder}
        samplelist="Radion Graviton Data"
        if [ ${doTheStrictMinimum} == 0 ]
        then
            samplelist="${samplelist} DataCS diphojet_sherpa_8TeV"
        fi
        for sample in `echo "${samplelist}"`
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
#            extraline="--applyMjjCut 0  --applyMggCut 0 --applyMtotCut 0"
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
                line[${i}]="${line[${i}]} --cutLevel ${cutLevel}"
                line[${i}]="${line[${i}]} --removeUndefinedBtagSF ${removeUndefinedBtagSF}"
                line[${i}]="${line[${i}]} --massCutVersion ${massCutVersion}"
                line[${i}]="${line[${i}]} --applyPhotonIDControlSample ${applyPhotonIDControlSample}"
                line[${i}]="${line[${i}]} --controlSampleWeights ${controlSampleWeights}"
                line[${i}]="${line[${i}]} ${extraline}"
                log[${i}]="${outfolder}/${outtree}_m${mass}.eo"
            #    echo -e "i= ${i}\tline= ${line[${i}]}"
            done
        done
    done
fi
##### END OF MGGJJ-FIT TREES
    
if [ ${doResonantLowMass} == 1 ]
then
    ##### PREPARE LOW-MASS RESONANCE MGG-FIT TREES
    kinfitlabel[0]="noKinFit"
    kinfitlabel[1]="withKinFit"
    kinfitjet[0]="base"
    kinfitjet[1]="kin"
    #kinfitlabel[2]="withReg"
    #kinfitlabel[3]="withRegKinFit"
    #kinfitjet[2]="reg"
    #kinfitjet[3]="regkin"
    for ikin in `seq 0 1`
    do
        outfolder="${version}_fitToMgg_${kinfitlabel[${ikin}]}"
        mkdir -p ${outfolder}
        samplelist="Radion Graviton MSSM ggh_m125_powheg_8TeV vbf_m125_8TeV wzh_m125_8TeV_wh wzh_m125_8TeV_zh tth_m125_8TeV bbh_m125_8TeV Data"
        if [ ${doTheStrictMinimum} == 0 ]
        then
            samplelist="${samplelist} DataCS diphojet_sherpa_8TeV"
        fi
        for sample in `echo "${samplelist}"`
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
                elif [ "${sample}" == "DataCS" ]
                then
                    itype="0"
                    applyPhotonIDControlSample=1
                    intree="Data"
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
                line[${i}]="${line[${i}]} --fitStrategy mgg"
                line[${i}]="${line[${i}]} --cutLevel ${cutLevel}"
                line[${i}]="${line[${i}]} --mass ${mass}"
                line[${i}]="${line[${i}]} --removeUndefinedBtagSF ${removeUndefinedBtagSF}"
                line[${i}]="${line[${i}]} --massCutVersion ${massCutVersion}"
                line[${i}]="${line[${i}]} --applyPhotonIDControlSample ${applyPhotonIDControlSample}"
                line[${i}]="${line[${i}]} --controlSampleWeights ${controlSampleWeights}"
                line[${i}]="${line[${i}]} ${extraline}"
                log[${i}]="${outfolder}/${outtree}_m${mass}.eo"
            #    echo -e "i= ${i}\tline= ${line[${i}]}"
            done
        done
    done # end of loop on kinfit scenarii
fi
##### END OF LOW-MASS RESONANCE TREES

#####
if [ ${doNonResonant} == 1 ]
then
    ##### PREPARE NON-RESONANT MGG-FIT TREES
    kinfitlabel[0]="noKinFit"
    kinfitlabel[1]="withKinFit"
    kinfitjet[0]="base"
    kinfitjet[1]="kin"
    for ikin in `seq 0 1`
    do
        outfolder="${version}_fitToMgg_${kinfitlabel[${ikin}]}"
        mkdir -p ${outfolder}
        samplelist="ggHH_8TeV"
        if [ ${doTheStrictMinimum} == 0 ]
        then
            samplelist="${samplelist} DataCS"
        fi
        for sample in `echo "${samplelist}"`
        do
            strategylist="mgg 2D"
            if [ ${applyFTR14001} == 1 ]
            then
                strategylist="FTR14001"
            fi
            for fitStrategy in `echo "${strategylist}"`
            do
                mass=0
                intree=${sample}
                outtree=${sample}
                itype="1"
                removeUndefinedBtagSF=0
                applyPhotonIDControlSample=0
                suffix=""
                extraline=""
                if [ "${sample}" == "Data" ]
                then
                    itype="0"
                elif [ "${sample}" == "DataCS" ]
                then
                    itype="0"
                    applyPhotonIDControlSample=1
                    intree="Data"
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
                line[${i}]="${line[${i}]} --fitStrategy ${fitStrategy}"
                line[${i}]="${line[${i}]} --cutLevel ${cutLevel}"
                line[${i}]="${line[${i}]} --mass 0"
                line[${i}]="${line[${i}]} --removeUndefinedBtagSF ${removeUndefinedBtagSF}"
                line[${i}]="${line[${i}]} --massCutVersion ${massCutVersion}"
                line[${i}]="${line[${i}]} --applyPhotonIDControlSample ${applyPhotonIDControlSample}"
                line[${i}]="${line[${i}]} --controlSampleWeights ${controlSampleWeights}"
                line[${i}]="${line[${i}]} ${extraline}"
                log[${i}]="${outfolder}/${outtree}_m${mass}.eo"
            #    echo -e "i= ${i}\tline= ${line[${i}]}"
            done
        done
    done # end of loop on kinfit scenarii
fi


#### PRODUCE EVERYTHING    
itot=${i}
for iline in `seq 0 ${itot}`
do
    echo -e "iline= ${iline} / ${itot}\t\t${line[${iline}]}"
    ./quickTrees.exe ${line[${iline}]} &> ${log[${iline}]}
done
