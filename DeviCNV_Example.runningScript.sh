#! /bin/bash
module load python
module load R

### Calculate Read Depth for each amps #####################################################################################
run_calculateReadDepthOfAmp=true
run_mergeReadDepthStatistics=true
## Caculate read Depth ratio ###############################################################################################
run_chrXNormalizeReadsDepth=true
run_filterLowQualSample=true
run_buildLinearRegression=true
### Call CNVs ################################################################################################################
run_plottingCNV=true
run_getCNVPerSample=true
run_mergeCNVFiles=true
run_scoreCNV=true

##############################################################################################################################
txtdir=./
export inSampleInfoTxt=${txtdir}/DeviCNV_Example.sampleInfo.txt  ### A txt file that contains sample’s sex information.
export inAmpliconTxt=${txtdir}/DeviCNV_Example.probeInformation.txt  ### A txt file that contains the genomic position and pool information about target capture probes.
export inScoringThTxt=${txtdir}/DeviCNV_scoringSystemThresholds.txt  ### A txt file that contains thresholds for DeviCNV's scoring system.
export inBamDir=./ExampleBams/  ### A directory that contains input bam files.
export codedir=./Code_v1.5/  ### A directory that contains code files.

export outdir=./ExampleOutputs_v1.5/  ### A name of your project directory.
export batchTag=Example ### A name of input batch of your project.
export datType="HYB"  ### "HYB" or "PCR".
export dedupOp=true   ### "true" or "false". "true" is recommended. 
export PoolList="Pool1" ### "Pool1,Pool2,Pool3" if there are 3 pools.
export MQList="MQ0,MQ20" ### "MQ0" or "MQ20" is recommended.
export dupdelList="1.1_0.9,1.3_0.7" ### "1.2_0.8,1.3_0.7". "1.3_0.7" is recommended. "1.3" is for duplication and "0.7" is for deletion.

#### make out directory #######################################################################################################
export readDepthDir=${outdir}/01.readDepthPerAmplicon/
export readDepthStatDir=${outdir}/02.readDepthStatistics/
export norReadDepthDir=${outdir}/03.norReadDepth/
export lqSampleDir=${outdir}/04.lqSampleFilter/
export RDRatioDir=${outdir}/05.readDepthRatio/
export plotDir=${outdir}/06.CNVPlot/
export CBSDir=${outdir}/07.segmentByCBS/
export CNVDir=${outdir}/08.CNV/
if [[ ! -d ${outdir} ]] ; then
   mkdir -p ${outdir}
fi
if [[ ! -d ${readDepthDir} ]] ; then
   mkdir -p ${readDepthDir}
fi
if [[ ! -d ${readDepthStatDir} ]] ; then
   mkdir -p ${readDepthStatDir}
fi
if [[ ! -d ${norReadDepthDir} ]] ; then
   mkdir -p ${norReadDepthDir}
fi
if [[ ! -d ${lqSampleDir} ]] ; then
   mkdir -p ${lqSampleDir}
fi
if [[ ! -d ${RDRatioDir} ]] ; then
   mkdir -p ${RDRatioDir}
fi
if [[ ! -d ${plotDir} ]] ; then
   mkdir -p ${plotDir}
fi
if [[ ! -d ${CBSDir} ]] ; then
   mkdir -p ${CBSDir}
fi
if [[ ! -d ${CNVDir} ]] ; then
   mkdir -p ${CNVDir}
fi

#### run DeviCNV ####################################################################################################################
#### Example of script ##############################################################################################################
for inSample in `sed -n '2,$p' ${inSampleInfoTxt} |awk '{print $1}'` ; do
  if [[ ${run_calculateReadDepthOfAmp} == true ]] ; then
    echo STEP1_calculateReadsDepthOfAmp      
    srun -s -p all.q python ${codedir}/python.calculateReadsDepthOfAmp.py ${inSample} ${inBamDir} ${inAmpliconTxt} ${readDepthDir} ${readDepthStatDir} ${dedupOp} ${datType} ${MQList}
    sleep 1
  fi
done

if [[ ${run_mergeReadDepthStatistics} == true ]] ; then
    echo STEP2_mergeReadDepthStatistics
    srun -s -p all.q python ${codedir}/python.mergeReadDepthStatistics.py ${batchTag} ${inSampleInfoTxt} ${readDepthStatDir}
    sleep 1
fi  
if [[ ${run_chrXNormalizeReadsDepth} == true ]] ; then
  echo STEP3_chrXNormalizeReadsDepth
  srun -s -p all.q python ${codedir}/python.chrXNormalizeReadDepth.py ${batchTag} ${inSampleInfoTxt} ${readDepthDir} ${norReadDepthDir} ${PoolList} ${MQList} 
  sleep 1
fi  
if [[ ${run_filterLowQualSample} == true ]] ; then
  echo STEP4_filterLowQualSample
  srun -s -p all.q Rscript ${codedir}/r.filterLowQualSample.r ${batchTag} ${inSampleInfoTxt} ${norReadDepthDir} ${lqSampleDir} ${MQList}
  sleep 1
fi   
if [[ ${run_buildLinearRegression} == true ]] ; then
  echo STEP5_buildLinearRegression      
  srun -s -p all.q python ${codedir}/python.buildLinearRegression.py ${batchTag} ${readDepthStatDir} ${norReadDepthDir} ${lqSampleDir} ${RDRatioDir} ${MQList} ${dupdelList} 
  sleep 1
fi

for inSample in `sed -n '2,$p' ${inSampleInfoTxt} |awk '{print $1}'` ; do
  if [[ ${run_plottingCNV} == true ]] ; then
    echo STEP6_plottingCNV
    Rscript ${codedir}/r.plotPerSample.r ${batchTag} ${inSample} ${RDRatioDir} ${plotDir} ${CBSDir} ${PoolList} ${MQList} ${dupdelList}
    sleep 1
  fi 
  if [[ ${run_getCNVPerSample} == true ]] ; then 
    echo STEP7_getCNVfromCBS
    python ${codedir}/python.getCNVPerSample.py ${batchTag} ${inSample} ${RDRatioDir} ${CBSDir} ${CNVDir} ${MQList} ${dupdelList} 
    sleep 1
  fi    
done
    
if [[ ${run_mergeCNVFiles} == true ]] ; then 
  echo STEP8_mergeCNVFiles
  python ${codedir}/python.mergeCNVFiles.py ${batchTag} ${inSampleInfoTxt} ${CNVDir} ${MQList} ${dupdelList}
  sleep 1
fi  
if [[ ${run_scoreCNV} == true ]] ; then 
  echo STEP9_scoreCNV
  python ${codedir}/python.scoreCNV.py ${batchTag} ${CNVDir} ${lqSampleDir} ${MQList} ${dupdelList} ${inScoringThTxt}
  sleep 1
fi     
