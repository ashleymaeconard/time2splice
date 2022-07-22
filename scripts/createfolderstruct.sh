#!/bin/bash
# createfolderstruct.sh
# Pranav Mahableshwarkar
# Last Mod. 7/15/22
# Purpose: Create a folder structure that can be used in the time2splice scripts. 

if [ $# -ne 1 ]; then
	echo $0: "Usage: 'bash createfolderstruct.sh'
            1) Directory to store all data and results. (i.e. ~/CLAMP_RNAi_Experiment)"
        exit 1
fi

INDIR=$1

echo "Making time2splice folder within the directory."
mkdir "${INDIR}/time2splice"

echo "Making data and results folders."
echo "Please put all input fastq files into inputdir/time2splice/data/fastq."
echo "These fastq files should be in folders labeled with replicate name (no _R* in naming)."
mkdir "${INDIR}/time2splice/data"
mkdir "${INDIR}/time2splice/results"
mkdir "${INDIR}/time2splice/data/fastq"
mkdir "${INDIR}/time2splice/data/count"
mkdir "${INDIR}/time2splice/results/preprocess"
mkdir "${INDIR}/time2splice/results/analysis"

mkdir "${INDIR}/time2splice/results/preprocess/fastqc"
mkdir "${INDIR}/time2splice/results/preprocess/trim_galore_fastqc"
mkdir "${INDIR}/time2splice/results/preprocess/alignment"
mkdir "${INDIR}/time2splice/results/analysis/salmon"
mkdir "${INDIR}/time2splice/results/analysis/suppa"

echo "Prior to running the scripts, make sure you have all packages downloaded and the
      time2splice conda environment activated. Please place all reference files in a reference
      folder. The working file path for that directory is needed for the scripts."

echo "Done."
