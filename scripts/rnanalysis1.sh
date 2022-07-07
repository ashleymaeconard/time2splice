#!/bin/bash

#SBATCH -J RNAAnalysis

#SBATCH -n 24
#SBATCH --mem=8G
#SBATCH -t 12:00:00

# Create a SBATCH Job for all files
#SBATCH -o output/output.out

#SBATCH --mail-type=END
#SBATCH --mail-user=pranav_mahableshwarkar@brown.edu

if [ $# -ne 4 ]; then
	echo $0: "Usage: 'sbatch preprocess.sh' or 'bash preprocess.sh'
            1) Data Input Directory (data in time2splice folder structure)
            2) Results Output Directory (results in time2splice folder structure)
            3) Reference Directory (look below to see specific files)
            4) SUPPA PATH (install from GitHub)
            This is ONLY part 1 of the RNA Analysis Scripts."
        exit 1
fi

echo "The RNA Analysis is done in two parts. This is ONLY steps 1 and 2."

# Activate in the time2splice conda environment. 
module load anaconda/2020.02
source /gpfs/runtime/opt/anaconda/2020.02/etc/profile.d/conda.sh
conda activate time2splice

DATA_DIR=$1
RES_DIR=$2 
REF_DIR=$3
SUPPA_PATH=$4

# Load the appropriate modules!
module load samtools/1.13
module load salmon/1.3.0
module load gcc/8.3

# Establish File Paths for all reference files. 
DMELALLTRANSCRIPT="${REF_DIR}/dmel-all-transcript-r6.46.fasta"

# Checks to make sure that ALL reference files exist. 
if [ ! -e $DMELALLTRANSCRIPT ] 
then
    echo $0:"${DMELALLTRANSCRIPT} not found. Please make sure file path is correct."
    exit 1
fi
echo "Reference Files Found. Starting..."

# IMPORTANT: Step 1 needs to be modified based on the formatting of the trim galore fastq file (Line 89).
./rna/1_run_salmon.sh ${RES_DIR}/preprocess/trim_galore_fastqc ${RES_DIR}/analysis/salmon ${REF_DIR}/reference ${DMELALLTRANSCRIPT} 16 4 1

# 2. Run SUPPA to isolate splicing differences within a sample type.
./rna/2_run_suppa.sh ${RES_DIR}/analysis/salmon ${RES_DIR}/analysis/suppa ${DMELALLTRANSCRIPT}  1

echo "Finished! Please read the instructions in rnanlysis2.sh before running the next script."

