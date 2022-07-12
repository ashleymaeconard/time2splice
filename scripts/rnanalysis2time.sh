#!/bin/bash

#SBATCH -J RNAAnalysis

#SBATCH -n 16
#SBATCH --mem=8G
#SBATCH -t 4:00:00

# Create a SBATCH Job for all files
#SBATCH -o output/output.out

#SBATCH --mail-type=END
#SBATCH --mail-user=pranav_mahableshwarkar@brown.edu

echo "##### INSTRUCTIONS FOR PART TWO - TIME SPECIFIC ANALYSIS (SCRIPT THREE) ######
This should be run AFTER rnanlysis2.sh has been run for BOTH time points. This script will only run SUPPA and
suppa_formatting to compare across time points to allow you to get dpsi files. The general naming conventions from
rnanalysis.sh STILL apply. 

Output_Name -> The second four should be the different comparisons you want to make.
            -> C Female Time 1 v 2, C Male Time 1 v 2, M Female Time 1 v 2, M Male Time 1 v 2
 
Filter -> To recognize which  folder to get input from, the scripts looks at the folder names.
       -> So to differentiate you need to make four folders. (One for each sample type.)
       -> 1. m_control 2. f_control 3. m_mutant 4. f_mutant (mutant can be replaced with your specific mutant name.)
       -> If you have multiple time points, you can repeat this process to have DIFFERENT names for the second time point.
       -> eg: m_time2_control, etc. 

Converted -> If you use the files listed in this script for salmon, you should automatically be in Fbgn Names
             so your value for converted is listed as 1 in the script below (in 3_suppa_formatting). 

Once you have completed this, comment out the following 'exit 1' and this instruction and run! 
#####          COMPLETED          ######"
# exit 1

if [ $# -ne 4 ]; then
	echo $0: "Usage: 'sbatch preprocess.sh' or 'bash preprocess.sh'
            1) Data Input Directory (data in time2splice folder structure)
            2) Results Output Directory (results in time2splice folder structure)
            3) Reference Directory (look below to see specific files)
            4) SUPPA PATH (install from GitHub)

            This is Part Two of the RNA Analysis scripts. Please read the following before running."
        exit 1
fi

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
REFSEQ="${REF_DIR}/fbtr_refseq.tsv" 
FBIOE="${REF_DIR}/flybase_events/flybase.events.ioe" # You might have to look in to acquiring this file. 
GENEGTF="${REF_DIR}/genes.gtf"
DMELGTF="${REF_DIR}/dmel-all-r6.46.gtf"
TRANSIOI="${REF_DIR}/fb_trans.isoforms.ioi"

# Checks to make sure that ALL reference files exist. 
if [ ! -e $REFSEQ ] 
then
    echo $0:"${REFSEQ} not found. Please make sure file path is correct."
    exit 1
elif [ ! -e $FBIOE ] 
then 
    echo $0:"${FBIOE} not found. Please make sure file path is correct."
    exit 1
elif [ ! -e $GENEGTF ] 
then 
    echo $0:"${GENEGTF} not found. Please make sure file path is correct."
    exit 1
elif [ ! -e $DMELGTF ] 
then 
    echo $0:"${DMELGTF} not found. Please make sure file path is correct."
    exit 1
elif [ ! -e $TRANSIOI ] 
then 
    echo $0:"${TRANSIOI} not found. Please make sure file path is correct."
    exit 1
fi

echo "All Reference Files Found. Starting..."

SUPPARESULTS="${RES_DIR}/analysis/suppa"

# USER INPUT: My respective folders are f_L3_BrControl, m_L3_BrControl, etc. etc.
#             Your input MUST match the information included below. Mismatches will disable the scripts.
TIME1CON="L3_BrControl"
TIME1MUT="L3_BrCLAMPRNAi"
TIME2CON="Ad_BrControl"
TIME2MUT="Ad_BrCLAMPRNAi"

# These are the four comparisons in 3. 
CAT1="FL3vAdControl"
CAT2="ML3vAdControl"
CAT3="FL3vAdCLAMPRNAi"
CAT4="ML3vAdCLAMPRNAi"

# These will identify the two respective samples for CAT1, 2, 3, and 4 in 4_suppa.sh
IDENT1="f*BrControl"
IDENT2="m*BrControl"
IDENT3="f*BrCLAMPRNAi"
IDENT4="m*BrCLAMPRNAi"

# This is used in 6_get_bias_genes to name the Male and Female samples. 
FFLYTYPE1="L3F"
MFLYTYPE1="L3M"
FFLYTYPE2="AdF"
MFLYTYPE2="AdM"

# 3. A Merge must be run for EVERY 4_suppa.sh data comparison that is run (folder required to run SUPPA).
# a) Category 1.
python3 ./rna/3_suppa_formatting.py ${SUPPARESULTS}/ ${SUPPARESULTS}/${CAT1} ${IDENT1} ${REFSEQ} 1
# b) Category 2.
python3 ./rna/3_suppa_formatting.py ${SUPPARESULTS}/ ${SUPPARESULTS}/${CAT2} ${IDENT2} ${REFSEQ} 1
# c) Category 3.
python3 ./rna/3_suppa_formatting.py ${SUPPARESULTS}/ ${SUPPARESULTS}/${CAT3} ${IDENT3} ${REFSEQ} 1
# d) Category 4.
python3 ./rna/3_suppa_formatting.py ${SUPPARESULTS}/ ${SUPPARESULTS}/${CAT4} ${IDENT4} ${REFSEQ} 1

# 4. Run SUPPA Differential Splicing Analysis. For each run, there must be a corresponding merge call made. 
#    Need to make sure that input merged_iso_tpm file ONLY has replicate names in first row. Check the if group at 229 of this script if errors arise.
module load R/4.2.0 gcc/10.2 pcre2/10.35 intel/2020.2 texlive/2018
# a) Control Males v Mutant Males. 
./rna/4_suppa.sh ${SUPPARESULTS}/${CAT1} ${FBIOE} $SUPPA_PATH ${GENEGTF} f_${TIME1CON} f_${TIME2CON} $TRANSIOI
# b) Control Females v Mutant Females.
./rna/4_suppa.sh ${SUPPARESULTS}/${CAT2} ${FBIOE} $SUPPA_PATH ${GENEGTF} m_${TIME1CON} m_${TIME2CON} $TRANSIOI
# c) Control Males v Control Females
./rna/4_suppa.sh ${SUPPARESULTS}/${CAT3} ${FBIOE} $SUPPA_PATH ${GENEGTF} f_${TIME1MUT} f_${TIME2MUT} $TRANSIOI
# d) Mutant Males v Mutant Females
./rna/4_suppa.sh ${SUPPARESULTS}/${CAT4} ${FBIOE} $SUPPA_PATH ${GENEGTF} m_${TIME1MUT} m_${TIME2MUT} $TRANSIOI

# 7. Plots Splicing: You HAVE to run Jupyter Notebook 7_plots_splicing on both time points before running Step 8. Step 8 will require the
#    two .tsv files that the individual timepoint analysis produces. 

# 8. If you are using time-specific samples, run this step. Read through the notebook for debugging purposes for input formatting.
echo "Done"