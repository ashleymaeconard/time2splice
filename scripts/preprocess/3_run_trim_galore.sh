#!/bin/bash
# run_trim_Galore.sh
# Ashley Mae Conard
# Last Mod. April 6, 2020
# Runs trim Galore THEN FastQC in for all .fastq.gz files in a given directory

# Check to make sure input is correct
if [ $# -ne 6 ]; then
	echo $0: "Usage: ./run_trim_Galore.sh 
		1) NUM_COMMANDS 
		2) /PATH/TO/DIR/ (with 2 subdirs (SAMPLE_NAME/REP_NAME/fastq.gz)) 
		3) /PATH/TO/OUTPUT_DIR/ (suggest placing in /time2splice/results/trim_galore_fastqc) 
		4) Q_SCORE (recommend 30) 
		5) PATH_TO_CUTADAPT (e.g. /time2splice_env/bin/cutadapt) 
        6) FOLDERS BETWEEN INDIR AND FASTQ FILES (ex: Indir/*/.fastq.gz would be 1)" #MODIFICATION: Added #6 for increased compatibility. 
	exit 1
fi

echo "If 3_trim_galore fails, this could be due to naming assumptions made in the script."

# CORES_PER_COMMAND_RUN
NUM_COMMANDS=$1
INDIR=$2
OUTDIR=$3
Q_SCORE=$4
PATH_CUTADAPT=$5
FOLDERS_BETWEEN=$6 #MOD: Loading the modification for #6. 

#Make output directory if necessary
mkdir -p $OUTDIR

# Create command script for Bowtie2
COMMAND_SCRIPT=${OUTDIR}/"commands_trim_galore.txt" 

# Remove command script if already made
if [ -f ${COMMAND_SCRIPT} ]; then
	echo -e "REPLACING ${COMMAND_SCRIPT}\n"
	rm -rf ${COMMAND_SCRIPT}
fi

START=0
i=$START

# Determine how many "*/" are needed in the folder structure. 
chr="*/"
folders=
for ((i = 0; i < FOLDERS_BETWEEN; i++)); do 
    folders="$folders$chr"
done

# Iterate through only the R1 replicates and match the second replicate
echo "Initiating paired-end RNA-seq data processing.";
for R1 in "${INDIR}/${folders}*_R1*f*q*" #MOD.
    do
        echo $R1
        echo "Getting paired .fastq for $R1"
        # iterate through each R1 to determine when to add 'wait' to script
        ((i = i + 1))

        # if number of processors reached, add wait to form a batch that will finish, and then process the next batch
        if (( ${i}%${NUM_COMMANDS}==0 )); then
            echo "wait" >> $COMMAND_SCRIPT
        fi 

        R2=${R1/"_R1"/"_R2"}
        echo $R2

        fileName=$(echo `basename $R1`) # get filename from .fq.gz
        fileParent=$(echo `basename $(dirname $R1)`)
        fileName2=$(echo `basename $R2`) #MOD: Need to find the second file path. 

        # create folder for all ooutputs per file by removing R1 and R2 (e.g. MTb8-8)
        # OUTDIR MUST be Sample specific. 
        folderName=${OUTDIR}/ #${fileParent}/ #"/"${fileName}
        touch "${folderName}/${fileName}_trimming_report.txt" #MOD: Added these two lines so the report can go somewhere. 
        touch "${folderName}/${fileName2}_trimming_report.txt"

        # output directory for BOWTIE2/sample_name
        #mkdir -p ${folderName}

        #MOD: Removed () and &           
        echo "trim_galore -q ${Q_SCORE} --phred33 --gzip --fastqc --paired --retain_unpaired --path_to_cutadapt ${PATH_CUTADAPT} --output_dir ${folderName} ${R1} ${R2}">> $COMMAND_SCRIPT   
done
echo "running command script"
bash $COMMAND_SCRIPT      
echo "done"
