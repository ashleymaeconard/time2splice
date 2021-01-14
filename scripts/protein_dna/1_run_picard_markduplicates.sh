#!/bin/bash
# run_picard_markduplicates.sh
# Ashley Mae Conard
# Last Mod. June 23, 2019
# Purpose: Runs Picard's MarkDuplicates in for all .sorted.bam files in a given directory

PICARD="/data/compbio/aconard/tools/gatk-4.1.7.0/"

# Check to make sure input is correct
if [ $# -ne 2 ]; then
	echo $0: Usage: "./run_picard_markduplicates.sh 
			1) /PATH/TO/SORTED.BAMS/FOLDER(s)/ (containing folders of .sorted.bam files, assuming at most 2 subdirectories (e.g. only write up to bowtie2 here: .../results/bowtie2/clamp_sg_fm/MR_1-1/out.sorted.bam) 
			2) PROCESSORS"
	exit 1
fi
#e.g. ./run_picard_markduplicates.sh /data/compbio/aconard/cut-n-run/results/bowtie2/ 25

# Assign input to variable names
INPUT_DIR=$1
PROCESSORS=$2

# Finding /results/ parent dir for results and creating output_dir
RESULTS_DIR=$(echo `dirname $INPUT_DIR`)
OUTPUT_DIR=${RESULTS_DIR}"/markduplicates/"

# Check to see if results/markduplicates/ directory exists
if [ -d "$RESULTS_DIR/markduplicates/" ]; then
	echo "${RESULTS_DIR}/markduplicates/ already exists."
    OUTPUT_DIR=${RESULTS_DIR}"/markduplicates/"
else
    OUTPUT_DIR=${RESULTS_DIR}"/markduplicates/"
    mkdir $OUTPUT_DIR
fi

# Create command script for Picard's MarkDuplicates
COMMAND_SCRIPT=${OUTPUT_DIR}"/run_picard_markduplicates.txt" 
echo $COMMAND_SCRIPT

# Remove run_picard_markduplicates.txt command file if already made
if [ -f ${COMMAND_SCRIPT} ]; then
	echo "REPLACING ${COMMAND_SCRIPT}\n"
	rm -rf ${COMMAND_SCRIPT}
fi

echo "Removing duplicates files with Picard's MarkDuplicates on $PROCESSORS processors."

START=0
i=$START

# Add commands to a script.txt file
for dir in $INPUT_DIR/*/
    do
    echo ${dir}
    # iterate through input directory
    if [ -d ${dir} ] ; then
    
        # iterate through only all sorted.bam files
        for file in ${dir}/*/*sorted.bam
            do
                echo "Getting $file to remove duplicates"
                # iterate through each sorted.bam to determine when to add 'wait' to script
                ((i = i + 1))

                # if number of processors reached, add wait to form a batch that will finish, and then process the next batch
                if (( ${i}%${PROCESSORS}==0 )); then
                    echo "wait" >> $COMMAND_SCRIPT
                fi
                
                # generate the MarkDuplicates command for the sorted.bam
                #fileName=$(echo `basename $file`) # get filename from .fq.gz
                fileReplicate=$(echo `basename $(dirname $file)`)
                fileSample=$(echo `basename $(dirname $(dirname $file))`)
                
                # create folder for all MarkDuplicates outputs per file by removing R1 and R2 (e.g. MTb8-8)
                #fName=$(echo ${fileName::$NUM_CHAR_KEEP_NAME}) 
                folderName=${OUTPUT_DIR}${fileSample}/${fileReplicate}

                #output directory for markDuplicates/sample_name
                mkdir -p ${folderName}

                # write Picard's MarkDuplicates command to a file
                echo "($PICARD/gatk --java-options '-Djava.io.tmpdir=/ltmp/ -Xmx12G' MarkDuplicates -I $file -O $folderName/out.sorted_nodup.bam -M ${folderName}/marked_dup_metrics.txt) &" >> $COMMAND_SCRIPT

        done
    fi
done

# run command_script (.txt file saved in $OUTPUT_DIR)
bash $COMMAND_SCRIPT
