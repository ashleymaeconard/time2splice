#!/bin/bash
# run_macs2.sh
# Ashley Mae Conard
# Last Mod. June 24, 2019
# Purpose: Runs MACs2 for all .sorted.bam files in a given directory
# Resource: http://liulab.dfci.harvard.edu/MACS/FAQ.html

# Check to make sure input is correct
if [ $# -ne 4 ]; then
	echo $0: "Usage: ./run_macs2.sh 
			1) /PATH/TO/SORTED.BAMS/FOLDER(s)/ (containing folders of .sorted.bam files, assuming at most 2 subdirectories (e.g. only write up to /bowtie2/ here: .../results/bowtie2/clamp_sg_fm/MR_1-1/out.sorted.bam)
			2) PROCESSORS
			3) MFOLD_LOWER 
			4) MFOLD_UPPER (default is 5 50) \n" # NOTE: MAKE SURE source ~/venv/bin/activate FOR PROPER PYTHON ENVIRONMENT.
	exit 1
fi

# Assign input to variable names
INPUT_DIR=$1
PROCESSORS=$2
MFOLD_LOWER=$3
MFOLD_UPPER=$4

RESULTS_DIR=$(echo `dirname $INPUT_DIR`)

# Check to see if results/bowtie2/ directory exists
if [ -d "$RESULTS_DIR/macs2/" ]; then
	echo "${RESULTS_DIR}/macs2/ already exists."
    OUTPUT_DIR=${RESULTS_DIR}"/macs2/"
else
    OUTPUT_DIR=${RESULTS_DIR}"/macs2/"
    mkdir $OUTPUT_DIR
fi

# Create command script for Bowtie2
COMMAND_SCRIPT=${OUTPUT_DIR}"/run_macs2.txt" 
echo $COMMAND_SCRIPT

# Remove command script if already made
if [ -f ${COMMAND_SCRIPT} ]; then
	echo -e "REPLACING ${COMMAND_SCRIPT}\n"
	rm -rf ${COMMAND_SCRIPT}
fi

START=0
i=$START

echo "Processing .sorted.bam files with MACS2 on $PROCESSORS processors."

# Add commands to a script.txt file
for dir in $INPUT_DIR/*/
    do
    echo ${dir}
    # iterate through input directory (with sample names)
    if [ -d ${dir} ] ; then
        declare -a my_array
        unset my_array
        
        # iterate through only replicates to add them to list to be merged
        for file in ${dir}/*/*.sorted*.bam
            do       
  	    echo ${file} 
            # iterate through each sample to determine when to add 'wait' to script
            ((i = i + 1))

            # if number of processors reached, add wait to form a batch that will finish, and then process the next batch
            if (( ${i}%${PROCESSORS}==0 )); then
                echo "wait" >> $COMMAND_SCRIPT
            fi 

            # generate folder structure 
            fileName=$(echo `basename $file`) # get filename from .fq.gz
            fileReplicate=$(echo `basename $(dirname $file)`)
            fileSample=$(echo `basename $(dirname $(dirname $file))`)

            # create folder for all outputs per file by removing replicate identifer
            folderName=${OUTPUT_DIR}${fileSample}/${fileReplicate}/

            # output directory for BOWTIE2/sample_name
            mkdir -p ${folderName}

            #echo "Adding ${my_array[@]} to run_macs2.txt script."
            echo " "
            # ${my_array[@]} instead of $file
            echo "(macs2 callpeak -t $file -B -f AUTO --nomodel --SPMR --keep-dup all -g dm --trackline -n $fileName --cutoff-analysis --call-summits -p 0.05 -m ${MFOLD_LOWER} ${MFOLD_UPPER} --outdir $folderName 2> $folderName/$fileName.log) &" >> $COMMAND_SCRIPT
        done
    fi
done

# run command_script (.txt file saved in /results/bowtie2/
echo "RUNNING ${COMMAND_SCRIPT}"
bash ${COMMAND_SCRIPT}
