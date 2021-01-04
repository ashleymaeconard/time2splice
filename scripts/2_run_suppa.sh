#!/bin/bash
# run_suppa.sh
# Ashley Mae Conard
# Last Mod: 7/5/2019
# Purpose: Run Suppa for treatment and control samples.

if [ $# -lt 4 ]; then
	echo $0: "Usage: ./1_run_suppa.sh /PATH/TO/INPUT_DIR/ containing inputs /PATH/TO/OUTPUT_DIR/ (suggest placing in /results/suppa/) TRANSCRIPTOME_LOC PROCESSORS (per batch to run on cluster)"
	exit 1
fi

# Assign inputs to variable names
INPUT_DIR=$1 # 
OUTPUT_DIR=$2
TRANS_LOC=$3
PROCESSORS=$4

# Make output directory if does not exist
mkdir -p $OUTPUT_DIR

# Create command script for Bowtie2
COMMAND_SCRIPT=${OUTPUT_DIR}"/run_suppa.txt" 

# Remove command script if already made
if [ -f ${COMMAND_SCRIPT} ]; then
	echo -e "REPLACING ${COMMAND_SCRIPT}\n"
	rm -rf ${COMMAND_SCRIPT}
fi

START=0
i=$START

for dir in $INPUT_DIR/*/
    do
    echo "Subdir:" ${dir}
    
    # iterate through input directory
    if [ -d ${dir} ] ; then  
        echo "Initiating suppa ID of alternative splicing events.";
        # iterate through only the R1 replicates
        for file in ${dir}/*/quant.sf
            do
                echo "Getting quantification information from $file"
                # iterate through each R1 to determine when to add 'wait' to script
                ((i = i + 1))

                # if number of processors reached, add wait to form a batch that will finish, and then process the next batch
                if (( ${i}%${PROCESSORS}==0 )); then
                    echo "wait" >> $COMMAND_SCRIPT
                fi 

                # generate the bowtie2 command for the next R1 and R2
                fileName=$(echo `basename $file`) # get filename from .fq.gz
                fileParent=$(echo `basename $(dirname $file)`)
                fileParent2=$(echo `basename $(dirname $(dirname $file))`)
                
                # create folder for all Bowtie2 outputs per file by removing R1 and R2 (e.g. MTb8-8)             
                folderName=${OUTPUT_DIR}'/'${fileParent2}'/'${fileParent}'/'
                echo "folderName $folderName"
                # output directory for BOWTIE2/sample_name
                mkdir -p ${folderName}

                # write to screen what file is being added to suppa script
                echo "Adding ${fileName} files to run_suppa.txt script"
                echo " "

                # QUANTIFICATION
                echo "(multipleFieldSelection.py -i $file -k 1 -f 4 -o $folderName/iso_tmp.txt ) &">> $COMMAND_SCRIPT
        done
    fi
done

# run command_script (.txt file saved in $OUTPUT_DIR)
bash $COMMAND_SCRIPT

