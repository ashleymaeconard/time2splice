#!/bin/bash
# run_process_macs2_output.sh
# Ashley Mae Conard
# Last Mod. June 25, 2019
# Runs MACS2 bdgcmp in for all macs2 callpeak files where there is a treatment and control

# Check to make sure input is correct
if [ $# -ne 2 ]; then
	echo $0: "Usage: ./macs2 
		1) /PATH/TO/MACS_OUTPUT/FOLDER/ (containing folders of MACS2 callpeak output files, (e.g. only need to write up to macs2 here: .../results/macs2/clamp_sg_fm/out.sorted.bam_peaks.xls) 
		2) PROCESSORS \n."
	exit 1
fi

# Assign input to variable names
INPUT_DIR=$1
PROCESSORS=$2

RESULTS_DIR=$(echo `dirname $INPUT_DIR`)
OUTPUT_DIR=${RESULTS_DIR}"/macs2/"

# Create command script for Bowtie2
COMMAND_SCRIPT=${OUTPUT_DIR}"/run_process_macs2_output.txt" 
echo $COMMAND_SCRIPT

# Remove command script if already made
if [ -f ${COMMAND_SCRIPT} ]; then
	echo -e "REPLACING ${COMMAND_SCRIPT}\n"
	rm -rf ${COMMAND_SCRIPT}
fi

START=0
i=$START

echo "Processing macs2 callpeak files with MACS2 bdgcmp on $PROCESSORS processors."

# Add commands to a script.txt file
for dir in $INPUT_DIR/*/
    do
    echo "Directory:" ${dir}
    # iterate through input directory (with sample names)
    if [ -d ${dir} ] ; then  
        for treat in ${dir}/*_treat_pileup.bdg # script will only proceed if there is a treatment and control
            do
                control="out.sorted.bam_control_lambda.bdg"
                # iterate through each sample to determine when to add 'wait' to script
                ((i = i + 1))
                # if number of processors reached, add wait to form a batch that will finish, and then process the next batch
                if (( ${i}%${PROCESSORS}==0 )); then
                    echo "wait" >> $COMMAND_SCRIPT
                fi 

                # generate folder structure 
                fileName=$(echo `basename $treat`) # get filename from .fq.gz
                fileSample=$(echo `basename $(dirname $treat)`)

                # create folder for all outputs per file by removing replicate identifer
                folderName=${OUTPUT_DIR}${fileSample}

                # write Bowtie2 command to a file, followed by .bam creation and sorting
                echo "Adding ${my_array[@]} to run_macs2.txt script."
                echo " "
                echo "(macs2 bdgcmp -t $treat -c $folderName/$control -o $folderName/out.sorted.bam_FE.bdg -m FE 2> $folderName/out.sorted.bam_FE.log; macs2 bdgcmp -t $treat -c $folderName/$control -o $folderName/out.sorted.bam_logLR.bdg -m logLR -p 0.00001 2> $folderName/out.sorted.bam_logLR.log; cd $folderName; Rscript out.sorted.bam_model.r) &" >> $COMMAND_SCRIPT  
        done
    fi
done

# run command_script (.txt file saved in /results/bowtie2/
bash $COMMAND_SCRIPT
