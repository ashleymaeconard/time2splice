#!/bin/bash
# run_BWA.sh
# Ashley Mae Conard
# Last Mod. May 22, 2020
# Purpose: Runs BWA in for all .fastq.gz files in a given directory

# Check to make sure input is correct
if [ $# -ne 7 ]; then
	echo $0: "Usage: ./run_BWA.sh 
			1) /PATH/TO/NAMED_FASTQ_DIR/ (NAMED folders containing replicate fastq.gz files, e.g. ../data/clampChip_f/a.fastq.gz, only go to ../data/) 
			2) /PATH/TO/RESULTS_DIR/ (desired location to add or create results directory) 
			3) NUM_CHAR_KEEP_NAME (i.e. MR-23_R1_001.fastq.gz, keep 5 so MR-23) 
			4) PROCESSORS (per batch, recommend 4 with 7 threads) 
			5) THREADS (number algnment threads to launch) 
			6) PAIRED_OR_NOT (YES is 1, NO is 0) 
			7) GENOME_LOCATION (/PATH/TO/GENOME/genome)"
	exit 1
fi
# genome can be found here: /data/compbio/aconard/genomes/BDGP6/bwaGenome/genome.fa

# Assign input to variable names
INPUT_DIR=$1
RESULTS_DIR=$2
NUM_CHAR_KEEP_NAME=$3
PROCESSORS=$4
THREADS=$5
PAIRED_OR_NOT=$6
LOC_GENOME=$7 

# Check to see if results directory exists
if [ -d "$RESULTS_DIR" ]; then
	echo "${RESULTS_DIR} already exists."
else
    mkdir $RESULTS_DIR
fi

# Check to see if results/bwa/ directory exists
if [ -d "$RESULTS_DIR/bwa/" ]; then
	echo "${RESULTS_DIR}/bwa/ already exists."
    OUTPUT_DIR=${RESULTS_DIR}"/bwa/"
else
    OUTPUT_DIR=${RESULTS_DIR}"/bwa/"
    mkdir $OUTPUT_DIR
fi

# Create command script for BWA
COMMAND_SCRIPT=${OUTPUT_DIR}"/run_BWA.txt" 
echo $COMMAND_SCRIPT

# Remove command script if already made
if [ -f ${COMMAND_SCRIPT} ]; then
	echo -e "REPLACING ${COMMAND_SCRIPT}\n"
	rm -rf ${COMMAND_SCRIPT}
fi
    
START=0
i=$START

echo "Processing .fastq files with BWA on $PROCESSORS processors, keeping $NUM_CHAR_KEEP_NAME characters of each replicate name."
# Add commands to a script.txt file
for dir in $INPUT_DIR/*/
    do
    echo ${dir}
    # iterate through input directory
    if [ -d ${dir} ] ; then
        
        # Check if RNA-seq data is paired or not
        # SINGLE END READS
        if [ "$PAIRED_OR_NOT" -ne "1" ]; then
           echo "Initiating single-end RNA-seq data processing.";
           for fastq in ${dir}/*.fastq.gz*
               do
                   echo "Getting unpaired .fastq for $fastq"
                   ((i = i + 1))

                    # if number of processors reached, add wait to form a batch that will finish, and then process the next batch
                    if (( ${i}%${PROCESSORS}==0 )); then
                        echo "wait" >> $COMMAND_SCRIPT
                    fi 
                    fileName=$(echo `basename $fastq`) # get filename from .fq.gz
                    fileParent=$(echo `basename $(dirname $fastq)`)
                    
                    # create folder for all BWA outputs per file by removing R1
                    fName=$(echo ${fileName::$NUM_CHAR_KEEP_NAME}) 
                    folderName=${OUTPUT_DIR}${fileParent}"/"${fName}
                    
                    # output directory for BWA/sample_name
                    mkdir -p ${folderName}
                    
                    # write BWA command to a file, followed by .bam creation and sorting
                    echo "Adding ${fileName} to run_BWA.txt script"
                    echo " "
                    echo "(bwa mem -t ${THREADS} $LOC_GENOME ${fastq} > $folderName/out.sam; samtools view -bS $folderName/out.sam > $folderName/out.bam; samtools sort $folderName/out.bam -o $folderName/out.sorted.bam; rm -rf $folderName/out.bam; samtools flagstat $folderName/out.sorted.bam  > $folderName/out.sorted.stats.txt) &" >> $COMMAND_SCRIPT 

           done
        
        # PAIRED END READS
        else
            echo "Initiating paired-end RNA-seq data processing.";
            # iterate through only the R1 replicates
            for R1 in ${dir}/*R1*
                do
                    echo "Getting paired .fastq for $R1"
                    # iterate through each R1 to determine when to add 'wait' to script
                    ((i = i + 1))

                    # if number of processors reached, add wait to form a batch that will finish, and then process the next batch
                    if (( ${i}%${PROCESSORS}==0 )); then
                        echo "wait" >> $COMMAND_SCRIPT
                    fi 

                    # generate the bwa command for the next R1 and R2
                    fileName=$(echo `basename $R1`) # get filename from .fq.gz
                    fileParent=$(echo `basename $(dirname $R1)`)

                    # create folder for all BWA outputs per file by removing R1 and R2 (e.g. MTb8-8)
                    fName=$(echo ${fileName::$NUM_CHAR_KEEP_NAME}) 
                    folderName=${OUTPUT_DIR}${fileParent}"/"${fName}

                    # get 2nd read pair
                    R2=${R1//R1/R2}

                    # output directory for BWA/sample_name
                    mkdir -p ${folderName}

                    # write BWA2 command to a file, followed by .bam creation and sorting
                    echo "Adding ${fileName} to run_BWA.txt script"
                    echo " "
		    echo "(bwa mem -t ${THREADS} $LOC_GENOME $R1 $R2 > $folderName/out.sam; samtools view -bS $folderName/out.sam > $folderName/out.bam; samtools sort $folderName/out.bam -o $folderName/out.sorted.bam; rm -rf $folderName/out.bam; samtools flagstat $folderName/out.sorted.bam  > $folderName/out.sorted.stats.txt) &" >> $COMMAND_SCRIPT 

            done
       fi
    fi
done

# run command_script (.txt file saved in $OUTPUT_DIR)
bash $COMMAND_SCRIPT
