#!/bin/bash
# run_Bowtie2.sh
# Ashley Mae Conard
# Last Mod. June 21, 2019
# Purpose: Runs Bowtie2 in for all .fastq.gz files in a given directory

# Check to make sure input is correct
if [ $# -ne 8 ]; then
	echo $0: "Usage: ./run_Bowtie.sh /PATH/TO/NAMED_FASTQ_DIR/ (NAMED folders containing replicate fastq.gz files, e.g. ../data/clampChip_f/a.fastq.gz, only go to ../data/) /PATH/TO/RESULTS_DIR/ (desired location to add or create results directory) SUFFIX_TO_REMOVE (e.g. R1_001.fastq.gz) PROCESSORS (per batch, recommend 4 with 7 threads) THREADS (number algnment threads to launch) PAIRED_OR_NOT (YES is 1, NO is 0) GENOME_LOCATION (/PATH/TO/GENOME/genome) SEARCH_TERM (e.g. _val_1.fq.gz, or *.fastq.gz*)" 
	exit 1
fi
# genome can be found here: /data/compbio/aconard/rnaseq_larschan/BDGP6/Drosophila_melanogaster/Ensembl/BDGP6/Sequence/Bowtie2Index/genome

# Assign input to variable names
INPUT_DIR=$1
RESULTS_DIR=$2
SUFFIX_TO_REMOVE=$3
PROCESSORS=$4
THREADS=$5
PAIRED_OR_NOT=$6
LOC_GENOME=$7
SEARCH_TERM=$8

# Check to see if results directory exists
if [ -d "$RESULTS_DIR" ]; then
	echo "${RESULTS_DIR} already exists."
else
    mkdir $RESULTS_DIR
fi

# Check to see if results/bowtie2/ directory exists
if [ -d "$RESULTS_DIR/bowtie2/" ]; then
	echo "${RESULTS_DIR}/bowtie2/ already exists."
    OUTPUT_DIR=${RESULTS_DIR}"/bowtie2/"
else
    OUTPUT_DIR=${RESULTS_DIR}"/bowtie2/"
    mkdir $OUTPUT_DIR
fi

# Create command script for Bowtie2
COMMAND_SCRIPT=${OUTPUT_DIR}"/commands_Bowtie2.txt" 
echo $COMMAND_SCRIPT

# Remove command script if already made
if [ -f ${COMMAND_SCRIPT} ]; then
	echo -e "REPLACING ${COMMAND_SCRIPT}\n"
	rm -rf ${COMMAND_SCRIPT}
fi
    
START=0
i=$START

echo "Processing .fastq files with Bowtie2 on $PROCESSORS processors, removing $SUFFIX_TO_REMOVE from each replicate name."
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
           for fastq in ${dir}/${SEARCH_TERM}
               do
                   echo "Getting unpaired .fastq for $fastq"
                   ((i = i + 1))

                    # if number of processors reached, add wait to form a batch that will finish, and then process the next batch
                    if (( ${i}%${PROCESSORS}==0 )); then
                        echo "wait" >> $COMMAND_SCRIPT
                    fi 
                    fileName=$(echo `basename $fastq`) # get filename from .fq.gz
                    fileParent=$(echo `basename $(dirname $fastq)`)
                    
                    # create folder for all Bowtie2 outputs per file by removing R1 and R2 (e.g. MTb8-8) 
                    fName=$(echo ${fileName%$SUFFIX_TO_REMOVE}) 
                    folderName=${OUTPUT_DIR}${fileParent}"/"${fName}
                    
                    # output directory for BOWTIE2/sample_name
                    mkdir -p ${folderName}
                    
                    # write Bowtie2 command to a file, followed by .bam creation and sorting
                    echo "Adding ${fileName} to run_Bowtie.txt script"
                    echo " "
                    echo "(bowtie2 -p ${THREADS} -x $LOC_GENOME -U $fastq --un-gz $folderName/out_un.sam --al-gz $folderName/out_al.sam --met-file $folderName/out_met-file.tsv -S $folderName/out.sam 2> $folderName/alignment_info.txt; samtools view -bS $folderName/out.sam > $folderName/out.bam; rm -rf $folderName/out.sam; samtools sort $folderName/out.bam -o $folderName/out.sorted.bam) &" >> $COMMAND_SCRIPT  
           done
        
        # PAIRED END READS
        else
            echo "Initiating paired-end RNA-seq data processing.";
            # iterate through only the R1 replicates
            for R1 in ${dir}/*R1*${SEARCH_TERM}
                do
                    echo "Getting paired .fastq for $R1"
                    # iterate through each R1 to determine when to add 'wait' to script
                    ((i = i + 1))

                    # if number of processors reached, add wait to form a batch that will finish, and then process the next batch
                    if (( ${i}%${PROCESSORS}==0 )); then
                        echo "wait" >> $COMMAND_SCRIPT
                    fi 

                    # generate the bowtie2 command for the next R1 and R2
                    fileName=$(echo `basename $R1`) # get filename from .fq.gz
                    fileParent=$(echo `basename $(dirname $R1)`)

                    # create folder for all Bowtie2 outputs per file by removing R1 and R2 (e.g. MTb8-8)
                    fName=${fileName} # $(echo ${fileName%$SUFFIX_TO_REMOVE})  
                    folderName=${OUTPUT_DIR}${fileParent}"/"${fName}

                    # get 2nd read pair
                    R2=${R1//R1/R2}
                    R2v=${R2//val_1/val_2}

                    # output directory for BOWTIE2/sample_name
                    mkdir -p ${folderName}

                    # write Bowtie2 command to a file, followed by .bam creation and sorting
                    echo "Adding ${fileName} to commands_Bowtie.txt script"
                    echo " "
                    echo "(bowtie2 -p ${THREADS} -x $LOC_GENOME -1 $R1 -2 $R2v --al-conc-gz $folderName/out_al-conc.sam --un-conc-gz $folderName/out_unconc.sam --met-file $folderName/out_met-file.tsv -S $folderName/out.sam 2> $folderName/alignment_info.txt; samtools view -bS $folderName/out.sam > $folderName/out.bam; rm -rf $folderName/out.sam; samtools sort $folderName/out.bam -o $folderName/out.sorted.bam) &" >> $COMMAND_SCRIPT  
            done
       fi
    fi
done

# run command_script (.txt file saved in $OUTPUT_DIR)
bash $COMMAND_SCRIPT
tail -n 1 ${OUTPUT_DIR}*/*/align*.txt > ${OUTPUT_DIR}/all_percent_alignment_bowtie.txt
