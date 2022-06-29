#!/bin/bash
# run_salmon.sh
# Ashley Mae Conard
# Last Mod. 7/4/2019
# Purpose: Run salmon to quantify transcript expression for treatment and control samples.

if [ $# -lt 7 ]; then
	echo $0: "Usage: ./run_salmon.sh /PATH/TO/INPUT_DIR/ containing inputs /PATH/TO/OUTPUT_DIR/ (suggest placing in /results/salmon/) REFERENCE_DIR TRANS_LOC NUM_THREADS (for each software command) PROCESSORS (per batch to run on cluster) IS_IT_PAIRED (0=no, 1=yes))"
	exit 1
fi

# Assign inputs to variable names
INPUT_DIR=$1 # 
OUTPUT_DIR=$2
REFERENCE_DIR=$3
TRANS_LOC=$4
NUM_THREADS=$5
PROCESSORS=$6
IS_IT_PAIRED=$7

# Salmon env
eval "$(conda shell.bash hook)"
conda activate salmon

# Make output directory if does not exist
mkdir -p $OUTPUT_DIR

# Create command script for Bowtie2
COMMAND_SCRIPT=${OUTPUT_DIR}"/run_salmon.txt" 

# Remove command script if already made
if [ -f ${COMMAND_SCRIPT} ]; then
	echo -e "REPLACING ${COMMAND_SCRIPT}\n"
	rm -rf ${COMMAND_SCRIPT}
fi

# Build index
trans_index=${REFERENCE_DIR}"/hisat2_dm6_trans_salmon_index"
if [ -d ${trans_index} ]; then
    echo -e "Already have transcriptome index $trans_index"
else
    echo "Building transcriptome index $trans_index"
    salmon index -t $TRANS_LOC -i $trans_index 
fi

START=0
i=$START

for dir in $INPUT_DIR/*/
    do
    echo "Subdir:" ${dir}
    
    # iterate through input directory
    if [ -d ${dir} ] ; then
        
        # Check if RNA-seq data is paired or not
        # SINGLE END READS
        if [ "$IS_IT_PAIRED" -ne "1" ]; then
           echo "Initiating single-end RNA-seq data processing.";
           for fastq in ${dir}/*.fastq*
               do
                   echo "Getting unpaired .fastq for $fastq"
                   ((i = i + 1))

                    # if number of processors reached, add wait to form a batch that will finish, and then process the next batch
                    if (( ${i}%${PROCESSORS}==0 )); then
                        echo "wait" >> $COMMAND_SCRIPT
                    fi 
                    fileName=$(echo `basename $fastq`) # get filename from .fq.gz
                    fileParent=$(echo `basename $(dirname $fastq)`)
                    
                    # create folder for all outputs per file by removing R1 and R2 (e.g. MTb8-8)                
                    fName=$(echo ${fileName%"_val_*"}); #Remove suffix
                    folderName=${OUTPUT_DIR}${fileParent}"/"${fName}
                    
                    # output directory for BOWTIE2/sample_name
                    mkdir -p ${folderName}
                    
                    # write Bowtie2 command to a file, followed by .bam creation and sorting
                    echo "Adding ${fileName} files to run_salmon.txt script"
                    echo " "     
                    echo "salmon quant -i $trans_index -l ISF --gcBias -1 ${dir}${fileName} -2 ${dir}${fileName} -p $NUM_THREADS -o $folderName">> $COMMAND_SCRIPT
           done
        
        # PAIRED END READS
        else
            echo "Initiating paired-end RNA-seq data processing.";
            # iterate through only the R1 replicates
            for R1 in ${dir}*_R1_*val_*.f*q.gz #MOD: Potentially unncessary formatting changes.
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
                    outFolder=${OUTPUT_DIR}"/"${fileParent}

                    echo $folderName

                    # get 2nd read pair
                    R2=${R1//"_R1_"/"_R2_"}
                    R2=${R2//"val_1"/"val_2"}

                    echo $R2

                    # output directory for BOWTIE2/sample_name
                    mkdir -p ${outFolder}
                    
                    # write Bowtie2 command to a file, followed by .bam creation and sorting
                    echo "Adding ${fileName} files to run_salmon.txt script"
                    echo " "
                    
                    # QUANTIFICATION
                    echo "salmon quant -i $trans_index -l ISF --gcBias -1 $R1 -2 $R2 -p $NUM_THREADS -o $outFolder">> $COMMAND_SCRIPT
            done
       fi
    fi
done
echo "conda deactivate">> $COMMAND_SCRIPT

# run command_script (.txt file saved in $OUTPUT_DIR)
bash $COMMAND_SCRIPT
