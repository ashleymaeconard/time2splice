#!/bin/bash
# run_HISAT2.sh
# Ashley Mae Conard
# Last Mod. April 5, 2020
# Purpose: Runs HISAT2 in for all .fastq.gz files in a given directory

# Check to make sure input is correct
if [ $# -ne 8 ]; then
	echo $0: "Usage: ./run_HISAT2.sh /PATH/TO/NAMED_FASTQ_DIR/ (NAMED folders containing replicate fastq.gz files) /PATH/TO/RESULTS_DIR/ (desired location to add or create results directory) SUFFIX_TO_REMOVE (e.g. R1_001.fastq.gz) THREADS (per HISAT2 command) PROCESSORS (number HISAT2 commands to launch per batch) PAIRED_OR_NOT (YES is 1, NO is 0) REF_TYPE \(trans, genom\) SEARCH_TERM (e.g. _val_1.fq.gz, or *.fastq.gz*)" 
	exit 1
fi

# Assign input to variable names
INPUT_DIR=$1
RESULTS_DIR=$2
SUFFIX_TO_REMOVE=$3
THREADS=$4
PROCESSORS=$5
PAIRED_OR_NOT=$6
REF_TYPE=$7
SEARCH_TERM=$8

# Check to see if results directory exists
if [ -d "$RESULTS_DIR" ]; then
	echo "${RESULTS_DIR} already exists."
else
    mkdir $RESULTS_DIR
fi

# Check to see if results/bowtie2/ directory exists
if [ -d "$RESULTS_DIR/hisat2/" ]; then
	echo "${RESULTS_DIR}/hisat2/ already exists."
    OUTPUT_DIR=${RESULTS_DIR}"/hisat2/"
else
    OUTPUT_DIR=${RESULTS_DIR}"/hisat2/"
    mkdir $OUTPUT_DIR
fi

# Create command script for Bowtie2
COMMAND_SCRIPT=${OUTPUT_DIR}"commands_HISAT2.txt" 

# Remove command script if already made
if [ -f ${COMMAND_SCRIPT} ]; then
	echo -e "REPLACING ${COMMAND_SCRIPT}\n"
	rm -rf ${COMMAND_SCRIPT}
fi

if [ "$REF_TYPE" = "genom" ]; then
    genome_HISAT2="/data/compbio/aconard/BDGP6/hisat2Genome/bdgp6/genome"
                          
# elif [ "$REF_TYPE" = "trans" ]; then
#     genome_HISAT2="/data/compbio/aconard/BDGP6/transcriptome_dir/pub/infphilo/hisat2/data/bdgp6_tran/genome_tran"
#         genome_HTSeq="/data/compbio/aconard/BDGP6/transcriptome_dir/pub/infphilo/hisat2/data/bdgp6_tran/Drosophila_melanogaster.BDGP6.84.gtf"
else
    echo "ERROR - must input 'genom' or 'trans' for reference type"
    exit 1
fi 

START=0
i=$START

echo "Processing .fastq files with HISAT2 on $PROCESSORS processors, keeping $SUFFIX_TO_REMOVE characters of each replicate name."
# Add commands to a script.txt file
for dir in $INPUT_DIR/*/
    do
    echo "Subdir:" ${dir}
    
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
                    echo "Adding ${fileName} to run_HISAT2.txt script"
                    echo " "                    
                    echo "(hisat2 -p ${THREADS} --dta -x $genome_HISAT2 -U ${dir}/$fileName -t --un-gz $folderName/out_unalign --al-gz $folderName/out_al_atleastonce --known-splicesite-infile /data/compbio/aconard/BDGP6/splicesites.txt --novel-splicesite-outfile $folderName/novel_splicesite --summary-file $folderName/summaryfile.txt --met-file $folderName/met-file.txt -S $folderName/a.sam 2> $folderName/alignmentsummary.txt; samtools sort -O sam -T $folderName/a_sorted -o $folderName/a_sorted.sam $folderName/a.sam; samtools view -S -b $folderName/a_sorted.sam > $folderName/a_sorted.bam; rm $folderName/a.sam $folderName/a_sorted.sam) &">> $COMMAND_SCRIPT
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
                    fName=$(echo ${fileName%$SUFIX_TO_REMOVE}) 
                    folderName=${OUTPUT_DIR}${fileParent}"/"${fName}

                    # get 2nd read pair
                    R2=${R1//R1/R2}
                    R2v=${R2//val_1/val_2}
                    
                    # output directory for BOWTIE2/sample_name
                    mkdir -p ${folderName}

                    # write Bowtie2 command to a file, followed by .bam creation and sorting
                    echo "Adding ${fileName} to run_Bowtie.txt script"
                    echo " "
                    echo "(hisat2 -p ${THREADS} --dta -x $genome_HISAT2 -1 $R1 -2 $R2v -t --un-conc-gz $folderName/out_unconc --al-conc-gz $folderName/out_al_conc_atleastonce --known-splicesite-infile /data/compbio/aconard/BDGP6/splicesites.txt --novel-splicesite-outfile $folderName/novel_splicesite --summary-file $folderName/summaryfile.txt --met-file $folderName/met-file.txt -S $folderName/a.sam 2> $folderName/alignmentsummary.txt; samtools sort -O sam -T $folderName/a_sorted -o $folderName/a_sorted.sam $folderName/a.sam; samtools view -S -b $folderName/a_sorted.sam > $folderName/a_sorted.bam; rm $folderName/a.sam $folderName/a_sorted.sam) &">> $COMMAND_SCRIPT
            done
       fi
    fi
done

# run command_script (.txt file saved in $OUTPUT_DIR)
bash $COMMAND_SCRIPT
