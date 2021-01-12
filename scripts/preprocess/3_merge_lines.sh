#!/bin/bash
# merge_lines.sh
# Ashley Mae Conard
# Last Mod. April 16, 2020
# Purpose: Merges the fastq.gz files that result from different lanes of the same flow cell. 

# Check to make sure input is correct
if [ $# -ne 4 ]; then
	echo $0: Usage: ./merge_lines.sh NUM_COMMANDS \(to run at a time\) /PATH/TO/DIR/ \(with 2 subdirs - containing fastq.gz files\) /PATH/TO/OUTPUT_DIR/ \(suggest placing in /results/trim_galore_fastqc\) SEARCH_TERM \(e.g. _val_1.fq.gz, or *.fastq.gz*\) 
	exit 1
fi

# CORES_PER_COMMAND_RUN
NUM_COMMANDS=$1
INDIR=$2
OUTDIR=$3
SEARCH_TERM=$4

# Make NUM_COMMANDS even, because we process commands for R1 and R2 in the same loop
if (( ${NUM_COMMANDS}%2==1 )) && (( ! ${NUM_COMMANDS}==0 ))
then
    ((NUM_COMMANDS = NUM_COMMANDS - 1))
    echo "NUM_COMMANDS will be reduced by 1 to ${NUM_COMMANDS} to be even to process R1 and R2 in the same batch."
fi 
 
#Make output directory if necessary
mkdir -p $OUTDIR

# Create command script for Bowtie2
COMMAND_SCRIPT=${OUTDIR}/"commands_merge_lines.txt" 

# Remove command script if already made
if [ -f ${COMMAND_SCRIPT} ]; then
	echo -e "REPLACING ${COMMAND_SCRIPT}\n"
	rm -rf ${COMMAND_SCRIPT}
fi

# Remove odd folder delimeters 
# for i in *-ds*; do echo ${i};n=${i%-ds.*}; mv ${i} ${n}; done

START=0
i=$START

# Iterate through only the R1 replicates
echo "Initiating paired-end RNA-seq data line merging for R1 and R1 separately.";
FILES=$(find ${INDIR} -type f -name "*L001*R1*"${SEARCH_TERM}"*")
#echo ${FILES}

for R1 in ${FILES} 
    do
        #echo " "
        #echo "Getting R1 lines L2, L3, L3 matching paired .fastq $R1"
        # iterate through each R1 to determine when to add 'wait' to script
        ((i = i + 2)) # 2 because we have R1 and R2 separately

        # if number of processors reached, add wait to form a batch that will finish, and then process the next batch
        if (( ${i}%${NUM_COMMANDS}==0 )); then
            echo "wait" >> $COMMAND_SCRIPT
        fi 
        
        # Finding R1 L2,L3,L4
        #R1 for L2
        L2_folder=$(echo `dirname ${R1/"L001"/"L002"}`)
        to_find_R1L2=${L2_folder}"/*L002*R1*"${SEARCH_TERM}"*"
        R1L2=$(find ${L2_folder} -type f -wholename ${to_find_R1L2})
        
        #R1 for L3
        L3_folder=$(echo `dirname ${R1/"L001"/"L003"}`)
        to_find_R1L3=${L3_folder}"/*L003*R1*"${SEARCH_TERM}"*"
        R1L3=$(find ${L3_folder} -type f -wholename ${to_find_R1L3})
        
        # R1 for L4
        L4_folder=$(echo `dirname ${R1/"L001"/"L004"}`)
        to_find_R1L4=${L4_folder}"/*L004*R1*"${SEARCH_TERM}"*"
        R1L4=$(find ${L4_folder} -type f -wholename ${to_find_R1L4})
        
        if [[ ! -f ${R1} ]] || [[ ! -f ${R1L2} ]] || [[ ! -f ${R1L3} ]] || [[ ! -f ${R1L4} ]]
        then
          echo "One or more of R1 is missing."
          echo ${R1}, ${R1L2}, ${R1L3}, ${R1L4}
          exit
        fi
        
        # Finding R2
        #R2 for L2, L3, L4
        R2=${R1/"_R1_"/"_R2_"}
        R2v=${R2//val_1/val_2}
        SEARCH_TERM2=${SEARCH_TERM/"_1."/"_2."}

        #echo "Getting R2 lines L2, L3, L3 matching ${R2v}"
        #R2 for L2
        R2L2_folder=$(echo `dirname ${R2v/"L001"/"L002"}`)
        to_find_R2L2=${R2L2_folder}"/*L002*R2*"${SEARCH_TERM2}"*"
        R2L2=$(find ${R2L2_folder} -type f -wholename ${to_find_R2L2})
        
        #R2 for L3
        R2L3_folder=$(echo `dirname ${R2v/"L001"/"L003"}`)
        to_find_R2L3=${R2L3_folder}"/*L003*R2*"${SEARCH_TERM2}"*"
        R2L3=$(find ${R2L3_folder} -type f -wholename ${to_find_R2L3})
        
        #R2 for L4
        R2L4_folder=$(echo `dirname ${R2v/"L001"/"L004"}`)
        to_find_R2L4=${R2L4_folder}"/*L004*R2*"${SEARCH_TERM2}"*"
        R2L4=$(find ${R2L4_folder} -type f -wholename ${to_find_R2L4})
        
        if [[ ! -f ${R2v} ]] || [[ ! -f ${R2L2} ]] || [[ ! -f ${R2L3} ]] || [[ ! -f ${R2L4} ]]
        then
          echo "One or more of R2 is missing."
          echo ${R2v}, ${R2L2}, ${R2L3}, ${R2L4}
          exit
        fi
        
        # Output 
        final_file_name=$(echo `basename ${R1}` |cut -d$'_' -f1)

        # create folder for all outputs per file by removing R1 and R2 (e.g. MTb8-8)
        #fName=$(echo ${fileName%$SUFFIX_TO_REMOVE}) 
        outdirfolder=${OUTDIR}/${final_file_name}

        # output directory for BOWTIE2/sample_name
        mkdir -p ${outdirfolder}           
        
        echo "(zcat ${R1} ${R1L2} ${R1L3} ${R1L4} | gzip > ${outdirfolder}/${final_file_name}_R1.fq.gz) &">> $COMMAND_SCRIPT
        echo "(zcat ${R2v} ${R2L2} ${R2L3} ${R2L4} | gzip > ${outdirfolder}/${final_file_name}_R2.fq.gz) &">> $COMMAND_SCRIPT
done
echo "Check ${COMMAND_SCRIPT} for commands."
bash $COMMAND_SCRIPT      

