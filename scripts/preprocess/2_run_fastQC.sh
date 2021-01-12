#!/bin/bash
# run_fastQC.sh
# Ashley Mae Conard
# Last Mod. Nov 5, 2018
# Purpose: Runs FastQC for all .fastq.gz files in a given directory

# Check to make sure input is correct
if [ $# -ne 3 ]; then
	echo $0: Usage: ./run_fastQC.sh NUM_THREADS /PATH/TO/DIR/ \(containing fastq.gz files\) /PATH/TO/OUTPUT_DIR/ \(suggest placing in /results/fastQC\)
	exit 1
fi

NUM_THREADS=$1
INDIR=$2
OUTDIR=$3

#Make output directory if necessary
mkdir -p $OUTDIR

#CREATE FASTQC DIR

# Move into directory and create fastQC script

find $INDIR -type f -name "*.fastq*"> $INDIR/file.txt # grab all .fastq.gz files in DIR
sed -i "1d" $INDIR/file.txt # remove dir itself from file.txt
sed -i "$ d" $INDIR/file.txt # remove "file.txt" from file.txt
sed -i "1 i\--outdir=$OUTDIR :::" $INDIR/file.txt
sed -i '1 i\fastqc {1}\n' $INDIR/file.txt  
sed -i "1 i\parallel -j $NUM_THREADS" $INDIR/file.txt

echo $INDIR/file.txt
tr '\n' ' ' < $INDIR/file.txt > $OUTDIR/commands_fastqc.txt # replace new line with space
rm -rf $INDIR/file.txt
bash $OUTDIR/commands_fastqc.txt
#rm -rf $OUTDIR/commands_fastqc.txt
