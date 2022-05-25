# get_coord_run_meme.sh
# Ashley Mae Conard
# Last Mod. Dec 1, 2020
# Purpose: Get coordinates of bed file and run through MEME.

# Check to make sure input is correct
if [ $# -ne 4 ]; then
	echo $0: "Usage: ./get_fastas_from_bed.sh 
                1) /FULL/PATH/TO/TF_OVERLAP_FILE (req first 3 cols: ^chr\tstart\tend) 
                2) OUTPUT_NAME (e.g. mle.clampR.v.clamp.X.2.4.M )
                3) OUTDIR 
                4) MAXSIZE"
	exit 1
fi

# Assign input to variable names
TF_OVERLAP_FILE=$1
OUTPUT_NAME=$2
OUTDIR=$3
MAXSIZE=$4

FASTA_FILE="/data/compbio/aconard/genomes/BDGP6/Drosophila_melanogaster/Ensembl/BDGP6/Sequence/WholeGenomeFasta/genome.fa"

# Make directory if does not exist
mkdir -p $OUTDIR/$OUTPUT_NAME/

# Create bed file from FILE. Remove chr if exists and remove top of file where it says chr start end
cat ${TF_OVERLAP_FILE} | cut -d$'\t' -f1,2,3 | sed 's/chr//' | grep -v seq > $OUTDIR/$OUTPUT_NAME/$OUTPUT_NAME.bed

# Determine maxsize
set_maxSize=$(< $OUTDIR/$OUTPUT_NAME/$OUTPUT_NAME.bed wc -l)

# Get sequences for each bed file coordinate (one pre line)
bedtools getfasta -fi ${FASTA_FILE} -bed $OUTDIR/$OUTPUT_NAME/$OUTPUT_NAME.bed > $OUTDIR/$OUTPUT_NAME/$OUTPUT_NAME.seq

# Run MEME to find motifs de novo (could add -maxw 12)
(meme $OUTDIR/$OUTPUT_NAME/$OUTPUT_NAME.seq -dna -nmotifs 2 -maxsize ${MAXSIZE} -mod anr -oc $OUTDIR/$OUTPUT_NAME/) &

echo "FINISHED MEME. Results for $OUTPUT_NAME are here: $OUTDIR/$OUTPUT_NAME/"
