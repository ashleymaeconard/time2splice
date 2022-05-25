#!/bin/bash
# overlap_protein_DNA_peaks.sh
# Ashley Mae Conard
# Last Mod. 07/2019
# Purpose: Runs Intervene to view intersection of each narrowpeak file

# Check to make sure input is correct
if [ $# -eq 0 ]; then
	echo 'Usage: ./overlap_cutNrun_chip.sh OUTPUT_DIR PLOT_TITLE /PATH/TO/FIRST_narrowPeak/ /PATH/TO/SECOND_narrowPeak/ /PATH/TO/THIRD_narrowPeak/ NOTE: underscore in title'
	exit 1
fi

# Joining array function to join values into comma separated list (for --names param)
function join_by { local d=$1; shift; echo -n "$1"; shift; printf "%s" "${@/#/$d}"; }

echo "
Running intervene to produce 5 different kinds of outputs.
"
input=$@
START=0
i=$START
declare -a arr
declare -a name_arr

for file in $input
do
    # Get output directory
    ((i = i + 1))
    if (( ${i}==1 )); then
        OUTPUT_DIR=${file}
        echo "OUTPUT_DIR" $OUTPUT_DIR 
    elif (( ${i}==2 )); then
        TITLE=${file}
    else
        # Get input bed/narrowPeak files
        echo input file: "$file"
        arr+=($file)
        
        fileRep=$(echo `basename $(dirname $file)`)
        fileSample=$(echo `basename $(dirname $(dirname $file))`)
        echo name: $fileSample\_$fileRep 
	name_arr+=($fileRep)
	#name_arr+=($fileSample\_$fileRep)
    fi
done

names=$( IFS=$','; echo "${name_arr[*]}" )
echo "
Creating intersection plots for ${arr[@]}
"

# # Venn
mkdir -p "$OUTPUT_DIR/$TITLE/venn/"
(intervene venn -i ${arr[@]} --names=${names[@]} --save-overlaps --title $TITLE -o "$OUTPUT_DIR/$TITLE/venn/") &
if [ ! "(ls -A $OUTPUT_DIR/$TITLE/venn/)" ]; then
    rm -rf "$OUTPUT_DIR/$TITLE/venn/"
fi

# # Bar Plot
mkdir -p "$OUTPUT_DIR/$TITLE/barPlot/"
(intervene upset -i ${arr[@]} --names=${names[@]} --save-overlaps -o "$OUTPUT_DIR/$TITLE/barPlot/") &

# # Pairwise pie
mkdir -p "$OUTPUT_DIR/$TITLE/pieChart/"
(intervene pairwise -i ${arr[@]} --names=${names[@]} --title $TITLE --htype pie -o "$OUTPUT_DIR/$TITLE/pieChart/") &

# # Pairwise heatmap
mkdir -p "$OUTPUT_DIR/$TITLE/heatmap/"
(intervene pairwise -i ${arr[@]} --names=${names[@]} --title $TITLE --htype color -o "$OUTPUT_DIR/$TITLE/heatmap/") &

# # # Pairwise tribar
mkdir -p "$OUTPUT_DIR/$TITLE/tribar/"
(intervene pairwise -i ${arr[@]} --names=${names[@]} --title $TITLE --htype tribar -o "$OUTPUT_DIR/$TITLE/tribar/") &

echo "
Results in $OUTPUT_DIR/$TITLE/venn/, barPlot/, pieChart/, heatmap/, and tribar/.
"

# e.g. "./overlap_cutNrun_chip.sh /data/compbio/aconard/overlap_cut_chip/ "Clamp_CUTnRUN_vs_ChIP_S2" ../cut-n-run/results/macs2/clamp_cell_line_m/MR-1_/out.sorted.bam_peaks.narrowPeak ../cut-n-run/results/macs2/clamp_cell_line_m/MR-2_/out.sorted.bam_peaks.narrowPeak ../cut-n-run/results/macs2/clamp_cell_line_m/MR-6_/out.sorted.bam_peaks.narrowPeak"


