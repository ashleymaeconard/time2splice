#!/bin/bash

#SBATCH -J 4Suppa

#SBATCH -n 4
#SBATCH --mem=8G
#SBATCH -t 1:00:00

# Create an array SBATCH Job for all files
#SBATCH -o output/output-1.out

#SBATCH --mail-type=END
#SBATCH --mail-user=pranav_mahableshwarkar@brown.edu

# Activate in the time2splice conda environment. 
module load anaconda/2020.02
source /gpfs/runtime/opt/anaconda/2020.02/etc/profile.d/conda.sh
conda activate time2splice


DATA_DIR=$1
RES_DIR=$2 
SUPPA_PATH=/users/pmahable/data/pmahable/SUPPA


# Load the appropriate modules!
module load samtools/1.13
module load salmon/1.3.0
module load gcc/8.3

# IMPORTANT: Step 1 needs to be modified based on the formatting of the trim galore fastq file (Line 89).

# 1. Run Salmon using the .sam files from Bowtie2 (These reads are already aligned to sequences on the transcriptome.)
# ./rna/1_run_salmon.sh ${RES_DIR}/preprocess/trim_galore_fastqc ${RES_DIR}/analysis/salmon ${DATA_DIR}/reference ${DATA_DIR}/reference/dmel-all-transcript-r6.46.fasta 16 4 1

# 2. Run SUPPA to isolate splicing differences within a sample type.
# ./rna/2_run_suppa.sh ${RES_DIR}/analysis/salmon ${RES_DIR}/analysis/suppa2 ${DATA_DIR}/reference/dmel-all-transcript-r6.46.fasta  1

# 3. Format the data from SUPPA. Had to reformat the folder structure (look at directory within 2_run_suppa). Then ran to merge data sets. 
#    A Merge must be run for EVERY 4_suppa.sh data comparison that is run (folder required to run SUPPA).
# a) Control Males v Mutant Males. 
# python3 ./rna/3_suppa_formatting.py ${RES_DIR}/analysis/suppa/2_run_suppa_output ${RES_DIR}/analysis/suppa/merged_clampvgfp_m KC /users/pmahable/data/pmahable/data/reference/fbtr_refseq.tsv
# b) Control Females v Mutant Females.
# python3 ./rna/3_suppa_formatting.py ${RES_DIR}/analysis/suppa/2_run_suppa_output ${RES_DIR}/analysis/suppa/merged_clampvgfp_f S2 /users/pmahable/data/pmahable/data/reference/fbtr_refseq.tsv
# c) Control Males v Control Females
# python3 ./rna/3_suppa_formatting.py ${RES_DIR}/analysis/suppa/2_run_suppa_output/GFP ${RES_DIR}/analysis/suppa/merged_gfp_fm b /users/pmahable/data/pmahable/data/reference/fbtr_refseq.tsv
# d) Mutant Males v Mutant Females
# python3 ./rna/3_suppa_formatting.py ${RES_DIR}/analysis/suppa/2_run_suppa_output/CLAMP ${RES_DIR}/analysis/suppa/merged_clamp_fm b /users/pmahable/data/pmahable/data/reference/fbtr_refseq.tsv

# 4. Run SUPPA Differential Splicing Analysis. For each run, there must be a corresponding merge call made. 
#    Need to make sure that input merged_iso_tpm file ONLY has replicate names in first row. Get Control/Treatment names by looking at the merged column names (this changed for me halfway through).
# module load R/4.2.0 gcc/10.2 pcre2/10.35 intel/2020.2 texlive/2018
# a) Control Males v Mutant Males. 
# ./rna/4_suppa.sh ${RES_DIR}/analysis/suppa/merged_clampvgfp_m ${DATA_DIR}/reference/flybase_events/flybase.events.ioe $SUPPA_PATH ${DATA_DIR}/reference/genes.gtf mt_KC_CLAMP mc_KC_GFP
# # b) Control Females v Mutant Females.
# ./rna/4_suppa.sh ${RES_DIR}/analysis/suppa/merged_clampvgfp_f ${DATA_DIR}/reference/flybase_events/flybase.events.ioe $SUPPA_PATH ${DATA_DIR}/reference/genes.gtf ft_S2_CLAMP fc_S2_GFP
# # c) Control Males v Control Females
# ./rna/4_suppa.sh ${RES_DIR}/analysis/suppa/merged_gfp_fm ${DATA_DIR}/reference/flybase_events/flybase.events.ioe $SUPPA_PATH ${DATA_DIR}/reference/genes.gtf bf_S2_GFP bm_KC_GFP
# d) Mutant Males v Mutant Females
# ./rna/4_suppa.sh ${RES_DIR}/analysis/suppa/merged_clamp_fm ${DATA_DIR}/reference/flybase_events/flybase.events.ioe $SUPPA_PATH ${DATA_DIR}/reference/genes.gtf bf_S2_CLAMP bm_KC_CLAMP

# 5. Create Pie Charts for Alternate Splicing Pattern Recognition
#    Had to pipinstall seaborn and matplotlib_venn (Install these in your conda environment if you run into issues.)
# python3 rna/5_calc_total_alt_splicing_piechart.py ${RES_DIR}/analysis/suppa/merged_clamp_fm/KC_CLAMP_iso_events.psi ${RES_DIR}/analysis/suppa/5_PieCharts KC_CLAMP >> ${RES_DIR}/analysis/suppa/5_PieCharts/KC_CLAMP_pie.txt
# python3 rna/5_calc_total_alt_splicing_piechart.py ${RES_DIR}/analysis/suppa/merged_clamp_fm/S2_CLAMP_iso_events.psi ${RES_DIR}/analysis/suppa/5_PieCharts S2_CLAMP >> ${RES_DIR}/analysis/suppa/5_PieCharts/S2_CLAMP_pie.txt
# python3 rna/5_calc_total_alt_splicing_piechart.py ${RES_DIR}/analysis/suppa/merged_gfp_fm/KC_GFP_iso_events.psi ${RES_DIR}/analysis/suppa/5_PieCharts KC_GFP >> ${RES_DIR}/analysis/suppa/5_PieCharts/KC_GFP_pie.txt
# python3 rna/5_calc_total_alt_splicing_piechart.py ${RES_DIR}/analysis/suppa/merged_gfp_fm/S2_GFP_iso_events.psi ${RES_DIR}/analysis/suppa/5_PieCharts S2_GFP >> ${RES_DIR}/analysis/suppa/5_PieCharts/S2_GFP_pie.txt

# 6. Identify the Bias Genes. This is ONLY done using the merged iso_tpm files from Step 3 on Control Data. 
#    In this case it is ONLY Control Female v Control Male however this can change with multiple time points. 
#    This is used as a control in future experiments as these genes can be thought to be exclusive. 

# Needed to download dmel-all-r6.46.gtf from http://ftp.flybase.org/genomes/dmel/current/gtf/
# Automatic --threshold is 3. Use --threshold x to change the new threshold to x. 
python3 rna/6_get_bias_genes.py ${RES_DIR}/analysis/suppa/merged_gfp_fm/KC_GFP_iso.tpm ${RES_DIR}/analysis/suppa/merged_gfp_fm/S2_GFP_iso.tpm ${RES_DIR}/analysis/suppa/6_control_bias_genes KC S2 ${DATA_DIR}/reference/dmel-all-r6.46.gtf
echo "Done"

conda deactivate
