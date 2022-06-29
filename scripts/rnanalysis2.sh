#!/bin/bash

#SBATCH -J RNA Analysis

#SBATCH -n 16
#SBATCH --mem=8G
#SBATCH -t 12:00:00

# Create a SBATCH Job for all files
#SBATCH -o output/output.out

#SBATCH --mail-type=END
#SBATCH --mail-user=pranav_mahableshwarkar@brown.edu

echo "##### INSTRUCTIONS FOR PART TWO (SCRIPT THREE) ######
The format is python3 ./rna/3_suppa_formatting.py input_folder output_name filter fbtr_refseq.tsv
Step 2 produces iso_tmp.txt files for each replicate. These need to be merged in two different ways. 

The first set of 4 merges will just be merging individual samples with themselves. The second set of four
will merge x replicates of a sample with y replicates of the compared sample. 

Output_Name -> The first four should just be: m_control, f_control, m_delgi (for example)
            -> The second four should be the different comparisons you want to make.
            -> C Male v C Female, C Male v M Male, C Female v M Female, M Male v M Female
 
Filter -> To recognize which  folder to get input from, the scripts looks at the folder names.
       -> So to differentiate you need to make four folders. (One for each sample type.)
       -> 1. m_control 2. f_control 3. m_mutant 4. f_mutant (mutant can be replaced with your specific mutant name.)
       -> Make sure you put all replicates for each sample category in the respective folders.

Once you have completed this, comment out the following 'exit 1' and this instruction and run! 
#####          COMPLETED          ######"
exit 1

if [ $# -ne 4 ]; then
	echo $0: "Usage: 'sbatch preprocess.sh' or 'bash preprocess.sh'
            1) Data Input Directory (data in time2splice folder structure)
            2) Results Output Directory (results in time2splice folder structure)
            3) Reference Directory (look below to see specific files)
            4) SUPPA PATH (install from GitHub)

            This is Part Two of the RNA Analysis scripts. Please read the following before running."
        exit 1
fi

# Activate in the time2splice conda environment. 
module load anaconda/2020.02
source /gpfs/runtime/opt/anaconda/2020.02/etc/profile.d/conda.sh
conda activate time2splice

DATA_DIR=$1
RES_DIR=$2 
REF_DIR=$3
SUPPA_PATH=/users/pmahable/data/pmahable/SUPPA #$4

# Load the appropriate modules!
module load samtools/1.13
module load salmon/1.3.0
module load gcc/8.3

# Establish File Paths for all reference files. 
REFSEQ="${REF_DIR}/reference/fbtr_refseq.tsv" 
FBIOE="${REF_DIR}/reference/flybase_events/flybase.events.ioe"
GENEGTF="${REF_DIR}/reference/genes.gtf"
DMELGTF="${REF_DIR}/reference/dmel-all-r6.46.gtf"

# Checks to make sure that ALL reference files exist. 
if [ ! -e $REFSEQ ] 
then
    echo $0:"${REFSEQ} not found. Please make sure file path is correct."
    exit 1
elif [ ! -e $FBIOE ] 
then 
    echo $0:"${FBIOE} not found. Please make sure file path is correct."
    exit 1
elif [ ! -e $GENEGTF ] 
then 
    echo $0:"${GENEGTF} not found. Please make sure file path is correct."
    exit 1
elif [ ! -e $DMELGTF ] 
then 
    echo $0:"${DMELGTF} not found. Please make sure file path is correct."
    exit 1
fi
echo "All Reference Files Found. Starting..."

# USER INPUT: Highly recommended that m_ or f_ starts each one! Only need to edit Lines 90/91.
CONNAME="control" 
MUTNAME="delgi"
# These are the four comparisons in 3. 
CONMvF="control_m_f"
MCONvMUT="m_control_delgi"
FCONvMUT="f_control_delgi"
MUTMvF="delgi_m_f"

FFLYTYPE="L3F"
MFLYTYPE="L3M"

# 3. A Merge must be run for EVERY 4_suppa.sh data comparison that is run (folder required to run SUPPA).
# a) Control Males v Mutant Males. 
python3 ./rna/3_suppa_formatting.py ${SUPPARESULTS}/ ${SUPPARESULTS}/${MCONvMUT} m ${REFSEQ}
# b) Control Females v Mutant Females.
python3 ./rna/3_suppa_formatting.py ${SUPPARESULTS}/ ${SUPPARESULTS}/${FCONvMUT} f ${REFSEQ}
# c) Control Males v Control Females
python3 ./rna/3_suppa_formatting.py ${SUPPARESULTS}/ ${SUPPARESULTS}/${CONMvF} ${CONNAME} ${REFSEQ}
# d) Mutant Males v Mutant Females
python3 ./rna/3_suppa_formatting.py ${SUPPARESULTS}/ ${SUPPARESULTS}/${MUTMvF} ${MUTNAME} ${REFSEQ}

# 4. Run SUPPA Differential Splicing Analysis. For each run, there must be a corresponding merge call made. 
#    Need to make sure that input merged_iso_tpm file ONLY has replicate names in first row. Check the if group at 229 of this script if errors arise.
module load R/4.2.0 gcc/10.2 pcre2/10.35 intel/2020.2 texlive/2018
# a) Control Males v Mutant Males. 
./rna/4_suppa.sh ${SUPPARESULTS}/${MCONvMUT} ${FBIOE} $SUPPA_PATH ${GENEGTF} m_${CONNAME} m_${MUTNAME}
# b) Control Females v Mutant Females.
./rna/4_suppa.sh ${SUPPARESULTS}/${FCONvMUT} ${FBIOE} $SUPPA_PATH ${GENEGTF} f_${CONNAME} f_${MUTNAME}
# c) Control Males v Control Females
./rna/4_suppa.sh ${SUPPARESULTS}/${CONMvF} ${FBIOE} $SUPPA_PATH ${GENEGTF} f_${CONNAME} m_${CONNAME}
# d) Mutant Males v Mutant Females
./rna/4_suppa.sh ${SUPPARESULTS}/${MUTMvF} ${FBIOE} $SUPPA_PATH ${GENEGTF} f_${MUTNAME} m_${MUTNAME}

# 5. Create Pie Charts for Alternate Splicing Pattern Recognition
#    Had to pipinstall seaborn and matplotlib_venn (Install these in your conda environment if you run into issues.)
mkdir ${SUPPARESULTS}/5_PieCharts
python3 rna/5_calc_total_alt_splicing_piechart.py ${SUPPARESULTS}/${CONMvF}/f_${CONNAME}_iso_events.psi ${SUPPARESULTS}/5_PieCharts f_${CONNAME} >> ${SUPPARESULTS}/5_PieCharts/f_${CONNAME}_pie.txt
python3 rna/5_calc_total_alt_splicing_piechart.py ${SUPPARESULTS}/${CONMvF}/m_${CONNAME}_iso_events.psi ${SUPPARESULTS}/5_PieCharts m_${CONNAME} >> ${SUPPARESULTS}/5_PieCharts/m_${CONNAME}_pie.txt
python3 rna/5_calc_total_alt_splicing_piechart.py ${SUPPARESULTS}/${MUTMvF}/f_${MUTNAME}_iso_events.psi ${SUPPARESULTS}/5_PieCharts f_${MUTNAME} >> ${SUPPARESULTS}/5_PieCharts/f_${MUTNAME}_pie.txt
python3 rna/5_calc_total_alt_splicing_piechart.py ${SUPPARESULTS}/${MUTMvF}/m_${MUTNAME}_iso_events.psi ${SUPPARESULTS}/5_PieCharts m_${MUTNAME} >> ${SUPPARESULTS}/5_PieCharts/m_${MUTNAME}_pie.txt

# 6. Identify the Bias Genes. This is ONLY done using the merged iso_tpm files from Step 3 on Control Data. 
#    In this case it is ONLY Control Female v Control Male however this can change with multiple time points. 
#    This is used as a control in future experiments as these genes can be thought to be exclusive. 

# Needed to download dmel-all-r6.46.gtf from http://ftp.flybase.org/genomes/dmel/current/gtf/
# Automatic --threshold is 3. Use --threshold x to change the new threshold to x. 
# mkdir 6_control_bias_genes in the folder PRIOR to running if not done: mkdir ${SUPPARESULTS}/6_control_bias_genes
mkdir ${RES_DIR}/analysis/suppa/6_control_bias_genes
python3 rna/6_get_bias_genes.py ${SUPPARESULTS}/${CONMvF}/f_${CONNAME}_iso.tpm ${SUPPARESULTS}/${CONMvF}/m_${CONNAME}_iso.tpm ${SUPPARESULTS}/6_control_bias_genes ${FFLYTYPE} ${MFLYTYPE} ${DMELGTF}

# 7. Using the DPSI files produced from SUPPA, create plots to visualize the data. (i.e. volcano, violin, notch plots, etc.)
#    To do this part, you have to go through 7_plots_splicing.ipynb and manually edit everything and then run from there!
#    When I ran this, I had to manually go through the dPSI files and make sure that the column headers were identical across all samples. 

# 8. If you are using time-specific samples, run this step. I am not so I did not do this step. You will likely have to reformat the Notebook.

echo "Done"