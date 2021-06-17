<img src="https://github.com/ashleymaeconard/time2splice/blob/master/img/time2splice_logo.png" width="100%">

### A method to identify temporal and sex-specific alternative splicing from multi-omics data.

####  Author: Ashley Mae Conard

Motivation
==========
Alternative splicing can occur in at least 3/4th of human genes to encode two or more splice isoforms. These isoforms occur in different proportions over time, and between sexes. Thus, we present a method to characterize these isoforms, so to better understand gene regulation happening in normal and diseased states. Time2splice identifies temporal and sex-specific alternative splicing combinding multi-omic (i.e. both expression via RNA-seq, and protein-DNA interaction via CUT&RUN and ChIP-seq) data. Analysis is done in 3 parts. 1) Temporal expression analysis, 2) Temporal protein-DNA analysis, and 3) Temporal multi-omics integration.

**NOTE: Snakemake pipeline coming soon for time2splice!**

Paper
==========
Ray*, M., Conard*, A. M., Urban, J., & Larschan, E. (2021). Sex-specific transcript diversity is regulated by a maternal pioneer factor in early Drosophila embryos. bioRxiv.

*contributed equally

Outline
==========

  1. [Preprocess](#preprocess-scriptspreprocess)
  2. [Temporal expression analysis](#temporal-expression-analysis-scriptsrna)
  3. [Temporal protein-DNA analysis](#temporal-protein-dna-analysis-scriptsprotein_dna)
  4. [Temporal multi-omics integration](#temporal-multi-omics-integration-scriptsmultio_analysis)

Preprocess (scripts/preprocess)
==========

### Retrieve raw data, quality control, trimming, alignment. Perform steps as needed.

`1_parse_sraRunTable.sh`

Creates `time2splice/` folder structure, as well as `metadatafile.csv` and `SraAccList.txt` (which is needed for next command to get .fastq files).

`1_get_fastq_files.sh`

Retrieves .fastq files by passing in `SraAccList.txt` from aforementioned step.

`2_run_fastQC.sh`

Runs FastQC for all .fastq files in a given directory.

`3_run_trim_galore.sh`

Run Trim Galore! followed by FastQC to trim any reads below quality threshold.

`3_merge_lines.sh`

Merges all the different lanes of the same flow cell .fastq files. 

`4_run_Bowtie2.sh` or `preprocess/4_run_BWA.sh` or `preprocess/4_run_HISAT2.sh`. 

Runs one or more of these three aligners (Bowtie2, BWA, or HISAT2) on .fastq data in a given directory. 

`5_plot_alignment.py`

Plot the alignments from either one or two different aligners (Bowtie2 or HISAT2).

Temporal expression analysis (scripts/rna)
==========

`1_run_salmon.sh`

Run salmon to quantify transcript expression for treatment and control samples.

e.g. `./1_run_salmon.sh /nbu/compbio/aconard/larschan_data/sexed_embryo/ /data/compbio/aconard/splicing/results/salmon_results_ncbi_trans/ /data/compbio/aconard/BDGP6/transcriptome_dir/pub/infphilo/hisat2/data/bdgp6_tran/genome.fa 3 10 1 _001.fastq.gz` 

`2_run_suppa.sh`

Run Suppa for treatment and control samples.

e.g. `./2_run_suppa.sh /data/compbio/aconard/splicing/results/salmon_results/ /data/compbio/aconard/splicing/results/suppa_results_ncbi_trans/ /data/compbio/aconard/BDGP6/transcriptome_dir/pub/infphilo/hisat2/data/bdgp6_tran/genome.fa 20`

`3_suppa_formatting.py`

Converts NM_ gene names to flybase name, then merging outputs from run_suppa (NM_ gene names by 1 TPM value column for each replicate)

`4_suppa.sh`

Identifies various forms of differential splicing (e.g. using PSI and DTU)

`5_calc_total_alt_splicing_controls.py`

Calculate and plot the proportions of alternative splicing (in pie chart) in control samples.

`6_calc_total_alt_differential_splicing.py`

Calculate and plot the proportions of alternative splicing (in pie chart) in treatment samples.

`7_get_bias_genes.py`

Retrieve male and female biased genes and create bed files for average profile plotting.

`8_plots_splicing.ipynb`

Plotting transcript expression using PSI and DTU measures.

`8_alt_plots_splicing.ipynb`

Alternative code base to plot transcript expression using PSI and DTU measures.

`9_plots_splicing_time.ipynb`

Plot alternative splicing genes within categories (all females, all males, females sex specific, male sex specific, female all rest, male all rest, female non-sex specific, male non-sex specific, female new sex specific, male new sex specific) over time.

Temporal protein-DNA analysis (scripts/protein_dna)
==========

`1_run_picard_markduplicates.sh`

Run Picard's MarkDuplicates in for all .sorted.bam files in a given directory.

`2_run_macs2.sh`

Runs MACS2 to call peaks for all .sorted.bam files in a given directory.

`3_run_macs2_fold_enrich.sh`

Generate signal track using MACS2 to profile transcription factor modification enrichment levels genome-wide.

Temporal multi-omics integration (scripts/multio_analysis)
==========

Note, there is no order to these scripts. Each analysis / results exploration is independent. More analysis scripts to come.

`overlap_protein_DNA_peaks.sh`

Runs Intervene to view intersection of each narrowpeak file.

`histogram_peak_val_intensity.ipynb`

Plot peak intensity for a given narrow peak file.

`get_coord_run_meme.sh`

Get coordinates of bed file and run through MEME.

`alt_splicing_chi_squared.ipynb`

Perform chi-squared test on alternative splicing categories. Mutually Exclusive Exons (MXE) used in this example.



