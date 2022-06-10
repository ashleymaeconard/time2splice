# 7_get_bias_genes.py
# Ashley Mae Conard
# Last Mod. 6/17/2020
# Purpose: get category 1 (e.g. female) and category 2 (e.g. male) biased genes, plot Venn, and create bed files. 
#          .bed files used for average profile plotting downstream.

# Libraries
import argparse, warnings, glob, sys, os
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib_venn import venn2, venn2_circles
warnings.filterwarnings('ignore')
#search_word='gene_id'

def create_cols(row, searchw):
    string = row['allInfo']
    list_of_words = string.split()
    next_word = list_of_words[list_of_words.index(searchw) + 1]
    return(next_word)

def get_bias_genes(INPUTFILE_1, INPUTFILE_2, OUTDIR, OUT_FILE_NAME_1, OUT_FILE_NAME_2, DM_FILE, THRESH):

    # Import, filter, and create dataframe of gene and transcript info
    df_dm_gtf = pd.read_csv(DM_FILE, sep="\t",usecols=[0,2,3,4,8], names=['chrom','type_prot', 'start_chrom', 'end_chrom', 'allInfo'], engine='python')
    df_dm_gtf = df_dm_gtf[df_dm_gtf['allInfo'].str.contains("transcript")] # keep lines with transcript info
    #df_dm_gtf['allInfo'] = df_dm_gtf['allInfo']
    df_dm_gtf.columns.str.replace(' ', '')
    list_cols =["gene_id", "gene_symbol", "transcript_id"]

    # Iterate through all gene information (from allInfo in df_genes_gtf of genes.gtf) to keep
    for search_word in list_cols:
        df_dm_gtf[search_word] = df_dm_gtf.apply(create_cols,args=(search_word,), axis=1)
        df_dm_gtf[search_word] = df_dm_gtf[search_word].str.replace(';','')
        df_dm_gtf[search_word] = df_dm_gtf[search_word].str.replace('"', '')

    df_dm_gtf = df_dm_gtf.drop(columns=['allInfo'])
    df_dm_gtf = df_dm_gtf.drop_duplicates(subset="transcript_id")

    # Create dictionary of gtf for quick .bed info. retrieval 
    a_dict = {}
    for idx, row in df_dm_gtf.iterrows():
        a_dict[row["transcript_id"]] = (row['chrom'], row['start_chrom'], row['end_chrom'], row["gene_id"], row["gene_symbol"])

    # Read in .csvs
    df_fbg_0_2 = pd.read_csv(INPUTFILE_1, sep="\t")
    df_fbg_0_2 = df_fbg_0_2[df_fbg_0_2.apply(lambda x : np.mean(x.values) >= THRESH, axis=1)]

    df_mbg_0_2 = pd.read_csv(INPUTFILE_2, sep="\t")
    df_mbg_0_2 = df_mbg_0_2[df_mbg_0_2.apply(lambda x : np.mean(x.values) >= THRESH, axis=1)]


    # Get values to create .bed file (from allInfo in df_genes_gtf of .gtf)
    # chrom
    df_fbg_0_2["chrom"] = df_fbg_0_2.apply(lambda x : a_dict[x.name][0], axis=1)
    df_mbg_0_2["chrom"] = df_mbg_0_2.apply(lambda x : a_dict[x.name][0], axis=1)

    # start
    df_fbg_0_2["start"] = df_fbg_0_2.apply(lambda x : a_dict[x.name][1], axis=1)
    df_mbg_0_2["start"] = df_mbg_0_2.apply(lambda x : a_dict[x.name][1], axis=1)

    # end
    df_fbg_0_2["end"] = df_fbg_0_2.apply(lambda x : a_dict[x.name][2], axis=1)
    df_mbg_0_2["end"] = df_mbg_0_2.apply(lambda x : a_dict[x.name][2], axis=1)

    # gene if
    df_fbg_0_2["gene_id"] = df_fbg_0_2.apply(lambda x : a_dict[x.name][3], axis=1)
    df_mbg_0_2["gene_id"] = df_mbg_0_2.apply(lambda x : a_dict[x.name][3], axis=1)

    # gene name
    df_fbg_0_2["gene_name"] = df_fbg_0_2.apply(lambda x : a_dict[x.name][4], axis=1)
    df_mbg_0_2["gene_name"] = df_mbg_0_2.apply(lambda x : a_dict[x.name][4], axis=1)

    # strand
    df_fbg_0_2["strand"] = "-"
    df_mbg_0_2["strand"] = "-"


    # Plot category overlaps
    venn2([set(list(df_fbg_0_2.index)), set(list(df_mbg_0_2.index))],
    set_labels = (OUT_FILE_NAME_1, OUT_FILE_NAME_2))
    plt.title('Transcript overlaps between %s and %s\n'%(OUT_FILE_NAME_1, OUT_FILE_NAME_2))
    plt.savefig(OUTDIR+'/venn_cat1_%s_cat2_%s.png'%(OUT_FILE_NAME_1, OUT_FILE_NAME_2))
    
    # Get differences between category 1 and category 2 (which genes are only on in category 1 and category 2 - set difference)
    fb_0_2_trans = set.difference(set(df_fbg_0_2.index), set(df_mbg_0_2.index))
    print("Set difference category 1: ", len(fb_0_2_trans))
    mb_0_2_trans = set.difference(set(df_mbg_0_2.index), set(df_fbg_0_2.index))
    print("Set difference category 2: ", len(mb_0_2_trans))


    # Save bias gene lists as bed files
    bed_fbg_0_2 = df_fbg_0_2[df_fbg_0_2.index.isin(list(fb_0_2_trans))].drop_duplicates(subset="gene_id")[["chrom","start","end","gene_id","gene_name","strand"]]
    bed_fbg_0_2.to_csv(OUTDIR+"/"+ OUT_FILE_NAME_1 +".bed", sep="\t", index=False, header=False)

    bed_mbg_0_2 = df_mbg_0_2[df_mbg_0_2.index.isin(list(mb_0_2_trans))].drop_duplicates(subset="gene_id")[["chrom","start","end","gene_id","gene_name","strand"]]
    bed_mbg_0_2.to_csv(OUTDIR+"/"+ OUT_FILE_NAME_2 +".bed", sep="\t", index=False, header=False)

    
    # See the number of genes vs. transcripts when dropping gene ID duplicates
    print("Category 1 number of transcripts: ", len(bed_fbg_0_2))
    print("Category 1 number of genes: ", len(bed_fbg_0_2.drop_duplicates(subset="gene_id")))

    print("Category 1 number of transcripts: ", len(bed_mbg_0_2))
    print("Category 1 number of genes: ", len(bed_mbg_0_2.drop_duplicates(subset="gene_id")))


    # Plot to choose TPM threshold
    dfs = [df_fbg_0_2, df_mbg_0_2]
    df_list = [OUT_FILE_NAME_1, OUT_FILE_NAME_2]

    cnt = 0
    for d, dl in zip(dfs, df_list):
        cnt+=1
        hist_nums = list()
        thresh_list = list()

        for thresh_toggle in range(15):
            df = d.select_dtypes(exclude=['object'])
            num_genes = len(df[df.apply(lambda x : np.mean(x.values) >= thresh_toggle, axis=1)])
            thresh_list.append(thresh_toggle)
            hist_nums.append(num_genes)
            
        fig = plt.figure()
        ax = plt.axes()
        plt.title(dl)
        ax.plot(thresh_list, hist_nums,marker=".")
        plt.savefig(OUTDIR+"/"+OUT_FILE_NAME_1+"_"+OUT_FILE_NAME_2+"line_threshold_plot_%s.png"%(str(cnt)))
        plt.close()
        

def main(args):
    INPUTFILE_1 = args.input_file_cat_1
    INPUTFILE_2 = args.input_file_cat_2
    OUTDIR = args.output_dir 
    OUT_FILE_NAME_1 = args.output_file_name_1
    OUT_FILE_NAME_2 = args.output_file_name_2
    DM_FILE = args.dm_file
    THRESH = args.threshold

    get_bias_genes(INPUTFILE_1, INPUTFILE_2, OUTDIR, OUT_FILE_NAME_1, OUT_FILE_NAME_2, DM_FILE, THRESH)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find bias genes using control samples. <!> NOTE: Run this for every group pairs (2 only) comparison separately (e.g. females vs. males, females time 1 vs. females time 2, etc.) <!>")
    parser.add_argument("input_file_cat_1", type=str, help="Import control tpm files for each sex and timepoint. (e.g. ./results/suppa_results_ncbi_trans/merged_0-2_f/control_0-2_f_iso.tpm)") 
    parser.add_argument("input_file_cat_2", type=str, help="Import control tpm files for each sex and timepoint. (e.g. ./results/suppa_results_ncbi_trans/merged_0-2_m/control_0-2_m_iso.tpm)") 
    parser.add_argument("output_dir", type=str, help="Choose output directory location (e.g. ./results/gene_groups/)")
    parser.add_argument("output_file_name_1", type=str, help="Choose category 1 output file name.") 
    parser.add_argument("output_file_name_2", type=str, help="Choose category 2 output file name.") 
    parser.add_argument("dm_file", type=str, help="Use .gtf file for mapping transcript to gene id and name (along with importing info for .bed creation), (e.g. ./genomes/BDGP6/dmel-all-r6.29.gtf).")
    
    parser.add_argument("--threshold", type=int, default=3, help="Filter each sample to include transcripts where the mean value is less than or equal to X TPMs per gene. (recommend 3 or so)")
    
    args = parser.parse_args()
    
    main(args)

#### Notes ####
# Use .gtf file for mapping transcript to gene id and name (along with importing info for 
#  .bed creation)
# NOTE - used dmel-all-r6.29.gtf instead of genes.gtf because some transcripts not found 
#  in genes.gtf! e.g. df_dm_gtf[df_dm_gtf["transcript_id"]=="FBtr0445398"] not in genes.gtf

# Paper references:
# Import control tpm files for each sex and timepoint
#fbg_0_2 = "/data/compbio/aconard/splicing_pj/results/suppa_results_ncbi_trans/merged_0-2_f/control_0-2_f_iso.tpm"
#mbg_0_2 = "/data/compbio/aconard/splicing_pj/results/suppa_results_ncbi_trans/merged_0-2_m/control_0-2_m_iso.tpm"
#fbg_2_4 = "/data/compbio/aconard/splicing_pj/results/suppa_results_ncbi_trans/merged_2-4_f/control_iso.tpm"
#mbg_2_4 = "/data/compbio/aconard/splicing_pj/results/suppa_results_ncbi_trans/merged_2-4_m/control_iso.tpm"
#DM_GTF = "/data/compbio/aconard/genomes/BDGP6/dmel-all-r6.29.gtf"
#GENES_GTF = "/data/compbio/aconard/genomes/BDGP6/genes.gtf"