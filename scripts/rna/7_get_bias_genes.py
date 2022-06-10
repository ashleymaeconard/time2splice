# get_bias_genes.py
# Ashley Mae Conard
# Last Mod. 6/17/2020
# Purpose: get male and female biased genes and create bed files for average profile plotting 

# Libraries
#get_ipython().run_line_magic('matplotlib', 'inline')
import seaborn as sns
import pandas as pd
import numpy as np
import seaborn as sns
import glob, sys, os
import matplotlib.pyplot as plt 
#mpl.use('Agg')
from matplotlib_venn import venn2, venn2_circles
import warnings
warnings.filterwarnings('ignore')

# Import necessary files and folders
OUTDIR = "/data/compbio/aconard/"
#OUTDIR = "/data/compbio/aconard/cutnrun_v_chip/embryo_cut-n-run/results/gene_groups/"

# Import control tpm files for each sex and timepoint
fbg_0_2 = "/data/compbio/aconard/splicing_pj/results/suppa_results_ncbi_trans/merged_0-2_f/control_0-2_f_iso.tpm"
mbg_0_2 = "/data/compbio/aconard/splicing_pj/results/suppa_results_ncbi_trans/merged_0-2_m/control_0-2_m_iso.tpm"
fbg_2_4 = "/data/compbio/aconard/splicing_pj/results/suppa_results_ncbi_trans/merged_2-4_f/control_iso.tpm"
mbg_2_4 = "/data/compbio/aconard/splicing_pj/results/suppa_results_ncbi_trans/merged_2-4_m/control_iso.tpm"

# Use .gtf file for mapping transcript to gene id and name (along with importing info for 
#  .bed creation)
# NOTE - used dmel-all-r6.29.gtf instead of genes.gtf because some transcripts not found 
#  in genes.gtf! e.g. df_dm_gtf[df_dm_gtf["transcript_id"]=="FBtr0445398"] not in genes.gtf

DM_GTF = "/data/compbio/aconard/genomes/BDGP6/dmel-all-r6.29.gtf"
# GENES_GTF = "/data/compbio/aconard/genomes/BDGP6/genes.gtf"

# Import, filter, and create dataframe of gene and transcript info
df_dm_gtf = pd.read_csv(DM_GTF, sep="\t",usecols=[0,2,3,4,8], names=['chrom','type_prot', 'start_chrom', 'end_chrom', 'allInfo'], engine='python')
df_dm_gtf = df_dm_gtf[df_dm_gtf['allInfo'].str.contains("transcript")] # keep lines with transcript info
#df_dm_gtf['allInfo'] = df_dm_gtf['allInfo']
df_dm_gtf.columns.str.replace(' ', '')
list_cols =["gene_id", "gene_symbol", "transcript_id"]

def create_cols(row):
    string = row['allInfo']
    list_of_words = string.split()
    next_word = list_of_words[list_of_words.index(search_word) + 1]   
    return(next_word)

# iterate through all gene information (from allInfo in df_genes_gtf of genes.gtf) to keep
for search_word in list_cols:
    df_dm_gtf[search_word] = df_dm_gtf.apply(create_cols,axis=1)
    df_dm_gtf[search_word] = df_dm_gtf[search_word].str.replace(';','')
    df_dm_gtf[search_word] = df_dm_gtf[search_word].str.replace('"', '')

df_dm_gtf = df_dm_gtf.drop(columns=['allInfo'])
df_dm_gtf = df_dm_gtf.drop_duplicates(subset="transcript_id")

# Create dictionary of gtf for quick .bed info. retrieval 
a_dict = {}
for idx, row in df_dm_gtf.iterrows():
    a_dict[row["transcript_id"]] = (row['chrom'], row['start_chrom'], row['end_chrom'], row["gene_id"], row["gene_symbol"])

# Filter each sample to include transcripts where the mean value is less than or equal to 3 TPMs per gene.
threshold = 3

df_fbg_0_2 = pd.read_csv(fbg_0_2, sep="\t")
print(len(df_fbg_0_2))
df_fbg_0_2 = df_fbg_0_2[df_fbg_0_2.apply(lambda x : np.mean(x.values) >= threshold, axis=1)]
print(len(df_fbg_0_2))

df_mbg_0_2 = pd.read_csv(mbg_0_2, sep="\t")
print(len(df_mbg_0_2))
df_mbg_0_2 = df_mbg_0_2[df_mbg_0_2.apply(lambda x : np.mean(x.values) >= threshold, axis=1)]
print(len(df_mbg_0_2))

df_fbg_2_4 = pd.read_csv(fbg_2_4, sep="\t")
print(len(df_fbg_2_4))
df_fbg_2_4 = df_fbg_2_4[df_fbg_2_4.apply(lambda x : np.mean(x.values) >= threshold, axis=1)]
print(len(df_fbg_2_4))

df_mbg_2_4 = pd.read_csv(mbg_2_4, sep="\t")
print(len(df_mbg_2_4))
df_mbg_2_4 = df_mbg_2_4[df_mbg_2_4.apply(lambda x : np.mean(x.values) >= threshold, axis=1)]
print(len(df_mbg_2_4))


# Get values to create .bed file (from allInfo in df_genes_gtf of .gtf)

# chrom
df_fbg_0_2["chrom"] = df_fbg_0_2.apply(lambda x : a_dict[x.name][0], axis=1)
df_mbg_0_2["chrom"] = df_mbg_0_2.apply(lambda x : a_dict[x.name][0], axis=1)

df_fbg_2_4["chrom"] = df_fbg_2_4.apply(lambda x : a_dict[x.name][0], axis=1)
df_mbg_2_4["chrom"] = df_mbg_2_4.apply(lambda x : a_dict[x.name][0], axis=1)

# start
df_fbg_0_2["start"] = df_fbg_0_2.apply(lambda x : a_dict[x.name][1], axis=1)
df_mbg_0_2["start"] = df_mbg_0_2.apply(lambda x : a_dict[x.name][1], axis=1)

df_fbg_2_4["start"] = df_fbg_2_4.apply(lambda x : a_dict[x.name][1], axis=1)
df_mbg_2_4["start"] = df_mbg_2_4.apply(lambda x : a_dict[x.name][1], axis=1)

# end
df_fbg_0_2["end"] = df_fbg_0_2.apply(lambda x : a_dict[x.name][2], axis=1)
df_mbg_0_2["end"] = df_mbg_0_2.apply(lambda x : a_dict[x.name][2], axis=1)

df_fbg_2_4["end"] = df_fbg_2_4.apply(lambda x : a_dict[x.name][2], axis=1)
df_mbg_2_4["end"] = df_mbg_2_4.apply(lambda x : a_dict[x.name][2], axis=1)

# gene if
df_fbg_0_2["gene_id"] = df_fbg_0_2.apply(lambda x : a_dict[x.name][3], axis=1)
df_mbg_0_2["gene_id"] = df_mbg_0_2.apply(lambda x : a_dict[x.name][3], axis=1)

df_fbg_2_4["gene_id"] = df_fbg_2_4.apply(lambda x : a_dict[x.name][3], axis=1)
df_mbg_2_4["gene_id"] = df_mbg_2_4.apply(lambda x : a_dict[x.name][3], axis=1)

# gene name
df_fbg_0_2["gene_name"] = df_fbg_0_2.apply(lambda x : a_dict[x.name][4], axis=1)
df_mbg_0_2["gene_name"] = df_mbg_0_2.apply(lambda x : a_dict[x.name][4], axis=1)

df_fbg_2_4["gene_name"] = df_fbg_2_4.apply(lambda x : a_dict[x.name][4], axis=1)
df_mbg_2_4["gene_name"] = df_mbg_2_4.apply(lambda x : a_dict[x.name][4], axis=1)

# strand
df_fbg_0_2["strand"] = "-"
df_mbg_0_2["strand"] = "-"

df_fbg_2_4["strand"] = "-"
df_mbg_2_4["strand"] = "-"

# Plot overlap transcripts
dfs = [df_fbg_0_2, df_mbg_0_2, df_fbg_2_4, df_mbg_2_4]
df_list = ["df_fbg_0_2", "df_mbg_0_2", "df_fbg_2_4", "df_mbg_2_4"]

#0-2 h
venn2([set(list(df_fbg_0_2.index)), set(list(df_mbg_0_2.index))],
set_labels = ('f_trans_0_2', 'm_trans_0_2'))
plt.title('Transcript overlaps 0-2 hours\n')
plt.savefig(OUTDIR+'/test1.png')

#2-4 h
venn2([set(list(df_fbg_2_4.index)), set(list(df_mbg_2_4.index))],
set_labels = ('f_trans_2_4', 'm_trans_2_4'))
plt.title('Transcript overlaps 2-4 hours\n')
plt.savefig(OUTDIR+'/test2.png')


# Get differences between males and females 0 - 2 hours (which genes are only on in females and males - set difference)
fb_0_2_trans = set.difference(set(df_fbg_0_2.index), set(df_mbg_0_2.index))
print(len(fb_0_2_trans))
mb_0_2_trans = set.difference(set(df_mbg_0_2.index), set(df_fbg_0_2.index))
print(len(mb_0_2_trans))

fb_2_4_trans = set.difference(set(df_fbg_2_4.index), set(df_mbg_2_4.index))
print(len(fb_2_4_trans))
mb_2_4_trans = set.difference(set(df_mbg_2_4.index), set(df_fbg_2_4.index))
print(len(mb_2_4_trans))


# Save female and male bias gene lists as bed files

bed_fbg_0_2 = df_fbg_0_2[df_fbg_0_2.index.isin(list(fb_0_2_trans))].drop_duplicates(subset="gene_id")[["chrom","start","end","gene_id","gene_name","strand"]]
bed_fbg_0_2.to_csv(OUTDIR+"f_bias_genes_0_2.bed", sep="\t", index=False, header=False)

bed_mbg_0_2 = df_mbg_0_2[df_mbg_0_2.index.isin(list(mb_0_2_trans))].drop_duplicates(subset="gene_id")[["chrom","start","end","gene_id","gene_name","strand"]]
bed_mbg_0_2.to_csv(OUTDIR+"m_bias_genes_0_2.bed", sep="\t", index=False, header=False)

bed_fbg_2_4 = df_fbg_2_4[df_fbg_2_4.index.isin(list(fb_2_4_trans))].drop_duplicates(subset="gene_id")[["chrom","start","end","gene_id","gene_name","strand"]]
bed_fbg_2_4.to_csv(OUTDIR+"f_bias_genes_2_4.bed", sep="\t", index=False, header=False)

bed_mbg_2_4 = df_mbg_2_4[df_mbg_2_4.index.isin(list(mb_2_4_trans))].drop_duplicates(subset="gene_id")[["chrom","start","end","gene_id","gene_name","strand"]]
bed_mbg_2_4.to_csv(OUTDIR+"m_bias_genes_2_4.bed", sep="\t", index=False, header=False)

# See the number of genes vs. transcripts when dropping gene ID duplicates

print(len(bed_fbg_0_2))
print(len(bed_fbg_0_2.drop_duplicates(subset="gene_id")))

print(len(bed_mbg_0_2))
print(len(bed_mbg_0_2.drop_duplicates(subset="gene_id")))

print(len(bed_fbg_2_4))
print(len(bed_fbg_2_4.drop_duplicates(subset="gene_id")))

print(len(bed_mbg_2_4))
print(len(bed_mbg_2_4.drop_duplicates(subset="gene_id")))


# Plot to choose TPM threshold

dfs = [df_fbg_0_2, df_mbg_0_2, df_fbg_2_4, df_mbg_2_4]
df_list = ["df_fbg_0_2", "df_mbg_0_2", "df_fbg_2_4", "df_mbg_2_4"]

ctn = 0
for d, dl in zip(dfs, df_list):
    print(dl)
    ctn+=1
    hist_nums = list()
    thresh_list = list()

    for threshold in range(15):
        #print(d)
        df = d.select_dtypes(exclude=['object'])
        num_genes = len(df[df.apply(lambda x : np.mean(x.values) >= threshold, axis=1)])
        print(threshold, num_genes)
        dfT = df.T
        print(dfT[dfT > threshold].count())
        thresh_list.append(threshold)
        hist_nums.append(num_genes)
        
    fig = plt.figure()
    ax = plt.axes()
    plt.title(dl)
    ax.plot(thresh_list, hist_nums,marker=".")
    fig.savefig(OUTDIR+"/testa_plot_%s.png"%(str(ctn)))
    plt.close()