# suppa_formatting.py 
# Ashley Mae Conard
# Last Mod. 9/12/2019
# Purpose: Converts NM_ gene names to flybase name, then merging outputs from run_suppa (NM_ gene names by 1 TPM value column for each replicate)

# Libraries
import argparse
import pandas as pd
pd.set_option('display.max_colwidth',-1)
import mygene
import glob
import collections as c
import os
import sys
from io import StringIO
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
from biothings_client import get_client

def plot_NaN(filename, df3):
    # Total number of NaNs when converting gene IDs
    list_nans=[len(df3), df3["ensembl.gene"].isna().sum(), df3["reporter.DroGene-1_0"].isna().sum(), 
          df3["symbol"].isna().sum(),df3["uniprot.Swiss-Prot"].isna().sum(),
          df3["uniprot.TrEMBL"].isna().sum()]

    objects = ('total.genes','flybase.ID', 'reporter.DroGene','symbol', 'uniprot.Swiss-Prot', 
               'uniprot.TrEMBL')
    y_pos = np.arange(len(objects))
    conversion = list_nans

    fig = plt.figure()
    plt.bar(y_pos, conversion, align='center', alpha=0.5)
    plt.xticks(y_pos, objects, rotation=35)
    plt.grid(True)
    plt.ylabel('Usage')
    plt.title('How many NaNs when converting NM_ ID to other IDs?')
    fig.savefig(OUTPUTDIR+"/"+"gene_id_conversion_"+filename+".pdf")
    
def convert_nm_ids_to_flybase(df1):
    # Remove all rows summing to 0
    df11 = df1.loc[~(df1==0).all(axis=1)]
    
    # Use mygene to change refseq into flybase ID
    mg = mygene.MyGeneInfo()
    mg = get_client('gene')

    # Calling mygene to map NM_ IDs to Flybase names 
    print("Calling mygene.")
    refseq_list = df11.index.tolist()
    df_geneIDs = mg.querymany(refseq_list, scopes="refseq", 
    fields=["ensembl.gene","uniprot","symbol", "reporter"], 
    species="fruitfly", as_dataframe=True)
    new_index_list = df_geneIDs["ensembl.gene"].tolist()

    # Plotting loss of gene IDs per replicate
    plot_NaN("df_merged", df_geneIDs)

    df11['flybase_id']= new_index_list
    cols = list(df11.columns)
    cols = [cols[-1]]+cols[:-1]
    df11 = df11[cols]
    # Adding the flybase names to dataframe
    #df1 = df1.set_index([pd.Index(new_index_list)])
    #df_converted=df1.reset_index().dropna().set_index("gene_id")
    #print("Convertion complete: \n", df_converted.head())
    #print("df_merged lost ", len(df1)-len(df_converted)," thus ",
    #      1-(len(df_converted)/len(df1)),"% gene IDs.")
    return df11

def read_nm_fb_map(DM6_NM_FB):
    df = pd.read_csv(DM6_NM_FB, usecols=[0,3], names=["flybase_id","refseq"],sep="\t")
    dict_n_f = pd.Series(df.flybase_id.values,index=df.refseq).to_dict()
    print(dict(list(dict_n_f.items())[0:2]))
    return dict_n_f

def remapping_NM_to_FBtr(df_merged_iso_tpms, DM6_NM_FB_MAP, OUTPUTDIR):
    
    # Adding flybase ID to df
    dict_nm_fb = read_nm_fb_map(DM6_NM_FB_MAP)

    # Remove .num from NM_ name
    list_nm = df_merged_iso_tpms.index.tolist()
    list_nm_no_nums = [sub[ : -2] for sub in list_nm]
    df_merged_iso_tpms.index = list_nm_no_nums

    df_merged_iso_tpms['flybase_id']= df_merged_iso_tpms.index.map(dict_nm_fb)
    #print(df_merged_iso_tpms)
    
    # Check if any gene_ids were unable to be converted
    if df_merged_iso_tpms[df_merged_iso_tpms.isna().any(axis=1)].empty:
        print("ALL REFSEQ IDs were converted to Flybase transcript IDs!")
    else:
        print("WARNING: Lost transcript IDs were not converted from REFSEQ to Flybase tr. IDs:")
        print(df_merged_iso_tpms[df_merged_iso_tpms.isna().any(axis=1)].shape[0]/
             df_merged_iso_tpms.shape[0], " that is ", 
             df_merged_iso_tpms[df_merged_iso_tpms.isna().any(axis=1)].shape[0], "/",
             df_merged_iso_tpms.shape[0])
        df_nan = df_merged_iso_tpms[df_merged_iso_tpms.isna().any(axis=1)]
        df_nan.to_csv(OUTPUTDIR+"nan_IDs_iso_tpm_merged.txt",sep="\t")
    
   # Move FBtr id to front and remove NaN
    df_merged_iso_tpms = df_merged_iso_tpms.dropna()
    df_merged_iso_tpms = df_merged_iso_tpms.set_index(list(df_merged_iso_tpms)[-1])
    df_merged_iso_tpms = df_merged_iso_tpms.rename_axis(None)
    
    return df_merged_iso_tpms

def reformat_merge_iso_tpm(INPUTDIR, TIMEPOINT, SEX, CONTROLS_ONLY):
    dict_samples_reps={}  
    converted_first_replicate=0
    
    # Iterating through samples
    if CONTROLS_ONLY:
        sample_types_to_compare = INPUTDIR+"/"+TIMEPOINT+"*"+"control"
    else:
        sample_types_to_compare = INPUTDIR+"/"+TIMEPOINT+"_"+SEX+"_*"
    for sample_dir in glob.glob(sample_types_to_compare):
        print("sample_dir", sample_dir)
        if os.path.isdir(sample_dir): # check only folders
            sample_name = sample_dir.split("/")[-1]

            # Iterating through replicates
            rep_count=0
            for rep_dir in glob.glob(sample_dir+"/*"):
                print("rep_dir", rep_dir)
                rep_count+=1
                rep_name = rep_dir.split("/")[-1]
                sample_rep_name = sample_name+"."+str(rep_count)
                dict_samples_reps[rep_name]= sample_rep_name

                # Iterating through iso_tmp.txt files
                for file in glob.glob(rep_dir+"/iso_tmp.txt"):
                    print("file", file)
                    df = pd.read_csv(file, sep="\t", index_col=0)
                    print("Converting this dataframe to have sample.repnum and Flybase IDs.\n")
                    #print(df.head())

                    # Changing replicate_name column to sample_name.repnum column
                    print("Changing ", list(df)[0], " to ", dict_samples_reps[list(df)[0]])
                    df_sample_rep = df.rename(columns={list(df)[0]:dict_samples_reps[list(df)[0]]})
                    df_sample_rep.index.names=['gene_id']

                    # Merging df_converted into new dataframe
                    if not converted_first_replicate:
                        print("Creating merged df")
                        converted_first_replicate=1
                        df_merged_iso_tpms1 = df_sample_rep

                    else:
                        print("Adding", sample_rep_name, "df to merged_df.")
                        df_merged_iso_tpms1 = pd.merge(df_merged_iso_tpms1, df_sample_rep, on='gene_id')
                        #print("Merged dfs: \n", df_merged_iso_tpms.head())
    return df_merged_iso_tpms1

def adding_flybase_IDs(df_merged_iso_tpms, DM6_NM_FB_MAP):
    # Adding flybase ID to df
    dict_nm_fb = read_nm_fb_map(DM6_NM_FB_MAP)
    df_merged = convert_nm_ids_to_flybase(df_merged_iso_tpms)
    df_merged2 = df_merged.fillna(value="hi")

    # Remove .num from NM_ name
    list_nm = df_merged2.index.tolist()
    list_nm_no_nums = [sub[ : -2] for sub in list_nm]
    df_merged2.index = list_nm_no_nums

    # Splitting df_merged into those with converted flybase names and those NaN
    df_merged2_NAN = df_merged2[df_merged2.flybase_id=='hi'].copy()
    df_merged2_noNAN = df_merged2[~df_merged2.flybase_id.str.contains('hi')]

    # Map NaNs using file from Gil dos Santos at Flybase
    df_merged2_NAN['flybase_id']=df_merged2_NAN.index.map(dict_nm_fb)

    # Merge both dataframes back in!
    df_final = df_merged2_noNAN.append(df_merged2_NAN)

    # Have multiple indices because flybase_id are not unique
    # Have FBgn and FBtr
    df_final = df_final.set_index(list(df_final)[0:4])
    df_final = df_final.rename(columns={"flybase_id": ""})

    # Check if any gene_ids were unable to be converted
    if df_final[df_final.isna().any(axis=1)].empty:
        print("ALL REFSEQ IDs were converted to Flybase IDs!")
    else:
        print("WARNING: These REFSEQ IDs were not converted to Flybase IDs:")
        print(df_final)
    return df_final

def plot_NaN_fraction(df2):
    # Fraction of NaNs when converting gene IDs
    frac_nas = [x/len(df2) for x in list_nas]
    objects = ('flybase.ID', 'reporter.DroGene','symbol', 'uniprot.Swiss-Prot', 'uniprot.TrEMBL')
    y_pos = np.arange(len(objects))
    conversion = frac_nas[1:6]

    plt.bar(y_pos, conversion, align='center', alpha=0.5)
    plt.xticks(y_pos, objects, rotation=45)
    plt.grid(True)
    plt.ylabel('Usage')
    plt.title('What fraction of total genes are NaNs when converted from NM_ ID?')

    plt.show()

def main(args):
    """Formatting SUPPA arguments."""
    INPUTDIR = args.inputdir 
    TIMEPOINT = args.time_point
    SEX = args.sex
    OUTPUTDIR = args.outputdir
    CONTROLS_ONLY = args.controls_only
    DM6_NM_FB_MAP = args.map

    # Create output dir if needed
    if not os.path.exists(OUTPUTDIR):
        os.makedirs(OUTPUTDIR)

    df_final_1 = reformat_merge_iso_tpm(INPUTDIR, TIMEPOINT, SEX, CONTROLS_ONLY)
    df_final_2 = remapping_NM_to_FBtr(df_final_1, DM6_NM_FB_MAP, OUTPUTDIR)
    #df_final_2 = adding_flybase_IDs(df_final_1, DM6_NM_FB_MAP)
    df_final_2.to_csv(OUTPUTDIR+"iso_tpm_merged.txt",sep="\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Format SUPPA")
    parser.add_argument("-t", "--time-point", nargs='?', default="0", type=str) # experiment timepoint (could be 0, signifying none, or e.g. "2-4")
    parser.add_argument("--controls-only", default=0, type=int, help="1 if controls only, 0 otherwise. Default is 0.")

    parser.add_argument("inputdir", type=str) # e.g. "/data/compbio/aconard/splicing/results/suppa_results_ncbi_trans/"
    parser.add_argument("outputdir", type=str) #e.g. "/data/compbio/aconard/splicing/results/suppa_results_ncbi_trans/merged_2-4_cf_cm/"
    parser.add_argument("sex", type=str, help="m or f")
    parser.add_argument("map", type=str, help="map from DM6, NM and FB (e.g. see BDGP6's fbtr_refseq.tsv file, could be less or more recent than dm6)")
    
    args = parser.parse_args()
    
    main(args)
