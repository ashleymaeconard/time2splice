# 5_calc_total_alt_splicing_controls.py
# Ashley Mae Conard
# Last Mod. 6/15/2020
# Purpose: Calculate and plot the proportions of alternative splicing (in pie chart) in control samples (run for each control timepoint and each category (e.g. male and female)).

# Import libraries
import argparse
import pandas as pd
pd.set_option('display.max_colwidth',-1)
import glob, os, sys
import collections as c
from scipy import stats
import seaborn as sns
from io import StringIO
import matplotlib.pyplot as plt; plt.rcdefaults()
get_ipython().run_line_magic('matplotlib', 'inline')
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import numpy as np


def calc_tot_alt_splicing_controls(INPUTFILE, OUTDIR, OUTFILENAME):

    # Read csv
    df_events = pd.read_csv(INPUTFILE,sep="\t+|;", engine='python')

    # Import events
    df_1event_1cat = df_events.dropna().loc[(df_events.sum(axis=1) != 0)]
    df_1event_1cat = df_1event_1cat.reset_index()
    df_1event_1cat = df_1event_1cat.rename(columns={"level_0": "gene_id", "level_1": "AS"})
    df_1event_1cat['event_type'] = df_1event_1cat['AS'].apply(lambda x: x[0:2])
    df_1event_1cat_num_AS_events= df_1event_1cat.groupby(['gene_id', 'event_type']).size().reset_index(name='num_events')

    print(df_1event_1cat_num_AS_events.head(10))

    # Plotting fractions of alternative splicing
    events_list = ["SE","A5","A3","MX","RI","AF","AL"]
    ee_experiments = [("df_1event_1cat", df_1event_1cat_num_AS_events)]
    for df_name, df_time_sex_events in ee_experiments:
        print("df_time_sex: ", df_name)
        list_evs = []
        list_evs_name = []
        for ev in events_list:
            list_evs_name.append(ev) 
            list_evs.append(df_time_sex_events[df_time_sex_events["event_type"]== ev]["num_events"].sum())
        print(list_evs_name)
        print(list_evs, sum(list_evs), 66927-sum(list_evs))

    labels = ['SE', 'A5SS', 'A3SS', 'MXE','RI']
    sizes = [1703, 1918, 1701, 1519, 1060]
    colors = ['red', 'yellowgreen', 'lime', 'blue','purple','yellow','orange']
    patches, texts = plt.pie(sizes, colors=colors, shadow=True, startangle=90)
    plt.legend(patches, labels, loc="best")
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig(OUTDIR+"/"+OUTFILENAME+"_piechart.png")

def main(args):
    INPUTFILE = args.input_file
    OUTDIR = args.output_dir 
    OUTFILENAME = args.output_file_name

    calc_tot_alt_splicing_controls(INPUTFILE, OUTDIR, OUTFILENAME)
    #df_final_2.to_csv(OUTPUTDIR+"iso_tpm_merged.txt",sep="\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate and plot proportions of alternative splicing (in pie chart) in control samples. <!> NOTE: Run this for every timepoint and condition separately.<!>")
    parser.add_argument("input-file", type=str, help="Use events.psi generated in previous step. (e.g. ./results/suppa_results_ncbi_trans/merged_0-2_f/control_iso_events.psi)") 
    parser.add_argument("output-dir", type=str, help="Choose output directory location") 
    parser.add_argument("output-file-name", type=str, help="Choose output file name.") 
    args = parser.parse_args()
    
    main(args)