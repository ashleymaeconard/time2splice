# 5_calc_total_alt_splicing_controls.py
# Ashley Mae Conard
# Last Mod. 6/15/2020
# Purpose: Calculate and plot the proportions of alternative splicing (in pie chart) in control samples.

# Import libraries
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

# Import events
ee_events_0_2_cf_cm = "/gpfs/data/compbio/aconard/splicing/results/suppa_results_ncbi_trans/merged_0-2_cf_cm/earlyEmbryo_events.psi"
ee_events_2_4_cf_cm = "/gpfs/data/compbio/aconard/splicing/results/suppa_results_ncbi_trans/merged_2-4_cf_cm/earlyEmbryo_events.psi"

# Read csv
df_0_2_cf_cm = pd.read_csv(ee_events_0_2_cf_cm,sep="\t+|;", engine='python')
df_2_4_cf_cm = pd.read_csv(ee_events_2_4_cf_cm,sep="\t+|;", engine='python')

# 0-2
# females 0-2
df_0_2_cf = df_0_2_cf_cm.iloc[:, : 4].dropna().loc[(df_0_2_cf_cm.sum(axis=1) != 0)]
df_0_2_cf = df_0_2_cf.reset_index()
df_0_2_cf = df_0_2_cf.rename(columns={"level_0": "gene_id", "level_1": "AS"})
df_0_2_cf['event_type'] = df_0_2_cf['AS'].apply(lambda x: x[0:2])
df_0_2_cf_num_AS_events= df_0_2_cf.groupby(['gene_id', 'event_type']).size().reset_index(name='num_events')

# males 0-2
df_0_2_cm = df_0_2_cf_cm.iloc[:, 4 : ].dropna().loc[(df_0_2_cf_cm.sum(axis=1) != 0)]
df_0_2_cm = df_0_2_cm.reset_index()
df_0_2_cm = df_0_2_cm.rename(columns={"index": "AS"})
df_0_2_cm = df_0_2_cm.rename(columns={"level_0": "gene_id", "level_1": "AS"})
df_0_2_cm['event_type'] = df_0_2_cm['AS'].apply(lambda x: x[0:2])
df_0_2_cm_num_AS_events= df_0_2_cm.groupby(['gene_id', 'event_type']).size().reset_index(name='num_events')

# 2-4
# females 2-4
df_2_4_cf = df_2_4_cf_cm.iloc[:, : 4].dropna().loc[(df_2_4_cf_cm.sum(axis=1) != 0)]
df_2_4_cf = df_2_4_cf.reset_index()
df_2_4_cf = df_2_4_cf.rename(columns={"level_0": "gene_id", "level_1": "AS"})
df_2_4_cf['event_type'] = df_2_4_cf['AS'].apply(lambda x: x[0:2])
df_2_4_cf_num_AS_events= df_2_4_cf.groupby(['gene_id', 'event_type']).size().reset_index(name='num_events')

# males 2-4
df_2_4_cm = df_2_4_cf_cm.iloc[:, 4 : ].dropna().loc[(df_2_4_cf_cm.sum(axis=1) != 0)]
df_2_4_cm = df_2_4_cm.reset_index()
df_2_4_cm = df_2_4_cm.rename(columns={"index": "AS"})
df_2_4_cm = df_2_4_cm.rename(columns={"level_0": "gene_id", "level_1": "AS"})
df_2_4_cm['event_type'] = df_2_4_cm['AS'].apply(lambda x: x[0:2])
df_2_4_cm_num_AS_events= df_2_4_cm.groupby(['gene_id', 'event_type']).size().reset_index(name='num_events')

print(df_0_2_cf_num_AS_events.head(10))
print(df_0_2_cm_num_AS_events.head(10))
print(df_2_4_cf_num_AS_events.head(10))
print(df_2_4_cm_num_AS_events.head(10))

# Printing fractions of alternative splicing
events_list = ["SE","A5","A3","MX","RI","AF","AL"]
ee_experiments = [("df_0_2_cf", df_0_2_cf_num_AS_events), ("df_0_2_cm", df_0_2_cm_num_AS_events), 
                      ("df_2_4_cf",df_2_4_cf_num_AS_events), ("df_2_4_cm", df_2_4_cm_num_AS_events)]
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
plt.show()