# 6_calc_total_alt_differential_splicing.py
# Ashley Mae Conard
# Last Mod. 5/3/2020
# Purpose: Calculate and plot the proportions of alternative splicing (in pie chart) in treatment samples.

# Import libraries
import pandas as pd
pd.set_option('display.max_colwidth',-1)
import collections as c
import glob, os, sys
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
from io import StringIO
import matplotlib.pyplot as plt; plt.rcdefaults()
get_ipython().run_line_magic('matplotlib', 'inline')
from matplotlib_venn import venn2
import numpy as np

# Import events
clampRNAi_events_0_2_f = "/gpfs/data/compbio/aconard/splicing/results/suppa_results_ncbi_trans/merged_0-2_f/clampRNAi_diffSplice.dpsi"
clampRNAi_events_0_2_m = "/gpfs/data/compbio/aconard/splicing/results/suppa_results_ncbi_trans/merged_0-2_m/clampRNAi_diffSplice.dpsi"
clampRNAi_events_2_4_f = "/gpfs/data/compbio/aconard/splicing/results/suppa_results_ncbi_trans/merged_2-4_f/clampRNAi_diffSplice.dpsi"
clampRNAi_events_2_4_m = "/gpfs/data/compbio/aconard/splicing/results/suppa_results_ncbi_trans/merged_2-4_m/clampRNAi_diffSplice.dpsi"

# Read csv
df_0_2_f1 = pd.read_csv(clampRNAi_events_0_2_f,sep="\t+|;", engine='python')
df_0_2_m1 = pd.read_csv(clampRNAi_events_0_2_m,sep="\t+|;", engine='python')
df_2_4_f1 = pd.read_csv(clampRNAi_events_2_4_f,sep="\t+|;", engine='python')
df_2_4_m1 = pd.read_csv(clampRNAi_events_2_4_m,sep="\t+|;", engine='python')

# 0-2
# females 0-2
df_0_2_f = df_0_2_f1.dropna().loc[(df_0_2_f1.sum(axis=1) != 0)]
df_0_2_f = df_0_2_f.reset_index()
df_0_2_f = df_0_2_f.rename(columns={"level_0": "gene_id", "level_1": "AS"})
df_0_2_f['event_type'] = df_0_2_f['AS'].apply(lambda x: x[0:2])
df_0_2_f_num_AS_events= df_0_2_f.groupby(['gene_id', 'event_type']).size().reset_index(name='num_events')

# males 0-2
df_0_2_m = df_0_2_m1.dropna().loc[(df_0_2_m1.sum(axis=1) != 0)]
df_0_2_m = df_0_2_m.reset_index()
df_0_2_m = df_0_2_m.rename(columns={"index": "AS"})
df_0_2_m = df_0_2_m.rename(columns={"level_0": "gene_id", "level_1": "AS"})
df_0_2_m['event_type'] = df_0_2_m['AS'].apply(lambda x: x[0:2])
df_0_2_m_num_AS_events= df_0_2_m.groupby(['gene_id', 'event_type']).size().reset_index(name='num_events')

# 2-4
# females 2-4
df_2_4_f = df_2_4_f1.dropna().loc[(df_2_4_f1.sum(axis=1) != 0)]
df_2_4_f = df_2_4_f.reset_index()
df_2_4_f = df_2_4_f.rename(columns={"level_0": "gene_id", "level_1": "AS"})
df_2_4_f['event_type'] = df_2_4_f['AS'].apply(lambda x: x[0:2])
df_2_4_f_num_AS_events= df_2_4_f.groupby(['gene_id', 'event_type']).size().reset_index(name='num_events')

# males 2-4
df_2_4_m = df_2_4_m1.dropna().loc[(df_2_4_m1.sum(axis=1) != 0)]
df_2_4_m = df_2_4_m.reset_index()
df_2_4_m = df_2_4_m.rename(columns={"index": "AS"})
df_2_4_m = df_2_4_m.rename(columns={"level_0": "gene_id", "level_1": "AS"})
df_2_4_m['event_type'] = df_2_4_m['AS'].apply(lambda x: x[0:2])
df_2_4_m_num_AS_events= df_2_4_m.groupby(['gene_id', 'event_type']).size().reset_index(name='num_events')

print(df_0_2_f_num_AS_events.head(10))
print(df_0_2_m_num_AS_events.head(10))
print(df_2_4_f_num_AS_events.head(10))
print(df_2_4_m_num_AS_events.head(10))

events_list = ["SE","A5","A3","MX","RI","AF","AL"]
ee_experiments = [("df_0_2_f", df_0_2_f_num_AS_events), ("df_0_2_m", df_0_2_m_num_AS_events), 
                      ("df_2_4_f",df_2_4_f_num_AS_events), ("df_2_4_m", df_2_4_m_num_AS_events)]
for df_name, df_time_sex_events in ee_experiments:
    print("df_time_sex: ", df_name)
    list_evs = []
    list_evs_name = []
    for ev in events_list:
        list_evs_name.append(ev) 
        list_evs.append(df_time_sex_events[df_time_sex_events["event_type"]== ev]["num_events"].sum())
    print(list_evs_name)
    print(list_evs, sum(list_evs), 66927-sum(list_evs))

labels = ['SE', 'A5SS', 'A3SS', 'MXE','RI','AF', 'AL']
sizes = [1641, 1820, 1644, 535, 1014, 2543, 347]
colors = ['red', 'greenyellow', 'lime', 'blue','purple','yellow','orange']
patches, texts = plt.pie(sizes, colors=colors, shadow=True, startangle=90)
#plt.legend(patches, labels, loc="best")
plt.axis('equal')
plt.tight_layout()
plt.show()

labels = ['SE', 'A5SS', 'A3SS', 'MXE','RI']
colors = ['red', 'greenyellow', 'lime', 'blue','purple']
patches, texts = plt.pie(sizes[0:-2], colors=colors, shadow=True, startangle=90)
#plt.legend(patches, labels, loc="best")
plt.axis('equal')
plt.tight_layout()
plt.show()