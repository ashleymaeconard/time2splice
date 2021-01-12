# plot_alignment.py
# Ashley Conard
# Last Mod. June 23, 2019
# Purpose: Plot the alignments from either one or two different aligners (Bowtie2 or HISAT2).

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.lines import Line2D
import random
import seaborn as sns
import glob, sys
import re, os
import pandas as pd
from matplotlib import pyplot
import matplotlib
import matplotlib.pyplot as plt 

if len(sys.argv) != 6 :
    print("Usage: python plot_alignment.py /FULL/PATH/TO/hisat2_htseq/ (with subfolders e.g. ../hisat2_htseq/SAMPLE/REP/htseq_counts stop at hisat2_htseq/) PAIRED_OR_NOT (0-single, 1-paired) METHOD_COMPARE (0-no, 1-yes) METHOD_NAME_1 (e.g. HISAT2) METHOD_NAME_2 (e.g. Bowtie2 or NA)")
    sys.exit (1)    

# Arguments
INPUT_DIR=sys.argv[1]
PAIRED = int(sys.argv[2]) # change to 1 if paired-end sequencing
method_comp = int(sys.argv[3]) # change to 1 if plotting 1 alignment tool vs. another (e.g. Bowtie2 vs. HISAT2)
method_name = str(sys.argv[4])
method2_name = str(sys.argv[5])

def main(argv):
    once = []
    overall = []
    reps = []

    if not method_comp and method2_name=="NA":
        # get all alignment files to read from method (e.g. Bowtie2)
        for alignment_file in glob.glob(INPUT_DIR+"/*/*R1*/*alignment_info*.txt"):
            print("Reading in: ", alignment_file)

            # open alignment file and read all lines
            with open(alignment_file) as f:
                lines = f.readlines()

            # If paired-end, check different lines for conc and overall.
            if PAIRED:
                # get only concordant once line and overall alignment line
                once_line = lines[3:4:1]
                overall_align_line = lines[14:15:1]
            else:
                once_line = lines[3:4:1]
                overall_align_line = lines[5:6:1]

            # get concordant alignment exactly once
            if once_line:
                in_parenthesis=re.search('\(([^)]+)', once_line[0]).group(1).strip("%")
                #print(in_parenthesis)
                once.append(float(in_parenthesis))

            # get overall alignment
            if overall_align_line:
                val=overall_align_line[0].split("%")[0]
                #print(val)
                overall.append(float(val))

                # Get name of replicate and sample files
                name_rep = alignment_file.split("/")[-2]
                reps.append(name_rep)

        align_list = [once, overall]
        df1 = pd.DataFrame(columns=reps, data=align_list)
        if PAIRED:
            df1['d']=["Concordant reads map exactly once","Overall alignment rate"]
        else:
            df1['d']=["Reads map exactly once","Overall alignment rate"]
        df1=df1.reset_index()
        df1.index = df1['d']
        del df1['d']
        df2= df1
    #     # Create dataframe from two input lists of values
    #     if not method_comp:
    #         df = pd.DataFrame(
    #             {"Concordant reads map exactly once": conc, "Overall alignment rate": overall, "label": reps})
    #     elif method_comp:
    #         print("Not set at the moment to compare methods")

    #     df1 = pd.melt(df, id_vars=['label'], var_name='measures', value_name= 'percent')
    #     df2 = df1.set_index('label').T
    #     df2.index.name = "d"
    #     df2['index']=np.arange(len(df2))

        # Plot concordant read mapping and overall alignment score for method and save to file
        ax = plt.figure(figsize=(17.7, 8.27))
        colors = ["blue","grey","lightcoral","green","purple","red", "gold", "olive","peru","aqua",
                  "palegreen","mediumorchid","wheat","darkkhaki","m","saddlebrown","lightblue","teal",
                  "orange","plum","deeppink","crimson","yellowgreen","palevioletred","indigo","mediumaquamarine",
                 "blueviolet","y","indianred","lightskyblue","blue","grey","lightcoral","green","purple","red", "gold", "olive","peru","aqua",
                  "palegreen","mediumorchid","wheat","darkkhaki","m","saddlebrown","lightblue","teal",
                  "orange","plum","deeppink","crimson","yellowgreen","palevioletred","indigo","mediumaquamarine",
                 "blueviolet","y","indianred","lightskyblue"]
        columns = list(df2.columns)
        #ax.set_ylim(0,100)
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print(df2)
        
        for j in range(1, len(df2.columns)):
            color = colors[j-1]
            g = sns.stripplot(x=df2.index, y=columns[j], color=color, jitter=0.3, size=14, dodge=True, data=df2, linewidth=2)

        for j in range(1, len(df2.columns)):
            for i in range(len(df2)):
                g.text(x=i+0.1, y=df2[df2.columns[j]].values[i]+0.001, s=df2.columns[j], horizontalalignment='right', size='medium', color='black')
        elements = [Line2D([0], [0], color=colors[i]) for i in range(len(df2.columns)-1)]

        ax.legend(handles=elements, labels=list(df2.columns)[1:], bbox_to_anchor=(.55,.85), prop={'size': 6})
        #.legend(loc=2,)
        
        plt.ylabel("percent")
        plt.ylim(0,100)
        plt.xlabel("measures")
        #plt.savefig('svm_conf.png', dpi=400)

        #p1= sns.stripplot(x="measures", y="percent", hue="label", size=11, palette="Set1", ax = ax, data=df1, jitter=0.2, linewidth=1)
        if not method_comp:
            title = plt.title("Mapping Metrics for "+method_name)
        else:
            print("Not working with method comparisons yet.")
            #title = plt.title("Alignment Scores Between HISAT2 (w/OUT Splice Site Info. and genome reference) and TopHat")

        plt.savefig(INPUT_DIR+"/"+method_name+"_plot_alignment.pdf")
        return(plt.show())

if __name__ == "__main__":
    main(sys.argv[1:])
    
