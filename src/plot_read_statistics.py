#!/usr/bin/python

# argv[1] = file with original read counts
# Hu_1    5116476
# Hu_2    4962652
# Hu_3    4744324
# Hu_4    1573608

# argv[4] = cell <tab> technology

# argv[5] = output name prefix

bar_figsize = (10,6)
dpi=400
# import system library - takes care of reading arguments from command line
import sys
# pandas library helps with reading csv file and with statistics
import pandas as pd
# library matplotlib takes care of plotting the graphs
import matplotlib.pyplot as plt
# library saborn provides more pleasant color palletes for the graphs
import seaborn as sns


original_counts = pd.read_csv(sys.argv[1],sep="\t",skiprows=1, names = ["cell","original_count"])
cell_technology = pd.read_csv(sys.argv[4],sep="\t",skiprows=1, names = ["cell","technology"])

counts = cell_technology.merge(original_counts,left_on="cell", right_on="cell").merge(trimmed_counts,left_on="cell", right_on="cell").merge(alignment_stats,left_on="cell", right_on="cell")
counts.set_index("cell",inplace=True)

technologies = counts["technology"].unique()

print counts

yrange = max(counts['original_count'])
fig, axes = plt.subplots(2, 3,figsize=bar_figsize);
i=0
for t in technologies:
		# bar plot of 1 technology
	counts[counts["technology"]==t].plot.bar(
		#stacked=True,
		ax=axes[i/3,i%3],
		legend=False, 
		ylim = [0,yrange],
		color=sns.color_palette()
	)
	axes[i/3,i%3].set_title(t)
	i += 1
axes[0,2].legend(loc='upper left', bbox_to_anchor=(1, 1),labels=["sequenced","after trimming","mapped","uniqly mapped","concordant"])
fig.tight_layout()
fig.subplots_adjust(right=0.85)
fig_filename = sys.argv[5]+"_read_count_statistics.png"
fig.savefig(fig_filename, dpi=dpi)
