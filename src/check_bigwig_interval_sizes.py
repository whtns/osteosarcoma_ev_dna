#!/usr/bin/python

# import system library - takes care of reading arguments from command line
import sys
# pandas library helps with reading csv file and with statistics
import pandas as pd
# library matplotlib takes care of plotting the graphs
import matplotlib.pyplot as plt
# library saborn provides more pleasant color palletes for the graphs
import seaborn as sns
import numpy as np

# #bedGraph section chr1:14473-1226777
# chr1    14473   14489   0.459
# chr1    14556   14572   3.672
# chr1    21123   21124   0.918
# chr1    21124   21140   11.934

numfiles = len(sys.argv)-1
dpi=400
line_figsize = (18,40)
fig, ax = plt.subplots(numfiles, 2, figsize=line_figsize)

for file_no, filename in enumerate(sys.argv[1:]):
	filename_root= filename.split("/")[-1]
	w = pd.read_csv(filename,comment="#", delim_whitespace=True, header=None, names=["chr","from","to","value"])
	
	weighted_sizes = []
	it = w.itertuples()
	last_row = it.next()
	start_segment = last_row[2]
	weight_of_segment = 0
	
	
	
	for this_row in it:
		if((last_row[3] < this_row[2]) or (last_row[1] != this_row[1])): #if previous interval ended before this one starts
			segment_size = last_row[3] - start_segment
			weighted_sizes.append( (segment_size,weight_of_segment) )
			start_segment = this_row[2]
			weight_of_segment = 0
			
		weight_of_segment += this_row[4]
		last_row = this_row
		#else:
		#	print "Y", last_row[3], this_row[2]
	
	#print sizes
	ws = pd.DataFrame(weighted_sizes, columns= ["size","w"])
	
	
	ws["size"].hist(ax=ax[file_no,0], bins=100)
	ax[file_no,0].set_yscale("log")
	ax[file_no,0].set_title(filename+" (sizes)")
	ws["size"].hist(ax=ax[file_no,1], bins=100, weights=ws["w"])
	ax[file_no,1].set_yscale("log")
	ax[file_no,1].set_title(filename+" (weighted sizes)")
	
fig.savefig("test.sizes.png",dpi=dpi)
