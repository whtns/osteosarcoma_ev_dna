#!/usr/bin/python

import fnmatch
import os
import subprocess
import glob
import sys

stringtie_binary = 'stringtie'
gtf_location = "/home/skevin/storage/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"

output_directory_root = sys.argv[1]
bam_files = glob.glob(output_directory_root+'N2*.bam')

bam_files = []
for root, dirnames, filenames in os.walk(output_directory_root):
    for filename in fnmatch.filter(filenames, '*.bam'):
        bam_files.append(os.path.join(root, filename))
        
num_threads = "4"

for i in bam_files:
	cell_name = os.path.basename(i).replace(".sorted.bam", "") 
	#~ print(cell_name)
	stringtie_output_dir = output_directory_root + "/stringtie/"
	stringtie_cell_dir = stringtie_output_dir + cell_name+"/"
	stringtie_gtf = stringtie_cell_dir + cell_name+".gtf"
	if not os.path.exists(stringtie_output_dir):
		os.makedirs(stringtie_output_dir)
	if not os.path.exists(stringtie_cell_dir):
		os.makedirs(stringtie_cell_dir)
		print "processing stringtie"
		cmd = [stringtie_binary]
		cmd += ["-G", gtf_location] 
		cmd += ["-x", "MT"]
		cmd += ["-eB"]
		cmd += ["-p",num_threads]
		cmd += ["-o",stringtie_gtf]
		cmd += [i]
		out_file = stringtie_cell_dir + "stringtie.log"
		err_file = stringtie_cell_dir + "stringtie.err"
		print " ".join(cmd)
		with open(err_file,"w")as outerr:
			with open(out_file,"w")as outf:
				subprocess.call(cmd,stdout=outf,stderr=outerr)
		del cmd
