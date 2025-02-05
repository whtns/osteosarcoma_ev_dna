#!/usr/bin/python

import subprocess
import sys
import re
import os
import datetime

gtf_location = "/media/thor/storage/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
transcriptome_index_location = gtf_location.rsplit(".",1)[0]
bowtie2index_location = "/media/thor/storage/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"

tophat_output_dir = "./tophat"

num_threads = "4"

tophat_binary = "tophat2"

fastq = sys.argv[1]
basename = fastq.split("/")[-1].split(".")[0]


tophat_arguments = ["--read-mismatches", "2", "--read-gap-length", "2", "--read-edit-dist", "2", "--max-multihits", "10", "--library-type", "fr-unstranded"]
print bowtie2index_location
cmd = [tophat_binary]
cmd += tophat_arguments
cmd += ["--GTF",gtf_location]
cmd += ["--transcriptome-index",transcriptome_index_location]
cmd += ["--num-threads", num_threads]
cmd += ["--output-dir",tophat_output_dir]
cmd += [bowtie2index_location]
cmd += [fastq]

print " ".join(cmd)
out_file = "tophat/"+basename+".log"
with open(out_file,"w")as outf:
 	subprocess.call(cmd,stdout=outf)
 	print cmd
del cmd
