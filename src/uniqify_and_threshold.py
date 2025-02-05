#!/usr/bin/python

import sys
import os 
import subprocess

bam_file=sys.argv[1]
# ~ r1_base_filename=os.path.basename(bam_file)
r1_base_filename = bam_file.split("_")
r1_base_filename="_".join(r1_base_filename[:3])

output_dir = os.path.dirname(bam_file)
print(bam_file)
print(r1_base_filename)
print(output_dir)

exonic_bed_hg19 = os.path.expanduser("~/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes_exonic.bed")
exonic_bed = exonic_bed_hg19

steps_to_process = ["uniqify", "threshold_bam"]


# GET UNIQUE INTRONIC READS
# =====================================================================

# hisat2 -x /home/thor/storage/Homo_sapiens/grch38/genome -1 76/trimmomatic/Hu_76_S127_R1_001.trimmed.fastq -2 76/trimmomatic/Hu_76_S127_R2_001.trimmed.fastq -U 76/trimmomatic/Hu_76_S127_R1_001.trimmed.unpaired.fastq -U 76/trimmomatic/Hu_76_S127_R2_001.trimmed.unpaired.fastq -S Hu_76_hisat2.sam

intronic_bam = r1_base_filename + "_intronic.bam"

if("uniqify" in steps_to_process):
	print "getting intronic reads"
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	cmd = ["bedtools", "intersect", "-wa", "-abam"]
	cmd += [bam_file, "-b", exonic_bed, "-v"]
	print " ".join(cmd)
	err_file = r1_base_filename+"_uniqify.log"
	with open(intronic_bam,"wb")as outf, open(err_file, "wb") as outerr:
		uniq = subprocess.call(cmd, stdout=outf, stderr=outerr)
	print("writing bam to "+intronic_bam)
	del cmd

# THRESHOLD BAM BELOW DEPTH N PILEUPS
# =====================================================================

# hisat2 -x /home/thor/storage/Homo_sapiens/grch38/genome -1 76/trimmomatic/Hu_76_S127_R1_001.trimmed.fastq -2 76/trimmomatic/Hu_76_S127_R2_001.trimmed.fastq -U 76/trimmomatic/Hu_76_S127_R1_001.trimmed.unpaired.fastq -U 76/trimmomatic/Hu_76_S127_R2_001.trimmed.unpaired.fastq -S Hu_76_hisat2.sam

threshold_bam = r1_base_filename + "_depth_under_3.bam"
threshold_bed  = r1_base_filename + "_depth_over_3.bed"

if("threshold_bam" in steps_to_process):
	print "filtering bam at max pileup (default 3)"
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	cmd = ["samtools", "depth", intronic_bam]
	print " ".join(cmd)
	err_file = r1_base_filename+"_threshold_bam.log"
	with open(threshold_bam,"wb")as outbam, open(threshold_bed, "wb") as outbed, open(err_file, "wb") as outerr:
		print " ".join(["samtools", "depth", intronic_bam])
		samtools_depth = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=outerr)
		bed_above_N = subprocess.Popen(["/home/thor/single_cell_pipeline/src/bed_from_areas_covered_above_N_v2.py", "3"], stdin=samtools_depth.stdout, stdout=subprocess.PIPE, stderr=outerr)
		bam_below_N = subprocess.Popen(["bedtools", "intersect", "-wa", "-abam", intronic_bam, "-b", "-", "-v"], stdin=bed_above_N.stdout, stdout=outbam, stderr=outerr)
		samtools_depth.wait()
		bed_above_N.wait()
		bam_below_N.wait()
		bam_below_N.communicate()
	print("writing bed to "+threshold_bed)
	print("writing bam to "+threshold_bam)
	del cmd
