#!/usr/bin/python

# argv[1,2] = R1/R2.fastq files
# argv[3] = reference genome (*.fasta)
# argv[4] = step to process

import gzip
import subprocess
import sys
import re
import os
import datetime
import argparse

# list here all steps of the pipeline. Pipeline will run all these steps, if not requested otherwise in -s argument
steps_to_process_all = ["trimmomatic", "prinseq", "bowtie", "samformatconverter", "sortsam", "addrg", "buildbamindex"]

parser = argparse.ArgumentParser(description="runs single cell pipeline")
parser.add_argument("-1", "--fastq-r1", dest="f1", help="fastq R1 file REQUIRED", metavar="FILE[.gz]", required=True)
parser.add_argument("-r", "--reference", dest="ref", default= "/media/thor/storage/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa", help="reference genome")
parser.add_argument("-d", "--out", dest="out_dir", help="output directory REQUIRED", metavar="DIRECTORY", required=True)
parser.add_argument("-s", "--steps-to-run", dest="steps", help="steps of pipeline to run")
parser.add_argument("--overwrite", dest="overwrite")

options = parser.parse_args()
fastq_location = options.f1
reference_genome = options.ref
output_directory = options.out_dir
if(options.steps == None):
	steps_to_process = steps_to_process_all
else:
	steps_to_process = options.steps.split(",")

print "Will run steps:", steps_to_process

num_threads = "7"

# *********************************************************************
# DEFINITION OF PATHS

fastqc = ["fastqc"]
Clumpify = ["clumpify"]
Trimmomatic = ["java", "-jar", "/usr/share/java/trimmomatic.jar"]
illuminaclip = "/media/thor/storage/os_exosome_pipeline/contaminants.fa"
min_base_quality = "20"
sliding_window = "10:25"
min_read_length = "25"
Prinseq = ["prinseq-lite"]
BBMerge = ["bbmerge"]
BowTie = ["/home/thor/TOOLS/bowtie-1.2/bowtie"]
SamTools = ["samtools"]
SamFormatConverter = ["java",  "-jar", "/usr/share/java/picard.jar", "SamFormatConverter"]
SortSam  = ["java", "-jar", "/usr/share/java/picard.jar", "SortSam"]
BuildBamIndex  = ["java", "-jar", "/usr/share/java/picard.jar", "BuildBamIndex"]
RNASeQC = ["/usr/lib/jvm/java-7-oracle/bin/java", "-jar", "/usr/share/java/RNA-SeQC_v1.1.8.jar"]
AddOrReplaceReadGroups = ["java", "-jar", "/usr/share/java/picard.jar", "AddOrReplaceReadGroups"]
MarkDuplicates = ["java", "-jar", "/usr/share/java/picard.jar", "MarkDuplicates"]
#genome_file = sys.argv[4]
bowtie_index = reference_genome.replace("WholeGenomeFasta/genome.fa", "BowtieIndex/genome")

#bam_base_filename = bam_location.replace(".sorted.bam","").split("/")[-1]
fastq_base_filename = fastq_location.replace("_R1_001.fastq.gz", "").split("/")[-1]


# /DEFINITION OF PATHS
# *********************************************************************

# READ GROUP PARSING
# =====================================================================
#if("rg_parse" in steps_to_process) or ("all" in steps_to_process):
with gzip.open(fastq_location) as fasta:
	s = fasta.readline().split(":")
	RGID =  "_".join(s[2:5])
	RGSM = fastq_base_filename
	RGLB = "ZEN456A1LI5"
	RGPU = "400"
	RGPL = "Illumina"

# fastqc fastq quality control
# =====================================================================
fastqc_out = "/media/thor/storage/os_exosome_pipeline/FASTQC/"+fastq_base_filename+".fastqc.html"
if(not os.path.isfile(fastqc_out)) and (("fastqc" in steps_to_process)):
	print "running fastqc on the original file"
	cmd = ["fastqc","--nogroup",fastq_location]
	print " ".join(cmd)
	subprocess.call(cmd)
	del cmd

# CLUMPIFY REMOVE DUPLICATES
# =====================================================================
clumpify_output_dir = output_directory+"clumpify/"
fastq_clumped = clumpify_output_dir+fastq_base_filename+"_clumped.c1.fastq.gz"
max_subs = "2"
if(not os.path.isfile(fastq_clumped)) and(("clumpify" in steps_to_process) or ("all" in steps_to_process)):
	subprocess.call(["mkdir",clumpify_output_dir])
	print "running clumpify on the original file"
	cmd = ["clumpify"]
	cmd += ["-Xmx10g"]
	cmd += ["in="+fastq_location]
	cmd += ["out="+fastq_clumped]
	cmd += ["dedupe", "subs="+max_subs]
	out_file = clumpify_output_dir + fastq_base_filename + "clumpify.c1.log"
	err_file = clumpify_output_dir + fastq_base_filename + "clumpify.c1.err"
	print " ".join(cmd)
	with open(err_file,"w")as outerr, open(out_file,"w")as outf:
		errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
		if(errcode == 0):
			print "clumpify finished successfully"
		else:
			print "clumpify failed !!!!"
			del steps_to_process[:]
	del cmd
        
# fastq-mcf trim read due to adapters and some quality measures
# =====================================================================
fastqmcf_output_dir = output_directory + "fastqmcf/"
fastq_trimmed = fastqmcf_output_dir + fastq_base_filename + ".trimmed.fastq"
if(not os.path.isfile(fastq_trimmed)) and(("fastqmcf" in steps_to_process) or ("all" in steps_to_process)):
	subprocess.call(["mkdir",fastqmcf_output_dir])
	print "running fastqc-mcf on the original file"
	cmd = ["fastq-mcf"]
	cmd += ["-s", "0.0"]
	cmd += ["-x", "0"]
	cmd += ["-p", "4"]
	cmd += ["-k", "0"]
	cmd += ["-q", "0"]
	cmd += ["-l", "48"]
	cmd += ["-o", fastq_trimmed]
	cmd += [illuminaclip]
	cmd += [fastq_clumped]
	print " ".join(cmd)
	subprocess.call(cmd)
	del cmd  

# TRIMMOMATIC
# =====================================================================
trimmomatic_output_dir = output_directory + "trimmomatic/"

fastq_trimmomatic = trimmomatic_output_dir+ fastq_base_filename + ".quality.trimmed.fastq"

if(not os.path.isfile(fastq_trimmomatic)) and(("trimmomatic" in steps_to_process) or ("all" in steps_to_process)):
	subprocess.call(["mkdir",trimmomatic_output_dir])
	print "running trimmomatic"
	cmd = Trimmomatic[:]
	#print cmd
	cmd.append("SE")
	cmd.append("-phred33")
	cmd.append("-threads")
	cmd.append("4")
	#cmd.append("-trimlog")
	#cmd.append(log)
	cmd.append(fastq_trimmed)
	cmd.append(fastq_trimmomatic)
	cmd.append("LEADING:"+min_base_quality)
	cmd.append("TRAILING:"+min_base_quality)
	cmd.append("SLIDINGWINDOW:"+sliding_window)
	cmd.append("MINLEN:"+min_read_length)
	print " ".join(cmd)
	#trimmomatic = subprocess.Popen(cmd)
	#trimmomatic.wait()
	#print trimmomatic.returncode
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "trimmomatic for sample "+fastq_base_filename+" finished successfully"
	else:
		print "trimmomatic for sample"+fastq_base_filename+" failed !!!!"
		del steps_to_process[:]
	del cmd

# PRINSEQ-LITE
# =====================================================================
#Command: prinseq-lite.pl -verbose -fastq <Input_trimmed.fastq> -trim_tail_left 5 -trim_tail_right 5 -trim_ns_right 2 -trim_qual_left 28 -trim_qual_right 28 -min_qual_mean 28 -min_len 16 -out_good <Input_trimmed_score> -out_bad <Input_trimmed_removed>

prinseq_output_dir = output_directory + "prinseq/"

prinseq_good = prinseq_output_dir+ fastq_base_filename + "_prinseq_good.trimmed"
prinseq_good_fastq = prinseq_good+".fastq"
prinseq_bad = prinseq_output_dir+ fastq_base_filename + "_prinseq_removed"
prinseq_bad_fastq = prinseq_bad+".fastq"
if(not os.path.isfile(prinseq_good_fastq)) and(("prinseq" in steps_to_process) or ("all" in steps_to_process)):
	subprocess.call(["mkdir",prinseq_output_dir])
	print "running prinseq"
	cmd = Prinseq[:]
	#print cmd
	cmd.append("-verbose")
	cmd.append("-min_len")
	cmd.append("16")
	cmd.append("-trim_tail_left")
	cmd.append("5")
	cmd.append("-trim_tail_right")
	cmd.append("5")
	cmd.append("-lc_method")
	cmd.append("dust")
	cmd.append("-lc_threshold")
	cmd.append("7")
	cmd.append("-fastq")
	cmd.append(fastq_trimmomatic)
	cmd.append("-out_good")
	cmd.append(prinseq_good)
	cmd.append("-out_bad")
	cmd.append(prinseq_bad)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "prinseq for sample "+fastq_base_filename+" finished successfully"
	else:
		print "prinseq for sample"+fastq_base_filename+" failed !!!!"
		del steps_to_process[:]
	del cmd


# BOWTIE
# =====================================================================
#Bowtie
bowtie_output_dir = output_directory+"bowtie/"
bowtie_unmapped_fastq = bowtie_output_dir+fastq_base_filename+"_unmapped.fastq.gz"
bowtie_sam = bowtie_output_dir+fastq_base_filename + ".sam"
if not os.path.isdir(bowtie_output_dir):
	subprocess.call(["mkdir", bowtie_output_dir])
if(not os.path.isfile(bowtie_sam)) and(("bowtie" in steps_to_process) or ("all" in steps_to_process)):
	cmd = BowTie[:]
	cmd.append("-a")
	cmd.append("--best")
	cmd.append("--strata")
	cmd.append("-v")
	cmd.append("2")
	cmd.append("--chunkmbs")
	cmd.append("200")
	cmd.append("-p")
	cmd.append("4")
	cmd.append("-q")
	cmd.append("-S")
	cmd.append("--un")
	cmd.append(bowtie_unmapped_fastq)
	cmd.append(bowtie_index)
	cmd.append(prinseq_good_fastq)
	cmd.append(bowtie_sam)
	out_file = bowtie_output_dir + "bowtie.log"
	err_file = bowtie_output_dir + "bowtie.err"
	print " ".join(cmd)
	with open(err_file,"w")as outerr, open(out_file,"w")as outf:
		errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
		if(errcode == 0):
			print "bowtie finished successfully"
		else:
			print "bowtie failed !!!!"
			del steps_to_process[:]
	del cmd
	
	
# SAMFORMATCONVERTER
# =====================================================================
#SamFormatConverter I=bwa_output.sam O=converted_bwa_output.bamt
picard_output_dir = output_directory+"picard/"
picard_bam = picard_output_dir+fastq_base_filename + ".bam"
if not os.path.isdir(picard_output_dir):
	subprocess.call(["mkdir", picard_output_dir])
if(not os.path.isfile(picard_bam)) and(("samformatconverter" in steps_to_process) or ("all" in steps_to_process)):
	cmd = SamFormatConverter[:]
	cmd.append("I=")
	cmd.append(bowtie_sam)
	cmd.append("O=")
	cmd.append(picard_bam)
	cmd.append("MAX_RECORDS_IN_RAM=200000")
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "samformatconverter finished successfully"
	else:
		print "samformatconverter failed !!!!"
		del steps_to_process[:]
	del cmd

# SORTSAM
# =====================================================================
#SortSam I=input.bam O=sorted.bam SORT_ORDER=coordinate 
sortsam_bam = picard_output_dir+fastq_base_filename + "_sorted.bam"
if (not os.path.isfile(sortsam_bam)) and(("sortsam" in steps_to_process) or ("all" in steps_to_process)):
	print "sortsam_bam:",sortsam_bam
	cmd = SortSam[:]
	cmd.append("I=")
	cmd.append(picard_bam)
	cmd.append("O=")
	cmd.append(sortsam_bam)
	cmd.append("SORT_ORDER=coordinate")
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "sortsam finished successfully"
	else:
		print "sortsam failed !!!!"
		del steps_to_process[:]
	del cmd
	


# Picard AddorReplace Read Groups
# =====================================================================
#SamFormatConverter I=bwa_output.sam O=converted_bwa_output.bamt
rg_bam = picard_output_dir + fastq_base_filename+ "_with_rg.bam"
if not os.path.isdir(picard_output_dir):
	subprocess.call(["mkdir", picard_output_dir])
if(not os.path.isfile(rg_bam)) and(("addrg" in steps_to_process) or ("all" in steps_to_process)):
	cmd = AddOrReplaceReadGroups[:]
	cmd.append("I=")
	cmd.append(sortsam_bam)
	cmd.append("O=")
	cmd.append(rg_bam)
	cmd.append("RGID=")
	cmd.append(RGID)
	cmd.append("RGLB=")
	cmd.append(RGLB)
	cmd.append("RGPL=")
	cmd.append(RGPL)
	cmd.append("RGPU=")
	cmd.append(RGPU)
	cmd.append("RGSM=")
	cmd.append(RGSM)
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "add or replace read groups finished successfully"
	else:
		print "add or replace read groups failed !!!!"
		del steps_to_process[:]
	del cmd

# BUILDBAMINDEX
# =====================================================================
#BuildBamIndex I=input.bam 
bam_index= picard_output_dir+fastq_base_filename+"_with_rg.bam.bai"
if (not os.path.isfile(bam_index)) and(("buildbamindex" in steps_to_process) or ("all" in steps_to_process)):
	cmd = BuildBamIndex[:]
	cmd.append("I=")
	cmd.append(rg_bam)
	cmd.append("O=")
	cmd.append(bam_index)
	print " ".join(cmd)
	out_file = picard_output_dir+ fastq_base_filename+"_buildbamindex.log"
	err_file = picard_output_dir+ fastq_base_filename+"_buildbamindex.err"
	with open(err_file,"w")as outerr:
		with open(out_file,"w")as outf:
			errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
			if(errcode == 0):
				print "buildbamindex finished successfully"
			else:
				print "buildbamindex failed !!!!"
				del steps_to_process[:]
	del cmd

# RNA-SeQC
# =====================================================================
rnaseqc_output = picard_output_dir + fastq_base_filename+ "_rnaseqc/"
gtf_location = "/home/thor/storage/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
if("rnaseqc" in steps_to_process) or ("all" in steps_to_process):
	cmd = RNASeQC[:]
	cmd.append("-r")
	cmd.append(reference_genome)
	cmd.append("-o")
	cmd.append(rnaseqc_output)
	cmd.append("-s")
	cmd.append(fastq_base_filename+"|"+rg_bam+"|"+"notes")
	cmd.append("-t")
	cmd.append(gtf_location)
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "rnaseqc finished successfully"
	else:
		print "rnaseqc failed !!!!"
		del steps_to_process[:]
	del cmd

	

	


