#!/bin/bash

rep_fasta=$1
sample=`basename $rep_fasta | sed 's/_.*//'`
echo $sample
sample_fastq_gz=`find ~/os_exosome_pipeline/data/FASTQ/ -name "$sample*fastq.gz"`
echo $sample_fastq_gz
sample_fastq=`echo $sample_fastq_gz | sed 's/.gz//'`
gunzip -k $sample_fastq_gz
repeatmask_fastq=`echo $sample_fastq_gz | sed 's/.fastq.gz/_rmask.fastq/'`
out_reads=`basename $rep_fasta | sed 's/\..*//'`.reads

zcat $rep_fasta | grep L1P1 | grep @ | cut -d ' ' -f5 | cut -c 2- > $out_reads

filterbyname.sh in=$sample_fastq out=$repeatmask_fastq names=$out_reads include=t




