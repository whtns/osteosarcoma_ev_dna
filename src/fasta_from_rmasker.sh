#!/bin/bash


rep_fasta=$1
sample=`basename $rep_fasta | sed 's/_.*//'`
echo $sample
out_fasta=`echo $rep_fasta | sed 's/fasta.cat.gz/_cleaned.fasta/'`

zcat $rep_fasta | grep L1P1 | grep -v @ | awk '{ print ">"$2 $4 }' > $out_fasta
