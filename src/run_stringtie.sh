#!/bin/bash

input_bam=$1
file_pattern=`basename $input_bam | cut -d '_' -f1`
output_gtf=`dirname $input_bam`/$file_pattern/$file_pattern.gtf
echo $file_pattern
echo $output_gtf

stringtie -G ~/Homo_sapiens/grch37/genes.gtf \
	-x MT -eB -p 4 -o $output_gtf \
	$input_bam
