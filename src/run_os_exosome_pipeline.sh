#!/bin/bash

# argv[1] = directory with the input fastq files
# argv[2] = organism
# argv[3] = output dir
# argv[4] = step to process

if [[ $# -eq 0 ]] ; then
    echo 'need to provide additional arguments:
        argv[1] = directory with the input fastq files
        argv[2] = reference genome
	argv[3] = output directory
	argv[4] = step to process'
    exit 0
fi

directory_with_fastq_files=$1 #"input_fastq"
reference_genome=$2
outdir=$3
step_to_process=$4


for i in `find $directory_with_fastq_files -name "*.fastq.gz"`;do 
	echo $i 
	/home/thor/os_exosome_pipeline/os_exosome_pipeline.py -1 $i -r $reference_genome -d $outdir  

done


