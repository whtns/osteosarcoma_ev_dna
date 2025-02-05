#!/bin/bash

# argv[1] = directory with the input fastq files
# argv[2] = organism
# argv[3] = output dir
# argv[4] = step to process

if [[ $# -eq 0 ]] ; then
    echo 'need to provide additional arguments:
        argv[1] = first directory with fastq files
        ...argv[n] = additional directories with fastq files'
    exit 0
fi

fq_dir_1=$1
fq_dir_2=$2
fq_dir_3=$3

regex="[^/]*$"
echo $FP | grep -oP "$regex"
b=$(basename $fq_dir_1)+"/read_table_from_fastqc.csv"

find $fq_dir_1 -name "*_fastqc.zip" | xargs ./get_read_statistics_from_fasqc_zip.py > $b
  

