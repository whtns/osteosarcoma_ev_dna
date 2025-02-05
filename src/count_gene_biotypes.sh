#!/usr/bin/bash

htseq-count -f bam -r pos -i gene_type 10-OS19-HLOH_S43_L007_sorted.bam ~/Homo_sapiens/grch37/gencode.v28lift37.annotation.gtf > test.txt

