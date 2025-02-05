#! /bin/bash

directory_with_bam_files=$1 #"input_fastq"

# file header  
echo -e "bam_file \t total \t mapped \t unmapped"

# loop to count reads in all files and add results to a table  
for i in `find $directory_with_bam_files -name "*.bam"`;do   
	total=$(samtools view -c $i)  
	mapped=$(samtools view -c -F 4 $i)  
	unmapped=$(samtools view -c -f 4 $i)  
	echo -e $i '\t' $total '\t' $mapped '\t' $unmapped  
#	./mapped_unmapped.R $bam_file $total $mapped $unmapped
done  
# ultimately the data is saved in a space-separated file

# saved in a file called "CountMappedReads.sh"
