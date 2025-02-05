#!/usr/bin/Rscript

#filepaths
args <- commandArgs(trailingOnly = TRUE)
os_directory =  args[1] #path to osteosarcoma analysis group ex. "/media/thor/storage/os_exosome_pipeline/osteosarcoma_biomarkers_analysis/All_OS_vs_Normal"

#load required libraries
library("stringr") 
library("DESeq") #http://bioconductor.org/packages/release/bioc/html/DESeq.html
require("ggplot2")

#input read counts table and clean for normalization
os_counts = file.path(os_directory, paste(list.files(path = os_directory, pattern = "raw_count_table*"), sep=""))
CountTable = read.table(os_counts, header=TRUE, row.names=1)
CountTable$chrom <- NULL
CountTable$start <- NULL
CountTable$stop <- NULL
CountTable$type <- NULL

#assign read counts columns to control and affected groups for summary analysis
names <- names(CountTable)
cleannames <- word(names, 6, sep="_") #extract group
osDesign = data.frame(row.names = colnames( CountTable ), condition = cleannames)
condition = osDesign$condition
cds = newCountDataSet(CountTable, condition)

#estimate size factor for read count correction (normalization)
cds = estimateSizeFactors(cds)
sizeFactors(cds)
normalized_counts = data.frame(counts(cds, normalized=TRUE))
write.table(normalized_counts, file=file.path(os_directory, paste("normalized_counts.tsv", sep="")), quote=FALSE, sep='\t', col.names = NA)

#estimate disperions
cds = estimateDispersions(cds)

# generrate summary statistics for control and affected groups
res  = nbinomTest(cds, "OS", "normal")
write.table(res, file=file.path(os_directory, paste("res.tsv", sep="")), quote=FALSE, sep='\t', col.names=NA)

##Highlight genes that have an absolute fold change > 2 and a p-value < 0.05
res$threshold = as.factor(abs(res$log2FoldChange) > 2 & res$padj < 0.05)

##Construct the plot object
g = ggplot(data=res, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  list(
    geom_point(alpha=0.4, size=1.75), 
    labs(title = "New plot title"),
    xlim(c(-11, 11)),
    ylim(c(0, 15)),
    xlab("log2 fold change"),
    ylab("-log10 p-value"),
    NULL
  )
g

ggsave(filename = "ks_volcano_plot.pdf", plot = g, path = os_directory)
