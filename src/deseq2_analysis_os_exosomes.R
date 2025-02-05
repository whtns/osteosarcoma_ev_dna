#!/usr/bin/Rscript

library(Rsubread)
library(DESeq2)
library(tidyverse)
library(gtools)
library(biomaRt)
library(data.table)
library(regionReport)

gene_sym_from_trsid <- function(x,grp_comp=grp_comp){
  ensembl  <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl")
  biomaRt_result <- getBM(attributes=c('ensembl_transcript_id', 'ensembl_gene_id', 'hgnc_symbol'), filters =
                            'ensembl_transcript_id', values = rownames(x), mart = ensembl)
  biomaRt_result <- data.frame(biomaRt_result, stringsAsFactors = FALSE)
  output <- setDT(as.data.frame(x), keep.rownames = TRUE)[]
  output <- rename(output, ensembl_transcript_id = rn)
  output <- full_join(as.data.frame(biomaRt_result), output, by = "ensembl_transcript_id") 
  filename = paste0(grp_comp,"_deseq2",".csv")
  write.table(output,filename, sep="\t", quote=FALSE, row.names = FALSE)
}


GTF=path.expand("~/os_exosome_pipeline/bin/gencode.v26lift37.annotation.gtf")
EXPTNAME="exosome_rna"
CPUS=6
MAPQ=10

#Make the total gene-wise matrix
fls <- dir("~/os_exosome_pipeline/output/picard/", "_sorted.bam$", full.names = TRUE)
coldata_all <- dir("~/os_exosome_pipeline/output/picard/", "_sorted.bam$")

condition = c(rep("OS", 3), rep("ctrl", 7), rep("OS", 1), rep("ctrl", 5), rep("OS", 8))
hospital = c(rep("HLOH", 3), rep("CHLA", 1), rep("HLOH", 1), rep("CHLA", 4), rep("HLOH", 1), rep("CHLA", 1), rep("HLOH", 5), rep("CHLA", 5), rep("HLOH", 3))


coldata_all <- data.frame(coldata_all, condition, hospital, row.names = 1)
fc_all <- featureCounts(fls, annot.ext = GTF, 
                        isGTFAnnotationFile = TRUE,
                        GTF.featureType="transcript",
                        ignoreDup=TRUE,
                        GTF.attrType="transcript_id") 
counts <- fc_all$counts
anno <- data.frame(rownames(counts))


#DESEQ2 for all samples begins
#============================

all_matrix <- fc_all$counts
colnames(all_matrix) <- rownames(coldata_all)
#coldata <- coldata[colnames(counts),]
all(rownames(coldata_all) == colnames(all_matrix))

dds_all <- DESeqDataSetFromMatrix(countData = all_matrix,
                              colData = coldata_all,
                              design = ~ condition)

dds_all <- dds_all[ rowSums(DESeq2::counts(dds_all)) > 1, ]
dds_all <- DESeq(dds_all, parallel = TRUE)


#dds_all$group <- factor(paste0(dds_all$condition, ddsMF$hospital))
#design(ddsMF) <- ~ group
# dds_all <- DESeq(dds_all, parallel = TRUE)
resultsNames(dds_all)

#DESEQ2 for all samples age less than 20 begins
#============================
remove <- c("N25", "N21", "N28", "OS27", "OS11")
raw_matrix <- fc_all$counts
all_cleaned_matrix <- raw_matrix[, grep(paste(remove,collapse="|"), colnames(raw_matrix), invert = TRUE, value = TRUE)]
coldata_all_cleaned <- dir("/media/thor/storage/os_exosome_pipeline/dedupe_output/picard/", "_sorted.bam$")
coldata_all_cleaned <- grep(paste(remove,collapse="|"), coldata_all_cleaned, invert = TRUE, value = TRUE)
colnames(all_cleaned_matrix) <- rownames(coldata_all_cleaned)
#coldata <- coldata[colnames(counts),]
all(rownames(coldata_all_cleaned) == colnames(all_cleaned_matrix))
condition_all_cleaned = c(rep("OS", 2), rep("ctrl", 7), rep("OS", 1), rep("ctrl", 2), rep("OS", 7))
coldata_all_cleaned <- data.frame(coldata_all_cleaned, condition_all_cleaned, row.names = 1)

dds_all_cleaned <- DESeqDataSetFromMatrix(countData = all_cleaned_matrix,
                                  colData = coldata_all_cleaned,
                                  design = ~ condition_all_cleaned)

dds_all_cleaned <- dds_all_cleaned[ rowSums(DESeq2::counts(dds_all_cleaned)) > 1, ]
dds_all_cleaned <- DESeq(dds_all_cleaned, parallel = TRUE)


#dds_all$group <- factor(paste0(dds_all$condition, ddsMF$hospital))
#design(ddsMF) <- ~ group
# dds_all <- DESeq(dds_all, parallel = TRUE)
resultsNames(dds_all_cleaned)

#Make the CHLA gene-wise matrix
# ================================
chla_matrix <- as.data.frame(fc_all$counts) %>%
  dplyr::select(contains("CHLA"))
chla_matrix <- as.matrix(chla_matrix)
coldata_chla <- dir("/media/thor/storage/os_exosome_pipeline/dedupe_output/picard/", "chla.*_sorted.bam$")
#anno_chla <- data.frame(rownames(counts))
condition_chla = c(rep("ctrl", 5), rep("OS", 6))
coldata_chla <- data.frame(coldata_chla, condition_chla, row.names = 1)

colnames(chla_matrix) <- rownames(coldata_chla)
#coldata <- coldata[colnames(counts),]
all(rownames(coldata_chla) == colnames(chla_matrix))

dds_chla <- DESeqDataSetFromMatrix(countData = chla_matrix,
                                colData = coldata_chla,
                                design = ~ condition_chla)

dds_chla <- dds_chla[ rowSums(DESeq2::counts(dds_chla)) > 1, ]
dds_chla <- DESeq(dds_chla, parallel = TRUE)

#Make the HLOH gene-wise matrix
# ================================
hloh_matrix <- as.data.frame(fc_all$counts) %>%
  dplyr::select(contains("hloh"))
hloh_matrix <- as.matrix(hloh_matrix)
coldata_hloh <- dir("/media/thor/storage/os_exosome_pipeline/dedupe_output/picard/", "HLOH.*_sorted.bam$")
#anno_hloh <- data.frame(rownames(counts))
condition_hloh = c(rep("OS", 3), rep("ctrl", 7), rep("OS", 3))
coldata_hloh <- data.frame(coldata_hloh, condition_hloh, row.names = 1)

colnames(hloh_matrix) <- rownames(coldata_hloh)
#coldata <- coldata[colnames(counts),]
all(rownames(coldata_hloh) == colnames(hloh_matrix))

dds_hloh <- DESeqDataSetFromMatrix(countData = hloh_matrix,
                                   colData = coldata_hloh,
                                   design = ~ condition_hloh)

dds_hloh <- dds_hloh[ rowSums(DESeq2::counts(dds_hloh)) > 1, ]
dds_hloh <- DESeq(dds_hloh, parallel = TRUE)

#Make the HLOH cleaned gene-wise matrix
# ================================
remove <- c("N25", "N21", "N28", "OS27", "OS11")
hloh_cleaned_matrix <- as.data.frame(fc_all$counts) %>%
  dplyr::select(contains("hloh")) 
hloh_cleaned_matrix <- hloh_cleaned_matrix[, grep(paste(remove,collapse="|"), colnames(hloh_cleaned_matrix), invert = TRUE, value = TRUE)]
hloh_cleaned_matrix <- as.matrix(hloh_cleaned_matrix)
coldata_hloh_cleaned <- dir("/media/thor/storage/os_exosome_pipeline/dedupe_output/picard/", "HLOH.*_sorted.bam$")
coldata_hloh_cleaned <- grep(paste(remove,collapse="|"), coldata_hloh_cleaned, invert = TRUE, value = TRUE)
#anno_hloh_cleaned <- data.frame(rownames(counts))
condition_hloh_cleaned = c(rep("OS", 2), rep("ctrl", 4), rep("OS", 2))
coldata_hloh_cleaned <- data.frame(coldata_hloh_cleaned, condition_hloh_cleaned, row.names = 1)

colnames(hloh_cleaned_matrix) <- rownames(coldata_hloh_cleaned)
#coldata <- coldata[colnames(counts),]
all(rownames(coldata_hloh_cleaned) == colnames(hloh_cleaned_matrix))

dds_hloh_cleaned <- DESeqDataSetFromMatrix(countData = hloh_cleaned_matrix,
                                   colData = coldata_hloh_cleaned,
                                   design = ~ condition_hloh_cleaned)

dds_hloh_cleaned <- dds_hloh_cleaned[ rowSums(DESeq2::counts(dds_hloh_cleaned)) > 1, ]
dds_hloh_cleaned <- DESeq(dds_hloh_cleaned, parallel = TRUE)


res_condition <- results(ddsMF, contrast = c("condition", "OS", "ctrl"), parallel = TRUE)
res_hospital <- results(ddsMF, contrast = c("hospital", "HLOH", "CHLA"), parallel = TRUE)
res_all <- results(dds_all, contrast=c("condition", "OS", "ctrl"), parallel = TRUE)
res_all_cleaned <- results(dds_all_cleaned, contrast=c("condition_all_cleaned", "OS", "ctrl"), parallel = TRUE)
res_chla <- results(dds_chla, contrast=c("condition_chla", "OS", "ctrl"), parallel = TRUE)
res_hloh <- results(dds_hloh, contrast=c("condition_hloh", "OS", "ctrl"), parallel = TRUE)
res_hloh_cleaned <- results(dds_hloh_cleaned, contrast=c("condition_hloh_cleaned", "OS", "ctrl"), parallel = TRUE)

res <- res_all_cleaned #res placeholder!!!!!
dds <- dds_all_cleaned ############# dds placeholder!!!!
resOrdered <- res[order(res$padj),]
res_fold_chng_Ordered <- res[order(res$log2FoldChange),]
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)

res05 <- results(dds, alpha=0.05)
sum(res05$padj < 0.05, na.rm=TRUE)

library("IHW")
resIHW <- results(dds, filterFun=ihw)
summary(resIHW)

sum(resIHW$padj < 0.1, na.rm=TRUE)
metadata(resIHW)$ihwResult

resultsNames(dds)

d <- plotCounts(dds, gene="ENST00000262215.7_1", intgroup="condition_all_cleaned", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition_all_cleaned, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

#interactive report with DESeq2Report
report <- DESeq2Report(dds_all, res = res_all, 'DESeq2-all_OS_v_ctrl', 'condition_all',
                       outdir = 'DESeq2Report_all', 
                       customCode = "/media/thor/storage/os_exosome_pipeline/dedupe_output/regionreport_child.Rmd") 
                       

# write output to file
#grp_comp = "chla_os_v_ctrl"
#subset_genes <- subset(res_fold_chng_Ordered, log2FoldChange > 0)
#top_genes <- as.data.frame(subset_genes, row.names = NULL)
#gene_sym_from_trsid(top_genes, grp_comp)


