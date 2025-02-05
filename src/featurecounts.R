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
  output <-  merge(as.data.frame(biomaRt_result), output, by = "ensembl_transcript_id") 
  filename = paste0(grp_comp,"_deseq2",".csv")
  write.table(output,filename, sep="\t", quote=FALSE, row.names = FALSE)
}

GTF="/media/thor/storage/os_exosome_pipeline/gencode.v26lift37.annotation.gtf"
EXPTNAME=exosome_rna
CPUS=8
MAPQ=10

#Make the gene-wise matrix
fls <- dir("/media/thor/storage/os_exosome_pipeline/dedupe_output/picard/", "with_rg.bam$", full.names = TRUE)
coldata <- dir("/media/thor/storage/os_exosome_pipeline/dedupe_output/picard/", "with_rg.bam$")
anno <- data.frame(rownames(counts))
condition = c(rep("OS", 3), rep("ctrl", 7), rep("OS", 1), rep("ctrl", 5), rep("OS", 8))
hospital = c(rep("HLOH", 3), rep("CHLA", 1), rep("HLOH", 1), rep("CHLA", 4), rep("HLOH", 1), rep("CHLA", 1), rep("HLOH", 5), rep("CHLA", 5), rep("HLOH", 3))
coldata <- data.frame(coldata, condition, hospital, row.names = 1)
#fc_all <- featureCounts(fls, annot.ext = GTF, 
#                        isGTFAnnotationFile = TRUE,
#                        GTF.featureType="transcript",
#                        ignoreDup=TRUE,
#                        GTF.attrType="transcript_id") 




counts <- fc_all$counts
colnames(counts) <- rownames(coldata)
#coldata <- coldata[colnames(counts),]
all(rownames(coldata) == colnames(counts))

coldata$hospital <- factor(coldata$hospital, levels=c("CHLA"))
coldata$hospital <- droplevels(coldata$hospital)
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

ddsMF <- ddsMF[ rowSums(counts(dds)) > 1, ]

coldata$hospital <- relevel(coldata$hospital, "CHLA")
dds_chla <- DESeqDataSetFromMatrix(countData = counts,
                                colData = coldata,
                                design = ~ condition)

ddsMF_chla <- ddsMF[, ddsMF$hospital %in% c("CHLA")]
ddsMF_chla$hospital <- droplevels(ddsMF_chla$hospital)
ddsMF_chla <- DESeq(ddsMF_chla, contrast)


#dds$condition <- factor(dds$condition, levels=c("OS","ctrl"))
ddsMF <- DESeq(ddsMF, parallel = TRUE)

ddsMF_chla <- ddsMF
ddsMF_chla$hospital <- factor(ddsMF_chla$hospital, levels=c("CHLA"))
ddsMF_chla$hospital <- droplevels(ddsMF_chla$hospital)
ddsMF_chla$group <- factor(paste0(ddsMF_chla$condition, ddsMF_chla$hospital))
design(ddsMF_chla) <- ~ group
ddsMF_chla <- DESeq(ddsMF_chla, parallel = TRUE)
resultsNames(ddsMF_chla)
resMf_chla <- results(ddsMF_chla, contrast=c("group", "OSCHLA", "ctrlCHLA"))

res_condition <- results(ddsMF, contrast = c("condition", "OS", "ctrl"), parallel = TRUE)
res_hospital <- results(ddsMF, contrast = c("hospital", "HLOH", "CHLA"), parallel = TRUE)
res_condition_CHLA <- results(ddsMF_chla, contrast = c("condition", "OS", "ctrl"), parallel = TRUE)
res <- resMf_chla
resOrdered <- res[order(res$padj),]
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)

res05 <- results(ddsMF, alpha=0.05)
sum(res05$padj < 0.05, na.rm=TRUE)

library("IHW")
resIHW <- results(ddsMF, filterFun=ihw)
summary(resIHW)

sum(resIHW$padj < 0.1, na.rm=TRUE)
metadata(resIHW)$ihwResult

resultsNames(ddsMF)

d <- plotCounts(ddsMF, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

#interactive report with DESeq2Report
report <- DESeq2Report(ddsMF_chla, 'DESeq2-CHLA_OS_v_ctrl', 'group',
                       outdir = 'DESeq2Report-example')

#Interactive MA plot with Glimma
glMDPlot <- glMDPlot(res, counts, anno, groups, samples = NULL,
         status = rep(0, nrow(x)), transform = FALSE, xlab = "Mean Expression",
         ylab = "log-fold-change", side.xlab = "Group", side.ylab = "Expression",
         side.log = FALSE, side.gridstep = ifelse(!transform || side.log, FALSE,
                                                  0.5), jitter = 30, side.main = "GeneID", display.columns = NULL,
         cols = c("#00bfff", "#858585", "#ff3030"), sample.cols = rep("#1f77b4",
                                                                      ncol(counts)), path = getwd(), folder = "glimma-plots",
         html = "MD-Plot", launch = TRUE, ...)

