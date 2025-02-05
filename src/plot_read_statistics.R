#!/usr/local/bin/Rscript

library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

qc_files <- list.files("~/os_exosome_pipeline/results/qc_read_count_plots/", full.names = T, pattern = ".txt")

qc_tables <- lapply(qc_files, read.table, sep = "\t", header = T, stringsAsFactors = F)

techs <- basename(gsub(".txt", "", qc_files))

clean_seq_reads <- function(qc_table, tech){
  # browser()
  cell_ids <- qc_table[,1]
  unif_cell_ids <- gsub("_.*", "", cell_ids)
  qc_table[,1] <- unif_cell_ids
  
  names(qc_table) <- c("cell", tech)
  
  return(qc_table)
}

qc_reads <- purrr::map2(qc_tables, techs, clean_seq_reads)

data <- Reduce(function(x, y) merge(x, y, all=T, by="cell"), qc_reads, accumulate=F) %>% 
  # dplyr::select(c("cell", "vanilla", "clumpify_1", "clumpify_2", "prinseq", "trimmomatic")) %>% 
  gather("technology", "counts", 2:dim(.)[2]) %>%
	dplyr::mutate(disease_group = ifelse(grepl("OS", cell), "OS", "ctrl")) %>% 
	dplyr::arrange(disease_group) %>% 
	identity()

data$cell<- ordered(data$cell, levels = unique(data$cell))
data$technology <- ordered(data$technology, levels = c("vanilla", "clumpify_1", "clumpify_2", "prinseq", "trimmomatic"))

pdf("~/os_exosome_pipeline/results/read_statistics_by_qc_tech.pdf")
g <- ggplot(data, aes(x = factor(cell), y = counts, fill=technology)) +
         geom_bar(position = position_dodge(), stat = "identity", width = .7) + 
         theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
				 scale_y_log10()
print(g)
dev.off()



