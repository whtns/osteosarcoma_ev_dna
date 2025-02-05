#!/usr/bin/Rscript

library(ggplot2)
library(data.table)
library(dtplyr)
library(tidyr)

os_metrics <- fread("/home/thor/os_exosome_pipeline/dedupe_output/picard/total_metrics.tsv")
os_plot_subset <- os_metrics %>%
  select(`Sample`, `Exonic Rate`, `Intronic Rate`, `Intergenic Rate`, `Mapping Rate`, `Mapped`, `Total Purity Filtered Reads Sequenced`) %>%
  mutate(Mapped.Reads = `Mapped`, Total.Reads = `Total Purity Filtered Reads Sequenced`, Exonic.Reads = `Exonic Rate`*`Mapped`, Intronic.Reads = `Intronic Rate`*`Mapped`, Intergenic.Reads = `Intergenic Rate`*`Mapped`) %>%
  select(Sample, Total.Reads, Mapped.Reads, Intergenic.Reads, Intronic.Reads, Exonic.Reads) %>%
  arrange(desc(Total.Reads)) 
os_plot_subset_2 <- gather(os_plot_subset, 'Locus', 'Read.Counts', Total.Reads:Exonic.Reads)
os_plot_melt <- melt(os_plot_subset_2, id.var = c('Sample', 'Locus'))

desc_sample_list <- factor(os_plot_subset$Sample)
exp_sample_list <- factor(filter(os_plot_subset, grepl("OS", Sample))$Sample)
control_sample_list <- factor(filter(os_plot_subset, !grepl("OS", Sample))$Sample)
control_exp_sample_list <- factor(c(as.character(control_sample_list),as.character(exp_sample_list)))
sample_list <- desc_sample_list #change to adjust ordering of x axis
p <- ggplot(os_plot_melt, aes(x=factor(Sample), y=value, fill=factor(Locus)))+geom_bar(stat='identity', position='dodge')+
  scale_fill_discrete(name="Locus", breaks=c('Exonic.Reads','Intronic.Reads','Intergenic.Reads','Mapped.Reads','Total Purity Filtered Reads Sequenced'), labels=c('Exonic', 'Intronic', 'Intergenic', 'Mapped', 'Total'))+
  scale_x_discrete(limits = control_exp_sample_list)+
  xlab('OS Samples')+ylab('Reads')
p <- p + theme(axis.text.x=element_text(angle = -90, hjust = 0))
p
ggsave("./dedupe_output/read_count_plots_ctrl_exp_split.png", p)

