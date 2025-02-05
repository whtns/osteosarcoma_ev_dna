#!/usr/bin/Rscript

library(dplyr)
library(ggplot2)

test <- read.table("~/os_exosome_pipeline/data/FASTQ/vanilla_read_counts_from_fastqc.txt", header = TRUE) %>% 
  rename(read_count = "vanilla") %>% 
  mutate(sample_id = gsub("_.*", "", cell)) %>% 
  select(-cell)

id.vars <- c("c1", "c2", "fastqmcf", "trimmomatic", "prinseq")
measure.vars <- c("read_count", "read_count.1", "read_count.2", "read_count.3", "read_count.4")

test0 <- read.csv("~/os_exosome_pipeline/output/clumpify/summary_fastqc_read.txt") %>% 
  mutate(sample_id = gsub("_.*", "", c1)) %>% 
  dplyr::select_(.dots = c(measure.vars, "sample_id")) %>% 
  rename_at(vars(measure.vars), ~ id.vars) %>% 
  inner_join(test, by = "sample_id") %>% 
  gather("technology", "read_counts", -sample_id)

write.csv(test0, "~/os_exosome_pipeline/doc/filtered_read_counts.csv")

test0 <- read.csv("~/os_exosome_pipeline/doc/filtered_read_counts.csv")

ggplot(test0, aes(x=sample_id, y=read_counts, fill=technology)) + 
  geom_col(position = position_dodge()) +
  theme_set(theme_gray(base_size = 20)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave("os_exosome_read_filtering.png")


