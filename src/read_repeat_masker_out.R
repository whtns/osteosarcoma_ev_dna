#!/usr/local/bin/Rscript


# load required libraries -------------------------------------------------

library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(data.table)


# load required functions -------------------------------------------------

strip_plot_by_family <- function(rep_table){
	strip_plot <- ggplot(rep_table, aes(x=genome_proportion, y=disease_group)) +  
		geom_point(aes(shape=disease_group, color=factor(sample_by_group), group = disease_group)) + facet_wrap(~repeat_class_family, ncol =1, strip.position = "left") + 
		theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), strip.text.y = element_text(angle=180), strip.text.x = element_blank(), 
					panel.border = element_rect(colour = "black", fill=NA, size=0.25), strip.background = element_rect(color = "black", size = 0.25), 
					panel.spacing = unit(0, "lines")) +
		scale_colour_manual(name = "Exosome Sample", labels = unique(rep_table$sample), values = getPalette(colourLength)) 
	print(strip_plot)	
}

strip_plot_by_repeat <- function(rep_table){
	strip_plot <- ggplot(rep_table, aes(x=genome_proportion, y=disease_group)) +  
		geom_point(aes(shape=disease_group, color=factor(sample_by_group), group = disease_group)) + facet_wrap(~matching_repeat, ncol =1, strip.position = "left") + 
		theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), strip.text.y = element_text(angle=180), strip.text.x = element_blank(), 
					panel.border = element_rect(colour = "black", fill=NA, size=0.25), strip.background = element_rect(color = "black", size = 0.25), 
					panel.spacing = unit(0, "lines")) +
		scale_colour_manual(name = "Exosome Sample", labels = unique(rep_table$sample), values = getPalette(colourLength)) 
	print(strip_plot)	
}

full_join_NA <- function(x, y, ...) {
  full_join(x = x, y = y, by = ...) %>% 
    mutate_all(funs(replace(., which(is.na(.)), 0)))
}

right_join_NA <- function(x, y, ...) {
  dplyr::right_join(x = x, y = y, by = ...) %>% 
    mutate_all(funs(replace(., which(is.na(.)), 0)))
}

n_fun <- function(x){
  return(data.frame(y = median(x), label = paste0("n = ",length(x))))
}


# load data ---------------------------------------------------------------

repeatmask_out_files = list.files("~/os_exosome_pipeline/output/repeatmasker/", full.names = TRUE, pattern = ".*.trimmed.fasta.out$")
rmask_summs <- lapply(repeatmask_out_files, function(x)read.table(x, header = FALSE, skip =2, sep = "", stringsAsFactors = FALSE, fill=TRUE))
names(rmask_summs) <- gsub(".*//|_S.*", "", repeatmask_out_files)

genome_sizes <- read.table("~/os_exosome_pipeline/results/prinseq_sequence_size.csv", sep = ",", header = FALSE)
colnames(genome_sizes) <- c("sample", "genome_size")

#repeatmask_out <- read.table(repeatmask_out_file, header = FALSE, skip =2, sep = "", stringsAsFactors = FALSE, fill=TRUE)

rep_mask_colnames <- c("SW_score", "perc_div.", "perc_del.", "perc_ind.", "query_sequence", "query_begin", 
                             "query_end", "query_left", "num_comp","matching_repeat", "repeat_class_family", "repeat_begin", 
                             "repeat_end", "repeat_left", "ID")
rmask_summs <- lapply(rmask_summs, setNames, rep_mask_colnames)

# tidy data ---------------------------------------------------------------

tidy_rmask <- rmask_summs %>%
  data.table::rbindlist(use.names=TRUE, fill=TRUE, idcol="sample") %>%
  mutate(Size = query_end - query_begin) %>%
  mutate(disease_group = ifelse(grepl("OS", sample), "OS","ctrl")) %>%
  left_join(genome_sizes, by = "sample")


# prep data for plotting --------------------------------------------------

zero_template = read.table("~/os_exosome_pipeline/results/sarcoma_zero_template.csv", sep = ",")

strip_plot.grouped <- tidy_rmask %>%
	group_by(disease_group, sample, genome_size, repeat_class_family, matching_repeat) %>%
	summarise(Count = n(), Size = sum(Size)) %>%
	mutate(genome_proportion = Size/genome_size, 
				 sample_by_group = paste(sample, disease_group, sep="_"),
				 sample_by_family = paste(sample, repeat_class_family, sep="_")) 
#  filter(!is.na(Size))

strip_plot_zeroes <- tidy_rmask %>%
  group_by(disease_group, sample, genome_size, repeat_class_family) %>%
  summarise(Count = n(), Size = sum(Size)) %>%
  mutate(genome_proportion = Size/genome_size, 
         sample_by_group = paste(sample, disease_group, sep="_"),
         sample_by_family = paste(sample, repeat_class_family, sep="_")) %>% 
  right_join_NA(zero_template, by=c("sample_by_family", "sample_by_group", "repeat_class_family", "sample", "disease_group"))

box_plot_zeroes <- strip_plot.grouped %>%
  group_by(repeat_class_family) %>%
  summarize(maxSize = max(Size)) %>%
  mutate(maxSize = maxSize ) %>% 
  left_join(strip_plot.grouped, by="repeat_class_family") %>%
  merge(zero_template, by=c("sample_by_family", "sample", "disease_group", "repeat_class_family", "sample_by_group"), all = TRUE)

box_plot_zeroes[is.na(box_plot_zeroes)] <- 0

box_plot_zeroes <- box_plot_zeroes %>%
  group_by(repeat_class_family) %>%
  mutate(maxSize  = maxSize + 1) %>%
  mutate(class_proportion = Size/maxSize)
  
box_plot.grouped <- strip_plot.grouped %>%
  group_by(repeat_class_family) %>%
  summarize(maxSize = max(Size)) %>%
  mutate(maxSize = maxSize) %>%
  left_join(strip_plot.grouped, by="repeat_class_family") %>%
  group_by(repeat_class_family) %>%
  mutate(class_proportion = Size/maxSize) %>%
  filter(!is.na(maxSize))
# filter(Size >= 1000) 
# filter(repeat_class_family != "rRNA" & repeat_class_family != "LINE/L1")

sine_alu_prop <- strip_plot.grouped[strip_plot.grouped$repeat_class_family == "SINE/Alu",]
line_l1_prop <- strip_plot.grouped[strip_plot.grouped$repeat_class_family == "LINE/L1",]

# do plotting -------------------------------------------------------------

colourCount = unique(strip_plot.grouped$sample_by_group)
colourLength = length(colourCount)
shapeCount = c(rep(16, 12), rep(17, 12))
shapeCount <- gsub(".*OS.*", "OS", gsub(".*ctrl", "ctrl", colourCount))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

# list of plots for repeatmasker analysis
p_l = list()

fam_strip_plot <- strip_plot_by_family(strip_plot.grouped)

sine_rep_strip_plot <- strip_plot_by_repeat(sine_alu_prop)
line_rep_strip_plot <- strip_plot_by_repeat(line_l1_prop)

p_l[[1]] = sine_rep_strip_plot
p_l[[2]] = line_rep_strip_plot

pdf("~/os_exosome_pipeline/results/sine_and_line_proportions_os_exosome.pdf", height = 12)
p_l
dev.off()



colourCount = unique(strip_plot_zeroes$sample_by_group)
colourLength = length(colourCount)
shapeCount = c(rep(16, 12), rep(17, 12))
shapeCount <- gsub(".*OS.*", "OS", gsub(".*ctrl", "ctrl", colourCount))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

strip_plot.zeros <- ggplot(strip_plot_zeroes, aes(x=genome_proportion, y=disease_group)) +  
  geom_point(aes(shape=disease_group, color=factor(sample_by_group), group = disease_group)) + facet_wrap(~repeat_class_family, ncol =1, strip.position = "left") + 
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), strip.text.y = element_text(angle=180), strip.text.x = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=0.25), strip.background = element_rect(color = "black", size = 0.25), 
        panel.spacing = unit(0, "lines")) +
  scale_colour_manual(name = "Exosome Sample", labels = unique(strip_plot_zeroes$sample), values = getPalette(colourLength)) 
print(strip_plot.zeros)
p_l[[2]] = strip_plot.zeros

ggsave(file="./results/fixed_stip_plot.png", width=7, height =10, dpi = 200)

box_plot.g <- ggplot(box_plot.grouped, aes(x=repeat_class_family, y=genome_proportion, fill=disease_group)) +
      geom_boxplot(outlier.shape = NA) + coord_flip()  +
  stat_summary(fun.data = n_fun, aes(x=factor(repeat_class_family), fill = factor(disease_group)), position=position_dodge(.9),
                                                                     geom = "text", vjust = +1, size = 4)
print(box_plot.g)
p_l[["box_plot"]] = box_plot
p_l[["box_plot.g"]] = box_plot.g

box_plot.zeros.g <- ggplot(box_plot_zeroes, aes(x=repeat_class_family, y=genome_proportion, fill=disease_group)) +
  geom_boxplot(outlier.shape = NA) + coord_flip() 
#  stat_summary(fun.data = n_fun, aes(x=factor(repeat_class_family), fill = factor(disease_group)), position=position_dodge(.9),
 #              geom = "text", vjust = +1, size = 4)
print(box_plot.zeros)
p_l[["box_plot.zeros"]] = box_plot.zeros
p_l[["box_plot.zeros.g"]] = box_plot.zeros.g

lapply(names(p_l), function(x){ ggsave(file=paste("./results/",x, ".png"), p_l[[x]], width=7, height=10, dpi=200)})



# find housekeeping repetivie elements ------------------------------------

# find 'housekeeping' repetitive element, exhibiting high and uniform expression despite disease_group

hk_in <- strip_plot.grouped %>% 
	group_by(disease_group, repeat_class_family, matching_repeat) %>% 
	summarize(values = list(genome_proportion)) %>% 
	group_by(disease_group, matching_repeat) %>% 
	tidyr::spread(disease_group, values) %>% 
	filter(!is.null(unlist(OS))) %>% 
	filter(!is.null(unlist(ctrl))) %>% 
	filter(length(unlist(ctrl)) > 1) %>% 
	filter(length(unlist(OS)) > 1) %>% 
	filter(!any(is.na(unlist(ctrl)))) %>% 
	filter(!any(is.na(unlist(OS)))) %>% 
	mutate(p_value = equivalence::tost(unlist(OS), unlist(ctrl))$tost.p.value) %>%
	mutate(OS_prop = equivalence::tost(unlist(OS), unlist(ctrl))$estimate[[1]]) %>%
	mutate(ctrl_prop = equivalence::tost(unlist(OS), unlist(ctrl))$estimate[[2]]) %>%
	mutate(FC = OS_prop/ctrl_prop) %>% 
	arrange(desc(OS_prop)) %>% 
	dplyr::select(-c(ctrl, OS)) %>% 
	identity()	

write.csv(hk_in, "~/os_exosome_pipeline/results/houskeeping_repetitive_elements.csv")


# run pam -----------------------------------------------------------------



mr_pamr <- as.data.frame(strip_plot.grouped) %>%
  ungroup() %>%
  dplyr::select(sample, matching_repeat, genome_proportion) %>%
  spread(sample, genome_proportion)

rpcf_pamr <- as.data.frame(strip_plot.grouped) %>%
  group_by(sample, repeat_class_family) %>% 
  summarise(genome_proportion = sum(genome_proportion)) %>%
  spread(sample, genome_proportion)

pamr_df <- as.data.frame(rpcf_pamr)
pamr_df <- data.frame(pamr_df[,-1], row.names=pamr_df[,1])

pamr_factors <- as.factor(ifelse(grepl("^.*OS.*$", colnames(pamr_df)), "OS", "ctrl"))

pamr_df[is.na(pamr_df)] <- 0

pamr_df <- as.matrix(pamr_df)


x <- pamr_df
y <- pamr_factors
mydata <- list(x=x,y=factor(y), geneid=rownames(x),
               genenames=rownames(x))
mytrain <-   pamr.train(mydata)

gene_list <- pamr.listgenes(mytrain, mydata, threshold=1.6)

pamr.geneplot(mytrain, mydata, threshold=3)
mycv <- pamr.cv(mytrain,mydata)
pamr.plotcen(mytrain, mydata,threshold=1.6)
pamr.plotcvprob(mycv,mydata,threshold=1.6)
pamr.plotcv(mycv)
myfdr <-  pamr.fdr(mytrain, mydata)
pamr.plotfdr(myfdr)
