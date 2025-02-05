# count reads using SAM tools  
library(ggplot2)  
library(reshape2)  
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
bam_feed_script = args[0]
# function to call BAMtools and store reads count in a text file  
CountReadsBAM <- function()  
{  
  command= "~/bin/CountMappedReads.sh > MappedReads.txt"  
  try(system(command))  
  res=read.table("MappedReads.txt", header=T)  
  return(res)  
}

# calling the funtion will output read count in an object:  
ReadCount <- CountReadsBAM()

# a bit of polishing to remove total (not needed for plotting)  
ReadCountSmall <- data.frame(BAM = ReadCount$bam_file, Mapped = ReadCount$mapped, Unmapped = ReadCount$unmapped)

# ggplot needs data in a specific layout  
MeltedReadCount = melt(ReadCountSmall, id=c('BAM'))  
names(MeltedReadCount) <- c('BAM', 'Mapping', 'Reads')

# calculate the fraction of mapped reads  
ReadsFraction <- ddply(  
  MeltedReadCount,  
  .(BAM),  
  summarise,  
  Count.Fraction = Reads / sum(Reads)  
)

# sort the data frame and add fraction to the data frame  
to_graph <- cbind(arrange(MeltedReadCount, BAM), fraction = ReadsFraction$Count.Fraction)

# Now all we have to do is plot the data  
gp <- ggplot(data=to_graph, aes(x=BAM, y=Reads, fill=Mapping)) +  
  geom_bar(stat="identity",position = "stack") +  
  geom_text(aes(label=paste(round(fraction*100),"%", sep="")),size = 3,vjust=0,position="stack") +  
  opts(axis.text.x=theme_text(angle=90))

pdf("MappedReads.pdf")  
gp  
dev.off() 