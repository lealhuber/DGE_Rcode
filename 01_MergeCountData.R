#### Make merged count data files for DGE analysis
### Lea Huber, November 2024

setwd("/Users/Lea/Library/CloudStorage/OneDrive-Aarhusuniversitet/LeaH/ChamberTempExp/data/working/RNAseq_countdata/")
library(tidyverse)
library(readxl)

# prepare SampleInfo df
sampleInfo = read.csv("datRNAseq_63.csv", stringsAsFactors = TRUE)
sampleInfo$SampleNr <- paste0("L", sampleInfo$SampleNr) # add the L again in sample names to make it definitely not a number
sampleInfo <- arrange(sampleInfo, SampleNr) # order by sample computer-style
rownames(sampleInfo) <- sampleInfo$SampleNr # might come in handy
sampleInfo$chickno <- as.factor(sampleInfo$chickno)
## add extraction stuff
extraction <- read_excel("../../../docs/RNAextraction/RNAlater_66_overview.xlsx", range = "C3:G65", col_names = c("SampleNr","Conc","r260v280","r260v230","tot_ng"))
sampleInfo <- left_join(sampleInfo, extraction, by = "SampleNr")

## Add info whether it will be merged or not
# get new sample names
newsamples <- read.table("Dedup_resequenced/SampleNames_dedup.txt")
newsamples <- newsamples$V1
newsamples <- str_remove(newsamples, ".count_matrix.txt")
# Add to sample info
sampleInfo$merged <- ifelse(sampleInfo$SampleNr %in% newsamples, "yes", "no")

### Merge counts of same sample
# get the count files in order
countmatrixnames1 = read.table("Firstround/matrixFileNames.txt", col.names = c("files","sample"))
# need to reorder it so sample number matches filename
countmatrixnames1$sample = str_sub(countmatrixnames1$files, start = 2, end = 3)
countmatrixnames1$sample = as.numeric(str_remove(countmatrixnames1$sample, "[.]"))
countmatrixnames1 <- arrange(countmatrixnames1, sample)
# add path in front
countmatrixnames1$files <- paste0("Firstround/", countmatrixnames1$files)

countmatrixnames2 <- read.table("Dedup_resequenced/matrixFileNames_dedup.txt", col.names = c("files","sample"))
# need to reorder it so sample number matches filename
countmatrixnames2$sample = str_sub(countmatrixnames2$files, start = 2, end = 3)
countmatrixnames2$sample = as.numeric(str_remove(countmatrixnames2$sample, "[.]"))
countmatrixnames2 <- arrange(countmatrixnames2, sample)
# add path in front
countmatrixnames2$files <- paste0("Dedup_resequenced/", countmatrixnames2$files)

# combine filenames old df after the new, by sample number
countmatrixnames <- rbind(countmatrixnames2, countmatrixnames1)
countmatrixnames <- arrange(countmatrixnames, sample)
countmatrixnames$sample <- as.factor(paste0("L", countmatrixnames$sample))

newfilenames <- c()
for (sampleNr in levels(countmatrixnames$sample)){
  ### merge count data from old and new sequencing run of same library
  files <- unlist(subset(countmatrixnames, sample == sampleNr, select = files))
  if (length(files)==2){ # only if there even are 2
    old <- read.table(files[1],col.names = c("gene","countold"))
    new <- read.table(files[2],col.names = c("gene","countnew"))
    # are gene names all the same?
    if (all(old$gene == new$gene)){
      new$countold <- old$countold
      new <- mutate(new, count = countnew + countold)
      counts <- subset(new, select = c(gene, count))
      write.table(counts, file = paste0("Dedup_combined/", sampleNr, ".dedup.countmatrix.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
      newfilenames <- c(newfilenames, paste0("Dedup_combined/",sampleNr, ".dedup.countmatrix.txt"))
    }else{
      echo(paste0("Error: Sample ",sampleNr," has not all same genes."))
    }
  }else{
    counts <- read.table(files, col.names = c("gene","countold"))
    write.table(counts, file = paste0("Dedup_combined/",sampleNr, ".dedup.countmatrix.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    newfilenames <- c(newfilenames, paste0("Dedup_combined/",sampleNr, ".dedup.countmatrix.txt"))
  }
}

filenames <- data.frame(files = newfilenames, sample = levels(countmatrixnames$sample))

## add duplication levels
stats_old <- read.table("../Processing_stats/general_stats_old.tsv", header = TRUE)
stats_new <- read.table("../Processing_stats/general_stats_new.tsv", header = TRUE, na.strings = ".")
# collapse to one row per sample
stats_old$shortsmpl <- str_remove(str_sub(stats_old$Sample, start = 1, end = 3), "_")
stats_new$shortsmpl <- str_remove(str_sub(stats_new$Sample, start = 1, end = 3), "_")
# summarise duplication level (mean)
stats_old_dups <- summarise(group_by(stats_old, shortsmpl), avgdup_old = mean(mqc.generalstats.fastqc.percent_duplicates))
stats_new_dups <- summarise(group_by(filter(stats_new, !is.na(Dups)), shortsmpl), avgdup_new = mean(Dups))
# add to sampleInfo
sampleInfo <- left_join(sampleInfo, stats_old_dups, join_by(SampleNr == shortsmpl))
sampleInfo <- left_join(sampleInfo, stats_new_dups, join_by(SampleNr == shortsmpl))


# save files
write.csv(filenames, "filenames_dedup.csv")
write.csv(sampleInfo, "datRNAseq_new.csv")

## changed files structure on 12.12.24 so have to redo old filenames file
filenames <- read.csv("filenames.csv", row.names = 1)
filenames$files <- paste0("Combined/",filenames$files)
write.csv(filenames, "filenames.csv")
