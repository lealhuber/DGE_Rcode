### Preparation of DGE analysis with glmmSeq
### Lea Huber

setwd("/Users/Lea/Library/CloudStorage/OneDrive-Aarhusuniversitet/LeaH/ChamberTempExp/data/working/RNAseq_countdata/")
# Packages
library(edgeR)
library(tidyverse)
library(ggpubr)
library(glmmSeq)

# import sample data
sampleInfo = read.csv("datRNAseq_new.csv", stringsAsFactors = TRUE, row.names = 1)
sampleInfo$chickno <- as.factor(sampleInfo$chickno)
sampleInfo$trtage <- as.factor(paste0(sampleInfo$age, ".", sampleInfo$treatment)) # maybe I'll need this
sampleInfo$Samplingtime <- as.POSIXct(sampleInfo$Samplingtime)
sampleInfo$group <- as.factor(sampleInfo$group)
# scale and center chamberT
sampleInfo$chamberTscaled <- scale(sampleInfo$chamberT)

# Read in count data
filenames <- read.csv("filenames.csv", row.names = 1)
filenames <- filter(filenames, !(sample == "L83")) # know L83 is bad from previous analysis
dgel <- readDGE(filenames, group = sampleInfo$trtage)
colnames(dgel$counts) <- sampleInfo$SampleNr # add correct sample names
# remove those meta tags from counts, so the last 5 rows
counts_raw <- dgel$counts
counts_rownr <- nrow(counts_raw)
counts_notags <- counts_raw[-c(counts_rownr,counts_rownr-1,counts_rownr-2,counts_rownr-3,counts_rownr-4),]
# put back into dgel
dgel$counts <- counts_notags
nrow(dgel$counts) # check
tail(dgel$counts) # good

### filter out Globin 
top100 <- read.csv("../../../results/DESeq_results/top100.csv") # some are in top 100 genes (made before resequencing)
# remove globin genes (four are in top 5) and ferritin (top 3)
keep <- !rownames(dgel$counts) %in% c(top100[1:5,2],"g12308","g12444","g14358","g10395","g17038")
dgel_noG <- dgel[keep, ,keep.lib.sizes = FALSE] # filter and recalculate library sizes
counts_noG <- dgel_noG$counts

hist(dgel_noG$samples$lib.size)
## two libsize outliers, remove those, but leave chicks in
libsamples <- data.frame(libsize = dgel_noG$samples$lib.size, SampleNr = dgel_noG$samples$sample)
libsamples <- libsamples[order(libsamples$libsize, decreasing = TRUE),] # orders largest to smallest
badsamples <- as.character(libsamples$SampleNr[c(1,2)]) # "L38" "L53"
sampleInfo <- droplevels(filter(sampleInfo, !(SampleNr %in% badsamples))) # remove from sample info
counts_18ch <- counts_noG[,-which(colnames(counts_noG) %in% badsamples)] # remove from counts
dgel_noG <- DGEList(counts = counts_18ch, samples = sampleInfo$SampleNr, group = sampleInfo$trtage) # remake dgel object


## normal filtering
keep3 <- filterByExpr(dgel_noG$counts, group = dgel_noG$samples$group)
# only genes > 10 reads in 10 samples kept (or genes with CPM >= 0.2845027 in 10 samples), and > 15 reads over all
table(keep3) # 11082 genes kept
dgel_f3 <- dgel_noG[keep3, ,keep.lib.sizes = FALSE] # filter and recalculate library sizes
# visualise
AveLogCPM <- aveLogCPM(dgel_f3) 
hist(AveLogCPM, breaks = c(-7:19))
# stricter filtering makes a nicer histogram with just one peak

# save raw counts separately for easier handling
counts <- dgel_f3$counts

## normalisation
dgel_f3 <- normLibSizes(dgel_f3) # default is TMM (trimmed mean of M-values)
# save normalisation factors
sizeFactors <- dgel_f3$samples$norm.factors
sizeFactors <- sizeFactors/mean(sizeFactors)
names(sizeFactors) <- dgel_f3$samples$sample # better to add sample names
sampleInfo$sizeFactorsf3 <- sizeFactors
sampleInfo$libsizef3_new<- dgel_f3$samples$lib.size
plot(dgel_f3$samples$lib.size, dgel_f3$samples$norm.factors) # looks better
ggplot(sampleInfo, aes(x = treatment, y = libsizef3_new))+
  geom_boxplot(outliers = FALSE)+
  geom_point()+
  #geom_smooth(method = "lm")+
  theme_pubr()
# hot does have a bit smaller size factors


## calculate dispersion
# try out this outlier robust method
# it's less strict but handles outliers better, which we might have because of the globin thing and it being an unusual tissue
disp  <- setNames(estimateGLMRobustDisp(dgel_f3, verbose = TRUE)$tagwise.dispersion, rownames(dgel_f3))
# are dispersions biased somehow?
hist(disp) # looks ok
table(round(disp)) # 10318   687    36    17     7     8     6     1     1     1 

# biases?
plot(disp, rowMeans(counts))
avlogexp <- rowMeans(log2(counts+1))
plot(avlogexp, sqrt(disp))


## save stuff
save(dgel_f3, counts, sampleInfo, disp, sizeFactors, file = "glmmIngredients_withchicks_rob.RData")

