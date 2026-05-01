# Calculate gene length needed to calculate transcripts per million (tpm) and to check effects on pi and dxy
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)

setwd("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/LeaH/ChamberTempExp/results/GeneLength")

GTFfile = "~/Library/CloudStorage/OneDrive-Aarhusuniversitet/LeaH/Annotation/data/Struthio_camelus_HiC_augustus.gtf"
FASTAfile = "~/Library/CloudStorage/OneDrive-Aarhusuniversitet/LeaH/Annotation/data/Struthio_camelus_HiC.fasta"

#Load the annotation and reduce it
GTF <- import.gff(GTFfile, format="gtf", genome="Struthio_camelus_HiC", feature.type="exon")
grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
reducedGTF <- unlist(grl, use.names=T)
elementMetadata(reducedGTF)$gene_id <- rep(names(grl),  elementNROWS(grl))

#Open the fasta file
FASTA <- FaFile(FASTAfile)
open(FASTA)

#Add the GC numbers
elementMetadata(reducedGTF)$nGCs <- letterFrequency(getSeq(FASTA, reducedGTF), "GC")[,1]
elementMetadata(reducedGTF)$widths <- width(reducedGTF)

#Create a list of the ensembl_id/GC/length
calc_GC_length <- function(x) {
  nGCs = sum(elementMetadata(x)$nGCs)
  width = sum(elementMetadata(x)$widths)
  c(width, nGCs/width)
}
output <- t(sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_GC_length))
colnames(output) <- c("Length", "GC")

write.table(output, file="GC_lengths.tsv", sep="\t")


### examine output a little
hist(output[,1]) # some very long but plausible
hist(output[which(output[,1]<2000),]) # many very short ones ie < 100, almost no 100 < x < 200, then some in most bins

hist(output[,2]) # looks close to normal, all good


