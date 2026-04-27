### Gene Set Enrichment Analysis
### Lea Huber
### January 2025

setwd("/Users/Lea/Library/CloudStorage/OneDrive-Aarhusuniversitet/LeaH/ChamberTempExp/data/working/RNAseq_countdata/")

# packages
library(tidyverse)
library(ggpubr)
library(gplots)
library(ggrepel)
library(topGO)
library(KEGGREST)
library(clusterProfiler)
library(org.Gg.eg.db)
library(ggarchery)

figurepath = "../../../figures/GSEA/glmmSeq_Lund/"

tempcols = c(benign = "#4DAF4A", cold = "#377EB8", hot = "#E41A1C")
agecols = c("1week" = "#1B9E77" , "8week" = "#D95F02")
sexcols = c("M" = "#7570B3", "F" = "#E7298A")

### set up go analysis
# go annotation path (easier to use in topGO apparently)
GO_path <- "../Struthio_camelus_HiC.emapper.annotations"
GO_annotation <- read_tsv(GO_path, comment = "##", na = "-")
# remove transcript information since I am not analysing those separately
GO_annotation$query <- str_remove(GO_annotation$query, ".t\\d" )
GO_annotation <- GO_annotation[!duplicated(GO_annotation$query), ] # each query only once

#### Pathway Over-representation -----------------------------------------------------------------------------------------------

# load the gene sets
load("NewPHgenesForGO.RData")
# for age have to make them first
load("cat_eff_age_resp.RData")

## for age do the same categories as targeted Fst
gw1_hot <- filter(cat_eff2, age_eff == "just_1week" & condition == "hvb")$gene
gw8_hot <- filter(cat_eff2, age_eff == "just_8week" & condition == "hvb")$gene
gsamedir_hot <- filter(cat_eff2, age_eff == "same_direction" & condition == "hvb")$gene
goppdir_hot <- filter(cat_eff2, age_eff == "opposite_direction" & condition == "hvb")$gene

gw1_cold <- filter(cat_eff2, age_eff == "just_1week" & condition == "cvb")$gene
gw8_cold <- filter(cat_eff2, age_eff == "just_8week" & condition == "cvb")$gene
gsamedir_cold <- filter(cat_eff2, age_eff == "same_direction" & condition == "cvb")$gene
goppdir_cold <- filter(cat_eff2, age_eff == "opposite_direction" & condition == "cvb")$gene
# some will have a bit too little

## other option: separate by age and treatment
cat_eff2_age <- filter(cat_eff2, ageeff != "same_slope") # only the ones that change with age
gw1_hot_up <- filter(cat_eff2_age, condition == "hvb" & week1 < 0)$gene
gw1_hot_down <- filter(cat_eff2_age, condition == "hvb" & week1 > 0)$gene
gw1_cold_up <- filter(cat_eff2_age, condition == "cvb" & week1 < 0)$gene
gw1_cold_down <- filter(cat_eff2_age, condition == "cvb" & week1 > 0)$gene

gw8_hot_up <- filter(cat_eff2_age, condition == "hvb" & week8 < 0)$gene
gw8_hot_down <- filter(cat_eff2_age, condition == "hvb" & week8 > 0)$gene
gw8_cold_up <- filter(cat_eff2_age, condition == "cvb" & week8 < 0)$gene
gw8_cold_down <- filter(cat_eff2_age, condition == "cvb" & week8 > 0)$gene

# checks
length(unique(gw8_cold_up))
duplicated(c(gw8_hot_up,gw8_hot_down))
duplicated(c(gw8_cold_down, gw8_cold_up))
# all good

### use gene symbols (preferred_name) to link to chicken entrezIDs
Ggentrez = bitr(GO_annotation$Preferred_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Gg.eg.db")
# 11% failed to map but the rest did go well, that's good enough
GO_annotation <- left_join(GO_annotation, Ggentrez, by = join_by(Preferred_name == SYMBOL))
### this should now directly work for analysis

# set universe to only the genes I analysed
load("glmmSeq_TchangeMods_andr.RData")
count_genes <- dimnames(fitfull_gr@countdata)$Tags
gene_universe <- count_genes[!(count_genes %in% errgs_2cat2lin)]
enrez_universe <- filter(GO_annotation, query %in% gene_universe)$ENTREZID # get entrez ID of those


# make gene lists with entrezID
gene_list <- hotup$gene # set to the one you want

entrezIDs <- GO_annotation[which(GO_annotation$query %in% gene_list),"ENTREZID"]
entrezIDs <- unique(entrezIDs$ENTREZID)
entrezIDs <- entrezIDs[!is.na(entrezIDs)]

### analyses

## KEGG enrichment
kegg_enrich <- enrichKEGG(gene = entrezIDs, organism = "gga", pvalueCutoff = 0.1, universe = enrez_universe) # find enriched KO terms
kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Gg.eg.db,keyType="ENTREZID")
# View results
head(kegg_enrich@result, n = 15)
# visualization
dotplot(kegg_enrich, showCategory=10)
# view pathway
browseKEGG(kegg_enrich,"gga00190")

## GO enrichment
GO_enrich <- enrichGO(gene = entrezIDs, OrgDb = org.Gg.eg.db, ont = "BP",pvalueCutoff = 0.1, readable = TRUE,, universe = enrez_universe) # find biological processes
head(GO_enrich@result, n =15)
#terms <- GO_enrich@result$ID[1:10]
# visualise
dotplot(GO_enrich, showCategory = 22, color = "qvalue")
goplot(GO_enrich) # not good
cnetplot(GO_enrich,categorySize="geneNum")
enrichplot::upsetplot(GO_enrich, n = 10)



### many GO terms in hotup, what genes are mostly involved?
gene_sym <- filter(kegg_enrich@result, qvalue < 0.1)$geneID
gene_sym <- unlist(strsplit(gene_sym,"/"))
table(gene_sym) # so it's only 8 of the 26 in hotup that appear in the top go terms, mark these for targeted fst
gene_sym_bothdown <- c(gene_sym_bothdown, gene_sym)
gene_sym_bothup <- gene_sym
gene_sym_hotdowncoldup <- gene_sym
gene_sym_hotupcolddown <- c(gene_sym_hotupcolddown,gene_sym)
gene_sym_hotdown <- gene_sym
gene_sym_hotup <- gene_sym

# save(gene_sym_bothdown, gene_sym_bothup, gene_sym_hotdown, gene_sym_hotup, gene_sym_hotdowncoldup, gene_sym_hotupcolddown, GO_annotation, file = "../functional_genes.RData")




#### Make supplementary files with DEGs including functions, and with GO results -----------------------------------------------------------------------

# load results including p-values
load("PostHoc_catInt.RData")
load("PostHoc_NoInt.RData")

### DEGs in heat
cat_eff_hot <- filter(cat_eff2, condition == "hvb")
## add columns effect_no_age, effect_1week, effect_8week, pval_ageint, pval_no_age from cat_effects and cat_effects_noint
# with int
cat_effects <- pivot_wider(cat_effects, names_from = age, values_from = c(hvb,cvb, hvb_p, cvb_p))
cat_effects <- dplyr::rename(cat_effects, effect_1week_hot = hvb_week1, effect_8week_hot = hvb_week8, effect_1week_cold = cvb_week1, effect_8week_cold = cvb_week8,
                             pval_1week_hot = hvb_p_week1, pval_1week_cold = cvb_p_week1, pval_8week_hot = hvb_p_week8, pval_8week_cold = cvb_p_week8)
cat_eff_hot <- left_join(dplyr::select(cat_eff_hot, c(gene, age_eff, condresp)), dplyr::select(cat_effects, c(gene, pval_1week_hot, pval_8week_hot, effect_1week_hot, effect_8week_hot)))

# no int
cat_effects_noint <- dplyr::rename(cat_effects_noint, effect_no_age_hot = hvb, effect_no_age_cold = cvb, pval_no_age_hot = p_hvb, pval_no_age_cold = p_cvb)
cat_eff_hot <- left_join(cat_eff_hot, dplyr::select(cat_effects_noint, c(gene, effect_no_age_hot, pval_no_age_hot)))

## add columns Description, preferred name, GOs and ENTREZID from GO_annotation
cat_eff_hot <- left_join(cat_eff_hot, dplyr::select(GO_annotation, c(query, Description, Preferred_name,ENTREZID, GOs)), by = join_by(gene == query))
# one has two different ENTREZ IDs, googling confirmed that 124416952 is the right one
cat_eff_hot <- cat_eff_hot[-which(cat_eff_hot$ENTREZID == "416989" & cat_eff_hot$gene == "g3646"),]
# remove all _hot from the ends
cat_eff_hot <- cat_eff_hot %>% rename_with(~str_remove(., '_hot$'))
# save
write.csv(cat_eff_hot, file = "../../../results/DEGs_heat.csv")

### DEGs in cold
cat_eff_cold <- filter(cat_eff2, condition == "cvb")

cat_eff_cold <- left_join(dplyr::select(cat_eff_cold, c(gene, age_eff, condresp)), dplyr::select(cat_effects, c(gene, pval_1week_cold, pval_8week_cold, effect_1week_cold, effect_8week_cold)))
cat_eff_cold <- left_join(cat_eff_cold, dplyr::select(cat_effects_noint, c(gene, effect_no_age_cold, pval_no_age_cold)))
cat_eff_cold <- left_join(cat_eff_cold, dplyr::select(GO_annotation, c(query, Description, Preferred_name,ENTREZID, GOs)), by = join_by(gene == query))
# again double same gene
cat_eff_cold <- cat_eff_cold[-which(cat_eff_cold$ENTREZID == "416989" & cat_eff_cold$gene == "g3646"),]
# remove all _cold from the ends
cat_eff_cold <- cat_eff_cold %>% rename_with(~str_remove(., '_cold$'))
# save
write.csv(cat_eff_cold, file = "../../../results/DEGs_cold.csv")



# make and save a list of all annotated DE genes, their expression category and functions
cat_effects_noint <- read.csv("cat_effects_noint.csv", row.names = 1)
cat_effects_noint <- rename(cat_effects_noint, no_age_cat = colour)
load("age_dir_cats.RData")
age_dirs <- rename(age_dirs, cold_ageint = cvb, hot_ageint = hvb)
Gene_functions <- inner_join(age_dirs, dplyr::select(GO_annotation, c(query,Description,Preferred_name)), by = join_by(gene == query))
Gene_functions <- left_join(Gene_functions,dplyr::select(cat_effects_noint, c(gene,no_age_cat)))
# write.csv(Gene_functions, file = "../Genefunctions.csv") # save
