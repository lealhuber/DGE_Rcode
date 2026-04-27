### post hoc analysis of gene expression
## Oct 2025
## Lea Huber

library(tidyverse)
library(glmmSeq)
library(emmeans)
library(ggpubr)

tempcols = c(benign = "#6800a4", cold = "#377EB8", hot = "#E41A1C")
agecols = c("week1" = "#1B9E77" , "week8" = "#D95F02")


setwd("/Users/Lea/Library/CloudStorage/OneDrive-Aarhusuniversitet/LeaH/ChamberTempExp/data/working/RNAseq_countdata/")


#### gather ingredients -----------------------------------------------------------------------

## load models
load("glmmSeq_catmods_group.RData")

### get genes that did not converge in all models
errgs_rgr <- names(fitfull_gr@errors)
errgs_lin_noint <- names(fit_lin_noint@errors)
errgs_cat <- names(fit_cat@errors)
errgs_cat_noint <- names(fit_cat_noint_par@errors)
errgs_2cat2lin <- unique(c(errgs_rgr, errgs_cat, errgs_cat_noint,errgs_lin_noint))

### for adding gene information
GO_path <- "../Struthio_camelus_HiC.emapper.annotations" # go annotation path (easier to use in topGO apparently)
GO_annotation <- read_tsv(GO_path, comment = "##", na = "-")
# remove transript information
GO_annotation$query <- str_remove(GO_annotation$query, ".t\\d" )
GO_annotation <- GO_annotation[!duplicated(GO_annotation$query), ] # each query only once

### get significant genes
# categorical with interaction
fit_cat <- glmmQvals(fit_cat, cutoff = 0.01)
fit_cat_par <- glmmQvals(fit_cat_par)
fit_cat_noint <- glmmQvals(fit_cat_noint, cutoff = 0.01) # is it the same? No, but similar, and better according to LRT
fit_cat_noint_par <- glmmQvals(fit_cat_noint_par)
statscat_withint <- summary(fit_cat)
statscat_withint_sig <- statscat_withint[order(statscat_withint[, 'P_treatment:age']), ] # take interaction term!!
statscat_withint_sig <- as.data.frame(statscat_withint_sig[1:395,])
colnames(statscat_withint_sig) <- c(colnames(statscat_withint_sig)[1:30], "q_treatment", "q_age", "q_sex", "q_treatment:age")
statscat_withint_sig$gene <- rownames(statscat_withint_sig)
statscat_withint_sig <- filter(statscat_withint_sig, !(gene %in% errgs_2cat2lin)) # remove the ones that don't converge for other models I will use
statscat_withint_sig <- left_join(statscat_withint_sig, subset(GO_annotation, select = c("query","Preferred_name","Description")), by = join_by(gene == query))
rownames(statscat_withint_sig) <- statscat_withint_sig$gene
# categorical with no interaction
statsFull_cat <- summary(fit_cat_noint)
statsFull_cat_sig <- statsFull_cat[order(statsFull_cat[, 'P_treatment']), ]
statsFull_cat_sig <- as.data.frame(statsFull_cat_sig[1:472,])
colnames(statsFull_cat_sig) <- c(colnames(statsFull_cat_sig)[1:23], "q_treatment", "q_age", "q_sex")
statsFull_cat_sig$gene <- rownames(statsFull_cat_sig)
statsFull_cat_sig <- filter(statsFull_cat_sig, !(gene %in% errgs_2cat2lin)) # remove the ones that don't converge for other models I will use
statsFull_cat_sig <- left_join(statsFull_cat_sig, subset(GO_annotation, select = c("query","Preferred_name","Description")), by = join_by(gene == query))
rownames(statsFull_cat_sig) <- statsFull_cat_sig$gene

#### refit genes and do post hoc analyses -------------------------------------------------------------------------

### model with interaction
head(statscat_withint_sig)
statscat_withint_sig <- filter(statscat_withint_sig, !(gene %in% errgs_2cat2lin))
# for each gene table with effect size in 1 week and 8 week X hvc and hvb
## gene # age # hvb  # cvb #
## g1   # 1w  #  0.2 # 0.3 #
## g1   # 8w  #  0.3 # ns  #
## g2   # 1w  # -0.1 # 0.1 #
## g2   # 8w  #  ns  # ns  #

# here I need post-hoc
genelist <- as.list(rownames(statscat_withint_sig))
names(genelist) <- rownames(statscat_withint_sig)
GeneFits_catint <- lapply(genelist, function(x) suppressMessages(glmmRefit(fit_cat_par, gene = x)))  
genefit <- GeneFits_catint[[3]] # inspect
summary(genefit)
# here I want the estimates COMPARED TO BENIGN for hot in age 1, hot in age 8 etc
EMM <- emmeans(genefit, ~ treatment * age)
summary(EMM)
empairs <- pairs(EMM, simple = "treatment")    # compare treats for each age, this is exactly what I need
summary(empairs)
pairs(EMM, simple = "age")     # compare ages for each treat
pairs(EMM)                   # compares all combinations (also ones I am not interested in)
# do emmeans for all genes
EMres_cat <- lapply(GeneFits_catint, function(fit) emmeans(fit, ~ treatment * age))
EMpairs_cat <- lapply(EMres_cat, function(emm) pairs(emm, simple = "treatment"))
w1_est <- sapply(EMpairs_cat, function(emm) summary(emm)$estimate[c(1,2)])
w8_est <- sapply(EMpairs_cat, function(emm) summary(emm)$estimate[c(4,5)])
w1_est <- t(w1_est)
w8_est <- t(w8_est)
w1_pval <- sapply(EMpairs_cat, function(emm) summary(emm)$p.value[c(1,2)])
w8_pval <- sapply(EMpairs_cat, function(emm) summary(emm)$p.value[c(4,5)])
w1_pval <- t(w1_pval)
w8_pval <- t(w8_pval)
cat_effects <- data.frame(gene = rownames(w1_est), age = rep("week1",nrow(w1_est)),cvb = w1_est[,1], hvb = w1_est[,2],
                          cvb_p = w1_pval[,1],hvb_p = w1_pval[,2], row.names = NULL)
cat_effects8 <- data.frame(gene = rownames(w8_est),age = rep("week8",nrow(w8_est)), cvb = w8_est[,1], hvb = w8_est[,2],
                           cvb_p = w8_pval[,1],hvb_p = w8_pval[,2], row.names = NULL)
cat_effects <- rbind(cat_effects,cat_effects8)
summary(cat_effects)
# make one with all values wide
cat_eff_wide_all <- pivot_wider(select(cat_effects, gene:hvb), names_from = age, values_from = c(cvb,hvb))
# make non significant effects NA
cat_effects[which(cat_effects$cvb_p > 0.01),"cvb"] <- NA
cat_effects[which(cat_effects$hvb_p > 0.01),"hvb"] <- NA
nocateff <-  filter(cat_effects, is.na(cvb) & is.na(hvb))$gene
cat_effects <- filter(cat_effects, !(is.na(cvb) & is.na(hvb))) # remove rows with NA in both effects, removes a lot actually I guess thats hvc
cat_eff_wide_all <- filter(cat_eff_wide_all, !(gene %in% nocateff))
summary(cat_effects)
table(cat_effects$age, is.na(cat_effects$cvb)) # 135 for 1w and 99 for 8w
table(cat_effects$age, is.na(cat_effects$hvb)) # 201 for 1w and 121 for 8w
# how many hvb are significant in both ages?
table(duplicated(filter(cat_effects, !is.na(hvb))$gene)) # 109
# how many cvb are significant in both ages?
table(duplicated(filter(cat_effects, !is.na(cvb))$gene)) # 73

# check if ages are sig different in both heat and cold
test <- pairs(EMM, simple = "age") 
summary(test)[3,6]
EMpairs_cat_age <- lapply(EMres_cat, function(emm) pairs(emm, simple = "age"))
int_pval_cold <- sapply(EMpairs_cat_age, function(emm) summary(emm)[2,7])
int_pval_hot <- sapply(EMpairs_cat_age, function(emm) summary(emm)[3,7])
summary(int_pval_cold) # not all sig
summary(int_pval_hot) # not all sig
int_pval_cold <- int_pval_cold[which(int_pval_cold <= 0.01)]
int_pval_hot <- int_pval_hot[which(int_pval_hot <= 0.01)]

# save(GeneFits_catint, cat_effects, file = "PostHoc_catInt.RData")


## how many are there?
length(unique(c(rownames(statsct_agesig), rownames(statscat_withint_sig)))) # 547



### model without interaction:
head(statsFull_cat_sig)
statsFull_cat_sig <- filter(statsFull_cat_sig, !(gene %in% errgs_2cat2lin))
genelist <- as.list(rownames(statsFull_cat_sig))
names(genelist) <- rownames(statsFull_cat_sig)
GeneFits_catNoint <- lapply(genelist, function(x) suppressMessages(glmmRefit(fit_cat_noint_par, gene = x))) # takes about a minute
genefit <- GeneFits_catNoint[[1]] # inspect
summary(genefit)
# here I want the estimates COMPARED TO BENIGN for hot and cold (no matter age)
EMres_trt <- lapply(GeneFits_catNoint, function(fit) emmeans(fit, ~ treatment))
EMpairs_trt <- lapply(EMres_trt, function(emm) pairs(emm))
summary(EMpairs_trt[[1]])
cvb_est_noint <- sapply(EMpairs_trt, function(emm) summary(emm)$estimate[1])
hvb_est_noint <- sapply(EMpairs_trt, function(emm) summary(emm)$estimate[2])
cvb_p_noint <- sapply(EMpairs_trt, function(emm) summary(emm)$p.value[1])
hvb_p_noint <- sapply(EMpairs_trt, function(emm) summary(emm)$p.value[2])
gene <- names(EMpairs_trt)
cat_effects_noint <- data.frame(gene=gene,cvb=cvb_est_noint,hvb=hvb_est_noint, p_cvb=cvb_p_noint,p_hvb=hvb_p_noint)
# make df with all values
cat_effects_both_all <- full_join(cat_eff_wide_all,select(cat_effects_noint, gene:hvb))
# replace non-sig with NA
cat_effects_noint[which(cat_effects_noint$p_cvb>0.01),"cvb"] <- NA
cat_effects_noint[which(cat_effects_noint$p_hvb>0.01),"hvb"] <- NA
nocateff <- filter(cat_effects_noint, is.na(cvb) & is.na(hvb))$gene # genes that have both NA so no effect after all
cat_effects_noint <- filter(cat_effects_noint, !(is.na(cvb) & is.na(hvb))) # remove rows where both NA
cat_effects_both_all <- filter(cat_effects_both_all, !(gene %in% nocateff))
summary(cat_effects_noint) # 434 sig for hvb and 109 sig for cvb


### now also make one for all genes in the analysis
genelist <- as.list(unique(c(cat_effects$gene, cat_effects_noint$gene)))
names(genelist) <- unique(c(cat_effects$gene, cat_effects_noint$gene))
GeneFits_all <- lapply(genelist, function(x) suppressMessages(glmmRefit(fit_cat_noint_par, gene = x))) # takes a while
lingene <- GeneFits_all$g10974
emlingene <- emmeans(lingene, ~treatment)
emplingene <- pairs(emlingene)
summary(emplingene)
EMres_all <- lapply(GeneFits_all, function(fit) emmeans(fit, ~ treatment))
EMpairs_all <- lapply(EMres_all, function(emm) pairs(emm))
summary(EMpairs_all[[1]])
cvb_est_noint <- sapply(EMpairs_all, function(emm) summary(emm)$estimate[1])
hvb_est_noint <- sapply(EMpairs_all, function(emm) summary(emm)$estimate[2])
cvb_p_noint <- sapply(EMpairs_all, function(emm) summary(emm)$p.value[1])
hvb_p_noint <- sapply(EMpairs_all, function(emm) summary(emm)$p.value[2])
gene <- names(EMpairs_all)
cat_effects_all <- data.frame(gene=gene,cvb=cvb_est_noint,hvb=hvb_est_noint, p_cvb=cvb_p_noint,p_hvb=hvb_p_noint)


# save(cat_effects_noint, linear_effects_noint, GeneFits_catNoint, GeneFitsLinNoint,GeneFits_all, cat_effects_all, file = "PostHoc_NoInt.RData")



### sanity check: are fits for interaction and non-interaction models similar?
table(cat_effects$gene %in% cat_effects_noint$gene)
table(cat_effects_noint$gene %in% cat_effects$gene) # do not overlap a lot

plot(rowMeans(cbind(cat_effects_both_all$cvb_week8,cat_effects_both_all$cvb_week1)),cat_effects_both_all$cvb)
abline(0,1)
plot(rowMeans(cbind(cat_effects_both_all$hvb_week8,cat_effects_both_all$hvb_week1)),cat_effects_both_all$hvb)
abline(0,1)
# yes they are similar


#### get numbers and overlaps -----------------------------------------------------------------

# load results
load("PostHoc_catInt.RData")
load("PostHoc_NoInt.RData")

# get genes
genes_hvb1w <- filter(cat_effects, age == "week1" & !is.na(hvb))$gene
genes_hvb8w <- filter(cat_effects, age == "week8" & !is.na(hvb))$gene
genes_hvbNoInt <- filter(cat_effects_noint, !is.na(hvb))$gene
genes_cvb1w <- filter(cat_effects, age == "week1" & !is.na(cvb))$gene
genes_cvb8w <- filter(cat_effects, age == "week8" & !is.na(cvb))$gene
genes_cvbNoInt <- filter(cat_effects_noint, !is.na(cvb))$gene

tempgenes3 <- unique(c(cat_effects$gene, cat_effects_noint$gene))


# categorise hot, cold and shared
make_colvec <- function(x){
  names(x) <- x
  if (x %in% genes_hvbNoInt & x %in% genes_cvbNoInt){
    return("hotandcold")
  }else if (x %in% genes_hvbNoInt){
    return("hot")
  }else if (x %in% genes_cvbNoInt){
    return("cold")
  }else{
    return("out")
  }
}
colvec <- sapply(cat_effects_noint$gene, function(x) make_colvec(x)) # changed here from cat_effects_all to cat_effects_noint because we're not using linear anymore
table(colvec)
cat_effects_noint$colour <- colvec
# write.csv(cat_effects_noint, file = "cat_effects_noint.csv") # save


## up or downregulated
hotdown <- filter(cat_effects_noint, hvb > 0 & colour == "hot") 
hotup <- filter(cat_effects_noint, hvb < 0 & colour == "hot") 
colddown <- filter(cat_effects_noint, cvb > 0 & colour == "cold") 
coldup <- filter(cat_effects_noint, cvb < 0 & colour == "cold")

### Are the ones shared between hot and cold, same or opposite direction? (using also non-significant opposite genes)
hotdowncoldup <- filter(cat_effects_noint, hvb > 0 & cvb < 0 & colour == "hotandcold") 
hotupcolddown <- filter(cat_effects_noint, hvb < 0 & cvb > 0 & colour == "hotandcold") 
bothdown <- filter(cat_effects_noint, hvb > 0 & cvb > 0 & colour == "hotandcold") 
bothup <- filter(cat_effects_noint, hvb < 0 & cvb < 0 & colour == "hotandcold")


# and are same direction ones the ones that are also in linear?
table(c(bothdown$gene,bothup$gene) %in% unique(c(linear_effects$gene, linear_effects_noint$gene)))
table(c(hotupcolddown$gene,hotdowncoldup$gene) %in% unique(c(linear_effects$gene, linear_effects_noint$gene)))
# more than the smiley ones but still less than half...
# to count as shared they still have to be significant in both heat and cold
gplots::venn(list("cat opposite" = c(hotupcolddown$gene,hotdowncoldup$gene),
                  "hvb" = filter(cat_effects, !is.na(hvb))$gene,
                  "cvb" = filter(cat_effects, !is.na(cvb))$gene))


## compare intensity of age responses
# how many each age?
length(unique(c(genes_hvb1w,genes_cvb1w)))
length(unique(c(genes_hvb8w,genes_cvb8w)))
# effect intensity?
hist(filter(cat_effects, age == "week1")$hvb)
hist(filter(cat_effects, age == "week8")$hvb)
hist(filter(cat_effects, age == "week1")$cvb)
hist(filter(cat_effects, age == "week8")$cvb)
# higher in 1w a bit


### make wide data frame
# need to run PH_ageefftable script first to get the dfs
make_colvec2 <- function(x){
  names(x) <- x
  if (x %in% c(genes_hvbNoInt,genes_hvb1w,genes_hvb8w) & x %in% c(genes_cvbNoInt,genes_cvb1w,genes_cvb8w)){
    return("hotandcold")
  }else if (x %in% c(genes_hvbNoInt,genes_hvb1w,genes_hvb8w)){
    return("hot")
  }else if (x %in% c(genes_cvbNoInt,genes_cvb1w,genes_cvb8w)){
    return("cold")
  }else{
    return("wtf")
  }
}
colvecall <- sapply(cat_effects_all$gene, function(x) make_colvec2(x))
table(colvecall)
cat_effects_all$colour <- colvecall

cat_eff_wide <- pivot_wider(dplyr::select(cat_effects, gene:hvb), names_from = age, values_from = c(cvb,hvb))
cat_eff_wide[is.na(cat_eff_wide$cvb_week1),"cvb_week1"] <- 0
cat_eff_wide[is.na(cat_eff_wide$cvb_week8),"cvb_week8"] <- 0
cat_eff_wide[is.na(cat_eff_wide$hvb_week1),"hvb_week1"] <- 0
cat_eff_wide[is.na(cat_eff_wide$hvb_week8),"hvb_week8"] <- 0
# add info where expressed

cat_eff_wide <- left_join(cat_eff_wide, subset(cat_effects_all, select = c(gene, colour)))

# save(hotup,hotdown,coldup,colddown,hotdowncoldup,hotupcolddown,bothup,bothdown,errgs_2cat2lin, file = "NewPHgenesForGO.RData") # should only be no interaction genes!

# also make list of genes for Fst analysis
genes_of_interest_noint <- cat_effects_noint$gene
table(duplicated(genes_of_interest_noint)) # good no overlaps
# writeLines(genes_of_interest_noint, "../genes_of_interest_noint.txt")

#### venn plots -------------------------

gplots::venn(list("hvb" = filter(cat_effects, !is.na(hvb))$gene,
                  "cvb" = filter(cat_effects, !is.na(cvb))$gene))
gplots::venn(list("cvb1w" = filter(cat_effects, age == "week1" & !is.na(cvb))$gene,
                  "cvb8w" = filter(cat_effects, age == "week8" & !is.na(cvb))$gene,
                  "hvb1w" = filter(cat_effects, age == "week1" & !is.na(hvb))$gene,
                  "hvb8w" = filter(cat_effects, age == "week8" & !is.na(hvb))$gene))
gplots::venn(list("cat_int" = cat_effects$gene,
                  "cat_noint" = cat_effects_noint$gene))

# quite different
# shared and distinct 1w
gplots::venn(list("hvb1w" = filter(cat_effects, age == "week1" & !is.na(hvb))$gene,
                  "cvb1w" = filter(cat_effects, age == "week1" & !is.na(cvb))$gene))
# shared and distinct 8
gplots::venn(list("hvb8w" = filter(cat_effects, age == "week8" & !is.na(hvb))$gene,
                  "cvb8w" = filter(cat_effects, age == "week8" & !is.na(cvb))$gene))
#### plot results --------------------------

### plot showing amount of shared and distinct response
## new version of this with both models where distinct = genes ever only in heat or ever only in cold across ages & models
# and shared = the rest (if it's ever DE in both across ages and models)
# run PH_ageefftable.R first

make_distshared_new <- function(row_gene){
  hvb <- unname(unlist(row_gene["hvb"]))
  cvb <- unname(unlist(row_gene["cvb"]))
  if(hvb == "not_sig" | cvb == "not_sig"){
    return("distinct")
  }else{
    return("shared")
  }
}
distshared <- apply(age_dirs,1,make_distshared_new)
table(distshared) # seems plausible
age_dirs$distshared <- distshared
table(age_dirs$distshared,age_dirs$hvb) # seems plausible

# reshape so we count hot and cold independently of age and treatment (same gene can appear in heat and in cold) 
age_dirs_long <- pivot_longer(age_dirs, cols = cvb:hvb, names_to = "cond", values_to = "age_eff")
age_dirs_long <- filter(age_dirs_long, age_eff != "not_sig") # delete the not significant rows
# make labels
nr_cold <- table(age_dirs_long$cond)[1]
nr_hot <- table(age_dirs_long$cond)[2]
age_dirs_long$label <- NA
age_dirs_long[1,"label"] <- nr_cold
age_dirs_long[2,"label"] <- nr_hot
# flip distinct and shared and gray shared out because we are talking about distinct here
age_dirs_long$distshared <- factor(age_dirs_long$distshared, levels = c("shared","distinct"))
# plot
ggplot(age_dirs_long, aes(x = cond, fill = distshared,alpha = distshared))+
  geom_bar(position = "fill")+
  geom_text(aes(label = label, y = 1.05), alpha = 1)+
  scale_fill_manual(values = c("#00a08a","#f98400"), labels = c("Distinct","Shared"))+
  scale_alpha_manual(values = c(0.4,1))+
  scale_x_discrete(labels = c("cvb" = "Cold", "hvb" = "Heat"))+
  guides(fill = guide_legend(reverse = TRUE, override.aes = list(alpha = c(1,0.6))), alpha = "none")+
  xlab("")+ylab("Proportion of genes")+
  theme_pubr(base_size = 20)+
  theme(legend.title = element_blank())
# voilà

# now from the distinct, how many in cold and how many in heat
table(filter(age_dirs_long, distshared == "distinct")$cond)
23/386*100
363/386*100
# and how many of all cold genes are distinct, and how many of all heat genes are distinct?
table(age_dirs_long$cond,age_dirs_long$distshared)
23/224*100
363/564*100
# how many of all genes are consistently involved in the response to both heat and cold across development?
table(filter(age_dirs_long, distshared == "shared" & age_eff == "same_slope")$cond)


## bar distinct with up and down
# make df with effect sizes for distinct genes from both models (then I don't have to deal with ages here at all)
distgenes <- filter(age_dirs, distshared == "distinct")$gene
cat_eff_dist <- bind_rows(subset(cat_effects, gene %in% distgenes, select = c(gene,cvb,hvb)),
                      filter(cat_effects_noint, gene %in% distgenes & !(gene %in% cat_effects$gene)))
rownames(cat_eff_dist) <- NULL
# in there are some genes that are upregulated in one age but downregulated in another, thus some genes may appear twice in the bar plot
# but will make genes that are expressed in different intensities across ages but same direction just once
same_dir_hot <- filter(age_dirs, distshared == "distinct" & hvb == "same_direction")$gene
same_dir_cold <- filter(age_dirs, distshared == "distinct" & cvb == "same_direction")$gene # there are none
which(duplicated(cat_eff_dist$gene)) # 19
cat_eff_dist <- cat_eff_dist[-intersect(which(cat_eff_dist$gene %in% same_dir_hot),which(duplicated(cat_eff_dist$gene))),] # remove double genes with same direction but only one! should remove 16
rownames(cat_eff_dist) <- NULL
which(duplicated(cat_eff_dist$gene)) # so note in legend that 3 genes are in there double (one cold, two hot)
cat_eff_dist$colour <- ifelse(is.na(cat_eff_dist$hvb), "cold","hot")
make_dir_dist <- function(row){
  colour <- unname(unlist(row["colour"]))
  cvb <- as.numeric(row["cvb"])
  hvb <- as.numeric(row["hvb"])
  if(!(is.na(hvb)) & colour == "hot" & (hvb > 0)){
    return("down_regulated")
  }else if(!(is.na(cvb)) & colour == "cold" & (cvb > 0)){
    return("down_regulated")
  }else{return("up_regulated")}
}
dist_dir <- apply(cat_eff_dist, 1, make_dir_dist)
cat_eff_dist$dist_dir <- dist_dir
table(cat_eff_dist$dist_dir, cat_eff_dist$colour)




### plot with hot dist and shared and cold dist and shared showing behaviour in ages
# run PH_ageefftable first
colvec2 <- sapply(cat_eff2$gene, function(x) make_colvec2(x))
table(colvec2) # idk could be right I guess?
cat_eff2$colour <- colvec2
# add response type per condition
hotcoldgenes_age <- unname(unlist(cat_eff2[which(duplicated(cat_eff2$gene)),"gene"]))
hotcoldgenes_age %in% filter(cat_eff2, colour == "hotandcold")$gene # good
filter(cat_eff2, colour == "hotandcold")$gene %in% hotcoldgenes_age # two false, g18643 and g17406, both just1week
# some kind of mistake just change it
cat_eff2[which(cat_eff2$gene %in% c("g18643","g17406") & cat_eff2$condition == "hvb"), "colour"] <- "hot"
# now add response type but check numbers first so no overlap!
nshared <- nrow(cat_eff2[which(cat_eff2$gene %in% hotcoldgenes_age & cat_eff2$condition == "hvb"),])
nshared <- nrow(cat_eff2[which(cat_eff2$gene %in% hotcoldgenes_age & cat_eff2$condition == "cvb"),])
nhot <- nrow(cat_eff2[which(!(cat_eff2$gene %in% hotcoldgenes_age) & cat_eff2$condition == "hvb"),])
ncold <- nrow(cat_eff2[which(!(cat_eff2$gene %in% hotcoldgenes_age) & cat_eff2$condition == "cvb"),])
nshared + nhot + ncold # correct!

cat_eff2[which(cat_eff2$gene %in% hotcoldgenes_age & cat_eff2$condition == "hvb"),"condresp"] <- "hot_shared"
cat_eff2[which(cat_eff2$gene %in% hotcoldgenes_age & cat_eff2$condition == "cvb"),"condresp"] <- "cold_shared"
cat_eff2[which(!(cat_eff2$gene %in% hotcoldgenes_age) & cat_eff2$condition == "hvb"),"condresp"] <- "hot_distinct"
cat_eff2[which(!(cat_eff2$gene %in% hotcoldgenes_age) & cat_eff2$condition == "cvb"),"condresp"] <- "cold_distinct"
# add column for label (adjust numbers manually if needed)
labsdf <- data.frame(condresp = c("cold_distinct","hot_distinct","cold_shared","hot_shared"), label = c(paste0("n = ",ncold), paste0("n = ",nhot),
                                                                                                        paste0("n = ", nshared),paste0("n = ", nshared)))
cat_eff2 <- left_join(cat_eff2, labsdf)
# make nicer order
cat_eff2$ageeff <- factor(cat_eff2$age_eff, c("same_slope","same_direction","opposite_direction","just_1week","just_8week"))
cat_eff2$condresp <- factor(cat_eff2$condresp, c("cold_distinct","hot_distinct","cold_shared","hot_shared"))
# also save it for targeted fst analysis
# save(cat_eff2, file = "cat_eff_age_resp.RData")


### new bar plot with modulated by age yes/no and colour shared/distinct
# each gene only once
age_dirs$age_mod <- ifelse((age_dirs$hvb == "same_slope" | age_dirs$cvb == "same_slope"), "identical_response", "age_modulated")
table(age_dirs$age_mod) # adds up now, 458-96 = 362!
diffgenes <- cat_effects_noint[which(!(cat_effects_noint$gene %in% filter(age_dirs, age_mod == "identical_response")$gene)), "gene"] # I mean it is different models, I guess they just have slightly different results
# I am guessing they are mostly just different slope but let's check
diffgenes <- filter(age_dirs, gene %in% diffgenes) # no it's all types of things
table(diffgenes$cvb, diffgenes$hvb)
# anyway, add response type hotandcold if there's anything going on in both temps in any age!
cat_eff2$resptype <- ifelse(cat_eff2$colour == "hotandcold","shared","distinct")
# test that the same gene always has the same resptype, thus, distinct genes should only occur once!
unname(unlist(cat_eff2[which(duplicated(cat_eff2$gene)),"gene"])) %in% filter(cat_eff2, resptype == "shared")$gene # good
table(filter(cat_eff2, resptype == "shared")$gene %in% unname(unlist(cat_eff2[which(duplicated(cat_eff2$gene)),"gene"]))) # good!
table(filter(cat_eff2, resptype == "distinct")$gene %in% unname(unlist(cat_eff2[which(duplicated(cat_eff2$gene)),"gene"]))) # good!

age_dirs2 <- inner_join(age_dirs, subset(cat_eff2, select = c("gene","resptype")))
age_dirs2 <- distinct(age_dirs2) # it duplicates the rows, so remove again
table(age_dirs2$age_mod) #
table(age_dirs2$resptype)
# plot
ggplot(age_dirs2, aes(x = age_mod, fill = resptype))+
  geom_bar(position = position_dodge())+
  scale_fill_manual(values = c("#f98400", "#00a08a"))+
  ylab("Number of genes")+ xlab("")+
  theme_pubr(base_size = 20)+
  theme(legend.title = element_blank())
# plot other way around
ggplot(age_dirs2, aes(x = resptype, fill = age_mod))+
  geom_bar(position = position_dodge())+
  scale_fill_manual(values = c("#f98400", "#00a08a"), labels = c("Age modulated","Not age modulated"))+
  scale_x_discrete(labels = c("Distinct","Shared"))+
  ylab("Number of genes")+ xlab("")+
  theme_pubr(base_size = 20)+
  theme(legend.title = element_blank())+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))


# NEW: hot, cold and hotandcold separately
# make facet labeller
colour_labs <- c("Distinct: hot", "Distinct: cold", "Shared")
names(colour_labs) <- c("hot","cold","hotandcold")
ggplot(filter(cat_eff2, condition == "cvb" & ageeff != "same_slope"), aes (x = ageeff, fill = colour))+
  facet_wrap(facets = vars(colour), ncol = 2, labeller = labeller(colour = colour_labs))+
  geom_bar()+
  xlab("")+ylab("Number of genes")+
  theme_pubr(base_size = 20)+
  scale_fill_manual(values = unname(tempcols[c(2,1)]), labels = c("Distinct: cold","Shared"))+
  scale_x_discrete(labels = c("Same direction","Opposite direction","Only 1-week-old", "Only 8-week-old"))+
  ylim(c(0,77))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 18), legend.position = "inside", legend.position.inside = c(0.2,0.9), legend.title = element_blank())
ggplot(filter(cat_eff2, condition == "hvb" & ageeff != "same_slope"), aes (x = ageeff, fill = colour))+
  facet_wrap(facets = vars(colour), ncol = 2, labeller = labeller(colour = colour_labs))+
  geom_bar()+
  xlab("")+ylab("Number of genes")+
  scale_fill_manual(values = unname(tempcols[c(3,1)]), labels = c("Distinct: hot","Shared"))+
  scale_x_discrete(labels = c("Same direction","Opposite direction","Only 1-week-old", "Only 8-week-old"))+
  ylim(c(0,77))+
  theme_pubr(base_size = 20)+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 18), legend.position = "inside", legend.position.inside = c(0.2,0.9), legend.title = element_blank())

# save this info for targeted fst analysis



### dot plots
# turn signs around so it fits with the text
cat_eff_wide$cvb_week1 <- cat_eff_wide$cvb_week1*-1
cat_eff_wide$cvb_week8 <- cat_eff_wide$cvb_week8*-1
cat_eff_wide$hvb_week1 <- cat_eff_wide$hvb_week1*-1
cat_eff_wide$hvb_week8 <- cat_eff_wide$hvb_week8*-1


ggplot(filter(cat_eff_wide, !(hvb_week1 == 0 & hvb_week8 == 0)), aes(x = hvb_week1, y = hvb_week8, colour = colour))+
  geom_jitter(height=0.02, width = 0.02, alpha = 0.6)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  scale_colour_manual(values = unname(tempcols[c(3,1)]), labels = c("Distinct: hot", "Shared"))+
  scale_x_continuous(name="Expression in heat compared to benign \n (1-week-olds)", breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5), limits = c(-1.85,1.85))+
  scale_y_continuous(name="Expression in heat compared to benign \n (8-week-olds)", breaks = c(-1,-0.5,0,0.5,1), limits = c(-1.2,1.2))+
  theme_minimal(base_size = 20)+
  theme(legend.title = element_blank(), legend.position = "top")

ggplot(filter(cat_eff_wide, !(cvb_week1 == 0 & cvb_week8 == 0)), aes(x = cvb_week1, y = cvb_week8, colour = colour))+
  geom_jitter(height=0.02, width = 0.02, alpha = 0.6)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  scale_colour_manual(values = unname(tempcols[c(2,1)]), labels = c("Distinct: cold", "Shared"))+
  scale_x_continuous(name="Expression in cold compared to benign \n (1-week-olds)", breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5), limits = c(-1.85,1.85))+
  scale_y_continuous(name="Expression in cold compared to benign \n (8-week-olds)", breaks = c(-1,-0.5,0,0.5,1), limits = c(-1.2,1.2))+
  theme_minimal(base_size = 20)+
  theme(legend.title = element_blank(), legend.position = "top") # quite a lot different directions!


