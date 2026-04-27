# age effect combinations in genes that change expression in both heat and cold
# Nov 2025
# Lea Huber


library(tidyverse)
library(ggpubr)

setwd("/Users/Lea/Library/CloudStorage/OneDrive-Aarhusuniversitet/LeaH/ChamberTempExp/data/working/RNASeq_countdata")
# load results
load("PostHoc_catInt.RData")
load("PostHoc_NoInt.RData")



# turn cat effects around
cat_eff_long <- pivot_longer(cat_effects, cols = cvb:hvb, names_to = "condition", values_to = "effect")
cat_eff2 <- pivot_wider(dplyr::select(cat_eff_long, c(gene,age,condition,effect)), names_from = age, values_from = effect)

whichageeff3 <- function(row){
  week1 <- as.numeric(row["week1"])
  week8 <- as.numeric(row["week8"])
  if(is.na(week1) & is.na(week8)){
    return("no_effect")
  } else if(!is.na(week1) & is.na(week8)){
    return("just_1week")
  } else if(is.na(week1) & !is.na(week8)){
    return("just_8week")
  } else if(week1 * week8 > 0){
    return("same_direction")
  } else if(week1 * week8 < 0){
    return("opposite_direction")
  } else {
    return("wtf")
  }
}

age_eff <- apply(cat_eff2, 1, whichageeff3)
cat_eff2$age_eff <- age_eff
cat_eff2 <- filter(cat_eff2, age_eff != "no_effect")

cat_eff_noint_dummy <- pivot_longer(dplyr::select(cat_effects_noint,gene:hvb), cols = cvb:hvb, names_to = "condition", values_to = "week1")
cat_eff_noint_dummy$week8 <- cat_eff_noint_dummy$week1
cat_eff_noint_dummy$age_eff <- rep("same_slope",nrow(cat_eff_noint_dummy))
cat_eff_noint_dummy <- filter(cat_eff_noint_dummy, !(is.na(week1) & is.na(week8))) # remove rows of genes with no effect in either cvb or hvb
table(unique(cat_eff_noint_dummy$gene) %in% unique(cat_eff2$gene)) # 96 no age effect genes are are also in age interaction model
table(unique(cat_eff2$gene) %in% unique(cat_eff_noint_dummy$gene)) # 129 age int genes are not in no interaction, but 96 are, there the same slope will potentially be overwritten!
cat_eff2 <- rbind(cat_eff2, filter(cat_eff_noint_dummy, !(gene %in% cat_eff2$gene))) # so this will add 362 genes with same_slope, 225 + 362 = 587
length(unique(cat_eff2$gene)) # yep
table(cat_eff2$age_eff)

# reshape again
age_dirs <- pivot_wider(dplyr::select(cat_eff2, c(gene,condition,age_eff)), names_from = condition, values_from = age_eff)
table(duplicated(age_dirs$gene)) # fuck yeah!
# the ones with NA are just not significant in that condition
age_dirs[is.na(age_dirs$cvb),"cvb"] <- "not_sig"
age_dirs[is.na(age_dirs$hvb),"hvb"] <- "not_sig"

table(age_dirs$cvb,age_dirs$hvb)
# tadaa

# make it nice
library(knitr)
library(kableExtra)

tbl <- table(age_dirs$cvb,age_dirs$hvb)
tbl2 <- table(filter(age_dirs, !(cvb %in% c("not_sig","same_slope") | hvb %in% c("not_sig","same_slope")))$cvb,
             filter(age_dirs, !(cvb %in% c("not_sig","same_slope") | hvb %in% c("not_sig","same_slope")))$hvb)
# this removes quite a few genes also from other categories 

#kbl(tbl, caption = "In heat") %>%
  kable_classic(full_width = FALSE) %>%
  pack_rows("In cold",1,6) %>%
  save_kable("../../../figures/glmmSeq_PostHoc/table_ageeffs.pdf")

#kbl(tbl2, caption = "In heat") %>%
  kable_classic(full_width = FALSE) %>%
  pack_rows("In cold",1,4)%>%
  save_kable("../../../figures/glmmSeq_PostHoc/table_ageeffs_nosigsame.pdf")

## also save as txt for targeted fst analysis
#writeLines(age_dirs$gene, "../genesOI_ageint.txt")
# and just as R data
#save(age_dirs, file = "age_dir_cats.RData")

