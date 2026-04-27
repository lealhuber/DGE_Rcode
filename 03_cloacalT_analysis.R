## chamber experiment chick cloacal temperature ##

#### load packages ####
library("readxl")
library("tidyverse")
library("ggpubr")
library("RColorBrewer")
library("ggExtra")
library(lme4)
library(viridis)

#### preliminaries ####
setwd("/Users/Lea/Library/CloudStorage/OneDrive-Aarhusuniversitet/LeaH/ChamberTempExp/data/working")
figurepath="/Users/Lea/Library/CloudStorage/OneDrive-Aarhusuniversitet/LeaH/ChamberTempExp/figures/Physiology"
### colours
tempcols = c(benign = "#6800a4", cold = "#377EB8", hot = "#E41A1C") # or should I do purple for benign? #4DAF4A
agecols = c("#1B9E77" ,"#D95F02")
sexcols = c("#7570B3", "#E7298A")
shadecols = c("orange","gray")
shadecols2 = c("black","gray")

# in here there's everything
load("ChamberTempSI_new.RData")
# re-level treatment so benign in the middle
sampleInfo$treatment <- factor(sampleInfo$treatment, levels = c("cold","benign","hot")) # maybe this messes up some things watch out!

# one row per chick
AllChicks <- pivot_wider(datTall, id_cols = !c(Date, SampleNr, Samplingtime, chamberT), names_from = treatment, values_from = cloacat) # should have 34 rows
SIwide <- pivot_wider(sampleInfo, id_cols = !c(Date, SampleNr, Samplingtime, chamberT), names_from = treatment, values_from = cloacat) # should have 20 rows

# often I called it datRNAseq, just the same as sampleInfo 
datRNAseq <- sampleInfo


#### plotting cloacal temp data ####

ggplot(datTalive, aes(x=treatment, y=cloacat, fill= sex))+
  geom_boxplot()+
  geom_jitter(position = position_jitterdodge(jitter.height = 0), aes(shape = Age))+
  facet_wrap(vars(as.factor(Date)))+
  scale_fill_manual(values=sexcols)+
  theme_pubr()
  
#per group: group 1 is colder in benign than cold (that was day 2), no consistent age difference
#per date: first day best differenciated
# using only survivors cold variation on first day 1 week and benign variation on 2nd day 8 week smaller
# don't see any sex differences whatsoever

# include info wether chick died the next day
datTall$deadnextday =c(rep(NA))
for (i in 1:nrow(datTall)){
    if(is.na(datTall$mort_dat[i])){ 
    }else if(datTall$Date[i] == as.POSIXct("2023-10-25", tz = "UTC") & datTall$mort_dat[i] == as.POSIXct("2023-10-26", tz = "UTC") |
       (datTall$Date[i] == as.POSIXct("2023-10-26", tz = "UTC") & datTall$mort_dat[i] == as.POSIXct("2023-10-27", tz = "UTC"))){
      datTall$deadnextday[i] = "yes"
    }else{
      datTall$deadnextday[i] = "no_or_NA"
  }
}
datTall$deadnextday <- as.factor(datTall$deadnextday)
ggplot(datTall, aes(x=treatment, y=cloacat, fill= Age))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(position = position_jitterdodge(jitter.height = 0), aes(color = deadnextday))+
  scale_fill_manual(values=agecols)+
  scale_color_manual(values = c("black","red"))+
  theme_pubr()


### Plot with temperature during experiment each day and treatment

ggplot(sampleInfo, aes(x=chamberT, y=cloacat, color=treatment, shape = Age))+
  geom_jitter(height = 0.1, width = 0.1, size = 3)+
  scale_x_continuous(limits = c(10,47))+
  scale_color_manual(values= tempcols)+
  xlab("Ambient temperature [°C]")+
  ylab("Cloacal temperature [°C]")+
  theme_pubr(base_size = 25)
# ggsave("cloacaTvchamberT_120to-30_SIdata.pdf", path = figurepath)


# import temp data
hdrs = c("Reading","Date.Time","temp","humidity","dewpoint","TrueTime","ManTime")
tempdata = vector("list", length = 9)
names(tempdata) = c("hot1", "hot2", "hot3","cold1","cold2","cold3","benign1","benign2","benign3")
tempdata$hot1 = read_xlsx("Temperatures/Day1/1HotDay1.xlsx", skip = 11, col_names = hdrs)[,1:7]
tempdata$hot3 = read_xlsx("Temperatures/Day3/1HotDay3.xlsx", skip = 11, col_names = hdrs)[,1:7]
tempdata$hot2 = read_xlsx("Temperatures/Day2/1HotDay2.xlsx", skip = 11, col_names = hdrs)[,1:7]
tempdata$cold1 = read_xlsx("Temperatures/Day1/2ColdDay1.xlsx", skip = 11, col_names = hdrs)[,1:7]
tempdata$cold2 = read_xlsx("Temperatures/Day2/2ColdDay2.xlsx", skip = 11, col_names = hdrs)[,1:7]
tempdata$cold3 = read_xlsx("Temperatures/Day3/2ColdDay3.xlsx", skip = 11, col_names = hdrs)[,1:7]
tempdata$benign1 = read_xlsx("Temperatures/Day1/3BenignDay1.xlsx", skip = 11, col_names = hdrs)[,1:7]
tempdata$benign2 = read_xlsx("Temperatures/Day2/3BenignDay2.xlsx", skip = 11, col_names = hdrs)[,1:7]
tempdata$benign3 = read_xlsx("Temperatures/Day3/3BenignDay3.xlsx", skip = 11, col_names = hdrs)[,1:7]


## plot all temperatures on all experiment days together
ggplot()+
  geom_line(data=tempdata$hot1, aes(x=ManTime,y=temp, linetype = "solid"),color = tempcols[3])+
  geom_line(data=tempdata$cold1,aes(x=ManTime,y=temp, linetype = "solid"), color = tempcols[2])+
  geom_line(data=tempdata$benign1,aes(x=ManTime,y=temp, linetype = "solid"), color = tempcols[1])+
  geom_line(data=tempdata$hot2, aes(x=ManTime,y=temp, linetype = "dashed"), color = tempcols[3])+
  geom_line(data=tempdata$cold2,aes(x=ManTime,y=temp, linetype = "dashed"), color = tempcols[2])+
  geom_line(data=tempdata$benign2,aes(x=ManTime,y=temp, linetype = "dashed"), color = tempcols[1])+
  geom_line(data=tempdata$hot3, aes(x=ManTime,y=temp, linetype = "twodash"), color = tempcols[3])+
  geom_line(data=tempdata$cold3,aes(x=ManTime,y=temp, linetype = "twodash"), color = tempcols[2])+
  geom_line(data=tempdata$benign3,aes(x=ManTime,y=temp, linetype = "twodash"), color = tempcols[1])+
  scale_x_datetime(limits = c(tempdata$hot1$ManTime[1],tempdata$benign3$ManTime[431]), date_breaks = "1 hour",
                   date_labels = "%H:%M", name = "Time")+
  scale_y_continuous(limits = c(10,48))+
  scale_color_discrete(name="Chamber", labels = c("Hot","Benign","Cold"))+
  scale_linetype_discrete(name = "Day", labels = c("1","2","3"))+
  ylab("Air temperature [°C]")+
  theme_pubr(base_size = 20)
#ggsave("Chamber_temp_all.pdf", path = figurepath)


#### Cloaca data statistics ------------------------------------------------------------------------
### Only with sequenced chicks
# Question: Does treatment significantly influence cloaca temp, and does age significantly influence cloaca temp?
# response variable: Cloaca temp, assume normal distribution
# regression: cloacat ~ treatment*Age + (1|group) + (1|chickno)
# and then two-way repeated measures ANOVA -> is it different?
# and then post-hoc test -> where is the difference?
library(lme4)
library(GGally)
library(DHARMa)
library(MuMIn)
library(lmerTest)
library(rstatix)
library(emmeans)
library(visreg)

SIwide = pivot_wider(dplyr::select(sampleInfo,c(treatment,chickno, cloacat,group,Age,sex,Mass,massScaled)),
                     names_from = treatment, values_from = cloacat) # wide with one row per chick
SIwide = mutate(SIwide,TchangeHot = hot-benign)
SIwide = mutate(SIwide,TchangeCold = cold-benign)
SIwide <- mutate(SIwide, TchangeAv = (abs(TchangeHot)+abs(TchangeCold))/2)
SIwide$chickno <- as.factor(SIwide$chickno)
# add to sampleInfo
sampleInfo <- left_join(sampleInfo, select(SIwide, c(chickno, TchangeHot, TchangeCold, TchangeAv)))

datRNAseq <- sampleInfo

m1 <- lmer(cloacat ~ treatment*Age + treatment*sex + (1|group) + (1|chickno) + (1|Date), data = datRNAseq)
simulationOutput <- simulateResiduals(fittedModel = m1, n = 500) # simulate data from our model n times
# This method checks if your model is useful by seeing if it can even produce data that looks like the data you used to fit your model (i.e. your observations). 
plot(simulationOutput) # looks just ok

# Are residuals normally distributed at each time point?
ggqqplot(datRNAseq, "cloacat") # hm not sure
shapiro.test(filter(datRNAseq, as.factor(Date) == "2023-10-25")$cloacat) # p-value =  0.1156 -> normal
shapiro.test(filter(datRNAseq, as.factor(Date) == "2023-10-26")$cloacat) # p-value =  0.09247 -> normal
shapiro.test(filter(datRNAseq, as.factor(Date) == "2023-10-27")$cloacat) # p-value =  0.4842 -> normal


options(na.action = "na.fail") # needed for dredge() function to prevent illegal model comparisons
dredgeOut<-dredge(m1) # fit and compare a model set representing all possible predictor combinations
# gives singularity warning
dredgeOut # seems like only treatment has an effect, age is at the border (but w/o interaction)

### so lets do it with only treatment, check different random effects
m2 <- lmerTest::lmer(cloacat ~ treatment + (1|chickno) + (1|group)+ (1|Date), data = datRNAseq)
m3 <- lmerTest::lmer(cloacat ~ treatment + (1|chickno) + (0 + treatment|group)+ (0 + treatment|Date), data = datRNAseq)
m4 <- lmerTest::lmer(cloacat ~ treatment + (1|chickno) + (1 + treatment|group)+ (1 + treatment|Date), data = datRNAseq)
m5 <- lmerTest::lmer(cloacat ~ treatment + (1|chickno) + (1|Date), data = datRNAseq) # this ends up being the best one
simulationOutput <- simulateResiduals(fittedModel = m5, n = 500) # simulate data from our model n times
# This method checks if your model is useful by seeing if it can even produce data that looks like the data you used to fit your model (i.e. your observations). 
plot(simulationOutput) # looks ok 
AIC(m2,m3,m4,m5,m6) # m2 and m5 have almost the same AIC, lower than m3 and m4

ranova(m2)
ranova(m3)
ranova(m4)
ranova(m5)

step(m2)
get_model(step(m2)) # here it says m5 is the best


summ5 <- summary(m5) # here there are already p values saying all treatments are significantly different from each other
summ_co <- summ5$coefficients # extract model coefficients
# extract variance explained by random effects
v2 <- VarCorr(m5)
v2df <- as.data.frame(v2)
v2df <- subset(v2df, select = -var2)
v2df <- reshape::rename(v2df, c(grp = "Random effect", var1 = "Var1", vcov = "Variance", sdcor = "Std. Dev."))
v2df
# chickno explains more variance than Date but less than residual variance

summary(m5)
(aov_res5 <- anova(m5))

## now that we've determined that treatment has an effect, which treatments are different from which?
emmeans(m5, pairwise ~ treatment, type = "response") # this is what I wanted !!!!
plot(emmeans(m5, pairwise ~ treatment, type = "response"), comparisons = TRUE)+
  coord_flip()+
  theme_pubr()

# visualise
ggplot(sampleInfo, aes(x=treatment, y=cloacat, fill=Age))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(shape = sex, group = interaction(treatment, Age)), size = 2, 
              position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.75))+
  geom_signif(comparisons = list(c("hot","cold")),map_signif_level=TRUE)+
  geom_signif(comparisons = list(c("hot","benign")),map_signif_level=TRUE, y_position = c(42,1))+
  ylab("cloaca temperature [°C]")+
  scale_fill_manual(values=agecols)+
  theme_pubr(base_size = 22)
# ggsave("cloacatVtreatment_box_stats.pdf", path = figurepath)

# Plot without differentiating ages
ggplot(sampleInfo, aes(x=treatment, y=cloacat, colour=treatment))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter()+
  geom_signif(comparisons = list(c("hot","cold")),map_signif_level=TRUE,color = "black")+
  geom_signif(comparisons = list(c("hot","benign")),map_signif_level=TRUE, y_position = c(42,1),color = "black")+
  ylab("cloaca temperature [°C]")+xlab("")+
  scale_colour_manual(values=tempcols)+
  theme_pubr(base_size = 22)
#ggsave("cloacatVtreatment_box_stats2.pdf", path = figurepath)

## for presentation hot and cold graphs separately
ggplot(filter(sampleInfo, treatment != "hot"), aes(x=treatment, y=cloacat, fill=Age))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(shape = sex, group = interaction(treatment, Age)), size = 2, 
              position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.75))+
  #geom_signif(comparisons = list(c("hot","benign")),map_signif_level=TRUE, y_position = c(42,1))+
  ylab("cloaca temperature [°C]")+
  scale_fill_manual(values=agecols)+
  theme_pubr(base_size = 22)
# ggsave("cloacatVtreatment_box_justcold.pdf", path = figurepath)

# or even without age
ggplot(filter(sampleInfo, treatment != "hot"), aes(x=treatment, y=cloacat, fill = treatment))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(shape = sex), size = 3, 
              position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.75))+
  #geom_signif(comparisons = list(c("hot","benign")),map_signif_level=TRUE, y_position = c(42,1))+
  scale_fill_manual(values = tempcols)+
  ylab("cloaca temperature [°C]")+
  theme_pubr(base_size = 25)
ggsave("cloacatVtreatment_box_coldnoage.pdf", path = figurepath)




#### Temp deviation in heat vs cold
# filter and scale
SIwide <- filter(SIwide, !is.na(TchangeCold) & !is.na(TchangeHot))
SIwide$massScaled <- as.numeric(SIwide$massScaled)
# fit first with everything that could influence it
fit_bm_as <- lmer(TchangeCold ~ TchangeHot*massScaled + Age + sex + (1|group), data = SIwide)

# check normality of residuals
simulationOutput <- simulateResiduals(fittedModel = fit_bm_as, n = 500) # simulate data from our model n times
plot(simulationOutput) # looks ok

ranova(fit_bm_as)
step(fit_bm_as)
# group as random is apparently not necessary

# and then look whether reduced is better
options(na.action = "na.fail") # needed for dredge() function to prevent illegal model comparisons
dredgeOut<-dredge(fit_bm_as) # fit and compare a model set representing all possible predictor combinations
dredgeOut
# tchange definitely in, plus mass basically the same but age is also close


# based on this new model
fit_bm <- lm(TchangeCold ~ TchangeHot + massScaled, data = SIwide) # can only be included as fixed effect obv
simulationOutput <- simulateResiduals(fittedModel = fit_bm, n = 500) # simulate data from our model n times
plot(simulationOutput, asFactor=FALSE) # looking OK


summary(fit_bm) # tchange hot is significan but body mass not significant

# visualise
visreg(fit = fit_bm, xvar = "TchangeHot",scale = "response", gg = TRUE, line=list(col="black"), rug = FALSE)+
  geom_point(data = SIwide, mapping = aes(x = TchangeHot, y = TchangeCold, shape = Age))+
  #scale_color_gradient(colours = mako(n=100, begin = 0.2, end = 0.8))+
  scale_shape_discrete(labels = c("1-week-old","8-week-old"))+
  xlab("Cloaca temperature change in heat [°C]")+
  ylab("Cloaca temperature change in cold [°C]")+
  #labs(size = "Body mass [kg]")+
  theme_pubr(base_size = 19)+
  theme(legend.title = element_blank())


#### incorporating morphological measurements ---------------------------------------------------------------------

# overview of chick size (with already summarised data)
ggplot(datTall, aes(y = Mass, x = Age, fill = sex))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(position = position_jitterdodge())+
  scale_fill_manual(values = sexcols)+
  theme_pubr()

dplyr::summarise(group_by(datRNAseq, Age), sd(Mass))


## plot chick size vs. deviation from median temp per treatment
datTall1 = mutate(group_by(filter(datTall, Age=="1week"),treatment), devmed = cloacat-mean(cloacat)) #only 1 week olds
ggplot(subset(datTall1, !is.na(Mass)), aes(x=Mass,y=devmed, color=treatment))+
  geom_point(size=2)+
  scale_color_manual(values= c("#4DAF4A","#377EB8" ,"#E41A1C"))+
  scale_y_continuous(name = "cloaca temp - median cloaca temp per treatment",limits = c(-2.2,2))+
  theme_pubr()
#ggsave("dev_med_cloacat_v_mass_1wk.pdf", path = figurepath)

#same for 8 week olds
datTall8 = mutate(group_by(filter(datTall, Age=="8week"),treatment), devmed = cloacat-mean(cloacat))
ggplot(subset(datTall8, !is.na(Mass)), aes(x=Mass,y=devmed, color=treatment))+
  geom_point(size=2)+
  scale_y_continuous(name = "cloaca temp - median cloaca temp per treatment")+
  scale_color_manual(values= c("#4DAF4A","#377EB8" ,"#E41A1C"))+
  theme_pubr()
#ggsave("dev_med_cloacat_v_mass_8wk.pdf", path = figurepath)
# the cold ones are on top here (esp when plotting mean) which I don't quite get, shouldn't each treatment have same amount above
# and below 0 -> something wrong? I think graph doesn't show me what I think it does

ggplot(filter(datTall, Age%in%c("1week")), aes(x=Mass, y=TchangeAv))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values = unname(tempcols[2:3]))+
  theme_pubr()
# same with mass:neck length ratio
ggplot(filter(datTall, Age%in%c("8week")), aes(x=Mass/NeckLength, y=TchangeAv))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values = unname(tempcols[2:3]))+
  theme_pubr()

ggplot(filter(SIwide, !is.na(TchangeCold) & Age == "1week"), aes(x=Mass, y = TchangeCold))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_pubr()

ggplot(filter(SIwide, Age == "1week"), aes(x=Mass, y = TchangeCold))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_pubr()

## statistical tests
SIwide1w <- filter(SIwide, Age == "1week")

fit1wH <- lm(TchangeHot ~ Mass, data = SIwide1w) # intercept is significant but mass is not
fit1wC <- lm(TchangeCold ~ Mass, data = SIwide1w) # not significant

simulationOutput <- simulateResiduals(fittedModel = fit1wC, n = 500) # simulate data from our model n times
# This method checks if your model is useful by seeing if it can even produce data that looks like the data you used to fit your model (i.e. your observations). 
plot(simulationOutput) # looks ok

summary(fit1wH)
summary(fit1wC)
# not significant



#### size correlations
ggplot(sampleInfo, aes(x=NeckLength, y= Mass, colour = Age))+
  geom_point(size=3.5, aes(shape=Age))+
  scale_color_manual(values= agecols)+
  labs(x = "Neck length [cm]", y = "Weight [kg]", title = "Weight and neck length of the chicks")+
  theme_pubr(base_size = 24)
#ggsave("FrontVBackGirth.pdf", path = figurepath)

p1=ggplot(sampleInfo, aes(x=NeckLength, y= Mass))+
  geom_point(size=2, aes(shape=Age))+
  scale_color_manual(values= tempcols)+
  theme_pubr()
p2=ggplot(sampleInfo, aes(x=NeckLength, y= FrontGirth))+
  geom_point(size=2, aes(shape=Age))+
  scale_color_manual(values= tempcols)+
  theme_pubr()
p3=ggplot(sampleInfo, aes(x=BackGirth, y= Mass))+
  geom_point(size=2, aes(shape=Age))+
  scale_color_manual(values= tempcols)+
  theme_pubr()
p4=ggplot(sampleInfo, aes(x=FrontGirth, y= Mass))+
  geom_point(size=2, aes(shape=Age))+
  scale_color_manual(values= tempcols)+
  theme_pubr()
ggarrange(p1,p2,p3,p4, nrow = 2, ncol = 2)
#best between girths and age and back girth, neck seems to be a bit plateauing


### visualise sex effects
# on morphology
psex1 = ggplot(filter(SIwide, Age == "1week"), aes(x = sex, y = Mass))+
  geom_boxplot()+
  geom_jitter(height = 0)+
  theme_pubr()
psex8 = ggplot(filter(SIwide, Age == "8week"), aes(x = sex, y = Mass))+
  geom_boxplot()+
  geom_jitter(height = 0)+
  theme_pubr()
ggarrange(psex1,psex8, nrow = 1)
# 8 weeks: most females heavier with a few very small ones so more variation, also same with girths but not neck length
# 1 week: basically the same

## on cloaca T
ggplot(datTalive, aes(x = treatment, y = cloacat, fill=sex))+
  geom_boxplot()+
  geom_jitter(position = position_jitterdodge(jitter.height = 0, jitter.width = 0.1))+
  scale_fill_manual(values = sexcols)+
  theme_pubr()


