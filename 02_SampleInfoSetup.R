### Gathering all information about samples
### Lea Huber

#### load packages ####
library("readxl")
library("tidyverse")
library("ggpubr")
library("RColorBrewer")
library("ggExtra")
library(lme4)

#### preliminaries ####
setwd("/Users/Lea/Library/CloudStorage/OneDrive-Aarhusuniversitet/LeaH/ChamberTempExp/data/working")



#### import ####
infoday1 = read_excel("SampleIDs.xlsx", sheet = 2)
infoday2 = read_excel("SampleIDs.xlsx", sheet = 3)
infoday3 = read_excel("SampleIDs.xlsx", sheet = 4)

datTday1 = read_excel("SampleIDs.xlsx", sheet = 6)
datTday2 = read_excel("SampleIDs.xlsx", sheet = 7)
datTday3 = read_excel("SampleIDs.xlsx", sheet = 8)
# add sex
datTday1 = left_join(datTday1,select(infoday1, c(chickno,sex)), by = "chickno")
datTday2 = left_join(datTday2,select(infoday2, c(chickno,sex)), by = "chickno")
datTday3 = left_join(datTday3,select(infoday3, c(chickno,sex)), by = "chickno")

datTall = bind_rows(subset(datTday1,select = -Time),subset(datTday2,select = -Time))
datTall = bind_rows(datTall, subset(datTday3,select = -c(Time,NoDNA)))
datTall$chickno = as.factor(datTall$chickno)

# sample number info only complete
datSampleNrs <- read_excel("samplesize.xlsx", sheet = 3)
datSampleNrs$group = as.factor(datSampleNrs$group)

### add weight and size
datmass = read_excel("morph_measures.xlsx")
datmass = rename(datmass, chickno = Chickno, date = Date)
datmass$Mass = as.numeric(datmass$Mass)
datmass$date = as.POSIXct(datmass$date, format = "%d_%m_%y")
datmass$chickno = as.factor(datmass$chickno)
datTall = left_join(datTall, datmass, by = "chickno")
## calculate some mass stuff
datTall <- datTall %>% group_by(Age) %>% mutate(massMean = mean(na.omit(Mass)))
datTall$massMC <- datTall$Mass - datTall$massMean
# scale()
datTall <- datTall %>% group_by(Age) %>% mutate(massScaled = scale(Mass))
datTall <- ungroup(datTall)

### add death date
deathdat = read_excel("Chick Trial Growth Rate.xlsx",sheet = 6, range=cell_cols(c("K","Q")))
deathdat = select(deathdat, c("chickno","mort_dat"))
datTall <- plyr::join(datTall, deathdat, by = "chickno")
#write.csv(datTall, "datTall.csv")

## add chamber temperature
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
# should take same time for each treatment each day, and always same length of window 1.5 h so starting ~1h before sampling
# hot: starting 13:00 till 14:30, so 257-347, 232-322, 242
# cold: starting 13:30 till 15:00, so 287-377, 261-351, 272,
# benign: starting 13:45 till 15:15, so 302-392, 276, 287

# Now I think it should rather be from around 11:00, so I am subtracting 120 from all the start times

tempmeans = c()
tempmeans["hot1"] = mean(tempdata$hot1$temp[(257-120):(257+90)])
tempmeans["hot2"] = mean(tempdata$hot2$temp[(232-120):(232+90)])
tempmeans["hot3"] = mean(tempdata$hot3$temp[(242-120):(242+90)])
tempmeans["cold1"] = mean(tempdata$cold1$temp[(278-120):(278+90)])
tempmeans["cold2"] = mean(tempdata$cold2$temp[(261-120):(261+90)])
tempmeans["cold3"] = mean(tempdata$cold3$temp[(272-120):(272+90)])
tempmeans["benign1"] = mean(tempdata$benign1$temp[(302-120):(302)])
tempmeans["benign2"] = mean(tempdata$benign2$temp[(276-120):(276+90)])
tempmeans["benign3"] = mean(tempdata$benign3$temp[(287-120):(287+90)])

## correlate mean cloacal temp with actual temperature (mean or at sampling time)
temps= c(rep(tempmeans["hot1"],10),rep(tempmeans["cold1"],10),rep(tempmeans["benign1"],10),
         rep(tempmeans["hot2"],10),rep(tempmeans["cold2"],10),rep(tempmeans["benign2"],10),
         rep(tempmeans["hot3"],10),rep(tempmeans["cold3"],10),rep(tempmeans["benign3"],10))
datTall$chamberT= temps


#### calculating change in cloacal temperature ---------------------------------------------------------------------

## calculate dev from benign per chick
# calculate plasticity per individual for hot and cold so standardize on benign
datTwide = pivot_wider(dplyr::select(datTall,-c(Date, SampleNr, Samplingtime,chamberT,comments)), # might have to add chamberT and devmed and comments
                       names_from = treatment, values_from = cloacat) # wide with one row per chick
datTwide = mutate(datTwide,TchangeHot = hot-benign)
datTwide = mutate(datTwide,TchangeCold = cold-benign)
datTwide <- mutate(datTwide, TchangeAv = (abs(TchangeHot)+abs(TchangeCold))/2)
# add to data frame
datTall <- left_join(datTall, dplyr::select(datTwide, c(chickno, TchangeHot:TchangeAv)), by = "chickno")

datTlong = pivot_longer(datTwide,cols = c(TchangeHot,TchangeCold),names_to = "treatChange",values_to = "TempDev")
datTlong$treatChange <- as.factor(datTlong$treatChange)

# mass:length
datTlong = mutate(group_by(datTlong, by=chickno), massVlength = Mass/NeckLength)


#### Make subsets of datTall ------------------------------------------------------------------------
## some improvements to datTall
# make factors
datTall$treatment <- as.factor(datTall$treatment)
datTall$group <- as.factor(datTall$group)
datTall$Age <- as.factor(datTall$Age)
datTall$sex <- as.factor(datTall$sex)
datTall$chickno <- droplevels(datTall$chickno) # drop old factor levels

# remove comments
datTall <- subset(datTall, select = -comments)

wrongdate <- datTall[14,"mort_dat"] # remove weird date, change to NA
datTall[which(datTall$mort_dat == wrongdate),"mort_dat"] <- NA


# filter out chicks that died next day
deadday2 = c(8247,7772)
deadday3 = c(8207,7830)
datTalive = filter(datTall, !(chickno%in%deadday2))
datTalive$Date = as.character(datTalive$Date)
datTalive = filter(datTalive, !(Date == "2023-10-26" & chickno%in%deadday3))

## Filter out chicks for RNAseq analysis
SampleType = data.frame(chickno= c(7742,7755,7760,7761,7778,7785,7807,7808,7813,7831,8202,8204,8209,8215,
                                   8228,8231,8235,8238,8242,8244,8245),
                        sampletype = c("both","both","both","later","both","both","both","both","both","Tempus","both","both",
                                       "both","both","both","both","both","both","both","both","both"))
sampleInfo = filter(datTall, chickno%in%SampleType$chickno)
chicks <- c("7742", "7755", "7760", "7761", "7778", "7785" ,"7807", "7808" ,"7813" ,"7831", "8202", "8204", "8209", "8215","8228", 
            "8231" ,"8235" ,"8238" ,"8242", "8244", "8245")
sampleInfo = filter(datTall, chickno%in%chicks)
sampleInfo <- droplevels(sampleInfo)

#### Save data for SI and other purposes -------------------------

# save(datTall, datTjuv,growthdatafull, growthdatwide,tempdata, file = "ChamberTempSI.RData")
# load("ChamberTempSI.RData")
# save(datTall,datTalive,sampleInfo, file = "ChamberTempSI_new.RData")




