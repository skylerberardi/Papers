library(tidyverse)
library(reshape2)
library(binom)

setwd("C:/Users/jrhod/Desktop/GradSchool/ArtificialPigmentation/LandingPad2")

darklist <- c("3R_25139132")

sampdata <- read.csv("samps_10Nov2020.csv", header = TRUE)
sampdata <- sampdata %>% select(sampleId, nFlies)
colnames(sampdata) <- c("sampleNames","nFlies")

chroms = c("3L","3R","X")

pos3Llist <- c(9020823, 1904606, 1271873)
posXlist <- c(9167319, 9125023)
pos3Rlist <- c(25139132)

darklist <- c("3R_25139132")

for (mychrom in chroms){
  if (mychrom == "3L"){
    poslist <- pos3Llist
  } else if (mychrom == "3R"){
    poslist <- pos3Rlist
  } else if (mychrom == "X"){
    poslist <- posXlist
  }

RDdata <- read.table(str_c("DESTdata_",mychrom,"_pooledsampsRD_LinvillaSeason_correctorder_Rel5.txt"), header=TRUE)
AFdata <- read.table(str_c("DESTdata_",mychrom,"_pooledsampsfreq_LinvillaSeason_correctorder_Rel5.txt"), header=TRUE)

RDdata <- RDdata %>% filter(POS %in% poslist)
AFdata <- AFdata %>% filter(POS %in% poslist)

RDdata[RDdata == "."] <- NA
AFdata[AFdata == "."] <- NA

sites <- AFdata %>% select(CHROM, POS, REF, ALT)
AFvals <- AFdata %>% select(-c(CHROM, POS, REF, ALT))

AFvals <- mutate_all(AFvals, function(x) as.numeric(as.character(x)))

afrows <- as.data.frame(colnames(AFvals))
colnames(afrows) <- "sample"
afrows <- afrows %>% separate_wider_delim(sample,"_",names=c(NA,NA,"year","season"))
afrows$year <- as.numeric(afrows$year)
afrows <- afrows %>% mutate(year = ifelse(season=="spring",year,year+0.5))

AFlabelled <- AFvals %>% t %>% as.data.frame() %>% cbind(afrows)

colnames(AFlabelled)[1:nrow(sites)] <- str_c(sites$CHROM,"_",sites$POS)

RDvals <- RDdata %>% select(-c(CHROM, POS))

RDvals <- mutate_all(RDvals, function(x) as.numeric(as.character(x)))

sampleNames <- as.data.frame(colnames(RDvals))
names(sampleNames) <- "sampleNames"
RDrows <- as.data.frame(colnames(RDvals))
colnames(RDrows) <- "sample"
RDrows <- RDrows %>% separate_wider_delim(sample,"_",names=c(NA,NA,"year","season"))
RDrows$year <- as.numeric(RDrows$year)
RDrows <- RDrows %>% mutate(year = ifelse(season=="spring",year,year+0.5))

RDlabelled <- RDvals %>% t %>% as.data.frame() %>% cbind(RDrows,sampleNames)

colnames(RDlabelled)[1:nrow(sites)] <- str_c(sites$CHROM,"_",sites$POS)

RDlabelled <- left_join(RDlabelled, sampdata, by = "sampleNames")


AFlabelled <- AFlabelled %>% mutate(fullyear = floor(year))
AFfall <- AFlabelled %>% filter(season=="fall")
AFspring <- AFlabelled %>% filter(season=="spring")
names(AFspring)[names(AFspring) != 'fullyear'] <- str_c(names(AFspring)[names(AFspring) != 'fullyear'],"_spring")
names(AFfall)[names(AFfall) != 'fullyear'] <- str_c(names(AFfall)[names(AFfall) != 'fullyear'],"_fall")
AFsegments <- merge(AFspring,AFfall, by="fullyear")

for (posnum in poslist){
  fullname <- str_c(mychrom,"_",posnum)
  fullnamefall <- str_c(fullname,"_fall")
  fullnamespring <- str_c(fullname,"_spring")
  
  if (fullname %in% darklist){
    AFlabelled[fullname] <- 1 - AFlabelled[fullname]
    AFsegments[fullnamefall] <- 1 - AFsegments[fullnamefall]
    AFsegments[fullnamespring] <- 1 - AFsegments[fullnamespring]
  }
  
  AFsegGraph <- AFsegments %>% select(fullnamespring,fullnamefall,year_spring,year_fall)
  AFsegGraph <- na.omit(AFsegGraph)
  
  if (mychrom != "X"){
  RDlabelled$Neff <- unlist(round((RDlabelled[fullname]*(2*RDlabelled["nFlies"]))/(RDlabelled[fullname]+(2*RDlabelled["nFlies"]-1))))
  } else {
  RDlabelled$Neff <- unlist(round((RDlabelled[fullname]*RDlabelled["nFlies"])/(RDlabelled[fullname]+(RDlabelled["nFlies"]-1))))
  }
  
  RDuse <- RDlabelled %>% select(year, Neff)
  AFuse <- AFlabelled %>% select(year,fullname)
  
  fullgraph <- merge(AFuse,RDuse, by="year")
  fullgraph <- na.omit(fullgraph)
  
  confinterval <- binom.confint(x=unlist(round(fullgraph[fullname]*fullgraph["Neff"])),n=unlist(fullgraph["Neff"]), conf.level = 0.8 ,methods = "exact")
  fullgraph$lower <- confinterval$lower
  fullgraph$upper <- confinterval$upper
  
  
  lightplot <- ggplot(data = fullgraph) + 
    geom_point(aes(x=year,y=!!as.symbol(fullname)),size=2) +
    geom_errorbar(aes(x=year,ymin=lower,ymax=upper),width=0.1,linewidth=0.8) +
    geom_segment(data = AFsegGraph, aes(x=year_spring,y=!!as.symbol(fullnamespring),xend=year_fall,yend=!!as.symbol(fullnamefall))) +
    scale_x_continuous(name="Collection Season",breaks=c(9,9.5,10,10.5,11,11.5,12,12.5,14,14.5,15,15.5,16,16.5),labels=c("9"="Spring\n2009","9.5"="Fall\n2009","10"="Spring\n2010","10.5"="Fall\n2010","11"="Spring\n2011","11.5"="Fall\n2011","12"="Spring\n2012","12.5"="Fall\n2012","14"="Spring\n2014","14.5"="Fall\n2014","15"="Spring\n2015","15.5"="Fall\n2015","16"="Spring\n2016","16.5"="Fall\n2016")) + 
    labs(y="Frequency of light-associated allele",title=str_c(mychrom,":",posnum)) +
    theme_classic(base_size=18)
  
  ggsave(str_c("C:/Users/jrhod/Desktop/GradSchool/ArtificialPigmentation/LandingPadPlots2/Seasonal_LightSNP_FreqWithErrorBars_",fullname,"_biggerText.png"),plot=lightplot,device = png, height = 6, width = 12, units = "in")
  
  
  }}
