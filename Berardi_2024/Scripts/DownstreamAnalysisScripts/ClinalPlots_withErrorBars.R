library(tidyverse)
library(reshape2)
library(binom)

setwd("C:/Users/jrhod/Desktop/GradSchool/ArtificialPigmentation/LandingPad2")

sampdata <- read.csv("C:/Users/jrhod/Desktop/GradSchool/ArtificialPigmentation/DownstreamCode/SupplementaryFiles/samps_10Nov2020.csv", header = TRUE)
sampdata <- sampdata %>% select(sampleId, nFlies, lat)
colnames(sampdata) <- c("sampleNames","nFlies","lat")

chroms = c("2L","3L","3R", "X")

pos2Llist <- c(20628607)
pos3Llist <- c(8888688, 1120350, 8878124, 10151119, 1102200)
pos3Rlist <- c(25139342)
posXlist <- c(9113092, 9113083)

darklist <- c("2L_20628607","3L_1120350","3L_10151119","3L_1102200", "3R_25139342")

for (mychrom in chroms){
  if (mychrom == "3L"){
    poslist <- pos3Llist
  } else if (mychrom == "3R"){
    poslist <- pos3Rlist
  } else if (mychrom == "X"){
    poslist <- posXlist
  } else if (mychrom == "2L"){
    poslist <- pos2Llist
  }





RDdata <- read.csv(str_c("DESTdata_",mychrom,"_pooledsampsRD_fivesamps_Rel5.csv"), header=TRUE)
AFdata <- read.csv(str_c("DESTdata_",mychrom,"_pooledsampsfreq_fivesamps_Rel5.csv"), header=TRUE)

RDdata <- RDdata %>% filter(POS %in% poslist)
AFdata <- AFdata %>% filter(POS %in% poslist)

RDdata[RDdata == "."] <- NA
AFdata[AFdata == "."] <- NA

sites <- AFdata %>% select(CHROM, POS, REF, ALT)
AFvals <- AFdata %>% select(-c(CHROM, POS, REF, ALT))

AFvals <- mutate_all(AFvals, function(x) as.numeric(as.character(x)))


afrows <- as.data.frame(colnames(AFvals))
colnames(afrows) <- "sampleNames"
afrows <- merge(afrows, sampdata, by = "sampleNames")


AFlabelled <- AFvals %>% t %>% as.data.frame() %>% cbind(afrows)

colnames(AFlabelled)[1:nrow(sites)] <- str_c(sites$CHROM,"_",sites$POS)

RDvals <- RDdata %>% select(-c(CHROM, POS))

RDvals <- mutate_all(RDvals, function(x) as.numeric(as.character(x)))

RDlabelled <- RDvals %>% t %>% as.data.frame() %>% cbind(afrows)

colnames(RDlabelled)[1:nrow(sites)] <- str_c(sites$CHROM,"_",sites$POS)



for (posnum in poslist){
  fullname <- str_c(mychrom,"_",posnum)
  
  if (fullname %in% darklist){
    AFlabelled[fullname] <- 1 - AFlabelled[fullname]
  }
  
  
  if (mychrom != "X"){
  RDlabelled$Neff <- unlist(round((RDlabelled[fullname]*(2*RDlabelled["nFlies"]))/(RDlabelled[fullname]+(2*RDlabelled["nFlies"]-1))))
  } else {
  RDlabelled$Neff <- unlist(round((RDlabelled[fullname]*RDlabelled["nFlies"])/(RDlabelled[fullname]+(RDlabelled["nFlies"]-1))))
  }
  
  RDuse <- RDlabelled %>% select(lat, Neff)
  AFuse <- AFlabelled %>% select(lat,fullname)
  
  fullgraph <- merge(AFuse,RDuse, by="lat")
  fullgraph <- na.omit(fullgraph)
  
  confinterval <- binom.confint(x=unlist(round(fullgraph[fullname]*fullgraph["Neff"])),n=unlist(fullgraph["Neff"]), conf.level = 0.8 ,methods = "exact")
  fullgraph$lower <- confinterval$lower
  fullgraph$upper <- confinterval$upper
  
  
  lightplot <- ggplot(data = fullgraph) + 
    geom_line(aes(x=lat,y=!!as.symbol(fullname))) +
    geom_point(aes(x=lat,y=!!as.symbol(fullname)),size=2) +
    geom_errorbar(aes(x=lat,ymin=lower,ymax=upper),width=0.1,linewidth=0.8) +
    scale_x_continuous(name="Latitude", limits =c(25,45), breaks=c(25.50, 33.95, 38.00, 39.88, 42.45),labels=c("25.50"="25.50\nFL", "33.95"="33.95\nGA", "38.00" = "38.00\nVA", "39.88" = "39.88\nPA", "42.45" = "42.45\nMA")) + 
    labs(y="Frequency of light-associated allele",title=str_c(mychrom,":",posnum)) +
    theme_classic(base_size=20)
  
  ggsave(str_c("C:/Users/jrhod/Desktop/GradSchool/ArtificialPigmentation/LandingPadPlots2/Clinal_LightSNP_FreqWithErrorBars_",fullname,"_biggerText.png"),plot=lightplot,device = png, height = 6, width = 12, units = "in")
  
  
  }}
