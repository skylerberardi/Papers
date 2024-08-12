library(tidyverse)
library(reshape2)
setwd("C:/Users/jrhod/Desktop/GradSchool/ArtificialPigmentation")

chroms <- c("2L","2R","3L","3R","X")

load("orch16_spring.RData")

fulldf <- cbind(sites,afmat)

graphingSNPs <- read.table("DownstreamCode/SupplementaryFiles/GraphingInfoSkylerSNPs.tsv", header = TRUE)
graphingSNPs <- graphingSNPs %>% filter(ANALYSIS=="orchard") %>% separate_wider_delim(HIT,"_",names = c("chrom","pos"))
graphingSNPs$pos <- as.numeric(graphingSNPs$pos)

for (ii in seq(1,nrow(graphingSNPs))){
  graphdf <- as.data.frame(matrix(nrow=0,ncol=3))
  graphreps <- as.data.frame(matrix(nrow=0,ncol=3))
  snp <- graphingSNPs[ii,]
  mychrom <- unlist(snp$chrom)
  mypos <- unlist(snp$pos)
  freqsnp <- fulldf %>% filter(chrom==mychrom,pos==mypos)
  for (TPnum in seq(1,3)){
    TPname <- str_c("TP",TPnum)
    TPfreq <- freqsnp %>% select(starts_with(TPname))
    TPmean <- mean(unlist(TPfreq))
    TPsd <- sd(unlist(TPfreq))/sqrt(length(TPfreq))
    if (unlist(snp$ALT) != unlist(snp$LIGHT)){
      TPmean <- 1-TPmean
    }
    entry <- c(TPnum,TPmean,TPsd)
    graphdf <- rbind(graphdf,entry)
    
    for(cagenum in seq(1,10)){
      cagename <- str_c("S",cagenum)
      repfreq <- TPfreq %>% select(ends_with(cagename))
      if (ncol(repfreq)>0){
        repfreq <- repfreq %>% unlist() %>% unname()
        if (unlist(snp$ALT) != unlist(snp$LIGHT)){
          repfreq = 1-repfreq
        }
      repentry <- c(TPnum,repfreq,cagename)
      graphreps <- rbind(graphreps, repentry)}
    }
  }
  colnames(graphdf) <- c("TP","FREQ","SDERR")
  colnames(graphreps) <- c("TP", "FREQ", "CAGE")
  graphreps$FREQ <- as.numeric(graphreps$FREQ)
  graphreps$TP <- as.numeric(graphreps$TP)
  lighttrajplot <- ggplot(data = graphdf) + 
    geom_line(data=graphreps,aes(x=TP,y=FREQ, group=CAGE),color="grey70",show.legend = FALSE) +
    geom_line(aes(x=TP,y=FREQ),linewidth=1) + geom_point(aes(x=TP,y=FREQ),size=2) +
    geom_errorbar(aes(x=TP,ymin=FREQ-SDERR,ymax=FREQ+SDERR),width=0.1,linewidth=0.8) +
    scale_x_continuous(name="Timepoint",breaks=c(1,2,3),labels=c("1"="TP 1","2"="TP 2","3"="TP 3")) + 
    labs(y="Frequency of light-associated allele",title=str_c(mychrom,":",mypos)) +
    theme_classic()
  ggsave(str_c("LandingPad3/Orchard_LightSNP_Freq_IndTrajs_",mychrom,"_",mypos,"_stderr.png"),plot = lighttrajplot,device=png)
  lightplot <- ggplot(data = graphdf) + 
    geom_line(aes(x=TP,y=FREQ),linewidth=1) + geom_point(aes(x=TP,y=FREQ),size=2) +
    geom_errorbar(aes(x=TP,ymin=FREQ-SDERR,ymax=FREQ+SDERR),width=0.1,linewidth=0.8) +
    scale_x_continuous(name="Timepoint",breaks=c(1,2,3),labels=c("1"="TP 1","2"="TP 2","3"="TP 3")) + 
    labs(y="Frequency of light-associated allele",title=str_c(mychrom,":",mypos)) +
    theme_classic(base_size=20)
  ggsave(str_c("LandingPad3/Orchard_LightSNP_Freq_",mychrom,"_",mypos,"_stderr_biggerText.png"),plot = lightplot,device=png)
}