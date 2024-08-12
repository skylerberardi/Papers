.libPaths("/home/users/jarhodes/R/x86_64-pc-linux-gnu-library/4.2")
library(tidyverse)
library(reshape2)
library(coxed)
setwd("/scratch/groups/dpetrov/jess/pigmentselect")

chroms <- c("2L","2R","3L","3R","X")

for (chrom in chroms){
	data <- read.table(str_c("GroupPermutations_LinRegBetas_Orch16_",chrom,".txt"),header=TRUE)
	hitdata <- read.table(str_c("GroupHits_LinRegBetas_Orch16_",chrom,".txt"),header=TRUE)
	if (chrom=="2L"){
		fulldata <- data
		fullhit <- hitdata
		} else{
			fulldata<-rbind(fulldata,data)
			fullhit<-rbind(fullhit,hitdata)
		}
}

fulldata <- as.data.frame(fulldata)
fullhit <- as.data.frame(fullhit)

colnames(fulldata) <- c("GROUP","CHROM","POS","BETA")
colnames(fullhit) <- c("GROUP","CHROM","POS","BETA")

fulldata$BETAABS<-abs(fulldata$BETA)
fullhit$BETAABS<-abs(fullhit$BETA)

for (i in seq(1,1000)){
	group <- fulldata[which(fulldata$GROUP==i),]
	meanperm <- mean(group$BETAABS)
	medianperm <- median(group$BETAABS)

	if (i==1){
		meanvec <- c(meanperm)
		medianvec <- c(medianperm)
		} else{
		meanvec <- append(meanvec,meanperm)
		medianvec <- append(medianvec,medianperm)
		}


	}

hitmean <- mean(fullhit$BETAABS)
hitmedian <- median(fullhit$BETAABS)

#makes a null distribution of matched SNPs
medianCDF <- ecdf(medianvec)
meanCDF <- ecdf(meanvec)

medianpval <- 1-medianCDF(hitmedian)
meanpval <- 1-meanCDF(hitmean)

print("Median, mean")
print(medianpval)
print(meanpval)

kstestcage <- ks.test(as.numeric(fullhit$BETAABS),as.numeric(fulldata$BETAABS),alternative="less")
print(kstestcage)

#makes and saves a plot of the cdf of candidates and matched SNPs
cdforchplot <- ggplot() + stat_ecdf(data=fulldata, aes(fulldata$BETAABS),geom = "step") + stat_ecdf(data=fullhit, aes(fullhit$BETAABS),geom = "step")
ggsave("CdfOrchPlot.png",plot=cdforchplot,device=png)
print(quantile(as.numeric(fullhit$BETAABS),probs=c(0.8,0.9,0.95,0.96,0.97,0.98,0.99)))
print(quantile(as.numeric(fulldata$BETAABS),probs=c(0.8,0.9,0.95,0.96,0.97,0.98,0.99)))
