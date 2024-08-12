.libPaths("/home/users/jarhodes/R/x86_64-pc-linux-gnu-library/4.2")
library(tidyverse)
library(reshape2)
library(coxed)
setwd("/scratch/groups/dpetrov/jess/pigmentselect")

chroms <- c("2L","2R","3L","3R","X")

args = commandArgs(trailingOnly=TRUE)

#load in orchard data and calculating the mean frequency for every SNP at each timepoint across the cages
load("orch16_spring.RData")

TP1 <- samps[which(samps$tpt==1),1]
TP1data <- afmat[,TP1]

TP1mean <- as.data.frame(apply(TP1data,1,mean))

TP2 <- samps[which(samps$tpt==2),1]
TP2data <- afmat[,TP2]

TP2mean <- as.data.frame(apply(TP2data,1,mean))

TP3 <- samps[which(samps$tpt==3),1]
TP3data <- afmat[,TP3]

TP3mean <- as.data.frame(apply(TP3data,1,mean))

data <- cbind(TP1mean,TP2mean,TP3mean)

colnames(data) <- c("TP1","TP2","TP3")

#calculates linear regression of timepoint vs. mean frequency for each SNP
linreg <- function(l) {
	data = data.frame(y=l,x=c(1,2,3))
	model = lm(y ~ x, data=data)
	model$coefficients["x"]
}

a <- apply(data,1,linreg)
betas <- as.data.frame(t(as.data.frame(as.list(a))))

#writes the linear regression results & the base frequencies to text files for downstream code
fullbeta <- as.data.frame(cbind(sites,betas))

colnames(fullbeta) <- c("CHROM","POS","BETA")

write.table(fullbeta, file = "Orch16_ThreeTimepoint_LinRegBetas_Full.txt", quote=FALSE)

fulldata <- as.data.frame(cbind(sites,data))

colnames(fulldata) <- c("CHROM","POS","TP1","TP2","TP3")

write.table(fulldata, file = "Orch16_ThreeTimepoint_Freqs_Full.txt", quote=FALSE)

for (chrom in chroms){
	chrombeta <- fullbeta[which(fullbeta$CHROM==chrom),]
	chromdata <- fulldata[which(fulldata$CHROM==chrom),]
	write.table(chrombeta, file = str_c("Orch16_ThreeTimepoint_LinRegBetas_",chrom,".txt"), quote=FALSE)
	write.table(chromdata, file = str_c("Orch16_ThreeTimepoint_Freqs_",chrom,".txt"), quote=FALSE)

}

