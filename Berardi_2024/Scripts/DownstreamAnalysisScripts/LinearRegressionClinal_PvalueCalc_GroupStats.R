.libPaths("/home/users/jarhodes/R/x86_64-pc-linux-gnu-library/4.2")
library(tidyverse)
library(reshape2)
library(coxed)
setwd("/scratch/groups/dpetrov/jess/pigmentselect")

args = commandArgs(trailingOnly=TRUE)

printpval <- data.frame(HIT=character(),NUM=numeric(),SLOPE=numeric(),PVAL=numeric())

chroms <- c("2L","2R","3L","3R","X")


for (chrom in chroms){
data <- read.table(str_c("GroupPermutations_LinearRegressionResults_Clinal_fivesamps_",chrom,".txt"),header=TRUE,stringsAsFactors=FALSE)
datahits <- read.table(str_c("LinearRegressionResults_Clinal_fivesamps_",chrom,".txt"),header=TRUE,stringsAsFactors=FALSE)
hits <- datahits[which(datahits$TYPE=="hit"),]

if (chrom=="2L"){
	groupperms<-data
	hitsdata <- hits
} else{
	groupperms<-rbind(groupperms,data)
	hitsdata <- rbind(hitsdata,hits)
}

}






groupperms$ABSSLOPE <- abs(as.numeric(groupperms$SLOPE))
hitsdata$ABSSLOPE <- abs(as.numeric(hitsdata$SLOPE))

#uses all matched SNP groups to make a null distribution of their group stats
#then generates a pvalue for the candidate SNP group stats

for (i in seq(1,1000)){
	group <- groupperms[which(groupperms$SAMPNUM==i),]
	distribution <- group$ABSSLOPE
	med <- median(distribution)
	avg <- mean(distribution)
	if (i==1){
		medvector <- c(med)
		avgvector <- c(avg)
	} else{
		medvector <- append(medvector,med)
		avgvector <- append(avgvector,avg)
	}

}

hitmed <- median(hitsdata$ABSSLOPE)
hitavg <- mean(hitsdata$ABSSLOPE)

medCDF <- ecdf(medvector)
avgCDF <- ecdf(avgvector)

medpval <- 1-medCDF(hitmed)
avgpval <- 1-avgCDF(hitavg)

print(medpval)
print(avgpval)

#prints out the the kstest and the top quantiles for the candidate & matched SNP distributions

kstest <- ks.test(as.numeric(hitsdata$ABSSLOPE),as.numeric(groupperms$ABSSLOPE),alternative="less")
print(kstest)

print(quantile(as.numeric(hitsdata$ABSSLOPE),probs=c(0.8,0.9,0.95,0.96,0.97,0.98,0.99)))
print(quantile(as.numeric(groupperms$ABSSLOPE),probs=c(0.8,0.9,0.95,0.96,0.97,0.98,0.99)))

print(hitmed)
print(hitavg)

print(medvector)
print(avgvector)