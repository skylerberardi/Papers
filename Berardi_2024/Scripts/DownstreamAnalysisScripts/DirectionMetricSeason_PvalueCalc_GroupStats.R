.libPaths("/home/users/jarhodes/R/x86_64-pc-linux-gnu-library/4.2")
library(tidyverse)
library(reshape2)
library(coxed)
setwd("/scratch/groups/dpetrov/jess/pigmentselect")

args = commandArgs(trailingOnly=TRUE)

printpval <- data.frame(HIT=character(),NUM=numeric(),SEASON=numeric(),RATIO=numeric(),PVALSEASON=numeric(),PVALRATIO=numeric())

chroms <- c("2L","2R","3L","3R","X")

#read in matched snps & candidate SNPs
for (chrom in chroms){
data <- read.table(str_c("DirectionSeason_GroupPermutations_",chrom,".txt"),header=TRUE,stringsAsFactors=FALSE)
datahits <- read.table(str_c("DirectionSeason_MatchedSNPs_Individual_",chrom,".txt"),header=TRUE,stringsAsFactors=FALSE)
hits <- datahits[which(datahits$TYPE=="hit"),]

if (chrom=="2L"){
	groupperms<-data
	hitsdata <- hits
} else{
	groupperms<-rbind(groupperms,data)
	hitsdata <- rbind(hitsdata,hits)
}

}


groupperms$ABSSEASON <- abs(as.numeric(groupperms$TOTALSEASON))
hitsdata$ABSSEASON <- abs(as.numeric(hitsdata$TOTALSEASON))

groupperms$ABSNORMAL <- abs(as.numeric(groupperms$TOTALNORMAL))
hitsdata$ABSNORMAL <- abs(as.numeric(hitsdata$TOTALNORMAL))

for (i in seq(1,1000)){
	group <- groupperms[which(groupperms$SAMPNUM==i),]
	distributionseason <- group$ABSSEASON
	distributionratio <- group$ABSNORMAL
	medseason <- median(distributionseason)
	avgseason <- mean(distributionseason)
	medratio <- median(distributionratio)
	avgratio <- mean(distributionratio)
	if (i==1){
		medvectorseason <- c(medseason)
		avgvectorseason <- c(avgseason)
		medvectorratio <- c(medratio)
		avgvectorratio <- c(avgratio)
	} else{
		medvectorseason <- append(medvectorseason,medseason)
		avgvectorseason <- append(avgvectorseason,avgseason)
		medvectorratio <- append(medvectorratio,medratio)
		avgvectorratio <- append(avgvectorratio,avgratio)
	}

}

hitmedseason <- median(hitsdata$ABSSEASON)
hitavgseason <- mean(hitsdata$ABSSEASON)
hitmedratio <- median(hitsdata$ABSNORMAL)
hitavgratio <- mean(hitsdata$ABSNORMAL)

#create null distributions of the matched SNPs' group stats
medCDFseason <- ecdf(medvectorseason)
avgCDFseason <- ecdf(avgvectorseason)
medCDFratio <- ecdf(medvectorratio)
avgCDFratio <- ecdf(avgvectorratio)

#generate pvalue of candidate SNP group's group stats
medpvalseason <- 1-medCDFseason(hitmedseason)
avgpvalseason <- 1-avgCDFseason(hitavgseason)
medpvalratio <- 1-medCDFratio(hitmedratio)
avgpvalratio <- 1-avgCDFratio(hitavgratio)

print("season pval - med, avg")
print(medpvalseason)
print(avgpvalseason)

print("ratio pval - med, avg")
print(medpvalratio)
print(avgpvalratio)


#prints out ks test results using the total seasonal metric and the normalized (averaged) seasonal metric
kstestseason <- ks.test(as.numeric(hitsdata$ABSSEASON),as.numeric(groupperms$ABSSEASON),alternative="less")
print("ks season")
print(kstestseason)

kstestratio <- ks.test(as.numeric(hitsdata$ABSNORMAL),as.numeric(groupperms$ABSNORMAL),alternative="less")
print("ks ratio")
print(kstestratio)

#compares the tail end of distribution of matched SNPs vs candidate snps
print("percentiles season")
print(quantile(as.numeric(hitsdata$ABSSEASON),probs=c(0.8,0.9,0.95,0.96,0.97,0.98,0.99)))
print(quantile(as.numeric(groupperms$ABSSEASON),probs=c(0.8,0.9,0.95,0.96,0.97,0.98,0.99)))

print("percentiles ratio")
print(quantile(as.numeric(hitsdata$ABSNORMAL),probs=c(0.8,0.9,0.95,0.96,0.97,0.98,0.99)))
print(quantile(as.numeric(groupperms$ABSNORMAL),probs=c(0.8,0.9,0.95,0.96,0.97,0.98,0.99)))