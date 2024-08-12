.libPaths("/home/users/jarhodes/R/x86_64-pc-linux-gnu-library/4.2")
library(tidyverse)
library(reshape2)
library(coxed)
setwd("/scratch/groups/dpetrov/jess/pigmentselect")

args = commandArgs(trailingOnly=TRUE)

printpval <- data.frame(HIT=character(),NUM=numeric(),SEASON=numeric(),RATIO=numeric(),PVALSEASON=numeric(),PVALRATIO=numeric())

chrom = args[1]

#file created from "MetricCalcSeason_MatchedSNPs_bootstrap.py"
data <- read.table(str_c("DirectionSeason_MatchedSNPs_Individual_",chrom,".txt"),header=TRUE,stringsAsFactors=FALSE)

data$ABSSEASON <- abs(data$TOTALSEASON)
data$ABSNORMAL <- abs(data$TOTALNORMAL)

hits <- unique(data$HITNAME)

#create a null distribution for each hit and each of 100 matched SNP groups and determine p-value of hit 100 times
for (hit in hits){
	samp <- data[which(data$HITNAME==hit),]
	hitvalue <- samp[which(samp$TYPE=="hit"),]
	season <- hitvalue$TOTALSEASON
	ratio <- hitvalue$TOTALNORMAL
	for (num in seq(1,100)){
	numsamp <- samp[which(samp$SAMPNUM==num),]
	matched <- numsamp[which(numsamp$TYPE!="hit"),]
	permCDFseason <- ecdf(matched$ABSSEASON)
	pvalueseason <- 1-permCDFseason(hitvalue$ABSSEASON)
	permCDFratio <- ecdf(matched$ABSNORMAL)
	pvalueratio <- 1-permCDFratio(hitvalue$ABSNORMAL)

	entry <- c(hit,num,season,ratio,pvalueseason,pvalueratio)
	printpval <- rbind(printpval,entry)
}
}

colnames(printpval) <- c("HIT","NUM","SEASON","RATIO","PVALSEASON","PVALRATIO")

write.table(printpval, file = str_c("SeasonPvals_Unadjusted_",chrom,".txt"), quote=FALSE)

printpval <- as.data.frame(printpval)

printpvalQuan <- data.frame(HIT=character(),SEASON=numeric(),RATIO=numeric(),PVALSEASON10=numeric(),PVALSEASON50=numeric(),PVALSEASON90=numeric(),PVALRATIO10=numeric(),PVALRATIO50=numeric(),PVALRATIO90=numeric())

#get quantiles of p-values generated from null distributions of matched SNPs per candidate SNP
for (hit in hits){
	pvals <- printpval[which(printpval$HIT==hit),]
	colnames(pvals) <- c("HIT","NUM","SEASON","RATIO","PVALSEASON","PVALRATIO")
	eightyPercentileSeason <- quantile(as.numeric(pvals$PVALSEASON),probs=c(0.1,0.5,0.9))
	eightyPercentileRatio <- quantile(as.numeric(pvals$PVALRATIO),probs=c(0.1,0.5,0.9))
	season <- unique(pvals$SEASON)
	ratio <- unique(pvals$RATIO)
	entry <- c(hit,season,ratio,unname(eightyPercentileSeason[1]),unname(eightyPercentileSeason[2]),unname(eightyPercentileSeason[3]),unname(eightyPercentileRatio[1]),unname(eightyPercentileRatio[2]),unname(eightyPercentileRatio[3]))
	printpvalQuan <- rbind(printpvalQuan,entry)

}
print(dim(printpvalQuan))

colnames(printpvalQuan) <- c("HIT", "SEASON", "RATIO", "PVALSEASON10", "PVALSEASON50", "PVALSEASON90", "PVALRATIO10", "PVALRATIO50", "PVALRATIO90")

write.table(printpvalQuan, file = str_c("SeasonPvals_Quantiles_",chrom,".txt"), quote=FALSE)