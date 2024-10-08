.libPaths("/home/users/jarhodes/R/x86_64-pc-linux-gnu-library/4.2")
library(tidyverse)
library(reshape2)
library(coxed)
setwd("/scratch/groups/dpetrov/jess/pigmentselect")

args = commandArgs(trailingOnly=TRUE)

printpval <- data.frame(HIT=character(),NUM=numeric(),SLOPE=numeric(),PVAL=numeric())

#takes chromosome as an argument so I could run as separate jobs on the cluster
chrom = args[1]

#file generated by "LinearRegressionClinal_MatchedSNPs_bootstrapSamps_efficient.py"

data <- read.table(str_c("LinearRegressionResults_Clinal_fivesamps_",chrom,".txt"),header=TRUE,stringsAsFactors=FALSE)

data$ABSSLOPE <- abs(data$SLOPE)

hits <- unique(data$HITNAME)

#creates a distribution from matched snps associated with a candidate SNP and then using that to generate a pvalue for the candidate
for (hit in hits){
	samp <- data[which(data$HITNAME==hit),]
	hitvalue <- samp[which(samp$TYPE=="hit"),]
	slope <- hitvalue$SLOPE
	for (num in seq(1,100)){
	numsamp <- samp[which(samp$SAMPNUM==num),]
	matched <- numsamp[which(numsamp$TYPE!="hit"),]
	# matchedgraph <- numsamp[which(numsamp$TYPE=="matched"),]
	permCDF <- ecdf(matched$ABSSLOPE)
	pvalue <- 1-permCDF(hitvalue$ABSSLOPE)

	entry <- c(hit,num,slope,pvalue)
	printpval <- rbind(printpval,entry)

	# plot <- ggplot() + geom_density(data=matchedgraph,aes(x=matchedgraph$ABSSLOPE)) + geom_vline(data=hitvalue,aes(xintercept=hitvalue$ABSSLOPE))
	# ggsave(plot,filename=str_c("densityplot_skylerclinalPval_",hit,".jpg"),device=jpeg)
}
}

colnames(printpval) <- c("HIT","NUM","SLOPE","PVAL")

write.table(printpval, file = str_c("ClinalPvals_Unadjusted_",chrom,".txt"), quote=FALSE)

printpval <- as.data.frame(printpval)

printpvalQuan <- data.frame(HIT=character(),SLOPE=numeric(),PVAL10=numeric(),PVAL50=numeric(),PVAL90=numeric())

for (hit in hits){
	pvals <- printpval[which(printpval$HIT==hit),]
	colnames(pvals) <- c("HIT","NUM","SLOPE","PVAL")
	eightyPercentile <- quantile(as.numeric(pvals$PVAL),probs=c(0.1,0.5,0.9))
	slope <- unique(pvals$SLOPE)
	entry <- c(hit,slope,unname(eightyPercentile[1]),unname(eightyPercentile[2]),unname(eightyPercentile[3]))
	printpvalQuan <- rbind(printpvalQuan,entry)
	

}

colnames(printpvalQuan) <- c("HIT","SLOPE","PVAL10","PVAL50","PVAL90")

write.table(printpvalQuan, file = str_c("ClinalPvals_Quantiles_",chrom,".txt"), quote=FALSE)






