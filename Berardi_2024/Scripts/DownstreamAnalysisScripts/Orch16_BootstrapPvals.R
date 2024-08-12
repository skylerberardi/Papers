.libPaths("/home/users/jarhodes/R/x86_64-pc-linux-gnu-library/4.2")
library(tidyverse)
library(reshape2)
library(coxed)
setwd("/scratch/groups/dpetrov/jess/pigmentselect")

args = commandArgs(trailingOnly=TRUE)
chrom = args[1]

printpval <- data.frame(HIT=character(),BETA=numeric(),PVAL10=numeric(),PVAL50=numeric(),PVAL90=numeric())

data <- read.table(str_c("MatchedSNPs_100Groups_LinRegBetas_Orch16_",chrom,".txt"),header=TRUE)
data$HITCHROMPOS <- str_c(data$HITCHROM,"_",data$HITPOS)
data$BETAABS <- abs(data$BETA)

hitdata <- data[which(data$TYPE=="hit"),]
hitsnames <- unique(hitdata$HITCHROMPOS)

matches <- data[which(data$TYPE=="matched"),]

for (hit in hitsnames){
	subsetmatch <- matches[which(matches$HITCHROMPOS==hit),]
	hitinfo <- hitdata[which(hitdata$HITCHROMPOS==hit),]
	beta <- unique(hitinfo$BETAABS)
	betasign <- unique(hitinfo$BETA)
	for (i in seq(1,100)){
		group <- subsetmatch[which(subsetmatch$SAMP==i),]
		groupCDF <- ecdf(group$BETAABS)
		grouppval <- 1-groupCDF(beta)
		if (i==1){
			betalist <- c(grouppval)
		} else{
			betalist <- append(betalist,grouppval)
		}
	}
	print(betalist)
	eightyPercentile <- quantile(as.numeric(betalist),probs=c(0.1,0.5,0.9))
	entry <- c(hit,betasign,unname(eightyPercentile[1]),unname(eightyPercentile[2]),unname(eightyPercentile[3]))
	printpval <- rbind(printpval,entry)

}

printpval <- as.data.frame(printpval)
colnames(printpval) <- c("HIT","BETA","PVAL10","PVAL50","PVAL90")

write.table(printpval,file = str_c("Cage16Pvals_Quantiles_",chrom,".txt"), quote=FALSE)