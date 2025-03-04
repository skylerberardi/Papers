##### Code used to run correlations between significant pigmentation SNPs in 3 contexts
#Latitudinal, seasonal, and experimental populations

##### Skyler Berardi; contact: berardis@sas.upenn.edu

#Packages
library(tidyverse)
library(reshape2)
library(ggplot2)
library(nlme)
library(readxl)
library(emmeans)
library(plotrix)
library(ggpubr)

### Read dataset with latitudinal/seasonal/experimental allele frequency shift metric 
### for significant SNPs across contexts.
pig_snps <- read_excel("Pig_SNP_correlations.xlsx")

view(pig_snps)

## Subset data for three comparisons, remove rows with n/a

# Latitudinal vs. seasonal (wild populations)

pig_snps_lat.sea <- pig_snps %>% select(-Experimental_slope)
pig_snps_lat.sea <- pig_snps_lat.sea %>% na.omit()
view(pig_snps_lat.sea)

# Latitudinal vs. experimental orchard populations

pig_snps_lat.exp <- pig_snps %>% select(-Seasonal_season.metric)
pig_snps_lat.exp <- pig_snps_lat.exp %>% na.omit()
view(pig_snps_lat.exp)

# Seasonal (wild) vs. experimental orchard populations

pig_snps_sea.exp <- pig_snps %>% select(-Latitudinal_slope)
pig_snps_sea.exp <- pig_snps_sea.exp %>% na.omit()
view(pig_snps_sea.exp)


## Testing correlations and plotting: latitudinal vs. seasonal

lat.sea.cor <- cor.test(pig_snps_lat.sea$Latitudinal_slope, pig_snps_lat.sea$Seasonal_season.metric,
                        method = "pearson")
lat.sea.cor #cor = 0.2950665; p = 0.2502 

#plot
lat.sea.plot <- ggscatter(pig_snps_lat.sea, x="Latitudinal_slope", y="Seasonal_season.metric",
                          add="reg.line", conf.int=TRUE,
                          cor.coef=TRUE, cor.method="pearson",
                          xlab="Allele frequency change from high to low latitudes\n(latitudinal slope)", 
                          ylab="Allele frequency change from early\nto late season (season metric)")
lat.sea.plot


## Testing correlations and plotting: latitudinal vs. experimental

lat.exp.cor <- cor.test(pig_snps_lat.exp$Latitudinal_slope, pig_snps_lat.exp$Experimental_slope,
                        method = "pearson")
lat.exp.cor #cor = -0.1452188; p = 0.6204

#plot
lat.exp.plot <- ggscatter(pig_snps_lat.exp, x="Latitudinal_slope", y="Experimental_slope",
                          add="reg.line", conf.int=TRUE,
                          cor.coef=TRUE, cor.method="pearson",
                          xlab="Allele frequency change from high to low latitudes\n(latitudinal slope)", 
                          ylab="Allele frequency change from early\nto late season (experimental slope)")
lat.exp.plot


## Testing correlations and plotting: seasonal vs. experimental

sea.exp.cor <- cor.test(pig_snps_sea.exp$Seasonal_season.metric, pig_snps_sea.exp$Experimental_slope,
                        method = "pearson")
sea.exp.cor #cor = 0.2128327, p = 0.4651

#plot
sea.exp.plot <- ggscatter(pig_snps_sea.exp, x="Seasonal_season.metric", y="Experimental_slope",
                          add="reg.line", conf.int=TRUE,
                          cor.coef=TRUE, cor.method="pearson",
                          xlab="Allele frequency change from early to late season\n(season metric)", 
                          ylab="Allele frequency change from early\nto late season (experimental slope)")
sea.exp.plot



