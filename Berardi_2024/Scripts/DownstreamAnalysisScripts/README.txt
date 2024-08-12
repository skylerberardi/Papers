README

Analysis pipeline written by Jess Rhodes

Laid out here is the order & brief description of the scripts in this directory. They should be run in  I wrote all this code prior to learning how to make fully replicable pipelines, so to successfully run this code on your local machine, you will have to go into each script and change the paths to your own files/libraries. If you have any questions or concerns, please reach out to me at jarhodes@stanford.edu.

#############

CLINAL ANALYSIS

1. CandidateSNPs_RegionInversionFreq_ClinalAnalysis.py -- Python script, generates genomic region (exon, intron, intergenic region, etc.), inversion status, and average frequency across clinal samples.

2. AllDESTSNPs_RegionInversionFreq_ClinalAnalysis.py -- Python script, does the same for all SNPs in the DEST dataset

3. DESTSNPs_MatchingAlgorithm_RegionInversionFreq_bootstrapSamps.py (with "Clinal" analysis argument) -- Python script, the same one is used for the seasonal DEST analysis. Generates groups of matched SNPs for each candidate SNP, matched for region, inversion status, and +-0.05 starting frequency.

4. LinearRegressionClinal_MatchedSNPs_bootstrapSamps_efficient.py -- Python script, does linear regression of latitude vs. frequency for every candidate SNP and every matched SNP.

5. P-VALUE CALCULATIONS:
	5a. LinearRegressionClinal_PvalueCalc_bootstrapSamps.R -- R script, generates a p-value for each candidate SNP by creating a null distribution from each of 100 groups of 100 matched SNPs for each candidate.

	5b. LinearRegressionClinal_PvalueCalc_GroupStats.R -- R script, generates a p-value for the whole group of candidate SNPs using a null distribution of group stats of the matched SNPs.

6. ClinalPlots_withErrorBars.R -- R script, generates figures for each of the significant clinal SNPs


#############

SEASONAL ANALYSIS

1. CandidateSNPs_RegionInversionFreq_SeasonAnalysis.py -- Python script, generates genomic region (exon, intron, intergenic region, etc.), inversion status, and average frequency across Linvilla seasonal samples.

2. AllDESTSNPs_RegionInversionFreq_SeasonAnalysis.py -- Python script, does the same for all SNPs in the DEST dataset

3. DESTSNPs_MatchingAlgorithm_RegionInversionFreq_bootstrapSamps.py (with "Season" analysis argument) -- Python script, the same one is used for the clinal DEST analysis. Generates groups of matched SNPs for each candidate SNP, matched for region, inversion status, and +-0.05 starting frequency.

4. MetricCalcSeason_MatchedSNPs_bootstrap.py -- Python script, calculates the total & average directional allele frequency change between spring and fall for all years that have data for both seasons.

5. P-VALUE CALCULATIONS:
	5a. DirectionMetricSeason_PvalueCalc_bootstrapSamps.R -- R script, generates a p-value for each candidate SNP by creating a null distribution from each of 100 groups of 100 matched SNPs for each candidate.

	5b. DirectionMetricSeason_PvalueCalc_GroupStats.R -- R script, generates a p-value for the whole group of candidate SNPs using a null distribution of group stats of the matched SNPs.

6. LinvillaSeason_Plots_withErrorBars.R -- R script, generates figures for each of the significant seasonal SNPs


#############

ORCHARD ANALYSIS

1. Orch16_LinRegMetric.R -- R script, reads in 2016 orchard data and calculates the mean frequency for every SNP at each of the 3 timepoints across all of the cages. Then, it does a linear regression of timepoint vs. averaged frequencies.

2. Orch16_RegionInversionFreq.py --  Python script, generates genomic region (exon, intron, intergenic region, etc.), inversion status, and average frequency across timepoints for every SNP. Then, it generates groups of matched SNPs.

3. Orch16_LinRegFileMerge.py -- Python script, grabs linear regression coefficient for all candidate and matched SNPs.

4. P-VALUE CALCULATIONS:
	4a. Orch16_BootstrapPvals.R -- R script, generates a p-value for each candidate SNP by creating a null distribution from each of 100 groups of 100 matched SNPs for each candidate.

	4b. Orch16_GroupPval.R -- R script, generates a p-value for the whole group of candidate SNPs using a null distribution of group stats of the matched SNPs.

5. Orch16_Plots_withErrorBars.R -- R script, generates figures for each of the significant orchard SNPs


#############

All necessary data can be found in the directory called "SupplementaryFiles"
