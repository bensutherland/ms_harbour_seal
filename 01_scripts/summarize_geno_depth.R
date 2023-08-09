# Obtain stats on depth of coverage per locus from the output VCF
#  Sutherland Bioinformatics
#  2023-08-01

#install.packages("vcfR")
library("vcfR")

# Set working directory
setwd("~/Documents/00_sbio/harbour_seal/simple_pop_stats_v.0.2/")

# Set input filename
input.FN <- "../stacks_workflow_2023-02-27/05-stacks/popn_out_single_snp/bhs_p7_r0.7_maf0.01_2023-02-28.vcf"

# Load data
vcf <- read.vcfR(file = input.FN)

# Extract genotype depth per sample and locus
dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)
dp[1:4, 1:6]

# Overall stats
mean(x = dp, na.rm = T) # 35.1x coverage
sd(x = dp, na.rm = T)   # +/- 273.9x
min(dp, na.rm = T)
max(dp, na.rm = T)

# Per individual averages
mean_cov.samples <- colMeans(x = dp, na.rm = T)
mean_cov.samples <- as.data.frame(mean_cov.samples)
mean_cov.samples$indiv <- rownames(mean_cov.samples)
mean_cov.samples <- mean_cov.samples[,c("indiv", "mean_cov.samples")]
head(mean_cov.samples)

# Add population variable
mean_cov.samples$pop <- gsub(pattern = "_.*", replacement = "", x = mean_cov.samples$indiv)

# Summarize average by group
aggregate(x = mean_cov.samples$mean_cov.samples, list(mean_cov.samples$pop), FUN=mean)

# Population averages
# Plot per pop
pdf(file = "03_results/mean_geno_depth_per_pop.pdf", width = 6.5, height = 4.5)
boxplot(mean_cov.samples$mean_cov.samples ~ mean_cov.samples$pop
        , las = 1
        , xlab = "Populations"
        , ylab = "Mean genotype depth, per individual"
        )
dev.off()

# Individual plot
pdf(file = "03_results/mean_geno_depth_per_indiv.pdf", width = 19, height = 6)
barplot(height = mean_cov.samples$mean_cov.samples
        , axisnames = TRUE, names.arg = mean_cov.samples$indiv
        , las = 3
        , cex.names = 0.7
        , ylab = "Mean genotype depth per indiv"
        )

#text(x = 160, y = 120, paste0("n = ", nrow(dp), " loci"))
dev.off()
