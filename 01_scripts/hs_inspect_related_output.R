# Inspect the output of simple_pop_stats relatedness_calc function 
#  with the objective to see if there are any discrete breaks in the relatedness metrics 
#  within populations that can be reasonably used as cutoffs

# Clear space
gc()

# Set working dir
setwd("~/Documents/00_sbio/harbour_seal/ms_harbour_seal")

# Set variables
input.FN <- "../simple_pop_stats/03_results/pairwise_relatedness_output_all_2023-07-12.txt"

# Read in data
rel.df <- read.table(file = input.FN, header = T, sep = "\t")
head(rel.df)


### SOG ###
# How many inds are there? 
all_inds.vec <- c(rel.df$ind1.id, rel.df$ind2.id)
SOG_inds.vec <- all_inds.vec[grep(pattern = "SOG", x = all_inds.vec)]
SOG_inds.vec <- unique(SOG_inds.vec)
length(SOG_inds.vec) # 41 inds

colnames(rel.df)

# Select type
compare.group <- "SOSO"
statistic <- "ritland"
statistic <- "quellergt"
statistic <- "wang"


# Inspect same-on-same distribution of values
print(paste0("Identifying outliers using the ", statistic, " statistic"))

# How many?
length(rel.df[rel.df$group==compare.group, statistic]) # There are 820 combinations where order does not matter and without repetitions

obs_rel.vec <- rel.df[rel.df$group==compare.group, statistic]

# Calculate 
median(obs_rel.vec)
boxplot.stats(obs_rel.vec)$out
# upper outliers
boxplot.stats(obs_rel.vec)$out[boxplot.stats(obs_rel.vec)$out > median(obs_rel.vec)]
num_outliers <- length(boxplot.stats(obs_rel.vec)$out[boxplot.stats(obs_rel.vec)$out > median(obs_rel.vec)])
print(paste0("This relatedness statistic identifies ", num_outliers, " outlier pairs"))
cutoff <- min(boxplot.stats(obs_rel.vec)$out[boxplot.stats(obs_rel.vec)$out > median(obs_rel.vec)])
cutoff

pdf(file = paste0("03_results/related_dist_", compare.group, "_", statistic, ".pdf")
    , width = 9, height = 5)
par(mfrow=c(1,2))
boxplot(obs_rel.vec, las = 1, main = compare.group
        , ylab = statistic)
abline(h = cutoff, lty = 3)

plot(obs_rel.vec, las = 1, main = compare.group
     , ylab = statistic)
abline(h = cutoff, lty = 3)
dev.off()


# # Calculate the margin of error
# #margin <- qt(0.975,df=n-1) * s / sqrt(n) # where s = sd, n = sample size
# # lowerinterval <- xbar - margin # where xbar is the mean
# xbar <- mean(rel.df$ritland[rel.df$group=="SOSO"])
# 
# lowerinterval 
# [1] 195.5191
# upperinterval <- xbar + margin
# upperinterval 
# [1] 204.4809
# 
# margin <- qt(0.975, df=length(obs_rel.vec)-1) * sd(obs_rel.vec) / sqrt(length(obs_rel.vec))
# 
# lowerinterval <- xbar - margin
# upperinterval <- xbar + margin
