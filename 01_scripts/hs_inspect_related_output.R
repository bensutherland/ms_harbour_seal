# Inspect the output of simple_pop_stats relatedness_calc function 
#  with the objective to see if there are any discrete breaks in the relatedness metrics 
#  within populations that can be reasonably used as cutoffs

# Run after running sections of hs_popn_analysis (and therefore are in active simple_pop_stats repo)
#  if starting again from here, need to launch simple_pop_stats

# Clear space
gc()

# Set variables
input.FN <- "03_results/pairwise_relatedness_output_all_2023-07-15.txt" # Atlantic analysis

# Read in data
rel.df <- read.table(file = input.FN, header = T, sep = "\t")
head(rel.df)

# These are the contrasts in this data
unique(rel.df$group)

# Set variables of interest
#popn <- "EQB"
#popn <- "NFL"
popn <- "LAB"

compare.group <- paste0(substr(x = popn, 1,2), substr(x = popn, 1,2))

# How many unique inds are there in the selected pop? 
all_inds.vec  <- c(rel.df$ind1.id, rel.df$ind2.id)
uniq_inds.vec <- unique(all_inds.vec[grep(pattern = popn, x = all_inds.vec)])
length(uniq_inds.vec)

colnames(rel.df)

# Select related stat
statistic <- "ritland"
# statistic <- "quellergt"
# statistic <- "wang"

# Inspect same-on-same distribution of values
print(paste0("Identifying outliers using the ", statistic, " statistic"))

# How many?
length(rel.df[rel.df$group==compare.group, statistic])

obs_rel.vec <- rel.df[rel.df$group==compare.group, statistic]

# Calculate 
median(obs_rel.vec)
boxplot.stats(obs_rel.vec)$out

# Upper outliers
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

# Then use the cutoff value to inspect the excel document to ID the pairs, removing one of every two pairs until no outlier pairs remain.
