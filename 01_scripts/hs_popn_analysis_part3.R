# General analysis of harbour seal data (Part 3)
# 2020-10-13 initialized
#  Requires that '01_scripts/hs_popn_analysis_part2.R' was run

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

# Start by sourcing simple_pop_stats, then load the previous image
load(file = "03_results/output_coast-sp_pruned_analysis.Rdata")
rm(date) # Clear earlier date variable


#### 02. Identify top FST markers ####
# Calculate FST for SOG vs. NBC and SOG vs. ORE
obj_pacific_filt.sep <- seppop(x = obj_pacific) 

### SOG vs. NBC
obj_pacific_filt <- repool(obj_pacific_filt.sep$NBC, obj_pacific_filt.sep$SOG)

## MAF
maf_filt(data = obj_pacific_filt, maf = 0.01)
obj_pacific_filt <- obj_maf_filt
myFreq_filt_SOG_vs_NBC <- myFreq

# Calculate per-locus FST
per_locus_stats(data = obj_pacific_filt)
date <- format(Sys.time(), "%Y-%m-%d")
file.copy(from = paste0("03_results/per_locus_stats_", date, ".txt"), to = paste0("03_results/per_locus_stats_SOG_vs_NBC_", date, ".txt"), overwrite = T)
per_loc_stats_SOG_vs_NBC.df <- per_loc_stats.df

### SOG vs. ORE
obj_pacific_filt <- repool(obj_pacific_filt.sep$ORE, obj_pacific_filt.sep$SOG)

## MAF 
maf_filt(data = obj_pacific_filt, maf = 0.01)
obj_pacific_filt <- obj_maf_filt
myFreq_filt_SOG_vs_ORE <- myFreq

# Calculate per-locus FST
per_locus_stats(data = obj_pacific_filt)
date <- format(Sys.time(), "%Y-%m-%d")
file.copy(from = paste0("03_results/per_locus_stats_", date, ".txt"), to = paste0("03_results/per_locus_stats_SOG_vs_ORE_", date, ".txt"), overwrite = T)
per_loc_stats_SOG_vs_ORE.df <- per_loc_stats.df
rm(date)


#### 03. Identify top HOBS markers #### 
# Calculate per-locus HOBS for BC or US pops separately
names(obj_pacific_filt.sep)

### BC
obj <- repool(obj_pacific_filt.sep$SOG, obj_pacific_filt.sep$NBC)
maf_filt(data = obj, maf = 0.01)
obj <- obj_maf_filt
per_locus_stats(data = obj)
per_loc_stats_BC.df <- per_loc_stats.df

### US
obj <- repool(obj_pacific_filt.sep$ORE, obj_pacific_filt.sep$CAL)
maf_filt(data = obj, maf = 0.01)
obj <- obj_maf_filt
per_locus_stats(data = obj)
per_loc_stats_US.df <- per_loc_stats.df


#### 04. Combine population or region-sp data ####
head(per_loc_stats_SOG_vs_ORE.df)
head(per_loc_stats_SOG_vs_NBC.df)

head(per_loc_stats_BC.df)
head(per_loc_stats_US.df)

## Rename FST column before merge for downstream plotting
colnames(per_loc_stats_SOG_vs_ORE.df)[which(colnames(per_loc_stats_SOG_vs_ORE.df)=="Fst")] <- "Fst.SOG.ORE"
colnames(per_loc_stats_SOG_vs_NBC.df)[which(colnames(per_loc_stats_SOG_vs_NBC.df)=="Fst")] <- "Fst.SOG.NBC"

colnames(per_loc_stats_BC.df)[which(colnames(per_loc_stats_BC.df)=="Hobs")] <- "Hobs.BC"
colnames(per_loc_stats_US.df)[which(colnames(per_loc_stats_US.df)=="Hobs")] <- "Hobs.US"

## Merge
SOG_comparisons_FST.df <- merge(x = per_loc_stats_SOG_vs_ORE.df, per_loc_stats_SOG_vs_NBC.df, by = "mname")

BC_and_US_sp_HOBS.df <- merge(x = per_loc_stats_BC.df, y = per_loc_stats_US.df, by = "mname")

## Plot
pdf(file = "03_results/FST_and_HOBS_region-specific.pdf", width = 9, height = 5)
par(mfrow=c(1,2))
plot(x = SOG_comparisons_FST.df$Fst.SOG.NBC, y = SOG_comparisons_FST.df$Fst.SOG.ORE
     , xlab = expression(italic(F)[ST] ~ SOG ~ vs. ~ NBC)
     , ylab = expression(italic(F)[ST] ~ SOG ~ vs. ~ ORE)
     , las = 1
)
text(x = 0.6, y = 0.5, labels = paste0("n = ", nrow(SOG_comparisons_FST.df)))

plot(x = BC_and_US_sp_HOBS.df$Hobs.BC, y = BC_and_US_sp_HOBS.df$Hobs.US
     , xlab = expression(per ~ locus ~ H[OBS] ~ (SOG ~ and ~ NBC))
     , ylab = expression(per ~ locus ~ H[OBS] ~ (ORE ~ and ~ CAL))
     , las = 1
)

dev.off()


### 06. Multi-coast plotting ####
## MAF
# Remove the HWE-filtered variants
myFreq.pac <- myFreq.pac[names(myFreq.pac) %in% locNames(obj_pacific)]
myFreq.atl <- myFreq.atl[names(myFreq.atl) %in% locNames(obj_atlantic)]

table(myFreq.pac < 0.1)
table(myFreq.atl < 0.1)

length(myFreq.pac)
length(myFreq.atl)

table(myFreq.pac < 0.1)[2] / length(myFreq.pac) # 56.3%
table(myFreq.atl < 0.1)[2] / length(myFreq.atl) # 46.5%

table(myFreq.pac < 0.05)
table(myFreq.atl < 0.05)

table(myFreq.pac < 0.05)[2] / length(myFreq.pac) # 36.8%
table(myFreq.atl < 0.05)[2] / length(myFreq.atl) # 29.4%


# Plot
pdf(file = paste0("03_results/MAF_hist_pac_atl.pdf"), width = 7, height = 4)
par(mfrow=c(1,2))
hist(myFreq.pac
     #, proba=T # note: does not sum to 1, not worth using
     , col="grey", xlab = "MAF, Pacific"
     , main = ""
     #, ylim = c(0, 1250)
     , ylab = "Number of loci"
     , las = 1
     , breaks = 20
)
text(x = 0.4, y = 1000, labels = paste("n = ", length(myFreq.pac), " loci", sep = "" ))

hist(myFreq.atl
     #, proba=T # note: does not sum to 1, not worth using
     , col="grey", xlab = "MAF, Atlantic"
     , main = ""
     #, ylim = c(0, 1250)
     , ylab = "Number of loci"
     , las = 1
     , breaks = 20
)
text(x = 0.4, y = 250, labels = paste("n = ", length(myFreq.atl), " loci", sep = "" ))
dev.off()

# Save out the MAF calculation as a table
myFreq.pac <- round(myFreq.pac, digits = 3)
write.table(x = myFreq.pac, file = "03_results/allele_freq_retained_loci_pac.txt"
            , sep = "\t", quote = F
            , row.names = T, col.names = F
)

myFreq.atl <- round(myFreq.atl, digits = 3)
write.table(x = myFreq.atl, file = "03_results/allele_freq_retained_loci_atl.txt"
            , sep = "\t", quote = F
            , row.names = T, col.names = F
)


##### 07. Final population analyses ####
# Go back to original genepop
obj

## Remove related individuals (previously identified)
length(inds.to.drop)
keep_inds_all <- setdiff(x = indNames(obj), y = inds.to.drop)
obj_all <- obj[keep_inds_all]
obj_all

## MAF 
maf_filt(data = obj_all, maf = 0.01)
obj_all <- obj_maf_filt

## HWE
hwe_eval(data = obj_all, alpha = 0.01)

# Which col contains the pval? 
head(per_locus_hwe_NFL.df)
col.oi <- grep(pattern = "Pr", x = colnames(per_locus_hwe_NFL.df))

# Identify hwe outliers
hwe_outlier_mname_EQB.vec     <-   per_locus_hwe_EQB.df[per_locus_hwe_EQB.df[, col.oi] < 0.01, "mname"]
hwe_outlier_mname_NFL.vec    <-    per_locus_hwe_NFL.df[per_locus_hwe_NFL.df[, col.oi] < 0.01, "mname"]
hwe_outlier_mname_LAB.vec    <-    per_locus_hwe_LAB.df[per_locus_hwe_LAB.df[, col.oi] < 0.01, "mname"]
hwe_outlier_mname_CAL.vec    <-    per_locus_hwe_CAL.df[per_locus_hwe_CAL.df[, col.oi] < 0.01, "mname"]
hwe_outlier_mname_NBC.vec    <-    per_locus_hwe_NBC.df[per_locus_hwe_NBC.df[, col.oi] < 0.01, "mname"]
hwe_outlier_mname_ORE.vec    <-    per_locus_hwe_ORE.df[per_locus_hwe_ORE.df[, col.oi] < 0.01, "mname"]
hwe_outlier_mname_SOG.vec     <-   per_locus_hwe_SOG.df[per_locus_hwe_SOG.df[, col.oi] < 0.01, "mname"]

# Reporting: how many outliers (p < 0.01) per population
print(length(hwe_outlier_mname_EQB.vec))   # 123   markers out of HWE
print(length(hwe_outlier_mname_NFL.vec))   # 131   markers out of HWE
print(length(hwe_outlier_mname_LAB.vec))   #   0   markers out of HWE
print(length(hwe_outlier_mname_CAL.vec))   # 281   markers out of HWE
print(length(hwe_outlier_mname_NBC.vec))   #   0   markers out of HWE
print(length(hwe_outlier_mname_ORE.vec))   # 289   markers out of HWE
print(length(hwe_outlier_mname_SOG.vec))   # 360   markers out of HWE

# How many total, and how many unique? 
markers_to_drop_all <- c(hwe_outlier_mname_EQB.vec, hwe_outlier_mname_NFL.vec, hwe_outlier_mname_LAB.vec
                         , hwe_outlier_mname_CAL.vec, hwe_outlier_mname_NBC.vec, hwe_outlier_mname_ORE.vec, hwe_outlier_mname_SOG.vec
                         )
print(paste0("Identified ", length(markers_to_drop_all), " markers out of HWE in at least one population (incl. redundant names)"))
markers_to_drop_all <- unique(markers_to_drop_all)
print(paste0("Of these, there are ", length(markers_to_drop_all), " unique markers"))

## Remove HWE outliers
markers_to_keep_all <- setdiff(x = locNames(obj_all), y = markers_to_drop_all)
length(markers_to_keep_all) # 9152 markers to keep

obj_all <- obj_all[, loc=markers_to_keep_all]
obj_all


## PCA, all pops
pca_from_genind(data = obj_all
                , PCs_ret = 4
                , plot_eigen = TRUE
                , plot_allele_loadings = TRUE
                , retain_pca_obj = TRUE
                , colour_file = "00_archive/harbour_seal_pops_colours.csv"
)

# Manually run pca_from_genind to pull out the pca1 obj (note: should assign this out to the global enviro as default)
pca_scores_result <- pca.obj$scores
write.csv(x = pca_scores_result, file = "03_results/pca_scores_result_all_pruned.csv", quote = F, row.names = T)

# Rename so the PCA figures are not overwritten
file.copy(from = "03_results/pca_samples_PC1_v_PC2.pdf", to = "03_results/pca_samples_PC1_v_PC2_all_pruned.pdf", overwrite = TRUE)
file.copy(from = "03_results/pca_samples_PC3_v_PC4.pdf", to = "03_results/pca_samples_PC3_v_PC4_all_pruned.pdf", overwrite = TRUE)

## FST
calculate_FST(format = "genind", dat = obj_all, separated = FALSE, bootstrap = TRUE)
file.copy(from = "03_results/gen_diff_wcfst_booted.csv", to = "03_results/gen_diff_wcfst_booted_all_pruned.csv", overwrite = TRUE)

## Dendrogram
make_tree(boot_obj = obj_all, bootstrap = TRUE, nboots = 10000, dist_metric = "edwards.dist", separated = FALSE)

#### 0.4 Export ####
# Write out object
save.image(file = "03_results/completed_analysis.RData")
