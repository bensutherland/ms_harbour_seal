# General analysis of harbour seal data (Part 2)
# 2020-10-13 initialized
#  Requires that '01_scripts/hs_popn_analysis_part1.R' was run, and any blacklist indiv were ID'd from '01_scripts/hs_inspect_related_output.R'

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

# Start by sourcing simple_pop_stats, then load the previous image
load(file = "03_results/output_coast-sp_relatedness.Rdata")
rm(date)

#### 00. Inspect private alleles between Burrard outlier cluster and the rest of the SOG samples
## Private Alleles, Burrard Inlet
# Explore Burrard Inlet samples further to see if any private alleles are present
indNames(obj_pacific)
outlier_samples.id <- c("SOG_108", "SOG_131", "SOG_135", "SOG_115", "SOG_127", "SOG_141", "SOG_125", "SOG_136", "SOG_102")
nonoutlier_samples.id <- c("SOG_134", "SOG_114", " SOG_123", "SOG_124", "SOG_130", "SOG_117", "SOG_119", "SOG_129", "SOG_101") # haphazardly selected inds from non-Burrard group

# Keep only SOG
obj_pacific.sep <- seppop(obj_pacific)
obj.SOG <- obj_pacific.sep$SOG

# Create a df defining the strata for each individual
strata.df <- as.data.frame(as.character(pop(obj.SOG)), stringsAsFactors = F)
colnames(strata.df)[1] <- "indiv.pop"
head(strata.df)

# Add column 'repunit'
strata.df$repunit <- strata.df$indiv.pop

strata.df[which(indNames(obj.SOG) %in% outlier_samples.id), "repunit"] <- "Burrard"
#strata.df[which(indNames(obj.SOG) %in% nonoutlier_samples.id), "repunit"] <- "Burrard" # haphazardly selected, for comparison

# Confirm selection method works:
indNames(obj.SOG)[which(indNames(obj.SOG) %in% outlier_samples.id)]

# Add strat object to genind
strata(obj.SOG) <- strata.df
obj.SOG
table(strata(obj.SOG))

# Confirms the correct inds were labeled by the strata
cbind(indNames(obj.SOG), strata(obj.SOG)) # compare with selected individual names above (i.e., outlier_samples.id or nonoutlier_samples.id)

per_repunit.privallele <- private_alleles(obj.SOG, alleles ~ repunit)

#per_repunit.privallele[, 1:5]
dim(per_repunit.privallele)

table(per_repunit.privallele["SOG", ])
table(per_repunit.privallele["Burrard", ])

which(per_repunit.privallele["Burrard", ]==9)
per_repunit.privallele[,1991]
myFreq.pac[which(names(myFreq.pac)=="586731_19")] # What was the MAF of this variant? 0.046875
obj_pacific # how many inds?
# 97 * 2 = 194; 9 / 194 = 0.04639175 (remember, this assumes no missing data); so this is actually the number of views, not the number of inds
obj_pacific[loc="586731_19"]$tab # observe the counts of the geno

which(per_repunit.privallele["Burrard", ]==6)
which(per_repunit.privallele["Burrard", ]==5)
which(per_repunit.privallele["Burrard", ]==4)


# Get raw number of private alleles per locus
#pal <- private_alleles(obj.SOG, locus ~ repunit, count.alleles = FALSE) # Note: runs slow
#table(pal) # this only gives a 0 or 1, does not count alleles, allows one to see exact how many private alleles exist per repunit
# This shows the number of rows that have a private allele
#rowSums(pal) # n = 2289 in SOG, 190 in Burrard

save.image(file = "03_results/private_allele_analysis.Rdata")


#### 01. Remove relatedness outliers ####
# Specify individuals to drop
inds.to.drop <- c(  "EQB_103", "EQB_104", "EQB_116", "EQB_117", "EQB_118", "EQB_120"
                  , "NFL_103", "NFL_111", "NFL_114", "NFL_120", "NFL_123", "NFL_119", "NFL_124", "NFL_122"
                  , "LAB_105", "LAB_103"
                  , "CAL_103", "CAL_113", "CAL_114", "CAL_117", "CAL_122", "CAL_126", "CAL_121", "CAL_123"
                  , "NBC_112", "NBC_113"
                  , "ORE_119", "ORE_120"
                  , "SOG_137", "SOG_108", "SOG_115", "SOG_136", "SOG_141", "SOG_128", "SOG_130", "SOG_129", "SOG_132", "SOG_131", "SOG_135", "SOG_138"
)

# Drop indiv from each obj
keep_inds_atl <- setdiff(x = indNames(obj_atlantic), y = inds.to.drop)
obj_atlantic <- obj_atlantic[keep_inds_atl]
obj_atlantic

keep_inds_pac <- setdiff(x = indNames(obj_pacific), y = inds.to.drop)
obj_pacific  <- obj_pacific[keep_inds_pac]
obj_pacific


#### 02. Filter new dataset per locus ####
## MAF 
maf_filt(data = obj_atlantic, maf = 0.01)
obj_atlantic <- obj_maf_filt
myFreq.atl <- myFreq

maf_filt(data = obj_pacific, maf = 0.01)
obj_pacific <- obj_maf_filt
myFreq.pac <- myFreq


## HWE
hwe_eval(data = obj_atlantic, alpha = 0.01)
hwe_eval(data = obj_pacific, alpha = 0.01)

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
print(length(hwe_outlier_mname_EQB.vec))   #  129 markers out of HWE
print(length(hwe_outlier_mname_NFL.vec))   #  136 markers out of HWE
print(length(hwe_outlier_mname_LAB.vec))   #    0 markers out of HWE

print(length(hwe_outlier_mname_CAL.vec))   #  282 markers out of HWE
print(length(hwe_outlier_mname_NBC.vec))   #    0 markers out of HWE
print(length(hwe_outlier_mname_ORE.vec))   #  291 markers out of HWE
print(length(hwe_outlier_mname_SOG.vec))   #  369 markers out of HWE

# How many total per region, and how many unique per region? 
markers_to_drop_atl <- c(hwe_outlier_mname_EQB.vec, hwe_outlier_mname_NFL.vec, hwe_outlier_mname_LAB.vec)
print(paste0("Identified ", length(markers_to_drop_atl), " markers out of HWE in at least one population (incl. redundant names)"))
markers_to_drop_atl <- unique(markers_to_drop_atl)
print(paste0("Of these, there are ", length(markers_to_drop_atl), " unique markers"))

markers_to_drop_pac <- c(hwe_outlier_mname_CAL.vec, hwe_outlier_mname_NBC.vec, hwe_outlier_mname_ORE.vec, hwe_outlier_mname_SOG.vec)
print(paste0("Identified ", length(markers_to_drop_pac), " markers out of HWE in at least one population (incl. redundant names)"))
markers_to_drop_pac <- unique(markers_to_drop_pac)
print(paste0("Of these, there are ", length(markers_to_drop_pac), " unique markers"))

## Remove HWE outliers
markers_to_keep_atl <- setdiff(x = locNames(obj_atlantic), y = markers_to_drop_atl)
length(markers_to_keep_atl) # 3454 markers to keep

obj_atlantic <- obj_atlantic[, loc=markers_to_keep_atl]
obj_atlantic

markers_to_keep_pac <- setdiff(x = locNames(obj_pacific), y = markers_to_drop_pac)
length(markers_to_keep_pac) # 7979 markers to keep

obj_pacific <- obj_pacific[, loc=markers_to_keep_pac]
obj_pacific


## Per locus statistics
rm(date)
per_locus_stats(data = obj_atlantic)
date <- format(Sys.time(), "%Y-%m-%d")
file.copy(from = paste0("03_results/per_locus_stats_", date, ".txt"), to = paste0("03_results/per_locus_stats_atl_", date, ".txt"), overwrite = T)
per_loc_stats_atl.df <- per_loc_stats.df

rm(date)
per_locus_stats(data = obj_pacific)
date <- format(Sys.time(), "%Y-%m-%d")
file.copy(from = paste0("03_results/per_locus_stats_", date, ".txt"), to = paste0("03_results/per_locus_stats_pac_", date, ".txt"), overwrite = T)
per_loc_stats_pac.df <- per_loc_stats.df

# Markers with HOBS > 0.5? 
hobs.outliers_atl <- per_loc_stats_atl.df[per_loc_stats_atl.df$Hobs > 0.5, "mname"] 
length(hobs.outliers_atl) # 320 markers

keep_atl <- setdiff(x = locNames(obj_atlantic), y = hobs.outliers_atl)
length(keep_atl)
obj_atlantic <- obj_atlantic[, loc=keep_atl] 
obj_atlantic

hobs.outliers_pac <- per_loc_stats_pac.df[per_loc_stats_pac.df$Hobs > 0.5, "mname"] 
length(hobs.outliers_pac) # 280 markers

keep_pac <- setdiff(x = locNames(obj_pacific), y = hobs.outliers_pac)
length(keep_pac)
obj_pacific <- obj_pacific[, loc=keep_pac] 
obj_pacific


## Coast-specific PCA
pca_from_genind(data = obj_atlantic
                , PCs_ret = 4
                , plot_eigen = TRUE
                , plot_allele_loadings = FALSE
                , retain_pca_obj = TRUE
                , colour_file = "00_archive/harbour_seal_pops_colours.csv"
)

atl_pc1_v_pc2.plot <- pc1_v_pc2.plot # Use this one for the future
atl_pc3_v_pc4.plot <- pc3_v_pc4.plot


# Rename so the PCA figures are not overwritten
file.copy(from = "03_results/pca_samples_PC1_v_PC2.pdf", to = "03_results/pca_samples_PC1_v_PC2_atl_pruned.pdf", overwrite = TRUE)
file.copy(from = "03_results/pca_samples_PC3_v_PC4.pdf", to = "03_results/pca_samples_PC3_v_PC4_atl_pruned.pdf", overwrite = TRUE)

# Retain and save out the PCA scores
pca_scores_result <- pca.obj$scores
write.csv(x = pca_scores_result, file = "03_results/pca_scores_result_atlantic_pruned.csv", quote = F, row.names = T)

# Pacific
pca_from_genind(data = obj_pacific
                , PCs_ret = 4
                , plot_eigen = TRUE
                , plot_allele_loadings = FALSE
                , retain_pca_obj = TRUE
                , colour_file = "00_archive/harbour_seal_pops_colours.csv"
                , outline_samples = TRUE
)

pac_pc1_v_pc2.plot <- pc1_v_pc2.plot # Use this one for the future
pac_pc3_v_pc4.plot <- pc3_v_pc4.plot


# Rename so the PCA figures are not overwritten
file.copy(from = "03_results/pca_samples_PC1_v_PC2.pdf", to = "03_results/pca_samples_PC1_v_PC2_pac_pruned.pdf", overwrite = TRUE)
file.copy(from = "03_results/pca_samples_PC3_v_PC4.pdf", to = "03_results/pca_samples_PC3_v_PC4_pac_pruned.pdf", overwrite = TRUE)

# Retain and save out the PCA scores
pca_scores_result <- pca.obj$scores
write.csv(x = pca_scores_result, file = "03_results/pca_scores_result_pacific_pruned.csv", quote = F, row.names = T)



## FST
calculate_FST(format = "genind", dat = obj_atlantic, separated = FALSE, bootstrap = TRUE)
file.copy(from = "03_results/gen_diff_wcfst_booted.csv", to = "03_results/gen_diff_wcfst_booted_atl.csv", overwrite = T)
calculate_FST(format = "genind", dat = obj_pacific, separated = FALSE, bootstrap = TRUE)
file.copy(from = "03_results/gen_diff_wcfst_booted.csv", to = "03_results/gen_diff_wcfst_booted_pac.csv", overwrite = T)

save.image(file = "03_results/output_coast-sp_pruned_analysis.Rdata")

# Move to next script, hs_popn_analysis_part3.R
