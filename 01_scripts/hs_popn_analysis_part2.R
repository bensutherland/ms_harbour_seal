# General analysis of harbour seal data (Part 2)
# 2020-10-13 initialized
# start by running '01_scripts/hs_popn_analysis_part1.R' 
#  or loading the previous image
# Note: to obtain variables such as the individuals to remove based on relatedness, you will need to have interactively
#    run 01_scripts/hs_inspect_related_output.R


# Load previous image
load(file = "03_results/output_coast-sp_relatedness.Rdata")

# Drop selected individuals, Atlantic
inds.to.drop <- c("EQB_103", "EQB_104", "EQB_116", "EQB_117", "EQB_118", "EQB_120"
                  , "NFL_103", "NFL_111", "NFL_114", "NFL_120", "NFL_123", "NFL_119", "NFL_124", "NFL_122"
                  , "LAB_105", "LAB_103"
)

keep.inds <- setdiff(x = indNames(obj_atlantic), y = inds.to.drop)
#obj_atlantic.bck <- obj_atlantic
obj_atlantic <- obj_atlantic[keep.inds]

# Drop selected individuals, Pacific
# TODO


### Add MAF filter since indiv have been dropped

## HWE filter
hwe_eval(data = obj_atlantic, alpha = 0.01)
head(per_locus_hwe_NFL.df)

# Which col contains pval? 
col.oi <- grep(pattern = "Pr", x = colnames(per_locus_hwe_NFL.df))

# Identify hwe outliers based on the dataset
if(dataset=="all"){
  
  # All Atlantic pops 
  # Identify mnames of outliers
  hwe_outlier_mname_EQB.vec     <-   per_locus_hwe_EQB.df[per_locus_hwe_EQB.df[, col.oi] < 0.01, "mname"]
  hwe_outlier_mname_NFL.vec    <-    per_locus_hwe_NFL.df[per_locus_hwe_NFL.df[, col.oi] < 0.01, "mname"]
  hwe_outlier_mname_LAB.vec    <-    per_locus_hwe_LAB.df[per_locus_hwe_LAB.df[, col.oi] < 0.01, "mname"]
  
  # How many outliers (p < 0.01) per population
  print(length(hwe_outlier_mname_EQB.vec))   #  173 markers out of HWE
  print(length(hwe_outlier_mname_NFL.vec))   #  233 markers out of HWE
  print(length(hwe_outlier_mname_LAB.vec))   #   85 markers out of HWE
  
  # How many unique HWE deviating markers?  
  markers_to_drop <- c(hwe_outlier_mname_EQB.vec, hwe_outlier_mname_NFL.vec, hwe_outlier_mname_LAB.vec)
  
}else if(dataset=="balanced"){
  
  # 'Balanced' pops
  # Identify mnames of outliers
  hwe_outlier_mname_EQB.vec     <-   per_locus_hwe_EQB.df[per_locus_hwe_EQB.df[, col.oi] < 0.01, "mname"]
  hwe_outlier_mname_NFL.vec    <-    per_locus_hwe_NFL.df[per_locus_hwe_NFL.df[, col.oi] < 0.01, "mname"]
  
  # How many outliers (p < 0.01) per population
  print(length(hwe_outlier_mname_EQB.vec))   #  173 markers out of HWE
  print(length(hwe_outlier_mname_NFL.vec))   #  233 markers out of HWE
  
  # How many unique HWE deviating markers?  
  markers_to_drop <- c(hwe_outlier_mname_EQB.vec, hwe_outlier_mname_NFL.vec)
  
}

length(markers_to_drop)             # 491 total markers identified
markers_to_drop <- unique(markers_to_drop)
length(markers_to_drop)             #  unique markers identified

markers_to_keep <- setdiff(x = locNames(obj_atlantic), y = markers_to_drop)
length(markers_to_keep) # 3470 markers to keep

obj_atlantic <- obj_atlantic[, loc=markers_to_keep]
obj_atlantic

## Per locus statistics
per_locus_stats(data = obj_atlantic)

# Markers with HOBS > 0.5? 
hobs.outliers <- per_loc_stats.df[per_loc_stats.df$Hobs > 0.5, "mname"] 
length(hobs.outliers) # 266 markers

keep <- setdiff(x = locNames(obj_atlantic), y = hobs.outliers)
length(keep)
obj_atlantic <- obj_atlantic[, loc=keep] 
obj_atlantic

# Re-run per loc stats
per_locus_stats(data = obj_atlantic)

# Save for later
per_loc_stats_atl.df <- per_loc_stats.df


### Save your PCA output into its own folder, then re-run
## Global PCA
pca_from_genind(data = obj_atlantic
                , PCs_ret = 4
                , plot_eigen = TRUE
                , plot_allele_loadings = TRUE
                , retain_pca_obj = TRUE
                , colour_file = "00_archive/harbour_seal_pops_colours.csv"
)

# Manually run pca_from_genind to pull out the pca1 obj (note: should assign this out to the global enviro as default)
pca_scores_result <- pca.obj$scores
write.csv(x = pca_scores_result, file = "03_results/pca_scores_result_atlantic_relatives_rem.csv", quote = F, row.names = T)

# Rename so the PCA figures are not overwritten
file.copy(from = "03_results/pca_samples_PC1_v_PC2.pdf", to = "03_results/pca_samples_PC1_v_PC2_atl_relatives_rem.pdf", overwrite = TRUE)
file.copy(from = "03_results/pca_samples_PC3_v_PC4.pdf", to = "03_results/pca_samples_PC3_v_PC4_atl_relatives_rem.pdf", overwrite = TRUE)



## FST
calculate_FST(format = "genind", dat = obj_atlantic, separated = FALSE, bootstrap = TRUE)

# Move all results into an 'Atlantic' folder, then proceed to Pacific analysis



# ## HWE filter
# hwe_eval(data = obj_pacific, alpha = 0.01)
# head(per_locus_hwe_SOG.df)
# 
# # Identify column with the p-val
# col.oi <- grep(pattern = "Pr", x = colnames(per_locus_hwe_SOG.df))
# 
# # Identify hwe outliers based on the dataset
# if(dataset=="all"){
#   
#   # All Pacific pops 
#   # Identify mnames of outliers
#   hwe_outlier_mname_SOG.vec     <-   per_locus_hwe_SOG.df[per_locus_hwe_SOG.df[, col.oi] < 0.01, "mname"]
#   hwe_outlier_mname_NBC.vec    <-    per_locus_hwe_NBC.df[per_locus_hwe_NBC.df[, col.oi] < 0.01, "mname"]
#   hwe_outlier_mname_CAL.vec    <-    per_locus_hwe_CAL.df[per_locus_hwe_CAL.df[, col.oi] < 0.01, "mname"]
#   hwe_outlier_mname_ORE.vec    <-    per_locus_hwe_ORE.df[per_locus_hwe_ORE.df[, col.oi] < 0.01, "mname"]
#   
#   # How many outliers (p < 0.01) per population
#   print(length(hwe_outlier_mname_SOG.vec))   #  513 markers out of HWE
#   print(length(hwe_outlier_mname_NBC.vec))   #  155 markers out of HWE
#   print(length(hwe_outlier_mname_CAL.vec))   #  323 markers out of HWE
#   print(length(hwe_outlier_mname_ORE.vec))   #  306 markers out of HWE
#   
#   # How many unique HWE deviating markers?  
#   markers_to_drop <- c(hwe_outlier_mname_SOG.vec, hwe_outlier_mname_NBC.vec, hwe_outlier_mname_CAL.vec, hwe_outlier_mname_ORE.vec)
#   
# }else if(dataset=="balanced"){
#   
#   # 'Balanced' pops
#   # Identify mnames of outliers
#   hwe_outlier_mname_SOG.vec     <-   per_locus_hwe_SOG.df[per_locus_hwe_SOG.df[, col.oi] < 0.01, "mname"]
#   hwe_outlier_mname_ORE.vec    <-    per_locus_hwe_ORE.df[per_locus_hwe_ORE.df[, col.oi] < 0.01, "mname"]
#   
#   # How many outliers (p < 0.01) per population
#   print(length(hwe_outlier_mname_SOG.vec))   #  513 markers out of HWE
#   print(length(hwe_outlier_mname_ORE.vec))   #  306 markers out of HWE
#   
#   # How many unique HWE deviating markers?  
#   markers_to_drop <- c(hwe_outlier_mname_SOG.vec, hwe_outlier_mname_ORE.vec)
#   
# }
# 
# 
# # How many unique HWE deviating markers?  
# length(markers_to_drop)             # 1297 markers out of HWE in at least one population
# markers_to_drop <- unique(markers_to_drop)
# length(markers_to_drop)             #  996 unique markers out of HWE in at least one population
# markers_to_keep <- setdiff(x = locNames(obj_pacific), y = markers_to_drop)
# length(markers_to_keep) # 7937 markers to keep
# 
# obj_pacific <- obj_pacific[, loc=markers_to_keep]
# obj_pacific

## Per locus statistics
per_locus_stats(data = obj_pacific)

# Markers with HOBS > 0.5? 
hobs.outliers <- per_loc_stats.df[per_loc_stats.df$Hobs > 0.5, "mname"] 
length(hobs.outliers) # 230 markers

keep <- setdiff(x = locNames(obj_pacific), y = hobs.outliers)
length(keep)
obj_pacific <- obj_pacific[, loc=keep] 
obj_pacific

# Re-run per loc stats
per_locus_stats(data = obj_pacific)

# Save for later
per_loc_stats_pac.df <- per_loc_stats.df


## Global PCA
pca_from_genind(data = obj_pacific
                , PCs_ret = 4
                , plot_eigen = TRUE
                , plot_allele_loadings = TRUE
                , retain_pca_obj = TRUE
                , colour_file = "00_archive/harbour_seal_pops_colours.csv"
)

# Manually run pca_from_genind to pull out the pca1 obj (note: should assign this out to the global enviro as default)
pca_scores_result <- pca.obj$scores
write.csv(x = pca_scores_result, file = "03_results/pca_scores_result.csv", quote = F, row.names = T)

## FST
calculate_FST(format = "genind", dat = obj_pacific, separated = FALSE, bootstrap = TRUE)


save.image(file = "03_results/output_coast-sp_pruned_analysis.Rdata")

# Move to next script, hs_popn_analysis_part3.R

