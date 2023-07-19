# General analysis of harbour seal data
# 2020-10-13 initialized
# start by sourcing simple_pop_stats then quit the menu

# Set flag, using all data or 'balanced' (limited)?  
dataset <- "all"
#dataset <- "balanced" # only has two populations per region

#### 01. Load data and check missing per ind ####
# Load genepop and characterize
load_genepop(datatype = "SNP")

## file sources: 
# all populations, here:  "02_input_data/bhs_p7_r0.7_maf0.01_2023-02-28.gen"
# balanced and normalized here: "02_input_data/bhs_p4_r0.7_maf0.01_2023-03-04.gen" 

# Clean up pop names
pop(obj) <- gsub(pattern = "_.*", replacement = "", x = pop(obj))
unique(pop(obj))

# Characterize genepop
characterize_genepop(df = obj, pdf_width = 7, pdf_height = 5, cex_names = 0.8, N=30)

# Characterize missing data per individual
percent_missing_by_ind(df = obj)
head(missing_data.df)
write.csv(x = missing_data.df, file = "03_results/missing_data_per_indiv.csv", row.names = F)

summary(missing_data.df$ind.per.missing)  # max: 0.127 (12.7%); mean: 2.9%
sd(missing_data.df$ind.per.missing)       # sd: 3.0%


#### 02. Coast-specific analyses ####
obj.sep <- seppop(obj)

## Repool based on the dataset
if(dataset=="all"){
  
  # All Atlantic pops being repooled
  obj_atlantic <- repool(obj.sep$EQB, obj.sep$NFL, obj.sep$LAB)
  
  # All Pacific pops being repooled
  obj_pacific  <- repool(obj.sep$NBC, obj.sep$SOG, obj.sep$ORE, obj.sep$CAL)
  
}else if(dataset=="balanced"){
  
  # 'Balanced' Atlantic pops being repooled
  obj_atlantic <- repool(obj.sep$EQB, obj.sep$NFL)
  
  # 'Balanced' Pacific pops being repooled
  obj_pacific <- repool(obj.sep$SOG, obj.sep$ORE)
  
}

obj_atlantic
obj_pacific


## Re-calculate AF to remove low MAF variants
maf_filt(data = obj_atlantic, maf = 0.01)
obj_atlantic <- obj_maf_filt
myFreq.atl <- myFreq

maf_filt(data = obj_pacific, maf = 0.01)
obj_pacific <- obj_maf_filt
myFreq.pac <- myFreq


## Run PCA to view sample relationships in PCA prior to pruning putative relatives
pca_from_genind(data = obj_atlantic
                , PCs_ret = 4
                , plot_eigen = TRUE
                , plot_allele_loadings = TRUE
                , retain_pca_obj = TRUE
                , colour_file = "00_archive/harbour_seal_pops_colours.csv"
)

# Rename so the PCA figures are not overwritten
file.copy(from = "03_results/pca_samples_PC1_v_PC2.pdf", to = "03_results/pca_samples_PC1_v_PC2_atl_all_inds.pdf", overwrite = TRUE)
file.copy(from = "03_results/pca_samples_PC3_v_PC4.pdf", to = "03_results/pca_samples_PC3_v_PC4_atl_all_inds.pdf", overwrite = TRUE)

# Retain and save out the PCA scores
pca_scores_result <- pca.obj$scores
write.csv(x = pca_scores_result, file = "03_results/pca_scores_result_atlantic_w_relatives.csv", quote = F, row.names = T)


pca_from_genind(data = obj_pacific
                , PCs_ret = 4
                , plot_eigen = TRUE
                , plot_allele_loadings = TRUE
                , retain_pca_obj = TRUE
                , colour_file = "00_archive/harbour_seal_pops_colours.csv"
)

# Rename so the PCA figures are not overwritten
file.copy(from = "03_results/pca_samples_PC1_v_PC2.pdf", to = "03_results/pca_samples_PC1_v_PC2_pac_all_inds.pdf", overwrite = TRUE)
file.copy(from = "03_results/pca_samples_PC3_v_PC4.pdf", to = "03_results/pca_samples_PC3_v_PC4_pac_all_inds.pdf", overwrite = TRUE)

# Retain and save out the PCA scores
pca_scores_result <- pca.obj$scores
write.csv(x = pca_scores_result, file = "03_results/pca_scores_result_pacific_w_relatives.csv", quote = F, row.names = T)


## Identify putative close relatives for removal using relatedness
# Calculate inter-individual relatedness
relatedness_calc(data = obj_atlantic, datatype = "SNP") # will output as "03_results/kinship_analysis_<date>.Rdata"
date <- format(Sys.time(), "%Y-%m-%d")
file.copy(from = paste0("03_results/kinship_analysis_", date, ".Rdata"), to = paste0("03_results/kinship_analysis_atl_", date, ".Rdata"))

# Plot
relatedness_plot(file = paste0("03_results/kinship_analysis_atl_", date, ".Rdata"), same_pops = TRUE, plot_by = "codes", pdf_width = 7, pdf_height = 5)
file.copy(from = paste0("03_results/relatedness_ritland_", date, ".pdf"), to = paste0("03_results/relatedness_ritland_atl_", date, ".pdf"), overwrite = T)
file.copy(from = paste0("03_results/relatedness_quellergt_", date, ".pdf"), to = paste0("03_results/relatedness_quellergt_atl_", date, ".pdf"), overwrite = T)
file.copy(from = paste0("03_results/relatedness_wang_", date, ".pdf"), to = paste0("03_results/relatedness_wang_atl_", date, ".pdf"), overwrite = T)

# Calculate inter-individual relatedness
relatedness_calc(data = obj_pacific, datatype = "SNP") # will output as "03_results/kinship_analysis_<date>.Rdata"
date <- format(Sys.time(), "%Y-%m-%d")
file.copy(from = paste0("03_results/kinship_analysis_", date, ".Rdata"), to = paste0("03_results/kinship_analysis_pac_", date, ".Rdata"))

# Plot
relatedness_plot(file = paste0("03_results/kinship_analysis_pac_", date, ".Rdata"), same_pops = TRUE, plot_by = "codes", pdf_width = 7, pdf_height = 5)
file.copy(from = paste0("03_results/relatedness_ritland_", date, ".pdf"), to = paste0("03_results/relatedness_ritland_pac_", date, ".pdf"), overwrite = T)
file.copy(from = paste0("03_results/relatedness_quellergt_", date, ".pdf"), to = paste0("03_results/relatedness_quellergt_pac_", date, ".pdf"), overwrite = T)
file.copy(from = paste0("03_results/relatedness_wang_", date, ".pdf"), to = paste0("03_results/relatedness_wang_pac_", date, ".pdf"), overwrite = T)


save.image("03_results/output_coast-sp_relatedness.Rdata")

### Move to next script (hs_popn_analysis_part_2.R)
