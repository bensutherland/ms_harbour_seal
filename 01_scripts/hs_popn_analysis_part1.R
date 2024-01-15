# General analysis of harbour seal data (Part 1)
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
## Separate populations
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
# Atlantic
pca_from_genind(data = obj_atlantic
                , PCs_ret = 4
                , plot_eigen = TRUE
                , plot_allele_loadings = FALSE
                , retain_pca_obj = TRUE
                , colour_file = "00_archive/harbour_seal_pops_colours.csv"
)

# Rename so the PCA figures are not overwritten
file.copy(from = "03_results/pca_samples_PC1_v_PC2.pdf", to = "03_results/pca_samples_PC1_v_PC2_atl_all_inds.pdf", overwrite = TRUE)
file.copy(from = "03_results/pca_samples_PC3_v_PC4.pdf", to = "03_results/pca_samples_PC3_v_PC4_atl_all_inds.pdf", overwrite = TRUE)
pca_atl_w_sibs.obj <- pca.obj

# Retain and save out the PCA scores
pca_scores_result <- pca.obj$scores
write.csv(x = pca_scores_result, file = "03_results/pca_scores_result_atlantic_w_relatives.csv", quote = F, row.names = T)

# Pacific
pca_from_genind(data = obj_pacific
                , PCs_ret = 4
                , plot_eigen = TRUE
                , plot_allele_loadings = FALSE
                , retain_pca_obj = TRUE
                , colour_file = "00_archive/harbour_seal_pops_colours.csv"
)

# Rename so the PCA figures are not overwritten
file.copy(from = "03_results/pca_samples_PC1_v_PC2.pdf", to = "03_results/pca_samples_PC1_v_PC2_pac_all_inds.pdf", overwrite = TRUE)
file.copy(from = "03_results/pca_samples_PC3_v_PC4.pdf", to = "03_results/pca_samples_PC3_v_PC4_pac_all_inds.pdf", overwrite = TRUE)
pca_pac_w_sibs.obj <- pca.obj

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
file.copy(from = paste0("03_results/pairwise_relatedness_output_all_", date, ".txt"), to = paste0("03_results/pairwise_relatedness_output_all_atl_", date, ".txt"))

# Calculate inter-individual relatedness
relatedness_calc(data = obj_pacific, datatype = "SNP") # will output as "03_results/kinship_analysis_<date>.Rdata"
date <- format(Sys.time(), "%Y-%m-%d")
file.copy(from = paste0("03_results/kinship_analysis_", date, ".Rdata"), to = paste0("03_results/kinship_analysis_pac_", date, ".Rdata"))

# Plot
relatedness_plot(file = paste0("03_results/kinship_analysis_pac_", date, ".Rdata"), same_pops = TRUE, plot_by = "codes", pdf_width = 7, pdf_height = 5)
file.copy(from = paste0("03_results/relatedness_ritland_", date, ".pdf"), to = paste0("03_results/relatedness_ritland_pac_", date, ".pdf"), overwrite = T)
file.copy(from = paste0("03_results/relatedness_quellergt_", date, ".pdf"), to = paste0("03_results/relatedness_quellergt_pac_", date, ".pdf"), overwrite = T)
file.copy(from = paste0("03_results/relatedness_wang_", date, ".pdf"), to = paste0("03_results/relatedness_wang_pac_", date, ".pdf"), overwrite = T)
file.copy(from = paste0("03_results/pairwise_relatedness_output_all_", date, ".txt"), to = paste0("03_results/pairwise_relatedness_output_all_pac_", date, ".txt"))


### Additional analysis for reviewer request (2024-01-07): within PCA grouping relatedness evaluation to 
##    determine whether relatedness is overestimated due to coast-wide allele frequencies
# Isolate to BC only samples
obj_pacific.sep <- seppop(x = obj_pacific)
obj_BC <- repool(obj_pacific.sep$NBC, obj_pacific.sep$SOG)
maf_filt(data = obj_BC, maf = 0.01)
obj_BC <- obj_maf_filt

relatedness_calc(data = obj_BC, datatype = "SNP") # will output as "03_results/kinship_analysis_<date>.Rdata"
date <- format(Sys.time(), "%Y-%m-%d")
file.copy(from = paste0("03_results/kinship_analysis_", date, ".Rdata"), to = paste0("03_results/kinship_analysis_BC_", date, ".Rdata"))

# Plot
relatedness_plot(file = paste0("03_results/kinship_analysis_BC_", date, ".Rdata"), same_pops = TRUE, plot_by = "codes", pdf_width = 7, pdf_height = 5)
file.copy(from = paste0("03_results/relatedness_ritland_", date, ".pdf"), to = paste0("03_results/relatedness_ritland_BC_", date, ".pdf"), overwrite = T)
file.copy(from = paste0("03_results/relatedness_quellergt_", date, ".pdf"), to = paste0("03_results/relatedness_quellergt_BC_", date, ".pdf"), overwrite = T)
file.copy(from = paste0("03_results/relatedness_wang_", date, ".pdf"), to = paste0("03_results/relatedness_wang_BC_", date, ".pdf"), overwrite = T)
file.copy(from = paste0("03_results/pairwise_relatedness_output_all_", date, ".txt"), to = paste0("03_results/pairwise_relatedness_output_all_BC_", date, ".txt"))

# What are the filenames
input_all_Pacific.FN <- paste0("03_results/pairwise_relatedness_output_all_pac_", date, ".txt")
input_all_BC.FN <- paste0("03_results/pairwise_relatedness_output_all_BC_", date, ".txt")

# Read in data
rel_pac.df <- read.table(file = input_all_Pacific.FN, header = T, sep = "\t")
head(rel_pac.df)
rel_pac.df$pair.ids <- paste0(rel_pac.df$ind1.id, "__", rel_pac.df$ind2.id)

rel_BC.df <- read.table(file = input_all_BC.FN, header = T, sep = "\t")
head(rel_BC.df)
rel_BC.df$pair.ids <- paste0(rel_BC.df$ind1.id, "__", rel_BC.df$ind2.id)

# What is the level and number of outliers for each? 
# original analysis with all Pacific
length(boxplot.stats(rel_pac.df[rel_pac.df$group=="SOSO", "ritland"])$out) # 33 pairs
outlier_cutoff_all_pac.val <- min(boxplot.stats(rel_pac.df[rel_pac.df$group=="SOSO", "ritland"])$out)    
outlier_cutoff_all_pac.val # 0.103

vals <- boxplot.stats(rel_BC.df[rel_BC.df$group=="SOSO", "ritland"])$out
vals <- vals[which(vals > 0)]
length(vals)   # 29 pairs
outlier_cutoff_BC_only.val <- min(vals)      
outlier_cutoff_BC_only.val # 0.033 

both_pac_and_bc_rel.df <- merge(x = rel_pac.df, y = rel_BC.df, by = "pair.ids")
dim(both_pac_and_bc_rel.df)
head(both_pac_and_bc_rel.df)

pdf(file = "03_results/estim_relat_BC_only_vs_all_pacific_loci.pdf", width = 8, height = 5.5)
plot(x = both_pac_and_bc_rel.df$ritland.x, y = both_pac_and_bc_rel.df$ritland.y
     , ylab = "Estim. relat. (Ritland), BC only loci"
     , xlab = "Estim. relat. (Ritland), all Pacific loci"
     , las = 1
     , ylim = c(-0.1, 0.22)
     , xlim = c(-0.1, 0.22)
     )
abline(h = outlier_cutoff_BC_only.val, lty = 2)
abline(v = outlier_cutoff_all_pac.val, lty = 2)
dev.off()

save.image("03_results/output_coast-sp_relatedness.Rdata")
# load("03_results/output_coast-sp_relatedness.Rdata")
### Move to next script (hs_inspect_outlier_output.R) then following (hs_popn_analysis_part_2.R)
