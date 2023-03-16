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


#### 02. Global PCA, FST, dendrogram ####
# note: have not removed HWE outliers yet, ok?  (#TODO#)
## PCA
pca_from_genind(data = obj
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
calculate_FST(format = "genind", dat = obj, separated = FALSE, bootstrap = TRUE)


#### 03. Coast-specific, Atlantic ####
obj.sep <- seppop(obj)
#obj_atlantic <- repool(obj.sep$EQB, obj.sep$NFL, obj.sep$LAB) # all pops
obj_atlantic <- repool(obj.sep$EQB, obj.sep$NFL) # balanced, normalized
obj_atlantic

## Re-calculate AF to remove low MAF variants
obj.gl <- gi2gl(gi = obj_atlantic, parallel = T) # Convert to genlight

# Calculate frequency of second allele
myFreq <- glMean(obj.gl)

# Ensure each locus second allele is the minor allele
for(i in 1:length(myFreq)){
  
  # if the second allele is > 0.5, this would be considered the major allele, and so calculate the minor allele frequency
  if(myFreq[i] > 0.5){
    
    myFreq[i] <- 1-myFreq[i]
    
  }else{
    
    myFreq[i] <- myFreq[i]
    
  }
  
}

## MAF filter
MAF_rem_final <- names(myFreq[which(myFreq < 0.01)])
length(MAF_rem_final)
markers_to_keep <- setdiff(x = locNames(obj_atlantic), y = MAF_rem_final)
length(markers_to_keep)
obj_atlantic <- obj_atlantic[, loc=markers_to_keep]
obj_atlantic

# Keep AF of only the retained variants
myFreq <- myFreq[which(myFreq >= 0.01)]
length(myFreq) # should match the number of variants kept in obj_atlantic

# Save
myFreq.atl <- myFreq
rm(myFreq)


## HWE filter
hwe_eval(data = obj_atlantic, alpha = 0.01)
head(per_locus_hwe_NFL.df)

# Which col contains pval? 
col.oi <- grep(pattern = "Pr", x = colnames(per_locus_hwe_NFL.df))

# Identify mnames of outliers
hwe_outlier_mname_EQB.vec     <-   per_locus_hwe_EQB.df[per_locus_hwe_EQB.df[, col.oi] < 0.01, "mname"]
hwe_outlier_mname_NFL.vec    <-    per_locus_hwe_NFL.df[per_locus_hwe_NFL.df[, col.oi] < 0.01, "mname"]
hwe_outlier_mname_LAB.vec    <-    per_locus_hwe_LAB.df[per_locus_hwe_LAB.df[, col.oi] < 0.01, "mname"]

# How many outliers (p < 0.01) per population
length(hwe_outlier_mname_EQB.vec)   #  173 markers out of HWE
length(hwe_outlier_mname_NFL.vec)   #  233 markers out of HWE
length(hwe_outlier_mname_LAB.vec)   #   85 markers out of HWE


# How many unique HWE deviating markers?  
markers_to_drop <- c(hwe_outlier_mname_EQB.vec, hwe_outlier_mname_NFL.vec, hwe_outlier_mname_LAB.vec)
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
write.csv(x = pca_scores_result, file = "03_results/pca_scores_result_atlantic.csv", quote = F, row.names = T)


## FST
calculate_FST(format = "genind", dat = obj_atlantic, separated = FALSE, bootstrap = TRUE)

## Relatedness
obj_atlantic

# Calculate inter-individual relatedness
relatedness_calc(data = obj_atlantic, datatype = "SNP") # will output as "03_results/kinship_analysis_<date>.Rdata"

# Explore the data
head(output$relatedness)


# Plot
relatedness_plot(file = "03_results/kinship_analysis_2023-03-01.Rdata", same_pops = TRUE, plot_by = "codes", pdf_width = 7, pdf_height = 5)
# It looks like Ritland > 0.1 is outlier 


# Plot distribution of relatedness statistics
pdf(file = "03_results/hist_relatedness_Ritland.pdf", width = 6, height = 3.5)
par(mfrow=c(1,1))
hist(output$relatedness$ritland, main = "", xlab = "relatedness (Ritland)", las = 1)
abline(v = 0.1, lty = 3)
text(paste0("median = ", median(output$relatedness$ritland))
     , y = 600, x = 0.15)
text(paste0("mean = "  , round(mean(output$relatedness$ritland), digits = 3))
     , y = 450, x = 0.15)
dev.off()


# Which pairs have high relatedness using the ritland statistic? 
highly_related.df <- output$relatedness[output$relatedness$ritland >= 0.1, c("ind1.id", "ind2.id", "group", "ritland")]
nrow(highly_related.df) # 17
highly_related.df

# How many individuals does this involve that are highly related? 
length(unique(x = c(highly_related.df$ind1.id, highly_related.df$ind2.id)))

# How many individuals does this involve? 
length(unique(x = c(output$relatedness$ind1.id, output$relatedness$ind2.id)))

write.table(x = highly_related.df, file = "03_results/highly_related_2023-03-01.csv", quote = F, sep = ","
            , row.names = F
)


# Move all results into an 'Atlantic' folder, then proceed to Pacific analysis


#### 04. Coast-specific, Pacific ####
obj.sep <- seppop(obj)

# Repool based on the dataset
if(dataset=="all"){
 
  # All Pacific pops being repooled
  obj_pacific <- repool(obj.sep$NBC, obj.sep$SOG, obj.sep$ORE, obj.sep$CAL) 
  
}else if(dataset=="balanced"){
  
  # 'Balanced' pops being repooled
  obj_pacific <- repool(obj.sep$SOG, obj.sep$ORE)
    
}

## Re-calculate AF to remove low MAF variants
obj.gl <- gi2gl(gi = obj_pacific, parallel = T) # Convert to genlight

# Calculate frequency of second allele
myFreq <- glMean(obj.gl)

# Ensure each locus second allele is the minor allele
for(i in 1:length(myFreq)){
  
  if(myFreq[i] > 0.5){
    
    myFreq[i] <- 1-myFreq[i]
    
  }else{
    
    myFreq[i] <- myFreq[i]
    
  }
  
}

## Final MAF filter
MAF_rem_final <- names(myFreq[which(myFreq < 0.01)])
length(MAF_rem_final)
markers_to_keep <- setdiff(x = locNames(obj_pacific), y = MAF_rem_final)
obj_pacific <- obj_pacific[, loc=markers_to_keep]
obj_pacific

# Keep AF of only the retained variants
myFreq <- myFreq[which(myFreq >= 0.01)]
length(myFreq)

# Save
myFreq.pac <- myFreq
rm(myFreq)

## HWE filter
hwe_eval(data = obj_pacific, alpha = 0.01)
head(per_locus_hwe_SOG.df)

# Identify column with the p-val
col.oi <- grep(pattern = "Pr", x = colnames(per_locus_hwe_SOG.df))

# Identify mnames of outliers
hwe_outlier_mname_SOG.vec     <-   per_locus_hwe_SOG.df[per_locus_hwe_SOG.df[, col.oi] < 0.01, "mname"]
hwe_outlier_mname_NBC.vec    <-    per_locus_hwe_NBC.df[per_locus_hwe_NBC.df[, col.oi] < 0.01, "mname"]
hwe_outlier_mname_CAL.vec    <-    per_locus_hwe_CAL.df[per_locus_hwe_CAL.df[, col.oi] < 0.01, "mname"]
hwe_outlier_mname_ORE.vec    <-    per_locus_hwe_ORE.df[per_locus_hwe_ORE.df[, col.oi] < 0.01, "mname"]

# How many outliers (p < 0.01) per population
length(hwe_outlier_mname_SOG.vec)   #  513 markers out of HWE
length(hwe_outlier_mname_NBC.vec)   #  155 markers out of HWE
length(hwe_outlier_mname_CAL.vec)   #  323 markers out of HWE
length(hwe_outlier_mname_ORE.vec)   #  306 markers out of HWE


# How many unique HWE deviating markers?  
markers_to_drop <- c(hwe_outlier_mname_SOG.vec, hwe_outlier_mname_NBC.vec, hwe_outlier_mname_CAL.vec, hwe_outlier_mname_ORE.vec)
length(markers_to_drop)             # 1297 markers out of HWE in at least one population
markers_to_drop <- unique(markers_to_drop)
length(markers_to_drop)             #  996 unique markers out of HWE in at least one population
markers_to_keep <- setdiff(x = locNames(obj_pacific), y = markers_to_drop)
length(markers_to_keep) # 7937 markers to keep

obj_pacific <- obj_pacific[, loc=markers_to_keep]
obj_pacific


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

#### Relatedness ####
obj_pacific

# Calculate inter-individual relatedness
relatedness_calc(data = obj_pacific, datatype = "SNP") # will output as "03_results/kinship_analysis_<date>.Rdata"

# Explore the data
head(output$relatedness)


# Plot
relatedness_plot(file = "03_results/kinship_analysis_2023-03-01.Rdata", same_pops = TRUE, plot_by = "codes", pdf_width = 7, pdf_height = 5)
# It looks like Ritland > 0.1 is outlier for both CACA and SOSO


# Plot distribution of relatedness statistics
pdf(file = "03_results/hist_relatedness_Ritland.pdf", width = 6, height = 3.5)
par(mfrow=c(1,1))
hist(output$relatedness$ritland, main = "", xlab = "relatedness (Ritland)", las = 1)
abline(v = 0.1, lty = 3)
text(paste0("median = ", median(output$relatedness$ritland))
     , y = 800, x = 0.15)
text(paste0("mean = "  , round(mean(output$relatedness$ritland), digits = 3))
     , y = 650, x = 0.15)
dev.off()


# Which pairs have more than 0.25 relatedness using the wang statistic? 
highly_related.df <- output$relatedness[output$relatedness$ritland >= 0.1, c("ind1.id", "ind2.id", "group", "ritland")]
nrow(highly_related.df) # 54
highly_related.df

# How many individuals does this involve that are highly related? 
length(unique(x = c(highly_related.df$ind1.id, highly_related.df$ind2.id)))

# How many individuals does this involve? 
length(unique(x = c(output$relatedness$ind1.id, output$relatedness$ind2.id)))

write.table(x = highly_related.df, file = "03_results/highly_related_2023-03-01.csv", quote = F, sep = ","
            , row.names = F
)


#### Private Alleles, Burrard Inlet ####
# Create backup
#obj_pacific.bck <- obj_pacific 
#obj_pacific <- obj_pacific.bck

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

# Confirms it works:
indNames(obj.SOG)[which(indNames(obj.SOG) %in% outlier_samples.id)]

# Add strat object to genind
strata(obj.SOG) <- strata.df
obj.SOG
table(strata(obj.SOG))

per_repunit.privallele <- private_alleles(obj.SOG, alleles ~ repunit)

#per_repunit.privallele[, 1:5]
dim(per_repunit.privallele)

table(per_repunit.privallele["SOG", ])
table(per_repunit.privallele["Burrard", ])

which(per_repunit.privallele["Burrard", ]==9)
which(per_repunit.privallele["Burrard", ]==6)
which(per_repunit.privallele["Burrard", ]==5)


# # Get raw number of private alleles per locus
# pal <- private_alleles(obj.SOG, locus ~ repunit, count.alleles = FALSE)
# table(pal) # this only gives a 0 or 1, does not count alleles, allows one to see exact how many private alleles exist
# # per repunit
# # This shows the number of rows that have a private allele
# rowSums(pal)


#### 05. BC-specific marker test for top FST ####
obj_pacific_filt.sep <- seppop(x = obj_pacific)
obj_pacific_filt <- repool(obj_pacific_filt.sep$NBC, obj_pacific_filt.sep$SOG)

## Re-calculate AF to remove low MAF variants
obj.gl <- gi2gl(gi = obj_pacific_filt, parallel = T) # Convert to genlight

# Calculate frequency of second allele
myFreq <- glMean(obj.gl)

# Ensure each locus second allele is the minor allele
for(i in 1:length(myFreq)){
  
  if(myFreq[i] > 0.5){
    
    myFreq[i] <- 1-myFreq[i]
    
  }else{
    
    myFreq[i] <- myFreq[i]
    
  }
  
}

## Final MAF filter
MAF_rem_final <- names(myFreq[which(myFreq < 0.01)])
length(MAF_rem_final)
markers_to_keep <- setdiff(x = locNames(obj_pacific_filt), y = MAF_rem_final)
obj_pacific_filt <- obj_pacific_filt[, loc=markers_to_keep]
obj_pacific_filt


# Calculate per-locus FST
per_locus_stats(data = obj_pacific_filt)

# Save for later
per_loc_stats_BC.df <- per_loc_stats.df


## Compare this to the per-locus stats for Pacific-wide, to see how much overlap
head(per_loc_stats_BC.df)
nrow(per_loc_stats_BC.df)
head(per_loc_stats_pac.df)
nrow(per_loc_stats_pac.df)

## Merge the all Pacific and the BC-only pops to compare FST per locus
per_loc_stats_PAC_BC.df <- merge(x = per_loc_stats_pac.df, y = per_loc_stats_BC.df, by = "mname")
nrow(per_loc_stats_PAC_BC.df)
head(per_loc_stats_PAC_BC.df)


pdf(file = "03_results/per_locus_FST_HOBS_all_Pacific_or_BC.pdf", width = 8, height = 3.5)
par(mfrow=c(1,2))
plot(per_loc_stats_PAC_BC.df$Fst.x, per_loc_stats_PAC_BC.df$Fst.y
     , xlab = "per locus FST, all Pacific", ylab = "per locus FST, BC only"
     )
plot(per_loc_stats_PAC_BC.df$Hobs.x, per_loc_stats_PAC_BC.df$Hobs.y
     , xlab = "per locus HOBS, all Pacific", ylab = "per locus HOBS, BC only"
     )
dev.off()

#### Note: we could also do a US only one too... ###


### 06. Multi-coast plotting ####

## MAF
# Remove the HWE-filtered variants
myFreq.pac <- myFreq.pac[names(myFreq.pac) %in% locNames(obj_pacific)]
myFreq.atl <- myFreq.atl[names(myFreq.atl) %in% locNames(obj_atlantic)]

#### SHORTCUT PLOTTING HERE ####

table(myFreq.pac < 0.1)
table(myFreq.atl < 0.1)

length(myFreq.pac)
length(myFreq.atl)

table(myFreq.pac < 0.05)[2] / length(myFreq.pac)
table(myFreq.atl < 0.05)[2] / length(myFreq.atl)


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
text(x = 0.4, y = 1000, labels = paste("n = ", length(myFreq.atl), " loci", sep = "" ))
dev.off()

#### TODO: do not plot with set ylim to be able to see the distribution independent of total markers ####
### Add text of % markers with 0.01 < MAF < 0.05
### Get same plot from the denovo, with the values as well

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

### HOBS
# Plot marker HOB
pdf(file = paste0("03_results/HOBS_hist_pac_atl.pdf"), width = 7, height = 4)
par(mfrow=c(1,2))
hist(per_loc_stats_pac.df$Hobs
     , ylim = c(0,1750)
     , main = ""
     , las = 1
     , xlab = "Hobs, Pacific"
       )

hist(per_loc_stats_atl.df$Hobs
     , ylim = c(0,1750)
     , main = ""
     , las = 1
     , xlab = "Hobs, Atlantic"
     )

dev.off()

#### 0.4 Export ####
# Write out object
save.image(file = "03_results/completed_analysis.RData")




