# General analysis of harbour seal data (Part 3)
# 2020-10-13 initialized
# start by running '01_scripts/hs_popn_analysis_part2.R' 
#  or loading the previous image

# Load previous image
load(file = "03_results/output_coast-sp_pruned_analysis.Rdata")

## Private Alleles, Burrard Inlet
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


#save.image(file = "03_results/output_coast-sp_pruned_analysis.Rdata")

# Move to next script, hs_popn_analysis_part3.R


#### 05. Identify top FST markers for several specific BC datasets ####
##### NBC vs SOG, no Burrard ####
obj_pacific_filt.sep <- seppop(x = obj_pacific) 
obj_pacific_filt <- repool(obj_pacific_filt.sep$NBC, obj_pacific_filt.sep$SOG)

# Remove Burrard Inlet outlier individuals
outlier_samples.id <- c("SOG_108", "SOG_131", "SOG_135", "SOG_115", "SOG_127", "SOG_141", "SOG_125", "SOG_136", "SOG_102")
keep.inds <- setdiff(indNames(obj_pacific_filt), outlier_samples.id)
obj_pacific_filt <- obj_pacific_filt[keep.inds]

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
per_loc_stats_SOG_no_Burrard_v_NBC.df <- per_loc_stats.df



##### ORE vs SOG, no Burrard ####
obj_pacific_filt <- repool(obj_pacific_filt.sep$ORE, obj_pacific_filt.sep$SOG)

# Remove Burrard Inlet outlier individuals
outlier_samples.id <- c("SOG_108", "SOG_131", "SOG_135", "SOG_115", "SOG_127", "SOG_141", "SOG_125", "SOG_136", "SOG_102")
keep.inds <- setdiff(indNames(obj_pacific_filt), outlier_samples.id)
obj_pacific_filt <- obj_pacific_filt[keep.inds]

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
per_loc_stats_SOG_no_Burrard_v_ORE.df <- per_loc_stats.df


##### SOG vs. Burrard ####
obj_pacific_filt <- obj_pacific_filt.sep$SOG

# Create separate population for Burrard Inlet outlier individuals
outlier_samples.id <- c("SOG_108", "SOG_131", "SOG_135", "SOG_115", "SOG_127", "SOG_141", "SOG_125", "SOG_136", "SOG_102")

# Obtain the pop for this object
SOG_only.vec <- pop(obj_pacific_filt)
SOG_only.vec <- as.character(SOG_only.vec) # convert to character
SOG_only.vec[which(indNames(obj_pacific_filt) %in% outlier_samples.id)] <- "BUR" # define the outlier samples as 'BUR'
SOG_only.vec

pop(obj_pacific_filt) <- SOG_only.vec # Re-assign the population with the new BUR and SOG designations
indNames(obj_pacific_filt)

obj_pacific_filt

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
per_loc_stats_SOG_vs_Burrard.df <- per_loc_stats.df

##### Per region HOBS ##### 
names(obj_pacific_filt.sep)
# Inspect HOBS for BC inds
obj <- repool(obj_pacific_filt.sep$SOG, obj_pacific_filt.sep$NBC)
per_locus_stats(data = obj)
per_loc_stats_BC.df <- per_loc_stats.df

# Inspect HOBS for US inds
obj <- repool(obj_pacific_filt.sep$ORE, obj_pacific_filt.sep$CAL)
per_locus_stats(data = obj)
per_loc_stats_US.df <- per_loc_stats.df


##### Comparing data together ##### 
head(per_loc_stats_SOG_no_Burrard_v_ORE.df)
head(per_loc_stats_SOG_vs_Burrard.df)
head(per_loc_stats_SOG_no_Burrard_v_NBC.df)

SOG_versus_no_burrard.df <- merge(x = per_loc_stats_SOG_no_Burrard_v_NBC.df, per_loc_stats_SOG_no_Burrard_v_ORE.df, by = "mname")

BC_and_US_sp_HOBS.df <- merge(x = per_loc_stats_BC.df, y = per_loc_stats_US.df, by = "mname")

pdf(file = "03_results/FST_and_HOBS_region-specific.pdf", width = 9, height = 5)
par(mfrow=c(1,2))
plot(x = SOG_versus_no_burrard.df$Fst.x, y = SOG_versus_no_burrard.df$Fst.y
     , xlab = expression(italic(F)[ST] ~ SOG ~ vs. ~ NBC)
     , ylab = expression(italic(F)[ST] ~ SOG ~ vs. ~ ORE)
     , las = 1
)
text(x = 0.6, y = 0.5, labels = paste0("n = ", nrow(SOG_versus_no_burrard.df)))

plot(x = BC_and_US_sp_HOBS.df$Hobs.x, y = BC_and_US_sp_HOBS.df$Hobs.y
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

table(myFreq.pac < 0.1)[2] / length(myFreq.pac) # 57.9%
table(myFreq.atl < 0.1)[2] / length(myFreq.atl) # 47.9%

table(myFreq.pac < 0.05)[2] / length(myFreq.pac) # 37.0%
table(myFreq.atl < 0.05)[2] / length(myFreq.atl) # 31.0%


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


##### FINAL ALL POPN ANALYSES TO BE ADDED HERE #####
#### 02. Global PCA, FST, dendrogram ####

# #### Note: will do this after identifying HWE and related outliers ####
# # note: have not removed HWE outliers yet, ok?  (#TODO#)
# ## PCA
# pca_from_genind(data = obj
#                 , PCs_ret = 4
#                 , plot_eigen = TRUE
#                 , plot_allele_loadings = TRUE
#                 , retain_pca_obj = TRUE
#                 , colour_file = "00_archive/harbour_seal_pops_colours.csv"
# )
# 
# # Manually run pca_from_genind to pull out the pca1 obj (note: should assign this out to the global enviro as default)
# pca_scores_result <- pca.obj$scores
# write.csv(x = pca_scores_result, file = "03_results/pca_scores_result.csv", quote = F, row.names = T)
# 
# ## FST
# calculate_FST(format = "genind", dat = obj, separated = FALSE, bootstrap = TRUE)
# 
# ## Dendrogram
# make_tree(boot_obj = obj, bootstrap = TRUE, nboots = 10000, dist_metric = "edwards.dist", separated = FALSE)
# 
# #### /END/ Note: will do this after identifying HWE and related outliers ####





#### 0.4 Export ####
# Write out object
save.image(file = "03_results/completed_analysis.RData")
