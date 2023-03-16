# Analysis of de novo, balanced n, coast-specific data to explore MAF distn
# 2023-03-02 initialized
# start by sourcing simple_pop_stats then quit the menu

# Set up null lists (*ONLY RUN ONCE*)
per_loc_stats.list <- list()
myFreq.list        <- list()

#### 01. Load data ####
# Load genepop and characterize
load_genepop(datatype = "SNP")

dataset <- "pacific" # from here:  "02_input_data/bhs_pac_p2_r0.7_maf0.01_2023-03-02.gen"
#dataset <- "atlantic" # from here:  "02_input_data/bhs_atl_p2_r0.7_maf0.01_2023-03-02.gen"

# Clean up pop names
pop(obj) <- gsub(pattern = "_.*", replacement = "", x = pop(obj))
unique(pop(obj))



# #### 02. HWE filter ####
# ### NOTE: it is not clear whether we need to run HWE check; it would be easier if not. However, the reference-based approach did. 
# hwe_eval(data = obj, alpha = 0.01)
# 
# ## Filter out HWE deviators (depends on dataset)
# if(dataset=="pacific"){
#   
#   head(per_locus_hwe_ORE.df)
#   
#   # Which col contains pval? 
#   col.oi <- grep(pattern = "Pr", x = colnames(per_locus_hwe_ORE.df))
#   
#   # Identify mnames of outliers
#   hwe_outlier_mname_SOG.vec     <-   per_locus_hwe_SOG.df[per_locus_hwe_SOG.df[, col.oi] < 0.01, "mname"]
#   hwe_outlier_mname_ORE.vec    <-    per_locus_hwe_ORE.df[per_locus_hwe_ORE.df[, col.oi] < 0.01, "mname"]
#   
#   # How many outliers (p < 0.01) per population
#   length(hwe_outlier_mname_SOG.vec)   #  159 markers out of HWE
#   length(hwe_outlier_mname_ORE.vec)   #  134 markers out of HWE
#   
#   # How many unique HWE deviating markers?
#   markers_to_drop <- c(hwe_outlier_mname_SOG.vec, hwe_outlier_mname_ORE.vec)
#   length(markers_to_drop)             # 293 markers out of HWE in at least one population
#   markers_to_drop <- unique(markers_to_drop)
#   length(markers_to_drop)             #  256 unique markers out of HWE in at least one population
#   markers_to_keep <- setdiff(x = locNames(obj), y = markers_to_drop)
#   length(markers_to_keep) # 4371 markers to keep
#   
#   obj <- obj[, loc=markers_to_keep]
#   obj
#   
# }else if(dataset=="atlantic"){
#   
#   head(per_locus_hwe_NFL.df)
#   
#   # Which col contains pval? 
#   col.oi <- grep(pattern = "Pr", x = colnames(per_locus_hwe_NFL.df))
#   
#   # Identify mnames of outliers
#   hwe_outlier_mname_EQB.vec     <-   per_locus_hwe_EQB.df[per_locus_hwe_EQB.df[, col.oi] < 0.01, "mname"]
#   hwe_outlier_mname_NFL.vec    <-    per_locus_hwe_NFL.df[per_locus_hwe_NFL.df[, col.oi] < 0.01, "mname"]
#   
#   # How many outliers (p < 0.01) per population
#   length(hwe_outlier_mname_EQB.vec)   #  xx markers out of HWE
#   length(hwe_outlier_mname_NFL.vec)   #  xx markers out of HWE
#   
#   # How many unique HWE deviating markers?  
#   markers_to_drop <- c(hwe_outlier_mname_EQB.vec, hwe_outlier_mname_NFL.vec)
#   length(markers_to_drop)             # xx total markers identified
#   markers_to_drop <- unique(markers_to_drop)
#   length(markers_to_drop)             # xx unique markers identified
#   
#   markers_to_keep <- setdiff(x = locNames(obj), y = markers_to_drop)
#   length(markers_to_keep) # xx markers to keep
#   
#   obj <- obj[, loc=markers_to_keep]
#   obj
#   
# }
# 
# #### 03. Excess HOBS filter ####
# ## Per locus statistics
# per_locus_stats(data = obj)
# 
# # Markers with HOBS > 0.5? 
# hobs.outliers <- per_loc_stats.df[per_loc_stats.df$Hobs > 0.5, "mname"] 
# length(hobs.outliers) # xx markers
# 
# keep <- setdiff(x = locNames(obj), y = hobs.outliers)
# length(keep)
# obj <- obj[, loc=keep] 
# obj
# 
# # Re-run per loc stats
# per_locus_stats(data = obj)
# 
# # Save for later
# per_loc_stats.list[[dataset]] <- per_loc_stats.df



#### 04. AF calc ####
## Re-calculate AF to remove low MAF variants
obj.gl <- gi2gl(gi = obj, parallel = T) # Convert to genlight

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

myFreq.list[[dataset]] <- myFreq


#### RETURN AND RUN ATL


####  04. Plotting ####
# Plot
pdf(file = paste0("03_results/MAF_hist_pac_atl_denovo.pdf"), width = 7, height = 4)
par(mfrow=c(1,2))
hist(myFreq.list[["pacific"]]
     #, proba=T # note: does not sum to 1, not worth using
     , col="grey", xlab = "MAF, Pacific"
     , main = ""
     , ylim = c(0, 900)
     , ylab = "Number of loci"
     , las = 1
     , breaks = 20
)
text(x = 0.4, y = 700, labels = paste("n = ", length(myFreq.list[["pacific"]]), " loci", sep = "" ))

hist(myFreq.list[["atlantic"]]
     #, proba=T # note: does not sum to 1, not worth using
     , col="grey", xlab = "MAF, Atlantic"
     , main = ""
     , ylim = c(0, 900)
     , ylab = "Number of loci"
     , las = 1
     , breaks = 20
)
text(x = 0.4, y = 700, labels = paste("n = ", length(myFreq.list$atlantic), " loci", sep = "" ))
dev.off()
