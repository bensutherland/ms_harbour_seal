# General analysis of harbour seal data
# 2020-10-13 initialized
# start by sourcing simple_pop_stats then quit the menu

# Set custom output folder
set_output_folder(result_path= "/hdd/harbour_seal/simple_pop_stats_2020-10-13/03_results/")

# Load genepop and characterize
load_genepop(datatype = "SNP")
# from here: /hdd/harbour_seal/simple_pop_stat/02_input_data/bhs_single_snp_p5_r0.7_maf0.01_out_rem_2020-10-24.gen 

# Drop oddball samples (identified from PCA in prev iteration)
rownames(obj$tab)
removeInd <- c("NFL_20192991", "NFL_20042268")

# remove individuals from *genind object*
# note that in this step there's no longer a comma needed before the
# closing square bracket
# obj.bck <- obj
obj <- obj[!row.names(obj@tab) %in% removeInd]

unique(pop(obj))

# Clean up pop names to just the region
pop(obj) <- gsub(pattern = "_.*", replacement = "", x = pop(obj))

characterize_genepop(df = obj, pdf_width = 7, pdf_height = 5, cex_names = 0.8, N=30)

# Make a dendrogram
make_tree(bootstrap = TRUE, boot_obj = obj, nboots = 10000, dist_metric = "edwards.dist", separated = FALSE)

# PCA
pca_from_genind(data = obj, PCs_ret = 3, plot_eigen=TRUE, plot_allele_loadings=TRUE, colour_file = "00_archive/harbour_seal_pops_colours.csv")

# Manually run pca_from_genind to pull out the pca1 obj (note: should assign this out to the global enviro as default)
pca_scores_result <- pca.obj$scores
write.csv(x = pca_scores_result, file = "03_results/pca_scores_result.csv", quote = F, row.names = T)

# Calculate FST
calculate_FST(format="genind", dat = obj, separated = FALSE)

# Calculate FST with 95% CI

# Change from genind file to hfstat
data.hf <- genind2hierfstat(obj)
rownames(data.hf) <- indNames(obj)

dim(data.hf)
data.hf[1:10,1:6]

# Pairwise Fst w/ bootstrapping (hierfstat)
boot.fst <- boot.ppfst(dat = data.hf, nboot = 1000, quant = c(0.025,0.975))
boot.fst

# Collect output
lower.limit <- t(boot.fst$ll)
upper.limit <- boot.fst$ul
upper.limit[is.na(upper.limit)] <- 0
lower.limit[is.na(lower.limit)] <- 0
boot.fst.output <- upper.limit + lower.limit
boot.fst.output

filename <- paste(result.path, datatype, "_boot_fst_output.csv", sep = "")
write.csv(x = boot.fst.output, file = filename)
