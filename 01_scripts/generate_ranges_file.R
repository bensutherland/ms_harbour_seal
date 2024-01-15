# Produce bed file for selecting a specific whitelist of loci out of a VCF
# Ben Sutherland, 2024-01-11

library("vcfR")

setwd("~/Documents/00_sbio/harbour_seal/ms_harbour_seal/")


#### 00. Set and read in data ####
# Select the locus filename
#loc.FN  <- "05_admixture/pacific_whitelist_loci.txt"   # Pacific
#loc.FN  <- "05_admixture/atlantic_whitelist_loci.txt" # Atlantic
#loc.FN  <- "05_admixture/both_coasts_whitelist_loci.txt" # Both

# Set output filename
output.FN <- gsub(pattern = ".txt", replacement = "_bed.txt", x = loc.FN)

# Read in locus names
loc.df <- read.delim(file = loc.FN, header = F)
loc.df <- as.data.frame(loc.df)
colnames(loc.df) <- "locus_id"
head(loc.df)     # ID in format of '825_55'
nrow(loc.df)

# Read in VCF (from the original stacks analysis output, single SNP per locus)
my_vcf <- read.vcfR(file = "02_input_data/populations.snps.vcf") # 10,847 variants

# Obtain info section of VCF
vcf_info <- my_vcf@fix
vcf_info.df <- as.data.frame(vcf_info)
head(vcf_info.df) # ID in format of '822:55:-'

#### 01. Match the ID column to the retain loci names ####
head(loc.df)
head(vcf_info.df)

# Format the IDs of the VCF to match the target locus names, making a new column
vcf_info.df$mname <- gsub(pattern = ":", replacement = "_", x = vcf_info.df$ID)
head(vcf_info.df)
vcf_info.df$mname <- gsub(pattern = "_\\-|_\\+", replacement = "", x = vcf_info.df$mname)
head(vcf_info.df)

table(vcf_info.df$mname %in% loc.df$locus_id) # Does this match? 
nrow(loc.df)

# Keep only the rows that are in the select list
vcf_info_select_loci.df <- vcf_info.df[vcf_info.df$mname %in% loc.df$locus_id, ]
dim(vcf_info_select_loci.df)
head(vcf_info_select_loci.df)

#### 02. Prepare a BED file ####
## Now use the CHROM and POS variable to create a bed file
bed.df <- vcf_info_select_loci.df[, c("CHROM", "POS")]
bed.df$POS <- as.numeric(bed.df$POS)
bed.df$left <- bed.df$POS-2
bed.df$right <- bed.df$POS+2
bed.df <- bed.df[,c("CHROM", "left", "right")]
head(bed.df)

## Write out 
write.table(x = bed.df, file = output.FN, quote = F, sep = "\t", col.names = F, row.names = F)

# Go back to terminal and use the bed file to extract specific loci from the VCF using bedtools 



