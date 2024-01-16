# Admixture analysis
# B. Sutherland, 2024-01-13

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

# Install and load packages
#

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)
rm(current.path)

options(scipen=99999999)

# Datasets
all_dataset <- "pvit_all"
pac_dataset <- "pvit_pac"
atl_dataset <- "pvit_atl"

# Choose which dataset to use
dataset <- all_dataset

k.list <- list()
k.list[["pvit_all"]] <- 6
k.list[["pvit_pac"]] <- 5
k.list[["pvit_atl"]] <- 5

## Other filenames
metadata.FN  <- "02_input_data/additional_file_s1_sample_metadata_2024-01-13.txt"
input_fam.FN <- paste0("05_admixture/", dataset, ".fam")
CV_error.FN <- paste0("05_admixture/", dataset, "_CV_error.txt")

# What is the data? 
dataset_files <- list.files(path = "05_admixture/", pattern = dataset, full.names = T)
dataset_files <- dataset_files[grep(pattern = ".Q", x = dataset_files)]
dataset_files <- head(sort(dataset_files), n = k.list[[dataset]])

# Read in order of samples 
fam <- read.table(file = input_fam.FN)
head(fam) # second column is indivs, and they are in the same order of the admixture matrix

# Read in all data
tbl <- NULL; tbl.list <- list()
for(i in 1:length(dataset_files)){
  
  # Read in data
  tbl <- read.table(file = dataset_files[[i]])
  
  # Assign indiv names to the admixture matrix
  rownames(tbl) <- fam$V2
  
  # Save it to the list
  tbl.list[[i]] <- tbl
  
}

str(tbl.list)

# Read in CV error
CV_error.df <- read.table(file = CV_error.FN, sep = ",")
CV_error.df
CV_error.df <- head(CV_error.df, n = k.list[[dataset]])

# Set cols
k1_col <- c("darkorchid3")
k2_col <- c("darkorchid3", "darkorange")
k3_col <- c("yellow", "darkorchid3", "darkorange")
k4_col <- c("darkorange", "darkorchid3",  "green", "yellow")
k5_col <- c("darkorange", "darkorchid3",  "green", "yellow", "black")

colour.list <- list()
colour.list[[1]] <- k1_col
colour.list[[2]] <- k2_col
colour.list[[3]] <- k3_col
colour.list[[4]] <- k4_col
colour.list[[5]] <- k5_col




#k_colours <- c("darkorange", "darkorchid3","green", "yellow")


length(tbl.list)
# Plot
pdf(file = paste0("05_admixture/", dataset, "_multi_k_admixture_plot.pdf")
    , width = 7.5, height = 12)
par(mfrow = c(k.list[[dataset]], 1)) # number of plots match k
for(i in 1:5){
  
  barplot(t(as.matrix(tbl.list[[i]]))
          , col=colour.list[[i]]
          , xlab="Individual #", ylab="Ancestry"
          , border=NA
          , las = 2
          , cex.names = 0.5
          , main = CV_error.df$V1[i]
  )
  
}

dev.off()



#### FINAL SINGLE FIGURE ####
## User set variables (for manual run of selected)
input_Q.FN   <- "05_admixture/pvit_all.4.Q"
metadata.FN  <- "02_input_data/additional_file_s1_sample_metadata_2024-01-13.txt"
input_fam.FN <- "05_admixture/pvit_all.fam"

# Read in data
tbl <- read.table(file = input_Q.FN)

barplot(t(as.matrix(tbl)), col=rainbow(5),
          xlab="Individual #", ylab="Ancestry"
        , border=NA
        )

# Read in sample order information
fam <- read.table(file = input_fam.FN)
head(fam) # second column is indivs, and they are in the same order of the admixture matrix

# Assign indiv names to the admixture matrix
rownames(tbl) <- fam$V2
head(tbl)
str(tbl)

# Add metadata for sorting
meta.df <- read.delim(file = metadata.FN, header = T)
head(meta.df)
meta.df <- meta.df[, c("sample_name", "collection_region", "lat_lon")]
head(meta.df)
meta.df <- separate(data = meta.df, col = "lat_lon", into = c("lat", "dir1", "lon", "dir2"), sep = " ", remove = T)
head(meta.df)

meta.df$sample_name <- gsub(pattern = ".fq.gz", replacement = "", x = meta.df$sample_name)
head(meta.df)

meta.df$coast <- NA
meta.df$coast[grep(pattern = "CAL|ORE|SOG|NBC", x = meta.df$collection_region)] <- "west"
meta.df$coast[grep(pattern = "CAL|ORE|SOG|NBC", x = meta.df$collection_region, invert = T)] <- "east"
head(meta.df)


meta.df <- meta.df[with(meta.df, order(meta.df$coast, meta.df$lat, decreasing = T)), ]


head(meta.df)
head(tbl)

meta.df <- meta.df[meta.df$sample_name %in% rownames(tbl), ]
dim(meta.df)
dim(tbl)

tbl$sample_name <- rownames(tbl)

data.df <- merge(x = meta.df, y = tbl, by = "sample_name", sort = F)
data.df
head(data.df)

tbl_ordered.df <- data.df[, grep(pattern = "V", x = colnames(data.df))]
rownames(tbl_ordered.df) <- data.df$sample_name


barplot(t(as.matrix(tbl_ordered.df)), col=rainbow(5),
        xlab="Individual #", ylab="Ancestry", border=NA
        , las = 2
        , cex.names = 0.5
)


# Set output filename based on input file

output_fig.FN <- gsub(pattern = ".Q", replacement = "_admixture_barplot.pdf", x = input_Q.FN)

# Set colours
setK <- as.numeric(str_sub(string = str_sub(string = input_Q.FN, start = -3), start = 1, end = 1))
# Note: could also get this from the tbl_ordered.df
k_colours <- c("darkorange", "darkorchid3","green", "yellow")

# SAVE OUT
pdf(file = output_fig.FN, width = 16, height = 7)
barplot(t(as.matrix(tbl_ordered.df)), col=k_colours,
        xlab="Individual #", ylab="Ancestry", border=NA
        , las = 2
        , cex.names = 0.5
        )

dev.off()










