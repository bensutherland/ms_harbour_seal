# Admixture analysis of harbour seal data
# B. Sutherland, 2024-01-13

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

# Load libraries
library("tidyr")


## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)
rm(current.path)

## Set options
options(scipen=99999999)

## User-set variables
metadata.FN  <- "02_input_data/additional_file_s1_sample_metadata_2024-01-13.txt"


# Define dataset names (based on filename prefix)
all_dataset <- "pvit_all"
pac_dataset <- "pvit_pac"
atl_dataset <- "pvit_atl"

dataset <- all_dataset # User: choose which dataset to use

## Based on dataset, complete filenames
input_fam.FN <- paste0("05_admixture/", dataset, ".fam")
CV_error.FN <- paste0("05_admixture/", dataset, "_CV_error.txt")

## Create list of the numbers of clusters (i.e., range of K) to display per datatype 
k.list <- list()
k.list[["pvit_all"]] <- 6
k.list[["pvit_pac"]] <- 5
k.list[["pvit_atl"]] <- 5

## Obtain list of relevant dataset filenames 
dataset_files <- list.files(path = "05_admixture/", pattern = dataset, full.names = T)
dataset_files <- dataset_files[grep(pattern = ".Q", x = dataset_files)]       # only keep the Q files
dataset_files <- head(sort(dataset_files), n = k.list[[dataset]])             # only keep the defined range of K as set above

#### Prepare CV error data ####
# Read in CV error
CV_error.df <- read.table(file = CV_error.FN, sep = ",")
CV_error.df
CV_error.df <- head(CV_error.df, n = k.list[[dataset]])



### Prepare metadata ####
# Read in metadata for sorting order based on GPS or grouping
meta.df <- read.delim(file = metadata.FN, header = T)
meta.df <- meta.df[, c("sample_name", "collection_region", "lat_lon")] # retain only specific cols
meta.df <- separate(data = meta.df, col = "lat_lon", into = c("lat", "dir1", "lon", "dir2"), sep = " ", remove = T) # sep lat / lon
# note: expect warning, due to NAs
head(meta.df)

# Format indiv string in metadata
meta.df$sample_name <- gsub(pattern = ".fq.gz", replacement = "", x = meta.df$sample_name)
head(meta.df)

# Add variable for coast for sorting, based on names of pops
meta.df$coast <- NA
meta.df$coast[grep(pattern = "CAL|ORE|SOG|NBC", x = meta.df$collection_region)] <- "west"
meta.df$coast[grep(pattern = "CAL|ORE|SOG|NBC", x = meta.df$collection_region, invert = T)] <- "east"
head(meta.df)

# Order metadata by coast then lat
meta.df <- meta.df[with(meta.df, order(meta.df$coast, meta.df$lat, decreasing = T)), ]


### Read in sample order in results ####
# Read in order of samples (note: consistent order for different Ks)
fam <- read.table(file = input_fam.FN)
head(fam) # second column is indivs, and they are in the same order of the admixture matrix


##### Q data per K ####
# Read in Q ppn data
tbl <- NULL; tbl.list <- list(); tbl_ordered.df <- NULL; data.df <- NULL
for(i in 1:length(dataset_files)){
  
  #print(i)
  
  # Read in data
  tbl <- read.table(file = dataset_files[[i]])
  
  # Assign indiv names to the admixture matrix
  rownames(tbl) <- fam$V2
  
  # The two df to be merged
  head(meta.df)
  head(tbl)
  
  # Only keep the present sample rows in metadata
  meta.df <- meta.df[meta.df$sample_name %in% rownames(tbl), ]
  dim(meta.df)
  dim(tbl)
  
  # Add sample name as a vector of the df
  tbl$sample_name <- rownames(tbl)
  head(tbl)
  
  # Merge the metadata and results by sample name
  data.df <- merge(x = meta.df, y = tbl, by = "sample_name", sort = F)
  dim(data.df)
  head(data.df)
  
  # Drop all info except the result proportions
  tbl_ordered.df <- data.df[, grep(pattern = "V|sample_name", x = colnames(data.df))]
  rownames(tbl_ordered.df) <- tbl_ordered.df$sample_name
  head(tbl_ordered.df)
  
  if(i == 1){ 
    
    # save order to bring it back
    sample_order.vec <- tbl_ordered.df$sample_name
    
    }
  
  tbl_ordered.df <- tbl_ordered.df[, grep(pattern = "V", x = colnames(tbl_ordered.df))]
  
  tbl_ordered.df <- as.data.frame(tbl_ordered.df)
  head(tbl_ordered.df)
  
  if(i == 1){
    
    # add back in the order
    rownames(tbl_ordered.df) <- sample_order.vec
    
  }
  head(tbl_ordered.df)
  
  
  # Save it to the list
  tbl.list[[i]] <- tbl_ordered.df
  
}

str(tbl.list)


# Set cols
k1_col <- c("darkorchid3")
k2_col <- c("darkorchid3", "darkorange")
k3_col <- c("yellow", "darkorchid3", "darkorange")
k4_col <- c("darkorange", "darkorchid3",  "green", "yellow")
k5_col <- c("yellow", "green", "darkorchid3",  "black", "darkorange")

k6_col <- c("darkorange", "black",  "yellow", "darkorchid3","pink", "green" )

colour.list <- list()
colour.list[[1]] <- k1_col
colour.list[[2]] <- k2_col
colour.list[[3]] <- k3_col
colour.list[[4]] <- k4_col
colour.list[[5]] <- k5_col
colour.list[[6]] <- k6_col

#k_colours <- c("darkorange", "darkorchid3","green", "yellow")


# Plot
pdf(file = paste0("05_admixture/", dataset, "_multi_k_admixture_plot.pdf")
    , width = 7.5, height = 12)
par(mfrow = c(k.list[[dataset]], 1)) # number of plots match k
for(i in 1:length(tbl.list)){
  
  barplot(t(as.matrix(tbl.list[[i]]))
          , col=colour.list[[i]]
          , xlab="Individual #"
          , ylab="Ancestry"
          , border=NA
          , las = 2
          , cex.names = 0.5
          , main = CV_error.df$V1[i]
  )
  
}

dev.off()


#### FINAL SINGLE FIGURE ####

## Plot k = 4
pdf(file = paste0("05_admixture/optimal_k_admixture_plot.pdf")
    , width = 10.5, height = 5)
par(mfrow = c(1,1))
barplot(t(as.matrix(tbl.list[[4]]))
        , col=colour.list[[4]]
        , xlab="Individual #"
        , ylab="Ancestry"
        , border=NA
        , las = 2
        , cex.names = 0.5
        #, main = CV_error.df$V1[4]
)
dev.off()

