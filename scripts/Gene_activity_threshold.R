#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# This script calculates the per-sample gene-activity thresholds for each sequencing sample and uses it to define active genes per sample. 
# The basis of this methodology was originally described in the paper: Hart, Traver, et al. "Finding the active genes in deep RNA-seq gene expression studies." BMC genomics 14.1 (2013): 778.
# --------------------------------------------------------------------------

# Load data files (make sure working directory is osteocyte_signature directory) and remove trailing decimal from gene ids
osteocyte_enrichment <- read.csv('data/Osteocyte_enriched_RSEM_gene_FPKM.csv', row.names = 1, header = TRUE)
rownames(osteocyte_enrichment) <- sub("(.*?)\\..*", "\\1", rownames(osteocyte_enrichment))

bone_comparison <- read.csv('data/Bone_comparison_RSEM_gene_FPKM.csv', row.names = 1, header = TRUE)
rownames(bone_comparison) <- sub("(.*?)\\..*", "\\1", rownames(bone_comparison))

skeletal_maturation <- read.csv('data/Skeletal_maturation_RSEM_gene_FPKM.csv', row.names = 1, header = TRUE)
rownames(skeletal_maturation) <- sub("(.*?)\\..*", "\\1", rownames(skeletal_maturation))

# Load gene annotation info
gene_annotations <- read.csv('data/Gene_annotation.csv', row.names = 1, header = TRUE)

# Calculate activity threshold for each sample in each dataset
# --------------------------------------------------------------------------

# Make list of expression dataframes
FPKM_list <- list(osteocyte_enrichment, bone_comparison, skeletal_maturation)

# Function to caluclate activity threshold for each sample (column) in each dataframe 
df_activity <- function(x){
  activity_summary <- data.frame()
  for (i in 1:ncol(x)) {
    #Calculate kernal density of log2FPKM values (removing genes with FPKM = 0)
    sample_expression <- x[ , i, drop = FALSE]
    sample_expression <- sample_expression[! sample_expression[ , 1] == 0, , drop = FALSE]
    sample_name <- colnames(sample_expression)
    se_num_vec <- log2(sample_expression)[ , 1]
    c <- density(se_num_vec, bw = "nrd")

    # Retrieve the log2FPKM value corresponding to the maximum kernel density estimate (MKDE) and calculate the standard deviation, where where U is the mean of all log2(FPKM) values > MKDE
    MKDE <- subset(c$x, c$y == max(c$y)) 
    U <- mean(subset(se_num_vec, se_num_vec > MKDE))
    sd <- (U-MKDE)*(sqrt(pi/2))
    
    # Calculate zFPKM scores using a Gaussian distribution with parameters (MKDE, sd).
    zFPKM <- data.frame((se_num_vec - MKDE)/sd, row.names = rownames(sample_expression))
    colnames(zFPKM) <- "zFPKM_values"

    # Hart et al. indicate that a reasonable expression cutoff describing the active genes in a cell would be the point where the ratio of active to repressed promoters drops below 1 (zFPKM = -2.82 +/- 0.22).
    # We chose the conservative end of this range preferring to be stringent in our definition of gene activity, setting a threshold of -2.6.
    # This threshold is used to subset active genes. 
    active <- zFPKM[zFPKM$zFPKM_values >= -2.6, , drop = FALSE]
    
    # Identify genes with minimum zFPKM threshold
    min_active_genes <- rownames(active[active$zFPKM_values == min(active$zFPKM_values), , drop = FALSE])
    
    # Retrieve FPKM value of genes with minimum zFPKM threshold
    FPKM_threshold <- unique(sample_expression[min_active_genes,])
    activity_summary <- rbind(activity_summary, c(MKDE, sd, FPKM_threshold))
  }
  
  colnames(activity_summary) <- c("Mean of distribution", "SD of distribution", "FPKM_threshold")
  rownames(activity_summary) <- colnames(x)
  
return(activity_summary)
}

# Create list of sample thresholds for each dataset
threshold_list <- lapply(FPKM_list, df_activity)

# Identify genes >= threshold in each sample for each dataset
# --------------------------------------------------------------------------

# Make list equal to the number of dataset in FPKM_list object
dataset_list <- as.list(1:length(FPKM_list))

# Function to generate boolean matrix of gene activity in each sample
above_threshold <- function(x){
  active_matrix <- FPKM_list[[x]] >= threshold_list[[x]]$FPKM_threshold
  return(active_matrix)
}

# List of activity matricies
active_matricies <- lapply(dataset_list, above_threshold)

# Identify genes actively expressed in each sample type (>= sample activity threshold in all replicates)
# --------------------------------------------------------------------------

# Function to generate sample information table based on sample name (note naming convention of Type_replicate)
sampleInfo <- function(x){
  sampleName <- colnames(x)
  sampleType <- factor(c(gsub( "\\_[0-9]*$", "", sampleName)))
  sampleRep <- factor(c(gsub( "^.*?_","", sampleName)))
  sampleTable <- data.frame(sampleName, sampleType, sampleRep)
  return(sampleTable)
  }

# Create list of sample info tables for each dataset
sampleTables <- lapply(FPKM_list, sampleInfo)

# Update gene annotation with mean expression and gene activity in each sample type
for(x in 1:length(FPKM_list)) {
  
  for(i in unique(sampleTables[[x]][,"sampleType"])) {
    FPKM_df <- FPKM_list[[x]][,grepl(i, x = colnames(FPKM_list[[x]]))]
    FPKM_mean <- rowMeans(FPKM_df)
    idx <- match(rownames(gene_annotations), names(FPKM_mean))
    gene_annotations[, paste0(i, "_mean_FPKM")] <- FPKM_mean [idx]
    
    activity_df <- active_matricies[[x]][,grepl(i, x = colnames(active_matricies[[x]]))]
    Activity_sum <- rowSums(activity_df)
    idx <- match(rownames(gene_annotations), names(Activity_sum))
    gene_annotations[, paste0(i, "_activity")] <- Activity_sum [idx]
  }
  
}

# Write gene annotation csv with gene activity information
write.csv(gene_annotations, "Gene_annotation_w_activity.csv")


