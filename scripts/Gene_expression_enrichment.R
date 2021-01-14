#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# This script calculates the gene expression enrichment in osteocyte enriched samples 
# --------------------------------------------------------------------------

# Load libraries
library("limma")
library("edgeR")

# Load data files (make sure working directory is osteocyte_signature directory) and remove trailing decimal from gene ids
osteocyte_enrichment <- read.csv('data/Osteocyte_enriched_RSEM_gene_counts.csv', row.names = 1, header = TRUE)
rownames(osteocyte_enrichment) <- sub("(.*?)\\..*", "\\1", rownames(osteocyte_enrichment))

# Load gene annotation info with gene activity information
gene_annotations <- read.csv('Gene_annotation_w_activity.csv', row.names = 1, header = TRUE)

# Use limma to get an estimate of log-fold change between osteocyte enriched bone samples (clean) and marrow containing bone samples (marrow).
# --------------------------------------------------------------------------

# Subset counts to genes activty in either osteocyte-enriched or marrow-containing bone samples
active_filter <- gene_annotations$Clean_activity == 5 | gene_annotations$Marrow_activity == 5
active_osteocyte_enrichment <- osteocyte_enrichment[active_filter,]

# Make a sample table for use in experimental design
sampleNames <- colnames(active_osteocyte_enrichment)
sampleTypes <- factor(c(gsub( "\\_[0-9]*$", "", sampleNames)))
sampleTable <- data.frame(sampleName = sampleNames, type = sampleTypes)

# Make dataframe with gene symbol and description for each gene
genes <- data.frame(rownames(active_osteocyte_enrichment), row.names = 1)
idx <- match(rownames(genes), rownames(gene_annotations))
genes$GeneSymbol <- gene_annotations$GeneSymbol[ idx ]
genes$Description <- gene_annotations$Description[ idx ]

# Make DGEList object using edgeR
activeDGE <- DGEList(counts = as.matrix(active_osteocyte_enrichment), samples = sampleTable, group = sampleTable$type, genes = genes)

# Voom counts, fit linear model and calculate log-fold-change with confidence intervals
type <- factor(sampleTable$type, levels = c("Clean", "Marrow"))
design <- model.matrix(~0 + type)
rownames(design) <- colnames(active_osteocyte_enrichment)
activeDGE <- voom(activeDGE, design)
fit <- lmFit(activeDGE, design)
contrast.matrix <- makeContrasts(typeClean-typeMarrow, levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit2 <- treat(fit, lfc=0)
geneEnrichment <- topTreat(fit2, coef=1, number=Inf, confint = TRUE, lfc = 0) # Used to calculate log-fold-change and CIs

# Write gene annotation csv with gene enrichment information
idx <- match(rownames(gene_annotations), rownames(geneEnrichment))
gene_annotations$LFC <- geneEnrichment$logFC [idx]
gene_annotations$LFC_CI_Upper <- geneEnrichment$CI.R [idx]
gene_annotations$LFC_CI_Lower <- geneEnrichment$CI.L [idx]

# Generate Gaussian Mixture Model (GMM) of osteocyte enrichment and calculate osteocyte enrichment threshold
# --------------------------------------------------------------------------

# Create a vector of gene log2FC between osteocyte-enriched and marrow containing bone samples (removing NA values)
geneLFC <- geneEnrichment$logFC
geneLFC <- geneLFC[!is.na(geneLFC)]

# Estimate optimum number of components using the Bayesian Information Criterion (BIC) and unequal variance.
BIC <- mclustBIC(geneLFC, model = "V")
# pdf("Gene_enrichment_BIC_plot.pdf")
# plot(BIC)
# dev.off()
# capture.output(BIC, file = "Gene_enrichment_BIC_summary.txt")

# Perform k-means assignment
logFC.kmeans <- kmeans(geneLFC, 4) # 4 chosen based on optimum BIC model
logFC.kmeans.cluster <- logFC.kmeans$cluster
logFC.df <- data_frame(x = geneLFC, cluster = logFC.kmeans.cluster)
logFC.summary.df <- logFC.df %>%
  group_by(cluster) %>%
  summarize(mu = mean(x), std = sd(x), size = n())
logFC.summary.df <- logFC.summary.df %>%
  mutate(alpha = size / sum(size))

# Inititate EM with parameters from k-means to get GMM
LFCmixmdl <- normalmixEM(geneLFC, lambda = logFC.summary.df$alpha, mu = logFC.summary.df$mu, sigma = logFC.summary.df$std, maxit = 10000)

# Orders modules mean and SD based on descending module mean
means <- LFCmixmdl$mu
sds <- LFCmixmdl$sigma
sds <- sds[order(-means)]
means <- means[order(-means)]

# Calculate threshold based on mean + 2SD of the module with the second highest mean (the highest is the most osteocyte enriched, ie the one we're trying to isolate)
osteocyte_threshold <- means[2] + 2*sds[2]

# Annotate osteocyte enriched threshold column on gene activity dataframe
gene_annotations$Above_threshold <- (!is.na(gene_annotations$LFC_CI_Lower)) & gene_annotations$LFC_CI_Lower > osteocyte_threshold
sum(gene_annotations$Above_threshold)
write.csv(gene_annotations, file = "Gene_annotation_w_activity_enrichment.csv")
