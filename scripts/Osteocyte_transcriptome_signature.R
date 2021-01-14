#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# This script defines the osteocyte transcriptome signature 
# --------------------------------------------------------------------------

# Load gene annotation info with gene enrichment information
gene_annotations <- read.csv('Gene_annotation_w_activity_enrichment.csv', row.names = 1, header = TRUE)

# Identify genes active in any of the sample types in the osteocyte_enrichment, bone_comparison or skeletal_maturation datasets (i.e. in the osteocyte transcriptome)
# --------------------------------------------------------------------------

osteocyte_enrichment_filter <- gene_annotations$Clean_activity == 5 | gene_annotations$Marrow_activity == 5
bone_comparison_filter <- gene_annotations$Tibia_activity == 8 | gene_annotations$Femur_activity == 8 | gene_annotations$Humerus_activity == 8
skeletal_maturation_filter <- gene_annotations$Male4_activity == 5 | gene_annotations$Male10_activity == 5 | gene_annotations$Male16_activity == 5 | gene_annotations$Male26_activity == 5 | gene_annotations$Female4_activity == 5 | gene_annotations$Female10_activity == 5 | gene_annotations$Female16_activity == 5 | gene_annotations$Female26_activity == 5
gene_annotations$Osteocyte_transcriptome <- osteocyte_enrichment_filter | bone_comparison_filter | skeletal_maturation_filter
sum(gene_annotations$Osteocyte_transcriptome) 

# Define osteocyte transcriptome signature as genes in the osteocyte transcriptome that are enriched in osteocytes and not contaminating tissues (as defined by Ayturk, Ugur M., et al. "An RNA‐seq protocol to identify mRNA expression changes in mouse diaphyseal bone: applications in mice with bone property altering Lrp5 mutations." Journal of Bone and Mineral Research 28.10 (2013): 2081-2093.)
gene_annotations$Osteocyte_transcriptome_signature <- gene_annotations$Osteocyte_transcriptome & gene_annotations$Above_threshold & (!gene_annotations$Ayturk_contaminant)
sum(gene_annotations$Osteocyte_transcriptome_signature)
write.csv(gene_annotations, file = "Gene_annotation_w_activity_enrichment_signature.csv")
