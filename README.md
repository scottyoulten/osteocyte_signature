# Defining the osteocyte transcriptome signature
A code walk-through for the definition of genes actively expressed in the osteocyte transcriptome, gene expression enrichment in osteocytes and definition of the osteocyte transcriptome signature as described in the publication Osteocyte Transcriptome Mapping Identifies a Molecular Landscape Controlling Skeletal Homeostasis and Susceptibility to Skeletal Disease (https://www.biorxiv.org/content/10.1101/2020.04.20.051409v2).

## data
The data directory contains preprocessed expression data (counts and FPKM) used in the analysis. Processed data files (BAM) and links to raw data (fastq) hosted on the European Nucleotide Archive are available at ArrayExpress with following accession IDs:  
E-MTAB-5532 [https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5532/] - Profiling the transcriptome of the adult mouse skeleton: RNA-seq of osteocytes from tibiae, femora, humeri and calvariae  
E-MTAB-5533 [https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5533/] - Identifying genes enriched for expression in osteocytes: transcriptome sequencing mouse bone with and without bone marrow  
E-MTAB-7447 [https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7447/] - The osteocyte transcriptome during post-natal skeletal maturation in male and female mice

## scripts
The scripts directory contains R code used to caluclate gene activity information, gene expression enrichment in osteocytes and define the osteocyte transcriptome signature.

## analysis
To run this analysis, first clone this repository git clone https://github.com/syouligan/osteocyte_signature.git.
From the downloaded directory run the analysis scripts sequentially:  
1.Gene_activity_threshold  
2.Gene_expression_enrichment  
3.Osteocyte_transcriptome_signature

## output
These scripts will produce a geneActivity.csv file, with columns containing annotation information, gene activity, gene expression enrichment and osteocyte transcriptome signature information.
