# Load necessary libraries
library(DESeq2)
library(readr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)

# Load DESeq2 results for Kallisto
load("path/to/kallisto/deseq2/dds.RData")
kallisto_dds <- dds

# Load DESeq2 results for RSEM
load("path/to/rsem/deseq2/dds.RData")
rsem_dds <- dds

# Load DESeq2 results for Salmon
load("path/to/salmon/deseq2/dds.RData")
salmon_dds <- dds

# Create a condition variable in colData and convert to factor for all datasets
colData(kallisto_dds)$condition <- factor(paste(colData(kallisto_dds)$Group1, colData(kallisto_dds)$Group2, sep = "_"))
colData(rsem_dds)$condition <- factor(paste(colData(rsem_dds)$Group1, colData(rsem_dds)$Group2, sep = "_"))
colData(salmon_dds)$condition <- factor(paste(colData(salmon_dds)$Group1, colData(salmon_dds)$Group2, sep = "_"))

# Convert Group3 to factor for all datasets
colData(kallisto_dds)$Group3 <- factor(colData(kallisto_dds)$Group3)
colData(rsem_dds)$Group3 <- factor(colData(rsem_dds)$Group3)
colData(salmon_dds)$Group3 <- factor(colData(salmon_dds)$Group3)

# Update design formula to include the condition variable
design(kallisto_dds) <- ~ condition
design(rsem_dds) <- ~ condition
design(salmon_dds) <- ~ condition

# Run DESeq2 analysis for all datasets with the updated design
kallisto_dds <- DESeq(kallisto_dds)
rsem_dds <- DESeq(rsem_dds)
salmon_dds <- DESeq(salmon_dds)

# Subset to common genes
common_genes <- Reduce(intersect, list(rownames(counts(kallisto_dds)), rownames(counts(rsem_dds)), rownames(counts(salmon_dds))))
kallisto_dds <- kallisto_dds[common_genes, ]
rsem_dds <- rsem_dds[common_genes, ]
salmon_dds <- salmon_dds[common_genes, ]

# List of Ensembl IDs of interest, including additional PGR transcripts
ensembl_ids_of_interest <- c(
  'ENSG00000091831', 'ENSG00000082175', 'ENSG00000141736', 'ENSG00000012048', 
  'ENSG00000139618', 'ENSG00000148773', 'ENST00000325455', 'ENST00000534013', 
  'ENST00000533207', 'ENST00000263463', 'ENST00000526300', 'ENST00000534780', 
  'ENST00000528960', 'ENST00000530764', 'ENST00000619228', 'ENST00000632634'
)

# Function to get normalized counts for Ensembl IDs of interest
get_normalized_counts <- function(dds, ensembl_ids) {
  norm_counts <- counts(dds, normalized = TRUE)
  norm_counts <- norm_counts[rownames(norm_counts) %in% ensembl_ids, ]
  return(norm_counts)
}

# Get normalized counts for Ensembl IDs of interest for each dataset
kallisto_counts <- get_normalized_counts(kallisto_dds, ensembl_ids_of_interest)
rsem_counts <- get_normalized_counts(rsem_dds, ensembl_ids_of_interest)
salmon_counts <- get_normalized_counts(salmon_dds, ensembl_ids_of_interest)

# Combine the counts into one data frame
combined_counts <- data.frame(
  Kallisto_7d = rowMeans(kallisto_counts[, colData(kallisto_dds)$condition == 'D492_7d']),
  Kallisto_14d = rowMeans(kallisto_counts[, colData(kallisto_dds)$condition == 'D492_14d']),
  Kallisto_21d = rowMeans(kallisto_counts[, colData(kallisto_dds)$condition == 'D492_21d']),
  RSEM_7d = rowMeans(rsem_counts[, colData(rsem_dds)$condition == 'D492_7d']),
  RSEM_14d = rowMeans(rsem_counts[, colData(rsem_dds)$condition == 'D492_14d']),
  RSEM_21d = rowMeans(rsem_counts[, colData(rsem_dds)$condition == 'D492_21d']),
  Salmon_7d = rowMeans(salmon_counts[, colData(salmon_dds)$condition == 'D492_7d']),
  Salmon_14d = rowMeans(salmon_counts[, colData(salmon_dds)$condition == 'D492_14d']),
  Salmon_21d = rowMeans(salmon_counts[, colData(salmon_dds)$condition == 'D492_21d'])
)

# Only keep rows that are present in the combined_counts
ensembl_ids_present <- rownames(combined_counts)
combined_counts <- combined_counts[ensembl_ids_present, ]

# Set row names to Ensembl IDs
rownames(combined_counts) <- ensembl_ids_present

# Map Ensembl IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids_present, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Replace missing symbols with original Ensembl IDs
gene_symbols[is.na(gene_symbols)] <- ensembl_ids_present[is.na(gene_symbols)]

# Set the row names of combined_counts to gene symbols
rownames(combined_counts) <- gene_symbols

# Exclude rows with all zero values
combined_counts <- combined_counts[rowSums(combined_counts != 0) > 0, ]

# Check the combined counts before generating the heatmap
print(combined_counts)

# Generate the heatmap
pheatmap(
  combined_counts,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Heatmap Comparing Kallisto, RSEM, and Salmon for Breast Cancer Genes",
  fontsize_row = 10, # Adjust fontsize for row names
  fontsize_col = 10, # Adjust fontsize for column names
  annotation_names_col = TRUE,
  legend = TRUE
)

head(combined_counts)
