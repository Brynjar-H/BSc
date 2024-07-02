# Load necessary libraries
library(DESeq2)
library(readr)
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(apeglm)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)

# Load DESeq2 results for Kallisto
load("path/to/kallisto/deseq2/dds.RData")
kallisto_dds <- dds

# Load DESeq2 results for Salmon
load("path/to/salmon/deseq2/dds.RData")
salmon_dds <- dds

# Create a condition variable in colData and convert to factor for both datasets
colData(kallisto_dds)$condition <- factor(paste(colData(kallisto_dds)$Group1, colData(kallisto_dds)$Group2, sep = "_"))
colData(salmon_dds)$condition <- factor(paste(colData(salmon_dds)$Group1, colData(salmon_dds)$Group2, sep = "_"))

# Convert Group3 to factor for both datasets
colData(kallisto_dds)$Group3 <- factor(colData(kallisto_dds)$Group3)
colData(salmon_dds)$Group3 <- factor(colData(salmon_dds)$Group3)

# Update design formula to include the condition variable
design(kallisto_dds) <- ~ condition
design(salmon_dds) <- ~ condition

# Run DESeq2 analysis for both datasets with the updated design
kallisto_dds <- DESeq(kallisto_dds)
salmon_dds <- DESeq(salmon_dds)

# Subset to common genes
common_genes <- intersect(rownames(counts(kallisto_dds)), rownames(counts(salmon_dds)))
kallisto_dds <- kallisto_dds[common_genes, ]
salmon_dds <- salmon_dds[common_genes, ]

# Function to average replicates
average_replicates <- function(dds, condition) {
  norm_counts <- counts(dds, normalized = TRUE)
  samples <- colnames(norm_counts)[dds$condition == condition]
  avg_counts <- rowMeans(norm_counts[, samples])
  return(avg_counts)
}

# Get averaged counts for each condition and aligner
avg_counts_kallisto_7d <- average_replicates(kallisto_dds, "Condition_7d")
avg_counts_kallisto_14d <- average_replicates(kallisto_dds, "Condition_14d")
avg_counts_kallisto_21d <- average_replicates(kallisto_dds, "Condition_21d")

avg_counts_salmon_7d <- average_replicates(salmon_dds, "Condition_7d")
avg_counts_salmon_14d <- average_replicates(salmon_dds, "Condition_14d")
avg_counts_salmon_21d <- average_replicates(salmon_dds, "Condition_21d")

# Combine the averaged counts
combined_counts <- data.frame(
  Kallisto_7d = avg_counts_kallisto_7d,
  Kallisto_14d = avg_counts_kallisto_14d,
  Kallisto_21d = avg_counts_kallisto_21d,
  Salmon_7d = avg_counts_salmon_7d,
  Salmon_14d = avg_counts_salmon_14d,
  Salmon_21d = avg_counts_salmon_21d
)

# Calculate row variances
row_vars <- apply(combined_counts, 1, var)

# Select top variable genes
top_genes <- head(order(row_vars, decreasing = TRUE), 25)
combined_counts_top <- combined_counts[top_genes, ]

# Map Ensembl IDs to gene symbols, retaining original Ensembl IDs where symbols are missing
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = rownames(combined_counts_top),
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Replace missing symbols with original Ensembl IDs
gene_symbols[is.na(gene_symbols)] <- rownames(combined_counts_top)[is.na(gene_symbols)]

# Replace row names with gene symbols
rownames(combined_counts_top) <- gene_symbols

# Generate the heatmap without annotations for Aligner and Condition
pheatmap(
  combined_counts_top,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Heatmap Comparing Kallisto and Salmon for Different Time Points",
  fontsize_row = 5, # Adjust fontsize for row names
  fontsize_col = 8, # Adjust fontsize for column names
  annotation_names_col = TRUE,
  legend = TRUE
)

