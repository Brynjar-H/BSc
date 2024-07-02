# Load necessary libraries
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(data.table)

# Set directories
kallisto_dir <- "path/to/kallisto/directory/"
rsem_dir <- "path/to/rsem/directory/"

# Load Kallisto TPM data
kallisto_tpm <- fread(file.path(kallisto_dir, "kallisto.merged.gene_tpm.tsv"))

# Load RSEM TPM data
rsem_tpm <- fread(file.path(rsem_dir, "rsem.merged.gene_tpm.tsv"))

# Reshape Kallisto data from wide to long format
kallisto_long <- melt(kallisto_tpm, id.vars = c("gene_id", "gene_name"), variable.name = "sample", value.name = "TPM_kallisto")

# Reshape RSEM data from wide to long format
rsem_long <- melt(rsem_tpm, id.vars = c("gene_id", "transcript_id(s)"), variable.name = "sample", value.name = "TPM_rsem")

# Merge the reshaped data
merged_tpm <- merge(kallisto_long, rsem_long, by = c("gene_id", "sample"))

# Calculate Pearson correlation coefficient
cor_test <- cor.test(merged_tpm$TPM_kallisto, merged_tpm$TPM_rsem)
pearson_correlation <- cor_test$estimate
p_value <- cor_test$p.value

# Plot the TPM values
ggplot(merged_tpm, aes(x = TPM_kallisto, y = TPM_rsem)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Comparison of Kallisto and RSEM TPM",
       x = "Kallisto TPM",
       y = "RSEM TPM") +
  theme_minimal() +
  annotate("text", x = Inf, y = Inf, label = paste("r =", round(pearson_correlation, 3)), 
           hjust = 1.7, vjust = 1.1, size = 3.5, color = "black")
