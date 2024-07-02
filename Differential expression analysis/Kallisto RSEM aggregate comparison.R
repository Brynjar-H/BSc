# Load Necessary Libraries
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)

# Set directories
kallisto_dir <- "path/to/kallisto/directory/"
rsem_dir <- "path/to/rsem/directory/"

# Load Kallisto TPM data
kallisto_tpm <- read_tsv(file.path(kallisto_dir, "kallisto.merged.gene_tpm.tsv"))

# Load RSEM TPM data
rsem_tpm <- read_tsv(file.path(rsem_dir, "rsem.merged.gene_tpm.tsv"))

# Reshape Kallisto data from wide to long format
kallisto_long <- kallisto_tpm %>%
  pivot_longer(cols = -c(gene_id, gene_name), names_to = "sample", values_to = "TPM_kallisto")

# Reshape RSEM data from wide to long format
rsem_long <- rsem_tpm %>%
  pivot_longer(cols = -c(gene_id, `transcript_id(s)`), names_to = "sample", values_to = "TPM_rsem")

# Aggregate TPM values by gene_id and gene_name for Kallisto
kallisto_agg <- kallisto_long %>%
  group_by(gene_id, gene_name) %>%
  summarise(TPM_kallisto = mean(TPM_kallisto, na.rm = TRUE))

# Aggregate TPM values by gene_id for RSEM
rsem_agg <- rsem_long %>%
  group_by(gene_id) %>%
  summarise(TPM_rsem = mean(TPM_rsem, na.rm = TRUE))

# Merge the aggregated data
merged_agg_tpm <- inner_join(kallisto_agg, rsem_agg, by = "gene_id")

# Calculate Pearson correlation coefficient
cor_test <- cor.test(merged_agg_tpm$TPM_kallisto, merged_agg_tpm$TPM_rsem)
pearson_correlation <- cor_test$estimate
p_value <- cor_test$p.value

# Plot the aggregated TPM values
ggplot(merged_agg_tpm, aes(x = TPM_kallisto, y = TPM_rsem)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Comparison of Aggregated Kallisto and RSEM TPM",
       x = "Kallisto TPM",
       y = "RSEM TPM") +
  theme_minimal() +
  annotate("text", x = Inf, y = Inf, label = paste("r =", round(pearson_correlation, 3)), 
           hjust = 1.7, vjust = 1.1, size = 3.5, color = "black")

# Filter genes by TPM RSEM = 0 and Kallisto > 100
filtered_genes_rsem_zero <- merged_agg_tpm %>%
  filter(TPM_rsem == 0 & TPM_kallisto > 100)

# Inspect the filtered genes
print(filtered_genes_rsem_zero)

# Filter genes by TPM Kallisto = 0 and RSEM > 100
filtered_genes_kallisto_zero <- merged_agg_tpm %>%
  filter(TPM_rsem > 100 & TPM_kallisto == 0)

# Inspect the filtered genes
print(filtered_genes_kallisto_zero)

# Count genes with TPM = 0 in RSEM but non-zero in Kallisto
count_rsem_zero_kallisto_nonzero <- merged_agg_tpm %>%
  filter(TPM_rsem == 0 & TPM_kallisto > 0) %>%
  nrow()

# Count genes with TPM = 0 in Kallisto but non-zero in RSEM
count_kallisto_zero_rsem_nonzero <- merged_agg_tpm %>%
  filter(TPM_kallisto == 0 & TPM_rsem > 0) %>%
  nrow()

# Print the counts
cat("Number of genes with TPM = 0 in RSEM but non-zero in Kallisto:", count_rsem_zero_kallisto_nonzero, "\n")
cat("Number of genes with TPM = 0 in Kallisto but non-zero in RSEM:", count_kallisto_zero_rsem_nonzero, "\n")

# Identify outliers where RSEM TPM is significantly higher than Kallisto TPM
outliers <- merged_agg_tpm %>%
  filter(TPM_rsem > TPM_kallisto + 5000) # You can adjust this factor as needed

# Inspect the outliers
print(outliers)

# Filter genes with TPM >= 100 in Kallisto but < 1 in RSEM
kallisto_high_rsem_low <- merged_agg_tpm %>%
  filter(TPM_kallisto >= 100 & TPM_rsem < 1)

# Filter genes with TPM >= 100 in RSEM but < 1 in Kallisto
rsem_high_kallisto_low <- merged_agg_tpm %>%
  filter(TPM_rsem >= 100 & TPM_kallisto < 1)

# Create a summary table
summary_table <- bind_rows(
  kallisto_high_rsem_low %>% mutate(Status = "High in Kallisto, Low in RSEM"),
  rsem_high_kallisto_low %>% mutate(Status = "High in RSEM, Low in Kallisto")
)

# Print the summary table
print(summary_table)
