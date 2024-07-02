# Kaan Okay - HI  aug 2023
# Snaevar Sigurdsson okt 2023
# Brynjar Halldórsson 1 jun 2024

# Load required libraries
library(data.table)
library(bsseq)
library(tools) # for resaveRdaFiles below
library(BiocParallel)
library(bsseqData)

# Set working directory for input files
input_dir <- "path/to/input_directory"
setwd(input_dir)

# Read BED files
Ctr1 <- fread(file.path(input_dir, "control_1.bed"), nThread = 10)
Ctr2 <- fread(file.path(input_dir, "control_2.bed"), nThread = 10)
test1 <- fread(file.path(input_dir, "test_1.bed"), nThread = 10)
test2 <- fread(file.path(input_dir, "test_2.bed"), nThread = 10)

# Add "chr" prefix to chromosome column
Ctr1$V1 <- paste0("chr", Ctr1$V1)
Ctr2$V1 <- paste0("chr", Ctr2$V1)
test1$V1 <- paste0("chr", test1$V1)
test2$V1 <- paste0("chr", test2$V1)

# Add column names
colnames(Ctr1) <- colnames(Ctr2) <- colnames(test1) <- colnames(test2) <- c(
  "chrom", "start", "end", "name", "score", "strand", "tstart", "tend",
  "color", "coverage", "freq", "canon", "mod", "filt"
)

# Create a new column combining mod and canon
Ctr1$mod_and_canon <- Ctr1$mod + Ctr1$canon
Ctr2$mod_and_canon <- Ctr2$mod + Ctr2$canon
test1$mod_and_canon <- test1$mod + test1$canon
test2$mod_and_canon <- test2$mod + test2$canon

# Filter out rows where mod is zero
Ctr1 <- Ctr1[Ctr1$mod > 0, ]
Ctr2 <- Ctr2[Ctr2$mod > 0, ]
test1 <- test1[test1$mod > 0, ]
test2 <- test2[test2$mod > 0, ]

# Identify common chromosomes across all datasets
all_chromosomes <- Reduce(intersect, list(Ctr1$chrom, Ctr2$chrom, test1$chrom, test2$chrom))

# Function to create BSseq object and standardize chromosomes
create_bsseq_object <- function(data, sample_name, common_chromosomes) {
  data <- data[data$chrom %in% common_chromosomes, ]
  tmp_f <- data[data$strand == "+", ]
  tmp_r <- data[data$strand == "-", ]
  BSseq_f <- BSseq(
    chr = as.character(tmp_f$chrom), pos = tmp_f$start + 1,
    M = matrix(tmp_f$mod), Cov = matrix(tmp_f$mod_and_canon),
    sampleNames = paste0(sample_name, "f")
  )
  BSseq_r <- BSseq(
    chr = as.character(tmp_r$chrom), pos = tmp_r$start + 1,
    M = matrix(tmp_r$mod), Cov = matrix(tmp_r$mod_and_canon),
    sampleNames = paste0(sample_name, "r")
  )
  BSseq_combined <- combine(BSseq_f, BSseq_r)
  BSseq_combined <- collapseBSseq(BSseq_combined, group = c(sample_name, sample_name), type = "integer")
  return(BSseq_combined)
}

# Create BSseq objects with standardized chromosomes
Ctr1_BSseq <- create_bsseq_object(Ctr1, "Ctr1", all_chromosomes)
Ctr2_BSseq <- create_bsseq_object(Ctr2, "Ctr2", all_chromosomes)
test1_BSseq <- create_bsseq_object(test1, "test1", all_chromosomes)
test2_BSseq <- create_bsseq_object(test2, "test2", all_chromosomes)

# Combine all samples into one object
BS_ctr <- combine(Ctr1_BSseq, Ctr2_BSseq, test1_BSseq, test2_BSseq)
pData(BS_ctr)$type <- c("control", "control", "test", "test")
BS_ctr <- collapseBSseq(BS_ctr, group = c("Ctr1", "Ctr2", "test1", "test2"))

# Smooth the data
BSc.fit <- BSmooth(BS_ctr, BPPARAM = MulticoreParam(workers = 22, progressbar = TRUE), verbose = TRUE)

# Set working directory for output files
output_dir <- "path/to/output_directory"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
setwd(output_dir)

# Save the BSmooth object
save(BSc.fit, file = file.path(output_dir, "BSfit.rda"))
load(file.path(output_dir, "BSfit.rda"))

# Plotting and calculating statistics
BS.cov <- getCoverage(BSc.fit)
keepLoci.ex <- which(rowSums(BS.cov[, BSc.fit$type == "test"] >= 2) >= 2 &
                       rowSums(BS.cov[, BSc.fit$type == "control"] >= 2) >= 2)
BSc.fit <- BSc.fit[keepLoci.ex,]

# Compute t-statistics
BS_ctr.tstat <- BSmooth.tstat(BSc.fit, 
                              group1 = c("Ctr1", "Ctr2"),
                              group2 = c("test1", "test2"), 
                              estimate.var = "group2",
                              local.correct = TRUE,
                              verbose = TRUE)

# Find differentially methylated regions (DMRs)
dmrs0 <- dmrFinder(BS_ctr.tstat, cutoff = c(-4.6, 4.6))
dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
nrow(dmrs)

# Save the DMRs to a file
write.table(dmrs, file = "DMRS.bed", row.names = FALSE, sep = "\t")

# Plotting regions
pData <- pData(BSc.fit)
pData$col <- rep(c("red", "blue"), each = 2)
pData(BSc.fit) <- pData

# Define the region for plotting a DMR (e.g., the first DMR in the list)
if (nrow(dmrs) > 0) {
  pdf("plot_region_DMR.pdf")
  plotRegion(BSc.fit, dmrs[1,], extend = 2000, addRegions = dmrs)
  dev.off()
}

# Define the THOR region based on the article
THOR_region <- GRanges(seqnames = "chr5", ranges = IRanges(start = 1295228, end = 1295991))

# Plot the THOR region
pdf("plot_THOR_region.pdf")
plotRegion(BSc.fit, THOR_region, extend = 0)
dev.off()
