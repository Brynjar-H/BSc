# Find DMRs
dmrs0 <- dmrFinder(BS_ctr.tstat, cutoff = c(-4.6, 4.6))
dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)

# Identify DMRs within the SERPINB5 region
serpinb5_region <- dmrs[dmrs$chr == "chr18" & dmrs$start >= 63474958 & dmrs$end <= 63507085, ]

# Extract and print the CpG site locations within the identified DMRs
if (nrow(serpinb5_region) > 0) {
  serpinb5_cpg_sites <- data.table(chrom = serpinb5_region$chr, 
                                   start = serpinb5_region$start, 
                                   end = serpinb5_region$end)
  print(serpinb5_cpg_sites)
} else {
  message("No DMRs found within the SERPINB5 region.")
}