#!/usr/bin/env Rscript
# ==============================================================================
# Population Genetics VCF Analysis Reference Solution
#
# This script reads a multi-sample VCF, merges metadata, performs population genetic
# calculations (Heterozygosity, Fis), constructs a Minimum Spanning Network (MSN),
# and conducts a PCA to identify population/pedigree structure.
# ==============================================================================

# Determine paths relative to the script's own location
script_dir <- tryCatch({
  dirname(rstudioapi::getActiveDocumentContext()$path)
}, error = function(e) {
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  if (length(script.name) > 0) dirname(script.name) else getwd()
})

# Define data paths
repo_root <- file.path(script_dir, "..", "..")
data_dir <- file.path(repo_root, "data", "session4")
plots_dir <- file.path(data_dir, "plots")

# Ensure plots directory exists
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

# 1. Load required libraries
library(vcfR)
library(poppr)
library(adegenet)
library(ggplot2)
library(ape)
library(magrittr)

message("Libraries loaded successfully.")

# 2. Read VCF and metadata
vcf_file <- file.path(data_dir, "population_structure.vcf")
meta_file <- file.path(data_dir, "sample_metadata.tsv")

if (!file.exists(vcf_file) || !file.exists(meta_file)) {
  stop("Missing Session 4 input data files. Please run the data generation script first.")
}

message("Reading VCF file...")
vcf <- read.vcfR(vcf_file, verbose = FALSE)

message("Reading metadata...")
metadata <- read.table(meta_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 3. Convert vcfR to genind object
message("Converting VCF to genind...")
gen_data <- vcfR2genind(vcf)

# Ensure sample ordering aligns with metadata
metadata <- metadata[match(indNames(gen_data), metadata$sample_id), ]

# Assign population factor based on family
pop(gen_data) <- factor(metadata$family)

# 4. Calculate heterozygosity and stats per locus
message("Calculating basic genetic statistics...")
# Load hierfstat package internally since it is a dependency of adegenet/poppr
library(hierfstat)
stats <- basic.stats(gen_data)

# Extract and print average observed/expected heterozygosity and Fis
Ho_avg <- mean(stats$Ho, na.rm = TRUE)
He_avg <- mean(stats$Hs, na.rm = TRUE)
Fis_avg <- mean(stats$Fis, na.rm = TRUE)

message("==================================================")
message("Mean Observed Heterozygosity (Ho): ", round(Ho_avg, 4))
message("Mean Expected Heterozygosity (He): ", round(He_avg, 4))
message("Mean Inbreeding Coefficient (Fis): ", round(Fis_avg, 4))
message("==================================================")

# 5. Compute Bruvo's distance and Minimum Spanning Network (MSN)
message("Constructing Minimum Spanning Network (MSN)...")
# Convert genind to genclone for MSN plotting in poppr
clone_data <- as.snpclone(gen_data)

# Since we don't have repeat length definitions for SNP loci, we set it to 1
png(file.path(plots_dir, "msn_plot.png"), width = 800, height = 800)
msn <- poppr.msn(
  clone_data,
  dist = bitwise.dist(clone_data),
  showplot = TRUE,
  gscale = FALSE,
  palette = colorRampPalette(brewer.pal(4, "Set1")),
  vertex.label = indNames(clone_data)
)
dev.off()

# 6. Run Principal Component Analysis (PCA) on allele frequencies
message("Performing PCA on genotypes...")
# Replace missing data with mean allele frequency (VCF is simulated complete, but it's good practice)
X <- tab(gen_data, NA.method = "mean")

# Run PCA
pca <- dudi.pca(X, cent = TRUE, scale = FALSE, scannf = FALSE, nf = 2)

# Create a data frame for plotting with ggplot2
pca_df <- data.frame(
  Sample = indNames(gen_data),
  PC1 = pca$li$Axis1,
  PC2 = pca$li$Axis2,
  Family = pop(gen_data)
)

# Plot PCA
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Family, label = Sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3, max.overlaps = 15) +
  theme_minimal() +
  labs(
    title = "PCA of Genomic Variants (Diploid)",
    x = paste0("PC1 (", round(pca$eig[1] / sum(pca$eig) * 100, 2), "%)"),
    y = paste0("PC2 (", round(pca$eig[2] / sum(pca$eig) * 100, 2), "%)")
  ) +
  scale_color_brewer(palette = "Set1")

ggsave(file.path(plots_dir, "pca_plot.png"), plot = pca_plot, width = 7, height = 6)
message("Plots successfully written to data/session4/plots/")
