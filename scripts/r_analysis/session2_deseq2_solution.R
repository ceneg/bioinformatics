#!/usr/bin/env Rscript
# ==============================================================================
# DESeq2 Differential Gene Expression Analysis Reference Solution
#
# This script performs transcript quantification import via tximport and runs
# the DESeq2 package to identify genes differentially expressed between Control (CTL)
# and Treated (TRT) conditions in E. coli.
# ==============================================================================

# Determine paths relative to the script's own location
script_dir <- tryCatch({
  dirname(rstudioapi::getActiveDocumentContext()$path)
}, error = function(e) {
  # Fallback to system frame of file or getwd()
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  if (length(script.name) > 0) dirname(script.name) else getwd()
})

# Define data paths relative to script location
repo_root <- file.path(script_dir, "..", "..")
data_dir <- file.path(repo_root, "data", "session2")
plots_dir <- file.path(data_dir, "plots")

# Ensure plots directory exists
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

# 1. Load required libraries
library(tximport)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(dplyr)
library(readr)
library(RColorBrewer)

message("Libraries loaded successfully.")

# 2. Define experimental metadata (coldata)
sample_ids <- c("CTL_rep1", "CTL_rep2", "CTL_rep3", "TRT_rep1", "TRT_rep2", "TRT_rep3")
conditions <- c("CTL", "CTL", "CTL", "TRT", "TRT", "TRT")

coldata <- data.frame(
  sample_id = sample_ids,
  condition = factor(conditions, levels = c("CTL", "TRT")),
  row.names = sample_ids
)

# 3. Locate and import Salmon quant.sf output files
# Students will run Salmon in their vaje to generate these folders inside data/session2/
quant_files <- file.path(data_dir, coldata$sample_id, "quant.sf")
names(quant_files) <- coldata$sample_id

# Verify files exist before running (avoids R runtime crashes)
if (!all(file.exists(quant_files))) {
  stop("Missing Salmon quant.sf output files. Make sure to run the Salmon quantification step first.")
}

# We use the transcript-to-gene mapping file (simulated data uses transcript IDs as gene IDs directly)
# Create a dummy tx2gene map linking each transcript to itself
tx_names <- read_tsv(quant_files[1], col_types = cols())$Name
tx2gene <- data.frame(TXNAME = tx_names, GENENAME = tx_names)

message("Importing Salmon quantifications with tximport...")
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene)

# 4. Construct DESeqDataSet and run differential expression
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)

message("Running DESeq2 pipeline...")
dds <- DESeq(dds)

# 5. Extract differential expression results
res <- results(dds, contrast = c("condition", "TRT", "CTL"), alpha = 0.05)

# Shrink fold changes using 'ashr' or 'apeglm' for better visualization
# Note: 'apeglm' requires coefficient specification rather than contrast
res_shrunk <- lfcShrink(dds, coef = "condition_TRT_vs_CTL", type = "apeglm")

# Save results table to TSV
res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column(var = "gene_id") %>%
  arrange(padj)

write_tsv(res_df, file.path(data_dir, "deseq2_results.tsv"))
message("Results saved to deseq2_results.tsv.")

# 6. Generate and save diagnostic visualizations
# 6.1. PCA Plot
message("Generating PCA plot...")
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA - Principal Component Analysis")

ggsave(file.path(plots_dir, "pca_plot.png"), plot = pca_plot, width = 6, height = 5)

# 6.2. Volcano Plot
message("Generating Volcano plot...")
res_df <- res_df %>%
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) > 1.0)

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("black", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  theme_minimal() +
  labs(title = "Volcano Plot (TRT vs CTL)", x = "log2 Fold Change", y = "-log10 adjusted P-value") +
  theme(legend.position = "none")

# Label the top 10 most significant genes
top_10_genes <- head(res_df, 10)
volcano_plot <- volcano_plot +
  geom_text_repel(data = top_10_genes, aes(label = gene_id), size = 3, max.overlaps = 15)

ggsave(file.path(plots_dir, "volcano_plot.png"), plot = volcano_plot, width = 7, height = 6)

# 6.3. MA Plot
message("Generating MA plot...")
png(file.path(plots_dir, "ma_plot.png"), width = 700, height = 600)
plotMA(res_shrunk, ylim = c(-5, 5), main = "MA Plot (Shrunk LFC)")
dev.off()

# 6.4. Heatmap of top 30 differentially expressed genes
message("Generating Heatmap...")
top_30_sig <- head(res_df$gene_id, 30)
vst_counts <- assay(vsd)[top_30_sig, ]

# Annotate heatmap columns
annotation_col <- data.frame(Condition = coldata$condition)
rownames(annotation_col) <- rownames(coldata)

png(file.path(plots_dir, "heatmap_top30.png"), width = 800, height = 800)
pheatmap(
  vst_counts,
  scale = "row",
  clustering_distance_cols = "correlation",
  annotation_col = annotation_col,
  show_rownames = TRUE,
  main = "Top 30 Differentially Expressed Genes"
)
dev.off()

# 7. Print summary statistics
sig_genes_count <- sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1.0, na.rm = TRUE)
message("======================================================================")
message("Analysis completed successfully.")
message("Total significant DE genes (padj < 0.05 & |log2FC| > 1): ", sig_genes_count)
message("======================================================================")
