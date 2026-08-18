# 1. Nalaganje knjižnic in priprava metapodatkov
library(tximport)
library(DESeq2)
library(readr)

sample_ids <- c("CTL_rep1", "CTL_rep2", "CTL_rep3", "TRT_rep1", "TRT_rep2", "TRT_rep3")
conditions <- c("CTL", "CTL", "CTL", "TRT", "TRT", "TRT")
coldata <- data.frame(sample_id = sample_ids, condition = factor(conditions, levels = c("CTL", "TRT")), row.names = sample_ids)

# Priprava poti do datotek programa Salmon
quant_files <- file.path("data/session2", coldata$sample_id, "quant.sf")
names(quant_files) <- coldata$sample_id

# Uvoz podatkov s paketom tximport (združevanje transkriptov v gene)
tx_names <- read_tsv(quant_files[1], col_types = cols())$Name
tx2gene <- data.frame(TXNAME = tx_names, GENENAME = tx_names)
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps = TRUE)

# 2. Statistično modeliranje z DESeq2
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)
dds <- DESeq(dds)

# Izpis rezultatov in shranjevanje v datoteko
res_shrunk <- lfcShrink(dds, coef = "condition_TRT_vs_CTL", type = "normal")
res_df <- as.data.frame(res_shrunk)
res_df <- res_df[order(res_df$padj), ]
write.table(res_df, "data/session2/deseq2_results.tsv", sep = "\t", quote = FALSE, col.names = NA)

# 3. Vizualizacija - Volcano plot
library(ggplot2)
res_plot <- res_df
res_plot$threshold <- res_plot$padj < 0.05 & abs(res_plot$log2FoldChange) > 1.0

dir.create("session2_plots", showWarnings = FALSE)
png("session2_plots/volcano_plot.png", width = 800, height = 600)
ggplot(res_plot, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot (TRT vs CTL)")
dev.off()

# 4. Analiza glavnih komponent (PCA) vzorcev
# Podatke najprej stabiliziramo glede variance (VST), da visoko izraženi geni
# ne prevladajo nad ostalimi.
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pca_var <- round(100 * attr(pca_data, "percentVar"))

png("session2_plots/pca_plot.png", width = 800, height = 600)
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 4, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA vzorcev (VST)",
       x = paste0("PC1: ", pca_var[1], " % variance"),
       y = paste0("PC2: ", pca_var[2], " % variance"))
dev.off()
