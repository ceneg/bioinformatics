# 1. Nalaganje knjižnic in uvoz podatkov
library(vcfR)
library(adegenet)
library(hierfstat)
library(poppr)
library(ggplot2)

# Branje VCF datoteke in metapodatkov
vcf <- read.vcfR("data/session4/population_structure.vcf", verbose = FALSE)
metadata <- read.table("data/session4/sample_metadata.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Pretvorba VCF v genind objekt za analizo v adegenet
gen_data <- vcfR2genind(vcf)
pop(gen_data) <- factor(metadata$family)

# 2. Izračun osnovnih genskih statistik
stats <- basic.stats(gen_data)

Ho_avg <- mean(stats$Ho, na.rm = TRUE)
He_avg <- mean(stats$Hs, na.rm = TRUE)
Fis_avg <- mean(stats$Fis, na.rm = TRUE)

print(paste("Povprečna opažena heterozigotnost (Ho):", round(Ho_avg, 4)))
print(paste("Povprečna pričakovana heterozigotnost (He):", round(He_avg, 4)))
print(paste("Povprečni koeficient inbridinga (Fis):", round(Fis_avg, 4)))

# 3. Analiza glavnih komponent (PCA)
# Zamenjava manjkajočih podatkov s povprečno frekvenco alela
X <- tab(gen_data, NA.method = "mean")
pca <- dudi.pca(X, cent = TRUE, scale = FALSE, scannf = FALSE, nf = 2)

pca_df <- data.frame(
  Sample = indNames(gen_data),
  PC1 = pca$li$Axis1,
  PC2 = pca$li$Axis2,
  Family = pop(gen_data)
)

dir.create("session4_plots", showWarnings = FALSE)
png("session4_plots/pca_plot.png", width = 800, height = 600)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Family, label = Sample)) +
  geom_point(size = 4, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA of Genomic Variants")
dev.off()

# 4. Minimalno vpeto omrežje (MSN)
# Osebke povežemo glede na gensko razdaljo med njimi. Debelina in dolžina
# povezav odražata sorodnost, barva pa pripadnost družini.
msn <- poppr.msn(gen_data, diss.dist(gen_data), showplot = FALSE)

png("session4_plots/msn_plot.png", width = 800, height = 600)
set.seed(42)
plot_poppr_msn(gen_data, msn, palette = rainbow, gadj = 15, nodescale = 8)
dev.off()
