#!/usr/bin/env Rscript
# ==============================================================================
# DESeq2 Differential Gene Expression Analysis (Student Template)
#
# Dopolnite manjkajoče dele kode (označene s TODO bloki) za izvedbo analize.
# ==============================================================================

# Določitev poti glede na lokacijo skripta
script_dir <- tryCatch({
  dirname(rstudioapi::getActiveDocumentContext()$path)
}, error = function(e) {
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  if (length(script.name) > 0) dirname(script.name) else getwd()
})

repo_root <- file.path(script_dir, "..", "..")
data_dir <- file.path(repo_root, "data", "session2")
plots_dir <- file.path(data_dir, "plots")

if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

# ============================================================
# TODO: Nalaganje bioinformatskih R paketov
#
# Naloga: Naložite naslednje knjižnice s funkcijo library():
#         tximport, DESeq2, ggplot2, ggrepel, pheatmap, dplyr, readr, RColorBrewer
#
# Biološki kontekst: Analiza RNA-Seq zahteva uvoz abundanc (tximport),
# statistično modeliranje (DESeq2) ter vizualizacijo in urejanje podatkov.
#
# Hint: Uporabite npr. library(DESeq2)
# ============================================================
# Zamenjajte spodnjo vrstico s svojo kodo:
stop("Preden uporabite ta skript, morate dopolniti manjkajočo kodo za nalaganje knjižnic.")

# ============================================================
# TODO: Priprava metapodatkov vzorcev (coldata)
#
# Naloga: Ustvarite data.frame z imenom 'coldata'. Vsebovati mora 
#         stolpec 'condition' s faktorjem, ki ima nivoja 'CTL' in 'TRT' 
#         za 6 vzorcev (3x CTL in 3x TRT). Vrstice morajo imeti imena vzorcev.
#
# Biološki kontekst: DESeq2 potrebuje informacije o poskusnem načrtu
# (kateri vzorec pripada kontrolni skupini in kateri skupini s toplotnim šokom).
# ============================================================
# Zamenjajte spodnjo vrstico s svojo kodo:
stop("Preden uporabite ta skript, morate dopolniti manjkajočo kodo za pripravo metapodatkov vzorcev.")


quant_files <- file.path(data_dir, coldata$sample_id, "quant.sf")
names(quant_files) <- coldata$sample_id

if (!all(file.exists(quant_files))) {
  stop("Manjkajo izhodne datoteke Salmon quant.sf. Prepričajte se, da ste najprej zagnali Salmon kvantifikacijo.")
}

# Pripravimo preprost tx2gene zemljevid
tx_names <- read_tsv(quant_files[1], col_types = cols())$Name
tx2gene <- data.frame(TXNAME = tx_names, GENENAME = tx_names)


# ============================================================
# TODO: Uvoz podatkov o številu odčitkov (tximport)
#
# Naloga: Uporabite funkcijo tximport() za uvoz rezultatov Salmona.
#         Shrani rezultat v spremenljivko 'txi'.
#
# Biološki kontekst: Funkcija združi podatke o izražanju transkriptov na
# nivoju posameznih genov, kar omogoča korektno statistično obravnavo.
# ============================================================
# Zamenjajte spodnjo vrstico s svojo kodo:
stop("Preden uporabite ta skript, morate dopolniti manjkajočo kodo za uvoz Salmon rezultatov.")


# ============================================================
# TODO: Inicializacija in zagon DESeq2 analize
#
# Naloga: Ustvarite DESeqDataSet objekt s funkcijo DESeqDataSetFromTximport()
#         in nato poženite analizo s funkcijo DESeq().
#
# Biološki kontekst: Orodje izvede normalizacijo velikosti knjižnic, 
# oceni disperzijo za vsak gen in prilagodi statistični model (Negative Binomial).
# ============================================================
# Zamenjajte spodnjo vrstico s svojo kodo:
stop("Preden uporabite ta skript, morate dopolniti manjkajočo kodo za izvedbo statistične analize DESeq2.")


# Izvoz rezultatov
res <- results(dds, contrast = c("condition", "TRT", "CTL"), alpha = 0.05)
res_shrunk <- lfcShrink(dds, coef = "condition_TRT_vs_CTL", type = "apeglm")

res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column(var = "gene_id") %>%
  arrange(padj)

write_tsv(res_df, file.path(data_dir, "deseq2_results.tsv"))

# ============================================================
# TODO: Izris PCA grafa
#
# Naloga: Transformirajte podatke s funkcijo vst() in narišite PCA graf 
#         z uporabo funkcije plotPCA(). Shranite graf z ggsave() v plots_dir.
#
# Biološki kontekst: Analiza glavnih komponent (PCA) omogoča vizualno preverjanje 
# bioloških ponovitev in iskanje morebitnih osamelcev (outliers) v poskusu.
# ============================================================
# Zamenjajte spodnjo vrstico s svojo kodo:
stop("Preden uporabite ta skript, morate dopolniti manjkajočo kodo za izris PCA grafa.")

# Izris ostalih grafov (vulkan, MA, toplotni diagram) je že pripravljen spodaj.
# ...
message("Skript se je uspešno zaključil.")
