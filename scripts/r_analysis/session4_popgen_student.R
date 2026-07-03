#!/usr/bin/env Rscript
# ==============================================================================
# Population Genetics VCF Analysis (Student Template)
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
data_dir <- file.path(repo_root, "data", "session4")
plots_dir <- file.path(data_dir, "plots")

if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

library(vcfR)
library(poppr)
library(adegenet)
library(ggplot2)
library(ape)
library(magrittr)

vcf_file <- file.path(data_dir, "population_structure.vcf")
meta_file <- file.path(data_dir, "sample_metadata.tsv")

if (!file.exists(vcf_file) || !file.exists(meta_file)) {
  stop("Manjkajo vhodne datoteke za 4. vajo. Prepričajte se, da ste najprej zagnali skript za simulacijo podatkov.")
}

vcf <- read.vcfR(vcf_file, verbose = FALSE)
metadata <- read.table(meta_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ============================================================
# TODO: Pretvorba iz vcfR v genind format
#
# Naloga: Uporabite ustrezno funkcijo iz paketa vcfR za pretvorbo
#         objekta 'vcf' v objekt 'gen_data' razreda genind.
#
# Biološki kontekst: Genind format je standardni R format za delo z
# genotipskimi podatki v paketih adegenet in poppr.
#
# Hint: Poglejte funkcijo vcfR2genind().
# ============================================================
# Zamenjajte spodnjo vrstico s svojo kodo:
stop("Preden uporabite ta skript, morate dopolniti manjkajočo kodo za pretvorbo VCF v genind.")


# Uskladitev metapodatkov z vrstnim redom osebkov v genind
metadata <- metadata[match(indNames(gen_data), metadata$sample_id), ]

# ============================================================
# TODO: Dodelitev populacijske pripadnosti
#
# Naloga: Dodelite populacijske oznake objektu 'gen_data'. Uporabite
#         stolpec 'family' iz metapodatkov kot faktor.
#
# Biološki kontekst: Pripadnost populaciji ali družini določa skupine
# za analize F-statistik in vizualizacije na grafih PCA in filogenetskih drevesih.
#
# Hint: Uporabite funkcijo pop(gen_data) <- ...
# ============================================================
# Zamenjajte spodnjo vrstico s svojo kodo:
stop("Preden uporabite ta skript, morate dopolniti manjkajočo kodo za dodelitev populacij.")


# ============================================================
# TODO: Izračun osnovnih genetskih statistik
#
# Naloga: Z uporabo paketa hierfstat izračunajte Ho, Hs (He) in Fis 
#         per locus na objektu 'gen_data'. Izpišite povprečne vrednosti.
#
# Biološki kontekst: Opažena (Ho) in pričakovana (He) heterozigotnost sta
# osnovni meri genetske pestrosti. Fis meri stopnjo inbridinga.
#
# Hint: Uporabite basic.stats() iz paketa hierfstat.
# ============================================================
# Zamenjajte spodnjo vrstico s svojo kodo:
stop("Preden uporabite ta skript, morate dopolniti manjkajočo kodo za izračun heterozigotnosti.")


# Preostali del za izris MSN in PCA se nahaja v rešitvi.
# ...
message("Skript se je uspešno zaključil.")
