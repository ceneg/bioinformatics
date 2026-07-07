#!/usr/bin/env bash
set -e
set -o pipefail

# Setup exit codes tracker
s1_status="PENDING"
s2_status="PENDING"
s2_r_status="PENDING"
s3_status="PENDING"
s4_r_status="PENDING"
s5_status="PENDING"

# Ensure we're in the repository root directory
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

echo "=========================================================="
echo "          STARTING EXTENDED WORKFLOW INTEGRATION TEST      "
echo "=========================================================="

# Create snapshot of baseline files
echo "Taking baseline snapshot of directory tree..."
find . -not -path '*/.pixi/*' -not -path '*/.git/*' | sort > baseline_files.txt

# Configure Pixi to run post-link scripts natively for future installations
pixi config set --local run-post-link-scripts insecure

# Install and instantiate all environments securely from the frozen lockfile
pixi install --frozen --all

# Manually execute post-link scripts for R packages if Pixi skipped them previously
echo "Checking and executing any pending R post-link scripts..."
shopt -s dotglob
for env in transkriptomika polimorfizmi diverziteta; do
    if [ -d ".pixi/envs/$env/bin" ]; then
        for script in .pixi/envs/$env/bin/*post-link.sh; do
            if [ -f "$script" ]; then
                echo "Running post-link script: $script"
                # Export PREFIX as CONDA_PREFIX for compatibility with Bioconductor scripts
                pixi run -e "$env" bash -c 'export PREFIX=$CONDA_PREFIX; bash "$1"' _ "$script" || echo "[WARNING] Post-link script failed. Continuing..."
            fi
        done
    fi
done
shopt -u dotglob

# Function to clean up safely and print the final report
cleanup_and_report() {
    # Disable strict error checking for the trap so we ensure everything prints
    set +e
    set +o pipefail

    # Extract error tails FIRST, before cleanup loop deletes the log files
    ERRORS=""
    if [[ "$s1_status" == FAILED* ]] && [ -f tmp_s1.log ]; then ERRORS+="\n--- Session 1 Error Log (Last 30 lines) ---\n$(tail -n 30 tmp_s1.log)\n"; fi
    if [[ "$s2_status" == FAILED* ]] && [ -f tmp_s2.log ]; then ERRORS+="\n--- Session 2 Error Log (Last 30 lines) ---\n$(tail -n 30 tmp_s2.log)\n"; fi
    if [[ "$s2_r_status" == FAILED* ]] && [ -f tmp_s2_r.log ]; then ERRORS+="\n--- Session 2 R Error Log (Last 30 lines) ---\n$(tail -n 30 tmp_s2_r.log)\n"; fi
    if [[ "$s3_status" == FAILED* ]] && [ -f tmp_s3.log ]; then ERRORS+="\n--- Session 3 Error Log (Last 30 lines) ---\n$(tail -n 30 tmp_s3.log)\n"; fi
    if [[ "$s4_r_status" == FAILED* ]] && [ -f tmp_s4.log ]; then ERRORS+="\n--- Session 4 Error Log (Last 30 lines) ---\n$(tail -n 30 tmp_s4.log)\n"; fi
    if [[ "$s5_status" == FAILED* ]] && [ -f tmp_s5.log ]; then ERRORS+="\n--- Session 5 Error Log (Last 30 lines) ---\n$(tail -n 30 tmp_s5.log)\n"; fi

    echo -e "\n=========================================================="
    echo "          CLEANING UP GENERATED WORKFLOW ARTIFACTS        "
    echo "=========================================================="

    if [ -f baseline_files.txt ]; then
        find . -not -path '*/.pixi/*' -not -path '*/.git/*' | sort > post_files.txt
        
        deleted_count=0
        # Find new files/folders created during the run
        while read -r path; do
            # Remove leading './'
            norm_path="${path#./}"
            
            if [ -z "$norm_path" ] || [ "$norm_path" = "." ] || [ "$norm_path" = "baseline_files.txt" ] || [ "$norm_path" = "post_files.txt" ]; then
                continue
            fi

            # Check whitelist of safe paths to delete
            is_safe=0
            if [[ "$norm_path" == spades_out* ]] || \
               [[ "$norm_path" == "genes.gff3" ]] || \
               [[ "$norm_path" == salmon_index* ]] || \
               [[ "$norm_path" == data/session2/CTL_rep1* ]] || \
               [[ "$norm_path" == data/session2/CTL_rep2* ]] || \
               [[ "$norm_path" == data/session2/CTL_rep3* ]] || \
               [[ "$norm_path" == data/session2/TRT_rep1* ]] || \
               [[ "$norm_path" == data/session2/TRT_rep2* ]] || \
               [[ "$norm_path" == data/session2/TRT_rep3* ]] || \
               [[ "$norm_path" == session2_plots* ]] || \
               [[ "$norm_path" == session4_plots* ]] || \
               [[ "$norm_path" == data/session2/deseq2_results.tsv ]] || \
               [[ "$norm_path" == data/session3/*.sorted.bam* ]] || \
               [[ "$norm_path" == data/session3/ecoli_ref_region.fasta.* ]] || \
               [[ "$norm_path" == data/session3/variants.vcf ]] || \
               [[ "$norm_path" == data/session3/vcf_stats.* ]] || \
               [[ "$norm_path" == data/session3/variants_filtered.recode.vcf ]] || \
               [[ "$norm_path" == rep-seqs.qza ]] || \
               [[ "$norm_path" == table.qza ]] || \
               [[ "$norm_path" == stats-dada2.qza ]] || \
               [[ "$norm_path" == base-transition-stats.qza ]] || \
               [[ "$norm_path" == aligned-rep-seqs.qza ]] || \
               [[ "$norm_path" == masked-aligned-rep-seqs.qza ]] || \
               [[ "$norm_path" == unrooted-tree.qza ]] || \
               [[ "$norm_path" == rooted-tree.qza ]] || \
               [[ "$norm_path" == taxonomy.qza ]] || \
               [[ "$norm_path" == core-metrics-results* ]] || \
               [[ "$norm_path" == taxa-bar-plots.qzv ]] || \
               [[ "$norm_path" == tmp_s*.log ]] || \
               [[ "$norm_path" == fastp_report.html ]] || \
               [[ "$norm_path" == fastp.json ]] || \
               [[ "$norm_path" == quast_out* ]] || \
               [[ "$norm_path" == "scratch_trans_list.txt" ]]; then
                is_safe=1
            fi

            if [ "$is_safe" -eq 1 ]; then
                if [ -d "$norm_path" ]; then
                    # Only remove directory if it's empty or we are removing spades_out/salmon_index/core-metrics-results/quast_out
                    if [[ "$norm_path" == "spades_out" ]] || [[ "$norm_path" == "salmon_index" ]] || [[ "$norm_path" == "core-metrics-results" ]] || [[ "$norm_path" == "quast_out" ]] || [[ "$norm_path" == "data/session2/"* ]] || [[ "$norm_path" == "session2_plots" ]] || [[ "$norm_path" == "session4_plots" ]]; then
                        rm -rf "$norm_path"
                        deleted_count=$((deleted_count + 1))
                    else
                        rmdir "$norm_path" 2>/dev/null && deleted_count=$((deleted_count + 1)) || true
                    fi
                elif [ -f "$norm_path" ] || [ -L "$norm_path" ]; then
                    rm -f "$norm_path"
                    deleted_count=$((deleted_count + 1))
                fi
            fi
        done < <(comm -13 baseline_files.txt post_files.txt | sort -r)
        rm -f baseline_files.txt post_files.txt
        echo "[CLEAN] Removed $deleted_count files and directories."
    fi

    rm -f tmp_s1.log tmp_s2.log tmp_s2_r.log tmp_s3.log tmp_s4.log tmp_s5.log

    # Print final test report
    echo -e "\n=========================================================="
    echo "                  EXTENDED TEST SUMMARY REPORT            "
    echo "=========================================================="
    echo "Session 1 (Genomics - SPAdes & Augustus):   $s1_status"
    echo "Session 2 (Transcriptomics - Salmon):       $s2_status"
    echo "Session 2 (R libraries & DESeq2 check):     $s2_r_status"
    echo "Session 3 (Polymorphisms - BWA & bcftools): $s3_status"
    echo "Session 4 (R libraries check):              $s4_r_status"
    echo "Session 5 (Microbiome - QIIME 2):           $s5_status"
    echo "=========================================================="

    if [ -n "$ERRORS" ]; then
        echo -e "\n=========================================================="
        echo "            FAILURE REASONS (LOG EXTRACTS)                "
        echo "=========================================================="
        echo -e "$ERRORS"
        echo "=========================================================="
    fi

    # Exit with failure code if any session failed
    if [ "$s1_status" = "SUCCESS" ] && \
       [ "$s2_status" = "SUCCESS" ] && \
       [ "$s2_r_status" = "SUCCESS" ] && \
       [ "$s3_status" = "SUCCESS" ] && \
       [ "$s4_r_status" = "SUCCESS" ] && \
       [ "$s5_status" = "SUCCESS" ]; then
        echo -e "\n[ALL PASSED] Extended workflow tests completed successfully!"
        exit 0
    else
        echo -e "\n[FAILED] One or more extended workflow tests failed!"
        exit 1
    fi
}

trap cleanup_and_report EXIT

# ----------------------------------------------------
# SESSION 1: De Novo Assembly & Gene Prediction
# ----------------------------------------------------
echo -e "\n--- [Running Session 1: Genomics] ---"
if (pixi run -e genomika fastp -i data/session1/yeast_chrom1_R1.fastq.gz -I data/session1/yeast_chrom1_R2.fastq.gz -h fastp_report.html && \
    pixi run -e genomika spades.py -1 data/session1/yeast_chrom1_R1.fastq.gz -2 data/session1/yeast_chrom1_R2.fastq.gz -o spades_out/ --only-assembler && \
    pixi run -e genomika quast spades_out/scaffolds.fasta -o quast_out && \
    pixi run -e genomika augustus --species=saccharomyces --gff3=on spades_out/scaffolds.fasta > genes.gff3) 2>&1 | tee tmp_s1.log; then
    
    # Validation
    if [ -s spades_out/scaffolds.fasta ] && [ -s genes.gff3 ] && [ -s fastp_report.html ] && [ -d quast_out ] && grep -q ">" spades_out/scaffolds.fasta && grep -q "augustus" genes.gff3; then
        s1_status="SUCCESS"
        echo "[OK] Session 1 tools and outputs successfully verified."
    else
        s1_status="FAILED (invalid output files)"
        echo "[ERROR] Session 1 output validation failed."
    fi
else
    s1_status="FAILED"
    echo "[ERROR] Session 1 execution failed."
fi

# ----------------------------------------------------
# SESSION 2: Transcript Quantification & DGE
# ----------------------------------------------------
echo -e "\n--- [Running Session 2: Transcriptomics (Salmon)] ---"
if (pixi run -e transkriptomika salmon index -t data/session2/escherichia_coli_transcriptome.fasta -i salmon_index/ && \
   (
     success=1
     for SAMPLE in CTL_rep1 CTL_rep2 CTL_rep3 TRT_rep1 TRT_rep2 TRT_rep3; do
         pixi run -e transkriptomika salmon quant -i salmon_index/ -l A -1 data/session2/${SAMPLE}_R1.fastq.gz -2 data/session2/${SAMPLE}_R2.fastq.gz -o data/session2/${SAMPLE} || success=0
     done
     [ $success -eq 1 ]
   )) 2>&1 | tee tmp_s2.log; then
    
    # Validation
    if [ -d salmon_index ] && [ -s data/session2/CTL_rep1/quant.sf ]; then
        s2_status="SUCCESS"
        echo "[OK] Session 2 Salmon quantification successfully verified."
    else
        s2_status="FAILED (missing quant.sf)"
        echo "[ERROR] Session 2 Salmon validation failed."
    fi
else
    s2_status="FAILED"
    echo "[ERROR] Session 2 Salmon execution failed."
fi

echo -e "\n--- [Running Session 2 R Modules & DESeq2 Functional Test] ---"
if pixi run -e transkriptomika Rscript -e '
message("Loading transkriptomika R libraries...")
libs <- c("tximport", "DESeq2", "ggplot2", "ggrepel", "pheatmap", "readr", "dplyr", "tibble", "RColorBrewer")
for (lib in libs) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    stop("Package ", lib, " is missing!")
  }
}
message("Testing DESeq2 workflow execution...")
library(tximport)
library(DESeq2)
library(readr)
library(ggplot2)

sample_ids <- c("CTL_rep1", "CTL_rep2", "CTL_rep3", "TRT_rep1", "TRT_rep2", "TRT_rep3")
conditions <- c("CTL", "CTL", "CTL", "TRT", "TRT", "TRT")
coldata <- data.frame(sample_id = sample_ids, condition = factor(conditions, levels = c("CTL", "TRT")), row.names = sample_ids)

quant_files <- file.path("data/session2", coldata$sample_id, "quant.sf")
names(quant_files) <- coldata$sample_id

tx_names <- read_tsv(quant_files[1], col_types = cols())$Name
tx2gene <- data.frame(TXNAME = tx_names, GENENAME = tx_names)
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps = TRUE)

dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)
dds <- DESeq(dds)

res_shrunk <- lfcShrink(dds, coef = "condition_TRT_vs_CTL", type = "normal")
res_df <- as.data.frame(res_shrunk)
res_df <- res_df[order(res_df$padj), ]
write.table(res_df, "data/session2/deseq2_results.tsv", sep = "\t", quote = FALSE, col.names = NA)

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
message("[OK] DESeq2 full workflow executed successfully.")
' 2>&1 | tee tmp_s2_r.log; then
    s2_r_status="SUCCESS"
else
    s2_r_status="FAILED"
    echo "[ERROR] Session 2 R validation failed."
fi

# ----------------------------------------------------
# SESSION 3: Read Mapping & Variant Calling
# ----------------------------------------------------
echo -e "\n--- [Running Session 3: Polymorphisms (BWA & bcftools)] ---"
if (pixi run -e polimorfizmi bwa index data/session3/ecoli_ref_region.fasta && \
   (
     success=1
     for IND in ind1 ind2 ind3 ind4 ind5; do
         pixi run -e polimorfizmi bwa mem data/session3/ecoli_ref_region.fasta data/session3/${IND}_R1.fastq.gz data/session3/${IND}_R2.fastq.gz | \
         pixi run -e polimorfizmi samtools view -bS - | \
         pixi run -e polimorfizmi samtools sort -o data/session3/${IND}.sorted.bam - || success=0
         pixi run -e polimorfizmi samtools index data/session3/${IND}.sorted.bam || success=0
     done
     [ $success -eq 1 ]
   ) && \
   pixi run -e polimorfizmi bcftools mpileup -f data/session3/ecoli_ref_region.fasta data/session3/*.sorted.bam | \
   pixi run -e polimorfizmi bcftools call -mv -o data/session3/variants.vcf && \
   pixi run -e polimorfizmi vcftools --vcf data/session3/variants.vcf --freq --out data/session3/vcf_stats && \
   pixi run -e polimorfizmi vcftools --vcf data/session3/variants.vcf --depth --out data/session3/vcf_stats && \
   pixi run -e polimorfizmi vcftools --vcf data/session3/variants.vcf --indv data/session3/ind1.sorted.bam --indv data/session3/ind2.sorted.bam --recode --out data/session3/variants_filtered) 2>&1 | tee tmp_s3.log; then
    
    # Validation
    if [ -s data/session3/variants.vcf ] && grep -q "#CHROM" data/session3/variants.vcf && [ -s data/session3/vcf_stats.frq ] && [ -s data/session3/vcf_stats.idepth ] && [ -s data/session3/variants_filtered.recode.vcf ]; then
        s3_status="SUCCESS"
        echo "[OK] Session 3 tools and outputs successfully verified."
    else
        s3_status="FAILED (invalid VCF)"
        echo "[ERROR] Session 3 variant calling validation failed."
    fi
else
    s3_status="FAILED"
    echo "[ERROR] Session 3 execution failed."
fi

# ----------------------------------------------------
# SESSION 4: Population Genetics (R checks)
# ----------------------------------------------------
echo -e "\n--- [Running Session 4: Population Structure R Workflow] ---"
if pixi run -e polimorfizmi Rscript -e '
message("Loading polimorfizmi R libraries...")
libs <- c("vcfR", "adegenet", "hierfstat", "ggplot2")
for (lib in libs) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    stop("Package ", lib, " is missing!")
  }
}
message("Testing Population Genetics workflow execution...")
library(vcfR)
library(adegenet)
library(hierfstat)
library(ggplot2)

vcf <- read.vcfR("data/session4/population_structure.vcf", verbose = FALSE)
metadata <- read.table("data/session4/sample_metadata.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

gen_data <- vcfR2genind(vcf)
pop(gen_data) <- factor(metadata$family)

stats <- basic.stats(gen_data)

Ho_avg <- mean(stats$Ho, na.rm = TRUE)
He_avg <- mean(stats$Hs, na.rm = TRUE)
Fis_avg <- mean(stats$Fis, na.rm = TRUE)

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
message("[OK] Polimorfizmi R workflow executed successfully.")
' 2>&1 | tee tmp_s4.log; then
    s4_r_status="SUCCESS"
else
    s4_r_status="FAILED"
    echo "[ERROR] Session 4 R validation failed."
fi

# ----------------------------------------------------
# SESSION 5: 16S Amplicon & Microbiome Diversity (QIIME 2)
# ----------------------------------------------------
echo -e "\n--- [Running Session 5: QIIME 2 Workflows] ---"
if (pixi run -e qiime2 qiime dada2 denoise-single \
     --i-demultiplexed-seqs data/session5/amplicon_demux.qza \
     --p-trim-left 0 \
     --p-trunc-len 120 \
     --o-representative-sequences rep-seqs.qza \
     --o-table table.qza \
     --o-denoising-stats stats-dada2.qza \
     --o-base-transition-stats base-transition-stats.qza && \
   pixi run -e qiime2 qiime phylogeny align-to-tree-mafft-fasttree \
     --i-sequences rep-seqs.qza \
     --o-alignment aligned-rep-seqs.qza \
     --o-masked-alignment masked-aligned-rep-seqs.qza \
     --o-tree unrooted-tree.qza \
     --o-rooted-tree rooted-tree.qza && \
   pixi run -e qiime2 qiime feature-classifier classify-sklearn \
     --i-classifier data/session5/silva_classifier_excerpt.qza \
     --i-reads rep-seqs.qza \
     --o-classification taxonomy.qza && \
   pixi run -e qiime2 qiime diversity core-metrics-phylogenetic \
     --i-phylogeny rooted-tree.qza \
     --i-table table.qza \
     --p-sampling-depth 1000 \
     --m-metadata-file data/session5/sample_metadata.tsv \
     --output-dir core-metrics-results && \
   pixi run -e qiime2 qiime taxa barplot \
     --i-table table.qza \
     --i-taxonomy taxonomy.qza \
     --m-metadata-file data/session5/sample_metadata.tsv \
     --o-visualization taxa-bar-plots.qzv) 2>&1 | tee tmp_s5.log; then
    
    # Validation
    if [ -s taxonomy.qza ] && [ -d core-metrics-results ] && [ -s taxa-bar-plots.qzv ]; then
        s5_status="SUCCESS"
        echo "[OK] Session 5 QIIME 2 tools and outputs successfully verified."
    else
        s5_status="FAILED (missing output QZA/QZV files)"
        echo "[ERROR] Session 5 output validation failed."
    fi
else
    s5_status="FAILED"
    echo "[ERROR] Session 5 execution failed."
fi
