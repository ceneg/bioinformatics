#!/usr/bin/env bash

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

# Function to clean up safely
cleanup_generated_files() {
    echo -e "\n=========================================================="
    echo "          CLEANING UP GENERATED WORKFLOW ARTIFACTS        "
    echo "=========================================================="

    if [ -f baseline_files.txt ]; then
        find . -not -path '*/.pixi/*' -not -path '*/.git/*' | sort > post_files.txt
        
        # Find new files/folders created during the run
        comm -13 baseline_files.txt post_files.txt | sort -r | while read -r path; do
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
               [[ "$norm_path" == data/session3/*.sorted.bam* ]] || \
               [[ "$norm_path" == data/session3/ecoli_ref_region.fasta.* ]] || \
               [[ "$norm_path" == data/session3/variants.vcf ]] || \
               [[ "$norm_path" == taxonomy.qza ]] || \
               [[ "$norm_path" == core-metrics-results* ]] || \
               [[ "$norm_path" == taxa-bar-plots.qzv ]] || \
               [[ "$norm_path" == "scratch_trans_list.txt" ]]; then
                is_safe=1
            fi

            if [ "$is_safe" -eq 1 ]; then
                if [ -d "$norm_path" ]; then
                    # Only remove directory if it's empty or we are removing spades_out/salmon_index/core-metrics-results
                    if [[ "$norm_path" == "spades_out" ]] || [[ "$norm_path" == "salmon_index" ]] || [[ "$norm_path" == "core-metrics-results" ]] || [[ "$norm_path" == "data/session2/"* ]]; then
                        rm -rf "$norm_path"
                        echo "[CLEAN] Removed directory: $norm_path"
                    else
                        rmdir "$norm_path" 2>/dev/null && echo "[CLEAN] Removed empty directory: $norm_path"
                    fi
                elif [ -f "$norm_path" ] || [ -L "$norm_path" ]; then
                    rm -f "$norm_path"
                    echo "[CLEAN] Removed file: $norm_path"
                fi
            else
                echo "[WARN] Not deleting untracked new path (not in whitelist): $norm_path"
            fi
        done
        rm -f baseline_files.txt post_files.txt
    fi
}

# Register cleanup trap to ensure cleanup happens even if user interrupts or test fails
trap cleanup_generated_files EXIT

# ----------------------------------------------------
# SESSION 1: De Novo Assembly & Gene Prediction
# ----------------------------------------------------
echo -e "\n--- [Running Session 1: Genomics] ---"
if pixi run -e genomika spades.py -1 data/session1/yeast_chrom1_R1.fastq.gz -2 data/session1/yeast_chrom1_R2.fastq.gz -o spades_out/ --only-assembler && \
   pixi run -e genomika augustus --species=saccharomyces --gff3=on spades_out/scaffolds.fasta > genes.gff3; then
    
    # Validation
    if [ -s spades_out/scaffolds.fasta ] && [ -s genes.gff3 ] && grep -q ">" spades_out/scaffolds.fasta && grep -q "augustus" genes.gff3; then
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
if pixi run -e transkriptomika salmon index -t data/original/escherichia_coli_transcriptome.fasta -i salmon_index/ && \
   (
     success=1
     for SAMPLE in CTL_rep1 CTL_rep2 CTL_rep3 TRT_rep1 TRT_rep2 TRT_rep3; do
         pixi run -e transkriptomika salmon quant -i salmon_index/ -l A -1 data/session2/${SAMPLE}_R1.fastq.gz -2 data/session2/${SAMPLE}_R2.fastq.gz -o data/session2/${SAMPLE} || success=0
     done
     [ $success -eq 1 ]
   ); then
    
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
message("Testing DESeq2 execution...")
library(DESeq2)
cnts <- matrix(rnbinom(n=100, mu=100, size=1/0.5), ncol=10)
cond <- factor(rep(c("A","B"), each=5))
coldata <- data.frame(row.names=colnames(cnts), cond)
dds <- DESeqDataSetFromMatrix(countData=cnts, colData=coldata, design=~cond)
dds <- DESeq(dds)
res <- results(dds)
if (nrow(res) != 10) stop("DESeq2 output length incorrect")
message("[OK] DESeq2 executes and libraries loaded successfully.")
'; then
    s2_r_status="SUCCESS"
else
    s2_r_status="FAILED"
    echo "[ERROR] Session 2 R validation failed."
fi

# ----------------------------------------------------
# SESSION 3: Read Mapping & Variant Calling
# ----------------------------------------------------
echo -e "\n--- [Running Session 3: Polymorphisms (BWA & bcftools)] ---"
if pixi run -e polimorfizmi bwa index data/session3/ecoli_ref_region.fasta && \
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
   pixi run -e polimorfizmi bcftools call -mv -o data/session3/variants.vcf; then
    
    # Validation
    if [ -s data/session3/variants.vcf ] && grep -q "#CHROM" data/session3/variants.vcf; then
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
echo -e "\n--- [Running Session 4: Population Structure R Modules Check] ---"
if pixi run -e polimorfizmi Rscript -e '
message("Loading polimorfizmi R libraries...")
libs <- c("vcfR", "adegenet", "hierfstat", "ggplot2")
for (lib in libs) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    stop("Package ", lib, " is missing!")
  }
}
message("[OK] Polimorfizmi R libraries loaded successfully.")
'; then
    s4_r_status="SUCCESS"
else
    s4_r_status="FAILED"
    echo "[ERROR] Session 4 R validation failed."
fi

# ----------------------------------------------------
# SESSION 5: 16S Amplicon & Microbiome Diversity (QIIME 2)
# ----------------------------------------------------
echo -e "\n--- [Running Session 5: QIIME 2 Workflows] ---"
if pixi run -e qiime2 qiime feature-classifier classify-sklearn \
     --i-classifier data/session5/silva_classifier_excerpt.qza \
     --i-reads data/session5/rep-seqs.qza \
     --o-classification taxonomy.qza && \
   pixi run -e qiime2 qiime diversity core-metrics-phylogenetic \
     --i-phylogeny data/session5/rooted-tree.qza \
     --i-table data/session5/table.qza \
     --p-sampling-depth 1000 \
     --m-metadata-file data/session5/sample_metadata.tsv \
     --output-dir core-metrics-results && \
   pixi run -e qiime2 qiime taxa barplot \
     --i-table data/session5/table.qza \
     --i-taxonomy taxonomy.qza \
     --m-metadata-file data/session5/sample_metadata.tsv \
     --o-visualization taxa-bar-plots.qzv; then
    
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
