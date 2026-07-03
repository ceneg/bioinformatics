#!/bin/bash
set -e

# --- Environment Tests ---
echo -e "\n--- [Running Basic Environment Checks] ---"
pixi run test

# ----------------------------------------------------
# EXTENDED WORKFLOW INTEGRATION TESTS
# ----------------------------------------------------
echo -e "\n=========================================================="
echo "          STARTING EXTENDED INTEGRATION TESTS"
echo "=========================================================="
echo "Note: This will generate real outputs and then clean them up."

# Store base directory
BASE_DIR=$(pwd)

# Trap to ensure cleanup happens even if script fails
cleanup() {
    echo -e "\n=========================================================="
    echo "          CLEANING UP GENERATED WORKFLOW ARTIFACTS"
    echo "=========================================================="
    
    # We will ONLY delete files that are explicitly known to be generated
    # by this script. To do this safely, we use a whitelist approach.
    
    # Run git clean but ONLY in dry-run mode to see what WOULD be deleted
    # compared to the tracked files and the .gitignore.
    
    # For safety, we just explicitly delete the specific output files/folders 
    # we expect our tools to have generated, IF they exist.
    
    # Find all untracked files/folders
    untracked=$(git ls-files --others --exclude-standard)
    
    for item in $untracked; do
        if [ -e "$item" ]; then
            # Normalize path
            norm_path=$(echo "$item" | sed "s|^$BASE_DIR/||")
            
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
               [[ "$norm_path" == "scratch_trans_list.txt" ]]; then
                is_safe=1
            fi
            
            if [ $is_safe -eq 1 ]; then
                if [ -d "$item" ]; then
                    rm -rf "$item"
                    echo "[CLEAN] Removed directory: $norm_path"
                else
                    rm -f "$item"
                    echo "[CLEAN] Removed file: $norm_path"
                fi
            fi
        fi
    done
    
    # Special cleanup for spades_out to ensure it is fully removed
    if [ -d "spades_out" ]; then
        rm -rf "spades_out"
        echo "[CLEAN] Removed directory: spades_out"
    fi
}
trap cleanup EXIT

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
if pixi run -e qiime2 qiime dada2 denoise-single \
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
