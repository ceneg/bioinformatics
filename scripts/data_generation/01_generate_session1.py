#!/usr/bin/env python3
"""
Simulate WGS reads for Saccharomyces cerevisiae Chromosome I using InSilicoSeq.

This script parses the complete yeast genome, extracts Chromosome I based on size
and name, and simulates paired-end Illumina-like sequencing reads at 30x coverage
using the pre-trained InSilicoSeq NovaSeq error model. The output reads are compressed
in .fastq.gz format.
"""

import sys
import logging
import pathlib
from Bio import SeqIO
import numpy as np

# Configure path so that it can find the local utils module
sys.path.append(str(pathlib.Path(__file__).resolve().parent))
import utils

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# Output execution guardrails: check python version
if sys.version_info < (3, 11):
    logger.error("Python 3.11 or higher is required.")
    sys.exit(1)

def main():
    repo_root = pathlib.Path(__file__).resolve().parent.parent.parent
    ref_path = repo_root / "data" / "original" / "saccharomyces_cerevisiae.fasta"
    out_dir = repo_root / "data" / "session1"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    out_prefix = out_dir / "yeast_chrom1"
    temp_fasta = out_dir / "yeast_chrom1_temp.fasta"

    # Verify input exists
    if not ref_path.exists():
        logger.error(f"Reference file not found: {ref_path}")
        logger.error("Please run the reference download script first.")
        sys.exit(1)

    # 1. Parse references
    logger.info(f"Parsing yeast genome reference from {ref_path}...")
    records = list(SeqIO.parse(ref_path, "fasta"))
    if not records:
        logger.error(f"No FASTA records found in {ref_path}.")
        sys.exit(1)

    # 2. Select Chromosome I (typically shortest chromosome ~230 kb)
    chrom1_record = None
    for rec in records:
        is_chrom_i = (
            rec.id == "I" 
            or rec.id == "chrI" 
            or "chromosome:I" in rec.description 
            or "chromosome I" in rec.description.lower()
            or (200000 <= len(rec) <= 250000)
        )
        if is_chrom_i:
            chrom1_record = rec
            break

    if chrom1_record is None:
        logger.error("Could not identify Chromosome I in the yeast reference file.")
        sys.exit(1)

    logger.info(f"Selected Chromosome I: ID={chrom1_record.id}, Length={len(chrom1_record)} bp")

    # Write selected chromosome to a temporary file for InSilicoSeq
    SeqIO.write(chrom1_record, temp_fasta, "fasta")
    logger.info(f"Wrote temporary chromosome FASTA to {temp_fasta}")

    # 3. Sequencing constants
    READ_LENGTH = 150
    TARGET_COVERAGE = 30

    # 4. Calculate number of fragments/reads (InSilicoSeq generates n_reads total read pairs)
    n_fragments = (len(chrom1_record) * TARGET_COVERAGE) // (2 * READ_LENGTH)
    logger.info(f"Simulating {n_fragments} read pairs to achieve {TARGET_COVERAGE}x coverage...")

    # 5. Run InSilicoSeq
    try:
        utils.run_insilicoseq(
            fasta_path=str(temp_fasta),
            output_prefix=str(out_prefix),
            n_reads=2 * n_fragments,
            model="NovaSeq",
            seed=42
        )
    finally:
        # Clean up temporary FASTA
        if temp_fasta.exists():
            temp_fasta.unlink()
            logger.info(f"Cleaned up temporary FASTA: {temp_fasta}")

    # Print summary of generated files
    out_r1 = out_dir / "yeast_chrom1_R1.fastq.gz"
    out_r2 = out_dir / "yeast_chrom1_R2.fastq.gz"
    
    if out_r1.exists() and out_r2.exists():
        r1_size = out_r1.stat().st_size / (1024 * 1024)
        r2_size = out_r2.stat().st_size / (1024 * 1024)
        logger.info("Session 1 WGS data generation completed.")
        logger.info(f"R1 file: {out_r1.name} ({r1_size:.2f} MB)")
        logger.info(f"R2 file: {out_r2.name} ({r2_size:.2f} MB)")
    else:
        logger.error("Expected output FASTQ.GZ files not found.")
        sys.exit(1)

if __name__ == "__main__":
    main()
