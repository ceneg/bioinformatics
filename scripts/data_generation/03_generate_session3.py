#!/usr/bin/env python3
"""
Simulate WGS reads for 5 clonal E. coli strains with phylogenetic structure using InSilicoSeq.

This script extracts a 500kb reference region from the E. coli genome, introduces
SNPs along a defined mutational tree to generate 5 distinct genomes, simulates
paired-end reads at 20x coverage for each using InSilicoSeq (NovaSeq model), and
outputs compressed .fastq.gz files.
"""

import sys
import logging
import pathlib
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np

# Configure path so that it can find the local utils module
sys.path.append(str(pathlib.Path(__file__).resolve().parent))
import utils

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# Output execution guardrails: check python version and seed global rng
if sys.version_info < (3, 11):
    logger.error("Python 3.11 or higher is required.")
    sys.exit(1)

GLOBAL_SEED = 42
rng = np.random.default_rng(GLOBAL_SEED)

def main():
    repo_root = pathlib.Path(__file__).resolve().parent.parent.parent
    ref_path = repo_root / "data" / "original" / "escherichia_coli.fasta"
    out_dir = repo_root / "data" / "session3"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    ref_region_path = out_dir / "ecoli_ref_region.fasta"

    if not ref_path.exists():
        logger.error(f"E. coli reference genome not found: {ref_path}")
        sys.exit(1)

    # 1. Parse complete genome and extract first 500,000 bp
    logger.info(f"Parsing E. coli reference genome from {ref_path}...")
    records = list(SeqIO.parse(ref_path, "fasta"))
    if not records:
        logger.error(f"No FASTA records found in {ref_path}.")
        sys.exit(1)
        
    ref_record = records[0]
    region_len = min(500000, len(ref_record))
    ref_region_seq = str(ref_record.seq[:region_len]).upper()
    logger.info(f"Extracted {region_len} bp from record: {ref_record.id}")

    # Write reference region to FASTA (for mapping reference)
    ref_region_record = SeqRecord(
        Seq(ref_region_seq),
        id=ref_record.id,
        description="E. coli reference 500kb region"
    )
    SeqIO.write(ref_region_record, ref_region_path, "fasta")
    logger.info(f"Wrote reference region to {ref_region_path}")

    # 2. Mutational tree variables
    SNP_DENSITY = 1 / 1500
    logger.info(f"Mutating genomes along tree branches (base SNP density: {SNP_DENSITY:.6f})...")

    # Build the tree branches
    # ind1: Ancestral branch from Reference
    ind1_seq, ind1_snps = utils.introduce_snps(ref_region_seq, SNP_DENSITY, rng)
    
    # ind2: Lineage A, branches off ind1
    ind2_seq, ind2_snps = utils.introduce_snps(ind1_seq, SNP_DENSITY, rng)
    
    # ind3: Lineage A.1, branches off ind2
    ind3_seq, ind3_snps = utils.introduce_snps(ind2_seq, SNP_DENSITY, rng)
    
    # ind4: Lineage B, branches off ind1 (independent from A/ind2/ind3)
    ind4_seq, ind4_snps = utils.introduce_snps(ind1_seq, SNP_DENSITY, rng)
    
    # ind5: Outgroup, highly divergent directly from Reference
    ind5_seq, ind5_snps = utils.introduce_snps(ref_region_seq, SNP_DENSITY * 3.0, rng)

    individuals = [
        ("ind1", ind1_seq, len(ind1_snps)),
        ("ind2", ind2_seq, len(ind2_snps)),
        ("ind3", ind3_seq, len(ind3_snps)),
        ("ind4", ind4_seq, len(ind4_snps)),
        ("ind5", ind5_seq, len(ind5_snps))
    ]

    # 3. Simulate and write paired-end reads for each individual using InSilicoSeq
    COVERAGE = 20
    READ_LENGTH = 150
    n_reads = (region_len * COVERAGE) // (2 * READ_LENGTH) # 33333 read pairs

    for name, seq, snp_count in individuals:
        logger.info(f"Simulating reads for {name} ({snp_count} branch SNPs)...")
        
        # Write clone sequence to temporary FASTA for InSilicoSeq
        temp_fasta = out_dir / f"temp_{name}.fasta"
        clone_record = SeqRecord(
            Seq(seq),
            id=name,
            description=f"E. coli clone {name}"
        )
        SeqIO.write(clone_record, temp_fasta, "fasta")
        
        try:
            utils.run_insilicoseq(
                fasta_path=str(temp_fasta),
                output_prefix=str(out_dir / name),
                n_reads=2 * n_reads,
                model="NovaSeq",
                seed=42
            )
        finally:
            # Clean up temporary FASTA
            if temp_fasta.exists():
                temp_fasta.unlink()

        out_r1 = out_dir / f"{name}_R1.fastq.gz"
        out_r2 = out_dir / f"{name}_R2.fastq.gz"
        logger.info(f"Finished {name}: generated gzipped read files: {out_r1.name}, {out_r2.name}")

    logger.info("Session 3 WGS data generation completed.")

if __name__ == "__main__":
    main()
