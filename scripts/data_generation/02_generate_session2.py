#!/usr/bin/env python3
"""
Simulate RNA-Seq reads for Escherichia coli comparing Control vs Treated using InSilicoSeq.

This script loads the E. coli transcriptome reference, assigns baseline expression
levels based on a log-normal distribution, introduces differential expression (DE)
for a subset of transcripts, simulates biological replication with Poisson noise,
and uses InSilicoSeq (NovaSeq model) to generate realistic gzipped FASTQ reads for
3 control replicates and 3 treated replicates.
"""

import sys
import logging
import pathlib
import math
from Bio import SeqIO
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
    ref_path = repo_root / "data" / "original" / "escherichia_coli_transcriptome.fasta"
    out_dir = repo_root / "data" / "session2"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    ground_truth_path = out_dir / "ground_truth_de_genes.tsv"

    if not ref_path.exists():
        logger.error(f"Reference transcriptome not found: {ref_path}")
        sys.exit(1)

    # 1. Load reference transcriptome
    logger.info(f"Parsing E. coli transcriptome reference from {ref_path}...")
    all_transcripts = list(SeqIO.parse(ref_path, "fasta"))
    logger.info(f"Total transcripts parsed: {len(all_transcripts)}")

    # 2. Filter transcripts with length >= 150 bp
    transcripts = [tx for tx in all_transcripts if len(tx.seq) >= 150]
    n_transcripts = len(transcripts)
    logger.info(f"Filtered transcripts (length >= 150bp): {n_transcripts}")

    # 3. Assign baseline expression
    baselines = rng.lognormal(mean=3.0, sigma=1.0, size=n_transcripts)
    baselines = (baselines / np.sum(baselines)) * 1_000_000

    # 4. Select DE transcripts (5% up, 5% down)
    de_count = math.floor(n_transcripts * 0.05)
    indices = rng.permutation(n_transcripts)
    up_indices = indices[:de_count]
    down_indices = indices[de_count : 2 * de_count]

    log2fcs = np.zeros(n_transcripts)
    control_abundances = baselines.copy()
    treated_abundances = baselines.copy()

    # Log ground truth
    with open(ground_truth_path, "w") as f_gt:
        f_gt.write("transcript_id\tdirection\tlog2fc\n")
        
        # Upregulated
        for idx in up_indices:
            tx_id = transcripts[idx].id
            l2fc = rng.uniform(1.5, 3.0)
            log2fcs[idx] = l2fc
            treated_abundances[idx] = control_abundances[idx] * (2 ** l2fc)
            f_gt.write(f"{tx_id}\tup\t{l2fc:.4f}\n")
            
        # Downregulated
        for idx in down_indices:
            tx_id = transcripts[idx].id
            l2fc = rng.uniform(-3.0, -1.5)
            log2fcs[idx] = l2fc
            treated_abundances[idx] = control_abundances[idx] * (2 ** l2fc)
            f_gt.write(f"{tx_id}\tdown\t{l2fc:.4f}\n")

    logger.info(f"Assigned {de_count} upregulated and {de_count} downregulated genes. Logged to {ground_truth_path}.")

    # Normalize abundances again to ensure they scale to CPM properly
    control_abundances = (control_abundances / np.sum(control_abundances)) * 1_000_000
    treated_abundances = (treated_abundances / np.sum(treated_abundances)) * 1_000_000

    # 5. Define samples and simulation parameters
    samples = [
        ("CTL_rep1", control_abundances),
        ("CTL_rep2", control_abundances),
        ("CTL_rep3", control_abundances),
        ("TRT_rep1", treated_abundances),
        ("TRT_rep2", treated_abundances),
        ("TRT_rep3", treated_abundances),
    ]

    TOTAL_READ_PAIRS_PER_SAMPLE = 30000

    # 6. Simulate sequencing reads per sample using InSilicoSeq
    for sample_name, abundance_profile in samples:
        logger.info(f"Simulating RNA-Seq reads for sample: {sample_name}...")
        
        # Apply Poisson noise to simulate biological replication variance
        sample_profile = rng.poisson(lam=abundance_profile)
        profile_sum = np.sum(sample_profile)
        if profile_sum == 0:
            sample_profile = abundance_profile
            profile_sum = np.sum(sample_profile)

        # Allocate reads proportionally
        reads_to_simulate = (sample_profile / profile_sum) * TOTAL_READ_PAIRS_PER_SAMPLE
        reads_to_simulate = np.round(reads_to_simulate).astype(int)

        # Write counts to temporary file for InSilicoSeq
        temp_counts_path = out_dir / f"{sample_name}_temp_counts.txt"
        with open(temp_counts_path, "w") as f_cnt:
            for idx, tx in enumerate(transcripts):
                n_reads = reads_to_simulate[idx]
                if n_reads > 0:
                    f_cnt.write(f"{tx.id}\t{2 * n_reads}\n")

        # Run InSilicoSeq with NovaSeq model
        try:
            utils.run_insilicoseq(
                fasta_path=str(ref_path),
                output_prefix=str(out_dir / sample_name),
                readcount_path=str(temp_counts_path),
                model="NovaSeq",
                seed=42
            )
        finally:
            # Clean up temp file
            if temp_counts_path.exists():
                temp_counts_path.unlink()

        out_r1 = out_dir / f"{sample_name}_R1.fastq.gz"
        out_r2 = out_dir / f"{sample_name}_R2.fastq.gz"
        logger.info(f"Finished {sample_name}: generated gzipped read files: {out_r1.name}, {out_r2.name}")

if __name__ == "__main__":
    main()
