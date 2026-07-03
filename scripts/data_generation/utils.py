"""
Shared utility functions for simulating genomic and transcriptomic sequencing datasets.

This module provides wrapper functions around InSilicoSeq for simulating realistic
Illumina reads and biological utilities like SNP introduction and reverse complementation.
"""

import sys
import logging
import pathlib
import subprocess
from typing import List, Dict, Tuple
import numpy as np

# Configure logging to ensure clear terminal feedback for instructors.
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# Base complement map for reverse complementation.
COMPLEMENT_MAP = str.maketrans("ATGCatgc", "TACGtacg")

def reverse_complement(sequence: str) -> str:
    """
    Generate the reverse complement of a DNA sequence.

    Reverses the strand direction (3'-to-5' to 5'-to-3') and swaps
    bases with their Watson-Crick partners (A<->T, G<->C).

    Parameters
    ----------
    sequence : str
        The input DNA sequence.

    Returns
    -------
    rev_comp_seq : str
        The reverse complement sequence.
    """
    return sequence.translate(COMPLEMENT_MAP)[::-1]


def introduce_snps(sequence: str, snp_density: float, rng: np.random.Generator) -> Tuple[str, List[Dict]]:
    """
    Introduce random single nucleotide polymorphisms (SNPs) into a DNA sequence.

    Simulates the accumulation of point mutations that distinguish a strain or
    individual from the reference genome. Each base is mutated with a probability
    defined by the snp_density.

    Parameters
    ----------
    sequence : str
        Reference DNA sequence (uppercase).
    snp_density : float
        Proportion of bases to mutate (e.g., 1/1500).
    rng : np.random.Generator
        A seeded NumPy random generator instance.

    Returns
    -------
    mutated_seq : str
        The sequence with SNPs introduced.
    snp_log : list of dict
        A list recording the details (position, ref_base, alt_base) of each SNP.
    """
    seq_list = list(sequence.upper())
    bases = ["A", "T", "G", "C"]
    snp_log = []

    for idx, base in enumerate(seq_list):
        if base in bases and rng.random() < snp_density:
            alt_bases = [b for b in bases if b != base]
            alt_base = rng.choice(alt_bases)
            seq_list[idx] = alt_base
            snp_log.append({
                "position": idx + 1,  # 1-indexed for biological/VCF convention
                "ref_base": base,
                "alt_base": alt_base
            })

    return "".join(seq_list), snp_log


def run_insilicoseq(
    fasta_path: str,
    output_prefix: str,
    n_reads: int = None,
    readcount_path: str = None,
    model: str = "NovaSeq",
    sequence_type: str = "metagenomics",
    seed: int = 42,
    mode: str = "kde"
) -> None:
    """
    Run InSilicoSeq simulation via subprocess.

    Invokes 'iss generate' to simulate reads using a pre-trained error model
    or basic/perfect modes, outputting gzipped FASTQ files.

    Parameters
    ----------
    fasta_path : str
        Path to the reference FASTA file.
    output_prefix : str
        Output file path prefix.
    n_reads : int, optional
        Number of read pairs to generate.
    readcount_path : str, optional
        Path to the abundance/readcount file.
    model : str, default 'NovaSeq'
        Pre-trained error model (HiSeq, NovaSeq, MiSeq, etc.).
    sequence_type : str, default 'metagenomics'
        Sequencing mode ('metagenomics' or 'amplicon').
    seed : int, default 42
        Seed for the random number generator.
    mode : str, default 'kde'
        Simulation mode ('kde', 'basic', or 'perfect').
    """
    repo_root = pathlib.Path(__file__).resolve().parent.parent.parent
    iss_bin = repo_root / ".venv" / "bin" / "iss"

    cmd = [
        str(iss_bin), "generate",
        "--genomes", fasta_path,
        "--output", output_prefix,
        "--compress",
        "--seed", str(seed)
    ]
    
    if mode != "kde":
        cmd += ["--mode", mode]
    else:
        cmd += ["--model", model]

    if sequence_type == "amplicon":
        cmd += ["--sequence_type", "amplicon"]

    if readcount_path:
        cmd += ["--readcount_file", readcount_path]
    elif n_reads:
        cmd += ["--n_reads", str(n_reads)]

    logger.info(f"Running InSilicoSeq: {' '.join(cmd)}")
    
    # Run the command and capture output for debugging
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        logger.error(f"InSilicoSeq simulation failed:\nStdout: {res.stdout}\nStderr: {res.stderr}")
        raise RuntimeError(f"InSilicoSeq failed with exit code {res.returncode}")
        
    logger.info(f"InSilicoSeq completed successfully. Outputs written with prefix: {output_prefix}")

    # Clean up temporary .vcf files created by InSilicoSeq
    for temp_suffix in [".iss.tmp.0.vcf", ".iss.tmp.1.vcf"]:
        vcf_path = pathlib.Path(f"{output_prefix}{temp_suffix}")
        if vcf_path.exists():
            try:
                vcf_path.unlink()
                logger.info(f"Cleaned up temporary InSilicoSeq file: {vcf_path.name}")
            except Exception as e:
                logger.warning(f"Could not delete temporary VCF file {vcf_path}: {e}")

