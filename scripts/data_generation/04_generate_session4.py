#!/usr/bin/env python3
"""
Simulate a multi-sample diploid VCF file for 20 individuals in 4 families.

This script generates synthetic variant call genotypes for 20 diploid samples
spanning 4 pedigree families. Point mutations are modeled with realistic allele
frequencies and Mendelian inheritance, and formatting complies with the VCF 4.2 spec.
"""

import sys
import logging
import pathlib
import numpy as np

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
    out_dir = repo_root / "data" / "session4"
    out_dir.mkdir(parents=True, exist_ok=True)

    vcf_path = out_dir / "population_structure.vcf"
    meta_path = out_dir / "sample_metadata.tsv"

    # 1. Define families and metadata
    # 20 samples: 4 families of 5 members each (Father, Mother, 3 Children)
    families = {
        "FamilyA": ["A_Father", "A_Mother", "A_Child1", "A_Child2", "A_Child3"],
        "FamilyB": ["B_Father", "B_Mother", "B_Child1", "B_Child2", "B_Child3"],
        "FamilyC": ["C_Father", "C_Mother", "C_Child1", "C_Child2", "C_Child3"],
        "FamilyD": ["D_Father", "D_Mother", "D_Child1", "D_Child2", "D_Child3"]
    }
    
    samples = []
    for fam_name, members in families.items():
        samples.extend(members)

    with open(meta_path, "w") as f_meta:
        f_meta.write("sample_id\tfamily\trole\tsex\n")
        for fam_name, members in families.items():
            for m in members:
                role = "Father" if "Father" in m else "Mother" if "Mother" in m else "Child"
                sex = "M" if "Father" in m or m[-1] in ["1", "3"] else "F"
                f_meta.write(f"{m}\t{fam_name}\t{role}\t{sex}\n")
    
    logger.info(f"Wrote metadata table for {len(samples)} samples to {meta_path}")

    # 2. Generate 500 polymorphic positions with minimum 5,000 bp gap
    logger.info("Sampling 500 polymorphic positions on chr1...")
    CHROM_LEN = 5_000_000
    MIN_GAP = 5000
    positions = []
    
    # Rejection sampling for well-spaced SNPs
    while len(positions) < 500:
        pos = int(rng.integers(10000, CHROM_LEN - 10000))
        if all(abs(pos - p) >= MIN_GAP for p in positions):
            positions.append(pos)
            
    positions.sort()
    logger.info("Finished sampling variant positions.")

    # 3. Create VCF file
    logger.info(f"Generating VCF records and writing to {vcf_path}...")
    with open(vcf_path, "w") as f_vcf:
        # Write headers
        f_vcf.write("##fileformat=VCFv4.2\n")
        f_vcf.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
        f_vcf.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n")
        f_vcf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f_vcf.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n")
        f_vcf.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n")
        f_vcf.write(f"##contig=<ID=chr1,length={CHROM_LEN}>\n")
        
        # Header line
        header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n"
        f_vcf.write(header_line)

        # Write records
        bases = ["T", "G", "C"]
        for pos in positions:
            ref = "A"
            alt = rng.choice(bases)
            
            # Allele frequency drawn from Beta(2,2)
            p = rng.beta(2.0, 2.0)
            
            # We will calculate genotypes for each sample
            genotypes = {}
            
            # Process each family
            for fam_name, members in families.items():
                father = members[0]
                mother = members[1]
                children = members[2:]
                
                # Father genotypes
                f_a1, f_a2 = rng.choice([0, 1], p=[1-p, p], size=2)
                genotypes[father] = (f_a1, f_a2)
                
                # Mother genotypes
                m_a1, m_a2 = rng.choice([0, 1], p=[1-p, p], size=2)
                genotypes[mother] = (m_a1, m_a2)
                
                # Children genotypes (Mendelian inheritance)
                for child in children:
                    c_a1 = rng.choice([f_a1, f_a2])
                    c_a2 = rng.choice([m_a1, m_a2])
                    genotypes[child] = (c_a1, c_a2)

            # Compute population Allele Frequency (AF) for the current variant
            total_alleles = 2 * len(samples)
            alt_alleles = sum(sum(gt) for gt in genotypes.values())
            computed_af = alt_alleles / total_alleles

            # Format the VCF row fields
            info = f"AF={computed_af:.4f}"
            fmt = "GT:DP:GQ"
            
            sample_fields = []
            for sample in samples:
                gt_tuple = genotypes[sample]
                gt_str = f"{gt_tuple[0]}/{gt_tuple[1]}"
                dp = rng.integers(15, 36)
                gq = rng.integers(40, 100)
                sample_fields.append(f"{gt_str}:{dp}:{gq}")

            row = f"chr1\t{pos}\t.\t{ref}\t{alt}\t60.0\tPASS\t{info}\t{fmt}\t" + "\t".join(sample_fields) + "\n"
            f_vcf.write(row)

    logger.info(f"Successfully wrote {len(positions)} VCF records to {vcf_path}.")

if __name__ == "__main__":
    main()
