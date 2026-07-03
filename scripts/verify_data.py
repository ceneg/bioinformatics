#!/usr/bin/env python3
import pathlib
import sys
import gzip

def main():
    repo_root = pathlib.Path(__file__).resolve().parent.parent
    data_dir = repo_root / "data"
    
    print("=== Verification of Simulated Data ===")

    # 1. Verify FASTQ / FASTQ.GZ files records count (>= 1000 records, which is 4000 lines)
    fastq_paths = list(data_dir.glob("session*/**/*.fastq")) + \
                  list(data_dir.glob("session*/**/*.fastq.gz")) + \
                  list(data_dir.glob("session*/*.fastq")) + \
                  list(data_dir.glob("session*/*.fastq.gz"))
    
    # Remove duplicates from multiple glob matches
    fastq_paths = sorted(list(set(fastq_paths)))
    fastq_ok = True
    
    for path in fastq_paths:
        # Check if the file is gzipped
        if path.name.endswith(".gz"):
            with gzip.open(path, "rt") as f:
                line_count = sum(1 for _ in f)
        else:
            with open(path, "r") as f:
                line_count = sum(1 for _ in f)
        
        record_count = line_count // 4
        size_mb = path.stat().st_size / (1024 * 1024)
        
        # Skip Low depth session5 samples (Sample 9 and 10 which have ~500-800 reads) from >= 1000 check, but make sure they are valid
        if "Sample9" in path.name or "Sample10" in path.name:
            expected_min = 500
        else:
            expected_min = 1000

        if record_count < expected_min:
            print(f"[FAIL] {path.relative_to(repo_root)}: only {record_count} records (Expected >= {expected_min})")
            fastq_ok = False
        else:
            print(f"[OK] {path.relative_to(repo_root)}: {record_count} records ({size_mb:.2f} MB)")

    # 2. Verify VCF
    vcf_path = data_dir / "session4" / "population_structure.vcf"
    if not vcf_path.exists():
        print("[FAIL] VCF file does not exist.")
        sys.exit(1)
        
    with open(vcf_path, "r") as f_vcf:
        vcf_lines = f_vcf.readlines()
        
    header_lines = [l for l in vcf_lines if l.startswith("#")]
    variant_lines = [l for l in vcf_lines if not l.startswith("#")]
    
    print("\n=== VCF Validation ===")
    print(f"Header lines count: {len(header_lines)}")
    print(f"Variant records count: {len(variant_lines)}")
    
    # Check 20 samples columns in the main header line
    main_header = header_lines[-1].strip().split("\t")
    samples = main_header[9:]
    print(f"Samples found in VCF ({len(samples)}): {', '.join(samples)}")
    
    vcf_ok = True
    if len(variant_lines) != 500:
        print(f"[FAIL] VCF variant count is {len(variant_lines)}, expected 500.")
        vcf_ok = False
    if len(samples) != 20:
        print(f"[FAIL] VCF sample count is {len(samples)}, expected 20.")
        vcf_ok = False
        
    if fastq_ok and vcf_ok:
        print("\n[SUCCESS] All datasets successfully validated.")
        sys.exit(0)
    else:
        print("\n[ERROR] Validation failed.")
        sys.exit(1)

if __name__ == "__main__":
    main()
