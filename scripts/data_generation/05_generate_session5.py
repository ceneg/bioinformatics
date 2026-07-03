#!/usr/bin/env python3
"""
Simulate a 16S rRNA amplicon sequencing dataset and import into QIIME2 using InSilicoSeq.

This script samples ASV representative sequences from the SILVA database, constructs
a sparse abundance matrix for 10 samples simulating distinct environment profiles,
runs InSilicoSeq in amplicon mode (NovaSeq model) to generate realistic reads,
renames files to Casava 1.8 format, writes a QIIME2 manifest, and imports the raw data
into a QIIME2 artifact (.qza).
"""

import sys
import logging
import pathlib
import subprocess
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

def check_qiime2():
    """Verify if QIIME2 is available in the current environment."""
    try:
        res = subprocess.run(["qiime", "--version"], capture_output=True, text=True)
        if res.returncode == 0:
            logger.info("QIIME2 detected successfully in the environment path.")
            return True
    except FileNotFoundError:
        pass
    
    # Print friendly help message and exit with 1 if missing
    print("\n" + "="*80)
    print("QIIME2 ENVIRONMENT NOT DETECTED")
    print("="*80)
    print("Session 5 data generation requires a configured QIIME2 environment.")
    print("Please follow these steps to proceed:\n")
    print("1. Install QIIME2 via Conda as described in the official docs:")
    print("   https://docs.qiime2.org/\n")
    print("2. If you are using the instructor's pre-configured environment:")
    print("   conda activate rachis-qiime2-2026.4\n")
    print("3. Alternatively, run this script inside the specific environment:")
    print("   conda run -n qiime2-2024 python scripts/data_generation/05_generate_session5.py")
    print("="*80 + "\n")
    
    sys.exit(1)

def main():
    # QIIME2 Check (Strict exit guardrail)
    check_qiime2()
    
    repo_root = pathlib.Path(__file__).resolve().parent.parent.parent
    ref_path = repo_root / "data" / "original" / "silva_16S.fasta"
    out_dir = repo_root / "data" / "session5"
    raw_reads_dir = out_dir / "raw_reads"
    raw_reads_dir.mkdir(parents=True, exist_ok=True)
    
    rep_seqs_path = out_dir / "representative_sequences.fasta"
    manifest_path = out_dir / "manifest.tsv"
    qza_path = out_dir / "amplicon_demux.qza"
    temp_asvs_path = out_dir / "temp_asvs.fasta"

    if not ref_path.exists():
        logger.error(f"SILVA 16S reference not found: {ref_path}")
        sys.exit(1)

    # 1. Parse SILVA and sample 50 distinct sequences as ASVs
    logger.info("Loading SILVA 16S reference sequences...")
    records = list(SeqIO.parse(ref_path, "fasta"))
    if len(records) < 50:
        logger.error(f"SILVA file only contains {len(records)} sequences; need at least 50.")
        sys.exit(1)
        
    sampled_indices = rng.choice(len(records), 50, replace=False)
    selected_records = [records[i] for i in sampled_indices]
    
    features = {}
    asv_records = []
    with open(rep_seqs_path, "w") as f_rep:
        for i, rec in enumerate(selected_records):
            asv_id = f"ASV{i+1:03d}"
            # Truncate to 253 bp to represent standard V3-V4 amplicon size
            asv_seq = str(rec.seq[:253]).upper().replace("U", "T")
            features[asv_id] = asv_seq
            f_rep.write(f">{asv_id}\n{asv_seq}\n")
            
            asv_records.append(SeqRecord(
                Seq(asv_seq),
                id=asv_id,
                description=f"16S ASV {asv_id}"
            ) )
            
    logger.info(f"Sampled 50 ASVs and wrote to {rep_seqs_path}")
    
    # Write temporary multifasta for InSilicoSeq
    SeqIO.write(asv_records, temp_asvs_path, "fasta")

    # 2. Define 10 sample profiles and depths
    # High: 3,000 - 4,000 reads; Medium: 1,000 - 1,500 reads; Low: 500 reads
    sample_depths = [
        int(rng.integers(3000, 4001)),  # Sample1 (High)
        int(rng.integers(3000, 4001)),  # Sample2 (High)
        int(rng.integers(3000, 4001)),  # Sample3 (High)
        int(rng.integers(1000, 1501)),  # Sample4 (Medium)
        int(rng.integers(1000, 1501)),  # Sample5 (Medium)
        int(rng.integers(3000, 4001)),  # Sample6 (High)
        int(rng.integers(3000, 4001)),  # Sample7 (High)
        int(rng.integers(3000, 4001)),  # Sample8 (High)
        int(rng.integers(500, 801)),    # Sample9 (Low)
        int(rng.integers(500, 801))     # Sample10 (Low)
    ]

    # 3. Construct sparse abundance matrix (50 features x 10 samples)
    n_samples = 10
    n_features = 50
    matrix = np.zeros((n_features, n_samples))

    # Features 1-15: High in Samples 1-5 only
    for i in range(0, 15):
        matrix[i, 0:5] = rng.lognormal(5.0, 0.5, size=5)
        
    # Features 16-30: High in Samples 6-10 only
    for i in range(15, 30):
        matrix[i, 5:10] = rng.lognormal(5.0, 0.5, size=5)
        
    # Features 31-45: Background present in all samples
    for i in range(30, 45):
        matrix[i, :] = rng.lognormal(2.0, 1.5, size=10)
        
    # Features 46-50: Rare biosphere (present in 1-2 samples at very low levels)
    for i in range(45, 50):
        target_s = rng.choice(10, size=rng.integers(1, 3), replace=False)
        matrix[i, target_s] = rng.integers(1, 5, size=len(target_s))

    # Normalize matrix to sample depths
    for col_idx in range(n_samples):
        depth = sample_depths[col_idx]
        col_vals = matrix[:, col_idx]
        col_sum = np.sum(col_vals)
        if col_sum > 0:
            normalized = np.round((col_vals / col_sum) * depth).astype(int)
            diff = depth - np.sum(normalized)
            max_idx = np.argmax(normalized)
            normalized[max_idx] = max(0, normalized[max_idx] + diff)
            matrix[:, col_idx] = normalized
        else:
            matrix[:, col_idx] = 0

    # 4. Simulate reads for each sample
    logger.info("Simulating reads for all samples...")
    
    with open(manifest_path, "w") as f_man:
        f_man.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
        
        for s_idx in range(n_samples):
            sample_id = f"Sample{s_idx+1}"
            
            # Write temporary readcounts for this sample
            temp_counts_path = out_dir / f"temp_{sample_id}_counts.txt"
            with open(temp_counts_path, "w") as f_cnt:
                for f_idx in range(n_features):
                    asv_id = f"ASV{f_idx+1:03d}"
                    n_reads = int(matrix[f_idx, s_idx])
                    if n_reads > 0:
                        f_cnt.write(f"{asv_id}\t{2 * n_reads}\n")
            
            # Run InSilicoSeq in amplicon mode
            out_prefix = raw_reads_dir / f"{sample_id}_S{s_idx+1}_L001"
            try:
                utils.run_insilicoseq(
                    fasta_path=str(temp_asvs_path),
                    output_prefix=str(out_prefix),
                    readcount_path=str(temp_counts_path),
                    model="NovaSeq",
                    sequence_type="amplicon",
                    seed=42
                )
            finally:
                if temp_counts_path.exists():
                    temp_counts_path.unlink()

            # Rename output files to Casava 1.8 convention: _R1.fastq.gz -> _R1_001.fastq.gz
            orig_r1 = raw_reads_dir / f"{sample_id}_S{s_idx+1}_L001_R1.fastq.gz"
            orig_r2 = raw_reads_dir / f"{sample_id}_S{s_idx+1}_L001_R2.fastq.gz"
            
            casava_r1 = raw_reads_dir / f"{sample_id}_S{s_idx+1}_L001_R1_001.fastq.gz"
            casava_r2 = raw_reads_dir / f"{sample_id}_S{s_idx+1}_L001_R2_001.fastq.gz"
            
            orig_r1.rename(casava_r1)
            orig_r2.rename(casava_r2)
            
            # Write absolute paths to manifest
            f_man.write(f"{sample_id}\t{casava_r1.resolve()}\t{casava_r2.resolve()}\n")
            logger.info(f"Simulated amplicon reads for {sample_id}.")

    # Clean up temporary ASV FASTA
    if temp_asvs_path.exists():
        temp_asvs_path.unlink()

    logger.info(f"Wrote QIIME2 manifest to {manifest_path}")

    # 5. Import into QIIME2
    logger.info("Importing FASTQ files into QIIME2 artifact...")
    try:
        subprocess.run([
            "qiime", "tools", "import",
            "--type", "SampleData[PairedEndSequencesWithQuality]",
            "--input-path", str(manifest_path),
            "--output-path", str(qza_path),
            "--input-format", "PairedEndFastqManifestPhred33V2"
        ], check=True)
        logger.info(f"Successfully created QIIME2 artifact: {qza_path}")
    except subprocess.CalledProcessError as e:
        logger.error(f"QIIME2 import failed with error: {e}")
        sys.exit(1)

    # 6. Generate sample_metadata.tsv
    sample_metadata_path = out_dir / "sample_metadata.tsv"
    logger.info(f"Generating sample metadata at {sample_metadata_path}...")
    with open(sample_metadata_path, "w") as f_meta:
        f_meta.write("sample-id\tgroup\n")
        for s_idx in range(n_samples):
            group = "Group1" if s_idx < 5 else "Group2"
            f_meta.write(f"Sample{s_idx+1}\t{group}\n")

    # 7. Generate table.qza
    table_tsv_path = out_dir / "table.tsv"
    table_biom_path = out_dir / "table.biom"
    table_qza_path = out_dir / "table.qza"
    
    logger.info("Writing feature table TSV...")
    with open(table_tsv_path, "w") as f_tab:
        sample_ids = [f"Sample{i+1}" for i in range(n_samples)]
        f_tab.write("#OTU ID\t" + "\t".join(sample_ids) + "\n")
        for f_idx in range(n_features):
            asv_id = f"ASV{f_idx+1:03d}"
            f_tab.write(f"{asv_id}\t" + "\t".join(str(int(matrix[f_idx, s_idx])) for s_idx in range(n_samples)) + "\n")
            
    logger.info("Converting feature table to BIOM format...")
    try:
        subprocess.run([
            "biom", "convert",
            "-i", str(table_tsv_path),
            "-o", str(table_biom_path),
            "--table-type", "OTU table",
            "--to-hdf5"
        ], check=True)
        
        logger.info("Importing BIOM table into QIIME2...")
        subprocess.run([
            "qiime", "tools", "import",
            "--type", "FeatureTable[Frequency]",
            "--input-path", str(table_biom_path),
            "--output-path", str(table_qza_path)
        ], check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to generate table.qza: {e}")
        sys.exit(1)
    finally:
        if table_tsv_path.exists():
            table_tsv_path.unlink()
        if table_biom_path.exists():
            table_biom_path.unlink()

    # 8. Generate rep-seqs.qza
    rep_seqs_qza_path = out_dir / "rep-seqs.qza"
    logger.info("Importing representative sequences into QIIME2...")
    try:
        subprocess.run([
            "qiime", "tools", "import",
            "--type", "FeatureData[Sequence]",
            "--input-path", str(rep_seqs_path),
            "--output-path", str(rep_seqs_qza_path)
        ], check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to import representative sequences: {e}")
        sys.exit(1)

    # 9. Generate rooted-tree.qza
    rooted_tree_path = out_dir / "rooted-tree.qza"
    aligned_rep_seqs = out_dir / "aligned-rep-seqs.qza"
    masked_aligned_rep_seqs = out_dir / "masked-aligned-rep-seqs.qza"
    unrooted_tree = out_dir / "unrooted-tree.qza"
    
    logger.info("Building phylogenetic tree with MAFFT and FastTree...")
    try:
        subprocess.run([
            "qiime", "phylogeny", "align-to-tree-mafft-fasttree",
            "--i-sequences", str(rep_seqs_qza_path),
            "--o-alignment", str(aligned_rep_seqs),
            "--o-masked-alignment", str(masked_aligned_rep_seqs),
            "--o-tree", str(unrooted_tree),
            "--o-rooted-tree", str(rooted_tree_path)
        ], check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to build phylogenetic tree: {e}")
        sys.exit(1)
    finally:
        for p in [aligned_rep_seqs, masked_aligned_rep_seqs, unrooted_tree]:
            if p.exists():
                p.unlink()

    # 10. Generate silva_classifier_excerpt.qza
    ref_taxonomy_tsv = out_dir / "ref_taxonomy.tsv"
    ref_taxonomy_qza = out_dir / "ref-taxonomy.qza"
    classifier_qza_path = out_dir / "silva_classifier_excerpt.qza"
    
    logger.info("Writing reference taxonomy mapping...")
    with open(ref_taxonomy_tsv, "w") as f_tax:
        for i, rec in enumerate(selected_records):
            asv_id = f"ASV{i+1:03d}"
            # Extract taxon path from SILVA header description
            parts = rec.description.split()
            taxon = parts[1] if len(parts) > 1 else "Unassigned"
            f_tax.write(f"{asv_id}\t{taxon}\n")
            
    logger.info("Importing reference taxonomy into QIIME2...")
    try:
        subprocess.run([
            "qiime", "tools", "import",
            "--type", "FeatureData[Taxonomy]",
            "--input-format", "HeaderlessTSVTaxonomyFormat",
            "--input-path", str(ref_taxonomy_tsv),
            "--output-path", str(ref_taxonomy_qza)
        ], check=True)
        
        logger.info("Training Naive Bayes classifier on excerpt...")
        subprocess.run([
            "qiime", "feature-classifier", "fit-classifier-naive-bayes",
            "--i-reference-reads", str(rep_seqs_qza_path),
            "--i-reference-taxonomy", str(ref_taxonomy_qza),
            "--o-classifier", str(classifier_qza_path)
        ], check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to train QIIME2 classifier: {e}")
        sys.exit(1)
    finally:
        if ref_taxonomy_tsv.exists():
            ref_taxonomy_tsv.unlink()
        if ref_taxonomy_qza.exists():
            ref_taxonomy_qza.unlink()

    logger.info("Session 5 data simulation and QIIME2 artifact generation completed successfully!")

if __name__ == "__main__":
    main()
