#!/usr/bin/env bash
set -euo pipefail

# Create output directory
OUT_DIR="data/original"
mkdir -p "$OUT_DIR"

echo "Downloading biological reference databases into $OUT_DIR..."

# Checksum mapping
YEAST_MD5="f4b28ca7cd1554c0bade6ca3ddb1aa8e"
ECOLI_MD5="182c4c48740aa62f1b6ca3c03afd663a"
ECOLI_TRANS_MD5="890645a43a4644a06cce026666934c91"
SILVA_MD5="1b10dc6f102a36fd203c0a3a9fdc67b8"

# Check if md5sum is available
if ! command -v md5sum >/dev/null 2>&1; then
    echo "Warning: md5sum command not found. Skipping checksum verification."
    HAS_MD5SUM=false
else
    HAS_MD5SUM=true
fi

check_md5() {
    local file=$1
    local expected_md5=$2
    if [ "$HAS_MD5SUM" = false ]; then
        [ -f "$file" ] && return 0 || return 1
    fi
    if [ ! -f "$file" ]; then
        return 1
    fi
    local actual_md5
    actual_md5=$(md5sum "$file" | cut -d' ' -f1)
    if [ "$actual_md5" = "$expected_md5" ]; then
        return 0
    else
        return 1
    fi
}

download_and_verify() {
    local label=$1
    local url=$2
    local dest=$3
    local expected_md5=$4
    local temp_dest="${dest}.tmp"
    
    if check_md5 "$dest" "$expected_md5"; then
        echo "[$label] Reference already exists and is verified: $dest"
        return 0
    fi
    
    echo "[$label] Downloading from $url..."
    
    # Trap to clean up temporary file on error/interrupt
    trap 'rm -f "$temp_dest"; echo "Download interrupted. Cleaned up temporary files."; exit 1' INT TERM EXIT
    
    if [[ "$dest" == *saccharomyces_cerevisiae.fasta ]]; then
        curl -sSL "$url" | gunzip -c > "$temp_dest"
    else
        curl -sSL -o "$temp_dest" "$url"
    fi
    
    # Remove trap after successful download
    trap - INT TERM EXIT
    
    if check_md5 "$temp_dest" "$expected_md5"; then
        mv "$temp_dest" "$dest"
        echo "[$label] Finished downloading and verified $dest."
        return 0
    else
        echo "[$label] ERROR: Checksum validation failed for downloaded file."
        rm -f "$temp_dest"
        return 1
    fi
}

# 1. Saccharomyces cerevisiae (Ensembl R64-1-1 dna.toplevel)
YEAST_URL="https://ftp.ensembl.org/pub/release-111/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"
download_and_verify "Yeast" "$YEAST_URL" "$OUT_DIR/saccharomyces_cerevisiae.fasta" "$YEAST_MD5"

# 2. Escherichia coli K-12 MG1655 complete genome (NC_000913.3)
ECOLI_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=fasta&retmode=text"
download_and_verify "E. coli" "$ECOLI_URL" "$OUT_DIR/escherichia_coli.fasta" "$ECOLI_MD5"

# 3. E. coli K-12 MG1655 Transcriptome (NC_000913.3 CDS sequences)
ECOLI_TRANS_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=fasta_cds_na&retmode=text"
download_and_verify "E. coli Transcriptome" "$ECOLI_TRANS_URL" "$OUT_DIR/escherichia_coli_transcriptome.fasta" "$ECOLI_TRANS_MD5"

# 4. SILVA 16S subset (SortMeRNA bac-16s-id85)
SILVA_URL="https://raw.githubusercontent.com/sortmerna/sortmerna/master/data/silva-bac-16s-database-id85.fasta"
download_and_verify "SILVA" "$SILVA_URL" "$OUT_DIR/silva_16S.fasta" "$SILVA_MD5"

echo "All reference files downloaded and verified successfully."
