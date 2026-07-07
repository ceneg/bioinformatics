#!/usr/bin/env bash
set -eu

# Ensure we're in the repository root directory
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

echo "=== Testing Data ==="
python scripts/verify_data.py

echo "=== Testing Environments ==="
# We disable pipefail temporarily because piping tool outputs to 'head'
# triggers a SIGPIPE in the tool when 'head' closes the pipe, which
# with pipefail would cause the script to fail.
set +o pipefail 2>/dev/null || true

echo "Checking genomika..."
pixi run -e genomika spades.py --version | head -n 1
pixi run -e genomika fastp --version 2>&1 | head -n 1
pixi run -e genomika quast --version | head -n 1

# Manually execute post-link scripts for R packages if Pixi skipped them
echo "Checking and executing any pending R post-link scripts..."
shopt -s dotglob
for env in transkriptomika polimorfizmi diverziteta; do
    if [ -d ".pixi/envs/$env/bin" ]; then
        for script in .pixi/envs/$env/bin/*post-link.sh; do
            if [ -f "$script" ]; then
                echo "Running post-link script: $script"
                # Export PREFIX as CONDA_PREFIX for compatibility with Bioconductor scripts
                pixi run -e "$env" bash -c 'export PREFIX=$CONDA_PREFIX; bash "$1"' _ "$script" || echo "[WARNING] Post-link script failed. Continuing..."
            fi
        done
    fi
done
shopt -u dotglob

echo "Checking transkriptomika..."
pixi run -e transkriptomika salmon --version 2>&1 | head -n 1
pixi run -e transkriptomika Rscript -e 'library(DESeq2); library(tximport)' 2>/dev/null && echo "R packages (DESeq2, tximport) OK" || { echo "[ERROR] R packages missing or broken in transkriptomika!"; exit 1; }
echo "Checking polimorfizmi..."
pixi run -e polimorfizmi samtools --version | head -n 1
pixi run -e polimorfizmi Rscript -e 'library(vcfR); library(adegenet)' 2>/dev/null && echo "R packages (vcfR, adegenet) OK" || { echo "[ERROR] R packages missing or broken in polimorfizmi!"; exit 1; }
echo "Checking diverziteta..."
pixi run -e diverziteta Rscript --version 2>&1 | head -n 1
echo "Checking qiime2..."
pixi run -e qiime2 qiime info | head -n 10

echo -e "\n[SUCCESS] All environments and data are verified successfully!"
