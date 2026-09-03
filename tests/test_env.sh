#!/usr/bin/env bash
set -eu

# Ensure we're in the repository root directory
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"


echo "=== Testing Data ==="
python scripts/verify_data.py

# ----------------------------------------------------------------------
# INSTALL EVERYTHING
# ----------------------------------------------------------------------
# This is the only script students run, and it must leave the machine ready
# for all five sessions on its own. Nothing later in the course may install
# anything — in particular not the extended test, which needs a key students
# do not have. If installation moved there, every lab machine would need the
# instructor present to finish its setup.
#
# The post-link setting has to be configured BEFORE the first environment is
# created. Bioconductor data packages (genomeinfodbdata) ship their metadata in
# a post-link script, which pixi skips by default; an environment created
# before this line would be built without it.
echo "=== Preparing Pixi ==="
pixi config set --local run-post-link-scripts insecure

# Install every environment from the frozen lockfile, so students get the
# instructor's exact pixi.lock with no local resolution drift. This is the slow
# step on a fresh machine (roughly 8 GB across five environments).
echo "=== Installing all environments (first run takes a while) ==="
pixi install --frozen --all

# Belt and braces: run any post-link script pixi still skipped.
echo "Checking and executing any pending R post-link scripts..."
shopt -s dotglob
for env in transkriptomika polimorfizmi; do
    if [ -d ".pixi/envs/$env/bin" ]; then
        for script in .pixi/envs/$env/bin/*post-link.sh; do
            if [ -f "$script" ]; then
                echo "Running post-link script: $script"
                # Export PREFIX as CONDA_PREFIX for compatibility with Bioconductor scripts
                if ! pixi run -e "$env" bash -c 'export PREFIX=$CONDA_PREFIX; bash "$1"' _ "$script"; then
                    echo "[WARNING] Post-link script failed: $script"
                    echo "          Bioconductor data packages install their metadata from"
                    echo "          this script, so the R library check below will fail."
                    echo "          The script calls 'yq', which needs 'jq' on PATH. 'jq' is"
                    echo "          declared in pixi.toml for exactly this reason -- if it is"
                    echo "          missing here, the environment was built from a stale lock."
                fi
            fi
        done
    fi
done
shopt -u dotglob

echo "=== Testing Environments ==="
# We disable pipefail temporarily because piping tool outputs to 'head'
# triggers a SIGPIPE in the tool when 'head' closes the pipe, which
# with pipefail would cause the script to fail.
set +o pipefail 2>/dev/null || true

echo "Checking genomika..."
pixi run -e genomika spades.py --version | head -n 1
pixi run -e genomika fastp --version 2>&1 | head -n 1
pixi run -e genomika quast --version | head -n 1
pixi run -e genomika augustus --version 2>&1 | head -n 1
pixi run -e genomika python -c 'import Bio; print("biopython", Bio.__version__)'
echo "Checking transkriptomika..."
pixi run -e transkriptomika salmon --version 2>&1 | head -n 1
pixi run -e transkriptomika Rscript -e 'library(DESeq2); library(tximport)' 2>/dev/null && echo "R packages (DESeq2, tximport) OK" || { echo "[ERROR] R packages missing or broken in transkriptomika!"; exit 1; }
echo "Checking polimorfizmi..."
pixi run -e polimorfizmi bwa 2>&1 | grep -i "^Version" | head -n 1
pixi run -e polimorfizmi samtools --version | head -n 1
pixi run -e polimorfizmi bcftools --version | head -n 1
pixi run -e polimorfizmi vcftools --version | head -n 1
pixi run -e polimorfizmi Rscript -e 'library(vcfR); library(adegenet); library(poppr); library(hierfstat); library(ape)' 2>/dev/null && echo "R packages (vcfR, adegenet, poppr, hierfstat, ape) OK" || { echo "[ERROR] R packages missing or broken in polimorfizmi!"; exit 1; }
echo "Checking qiime2..."
pixi run -e qiime2 qiime info | head -n 10

# ----------------------------------------------------------------------
# Confirm the machine is fully provisioned
# ----------------------------------------------------------------------
# A student must never have to run anything else, so state plainly whether
# every session environment is really on disk. `conda-meta` only exists once
# an environment has been installed completely.
set -o pipefail 2>/dev/null || true
missing=0
for env in genomika transkriptomika polimorfizmi qiime2; do
    if [ -d ".pixi/envs/$env/conda-meta" ]; then
        echo "[OK] Environment '$env' is installed."
    else
        echo "[ERROR] Environment '$env' is NOT installed." >&2
        missing=$((missing + 1))
    fi
done
if [ "$missing" -gt 0 ]; then
    echo "[FAILED] $missing environment(s) missing — rerun this script." >&2
    exit 1
fi
set +o pipefail 2>/dev/null || true

# ----------------------------------------------------------------------
# Protect the input data from accidental modification
# ----------------------------------------------------------------------
# Every file under data/session*/ is a read-only course input. Making them
# non-writable stops an accidental `>` redirect or in-place edit from
# destroying them. The DIRECTORIES stay writable, so the exercises can still
# write their own outputs (e.g. *.sorted.bam, CTL_rep1/) alongside them.
echo "Protecting input data files from accidental overwriting..."
protected=0
while IFS= read -r f; do
    if [ -f "$f" ]; then
        chmod a-w "$f" 2>/dev/null && protected=$((protected + 1))
    fi
done < <(git ls-files -- data 2>/dev/null)
if [ "$protected" -gt 0 ]; then
    echo "[OK] Write-protected $protected input files (directories remain writable)."
    echo "     If you ever need to modify one deliberately: chmod u+w <file>"
else
    echo "[INFO] Could not list tracked data files (not a git checkout); skipping."
fi

echo -e "\n[SUCCESS] All environments and data are verified successfully!"
echo "This computer is ready for all five sessions. Nothing else needs to be installed."
