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
echo "Checking transkriptomika..."
pixi run -e transkriptomika salmon --version 2>&1 | head -n 1
echo "Checking polimorfizmi..."
pixi run -e polimorfizmi samtools --version | head -n 1
echo "Checking diverziteta..."
pixi run -e diverziteta Rscript --version 2>&1 | head -n 1
echo "Checking qiime2..."
pixi run -e qiime2 qiime info | head -n 10

echo -e "\n[SUCCESS] All environments and data are verified successfully!"
