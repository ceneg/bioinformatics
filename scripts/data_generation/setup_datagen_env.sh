#!/usr/bin/env bash
set -euo pipefail

# Find a suitable python interpreter (prefer python3.11)
PYTHON_BIN=""
for cmd in python3.11 python3 python; do
    if command -v "$cmd" >/dev/null 2>&1; then
        VERSION_MAJOR_MINOR=$("$cmd" -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')
        if [ "$VERSION_MAJOR_MINOR" = "3.11" ] || [ "$(echo -e "3.11\n$VERSION_MAJOR_MINOR" | sort -V | head -n1)" = "3.11" ]; then
            PYTHON_BIN="$cmd"
            break
        fi
    fi
done

if [ -z "$PYTHON_BIN" ]; then
    echo "Error: Python 3.11 or higher is required." >&2
    exit 1
fi

echo "Using Python binary: $PYTHON_BIN"

# Create virtual environment if it doesn't exist
if [ ! -d ".venv" ]; then
    "$PYTHON_BIN" -m venv .venv
    echo "Created virtual environment in .venv/"
else
    echo "Virtual environment .venv/ already exists."
fi

# Activate and install dependencies
source .venv/bin/activate
pip install --upgrade pip
pip install biopython pandas numpy insilicoseq
echo "Dependencies installed successfully in .venv/"

