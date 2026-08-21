#!/usr/bin/env bash
#
# Loader for the extended integration test.
#
# The test runs all five session pipelines end to end, including the optional
# tasks, so it contains every solution. The repository is public, so the test
# itself is encrypted (`tests/test_extended.sh.enc`) and only the instructor,
# who has the key, can run it.
#
# Usage (instructor):
#   pixi run test-extended <KEY>
#   bash tests/test_extended.sh <KEY>
#
# This test is not meant for students: to verify their installation they run
# `pixi run test`, which is public and needs no key.
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

ENC="tests/test_extended.sh.enc"
KEY="${1:-${VAJE_KEY:-}}"

if [ -z "$KEY" ]; then
    cat >&2 <<'EOF'
==========================================================
  RAZŠIRJENI TEST JE ZAKLENJEN
==========================================================
Ta test je namenjen učitelju in zahteva ključ:

    pixi run test-extended <KLJUC>

Za preverjanje svoje namestitve uporabite javni test, ki
ključa ne potrebuje:

    pixi run test
==========================================================
EOF
    exit 2
fi

if [ ! -f "$ENC" ]; then
    echo "[NAPAKA] Šifrirane datoteke $ENC ni." >&2
    exit 1
fi

if ! command -v openssl >/dev/null 2>&1; then
    echo "[NAPAKA] Program openssl ni na voljo, odklepanje ni mogoče." >&2
    exit 1
fi

TMP="$(mktemp -t vaje_extended.XXXXXXXX)"
chmod 600 "$TMP"
cleanup() { rm -f "$TMP"; }
trap cleanup EXIT INT TERM

if ! openssl enc -aes-256-cbc -pbkdf2 -a -d -in "$ENC" -out "$TMP" -pass pass:"$KEY" 2>/dev/null; then
    echo "[NAPAKA] Napačen ključ — testa ni mogoče odkleniti." >&2
    exit 3
fi

# Make sure the decrypted content really is a script and not random bytes:
# a wrong key can occasionally produce output without reporting an error.
if ! head -n 1 "$TMP" | grep -q '^#!/usr/bin/env bash'; then
    echo "[NAPAKA] Napačen ključ — odšifrirana vsebina ni veljavna skripta." >&2
    exit 3
fi

# The script runs from a temporary file, so pass the repository path explicitly.
export VAJE_REPO_ROOT="$REPO_ROOT"
bash "$TMP"
