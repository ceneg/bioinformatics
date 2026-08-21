#!/usr/bin/env bash
#
# Seal the extended test.
#
# The extended test contains the solution to every task, including the ones
# students have to write themselves. The repository is public, so the plaintext
# source is never published: only the encrypted `test_extended.sh.enc` is.
#
#   tests/test_extended_src.sh   <- plaintext source (NOT in git; edit this one)
#   tests/test_extended.sh.enc   <- encrypted copy (in git)
#   tests/test_extended.sh       <- small loader that asks for the key
#
# Run this script after every change to the source, otherwise the encrypted
# copy goes stale.
#
# Usage:
#   bash tests/seal_extended.sh <KEY>
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

SRC="tests/test_extended_src.sh"
ENC="tests/test_extended.sh.enc"

KEY="${1:-${VAJE_KEY:-}}"
if [ -z "$KEY" ]; then
    echo "Uporaba: bash tests/seal_extended.sh <KLJUC>" >&2
    exit 2
fi

if [ ! -f "$SRC" ]; then
    echo "[NAPAKA] Izvorne datoteke $SRC ni." >&2
    echo "         Odšifrirajte jo iz $ENC ali jo obnovite iz varnostne kopije." >&2
    exit 1
fi

# Check the source is syntactically valid before sealing it: a broken script
# cannot easily be inspected once encrypted.
if ! bash -n "$SRC"; then
    echo "[NAPAKA] $SRC vsebuje sintaktično napako. Zaklepanje prekinjeno." >&2
    exit 1
fi

openssl enc -aes-256-cbc -pbkdf2 -a -salt -in "$SRC" -out "$ENC" -pass pass:"$KEY"

# Verify the result really decrypts back to the original.
TMP="$(mktemp)"
trap 'rm -f "$TMP"' EXIT
if ! openssl enc -aes-256-cbc -pbkdf2 -a -d -in "$ENC" -out "$TMP" -pass pass:"$KEY" 2>/dev/null; then
    echo "[NAPAKA] Šifrirane datoteke ni mogoče odšifrirati nazaj." >&2
    exit 1
fi
if ! cmp -s "$SRC" "$TMP"; then
    echo "[NAPAKA] Odšifrirana vsebina se ne ujema z izvirnikom." >&2
    exit 1
fi

echo "[OK] $SRC -> $ENC ($(wc -c < "$ENC") bajtov), preverjeno z odšifriranjem."
echo "     Ne pozabite v git dodati $ENC (izvorne datoteke NE)."
