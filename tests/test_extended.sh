#!/usr/bin/env bash
#
# Zaganjalnik razširjenega integracijskega testa.
#
# Ta test požene celotne poteke vseh petih vaj, vključno z dodatnimi nalogami,
# zato vsebuje tudi vse rešitve. Ker je repozitorij javen, je sama koda testa
# šifrirana (`tests/test_extended.sh.enc`) in jo lahko požene le učitelj s
# ključem.
#
# Uporaba (za učitelja):
#   pixi run test-extended <KLJUC>
#   bash tests/test_extended.sh <KLJUC>
#
# Študentom ta test ni namenjen — za preverjanje namestitve uporabite
# `pixi run test`, ki je javen in ne potrebuje ključa.
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

# Preverimo, da je odšifrirana vsebina res skripta in ne naključni bajti
# (napačen ključ lahko v redkih primerih vrne izhod brez napake).
if ! head -n 1 "$TMP" | grep -q '^#!/usr/bin/env bash'; then
    echo "[NAPAKA] Napačen ključ — odšifrirana vsebina ni veljavna skripta." >&2
    exit 3
fi

# Skripta teče iz začasne datoteke, zato ji pot do repozitorija podamo posebej.
export VAJE_REPO_ROOT="$REPO_ROOT"
bash "$TMP"
