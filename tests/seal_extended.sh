#!/usr/bin/env bash
#
# Zaklep razširjenega testa.
#
# Razširjeni test vsebuje rešitve vseh nalog — tudi tistih, ki jih morajo
# študenti napisati sami. Repozitorij je javen, zato izvorne skripte ne
# objavljamo: shranjena je le šifrirana različica `test_extended.sh.enc`.
#
#   tests/test_extended_src.sh   <- izvorna koda (NI v gitu, jo urejate vi)
#   tests/test_extended.sh.enc   <- šifrirana različica (je v gitu)
#   tests/test_extended.sh       <- majhen zaganjalnik, ki zahteva ključ
#
# Po vsaki spremembi izvorne skripte morate zagnati to skripto, sicer bo
# šifrirana različica zastarela.
#
# Uporaba:
#   bash tests/seal_extended.sh <KLJUC>
#
# Ključ je zapisan v README.private.md (ta datoteka ni na GitHubu).
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

SRC="tests/test_extended_src.sh"
ENC="tests/test_extended.sh.enc"

KEY="${1:-${VAJE_KEY:-}}"
if [ -z "$KEY" ]; then
    echo "Uporaba: bash tests/seal_extended.sh <KLJUC>" >&2
    echo "Ključ najdete v README.private.md." >&2
    exit 2
fi

if [ ! -f "$SRC" ]; then
    echo "[NAPAKA] Izvorne datoteke $SRC ni." >&2
    echo "         Odšifrirajte jo iz $ENC ali jo obnovite iz varnostne kopije." >&2
    exit 1
fi

# Preverimo, da je izvorna skripta sintaktično veljavna, preden jo zaklenemo —
# zaklenjene pokvarjene skripte se namreč ne da preprosto pregledati.
if ! bash -n "$SRC"; then
    echo "[NAPAKA] $SRC vsebuje sintaktično napako. Zaklepanje prekinjeno." >&2
    exit 1
fi

openssl enc -aes-256-cbc -pbkdf2 -a -salt -in "$SRC" -out "$ENC" -pass pass:"$KEY"

# Preverimo, da se da rezultat res odšifrirati nazaj v izvirnik.
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
