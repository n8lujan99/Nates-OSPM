#!/usr/bin/env bash
set -e

# Resolve directory of this script
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SRC="$HERE/Galaxies"
DST="$HERE/center"

mkdir -p "$DST"

tail -n +2 "$SRC" | while IFS=',' read -r name ra dec pa dist src; do
    d="$DST/$name"
    mkdir -p "$d"

    cat > "$d/center.txt" <<EOF
name = $name
ra_deg = $ra
dec_deg = $dec
pa_deg = $pa
distance_kpc = $dist
source = $src
EOF
done

