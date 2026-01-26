#!/usr/bin/env bash
# Utils/archive_daemon
#
# Archive daemon outputs for the active galaxy.
#
# Usage:
#   bash Utils/archive_daemon

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WHICH_FILE="$ROOT_DIR/which_galaxy"
DATA_ROOT="$ROOT_DIR/Data/Gal_Profiles"

if [ ! -f "$WHICH_FILE" ]; then
    echo "[ERROR] which_galaxy file not found"
    exit 1
fi

GALAXY="$(cat "$WHICH_FILE" | tr -d '[:space:]')"

if [ -z "$GALAXY" ]; then
    echo "[ERROR] which_galaxy is empty"
    exit 1
fi

PROFILE_DIR="$DATA_ROOT/$GALAXY"

if [ ! -d "$PROFILE_DIR" ]; then
    echo "[ERROR] Galaxy profile directory not found:"
    echo "        $PROFILE_DIR"
    exit 1
fi

ARCHIVE_ROOT="$PROFILE_DIR/archive"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
ARCHIVE_DIR="$ARCHIVE_ROOT/$TIMESTAMP"

mkdir -p "$ARCHIVE_DIR"

echo "[ARCHIVE] Galaxy: $GALAXY"
echo "[ARCHIVE] Destination: $ARCHIVE_DIR"
echo

shopt -s nullglob

FILES=(
    "$PROFILE_DIR"/daemon_deck*
    "$PROFILE_DIR"/daemon_log*
)

if [ ${#FILES[@]} -eq 0 ]; then
    echo "[ARCHIVE] No daemon files found to archive"
    exit 0
fi

for f in "${FILES[@]}"; do
    echo "[ARCHIVE] Moving $(basename "$f")"
    mv "$f" "$ARCHIVE_DIR/"
done

echo
echo "[ARCHIVE] Done"
