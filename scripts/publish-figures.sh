#!/usr/bin/env bash
# Copy accepted figures to the manuscript folder on OneDrive.
# Usage: scripts/publish-figures.sh hprc2-pca [more-figures ...]
# Names may be bare (hprc2-pca) or paths (Figures/hprc2-pca.pdf); either way the
# backing data table mrv_ggsave() wrote is copied alongside the figure.
set -euo pipefail

DEST="${PUBLISH_DEST:-/Users/mrvollger/Library/CloudStorage/OneDrive-UW/Stergachislab/Manuscripts/HPRC - Fiber-seq pilot/Figures/MRV-raw-figures}"

[ $# -gt 0 ] || { echo "usage: $0 <figure> [figure ...]" >&2; exit 2; }
[ -d "$DEST" ] || { echo "destination not found: $DEST" >&2; exit 1; }

for arg in "$@"; do
    name=$(basename "$arg")
    name="${name%.pdf}"
    fig="Figures/${name}.pdf"
    tbl="Figures/Tables/${name}.tbl.gz"

    [ -f "$fig" ] || { echo "no such figure: $fig" >&2; exit 1; }
    cp "$fig" "$DEST/"
    echo "copied $fig"

    if [ -f "$tbl" ]; then
        mkdir -p "$DEST/Tables"
        cp "$tbl" "$DEST/Tables/"
        echo "copied $tbl"
    else
        echo "  note: no data table for $name" >&2
    fi
done
