#!/usr/bin/env bash
# Genome-wide mean PCLAI coordinate per HPRC R2 haplotype.
# Streams the 463 per-haplotype BEDs from the public HPRC S3 bucket and averages
# the (PC1,PC2) in column 4 across windows. The per-window SD comes along too: it
# says whether the mean is a faithful summary (tight, unadmixed) or a centroid
# sitting between ancestry clusters (wide, admixed). Output is a ~463-row table
# small enough to commit, so the Rmd never touches S3.
set -euo pipefail

IDX_URL="https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/main/data_tables/annotation/pclai/pclai_v1.1_grch38_coord_local_hprc_r2.index.csv"
OUT="${1:-metadata/pclai_hap_means.tsv}"
JOBS="${2:-12}"

one() {
    local sample=$1 hap=$2 s3=$3
    local url="https://human-pangenomics.s3.amazonaws.com/${s3#s3://human-pangenomics/}"
    curl -sfL --retry 3 --max-time 600 "$url" |
        awk -F'\t' -v s="$sample" -v h="$hap" '
            { p = index($4, "("); split(substr($4, p + 1, length($4) - p - 1), c, ",");
              x += c[1]; y += c[2]; xx += c[1] ^ 2; yy += c[2] ^ 2; q += $5; n++ }
            END {
                if (!n) exit
                mx = x / n; my = y / n
                sx = sqrt((xx / n - mx ^ 2) > 0 ? xx / n - mx ^ 2 : 0)
                sy = sqrt((yy / n - my ^ 2) > 0 ? yy / n - my ^ 2 : 0)
                printf "%s\t%s\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.1f\n", s, h, n, mx, my, sx, sy, q / n
            }
        '
}
export -f one

printf 'sample_id\thaplotype\tn_windows\tx1\tx2\tsd_x1\tsd_x2\tmean_score\n' >"$OUT"
curl -sfL "$IDX_URL" | tr -d '\r' | tail -n +2 |
    awk -F, '{print $1, $2, $4}' |
    xargs -P "$JOBS" -n 3 bash -c 'one "$0" "$1" "$2"' >>"$OUT"

echo "wrote $(( $(wc -l <"$OUT") - 1 )) haplotypes to $OUT" >&2
