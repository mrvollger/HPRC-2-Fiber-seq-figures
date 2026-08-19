#!/usr/bin/env bash
# Per-sample PCA coordinates for the full 1000G + HGDP cohort.
# gnomAD v3.1 publishes the HGDP+1KG global PCA scores (4,117 samples, PC1-PC20)
# and a 198-column metadata table. We keep PC1-PC4 plus project/population/region,
# which is small enough to commit and covers 222 of the 232 HPRC2 samples directly
# (unlike the PCLAI reference space, which excludes HPRC2 samples by construction).
set -euo pipefail

BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2"
OUT="${1:-metadata/gnomad_hgdp_1kg_pca.tsv}"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

curl -sfL "$BASE/pca/pc_scores_with_outliers/GLOBAL_scores_with_outliers.txt.bgz" -o "$TMP/pcs.bgz"
curl -sfL "$BASE/metadata_and_qc/gnomad_meta_updated.tsv" -o "$TMP/meta.tsv"

# Column names are stable but not their positions; look them up from the header.
awk -F'\t' '
    NR == 1 {
        for (i = 1; i <= NF; i++) col[$i] = i
        print "sample_id\tproject\tpopulation\tregion"
        next
    }
    { print $col["s"] "\t" $col["hgdp_tgp_meta.Project"] "\t" \
            $col["hgdp_tgp_meta.Population"] "\t" $col["hgdp_tgp_meta.Genetic.region"] }
' "$TMP/meta.tsv" | sort -k1,1 >"$TMP/meta.small"

gunzip -c "$TMP/pcs.bgz" | awk -F'\t' 'NR > 1 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' |
    sort -k1,1 >"$TMP/pcs.small"

printf 'sample_id\tPC1\tPC2\tPC3\tPC4\tproject\tpopulation\tregion\n' >"$OUT"
join -t $'\t' -j 1 "$TMP/pcs.small" "$TMP/meta.small" >>"$OUT"

echo "wrote $(( $(wc -l <"$OUT") - 1 )) samples to $OUT" >&2
