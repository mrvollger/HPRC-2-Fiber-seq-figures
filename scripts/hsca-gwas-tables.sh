#!/usr/bin/env bash
# Tables behind Rscripts/hsca-gwas-enrichment.Rmd, reduced so the Rmd itself
# only reads small files.
#
#   Tables/PICS2-GWAScat-2024-09-20.txt.gz  PICS2 fine-mapped GWAS Catalog, the
#                                           release used for the FIRE paper
#                                           (664 MB, gitignored, downloaded)
#   Tables/pics2-gwas-sites.tsv.gz          one row per GRCh38 position linked to
#                                           a genome-wide significant association
#   metadata/all_peaks_hg38.tsv.gz          every consensus FIRE peak with a
#                                           GRCh38 interval
#   metadata/hg38.chrom.sizes               chr1-22, X, Y, M (for the shuffles)
#
# The driver calls (metadata/hsca_driver_peaks.tsv.gz) come from
# scripts/hsca-driver-peaks.sh.
set -euo pipefail

PICS_URL="https://pics2.ucsf.edu/Downloads/PICS2-GWAScat-2024-09-20.txt.gz"
PICS_RAW="Tables/PICS2-GWAScat-2024-09-20.txt.gz"
PICS_SHA256="c766e042b4f22b3156e6f21a50b1ac7fc7de5b0167851db2ecca4a73649dd97a"
UNION="${UNION_PEAKS:-Tables/union-peaks-asm.bed.gz}"
# genome-wide significance, p < 5e-8, as in utilities/intersect-with-gwas.py
MIN_MLOGP="${MIN_MLOGP:-7.30103}"
THREADS="${THREADS:-8}"

mkdir -p Tables metadata

# 1. PICS2 release. Checked against the hash recorded when it was first fetched.
[ -s "$PICS_RAW" ] || curl -sS -L -o "$PICS_RAW" "$PICS_URL"
echo "$PICS_SHA256  $PICS_RAW" | shasum -a 256 -c -

# 2. GWAS sites. PICS2 lists, for each index (lead) SNP, every SNP linked to it
# with a PICS probability. Keep associations with p < 5e-8 and collapse to one
# row per GRCh38 position. Columns used (tab-separated, 1-based):
#   2 IndexSNP  3 LinkedSNP  4 "chr:pos (hg19 chr:pos)"  11 PICS_probability
#   21 Index-mlogP
gunzip -c "$PICS_RAW" | awk -F'\t' -v OFS='\t' -v min="$MIN_MLOGP" '
    NR == 1 { next }
    $21 + 0 > min {
        split($4, hg38, " ")
        n = split(hg38[1], cp, ":")
        if (n != 2 || cp[2] !~ /^[0-9]+$/) next
        key = "chr" cp[1] "\t" cp[2]
        n_assoc[key]++
        if ($11 + 0 > pics[key]) pics[key] = $11 + 0
        if ($21 + 0 > mlogp[key]) mlogp[key] = $21 + 0
        if ($2 == $3) is_index[key] = 1
    }
    END {
        for (key in n_assoc)
            print key, pics[key] + 0, mlogp[key], n_assoc[key], is_index[key] + 0
    }' |
    sort -k1,1V -k2,2n |
    awk 'BEGIN { print "chrom\tpos\tmax_pics\tmax_mlogp\tn_assoc\tis_index" } 1' |
    bgzip -@ "$THREADS" >Tables/pics2-gwas-sites.tsv.gz

# 3. Every consensus FIRE peak with a GRCh38 interval: a peak in at least one
# primary sample. GRCh38 coordinates are constant across the rows of a peak, so
# the first row wins (the rule build_peak_coords.py uses in phylo-FIRE). hck
# returns columns in file order, hence the lookup by header name.
bgzip -dc -@ "$THREADS" "$UNION" |
    hck -d '\t' -F consensus_peak_id -F hg38_chr -F hg38_start -F hg38_end \
        -F exists_in_hg38 -F is_peak -F is_primary_sample |
    awk -F'\t' -v OFS='\t' '
        NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
        $c["is_peak"] == "TRUE" && $c["is_primary_sample"] == "TRUE" &&
            $c["exists_in_hg38"] == "TRUE" && !seen[$c["consensus_peak_id"]]++ {
            print $c["consensus_peak_id"], $c["hg38_chr"], $c["hg38_start"], $c["hg38_end"]
        }' |
    sort -k2,2V -k3,3n |
    awk 'BEGIN { print "consensus_peak_id\tchrom\tstart\tend" } 1' |
    bgzip -@ "$THREADS" >metadata/all_peaks_hg38.tsv.gz

# 4. GRCh38 chromosome sizes for bed_shuffle.
curl -sS https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes |
    awk -F'\t' '$1 ~ /^chr([0-9]+|X|Y|M)$/' | sort -k1,1V >metadata/hg38.chrom.sizes

echo "GWAS sites: $(($(bgzip -dc Tables/pics2-gwas-sites.tsv.gz | wc -l) - 1))"
echo "FIRE peaks on GRCh38: $(($(bgzip -dc metadata/all_peaks_hg38.tsv.gz | wc -l) - 1))"
