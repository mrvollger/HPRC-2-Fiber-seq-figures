#!/usr/bin/env bash
# GWAS risk alleles carried by each Fiber-seq haplotype, and each haplotype's
# FIRE peaks on GRCh38. Input for Rscripts/fiber-gwas-enrichment.Rmd.
#
#   Tables/fiber-gwas-index-sites.tsv.gz   PICS2 index SNPs and risk alleles
#   Tables/fiber-gwas-risk-alleles.tsv.gz  risk alleles carried, per haplotype
#   Tables/fiber-peaks-hg38.tsv.gz         FIRE peaks, per haplotype
#
# Risk alleles are matched by GRCh38 position and allele (REF or ALT), trying
# the complement strand except at A/T and C/G sites.
set -euo pipefail

S3="https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/release2/minigraph-cactus"
VCF="Tables/hprc-v2.0-mc-grch38.wave.vcf.gz"
PICS_RAW="Tables/PICS2-GWAScat-2024-09-20.txt.gz" # from scripts/hsca-gwas-tables.sh
MIN_MLOGP="${MIN_MLOGP:-0}"                       # -log10(p); 7.30103 is p < 5e-8
THREADS="${THREADS:-8}"
SITES="Tables/fiber-gwas-index-sites.tsv.gz"
OUT="Tables/fiber-gwas-risk-alleles.tsv.gz"

# HPRC v2.0 Minigraph-Cactus VCF
[ -s "$VCF" ] || curl -sS -L -o "$VCF" "$S3/hprc-v2.0-mc-grch38.wave.vcf.gz"
[ -s "$VCF.tbi" ] || curl -sS -L -o "$VCF.tbi" "$S3/hprc-v2.0-mc-grch38.wave.vcf.gz.tbi"

# PICS2 index SNPs. Columns: 2 IndexSNP, 3 LinkedSNP, 4 GRCh38 position,
# 19 SNP-riskAllele, 21 Index-mlogP
gunzip -c "$PICS_RAW" | awk -F'\t' -v OFS='\t' -v min="$MIN_MLOGP" '
    NR == 1 || $2 != $3 || $21 + 0 < min { next }
    {
        split($4, hg38, " ")
        if (split(hg38[1], cp, ":") != 2 || cp[2] !~ /^[0-9]+$/) next
        risk = toupper(substr($19, index($19, "-") + 1))
        if (risk !~ /^[ACGT]$/) next
        key = "chr" cp[1] OFS cp[2] OFS $2 OFS risk
        n[key]++
        if ($21 + 0 > mlogp[key]) mlogp[key] = $21 + 0
    }
    END { for (key in n) print key, mlogp[key], n[key] }' |
    sort -k1,1V -k2,2n |
    awk 'BEGIN { print "chrom\tpos\trsID\trisk_allele\tmax_mlogp\tn_assoc" } 1' |
    bgzip -@ "$THREADS" >"$SITES"

# Phased genotypes of the Fiber-seq individuals at those sites. GM12878 is not
# in the HPRC VCF. Top-level records only (LV=0).
samples=$(awk -F'\t' 'NR > 1 && $2 != "GM12878" { print $2 }' metadata/fiberseq.metadata.tbl |
    sort -u | paste -sd, -)

bcftools view --threads "$THREADS" -s "$samples" -e 'INFO/LV>0' \
    -T <(bgzip -dc "$SITES" | awk -F'\t' -v OFS='\t' 'NR > 1 { print $1, $2 }') \
    "$VCF" |
    bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' |
    awk -F'\t' -v OFS='\t' -v sites=<(bgzip -dc "$SITES") '
        function comp(b) { return b == "A" ? "T" : b == "T" ? "A" : b == "C" ? "G" : "C" }
        BEGIN {
            while ((getline line < sites) > 0) {
                split(line, f, "\t")
                if (f[1] == "chrom") continue
                k = f[1] ":" f[2]
                site[k] = site[k] f[3] "\t" f[4] "\n"
            }
        }
        NR == 1 {
            for (i = 5; i <= NF; i++) { s = $i; sub(/^\[[0-9]+\]/, "", s); sub(/:GT$/, "", s); name[i] = s }
            next
        }
        length($3) != 1 { next }
        {
            k = $1 ":" $2
            if (!(k in site)) next
            nalt = split($4, alt, ",")
            allele[0] = $3
            for (a = 1; a <= nalt; a++) allele[a] = substr(alt[a], 1, 1)
            m = split(site[k], rows, "\n")
            for (r = 1; r < m; r++) {
                split(rows[r], rr, "\t")
                risk = rr[2]; flipped = 0
                present = (risk == $3)
                for (a = 1; a <= nalt && !present; a++) present = (risk == allele[a])
                if (!present) {
                    c = comp(risk); hit = (c == $3); amb = 0
                    for (a = 1; a <= nalt; a++) {
                        if (allele[a] == c) hit = 1
                        if (allele[a] == comp($3)) amb = 1
                    }
                    if (!hit || amb) continue
                    risk = c; flipped = 1
                }
                for (i = 5; i <= NF; i++) {
                    nh = split($i, gt, /[|\/]/)
                    for (h = 1; h <= nh; h++) {
                        if (gt[h] == "." || allele[gt[h]] != risk) continue
                        print name[i], h, $1, $2, rr[1], risk, (risk == $3) + 0, flipped
                    }
                }
            }
        }' |
    sort -k1,1 -k3,3V -k4,4n -k2,2n -k5,5 | uniq |
    awk 'BEGIN { print "Individual_ID\thaplotype\tchrom\tpos\trsID\trisk_allele\trisk_is_ref\tstrand_flipped" } 1' |
    bgzip -@ "$THREADS" >"$OUT"

# FIRE peaks of each primary sample on GRCh38, chromosomes in hg38.chrom.sizes
bgzip -dc -@ "$THREADS" "${UNION_PEAKS:-Tables/union-peaks-asm.bed.gz}" |
    hck -d '\t' -F hg38_chr -F hg38_start -F hg38_end -F exists_in_hg38 -F is_peak \
        -F Individual_ID -F Haplotype -F Platform -F is_primary_sample |
    awk -F'\t' -v OFS='\t' '
        NR == FNR { keep[$1] = 1; next }
        FNR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
        $c["is_peak"] == "TRUE" && $c["is_primary_sample"] == "TRUE" &&
            $c["exists_in_hg38"] == "TRUE" && ($c["hg38_chr"] in keep) {
            print $c["Individual_ID"], $c["Platform"], $c["Haplotype"],
                $c["hg38_chr"], $c["hg38_start"], $c["hg38_end"]
        }' metadata/hg38.chrom.sizes - |
    sort -k1,1 -k3,3n -k4,4V -k5,5n -u |
    awk 'BEGIN { print "Individual_ID\tPlatform\thaplotype\tchrom\tstart\tend" } 1' |
    bgzip -@ "$THREADS" >Tables/fiber-peaks-hg38.tsv.gz

echo "PICS2 index SNP risk alleles: $(($(bgzip -dc "$SITES" | wc -l) - 1))"
echo "carried (haplotype x SNP): $(($(bgzip -dc "$OUT" | wc -l) - 1))"
