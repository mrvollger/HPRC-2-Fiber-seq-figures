#!/usr/bin/env bash
# Tables behind Rscripts/hsca-drivers-hla.Rmd, pulled from redwood and reduced to
# something small enough to commit, so the Rmd needs no cluster access.
#
#   metadata/hsca_driver_peaks.tsv.gz  one row per peak phylo-FIRE tested, with
#                                      has_driver from forward selection
#   metadata/all_peaks_chm13.tsv.gz    every consensus FIRE peak, the background
#
# Coordinates are CHM13 (the consensus/graph reference), matching the ideogram.
# The all-peaks set comes from the unique consensus_peak_id values, which encode
# chrom_start_end; IDs on per-sample contigs (HG00438#1#CM089167.1_...) are
# dropped, leaving the bare chr* CHM13 reference paths.
set -euo pipefail

PF="${PHYLO_FIRE_RUN:-/scratch/ucgd/lustre-labs/vollger/projects/HPRCv2/phylo-FIRE/results/hsca_hapsel/tables}"
FIG="${MRV_FIGURES_TABLES:-/scratch/ucgd/lustre-labs/vollger/projects/HPRCv2/mrv-figures/Tables}"
HOST="${PHYLO_FIRE_HOST:-redwood}"
CACHE="${PHYLO_FIRE_CACHE:-phylo-fire-tables}"
ANNA_HOST="${ANNA_HOST:-hyak}"
ANNA_META="${ANNA_META:-/mmfs1/gscratch/stergachislab/HPRC/fiber-seq-pilot/HLA_graphs/peak_calling_all_chrs/peak_calling_downstream_analyses/peak_metadata_file.tsv}"

mkdir -p "$CACHE" metadata

# The all-peaks reduction streams a 700 MB table, so do it on the cluster and
# only ever download the ~600k unique ids.
ssh "$HOST" "test -s '$FIG/all-peaks-cons-unique.txt.gz' || \
    gzip -dc '$FIG/peaks-cons.bed.gz' | awk -F'\t' 'NR>1{print \$5}' | sort -u | \
    gzip > '$FIG/all-peaks-cons-unique.txt.gz'"

rsync -a --no-motd "$HOST:$PF/{re_variants.tsv.gz,peak_coords.tsv}" "$CACHE/"
rsync -a --no-motd "$HOST:$FIG/all-peaks-cons-unique.txt.gz" "$CACHE/"

# Anchors for peaks that are insertions relative to CHM13: Anna's peak metadata
# carries last_chm13_peak, the preceding CHM13 peak. Lives on hyak, not redwood.
# Rows there are duplicated per is_SD value, hence the unique().
if [ ! -s metadata/peak_chm13_anchors.tsv.gz ]; then
    rsync -a --no-motd "$ANNA_HOST:$ANNA_META" "$CACHE/peak_metadata_file.tsv"
    Rscript -e '
        suppressMessages(library(data.table))
        m <- fread(file.path(Sys.getenv("PHYLO_FIRE_CACHE", "phylo-fire-tables"),
                             "peak_metadata_file.tsv"),
                   select = c("consensus_peak_id", "in_chm13_asm", "last_chm13_peak"))
        a <- unique(m[in_chm13_asm == FALSE, .(consensus_peak_id, last_chm13_peak)])
        stopifnot(uniqueN(a$consensus_peak_id) == nrow(a))
        p <- a[, tstrsplit(last_chm13_peak, "_", fixed = TRUE, type.convert = TRUE)]
        a[, `:=`(chrom = p$V1, start = p$V2, end = p$V3)]
        a <- a[grepl("^chr[0-9XY]+$", chrom)]
        setorder(a, chrom, start)
        fwrite(a[, .(consensus_peak_id, last_chm13_peak, chrom, start, end)],
               "metadata/peak_chm13_anchors.tsv.gz", sep = "\t")
        cat("anchors:", nrow(a), "\n")
    '
fi

Rscript scripts/hsca-driver-peaks.R
