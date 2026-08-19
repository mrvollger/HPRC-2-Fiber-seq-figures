#!/usr/bin/env Rscript
# Reduce the phylo-FIRE run tables to two committable tables in CHM13
# (consensus/graph) coordinates. Called by scripts/hsca-driver-peaks.sh.
suppressMessages(library(data.table))

cache <- Sys.getenv("PHYLO_FIRE_CACHE", "phylo-fire-tables")

# Peaks that are insertions relative to CHM13 have no CHM13 coordinate of their
# own. Rather than dropping them, place them at their last preceding CHM13 peak
# (Anna's `last_chm13_peak`), so insertions inside a region still count toward
# that region. Built by scripts/hsca-driver-peaks.sh.
anchors <- fread("metadata/peak_chm13_anchors.tsv.gz")

# Give every peak a CHM13 position: its own if the consensus path is a CHM13
# chromosome, otherwise its anchor.
place <- function(dt, id_col) {
    on_chm13 <- grepl("^chr[0-9XY]+$", dt$chrom)
    out <- rbind(
        dt[on_chm13],
        merge(
            dt[!on_chm13][, c("chrom", "start", "end") := NULL],
            anchors[, .(consensus_peak_id, chrom, start, end)],
            by.x = id_col, by.y = "consensus_peak_id"
        ),
        use.names = TRUE
    )
    out[grepl("^chr[0-9XY]+$", chrom)]
}

# --- tested peaks + driver calls ------------------------------------------
# has_driver is written lowercase true/false, which fread reads as logical.
rv <- fread(
    file.path(cache, "re_variants.tsv.gz"),
    select = c("peak_id", "has_driver", "driver_ids")
)
stopifnot(is.logical(rv$has_driver))

# A peak is a driver peak if ANY of its RE sequences carries an attributed driver.
pk <- rv[, .(
    has_driver = any(has_driver),
    n_driver_ids = uniqueN(unlist(
        strsplit(driver_ids[!is.na(driver_ids) & driver_ids != ""], ",")
    ))
), by = peak_id]

co <- fread(file.path(cache, "peak_coords.tsv"))
tested <- merge(pk, co, by = "peak_id")[
    , .(peak_id,
        chrom = cons_chr, start = cons_start, end = cons_end,
        has_driver, n_driver_ids)
]
n_before <- nrow(tested)
tested <- place(tested, "peak_id")
setorder(tested, chrom, start)
fwrite(tested, "metadata/hsca_driver_peaks.tsv.gz", sep = "\t")
cat("tested peaks:", nrow(tested), "of", n_before,
    "| driver peaks:", sum(tested$has_driver), "\n")

# --- all consensus peaks, the background ----------------------------------
# consensus_peak_id encodes chrom_start_end for CHM13 paths. tstrsplit rather
# than a regex backreference, which does not survive a shell heredoc.
ids <- fread(
    file.path(cache, "all-peaks-cons-unique.txt.gz"),
    header = FALSE, col.names = "consensus_peak_id"
)
ids[, chrom := sub("_[0-9]+_[0-9]+$", "", consensus_peak_id)]
on_chm13 <- grepl("^chr[0-9XY]+$", ids$chrom)
ids[on_chm13, c("chrom", "start", "end") := tstrsplit(
    consensus_peak_id, "_",
    fixed = TRUE, type.convert = TRUE
)]
ids[!on_chm13, c("start", "end") := NA_integer_]
n_before <- nrow(ids)
all_peaks <- place(ids, "consensus_peak_id")
stopifnot(!anyNA(all_peaks$start), all(all_peaks$end > all_peaks$start))
setorder(all_peaks, chrom, start)
fwrite(all_peaks[, .(chrom, start, end)], "metadata/all_peaks_chm13.tsv.gz", sep = "\t")
cat("all consensus peaks:", nrow(all_peaks), "of", n_before,
    sprintf("(%d placed by anchor)\n", nrow(all_peaks) - sum(on_chm13)))
