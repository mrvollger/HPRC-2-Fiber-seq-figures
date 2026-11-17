#!/usr/bin/env Rscript
# Are non-reference peaks (insertions relative to CHM13) enriched for segmental
# duplications? Reviewer point 4 on Fig. 6i. Runs on Anna's peak metadata, cached
# by scripts/hsca-driver-peaks.sh in phylo-fire-tables/.
suppressMessages(library(data.table))

m <- fread("phylo-fire-tables/peak_metadata_file.tsv")
# Rows are duplicated per is_SD value (README gotcha): a peak is SD if any row says so.
pk <- m[, .(
    is_SD = any(is_SD), nonref = !in_chm13_asm[1], asm = cons_peak_asm[1],
    anchor = last_chm13_peak[1], n_asm = num_of_asms_w_peak[1]
), by = consensus_peak_id]

fisher_row <- function(dt, label) {
    tab <- dt[, table(factor(nonref, c(TRUE, FALSE)), factor(is_SD, c(TRUE, FALSE)))]
    ft <- fisher.test(tab)
    data.table(
        test = label,
        nonref_SD = sprintf("%d/%d (%.1f%%)", tab[1, 1], sum(tab[1, ]), 100 * tab[1, 1] / sum(tab[1, ])),
        ref_SD = sprintf("%d/%d (%.1f%%)", tab[2, 1], sum(tab[2, ]), 100 * tab[2, 1] / sum(tab[2, ])),
        OR = round(unname(ft$estimate), 1),
        CI = sprintf("%.1f-%.1f", ft$conf.int[1], ft$conf.int[2]),
        p = signif(ft$p.value, 3)
    )
}

# 1. Every peak as an observation.
# 2. Non-reference peaks collapsed to their CHM13 anchor locus, so SD copies of
#    one insertion count once. The locus is SD if most of its peaks are.
loci <- rbind(
    pk[nonref == FALSE, .(nonref, is_SD)],
    pk[nonref == TRUE, .(nonref = TRUE, is_SD = mean(is_SD) > 0.5), by = anchor][, .(nonref, is_SD)]
)
res <- rbind(
    fisher_row(pk, "all peaks"),
    fisher_row(pk[nonref == FALSE | n_asm %in% 1], "single-assembly nonref only"),
    fisher_row(loci, "nonref collapsed to CHM13 anchor")
)
print(res)

# 3. Against the genome: SDs are ~7% of CHM13 (Vollger et al. 2022, Science).
nr <- pk[nonref == TRUE]
bt <- binom.test(sum(nr$is_SD), nrow(nr), 0.07, alternative = "greater")
cat(sprintf("\nnonref SD %.1f%% vs 7%% of genome: p = %.3g\n", 100 * mean(nr$is_SD), bt$p.value))

# 4. Consistency: SD fraction of non-reference peaks per assembly.
per_asm <- nr[!asm %in% "CHM13", .(n = .N, pct_SD = 100 * mean(is_SD)), by = asm][n >= 50]
cat(sprintf(
    "per assembly (n>=50 nonref peaks, %d assemblies): median %.1f%%, range %.1f-%.1f%%; all above reference %.1f%%: %s\n",
    nrow(per_asm), median(per_asm$pct_SD), min(per_asm$pct_SD), max(per_asm$pct_SD),
    100 * pk[nonref == FALSE, mean(is_SD)], all(per_asm$pct_SD > 100 * pk[nonref == FALSE, mean(is_SD)])
))

dir.create("Figures/Tables", showWarnings = FALSE, recursive = TRUE)
fwrite(res, "Figures/Tables/nonref-peaks-sd-test.tbl.gz", sep = "\t")
fwrite(per_asm, "Figures/Tables/nonref-peaks-sd-per-asm.tbl.gz", sep = "\t")
