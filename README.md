# HPRC-2-Fiber-seq-figures

Figures for the HPRC release 2 Fiber-seq pilot.

## Setup

```bash
pixi install     # builds the R environment (linux-64 and osx-64)
pixi run setup   # installs GenomeInfoDbData and mrvplot
```

Run `pixi run setup` once after `pixi install`. It does two things that
`pixi install` cannot:

1. Bioconductor *data* packages ship their R payload in a conda post-link
   script. Pixi does not run post-link scripts. Without this step
   `GenomeInfoDbData` is missing, and `GenomeInfoDb`, `GenomicRanges`, and
   `karyoploteR` all fail to load.
2. `mrvplot` is a GitHub package, so conda cannot provide it. It is installed
   with `dependencies = FALSE`, because every import it needs is already in the
   environment. Letting devtools resolve them makes it rebuild conda-managed
   packages from source, which fails.

The script is safe to re-run. It skips whatever already works.

## Figures

| Task | Script | Output |
|---|---|---|
| `pixi run hprc2-pca` | `Rscripts/hprc2-pca.Rmd` | `Figures/hprc2-pca.pdf` |
| `pixi run hsca-drivers` | `Rscripts/hsca-drivers-hla.Rmd` | `Figures/hsca-drivers-ideogram.pdf`, `Figures/hsca-drivers-denominators.pdf` |

Both read cached tables from `metadata/`, so neither needs cluster access.

Each figure writes the data behind it to `Figures/Tables/*.tbl.gz`.

Copy accepted figures to the manuscript folder on OneDrive:

```bash
pixi run publish hprc2-pca hsca-drivers-ideogram
```

## hprc2-pca.Rmd

Places the HPRC release 2 samples in the PCA space of HPRC2 Fig. 1d, and marks
the samples we have Fiber-seq data for.

The space is PCLAI's reference PCA: 1,561 low-admixture 1000G samples, two
haplotypes each. HPRC2 samples are **not** in that reference. PCLAI built it
from 1000G samples that were deliberately left out of HPRC2, so there is no
coordinate to look up for our samples. A haplotype's position here is the
genome-wide mean of PCLAI's per-window predictions.

That average is sound. It tracks an independent genotype PCA closely
(Procrustes R² = 0.97 against gnomAD's HGDP+1KG PCA over 222 shared samples),
and it does no worse for admixed samples than for single-ancestry ones. For a
recently admixed person the mean is a centroid between ancestry clusters, not a
place any window sits. `sd_x1` and `sd_x2` in `metadata/pclai_hap_means.tsv`
carry the per-window spread: 0.40 for admixed samples against 0.07 for
single-ancestry ones.

Haplotypes are drawn separately, joined by a segment. The segments are drawn
only if every Fiber-seq sample has both haplotypes. A partial set would read as
a real difference between samples rather than as missing data.

## hsca-drivers-hla.Rmd

Genome-wide CHM13 ideogram of phylo-FIRE driver peaks against all FIRE peaks.
It shows driver peaks piling up in the MHC, which motivates the HLA zoom-in.

**Coordinates are CHM13 v2.0**, with CHM13 cytoband stains.

**The MHC is GABBR1 through KIFC1**, chr6:29,430,312-33,231,258, 3.80 Mb. Gene
bounds come from the UCSC `hs1` `catLiftOffGenesV1` track. This matches how
driver elements are counted elsewhere in the project.

**Insertions are anchored, not dropped.** 50,862 peaks are insertions relative
to CHM13 and have no CHM13 coordinate of their own. Each one is placed at its
last preceding CHM13 peak, from Anna's `last_chm13_peak`. Insertions inside a
region therefore still count toward that region. Including them raises the MHC
odds ratio from 5.17 to 5.50.

Anchoring has one side effect. Where CHM13 has few peaks to anchor to, many
insertions collapse onto one anchor and make a spike in the background track.
This shows at the chr1 centromere, chr9, and the acrocentric chromosomes. Do
not read structure into those grey spikes.

**Densities slide.** The plotted density is a 5 Mb window that steps 500 kb.
Counts are binned at the step size, then summed over consecutive bins with
cumulative sums. Cumulative sums rather than `frollsum`, because `frollsum` has
to fill the leading bins, and filling them with zero makes the first 4.5 Mb of
every chromosome read as empty.

Each curve is scaled to sum to 1, so "all peaks" and "driver peaks" enclose the
same area. Height is then the share of that class found in a window, and the
two curves are comparable despite 617k against 3.2k peaks.

### Two tests

**A pre-specified Fisher test on the MHC.** This is what the claim rests on.

| denominator | MHC | rest | OR (95% CI) | p |
|---|---|---|---|---|
| HSCA peaks tested | 49/201 (24.4%) | 3175/22817 (13.9%) | 1.99 (1.41-2.77) | 8.6e-5 |
| all FIRE peaks | 49/1708 (2.9%) | 3175/593861 (0.5%) | 5.50 (4.04-7.32) | 2.1e-20 |
| HSCA tested / all peaks | 201/1708 (11.8%) | 22817/593861 (3.8%) | 3.34 (2.86-3.87) | 2.5e-43 |

The two denominators answer different questions. Against tested peaks: given a
haplotype-selective peak, is it more often sequence-driven in the MHC? Against
all peaks: does the MHC carry more sequence-driven HSCA than other regulatory
elements? The MHC is 3.3x enriched for haplotype-selective peaks to begin with,
and those peaks are then 2.0x more likely to be sequence-driven.

**A genome-wide scan in rolling windows of 200 tested peaks, stepping 50.** The
window is a fixed number of peaks, not a fixed number of base pairs. Peak
density varies by more than an order of magnitude along the genome, so fixed-bp
windows carry very unequal power. A 5 Mb window held anywhere from 20 to 200+
tested peaks. Equal-n windows make the p-values comparable.

This matters. Under a fixed-bp scan the MHC did not pass (FDR 0.29). Under
rolling peak-count windows it does, at FDR 1.1e-2.

The windows overlap, so BH here is a screen for where to look. It is not a
clean family-wise statement.

### Read this before using the chr19 result

The scan ranks two chr19 windows above the MHC. **Both are artifacts.**

chr19 as a whole has a 15.1% driver rate against a 13.9% genome baseline. Two
1 Mb bins carry the entire signal:

| bin (CHM13) | tested | drivers | rate |
|---|---|---|---|
| chr19:50-51 Mb | 28 | 19 | 67.9% |
| chr19:53-54 Mb | 30 | 15 | 50.0% |

Both are tandem repeat arrays. The gaps between consecutive driver peaks are
integer multiples of a ~5.3 kb unit (83% and 71% of gaps), and the peak widths
are nearly uniform. The MHC sits at 11% and chr1 at 7%, which is background.

Anna's `is_SD` flag agrees. Array 1 is 100% segmental duplication and array 2 is
76%, against 22% for the MHC and 10% for all tested peaks.

The scan treats each peak as independent. In those arrays that fails: 19 driver
calls in 199 kb are copies of one locus. Collapse each array to one observation
and chr19 has nothing left.

Full coordinates:

```
array 1  chr19:50,751,106-50,950,000  (199 kb)  19/21 drivers  unit ~5,323 bp
array 2  chr19:53,101,082-53,556,613  (456 kb)  15/30 drivers  unit ~5,338 bp
```

The hg38 liftover confirms this. Array 1 spans 199 kb in CHM13 but its peaks map
to 37 kb of hg38, and 17 of 21 peaks do not lift at all. That is an expanded
array against a reference that collapsed it.

## Data

Everything in `metadata/` is generated by a script in `scripts/` and committed
so the figures run without cluster access.

| File | Built by | Source |
|---|---|---|
| `pclai_hap_means.tsv` | `pclai-hap-means.sh` | PCLAI BEDs, public HPRC S3 |
| `pclai_reference_pca.tsv` | downloaded | `AI-sandbox/hprc-pclai` |
| `1kg_3202_samples_population.txt` | downloaded | 1000G FTP |
| `hsca_driver_peaks.tsv.gz` | `hsca-driver-peaks.sh` | phylo-FIRE `hsca_hapsel`, redwood |
| `all_peaks_chm13.tsv.gz` | `hsca-driver-peaks.sh` | consensus peaks, redwood |
| `peak_chm13_anchors.tsv.gz` | `hsca-driver-peaks.sh` | Anna's peak metadata, hyak |
| `chm13v2.0_cytobands.bed` | downloaded | T2T CHM13 S3 |
| `chm13v2.0.chrom.sizes` | downloaded | UCSC `hs1` |
| `hprc2.metadata.tbl` | provided | HPRC release 2 |
| `fiberseq.metadata.tbl` | provided | our Fiber-seq manifest |

`scripts/gnomad-1kg-hgdp-pcs.sh` builds the gnomAD HGDP+1KG PCA table. No
figure uses it. It is what produced the R² = 0.97 cross-check on the PCLAI
averages, and it regenerates that check on demand.

Downloaded cluster tables land in `phylo-fire-tables/`, which is gitignored.

### Going from DSA coordinates to a consensus peak

Query the union-peaks table. Column 5 is the consensus peak ID.

```bash
tabix union-peaks-asm.bed.gz "HG00438#1#CM089172.1:32480068-32480335"
```

A peak can exist in the consensus set while a given sample does not call it.
`is_peak` needs `fire_coverage >= 4` and `fire_coverage / coverage >= 0.2`.

## Gotchas

- **Anna's `peak_metadata_file.tsv` has 15,250 duplicated `consensus_peak_id`
  rows.** They differ only in `is_SD`. An unguarded join silently inflates
  counts. Deduplicate first.
- **`utils.R` loads IRanges, which masks `data.table::shift`.** Qualify the
  call.
- **`getCytobandColors()` already contains a `border` entry.** Appending a
  second one with `c()` does nothing, because `color.table["border"]` matches
  the first. Overwrite it in place. Black borders on 2 px bands render the whole
  ideogram solid black.
- **karyoploteR needs `chromosomes = "all"` with a custom genome.** Otherwise it
  filters heuristically and warns.
- **`pixi` uses `osx-64`, not `osx-arm64`.** `r-valr` and `r-weights` have no
  arm64 build. Under Rosetta the full dependency set solves unchanged.
