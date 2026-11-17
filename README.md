# HPRC-2-Fiber-seq-figures

Figures for the HPRC release 2 Fiber-seq pilot.

```bash
pixi install && pixi run setup   # R env; setup adds GenomeInfoDbData and mrvplot
pixi run <task>                  # tasks listed in pixi.toml
pixi run publish <figure> ...    # copy accepted figures to OneDrive
```

- `Rscripts/`: one `.Rmd` per figure; methods and caveats live there
- `scripts/`: builders for the cached tables in `metadata/`
- `metadata/`: committed inputs, so figures run without cluster access
- `Figures/`: PDFs, with backing data in `Figures/Tables/*.tbl.gz`
- `phylo-fire-tables/`: downloaded cluster tables (gitignored)
