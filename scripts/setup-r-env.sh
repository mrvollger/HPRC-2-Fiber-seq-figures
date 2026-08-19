#!/usr/bin/env bash
# R setup that `pixi install` cannot do on its own. Run via `pixi run setup`.
#
# 1. Bioconductor *data* packages (GenomeInfoDbData and friends) ship their R
#    payload in a conda post-link script. Pixi does not run post-link scripts
#    (run-post-link-scripts = false), so the R package never lands in the
#    library and GenomeInfoDb -- and therefore karyoploteR and GenomicRanges --
#    fail to load with "there is no package called 'GenomeInfoDbData'".
#    We run those specific scripts rather than setting run-post-link-scripts
#    globally, which would execute every package's script.
# 2. mrvplot is a GitHub package, so conda cannot provide it. dependencies =
#    FALSE because every Import is already in the env; letting devtools resolve
#    them makes it rebuild conda-managed packages from source, which fails.
#
# Idempotent: re-running skips whatever already works.
set -euo pipefail

has_r_pkg() {
    Rscript -e "quit(status = as.integer(!requireNamespace('$1', quietly = TRUE)))" >/dev/null 2>&1
}

if has_r_pkg GenomeInfoDbData; then
    echo "GenomeInfoDbData present"
else
    echo "installing bioconductor data packages via their post-link scripts"
    shopt -s nullglob
    for script in "$CONDA_PREFIX"/bin/.bioconductor-*-post-link.sh; do
        echo "  $(basename "$script")"
        PREFIX="$CONDA_PREFIX" bash "$script"
    done
    shopt -u nullglob
fi

if has_r_pkg mrvplot; then
    echo "mrvplot present"
else
    echo "installing mrvplot from GitHub"
    Rscript -e 'devtools::install_github("mrvollger/mrvplot", dependencies = FALSE, upgrade = "never")'
fi

# Fail loudly if the env still cannot do what the figures need.
Rscript -e 'suppressMessages({library(karyoploteR); library(mrvplot)});
    cat("OK: karyoploteR", as.character(packageVersion("karyoploteR")),
        "| mrvplot", as.character(packageVersion("mrvplot")), "\n")'
