# Development scripts

These scripts are not part of the built package (`dev/` is listed in
`.Rbuildignore`). Run them from the directory that contains the package
folder, for example `Rscript mglpa/dev/test_app.R`.

| Script | Purpose |
|---|---|
| `validate_engine.R` | Compares the EM engine with the Mplus reference fits of the MG-LPA Analyzer example project (`inst/extdata/mplus_reference_fits.csv`). |
| `validate_hard.R` | Solutions that need several thousand random starts to reach the global maximum. |
| `e2e.R` | End-to-end run of the whole workflow on the example data. |
| `test_app.R` | Headless test of the Shiny application (`shiny::testServer`). |
| `screenshots.R` | Drives the installed application in headless Chrome (chromote package) and writes the README screenshots to `man/figures/`; doubles as an end-to-end check of the user interface. The committed images were afterwards reduced to a 256-colour palette. |
| `test_parallel.R` | Timing of parallel random starts (needs the installed package). |
| `mplus_reference_extra.R` | Writes and runs the additional Mplus inputs in `mplus_extra/` (class-varying variances, three known groups, missing data) and updates the reference values; needs a local Mplus installation. |

`mplus_extra/` holds the Mplus 9.1.1 inputs and outputs behind those
reference values, including the `TECH11` runs (`t11_*`) used to verify
`vlmr_test()`.
