# Post-shrinkage estimation for spatial autoregressive models

R code for the simulations and the two applications in the paper. 
The estimators (OF, UF, PT, S, PS) and the benchmarks (adaptive Lasso and ridge) are coded directly in `sar_functions.R`.

## Files

| File | Content |
|---|---|
| `sar_functions.R` | SAR-QMLE, information matrix, adaptive Lasso and ridge tuning, PT/S/PS, ADR formulas, missing-response likelihood and prediction |
| `real_data_functions.R` | Design and weights helpers for Ames and Chicago |
| `simulation.R` | Section 5: both experiments, Tables 1–2, Figures 1–4 and appendix figures |
| `chicago.R` | Section 6.1: Chicago income, Table 3, Figure 5, validation maps |
| `ames.R` | Section 6.2: Ames housing, Table 4, Figure 6, validation maps |
| `demo.R` | Small example with one simulated data set |
| `data/` | `ames.rda` and `Chi-SDOH.geojson` |

## How to run

Set the working directory to this folder and run, for example, `source("demo.R")`. Results are written to `output/`.

| Script | Time (one core) |
|---|---|
| `demo.R` | a few seconds |
| `chicago.R` | about 10 minutes |
| `ames.R` | about 25 minutes |
| `simulation.R` | about 2.5 hours (`mode <- "small"` for a quick check) |

Simulation cells and validation fold fits are saved and reused when a script is run again. All seeds are fixed in the scripts.

Packages: `Matrix`, `ggplot2`, and `sf` (Chicago). The OpenStreetMap basemap in the maps is optional and needs `ggspatial`, `rosm`, `prettymapr`, `raster` and an internet connection; without them the maps are drawn without a basemap.

## Notes

- Validation: preprocessing, selection, tuning and refitting use training data only. Held-out responses are integrated out of the Gaussian likelihood (Suesse, 2018).
- Prediction rule for PT, S and PS: `ames.R` uses `shrink_theta()`, which applies the shrinkage weight to all coefficients and to rho. In `chicago.R`, `rule <- "original"` gives Table 3 as reported; `rule <- "fixed"` uses the same rule as Ames.

## Data

- Ames housing: processed release in the R package `modeldata` (De Cock, 2011, *Journal of Statistics Education* 19). The 648 sales from 2009 are used.
- Chicago: Chi-SDOH data from the GeoDa Center, https://geodacenter.github.io/data-and-lab/Chi-SDOH/ (Kolak et al., 2020, *JAMA Network Open*).

Please cite the original data sources when using these files.
