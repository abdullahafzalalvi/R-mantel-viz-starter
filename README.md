# R-mantel-viz-starter
R-mantel-viz-starter
# Mantel Test Visualization using linkET
# Demo script with synthetic data
 
# Author:  Abdullah Afzal Alvi
# Mantel Test Visualization with linkET

A reproducible R demo for performing pairwise Mantel tests between two variable blocks (physiological parameters and gene expression data) and visualizing the results using the [`linkET`](https://github.com/Hy4m/linkET) package.

---

## What This Script Does

1. Constructs a synthetic dataset with two treatment levels, four dose levels, and three biological replicates (24 total observations).
2. Defines two variable blocks:
   - **Physiological parameters** (`Param_1` to `Param_5`)
   - **Gene expression values** (`Gene_1` to `Gene_3`)
3. Performs pairwise Mantel tests (Pearson correlation on Euclidean distance matrices, 9999 permutations) for every parameter x gene combination.
4. Generates a combined visualization:
   - Lower-triangle **Pearson correlation heatmap** among physiological parameters
   - **Curved link lines** connecting gene expression variables to parameters, where line colour encodes Mantel *p*-value and line width encodes Mantel *r*

> **Note:** All data in this repository are synthetic and generated for demonstration purposes only. They do not represent real experimental measurements.

---

## Output

| File | Description |
|------|-------------|
| `outputs/mantel_linkET_demo_plot.png` | Main visualization (300 dpi) |
| `outputs/mantel_results_demo.csv` | Table of Mantel *r* and *p* values for all pairings |

### Example output figure

```
[Correlation heatmap of Param_1–Param_5] + [curved links to Gene_1, Gene_2, Gene_3]
```

---

## How to Use This for Your Own Data

Replace the synthetic data blocks in **Section 2** and **Section 3** of the script with your own measurements. The key requirements are:

- Both data frames must have the same number of rows (samples)
- A `Sample` column (or equivalent unique identifier) is used to merge the two blocks
- All variables passed to `env` should be on comparable scales, or standardize them with `scale()` as shown

---

## Installation

### R version

R >= 4.1.0 is recommended.

### CRAN packages

```r
install.packages(c("tidyverse", "RColorBrewer", "vegan"))
```

### linkET (GitHub)

```r
install.packages("devtools")
devtools::install_github("Hy4m/linkET")
```

---

## Usage

```bash
Rscript mantel_linkET_demo.R
```

Or open the script in RStudio and run it interactively.

---

## Dependencies

| Package | Version tested | Source |
|---------|---------------|--------|
| tidyverse | >= 1.3.0 | CRAN |
| linkET | >= 0.0.7 | GitHub (Hy4m/linkET) |
| RColorBrewer | >= 1.1-3 | CRAN |
| vegan | >= 2.6-4 | CRAN |

---

## Repository Structure

```
.
├── mantel_linkET_demo.R   # Main analysis script
├── README.md              # This file
└── outputs/               # Generated automatically on first run
    ├── mantel_linkET_demo_plot.png
    └── mantel_results_demo.csv
```

---

## Citation

If you adapt this script for your research, please cite the relevant packages:

- Oksanen J et al. (2022). *vegan: Community Ecology Package*. R package version 2.6-4. https://CRAN.R-project.org/package=vegan
- Hai Y (2023). *linkET: Everything is linkable*. https://github.com/Hy4m/linkET

---

## License

MIT License. See `LICENSE` for details.

---

## Contact

Contributions and issue reports are welcome via GitHub Issues.
