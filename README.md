# ProteomicsAnalysisShiny

An interactive R Shiny application for differential expression analysis and visualization of quantitative proteomics data (Spectronaut, MSFragger, and generic quantitative tables).

---

## Features

- **Flexible Ingestion:** Direct support for Spectronaut exports, MSFragger outputs, and custom long/wide tabular formats.
- **Data Preprocessing:** Robust quantile normalization options and hybrid imputation (KNN for MAR / MinDet for MNAR).
- **Statistical Modeling:** Moderated linear models and empirical Bayes contrast evaluation via `limma`.
- **Interactive & Publication Visualizations:** Real-time Volcano, MA, and Scatter plots with dynamic `ggrepel` labeling, customizable FDR/logFC thresholds, and synchronized Plotly exploration.
- **Reproducible Exports:** Publication-grade PDF/PNG figure generation, multi-tab Excel workbooks, and complete JSON audit trails capturing session parameters, R versions, and loaded package namespaces.

---

## Installation & Setup

### Prerequisites
- [R (>= 4.2.0)](https://cran.r-project.org/) installed on your machine.
- An active internet connection for first-time package installation.

---

### Option 1: Quick Start (Windows Standalone ZIP)

1. Navigate to the [Releases](https://github.com/weaversd/ProteomicsAnalysisShiny/releases) page.
2. Download the latest `Source code (zip)` archive and extract it to your desired folder.
3. Double-click **`launch_app.bat`**. 
   - *The launcher will automatically check for required R packages, install any missing dependencies, and open the app in your default browser.*

---

### Option 2: Clone for Development

For developers who want the latest updates and active branch builds:

```bash
# 1. Clone repository
git clone https://github.com/weaversd/ProteomicsAnalysisShiny.git

# 2. Enter project directory
cd ProteomicsAnalysisShiny
```

Open `ProteomicsAnalysisShiny.Rproj` in RStudio, or launch the app directly in R:

```r
# Inside R / RStudio console
source("launch.R")
```

---

## File Format Requirements

### Spectronaut
- Standard `.tsv` or `.txt` output containing `PG.ProteinAccessions`, `PG.Quantity`, `R.Condition`, and `R.Replicate`.

### MSFragger
- Standard combined protein report (`.tsv`) containing `Protein.ID` and sample intensity columns (e.g., `x11.Intensity`).

### Generic Tabular (CSV / TSV)
- **Wide Format:** A row per protein with an `Accession` or `Protein.ID` column and sample columns ending in `.Intensity`.
- **Long Format:** Columns for `Accession`, `Sample` (or `ID`), and `Intensity` (or `Quantity`).

---

## Running Unit Tests

The test suite covers file parsing fallbacks, hybrid imputation edge cases, and end-to-end Shiny reactivity.

Run tests directly from R:

```r
testthat::test_dir("tests/testthat")
```

---

## Project Structure

```text
ProteomicsAnalysisShiny/
├── app.R                  # Shiny entrypoint
├── global.R               # Global variables, packages, and helper functions
├── launch.R               # Headless R launcher script
├── launch_app.bat         # 1-click Windows batch runner
├── modules/               # Modular UI/server components
├── R/
│   └── parsers.R          # Format parsing routines
├── server.R               # Main server logic and reactives
├── ui.R                   # Layout and interface definitions
└── tests/
    ├── testthat.R         # Test runner script
    └── testthat/          # Test suites (parsers, imputation, server)
```

---

## License

This project is licensed under the MIT License - see the LICENSE file for details.
