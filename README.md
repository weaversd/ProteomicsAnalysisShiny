# ProteomicsAnalysisShiny

An interactive R Shiny application for end-to-end differential expression analysis, interactive visualization, and functional pathway enrichment (KEGG & GSEA) of quantitative proteomics data.

---

## Features

- **Flexible Ingestion:** Direct support for Spectronaut, MSFragger, DIA-NN, PEAKS, generic long/wide tables, and serialized app session states (.rds).
- **Pre-Normalization & Global Normalization:** Optional loading control/reference protein normalization (e.g., Gapdh, Actb) prior to global methods (median centering, mean centering, quantile normalization, or VSN).
- **Preprocessing & Imputation:** Hybrid imputation (KNN for MAR / MinDet for MNAR) alongside a dedicated 'No Imputation' mode projecting single-condition dropouts directly onto MA and Scatter margins.
- **Statistical Modeling:** Moderated linear models and empirical Bayes contrast evaluation via limma, with instant contrast inversion toggles.
- **Interactive Visualizations & Custom Highlighting:** Volcano, MA, and Scatter plots rendered synchronously via ggplot2 and Plotly. Dynamically highlight and label custom protein sets using regex patterns or pasted accession/gene lists with independent color and label toggles.
- **Functional Pathway Enrichment (KEGG & GSEA):** Over-Representation Analysis (ORA) and ranked Gene Set Enrichment Analysis (GSEA) via clusterProfiler for Mouse (Mus musculus) and Human (Homo sapiens). Includes interactive enrichment dot plots, running-score pathway curves (gseaplot2), and distribution overlays showing pathway proteins on Volcano/MA/Scatter plots.
- **Reproducible Exports & State Tracking:** Multi-tab Excel workbooks (.xlsx) with condition-specific headers, publication-grade PNG/PDF figure exports, standalone reproducible R analysis scripts, JSON audit trails capturing Git tags and R environment parameters, and RDS full-state saving/restoring.

---
## User Manual:
- Available as USER_MANUAL.pdf
---

---
## Installation & Setup

### Prerequisites
- R (>= 4.2.0) installed on your system.
- Bioconductor packages (QFeatures, limma, clusterProfiler, enrichplot, org.Mm.eg.db, org.Hs.eg.db).
- An active internet connection for first-time dependency resolution.

---

### Option 1: Quick Start (Windows Standalone ZIP)

1. Navigate to the Releases page.
2. Download the latest Source code (zip) archive and extract it to your desired location.
3. Double-click launch_app.bat.
   - The launcher verifies package dependencies, installs any missing CRAN or Bioconductor packages, and launches the app in your default browser.

---

### Option 2: Clone for Development

For developers working directly on features or running unit test suites:

```bash
git clone https://github.com/weaversd/ProteomicsAnalysisShiny.git
cd ProteomicsAnalysisShiny
```

Open ProteomicsAnalysisShiny.Rproj in RStudio, or execute the launch script from an R console:

```r
source("launch.R")
```

---

## File Format & Ingestion Requirements

Select the matching engine from the Select Source Format dropdown on the import panel:

| Platform | Expected File Name / Type | Key Columns Required | Notes |
| :--- | :--- | :--- | :--- |
| **Spectronaut** | *.tsv or *.txt (Export Report) | PG.ProteinAccessions, PG.Quantity, R.Condition, R.Replicate | Quantities <= 0 are converted to missing values (NA). |
| **MSFragger** | combined_protein.tsv | Protein.ID (or Protein), along with * Intensity or * MaxLFQ Intensity columns | Prioritizes standard intensity columns; falls back to MaxLFQ intensity if standard is absent. |
| **DIA-NN** | report.pg_matrix.tsv | Protein.Group, Genes, First.Protein.Description, and sample file path columns | Raw run file paths (e.g., *.wiff, *.raw) are stripped of directory structures and common run prefixes/suffixes to yield clean sample IDs. |
| **PEAKS** | db.proteins.csv or lfq.proteins.csv | Accession, Description, and Area Sample * or Sample * Area columns | Accessions with trailing pipes (e.g., P12273|PIP_HUMAN) are cleaned to the base accession. Gene symbols are extracted from the GN= description field. |
| **Generic Table** | *.csv or *.tsv (Wide or Long format) | Wide: Accession or Protein.ID with .Intensity numeric sample columns.<br>Long: Accession, Sample/ID, and Intensity/Quantity. | Supports custom workflows from MaxQuant, Skyline, or custom scripting pipelines. |
| **Saved State** | *_state.rds | Serialized application state list | Re-imports complete raw data, sample maps, normalized matrices, and DE calculations without re-running models. |

---
## Functional Enrichment Analysis (KEGG & GSEA)

The KEGG & GSEA Enrichment panel allows functional profiling directly from within the app:

- **Dual Ingestion:** Run enrichment on the active differential expression results calculated in the app or upload an exported Excel workbook (.xlsx).
- **Directional Awareness:** Inverting the contrast direction on the visualization tab propagates through to enrichment, reversing log2FC signs to align upregulated and downregulated subsets properly.
- **KEGG Over-Representation Analysis (ORA):** Filter proteins using user-defined adjusted p-value and log2FC cutoffs, choose protein directions (All Significant, Upregulated, or Downregulated), and render customizable enrichment dot plots.
- **Ranked GSEA:** Evaluates global distributions across ranked log2FC profiles, generates split-sign activation/suppression dot plots, and renders running enrichment score curves (gseaplot2) for individual selected pathways.
- **Pathway DE Distribution View:** Plots proteins belonging to any selected GSEA pathway onto Volcano, MA, or Scatter plots to inspect the specific fold-changes and significance levels of genes driving pathway enrichment.

---

## Running Unit Tests

The test suite covers file format parsing, reference protein adjustments, hybrid imputation edge cases, and contrast inversion handlers:

```r
testthat::test_dir("tests/testthat")
```

---

## Project Structure

```text
ProteomicsAnalysisShiny/
|-- app.R                  # Shiny entrypoint
|-- global.R               # Environment setup, package imports, and version tracking
|-- launch.R               # Standalone R launch script
|-- launch_app.bat         # 1-click Windows batch runner
|-- modules/               # Modular Shiny UI and server modules
|-- R/
|   `-- parsers.R          # Spectronaut, MSFragger, DIA-NN, PEAKS, and Generic parsers
|-- server.R               # Core reactive pipeline, modeling, visualization, and export logic
|-- ui.R                   # Layout, CSS styling, and tab definitions
`-- tests/
    |-- testthat.R         # Test runner script
    `-- testthat/          # Testthat test suites (parsers, normalization, stats)
```

---

## License

This project is licensed under the MIT License - see the LICENSE file for details.
