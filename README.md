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
Comprehensive documentation covering data ingestion formats, normalization methods, statistical modeling, custom highlighting, and functional enrichment:

* **[View User Manual (Markdown)](USER_MANUAL.md)**
* **[Download User Manual (PDF)](USER_MANUAL.pdf)**

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

## License

This project is licensed under the MIT License - see the LICENSE file for details.
