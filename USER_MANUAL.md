# Proteomics Explorer: User Manual & Technical Guide

**Version:** 1.2.0
**Author:** Simon Weaver (Notre Dame Mass Spectrometry and Proteomics Facility)
**License:** MIT

---

## Table of Contents
1. [Overview & Architecture](#1-overview--architecture)
2. [Prerequisites & Launch Options](#2-prerequisites--launch-options)
3. [Tab 1: Data Import & Metadata Setup](#3-tab-1-data-import--metadata-setup)
4. [Tab 2: Normalization & Imputation](#4-tab-2-normalization--imputation)
5. [Tab 3: Differential Expression & Visualization](#5-tab-3-differential-expression--visualization)
6. [Tab 4: Functional Pathway Enrichment (KEGG & GSEA)](#6-tab-4-functional-pathway-enrichment-kegg--gsea)
7. [Tab 5: Export & Reproducibility Audit](#7-tab-5-export--reproducibility-audit)
8. [Troubleshooting & FAQs](#8-troubleshooting--faqs)

---

## 1. Overview & Architecture

Proteomics Explorer is an interactive R/Shiny platform built for quantitative bottom-up and label-free mass spectrometry (LFQ/DIA) analysis. It integrates data ingestion across major search platforms, customizable normalization pipelines, missing value imputation, limma-powered moderated linear modeling, and functional pathway enrichment via clusterProfiler.

### Pipeline Flow

    Raw Search Engine / Matrix Output
                   |
                   v
    [Tab 1] Ingestion & Parsing ---> Sample Mapping (Bio/Tech Reps) ---> Technical Replicate Averaging (Log2)
                                                                                |
                   +------------------------------------------------------------+
                   v
    [Tab 2] Reference Alignment (Optional) ---> Global Normalization (QFeatures) ---> Imputation
                                                                                          |
                   +----------------------------------------------------------------------+
                   v
    [Tab 3] Moderated Linear Fit (limma) ---> Pairwise Contrasts ---> Volcano / MA / Scatter Visualizations
                                                                                |
                   +------------------------------------------------------------+---------------------+
                   v                                                                                  v
    [Tab 4] Functional Enrichment (KEGG ORA & GSEA)                                  [Tab 5] Export Bundle & Audit
            * Dot Plots                                                                      * Multi-Tab Excel (.xlsx)
            * Running-Score Plots (gseaplot2)                                                * Standalone R Script
            * Pathway DE Distribution Overlays                                               * Full State RDS (.rds)
                                                                                             * Audit Trail JSON

---

## 2. Prerequisites & Launch Options

### Prerequisites
* R (>= 4.2.0) installed.
* Active internet connection on initial launch to retrieve CRAN and Bioconductor packages.
* Dependencies are sandboxed into a local app_lib/ folder to prevent library permission conflicts.

### Option A: 1-Click Launch (Windows)
Double-click launch_app.bat in the repository root directory. The batch script will automatically:
1. Detect your local 64-bit R installation from PATH, Program Files, or Windows Registry.
2. Set the user library variable to ./app_lib.
3. Run launch.R to check and install missing packages, then open the Shiny app in your default browser.

### Option B: Terminal / RStudio (macOS, Linux, Windows)
From an active R session in the project directory, run:

    source("launch.R")

---

## 3. Tab 1: Data Import & Metadata Setup

### 3.1 Supported File Formats
Select the appropriate format from the Select Source Format dropdown:

| Engine / Platform | Expected Input File | Key Identifier Columns | Quantitative Parsing Logic |
| :--- | :--- | :--- | :--- |
| Spectronaut | *.tsv / *.txt Export Report | PG.ProteinAccessions, PG.Genes, PG.ProteinDescriptions | Extracts PG.Quantity per R.Condition + R.Replicate; converts values <= 0 to NA and applies log2. |
| MSFragger | combined_protein.tsv | Protein.ID (or Protein), Gene, Description | Reads *.Intensity columns (prioritizes standard intensity, falls back to MaxLFQ Intensity if standard is missing). |
| DIA-NN | report.pg_matrix.tsv | Protein.Group, Genes, First.Protein.Description | Parses sample raw data paths, automatically strips long directory paths and file extensions (.raw, .wiff, .d). |
| PEAKS | db.proteins.csv or lfq.proteins.csv | Accession, Description | Extracts gene symbols using regex (GN=([^\\s]+)); strips trailing pipe notations (|) from accessions; extracts Area Sample *. |
| Generic Table | *.tsv / *.csv | Accession, Protein.ID, or Protein | Auto-detects Long vs. Wide format. In Wide format, quant columns must have .Intensity or numeric values. |
| Saved State | *_state.rds | Serialized RDS list | Bypasses raw ingestion and directly restores entire analysis state. |

### 3.2 Sample Metadata Mapping
Once the file is uploaded, the app auto-populates sample mapping inputs:
* Condition: Biological experimental group (e.g., WT, KO, Treated).
* Bio Rep: Biological replicate number (1, 2, 3).

### 3.3 Technical Replicate (TR) Averaging
If multiple injections were performed on the same biological samples:
1. Check the box "Dataset contains Technical Replicates (Injections)".
2. A third input column (Tech Rep) will appear.
3. For each unique biological unit (Condition + BioRep), technical replicates are averaged on the log2 scale:
   LogInt_Consolidated = mean(LogInt_valid)
* Missing Values (NAs): If a protein is detected in 1 out of 2 technical injections, that single valid value is retained. If all technical injections are missing, the biological sample is marked as NA.
4. Click "Process & Map Samples" to register the design.

---

## 4. Tab 2: Normalization & Imputation

### 4.1 Step 1: Reference Protein Normalization (Optional)
Used for spike-ins, biological internal standards, or cell lysis controls (e.g., GAPDH, ACTB, or a specific accession):
1. Check "Normalize to Reference Protein / Loading Control".
2. Select "Gene Symbol" or "Protein Accession", and type the query.
3. For each sample s, the intensity offset is calculated as:
   Offset_s = Ref_s - mean(Ref_global)
   Intensity_corrected_s = Intensity_s - Offset_s

### 4.2 Step 2: Global Normalization (QFeatures)
Aligns whole-proteome sample distributions:
* None: Retains raw log2 measurements.
* Median Centering (center.median): Subtracts each sample median intensity to align medians at zero.
* Mean Centering (center.mean): Subtracts each sample mean intensity to align means at zero.
* Quantile Normalization (quantiles): Enforces identical distribution curves across all runs.
* Variance Stabilizing (vsn): Calibrates affine transformations to stabilize variance across the dynamic range.

### 4.3 Step 3: Missing Value Imputation
* Hybrid (MAR: KNN / MNAR: MinDet):
  - Missing at Random (MAR): Proteins present in >= 50% of replicates within a condition are imputed using global neighbor averaging (nbavg/KNN).
  - Missing Not at Random (MNAR): Proteins present in < 50% of replicates represent dropouts below detection limits and are imputed condition-wise via deterministic minimum substitution (MinDet, positioned at Min_observed - 0.5).
* No Imputation (Show 1-Condition Dropouts on Margins): No synthetic values are generated. Proteins completely missing in one condition but observed in another are retained and projected along plot margins at fixed offsets.
* Standard Methods: Full matrix imputation using MinDet, knn, nbavg, or None.

Click "Apply Transformation" to execute the pipeline and update the interactive QC density plots.

---

## 5. Tab 3: Differential Expression & Visualization

### 5.1 Model Fitting & Contrast Evaluation
Differential expression is computed automatically using limma linear modeling:
* The model matrix is fit without an intercept (~ 0 + Condition).
* Pairwise contrasts (e.g., CondA - CondB) are built and tested using Empirical Bayes moderation (eBayes).
* Invert Contrast Direction: Instantly flips log2 fold-changes (logFC -> -logFC) and swaps experimental/reference sample assignments across plots and tables without refitting the linear model.

### 5.2 Plot Layouts
* Volcano Plot: Plots Log2(Fold Change) vs. -Log10(Adjusted p-value).
* MA Plot: Plots average intensity (AveExpr) vs. Log2(Fold Change). 1-condition dropouts are pinned along visual dashed margins (y = +/- y_cap).
* Scatter Plot: Plots Log2 condition means directly (Ref vs. Exp), with dropouts positioned along the lower baseline.

### 5.3 Collapsible Accordion Views
* Interactive Plotly Explorer (Top): Pan, zoom, and mouse-over points to inspect gene symbols, accessions, and statistics.
* Publication-Ready ggplot (Bottom): Rendered with ggrepel text labeling and publication-quality styling.
* Export individual plots as publication-ready PNG (300 DPI) or vector PDF files with customizable dimensions.

### 5.4 Custom Protein Highlighting
Highlight biological pathways, complexes, or candidate genes:
1. Click "+ Add Protein Set".
2. Choose "List" (comma-, space-, or newline-separated entries) or "Regex Match" (e.g., ^HIST.*, ^RPL.*).
3. Assign custom colors and toggle highlighting and labels independently.

---

## 6. Tab 4: Functional Pathway Enrichment (KEGG & GSEA)

Analyze functional enrichment using clusterProfiler for Mouse (Mus musculus) or Human (Homo sapiens).

### 6.1 Dual Ingestion Modes
* Use Active App Analysis: Runs enrichment on the contrast currently selected in Tab 3 (respecting contrast direction inversion).
* Upload Excel Export (.xlsx): Ingests previously generated Excel workbooks to run enrichment without reloading raw MS runs.

### 6.2 Analysis Methods
* KEGG Over-Representation Analysis (ORA): Evaluates significant proteins against the measured proteome background using Fisher exact test. Filter by direction: All Significant, Upregulated Only, or Downregulated Only.
* Ranked GSEA: Uses a ranked list of all measured proteins sorted by Log2(Fold Change). Generates split-sign dot plots showing activated vs. suppressed pathways.

### 6.3 Pathway Distribution Inspection
* Running Score Curve: View gseaplot2 profiles for any selected pathway.
* Pathway DE Distribution Plot: Overlays all proteins from a specific pathway directly onto the Volcano, MA, or Scatter plot, allowing you to see the individual protein fold changes driving that pathway score.

---

## 7. Tab 5: Export & Reproducibility Audit

### 7.1 Export Deliverables
* Excel Results Workbook (.xlsx): Multi-sheet workbook with one tab per comparison. Column names dynamically incorporate condition labels (e.g., Log2_Treated_Mean, Log2_Control_Mean, Treated_Valid_Count, Control_Valid_Count).
* Reproducible R Pipeline Script (.R): Generates a self-contained R script recreating all file parsing, metadata mappings, normalization transformations, limma designs, and publication figures without needing the Shiny UI.
* Full Session State (.rds): A serialized snapshot containing all raw data, sample maps, normalized matrices, and DE tables. Can be reloaded into Tab 1 under "Re-import Exported RDS/State" to resume analysis instantly.
* JSON Audit Trail (.json): Tracks the application version, session timestamp, R version, loaded package versions, user mapping dictionary, normalization settings, and imputation counts.

---

## 8. Troubleshooting & FAQs

* Error: "No valid intensity columns found"
  - MSFragger: Ensure your table has columns ending with Intensity or MaxLFQ Intensity.
  - Spectronaut: Ensure PG.Quantity is present in your exported report.
  - Generic Table: Select Generic Table (Long or Wide) and verify your intensity columns end with .Intensity or are numeric.

* Why are my fold changes backwards?
  - Contrasts are built as Condition1 - Condition2 based on alphanumeric ordering. To swap the direction without refitting your model, check "Invert Contrast Direction" in Tab 3.

* How do I prevent Git from tracking large package downloads?
  - Verify that app_lib/ is listed in your .gitignore. If it was previously staged, run `git rm -r --cached app_lib` in your terminal.

* Can I run this without the GUI?
  - Yes. Download the standalone script via "Download Reproducible R Script" in Tab 5. It runs headless from any standard R environment.
