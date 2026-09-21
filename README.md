# Elk-1-Potential-Transcriptional-Regulator-in-Neuroblastoma
R pipeline for investigating ELK1 and PEA3 family transcription factors as regulators of cell surface markers CD44, CD133, and CD114 in neuroblastoma using transcriptomic and survival data.

Two R scripts that analyse the expression of the ETS transcription factors **ELK1** and the **PEA3 subfamily (ETV1, ETV4, ETV5)** in neuroblastoma, using public RNA-seq data from cell lines and patients. It also relates their expression to three focus markers (**CSF3R, CD44, PROM1**), to surface/CD-antigen marker lists, and to MYCN status.

## Scripts

| Script | Survival endpoint |
|---|---|
| `neuroblastoma_efs.R` | Event-free survival (EFS) |
| `neuroblastoma_os.R` | Overall survival (OS) |

The scripts differ in the survival endpoint, which is used for the optimal cutpoint, the Kaplan-Meier plots and the time-dependent ROC, and appears in the output file names (`EFS` / `OS`). They also differ in one further respect: in the OS script the multi-page violin PDFs for the CD antigen / surface / neuroblastoma marker lists (Part B, files `25`-`28`) are disabled; uncomment the `save_violin_pdf(...)` lines to generate them. Both scripts run top to bottom and are organised in three parts.

## What the scripts do

**Part A - Cell lines (GSE89413), ELK1**
- ELK1 expression across neuroblastoma cell lines (bar plot; SKNBE2, SHSY5Y and KELLY highlighted)
- ELK1 together with CSF3R / CD44 / PROM1 across cell lines
- Compact heatmap of the 500 most variable surface markers

**Part B - Patients (GSE62564), ELK1**
- Downloads clinical data and the RNA-seq matrix from GEO, maps RefSeq IDs to gene symbols, and matches clinical records to expression samples
- Splits patients into ELK1-High / ELK1-Low (optimal cutpoint, median or tertile)
- Kaplan-Meier survival curves and time-dependent ROC (12 / 36 / 60 months) for the chosen endpoint (EFS or OS)
- ELK1 correlation with CSF3R / CD44 / PROM1 and with a wider marker panel
- Violin plots (Wilcoxon and ANOVA) for the focus markers and, in the EFS script, for the CD antigen, surface and neuroblastoma marker lists
- MYCN expression across cell lines, colored by MYCN amplification status

**Part C - PEA3 family (ETV1, ETV4, ETV5)**
- Each gene is analysed **independently** (no averaged score): the full pipeline is run once per gene
- Cell line bar plots, surface-marker heatmap, and MYCN bar plot
- Per-gene patient classification, Kaplan-Meier survival, time-dependent ROC, correlation with focus markers, and violin plots (one PNG per marker)

## Requirements

- R (a recent 4.x release is recommended)
- An internet connection (packages and GEO data are downloaded on the first run)

All packages are installed (if missing) and loaded in a single block at the top of the script:

| Source | Packages |
|---|---|
| CRAN | ggplot2, dplyr, tidyr, tibble, stringr, readxl, ggrepel, pheatmap, ggpubr, writexl, igraph, scales, survival, survminer, timeROC, cmprsk |
| Bioconductor | GEOquery, limma, org.Hs.eg.db, Biobase |

## Input data

Put the following files in the data directories you configure (see [Setup](#setup)).

| File | Description |
|---|---|
| `GSE89413_2016-10-30-NBL-cell-line-STAR-fpkm.txt` | Neuroblastoma cell line FPKM matrix (GEO: GSE89413) |
| `2016-11-17-CellLine-MYCN-status.txt` | Tab-delimited table with two columns: cell line name, MYCN status (`Amplified` / `Nonamplified`) |
| `cd_antigen_markers.csv` | CD antigen marker list; must contain a `Gene Names` column (whitespace-separated gene symbols) |
| `surface_marker_candidates.csv` | Surface marker candidates; must contain a `Gene Names` column (whitespace-separated gene symbols) |
| `neuroblastoma_markers.csv` | Neuroblastoma markers; must contain `Marker` and `Alias` columns |

The patient dataset **GSE62564** (RNA-seq, ~498 samples) is downloaded automatically from GEO with `GEOquery`.

## Setup

Edit the placeholder paths in the script before running (both scripts use the same variables):

| Location in script | Variable(s) | Placeholder | Meaning |
|---|---|---|---|
| Part A `config` | `DATA_DIR`, `OUT_DIR` | `path/to/your/data`, `path/to/your/output` | Input files / output root |
| Part B `config` | `DATA_DIR`, `OUT_DIR` | `path/to/your/data`, `path/to/your/output` | Input files / output root |
| MYCN supplementary block | `DATA_DIR_CELL` | `path/to/your/cell_line_data` | Directory containing the MYCN status file |
| PEA3 `SETTINGS` | `DATA_DIR_CELL`, `DATA_DIR_PAT` | `path/to/your/pea3_data` | Input files for Part C |
| PEA3 `SETTINGS` | `OUT_DIR` | `path/to/your/pea3_output` | Output root for Part C |

All of these can point to the same folder if you keep every input file together.

## Usage

Run the whole script in order, for example:

```r
source("neuroblastoma_ets_analysis_efs.R")   # or neuroblastoma_ets_analysis_os.R
```

or from a terminal:

```bash
Rscript neuroblastoma_ets_analysis_efs.R
```

Packages are loaded once at the top, and Part C relies on them, so run the file from the top rather than executing Part C on its own. If you run both scripts, use different `OUT_DIR` / `pea3_output` folders (or move the results in between), otherwise the ELK1 files with identical names (e.g. `07_TimeROC`, `01_ELK1_distribution`) will be overwritten.

Key settings (defined in the config/settings block of the relevant part):

- `FOCUS_CELLS` - cell lines highlighted in the bar plots (`SKNBE2`, `SHSY5Y`, `KELLY`)
- `FOCUS3` - focus markers (`CSF3R`, `CD44`, `PROM1`)
- `ELK1_SPLIT` - patient grouping method: `"optimal"` (survminer cutpoint on EFS, minimum group proportion 0.2), `"median"` or `"tertile"`. In Part C the same setting is applied to each PEA3 gene.
- `PEA3_GENES` - genes analysed separately in Part C
- `neuro_genes`, `stem_sigs` - gene sets used to group genes (Neuroblastoma / Stemness / Neuro+Stem) and to build the correlation panels; edit them to add or remove genes from the literature

## Output

Every part writes to a `Figures/` folder (PDF and/or PNG) and a `Results/` folder (Excel tables).

| Part | Output location |
|---|---|
| A | `<OUT_DIR>/GSE89413_outputs/{Figures,Results}` |
| B (including the MYCN supplementary plot) | `<OUT_DIR>/GSE62564_outputs/{Figures,Results}`, plus `geo_cache/` and `GSE62564_Supp/` for downloaded GEO files |
| C | `<pea3_output>/{Figures,Results}`, plus `geo_cache/` and `GSE62564_Supp/` |

Typical files:

- Cell line expression bar plots, e.g. `02_ELK1_across_celllines`, `01_<GENE>_across_celllines`
- Surface marker heatmaps, e.g. `05_Surface_heatmap`, `02_Surface_heatmap`
- MYCN bar plots, e.g. `30_MYCN_status`, `04_MYCN_status`
- Patient groups: `01_ELK1_distribution`, `06_<GENE>_distribution`, `*_optimal_cutpoint`, `*_patient_classification.xlsx`
- Survival (EFS script): `03_KM_EFS_ELK1`, `07_KM_EFS_<GENE>`, `07_TimeROC`, `08_TimeROC_EFS_<GENE>`
- Survival (OS script): `03_KM_OS_ELK1`, `07_KM_OS_<GENE>`, `07_TimeROC`, `08_TimeROC_OS_<GENE>`
- Correlation: `*_focus_scatter`, `*_GeneCorr*`, `*_correlation_panel.xlsx`
- Violin plots: `11_<GENE>_Focus3_violin_<MARKER>.png` and `*_stats.xlsx` tables

`<GENE>` is one of `ETV1`, `ETV4`, `ETV5`.

## Data sources

- Cell lines: GEO accession [GSE89413](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89413)
- Patients: GEO accession [GSE62564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62564)

