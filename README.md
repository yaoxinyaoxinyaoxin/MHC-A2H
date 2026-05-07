# Code Availability for: "A shared MHC immunogenetic signal links aging, rheumatoid arthritis, and herpes zoster through high inflammatory burden–compensatory immune tolerance dysregulation"

This repository contains the complete, standardized analytical pipeline utilized to generate the results presented in the manuscript. By integrating large-scale genome-wide association studies (GWAS) with cell-resolved single-cell expression quantitative trait loci (sc-eQTLs) and high-resolution plasma protein QTLs (pQTLs), this framework dissects the shared genetic and molecular architecture linking aging, rheumatoid arthritis (RA), and herpes zoster (HZ).

## Overview of Methodology

The analytical pipeline is designed to overcome the extreme linkage disequilibrium (LD) and complex haplotype structures inherent to the major histocompatibility complex (MHC). It sequentially translates macroscopic phenotypic associations into specific molecular and cellular regulatory networks:

1. **LAVA (Local Analysis of [co]Variant Associations):** Maps local genetic correlations to demonstrate non-random regional aggregation of shared genetic effects across the extended MHC region.
2. **SuSiE & Colocalization:** Employs an LD-aware colocalization framework, integrating Sum of Single Effects (SuSiE) fine-mapping across dynamic genomic windows, to resolve complex LD and map shared pleiotropic signals.
3. **Colocalization-constrained Mendelian Randomization (MR):** Uses multi-omics data (bulk eQTLs, sc-eQTLs, and plasma pQTLs) as instrumental variables to infer causal directions between genetic signals, immune cell phenotypes, and disease outcomes while mitigating horizontal pleiotropy.
4. **LD-aware Signal Clustering & 3D Spatial Mapping:** Translates complex genetic pleiotropy into a unified principal component axis (Longitudinal Time Axis), polarizing signals into distinct functional trajectories.
5. **NicheNet & Functional Networks:** Reconstructs cell-type-specific intercellular communication (ligand-target regulatory networks) and protein-protein interaction (PPI) networks, validating downstream molecular consequences such as elevated pro-inflammatory mediators (IL15, CXCL9) and compensatory immune tolerance (TIGIT, PDCD1).

## Repository Structure

The directory structure strictly mirrors the sequential workflow detailed in the Methods section of the manuscript:

* **`1_Local Analysis of [co]Variant Associations (LAVA)`**
  Scripts for evaluating local genetic correlations across the MHC region.
* **`2_Colocalization-constrained MR analysis`**
  Integration of SuSiE fine-mapping, multi-window colocalization, and robust MR using OpenGWAS bulk eQTLs and OneK1K sc-eQTLs.
* **`3_Cross-trait genetic correlation and cell-type-specific pleiotropy extension`**
  Extraction of cell-gene-SNP triplets and evaluation of macroscopic consistency of genetic effects.
* **`4_LD-aware signal clustering and multidimensional spatial mapping`**
  LD calculations, independent signal clustering, and 3D PCA scatter plot visualizations.
* **`5_Cross-trait functional and network analyses`**
  Cross-trait GO/KEGG functional enrichment, PPI network integration, and NicheNet intercellular communication prediction.
* **`6_Cis-MR based on UKB pQTL`**
  Cis-MR analysis prioritizing plasma proteins (e.g., TIGIT) associated with RA and HZ.
* **`7_PheWAS and downstream analyses`**
  Deep profiling of the shared signal tagged by rs1800628, including systemic immune cell MR, cross-cohort validation, extended pleiotropy mapping, and hierarchical edge-bundling visualization of cis-MR networks.
* **`8_Sensitivity and robustness analyses`**
  Validation using independent GEO transcriptomic cohorts (e.g., GSE242252) to exclude reverse causality from acute infection.
* **`9_Genome_build_conversion`**
  Utility scripts for harmonizing genetic coordinates (hg38 to hg19) via `liftOver`.

## System Requirements

### R Environment
The analytical pipeline was executed in **R version 4.x**. The following core R packages must be installed:
* **Statistical Genetics & MR:** `TwoSampleMR`, `LAVA`, `susieR`, `coloc`, `MRPRESSO`
* **Single-Cell & Network Analysis:** `nichenetr`, `Seurat`, `igraph`, `ggraph`, `clusterProfiler`
* **Data Manipulation & Visualization:** `tidyverse`, `data.table`, `ggplot2`, `plotly`

### Command-line Tools
* **PLINK (v1.9 / v2.0):** Required for LD matrix generation and LD clumping. Ensure the `plink` executable is in your system PATH or correctly specified via command-line arguments.

## Usage

To ensure environmental portability and reproducibility across different computational clusters, all R scripts have been refactored into parameterized command-line tools utilizing the `optparse` package. All hardcoded absolute paths have been removed.

**Execution Example:**
```bash
Rscript 5.3.1.1_run_nichenet_analysis_template.R \
  --input_mr ./data/mr_results.csv \
  --cell_map ./data/cell_type_mapping.csv \
  --bg_genes ./data/background_genes.txt \
  --nichenet_db ./db/nichenet_database \
  --out_dir ./results/nichenet_output
```
*Note: Execute `Rscript <script_name.R> --help` to review all required and optional parameters for any given script.*


## License
This code is distributed under the MIT License. See `LICENSE` for more information.

---
### 🌐 Interactive 3D Visualizations
Interactive 3D causal networks and scatter plots related to this study are hosted on GitHub Pages:
[Explore the 3D Interactive Plots](https://yaoxinyaoxinyaoxin.github.io/MHC-A2H/)

---
### 🌐 Interactive 3D Visualizations
Interactive 3D causal networks and scatter plots related to this study are hosted on GitHub Pages:
[Explore the 3D Interactive Plots](https://yaoxinyaoxinyaoxin.github.io/MHC-A2H/)
