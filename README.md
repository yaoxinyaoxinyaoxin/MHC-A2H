# Code Availability for: "A shared MHC immunogenetic signal connects aging, rheumatoid arthritis, and herpes zoster through chronic high inflammatory burden–compensatory immune tolerance dysregulation"
This repository contains the complete, standardized analytical pipeline utilized to generate the results presented in this study. By integrating large-scale genome-wide association studies (GWAS) of a multivariate aging latent factor (mvAge), rheumatoid arthritis (RA), and herpes zoster (HZ) with multi-omics quantitative trait loci (including bulk eQTLs, single-cell eQTLs, and plasma pQTLs), this framework dissects the shared genetic and molecular architecture linking aging, autoimmunity, and viral reactivation.

## Overview of Methodology

The analytical pipeline is designed to overcome the extreme linkage disequilibrium (LD) and complex haplotype structures inherent to the major histocompatibility complex (MHC). It sequentially translates macroscopic phenotypic associations into specific molecular and cellular regulatory networks, providing a mechanistic foundation for the proposed life-course "chronic high inflammatory burden–compensatory immune tolerance dysregulation" model:

1. **LAVA (Local Analysis of [co]Variant Associations):** Maps local genetic correlations to demonstrate non-random regional aggregation of shared genetic effects across the extended MHC region among mvAge, RA, and HZ.
2. **SuSiE & Colocalization:** Employs an LD-aware colocalization framework, integrating Sum of Single Effects (SuSiE) fine-mapping across dynamic genomic windows, to resolve complex LD and map shared pleiotropic signals.
3. **Colocalization-constrained Mendelian Randomization (MR):** Employs eQTLs (bulk and single-cell) as instrumental variables to infer causal directions towards the three primary GWAS traits (mvAge, RA, and HZ). To mitigate horizontal pleiotropy and LD confounding, MR estimations are strictly constrained to single-instrument analyses anchored by strongly colocalized variants (e.g., PP.H4 > 0.8).
4. **LD-aware Signal Clustering & 3D Spatial Mapping:** Translates complex genetic pleiotropy into a unified principal component axis via LD-aware clustering (r² > 0.6). This maps the independent shared genetic signals into a three-dimensional spatial coordinate system, polarizing them into distinct functional trajectories (e.g., concordant protection representing healthier aging and reduced disease risk, versus shared pathogenesis).
5. **Cross-trait Functional and Network Analyses:** Reconstructs cell-type-specific intercellular communication (NicheNet) and protein-protein interaction (PPI) union networks across the three primary traits (mvAge, RA, and HZ) to evaluate shared downstream molecular consequences and key immune communication hubs.
6. **Cis-MR based on UKB pQTL:** Infers the causal effects of plasma proteins on the risk of RA and HZ using cis-acting instrumental variables (±250 kb). This includes cross-trait intersection analyses to identify whether their causal effect directions are concordant (driving shared pathogenic risks) or antagonistic (exerting protective trade-offs).
7. **PheWAS and Downstream Analyses:** Deep profiles the core shared signal tagged by rs1800628, systematically delineating the causal hierarchy of immune alterations (e.g., elevated pro-inflammatory mediators like IL15 and CXCL9, alongside compensatory T-cell tolerance markers like TIGIT and PDCD1) across multiple independent cohorts.
8. **Sensitivity and Robustness Analyses:** Validates the identified genetic associations using independent transcriptomic cohorts (e.g., GSE242252 RNA-seq) to confirm that the GWAS risk genes maintain stable baseline expression and to exclude reverse causality from acute infection.

## Repository Structure

The directory structure strictly mirrors the sequential workflow detailed in the Methods section:

* **`1_Local Analysis of [co]Variant Associations (LAVA)`**
  Scripts for evaluating local genetic correlations across the MHC region.
* **`2_Colocalization-constrained MR analysis`**
  Integration of SuSiE fine-mapping, multi-window colocalization, and robust single-instrument MR using OpenGWAS bulk eQTLs and OneK1K sc-eQTLs against the three primary traits.
* **`3_Cross-trait genetic correlation and cell-type-specific pleiotropy extension`**
  Extraction of cell-gene-SNP triplets and evaluation of macroscopic consistency of genetic effects.
* **`4_LD-aware signal clustering and multidimensional spatial mapping`**
  LD calculations (r² > 0.6), independent signal clustering, and 3D PCA spatial mapping to polarize shared genetic etiology.
* **`5_Cross-trait functional and network analyses`**
  Cross-trait GO/KEGG functional enrichment, PPI union network integration, and NicheNet intercellular communication prediction across mvAge, RA, and HZ.
* **`6_Cis-MR based on UKB pQTL`**
  Cis-MR analysis (±250 kb) inferring the causal effects of plasma proteins on RA and HZ to identify concordant or antagonistic trade-offs.
* **`7_PheWAS and downstream analyses`**
  Deep profiling of the shared signal tagged by rs1800628, including systemic immune cell MR, cross-cohort validation, and reconstruction of the causal hierarchy of immune alterations (e.g., IL15, CXCL9, TIGIT, PDCD1).
* **`8_Sensitivity and robustness analyses`**
  Scripts for external validation using independent transcriptomic cohorts (e.g., GSE242252 RNA-seq) to evaluate the robustness of the identified genetic associations against infection-induced transcriptional perturbations.
* **`9_Genome_build_conversion`**
  Utility scripts for harmonizing genetic coordinates (hg38 to hg19) via `liftOver`.

## System Requirements

### R Environment
The analytical pipeline was executed in **R version 4.5**. The following core R packages must be installed:
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
