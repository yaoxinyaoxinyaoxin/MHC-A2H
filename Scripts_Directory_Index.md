# Project Scripts Directory Index

## 1_Local Analysis of [co]Variant Associations (LAVA)
### 1.1_LAVA
- `1.1.1.1_prepare_GWAS_LAVA_inputs_template.py`
- `1.1.2.1_run_lava_mhc_refined_template.R`
- `1.1.3.1_Manhattan_Plots_GWAS+LAVA_template.R`

## 2_Colocalization-constrained MR analysis
### 2.1_Genome-wide_multi-window_colocalization
- `2.1.1.1_GWAS_preprocessing_for_SuSiE_coloc_template.R`
- `2.1.1.2_eQTL_preprocessing_for_SuSiE_coloc_template.R`
- `2.1.1.4_sc_eQTL_preprocessing_for_SuSiE_template.R`
- `2.1.2.1_OpenGWAS_bulk_eQTL_susie_coloc_blocks_4windows_template.R`
- `2.1.3.1_OneK1K_sceQTL_susie_coloc_blocks_4windows_template.R`
### 2.2_Colocalization_constrained_MR_analysis
- `2.2.1.1_MR_analysis_template.R`
- `2.2.2.1_Integrated_MR_Plots_Analysis_template.R`

## 3_Cross-trait genetic correlation and cell-type-specific pleiotropy extension
### 3.1_Integrate_cell_gene_snp
- `3.1.1.1_forest_plot_script_template.R`
### 3.2_Query_pleiotropy_of_shared_cell_gene_snp
- `3.2.1.1_extract_onek1k_snps_template.R`

## 4_LD-aware signal clustering and multidimensional spatial mapping
### 4.1_LD_aware_clustering_and_3D_scatter_plot
- `4.1.1.1_calculate_LD_clustering_all_template.R`
- `4.1.2.1_calculate_LD_and_filter_template.R`
- `4.1.3.1_batch_single_snp_mr_template.R`
- `4.1.5.1_visualize_3D_signals_template.R`
- `4.1.4.1_plot_3d_scatter_with_PCA_template.R`

## 5_Cross-trait functional and network analyses
### 5.1_Cross-trait_GO_and_KEGG_enrichment
- `5.1.1.1_Enrichment_Analysis_Comparison_template.R`
### 5.2_Cross-trait_PPI_network_analysis
- `5.2.1.1_visualize_integrated_PPI_template.R`
### 5.3_Cross-trait_NicheNet_prediction
- `5.3.1.1_run_nichenet_analysis_template.R`
- `5.3.2.1_NicheNet_Combined_Visualization_template.R`

## 6_Cis-MR based on UKB pQTL
### 6.1_Cis-MR_analysis_based_on_UKB_pQTL
- `6.1.1.1_get_gene_info_GRCh37_template.R`
- `6.1.1.2_extract_snps_by_gene_region_template.R`
- `6.1.1.3_IV_Filtering_250kb_template.R`
- `6.1.1.4_Local_LD_Clumping_template.R`
- `6.1.1.5_Fallback_LD_Clumping_template.R`
- `6.1.2.1_cisMR_Analysis_template.R`
- `6.1.3.1_MR_Visualization_Compare_template.R`

## 7_PheWAS and downstream analyses
### 7.1_PheWAS_based_on_rs1800628
- `7.1.1.1_extract_target_snp_template.R`
- `7.1.2.1_extract_target_snp_template.R`
- `7.1.3.1_extract_target_snp_template.R`
- `7.1.4.1_extract_target_snp_template.R`
### 7.2_PPI_analysis_of_candidate_plasma_proteins_P_5e-8
- `7.2.1.1_visualize_PPI_template.R`
### 7.3_NicheNet_analysis_of_candidate_plasma_proteins_P_5e-8
- `7.3.1.1_nichenet_plasma_proteins_template.R`
### 7.4_PPI_analysis_of_candidate_plasma_proteins_P_1e-11
- `7.4.1.1_visualize_PPI_template.R`
### 7.5_NicheNet_analysis_of_candidate_plasma_proteins_P_1e-11
- `7.5.1.1_nichenet_plasma_proteins_template.R`
### 7.6_Cis-MR_of_plasma_proteins_and_immune_cells
- `7.6.1.1_get_gene_info_GRCh37_template.R`
- `7.6.1.2_extract_snps_by_gene_region_template.R`
- `7.6.1.3_IV_Filtering_250kb_template.R`
- `7.6.1.4_Local_LD_Clumping_template.R`
- `7.6.1.5_Fallback_LD_Clumping_template.R`
- `7.6.2.1_cisMR_Analysis_template.R`
### 7.7_Colocalization_of_plasma_proteins_and_immune_cells
- `7.7.1.1_GWAS_preprocessing_template.R`
- `7.7.1.2_pQTL_preprocessing_template.R`
- `7.7.2.1_susie_coloc_pairwise_template.R`
### 7.8_Cross-trait_cross-cohort_colocalization_and_Manhattan_visual_alignment
- `7.8.1.1_UKB_GWAS_preprocessing_template.R`
- `7.8.1.2_FinnGen_GWAS_preprocessing_template.R`
- `7.8.2.1_susie_coloc_template.R`
- `7.8.3.1_plot_manhattan_template.R`

## 8_Sensitivity and robustness analyses
### 8.1_GEO_data_validation_GSE242252
- `8.1.1.1_process_GEO_data_template.R`
- `8.1.2.1_split_expression_matrix_template.R`
- `8.1.3.1_verify_target_genes_template.R`

## 9_Genome_build_conversion
### 9.1_Genome_build_conversion
- `9.1.1.1_liftOver_template.R`

