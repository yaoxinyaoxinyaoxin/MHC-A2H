# 项目代码库说明 (Bilingual README)

*For the English version, please refer to [README.md](./README.md).*

本项目包含了用于生成论文 **《A shared MHC immunogenetic signal links aging, rheumatoid arthritis, and herpes zoster through high inflammatory burden–compensatory immune tolerance dysregulation》** 结果的完整标准化分析流水线。

通过整合大规模全基因组关联研究 (GWAS)、单细胞表达数量性状基因座 (sc-eQTL) 以及高分辨率血浆蛋白 QTL (pQTL) 数据，本研究系统性地揭示了衰老、类风湿关节炎（RA）和带状疱疹（HZ）之间共享的遗传与分子架构。

## 核心方法学概述

由于主要组织相容性复合体（MHC）区域存在极端的连锁不平衡（LD）和复杂的单倍型结构，本分析流水线被设计用于克服这些统计学障碍，将宏观的表型关联逐步转化为微观的分子细胞调控网络：

1. **LAVA (局部遗传相关性分析):** 绘制局部遗传相关性图谱，证明共享遗传效应在 MHC 扩展区域内的非随机区域性聚集。
2. **SuSiE 与共定位分析:** 采用考虑 LD 结构的共定位框架，结合多窗口动态 SuSiE (Sum of Single Effects) 精细映射，解析复杂的 LD 并定位共享的多效性因果信号。
3. **共定位约束的孟德尔随机化 (MR):** 利用多组学数据（bulk eQTL、sc-eQTL 和血浆 pQTL）作为工具变量，推断遗传信号、免疫细胞表型与疾病结局之间的因果方向，同时严格控制水平多效性偏倚。
4. **LD 感知信号聚类与 3D 空间映射:** 将复杂的多效性遗传信号降维并映射到统一的主成分轴（纵向时间轴）上，将信号极化为具有明确方向性的功能轨迹。
5. **NicheNet 与功能网络分析:** 重建细胞类型特异性的细胞间通讯（配体-靶基因调控网络）及蛋白质-蛋白质相互作用（PPI）网络。这进一步验证了由遗传变异驱动的下游分子效应，例如促炎介质（IL15, CXCL9）的升高以及代偿性免疫耐受（TIGIT, PDCD1）机制的激活。

## 目录结构

本仓库的目录结构与手稿中的 Methods（方法学）部分严格对应，反映了分析的先后顺序：

* **`1_Local Analysis of [co]Variant Associations (LAVA)`**
  MHC 区域局部遗传相关性分析脚本。
* **`2_Colocalization-constrained MR analysis`**
  整合 OpenGWAS bulk eQTL 与 OneK1K sc-eQTL，进行 SuSiE 精细映射、多窗口共定位与稳健的 MR 分析。
* **`3_Cross-trait genetic correlation and cell-type-specific pleiotropy extension`**
  提取共享的 cell-gene-SNP 三元组，并评估跨表型遗传效应的宏观一致性。
* **`4_LD-aware signal clustering and multidimensional spatial mapping`**
  LD 计算、独立信号聚类及 3D PCA 散点图可视化（绘制纵向时间轴）。
* **`5_Cross-trait functional and network analyses`**
  跨性状 GO/KEGG 富集分析、PPI 网络整合及 NicheNet 细胞间通讯预测。
* **`6_Cis-MR based on UKB pQTL`**
  基于 UKB pQTL 的顺式孟德尔随机化分析，筛选潜在药物靶点（如 TIGIT）。
* **`7_PheWAS and downstream analyses`**
  对 rs1800628 标记的共享信号进行全表型关联深度剖析，包含系统性免疫细胞 MR、跨队列验证、多效性扩展以及 cis-MR 因果网络的分层边缘捆绑 (Hierarchical Edge Bundling) 可视化。
* **`8_Sensitivity and robustness analyses`**
  使用独立的 GEO 转录组队列（如 GSE242252）进行验证，以排除急性感染引发的反向因果干扰。
* **`9_Genome_build_conversion`**
  用于将遗传坐标（hg38 到 hg19）进行转换的 `liftOver` 实用脚本。

## 运行环境与依赖

### R 语言环境
分析流水线在 **R 4.x** 版本下执行。请确保已安装以下核心 R 包：
* **统计遗传学与 MR:** `TwoSampleMR`, `LAVA`, `susieR`, `coloc`, `MRPRESSO`
* **单细胞与网络分析:** `nichenetr`, `Seurat`, `igraph`, `ggraph`, `clusterProfiler`
* **数据处理与可视化:** `tidyverse`, `data.table`, `ggplot2`, `plotly`

### 命令行工具
* **PLINK (v1.9 / v2.0):** 脚本中涉及的 LD 矩阵生成及 LD clumping 强依赖于 `plink`。请确保该可执行文件已加入系统 PATH，或通过命令行参数明确指定绝对路径。

## 使用说明

为了保证代码的环境可移植性与跨集群复现性，本仓库中的所有 R 脚本均已利用 `optparse` 包重构为参数化的命令行工具，**清除了所有硬编码的绝对路径**。

**执行示例:**
```bash
Rscript 5.3.1.1_run_nichenet_analysis_template.R \
  --input_mr ./data/mr_results.csv \
  --cell_map ./data/cell_type_mapping.csv \
  --bg_genes ./data/background_genes.txt \
  --nichenet_db ./db/nichenet_database \
  --out_dir ./results/nichenet_output
```
*提示：执行 `Rscript <脚本名称.R> --help` 即可查看该脚本支持的所有必需及可选参数。*


## 许可证
本项目采用 MIT License 开源许可证。详情请参阅 `LICENSE` 文件。

---
### 🌐 Interactive 3D Visualizations (交互式 3D 可视化)
Interactive 3D causal networks and scatter plots related to this study are hosted on GitHub Pages:
本研究相关的交互式 3D 因果网络图与散点图已托管至 GitHub Pages，点击下方链接即可在线探索：
[Explore the 3D Interactive Plots (探索 3D 交互图)](https://yaoxinyaoxinyaoxin.github.io/MHC-A2H/)

---
### 🌐 Interactive 3D Visualizations (交互式 3D 可视化)
Interactive 3D causal networks and scatter plots related to this study are hosted on GitHub Pages:
本研究相关的交互式 3D 因果网络图与散点图已托管至 GitHub Pages，点击下方链接即可在线探索：
[Explore the 3D Interactive Plots (探索 3D 交互图)](https://yaoxinyaoxinyaoxin.github.io/MHC-A2H/)
