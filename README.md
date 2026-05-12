# MultiOmicsBridge

<!-- badges -->
[![Bioc-check](https://github.com/SubhadipJana1409/MultiOmicsBridge/actions/workflows/bioc-check.yml/badge.svg)](https://github.com/SubhadipJana1409/MultiOmicsBridge/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Bioconductor devel](https://img.shields.io/badge/Bioconductor-devel-blue)](https://bioconductor.org/packages/devel/bioc/html/MultiOmicsBridge.html)

## Overview

**MultiOmicsBridge** is an R/Bioconductor package providing an end-to-end,
reproducible computational framework for integrative analysis of paired host
transcriptomics (bulk RNA-seq) and gut microbiome (16S rRNA or shotgun
metagenomics) data.

No single Bioconductor package currently provides all four critical
capabilities for host-microbiome integration in a unified, microbiome-aware
workflow. MultiOmicsBridge fills this gap:

| Gap | What MultiOmicsBridge provides |
|---|---|
| Microbiome-aware normalization | CLR transformation (compositionally correct) alongside TMM/voom for RNA-seq |
| Joint dimensionality reduction | Sparse multi-block PLS-DA (DIABLO) across both data layers |
| Cross-omics biomarker selection | Loading scores + Spearman correlation network between host genes and microbial taxa |
| Multi-omics diagnostic value | Host-only vs microbiome-only vs joint classifier comparison with nested CV |

## The Five Modules

```
Module 1: Data Harmonization
  loadHostData()          TMM + voom normalization → SummarizedExperiment
  loadMicrobiomeData()    CLR/TSS transformation  → SummarizedExperiment
  matchSamples()          Paired sample matching  → MultiAssayExperiment

Module 2: Joint Dimensionality Reduction
  jointDimReduction()     DIABLO sparse multi-block PLS-DA

Module 3: Biomarker Discovery
  biomarkerDiscovery()    Loading scores + cross-omics correlation network

Module 4: Diagnostic Classification
  diagnosticClassifier()  Host-only vs MB-only vs joint Random Forest (CV)

Module 5: Visualization & Reporting
  plotIntegration()           Joint biplot with loading arrows
  plotBiomarkerNetwork()      Cross-omics correlation heatmap
  plotClassifierComparison()  ROC curves / AUC bar chart
  plotSankey()                Feature flow diagram
  generateReport()            Structured text summary
```

## Installation

```r
# From Bioconductor (once accepted)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("MultiOmicsBridge")

# Development version from GitHub
BiocManager::install("SubhadipJana1409/MultiOmicsBridge")
```

## Quick Start

```r
library(MultiOmicsBridge)

# Step 1: Load and normalize each data type
host_se <- loadHostData(host_count_matrix, col_data = sample_metadata)
mb_se   <- loadMicrobiomeData(taxa_count_matrix, normalization = "CLR")

# Step 2: Match paired samples
mae <- matchSamples(host_se, mb_se)

# Step 3: Run the full pipeline in one call
result <- MultiOmicsBridgeAnalysis(
    mae,
    outcome         = sample_metadata$condition,
    n_components    = 2,
    n_features_host = 50,
    n_features_mb   = 20,
    n_biomarkers    = 50,
    cv_folds        = 5
)

result
# MOBResult
#   Integration    : DIABLO
#   Samples        : 106
#   Components     : 2
#   Biomarkers     : 70 (50 host, 20 microbiome)
#   -- Classifier AUC (mean ± SD) -----------
#   host_only:       0.847 ± 0.042
#   microbiome_only: 0.793 ± 0.055
#   joint:           0.931 ± 0.029

# Access results
biomarkers(result)          # ranked biomarker DataFrame
performance(result)         # classifier AUC list
integrationScores(result)   # DIABLO sample score matrix
featureLoadings(result)     # host + microbiome loading matrices

# Visualize
plotIntegration(result, outcome = outcome)
plotBiomarkerNetwork(result, mae)
plotClassifierComparison(result, type = "bar")
plotSankey(result, n_features = 8)
generateReport(result)
```

## African Health Context

MultiOmicsBridge is designed with African health priorities in mind and
validated on three disease contexts:

| Disease | Dataset | Host | Microbiome |
|---|---|---|---|
| IBD | HMP_2019_ibdmdb (IBDMDB) | Simulated (IBD signatures) | Real (curatedMetagenomicData) |
| Tuberculosis | GEO: GSE79362 | Real (blood RNA-seq) | Simulated (TB dysbiosis) |
| HIV/ART | Published signatures | Simulated | Simulated |

Africa carries ~25% of the global TB burden, and gut-lung axis interactions
are increasingly recognized as relevant to TB progression and treatment
response. For HIV/ART, gut microbiome composition is a significant predictor
of immune reconstitution. MultiOmicsBridge provides a standardized,
accessible tool for these analyses in resource-limited research settings.

## Why CLR for microbiome data?

Microbiome count data is **compositional** — only relative abundances are
observed. Standard correlation and distance measures applied to compositional
data produce spurious results (the Aitchison problem). The centered log-ratio
(CLR) transformation maps compositional data to real Euclidean space:

```
clr(x_j) = log(x_j + δ) − mean_k[log(x_k + δ)]
```

where δ is a small pseudocount. This removes the unit-sum constraint and
enables valid correlation analysis between host genes and microbial taxa.

## Key functions

| Function | Input | Output |
|---|---|---|
| `loadHostData()` | count matrix | SummarizedExperiment |
| `loadMicrobiomeData()` | taxa table | SummarizedExperiment |
| `matchSamples()` | 2 × SE | MultiAssayExperiment |
| `jointDimReduction()` | MAE + outcome | list (scores, loadings) |
| `biomarkerDiscovery()` | MAE + DR result | DataFrame |
| `diagnosticClassifier()` | MAE + outcome | list (AUC per model) |
| `MultiOmicsBridgeAnalysis()` | MAE + outcome | MOBResult |

## Validation scripts

```r
# Real IBD microbiome data (curatedMetagenomicData, ~10-15 min)
source(system.file("scripts", "run_ibdmdb_demo.R",
                    package = "MultiOmicsBridge"))

# African disease contexts (TB + HIV, ~20-30 min)
source(system.file("scripts", "run_african_context.R",
                    package = "MultiOmicsBridge"))

# Performance benchmark
source(system.file("scripts", "run_benchmark.R",
                    package = "MultiOmicsBridge"))
```

## MOBResult S4 class

```r
result                     # print compact summary
integrationScores(result)  # matrix: samples × components
featureLoadings(result)    # list: $host and $microbiome loading matrices
biomarkers(result)         # DataFrame: feature, layer, loading_score, cross-cor
performance(result)        # list: host_only, microbiome_only, joint AUC
```

## References

- Rohart F et al. (2017). mixOmics. *PLoS Comput Biol* 13(11):e1005752.
- Franzosa EA et al. (2019). *Nature Microbiology* 4:293–305.
- Wright MN & Ziegler A (2017). ranger. *J Stat Softw* 77(1):1–17.
- Aitchison J (1982). The statistical analysis of compositional data. *JRSS-B* 44(2):139–177.

## Citation

> Jana S (2026). MultiOmicsBridge: Integrative Multi-Omics Analysis of
> Host Transcriptomics and Gut Microbiome Data. R package version 0.99.0.
> https://github.com/SubhadipJana1409/MultiOmicsBridge

## License

MIT © Subhadip Jana
