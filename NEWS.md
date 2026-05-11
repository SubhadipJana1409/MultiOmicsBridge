# MultiOmicsBridge 0.99.0

## New features

* Initial release submitted to Bioconductor.

* `loadHostData()`: Import and normalize bulk RNA-seq count data using
  TMM normalization and limma-voom precision weights. Returns a
  `SummarizedExperiment` with both raw counts and voom-transformed values.

* `loadMicrobiomeData()`: Import and normalize microbiome taxa tables using
  centered log-ratio (CLR) or total sum scaling (TSS) transformation,
  with automatic zero-inflation handling via pseudocounts. Returns a
  `SummarizedExperiment`.

* `matchSamples()`: Identify and retain paired samples present in both
  host and microbiome datasets. Returns a `MultiAssayExperiment`.

* `jointDimReduction()`: Run sparse multi-block PLS-DA (DIABLO from
  mixOmics) for joint dimensionality reduction across host and microbiome
  data layers. Returns integrated sample scores and feature loadings.

* `biomarkerDiscovery()`: Identify top host genes and microbial taxa using
  DIABLO sparse loadings and cross-omics Spearman correlation networks.
  Returns a ranked `DataFrame` of multi-omics biomarkers.

* `diagnosticClassifier()`: Train and evaluate host-only, microbiome-only,
  and joint Random Forest classifiers using nested cross-validation.
  Provides AUC-ROC comparison quantifying the multi-omics advantage.

* `MultiOmicsBridgeAnalysis()`: One-call wrapper executing the full
  pipeline from matched `MultiAssayExperiment` through classification.

* `plotIntegration()`: Joint biplot of sample scores in integrated
  component space with feature loading vectors.

* `plotBiomarkerNetwork()`: Clustered heatmap of host gene to microbial
  taxon cross-omics Spearman correlations.

* `plotClassifierComparison()`: ROC curves comparing host-only,
  microbiome-only, and joint classifier performance.

* `plotSankey()`: Feature-flow diagram from omics layer through selected
  biomarkers to outcome classes.

* `generateReport()`: Print or save a structured text summary of all
  analysis results.

* `MOBResult` S4 class with slots for integrated scores, feature loadings,
  biomarker table, classifier performance, and parameters. Accessor
  methods: `biomarkers()`, `performance()`, `integrationScores()`,
  `featureLoadings()`.
