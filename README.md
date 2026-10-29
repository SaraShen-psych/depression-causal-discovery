# Causal Discovery in Depressive Symptoms

This repository contains R code for cross-sectional network analysis and causal discovery of depressive symptoms using data from Wave 2 of the Survey of Health, Ageing and Retirement in Europe (SHARE).

## Methods

The analysis includes:

- **undirected network estimation**:  Estimation of an Ising network for binary depressive symptoms using the enhanced LASSO (eLASSO) method.
- **causal network estimation**: Application of the Fast Causal Inference (FCI) algorithm to estimate a partial ancestral graph (PAG).
- **Centrality comparison**: Comparison of strength centrality in the undirected network with in-degree and out-degree in the causal network.
- **Correlation analysis**: Spearman correlation between strength centrality and out-degree across symptoms.
- **Generalized Adjustment Criterion (GAC)** to explore covariate adjustment identifiability in the full-sample FCI PAG.

## Data availability

The SHARE dataset is not included in this repository due to data access restrictions.

Researchers can obtain access to SHARE data through the [SHARE website](https://share-eric.eu/data/).

## Requirements

The analyses were conducted in R using packages including `haven`, `dplyr`, `IsingFit`, `qgraph`, `bootnet`, `igraph`, `pcalg`, `ggplot2`, `reshape2`, and `ggrepel`.

## Reproducing the analysis

1. Obtain authorized access to SHARE Wave 2 data.
2. Place the required dataset in a local `data/` directory.
3. Set the R working directory to the repository root.
4. Run `network_causal discovery_analysis code.R`.

The script generates network estimates, stability analyses, and visualizations.
