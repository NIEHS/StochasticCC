# Tracking dynamic cell-cell communication using a stochastic ordering framework for spatially resolved transcriptomics
Intercellular communication is fundamental to tissue homeostasis and disease progression, yet quantitatively resolving its dynamics across biological systems remains challenging. Although single-cell and spatial transcriptomic technologies now enable high-resolution profiling of cell-cell interactions, most computational approaches provide static summaries that overlook temporal and contextual variability. Here, we introduce StochasticCC, a novel probabilistic framework that models ligand-receptor signaling using full bivariate expression distributions. By leveraging stochastic ordering and Bayesian network modeling, StochasticCC accounts for correlation structure between ligands and receptors, adjusts for spatial and technical variation, and reduces false discoveries in heterogeneous or multi-modal datasets. In simulations, controlled biological studies, and benchmarking against four major communication frameworks, it consistently outperformed mean-based and product scoring methods. Applying StochasticCC to induced ulcerative colitis and tamoxifen-induced lung injury, we identify conserved fibroblast activation programs converging on FAK/Src-centered signaling modules, mediated by distinct ligand-receptor pairs such as Cxcl12-Itga5 and Spp1-Itgav/Itga8. This cross-organ convergence reveals a shared mechanotransduction axis underlying fibrotic remodeling, demonstrating that stochastic communication dynamics encode conserved signaling logic across inflammatory diseases. Collectively,  StochasticCC provides a generalizable framework for uncovering conserved signaling motifs across tissues and disease states.

# Demo (Dynamics of communication)
In the following, we demonstrate how the StochasticCC is used to identify dynamic cell-cell communication


# Demo (Static communication)
In the following, we demonstrate how the StochasticCC is used to identify dynamic cell-cell communication

## Load Demo data

```{R}
load("~/Downloads/datExpr.Rdata")
load("~/Downloads/metaData.Rdata")
source("~/Downloads/DSFMix-main/EMTVenosa/StochasticOdering/StochatCCGuthub.R")
```

## Estimate communication parameters
```{R}
R = StochasticCC(metaData,
                 datExpr,
                 day_var ="Sample_type",
                 CTn = CTn
                )
```
## Plot Results
