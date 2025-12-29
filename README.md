Testgit
================
Osafu Egbon
2025-12-29

# Tracking dynamic cell-cell communication using a stochastic ordering framework for spatially resolved transcriptomics

Intercellular communication is fundamental to tissue homeostasis and
disease progression, yet quantitatively resolving its dynamics across
biological systems remains challenging. Although single-cell and spatial
transcriptomic technologies now enable high-resolution profiling of
cell-cell interactions, most computational approaches provide static
summaries that overlook temporal and contextual variability. Here, we
introduce StochasticCC, a novel probabilistic framework that models
ligand-receptor signaling using full bivariate expression distributions.
By leveraging stochastic ordering and Bayesian network modeling,
StochasticCC accounts for correlation structure between ligands and
receptors, adjusts for spatial and technical variation, and reduces
false discoveries in heterogeneous or multi-modal datasets. In
simulations, controlled biological studies, and benchmarking against
four major communication frameworks, it consistently outperformed
mean-based and product scoring methods. Applying StochasticCC to induced
ulcerative colitis and tamoxifen-induced lung injury, we identify
conserved fibroblast activation programs converging on FAK/Src-centered
signaling modules, mediated by distinct ligand-receptor pairs such as
Cxcl12-Itga5 and Spp1-Itgav/Itga8. This cross-organ convergence reveals
a shared mechanotransduction axis underlying fibrotic remodeling,
demonstrating that stochastic communication dynamics encode conserved
signaling logic across inflammatory diseases. Collectively, StochasticCC
provides a generalizable framework for uncovering conserved signaling
motifs across tissues and disease states.

# Demo (Dynamics of communication)

In the following, we demonstrate how the StochasticCC is used to
identify dynamic cell-cell communication

# Demo (Static communication)

In the following, we demonstrate how the StochasticCC is used to
identify dynamic cell-cell communication

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax
for authoring HTML, PDF, and MS Word documents. For more details on
using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that
includes both content as well as the output of any embedded R code
chunks within the document. You can embed an R code chunk like this:

``` r
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

## Including Plots

You can also embed plots, for example:

![](Untitled_files/figure-gfm/pressure-1.png)<!-- -->

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
