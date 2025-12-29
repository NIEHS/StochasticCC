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
HH = NetworkGet(metaData,datExpr)

# Extract Processed data

datExpr = HH$datExpr_withNeigborhood
Centers2 = HH$Centers
ax_reduce_META = HH$ax_reduce_META
datExpr$day.harvested = factor(ax_reduce_META$Sample_type,levels = c("Healthy","DSS3", "DSS9", "DSS21"))%>%as.numeric()

# Estimate network
set.seed(12345)
mst = ClusterToTree(Centers = Centers2)


# Estimate phi for day, say 3

DEGs = FindDEGs(datExpr,
                mst,
                day=3 
               )

# Collect outputs
Result  = DEGs$Result_Cancer_day
subdata = DEGs$subdata 
ctype_day1  =  DEGs$ctype_day1
LRInData = DEGs$LRInData
META_subb = DEGs$META_subb


# Extract effect

TreeTemp_data_effect = Teffect(Result,
                               mst,
                               subdata,
                               ctype_day1)

TreeTemp_data_effect = TreeTemp_data_effect$TreeTemp_data_effect
TreeTemp_data = TreeTemp_data_effect$TreeTemp_data


# Get Null distribution of communication scores

AA = getNull(TreeTemp_data_effect,LRInData)

DN = AA$DN
LR0_effect = AA$LR0_effect

 
# Estimate cell-cell communication

CTn = c("Stem cells","Goblet 1", "Colonocytes", "T (Cd4+ Ccr7+)", "Treg","Fibro 1")

BB = ComputeCCS(TreeTemp_data_effect,
                LRInData,
                CTn,
                DN,
                LR0_effect)


ResultMatrix     =     BB$ResultMatrix
ResultMatrixPval = BB$ResultMatrixPval
```
