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

# Tracking dynamic cell-cell communication using a stochastic ordering framework for spatially resolved transcriptomics
Intercellular communication is fundamental to tissue homeostasis and disease progression, yet quantitatively resolving its dynamics across biological systems remains challenging. Although single-cell and spatial transcriptomic technologies now enable high-resolution profiling of cell-cell interactions, most computational approaches provide static summaries that overlook temporal and contextual variability. Here, we introduce StochasticCC, a novel probabilistic framework that models ligand-receptor signaling using full bivariate expression distributions. By leveraging stochastic ordering and Bayesian network modeling, StochasticCC accounts for correlation structure between ligands and receptors, adjusts for spatial and technical variation, and reduces false discoveries in heterogeneous or multi-modal datasets. In simulations, controlled biological studies, and benchmarking against four major communication frameworks, it consistently outperformed mean-based and product scoring methods. Applying StochasticCC to induced ulcerative colitis and tamoxifen-induced lung injury, we identify conserved fibroblast activation programs converging on FAK/Src-centered signaling modules, mediated by distinct ligand-receptor pairs such as Cxcl12-Itga5 and Spp1-Itgav/Itga8. This cross-organ convergence reveals a shared mechanotransduction axis underlying fibrotic remodeling, demonstrating that stochastic communication dynamics encode conserved signaling logic across inflammatory diseases. Collectively,  StochasticCC provides a generalizable framework for uncovering conserved signaling motifs across tissues and disease states.

# Demo (Dynamics of communication)
In the following, we demonstrate how the StochasticCC is used to identify dynamic cell-cell communication


# Demo (Static communication)
In the following, we demonstrate how the StochasticCC is used to identify dynamic cell-cell communication

```{R}
load("~/Downloads/datExpr.Rdata")
load("~/Downloads/metaData.Rdata")
source("~/Downloads/DSFMix-main/EMTVenosa/StochasticOdering/StochatCCGuthub.R")

# Head of metadata
head(metaData)

# First few rows and column of expression data
datExpr[1:5,1:5]
```
<img width="2688" height="524" alt="image" src="https://github.com/user-attachments/assets/d0521e36-28f7-4437-afb6-5ac2bdec7521" />
<img width="802" height="304" alt="image" src="https://github.com/user-attachments/assets/ee2a684d-8400-45f2-8ec5-b9fc33bb62b0" />

```{R}

#########

HH = NetworkGet(metaData,datExpr)

### Extract Processed data


datExpr = HH$datExpr_withNeigborhood
Centers2 = HH$Centers
ax_reduce_META = HH$ax_reduce_META

```


```{R}
## Convert time to numeric id and estimate network

datExpr$day.harvested = factor(ax_reduce_META$Sample_type,levels = c("Healthy","DSS3", "DSS9", "DSS21"))%>%as.numeric()

########


# Estimate network
set.seed(12345)
mst = ClusterToTree(Centers = Centers2)

# NeighXCell_id is the node id

# Count number of cells per node
Ref = data.frame(NeighXCell_id = datExpr$NeighXCell_id) %>% 
                group_by(NeighXCell_id) %>%
                summarise(nn=n())

####################
# Plot nodes according to ther cell types
####################

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
Centers =  data.frame (ax_reduce_META,
                       NeighXCell_id = datExpr$NeighXCell_id ) %>% 
  group_by(NeighXCell_id) %>% summarise_all(Mode) %>% ungroup() %>% 
  dplyr::select(-NeighXCell_id)

# Extr
library(forstringr)
Centers$Slice_ID_day = Centers$Slice_ID %>%str_extract_part("_D",before = F)%>%str_extract_part("_m",before = T)
Centers$Slice_ID_day = factor(Centers$Slice_ID_day,levels = c(0,3,9,21))
Centers$Slice_ID = Centers$Slice_ID %>%str_extract_part("1_",before = F) %>%as.factor()



plotTree(mst,Centers$Tier3,vertex.size =Ref$nn, 
         main = "",Lab = F,noLegend = F,
         legend.size = 5,
         edge_color="grey",
         edge_alpha=.1,
         cols = c25)

plotTree(mst,Centers$Slice_ID,vertex.size =Ref$nn, 
         main = "",Lab = F,noLegend = F,
         legend.size = 5,
         edge_color="grey",
         edge_alpha=.1,
         cols = c25)

```

<img width="1542" height="800" alt="image" src="https://github.com/user-attachments/assets/74d9f5df-a73e-4af4-855c-68a0e73a0b39" />
<img width="1542" height="800" alt="image" src="https://github.com/user-attachments/assets/00a8cdbe-32d1-4d85-a7ad-06f40433dd20" />


