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
```{R}
#########################
######## Results & Plots 
#########################
day=1
# Plot communication network 
CC = GetNetWk(R$ResultMatrix[[day]],
              R$ResultMatrixPval[[day]],
              R$subdata[[day]],
              mst,
              cutoff=20)


```
<img width="5023" height="1902" alt="Picture1" src="https://github.com/user-attachments/assets/6b1ec824-6bae-4735-821d-79932b55e03d" />


```{R}
# Show top Ligand-receptor pairs
grid.arrange(CC$tg)
df = CC$df
sa = union(df$from,df$to)

```
<img width="5234" height="898" alt="Picture2" src="https://github.com/user-attachments/assets/f2f4af6c-3aa0-4e7b-ba58-b3ea34df5f8e" />

```{R}
# Plot cell types on tissue

Graphh(R$Result[[day]],
       R$subdata[[day]],
       R$META_subb[[day]],
       mst,
       tissue =tissue,
       sa, # Cell types to view
       #antibody = "Itgb2"
       antibody = NULL
       
)

```
<img width="4565" height="1127" alt="Picture3" src="https://github.com/user-attachments/assets/80998ef8-68ba-4d06-95d3-5789663dbc02" />

```{R}
# Plot effects on tissue

Graphh(R$Result[[day]],
       R$subdata[[day]],
       R$META_subb[[day]],
       mst,
       tissue =tissue,
       sa, # Cell types to view
       #antibody = "Itgb2"
       antibody = NULL
       
)

```

<img width="4358" height="1141" alt="Picture4" src="https://github.com/user-attachments/assets/c05b3837-3417-409a-b834-adaf03d3e770" />


```{R}

#######
HeatMap(
  ResultMatrix =R$ResultMatrix[[day]],
  ResultMatrixPval =R$ResultMatrixPval[[day]] ,
  sa=sa,
  fontsize_row = 3.5,
  fontsize_col = 7,
  zoom=1,
  show.onlySig = FALSE,
  toCell = NULL)

```


```{R}
########

# NeighXCell_id is the node id

# Count number of cells per node
Ref = data.frame(NeighXCell_id = R$datExpr$NeighXCell_id) %>% 
  group_by(NeighXCell_id) %>%
  summarise(nn=n())

# Plot nodes on the network according to their cell types


Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
Centers =  data.frame (R$ax_reduce_META,
                       NeighXCell_id = R$datExpr$NeighXCell_id ) %>% 
  group_by(NeighXCell_id) %>% summarise_all(Mode) %>% ungroup() %>% 
  dplyr::select(-NeighXCell_id)

# Extract
library(forstringr)
Centers$Slice_ID_day = Centers$Slice_ID %>%str_extract_part("_D",before = F)%>%str_extract_part("_m",before = T)
Centers$Slice_ID_day = factor(Centers$Slice_ID_day,levels = c(0,3,9,21))
Centers$Slice_ID = Centers$Slice_ID %>%str_extract_part("1_",before = F) %>%as.factor()

```

```{R}

plotTree(mst,Centers$Tier3,vertex.size =Ref$nn, 
         main = "",Lab = F,noLegend = F,
         legend.size = 5,
         edge_color="grey",
         edge_alpha=.1,
         cols = c25)
```

```{R}
plotTree(mst,Centers$Slice_ID,vertex.size =Ref$nn, 
         main = "",Lab = F,noLegend = F,
         legend.size = 5,
         edge_color="grey",
         edge_alpha=.1,
         cols = c25)
```

```{R}
## Sort genes based on their SNR

Significantgene = R$Result[[1]]$SNR[1,] %>% sort(decreasing = T) %>% names
ge = R$Result[[1]]$treeEffect[,Significantgene[1]]

#################
# Plot phi on network
#################

plotTree(mst,ge,vertex.size =Ref$nn, 
         main = Significantgene[1],Lab = F,noLegend = F,
         legend.size = 5,
         edge_color="grey",
         edge_alpha=.1,
         cols = c25)

```
