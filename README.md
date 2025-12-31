# Tracking dynamic cell-cell communication using a stochastic ordering framework for spatially resolved transcriptomics
We introduce StochasticCC, a probabilistic framework that models ligand-receptor cell-cell signaling using full bivariate distributions of expression. The framework accounts for the correlation structure between ligands and receptors and adjusts for spatial and technical variations in spatial transcriptomics/sc-RNA seq datasets. In this page, we demonstrate how to use the framework to analyze a MERFISH spatial transcriptomics dataset from a mouse model that involved the developmental stages of Ulcerative Colitis.
# Demo (Static communication)
In the following, we demonstrate how the StochasticCC is used to identify dynamic cell-cell communication

## Load Demo data

```{R}
load("~/Downloads/datExpr.Rdata")
load("~/Downloads/metaData.Rdata")
source("~/Downloads/DSFMix-main/EMTVenosa/StochasticOdering/StochatCCGuthub.R")

# Convert string to numeric days
metaData$Sample_type  = factor(metaData$Sample_type,levels = c("Healthy","DSS3", "DSS9", "DSS21"))%>%as.numeric()

# Cell types of interest
CTn = c("Stem cells", "Goblet 1", "Colonocytes", "T (Cd4+ Ccr7+)", "Treg", "Fibro 1")
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
<img width="4491" height="2323" alt="Picture5" src="https://github.com/user-attachments/assets/60c1e84b-f1ce-444d-8656-775764b66d5b" />

## Descriptive Plots on Network
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

plotTree(mst,Centers$Slice_ID,vertex.size =Ref$nn, 
         main = "",Lab = F,noLegend = F,
         legend.size = 5,
         edge_color="grey",
         edge_alpha=.1,
         cols = c25)
```

<img width="4624" height="1292" alt="Picture6" src="https://github.com/user-attachments/assets/05ddbc85-b250-4c03-aca0-196e257f828a" />

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


# Demo (Dynamics of communication)
In the following, we demonstrate how the StochasticCC is used to identify dynamic cell-cell communication

```{R}
##################
# Test ordering
################

Test = list(
  c(1,2,3), # Time T1 -> T3
  c(3,2,1), # Time T3 -> T2
  c(2,1,3) # Time T2  -> T3
)
```
### Estimate communication parameters
```{R}
RRe  = GetDynamicNtwk(Re = R,
                      LRInData = LRInData[1:10,],
                      ctype_from = "Colonocytes",
                      ctype_to = "Stem cells",
                      uniq_day = Test,
                      nCores = 9,
                      usePvalue =TRUE,
                      NoG = 5
)
```
### Plot dynamics 

```{R}
#################
#### Plot Network 
##################

OR  = RRe$OR
df <- data.frame(
  from = OR[,1],
  to   =  OR[,3],
  ligand_receptor = rownames(OR)
)

# Create graph

g <- graph_from_data_frame(df, directed = TRUE)
E(g)$width <- 1/as.numeric(OR[,4])

pairwise_colors <- c(
  "T1T1" = "#000000",   # Black
  "T1T2" = "#E69F00",   # Orange
  "T1T3" = "#56B4E9",   # Sky Blue
  "T2T1" = "#009E73",   # Bluish Green
  "T2T2" = "#F0E442",   # Yellow
  "T2T3" = "#0072B2",   # Blue
  "T3T1" = "#D55E00",   # Vermillion
  "T3T2" = "#CC79A7",   # Reddish Purple
  "T3T3" = "#999999"    # Gray
)

# Assign colors to edges
edge_groups <- as.factor(paste0(df$from, df$to))
E(g)$color <- pairwise_colors[as.character(edge_groups)]

# Plot with edge labels
plot(g,
     edge.label = E(g)$ligand_receptor,  # assign edge labels
     edge.arrow.size = 0.5,              # arrow size
     vertex.color = "lightblue",
     vertex.size = 30,
     edge.label.cex = 0.9,
     vertex.label.cex = 1.5)

```

<img width="1688" height="1738" alt="image" src="https://github.com/user-attachments/assets/98a58d81-fa18-4f01-bdcc-9713c4afed08" />

### Plot simplex for c(1,2,3) and Ligand-recpeptor pair `Adam12^Itga9`

```{R}
# Plot simplex for c(1,2,3) and Ligand-recpeptor pair `Adam12^Itga9`

hk = paste0("T",c(1,2,3),collapse = "")

df = RRe$Simplex[[hk]]$`Adam12^Itga9`

ggtern(df %>% unique(), aes(x = T1, y = T2, z = T3)) +
  geom_point(size = 3, color = "blue", alpha = 0.6) +
  theme_bw() +
  labs(
    title = "Simplex Plot of T1, T2, T3",
    T = "T2", L = "T1", R = "T3"
  )
```

<img width="1988" height="1692" alt="image" src="https://github.com/user-attachments/assets/164f48eb-9273-4bd4-90e4-4b28538d449f" />
