
########################
# Packages Required
########################
library(quadprog)
library(INLA) # For Imputation
library(gamlss.spatial) # For SNR
library(tidyverse) # Data manipulation
library(igraph)
library(doParallel)
library(scales)
library(genie)
library(ggraph)
library(ar.matrix)
library(viridis)
library(forstringr)
library(reshape2)
library(Iso)
library(ggtern)

Hcluster <- function(DataTocluster,thresholdGini=0.2,k,ClusterName){
  # Goal: Function to cluster a data using hierarchical clustering
  # INPUT: 
  # DataTocluster - > data matrix to cluster
  # thresholdGini - > Gini hyperparameter
  # OUTPUT: Cluste labled class
  
  cluster  <- genie::hclust2(objects=as.matrix(DataTocluster), thresholdGini=thresholdGini)
  clust    <- cutree(cluster, k = k)
  
  clust = as.matrix(clust) 
  colnames(clust) = ClusterName
  ClusteredData = cbind(clust,DataTocluster)
  
  return(ClusteredData)
}

AddColumnToFCS = function(InputDirOfData, ColumToadd,ColName,OutputDirOfData){
  # Goal: Function to add a column to an FCS file
  
  h <- spade:::SPADE.read.FCS(InputDirOfData)
  
  params <- flowCore::parameters(h)
  params1 <- params
  desc <- description(h)
  pd <- pData(params1)
  
  in_data <-  exprs(h)
  channel_number <- ncol(in_data) + 1
  channel_id <- paste("$P", channel_number, sep = "")
  channel_name <- ColName
  channel_range <- k + 1
  plist <- matrix(c(channel_name, channel_name, channel_range,
                    0, channel_range - 1))
  rownames(plist) <- c("name", "desc", "range", "minRange",
                       "maxRange")
  colnames(plist) <- c(channel_id)
  pd <- rbind(pd, t(plist))
  pData(params1) <- pd
  out_data <- cbind(in_data, ColumToadd)
  colnames(out_data) <- c(colnames(in_data),ColName)
  out_frame <- flowFrame(out_data, params1, description = desc)
  keyval <- list()
  keyval[[paste("$P", channel_number, "B", sep = "")]] <- "32"
  keyval[[paste("$P", channel_number, "R", sep = "")]] <- toString(channel_range)
  keyval[[paste("$P", channel_number, "E", sep = "")]] <- "0,0"
  keyval[[paste("$P", channel_number, "N", sep = "")]] <- channel_name
  keyval[[paste("$P", channel_number, "S", sep = "")]] <- channel_name
  keyword(out_frame) <- keyval
  write.FCS(out_frame, OutputDirOfData)
}


IntegrateTreeWithPriorTree <- function(ExprsData,ClusterCol,PriorColumn,  weighted.prob = 0.5,scale=0.1){
  
  # Determine the clusters with maximum similar cells between the data cluster
  # and manual clusters
  in_data   = data.frame(ExprsData)
  cluster   = ClusterCol
  spatialid = PriorColumn
  grp = c(cluster,spatialid)
  ######################
  in_data_grp = in_data %>% group_by_at(grp)
  in_data_grp_crossTab = in_data_grp %>% summarise(n = n()) %>% 
    ungroup() %>% group_by_at(cluster)%>%
    filter(n==max(n))
  
  n2=runif(nrow(in_data_grp_crossTab),0,1)
  in_data_grp_crossTab = in_data_grp_crossTab %>% ungroup%>%
    mutate(n2=n2)%>% group_by_at(cluster)%>%
    filter(n2==max(n2))
  
  
  # Compute the weighted average between centers of data and prior centers
  
  AdjustedCenters = matrix(0,nrow=nrow(in_data_grp_crossTab),ncol=ncol(in_data)-1)
  
  colnames(AdjustedCenters) = colnames(in_data[,!(colnames(in_data) %in% cluster)])
  
  auxA = in_data_grp_crossTab[,cluster] %>% as.matrix() %>% as.vector()
  auxB = in_data_grp_crossTab[,spatialid] %>% as.matrix() %>% as.vector()
  
  for(i in 1:nrow(in_data_grp_crossTab)){
    Aux_data_data = in_data %>% filter_at(vars(cluster),all_vars(.==auxA[i])) %>% select_at(vars(-cluster) )%>% mutate_at(vars(spatialid),function(x)x=0.001) %>% colMeans()
    Aux_data_prior = in_data %>% filter_at(vars(spatialid),all_vars(.==auxB[i])) %>% select_at(vars(-cluster) )%>% mutate_at(vars(spatialid),function(x)x=x*scale) %>% colMeans()
    
    aAux_data_prior = Aux_data_prior
    Aux_data_prior  = rep(0,length(aAux_data_prior))
    names(Aux_data_prior) = names(aAux_data_prior)
    Aux_data_prior[spatialid] =aAux_data_prior[spatialid]
    AdjustedCenters[i,] = weighted.prob * Aux_data_prior  +(1-weighted.prob)*Aux_data_data
    
  }
  return(list(UpdatedCenters=AdjustedCenters,ClusterAndPriorID = in_data_grp_crossTab))
}

ClusterToTree <- function(Centers,weighted = TRUE){
  # Derive a tree from a center
  # INPUT: Centers
  # OUTPUT: mst (igraph) object.
  adjacency <- dist(Centers, method = "manhattan")
  full_graph <- graph_from_adjacency_matrix(as.matrix(adjacency), mode = "undirected", 
                                            weighted = weighted)
  mst <- minimum.spanning.tree(full_graph)
  
  return(mst)
}


ReconectDisconectedNetwk = function(mstdel){
  m2 = as_adjacency_matrix(mstdel)
  
  while(sum(rowSums(m2)==0)>0){
    
    Disc.Node =  which(rowSums(m2)==0)
    layoutCoord = layout_nicely(mstdel)
    dist = proxy::dist(layoutCoord[Disc.Node[1],,drop=F] ,layoutCoord)
    aux = which(as.vector(dist)==sort(as.vector(dist),decreasing =F)[2])
    mstdel= add_edges(mstdel, edges = c(Disc.Node[1], aux))
    m2 = as_adjacency_matrix(mstdel)
  }
  return(mstdel)
}

plotTree <- function(mst,PlotVar,
                     #community=FALSE,
                     main="Weighted",Lab=T,noLegend=F,
                     vertex.size=3,vertex.label.cex=1,
                     sizelegend="none",
                     limits = NULL,
                     cols=NULL,
                     legend.size = 0.2,
                     legend.text.size=10,
                     legend.size.alpha=1,
                     seed=11234){
  # Plot function
  # INPUT: mst -> graph object(igraph),community-> whether to plot community (logical),
  # main-> title,... Other parameters are igraph parameters. 
  # OUTPU: plot
  # Note, if PlotVar is a factor, plot is made from igraph base plot, if countitnous, ggplot is used.
  #if(community){
  #  a = cluster_fast_greedy(mst)
  #  V(mst)$id_com = membership(a)
  #  set.seed = seed
  #  plot(a,mst,vertex.size=8,vertex.label.cex=1,main=main)
  #}else{
  
  if(is.factor(PlotVar)){
    # a=colorRampPalette(c("blue","yellow","red"))
    #PlotVar = PlotVar %>%as.numeric()
    # aux_PlotVar =unique(PlotVar) %>% sort()
    # pal =  a(length(aux_PlotVar))[PlotVar]
    # set.seed = seed
    #  plot(mst,vertex.size=3,vertex.label.cex=1,vertex.label=NA,
    #       vertex.color =pal,main=main)
    #  legend("bottomright",bty = "n",
    #         legend=(aux_PlotVar),
    #        fill=(a(length(aux_PlotVar))), cex=.8,border=NA, horiz = TRUE,text.width = 0.001)
    p= ggraph(mst, layout = "stress") + 
      geom_edge_link() + 
      geom_node_point(aes(color=PlotVar,size=vertex.size))+ guides(size=sizelegend)+
      scale_color_manual(values = cols)+labs(color="",x="",y="",title=main,size="")+
      #  geom_node_text(label =V(mst),
      #                 colour = 'black', size=3,
      #                 show.legend = FALSE, family = "serif")+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            rect = element_blank(),
            legend.key.width= unit(0.2, 'cm'),
            legend.text = element_text(size=legend.text.size))+
      guides(colour = guide_legend(override.aes = list(size=legend.size)),
             size = guide_legend(override.aes = list(alpha = 0.5))
      )
    
    #   par(old.par)
    if(Lab){
      p = p +  geom_node_text(label =V(mst),
                              colour = 'black', size=3,
                              show.legend = FALSE, family = "serif")
    }
    if(noLegend){
      p=p+ guides(color="none")
    }
    if(sizelegend=="none"){
      p=p+ guides(size="none")
    }
    p
  }else{
    
    set.seed = seed
    pal <- gradient_n_pal(brewer_pal(palette = "Spectral", direction = -1)(7))
    p= ggraph(mst, layout = "stress") + 
      geom_edge_link() + 
      geom_node_point(aes(color=PlotVar,size=vertex.size))+ guides(size=sizelegend)+
      scale_color_distiller(palette = "Spectral",limits=limits)+labs(color="",x="",y="",title=main,size="")+
      #  geom_node_text(label =V(mst),
      #                 colour = 'black', size=3,
      #                 show.legend = FALSE, family = "serif")+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            rect = element_blank(),
            legend.key.width= unit(0.2, 'cm'),
            legend.text = element_text(size=legend.text.size))+
      guides(#colour = guide_legend(override.aes = list(size=legend.size)),
        
        size = guide_legend(override.aes = list(alpha = alpha)))
    
    #   par(old.par)
    if(Lab){
      p = p +  geom_node_text(label =V(mst),
                              colour = 'black', size=3,
                              show.legend = FALSE, family = "serif")
    }
    if(noLegend){
      p=p+ guides(color="none")
    }
    if(sizelegend=="none"){
      p=p+ guides(size="none")
    }
    p
  }
}


plotTree <- function(mst,PlotVar,
                     #community=FALSE,
                     main="Weighted",Lab=T,noLegend=F,
                     vertex.size=3,vertex.label.cex=1,
                     sizelegend="none",
                     edge_color="black",
                     edge_alpha=1,
                     limits = NULL,
                     cols=NULL,
                     legend.size = 0.2,
                     legend.text.size=10,
                     legend.size.alpha=1,
                     seed=11234){
  # Plot function
  # INPUT: mst -> graph object(igraph),community-> whether to plot community (logical),
  # main-> title,... Other parameters are igraph parameters. 
  # OUTPU: plot
  # Note, if PlotVar is a factor, plot is made from igraph base plot, if countitnous, ggplot is used.
  #if(community){
  #  a = cluster_fast_greedy(mst)
  #  V(mst)$id_com = membership(a)
  #  set.seed = seed
  #  plot(a,mst,vertex.size=8,vertex.label.cex=1,main=main)
  #}else{
  
  if(is.factor(PlotVar)){
    # a=colorRampPalette(c("blue","yellow","red"))
    #PlotVar = PlotVar %>%as.numeric()
    # aux_PlotVar =unique(PlotVar) %>% sort()
    # pal =  a(length(aux_PlotVar))[PlotVar]
    # set.seed = seed
    #  plot(mst,vertex.size=3,vertex.label.cex=1,vertex.label=NA,
    #       vertex.color =pal,main=main)
    #  legend("bottomright",bty = "n",
    #         legend=(aux_PlotVar),
    #        fill=(a(length(aux_PlotVar))), cex=.8,border=NA, horiz = TRUE,text.width = 0.001)
    p= ggraph(mst, layout = "stress") + 
      geom_edge_link(color=edge_color,alpha=edge_alpha) + 
      geom_node_point(aes(color=PlotVar,size=vertex.size))+ guides(size=sizelegend)+
      scale_color_manual(values = cols)+labs(color="",x="",y="",title=main,size="")+
      #  geom_node_text(label =V(mst),
      #                 colour = 'black', size=3,
      #                 show.legend = FALSE, family = "serif")+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            rect = element_blank(),
            legend.key.width= unit(0.2, 'cm'),
            legend.text = element_text(size=legend.text.size))+
      guides(colour = guide_legend(override.aes = list(size=legend.size)),
             size = guide_legend(override.aes = list(alpha = 0.5))
      )
    
    #   par(old.par)
    if(Lab){
      p = p +  geom_node_text(label =V(mst),
                              colour = 'black', size=3,
                              show.legend = FALSE, family = "serif")
    }
    if(noLegend){
      p=p+ guides(color="none")
    }
    if(sizelegend=="none"){
      p=p+ guides(size="none")
    }
    p
  }else{
    
    set.seed = seed
    pal <- gradient_n_pal(brewer_pal(palette = "Spectral", direction = -1)(7))
    p= ggraph(mst, layout = "stress") + 
      geom_edge_link(color=edge_color,alpha=edge_alpha) + 
      geom_node_point(aes(color=PlotVar,size=vertex.size))+ guides(size=sizelegend)+
      scale_color_distiller(palette = "Spectral",limits=limits)+labs(color="",x="",y="",title=main,size="")+
      #  geom_node_text(label =V(mst),
      #                 colour = 'black', size=3,
      #                 show.legend = FALSE, family = "serif")+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            rect = element_blank(),
            legend.key.width= unit(0.2, 'cm'),
            legend.text = element_text(size=legend.text.size))+
      guides(#colour = guide_legend(override.aes = list(size=legend.size)),
        
        size = guide_legend(override.aes = list(alpha = alpha)))
    
    #   par(old.par)
    if(Lab){
      p = p +  geom_node_text(label =V(mst),
                              colour = 'black', size=3,
                              show.legend = FALSE, family = "serif")
    }
    if(noLegend){
      p=p+ guides(color="none")
    }
    if(sizelegend=="none"){
      p=p+ guides(size="none")
    }
    p
  }
}


c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

c35 <- c(
  "dodgerblue2", "green4", "#6A3D9A", "#FF7F00", "black",
  "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F",
  "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1",
  "steelblue4", "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown",
  # Added 10 more
  "mediumpurple", "turquoise4", "darkseagreen", "tomato", "darkmagenta",
  "lightcoral", "midnightblue", "olivedrab", "sienna", "cyan4", "#E31A1C"
)
plotScatter = function(Xcord ,
                       Ycord,
                       Gene,
                       size=1,
                       main,
                       pal="Accent",
                       legend.size = 0.2,
                       legend.text.size=10,
                       noLegend=FALSE,
                       limits=NULL,
                       ManualColor =FALSE,
                       cols = NULL){
  
  if(!is.factor(Gene)){
    p = data.frame(Xcord,Ycord,Gene) %>% ggplot()+
      geom_point(aes(x=Xcord,y=Ycord,color=Gene),size=size)+
      scale_color_distiller(palette = "Spectral",limits=limits)+
      labs(color="",x="",y="",title=main)+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            rect = element_blank(),
            legend.key.width= unit(0.2, 'cm'),
            legend.text = element_text(size=legend.text.size))#+
    # +
    # guides(colour = guide_legend(override.aes = list(size=legend.size)))
  }else{
    
    p= data.frame(Xcord,Ycord,Gene) %>% ggplot()+
      geom_point(aes(x=Xcord,y=Ycord,color=Gene),size=size)+
      scale_color_brewer(palette=pal)+
      labs(color="",x="",y="",title=main)+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            rect = element_blank(),
            legend.key.width= unit(0.2, 'cm'),
            #legend.key.size= unit(legend.size, 'cm'),
            legend.text = element_text(size=legend.text.size)
      )+
      guides(colour = guide_legend(override.aes = list(size=legend.size)))
  }
  if(noLegend){
    p=p+ guides(color="none")
  }
  if(ManualColor){
    p= p+scale_color_manual(values =cols )
  }
  p
}

#cv = function(v) {return((sd(v)*100)/mean(v))} 

ComputeCenter <- function(ExprsData,Genes,ClusterCol="cluster",f = mean){
  # Functions to compute centers
  # INPUT: ExprsData -> Expression data (cellxgenes), Genes-> Genes to which computation should be made
  #       ClusterCol -> Column with cluster labels, f-> Function to summarise for each cluster
  # OUTPU: Data frame with each Genes summarised with function f for each unique ClusterCol.
  
  ExprsData = as.data.frame(ExprsData)
  if(is.null(Genes)) Genes = colnames(ExprsData)[!(colnames(ExprsData)%in%ClusterCol)]
  G = c(ClusterCol,Genes)
  ExprsData =  ExprsData %>% select_at(G)
  aux =  ExprsData %>% group_by_at(ClusterCol) %>% summarise(
    across(where(is.numeric), list(f))
  )%>% ungroup
  colnames(aux) = ExprsData %>% ungroup %>% select_at(G) %>% colnames()
  
  return(aux)
}

IntegrateTwoTrees = function(ExprsDataRef,ExprsData,ClusterColRef,ClusterCol){
  
  # Funtion to update a Tree based on a new expression data
  # INPUT: ExprsDataRef-> Reference expression data(cellxgenes), ExprsData-> new expr data, 
  #        ClusterColRef -> cluster column for  ExprsDataRef, ClusterCol-> cluster column for ExprsData.
  # OUTPUT: Updated centers, outlier clusters id, and outlier id on updated clusters
  #Note: The two impute expressions must have the same genes;
  
  in_data_A = as.data.frame(ExprsDataRef)
  in_data_B = as.data.frame(ExprsData)
  oo = sort(unique(in_data_A[,ClusterColRef]))
  in_data_grp_crossTab = matrix(0, nrow = length(oo),ncol = 2)
  
  for (i in 1:length(oo)) {
    
    A =  in_data_A %>% dplyr::filter_at(vars(ClusterColRef),all_vars(.==oo[i])) %>% select_at(all_of(vars(-ClusterColRef))) %>% colMeans()
    #A$id = 1
    o = in_data_B[,ClusterCol] %>% unique
    pv=vector("numeric",length = length(o))
    
    sim = pv = NULL
    for (x in 1:length(o)){ 
      
      xy = o[x]
      B =  in_data_B %>% dplyr::filter_at(vars(ClusterCol),all_vars(.==xy)) %>% select_at(vars(-ClusterCol)) %>% colMeans()
      
      C   =  rbind(A,B)
      
      sim = c(sim, as.vector(proxy::simil(C, method = "manhattan")))
      
      
      
    }
    
    max.table2 = o[which(sim==max(sim,na.rm = T))]
    in_data_grp_crossTab[i,] =c(oo[i],max.table2)
    
  }
  
  colnames(in_data_grp_crossTab) = c("cluster_A","cluster_B")
  in_data_grp_crossTab = as.data.frame(in_data_grp_crossTab)
  
  # order according to A
  in_data_grp_crossTab = in_data_grp_crossTab[order(in_data_grp_crossTab$cluster_A),]
  
  # Determine clusters in B that are not captured
  outlier_cluster = unique(in_data_B[,ClusterCol])[!(unique(in_data_B[,ClusterCol]) %in% in_data_grp_crossTab$cluster_B)]
  outlier_center = in_data_B %>% dplyr::filter_at(vars(ClusterCol),all_vars(. %in% outlier_cluster)) %>% group_by_at(ClusterCol) %>%
    summarise(
      across(where(is.numeric), list(mean))
    ) %>% ungroup %>%  select_at(vars(-ClusterCol))
  colnames(outlier_center) = gsub("_1","",colnames(outlier_center))
  outlier_spatialid = max(in_data_A[,ClusterColRef])+(1:length(outlier_cluster))
  
  
  # Compute the weighted average between centers of data and prior ceenters
  weighted.prob = 0.5
  AdjustedCenters = matrix(0,nrow=nrow(in_data_grp_crossTab),ncol=ncol(in_data_A)-1)
  colnames(AdjustedCenters) = colnames(in_data_A %>% select_at(vars(-ClusterColRef)))
  for(i in 1:nrow(in_data_grp_crossTab)){
    Aux_data_data = in_data_A %>% dplyr::filter_at(vars(ClusterColRef),all_vars(.==in_data_grp_crossTab$cluster_A[i])) %>% dplyr::select_at(vars(-ClusterColRef)) %>% colMeans()
    Aux_data_prior = in_data_B %>% dplyr::filter_at(vars(ClusterCol),all_vars(.==in_data_grp_crossTab$cluster_B[i]))%>% dplyr::select_at(vars(-ClusterCol)) %>% colMeans()
    
    AdjustedCenters[i,] = weighted.prob * Aux_data_prior  + (1-weighted.prob) * Aux_data_data
    
  }
  AdjustedCenters = rbind(AdjustedCenters,outlier_center)
  
  return(list(UpdatedCenters=AdjustedCenters,
              unmachedCluster = outlier_cluster,
              IdOfUnmachedClusterOnUpdatedCenters=outlier_spatialid ))
}


IntegrateTwoTrees.impute = function(ExprsDataRef,ExprsData,ClusterColRef,ClusterCol){
  # Function to update a Tree based on a new expression data if some Genes are different.
  # INPUT: ExprsDataRef-> Reference expression data(cellxgenes), ExprsData-> new expr data, 
  #        ClusterColRef -> cluster column for  ExprsDataRef, ClusterCol-> cluster column for ExprsData.
  # OUTPUT: Updated centers, outlier clusters id, outlier id on updated clusters, and imputed expressions
  LabRef  = colnames(ExprsDataRef)
  Lab     = colnames(ExprsData)
  InT     = intersect(Lab,LabRef)
  FLab    = LabRef[!(LabRef%in%InT)]
  FLabRef = Lab[!(Lab%in%InT)]
  for (i in seq_len(length(FLabRef))) ExprsDataRef[,FLabRef[i]] = NA
  for (i in seq_len(length(FLab))) ExprsData[,FLab[i]] = NA
  
  in_data_A.mis  = as.data.frame(ExprsDataRef)
  in_data_B.mis  = as.data.frame(ExprsData)
  
  in_data_A.mis$id = 1
  in_data_B.mis$id = 2
  in_data.mis = rbind(in_data_A.mis,in_data_B.mis)
  
  
  # Remove cols that are all Missing (that is, the variable is missing in both data)
  resB = colSums(in_data.mis,na.rm = TRUE)
  in_data.mis = in_data.mis[,names(resB)[resB !=0]]
  
  # Set Missing col with missing as respons variable 
  
  resA = colSums(in_data.mis)
  aux.namesOfMis = names(resA[is.na(resA)])
  FiTTed = list()
  for(i in seq_len(length(aux.namesOfMis))){
    Y = in_data.mis %>% dplyr::select_at(vars(aux.namesOfMis[1]))
    Y = Y+1e-5 # Computational reasons
    X = in_data.mis %>% select (-all_of(aux.namesOfMis))
    
    DATA_imputation = cbind(Y,X)
    
    formula = paste(colnames(X%>%select(-id)%>% select_at((vars(-ClusterColRef,-ClusterColRef)))),
                    collapse = "+")
    formula = paste0(colnames(Y),"~",formula)
    formula = paste0(formula,"+f(id,model='iid')")
    formula = paste0(formula,"+f(",ClusterColRef,",model='iid')")
    if(ClusterColRef!=ClusterCol) formula = paste0(formula,"+f(",ClusterCol,",model='iid')")
    formula = as.formula(formula)
    
    Result = inla(formula = formula,data = DATA_imputation,family = "gamma",control.predictor=list(compute=TRUE),
    )
    aux = Result$summary.fitted.values$mean[is.na(Y)] %>% exp
    Y[is.na(Y)] = aux
    FiTTed[colnames(Y)][[1]] = aux
    
    ####################################### Using semi-structured modeling
    #formula = paste(colnames(X),
    #                collapse = "+")
    #formula = paste0("nn(~",formula,",size=",5,",decay=",0.01,")")
    #formula = paste0(colnames(Y),"~",formula)
    #formula = as.formula(formula)
    #Result = gamlss(formula, data=DATA_imputation,family ="GA")
    #Y[is.na(Y)] = fitted(Result)[is.na(Y)]
    ########################################
    in_data.mis[,colnames(Y)] = Y
    # Remove already imputed variable
    
    aux.namesOfMis = aux.namesOfMis[aux.namesOfMis!=colnames(Y)]
    
  }
  ##### End of imputation #########
  
  in_data_A = in_data.mis %>% filter(id==1) %>% select(-id)
  in_data_B = in_data.mis %>% filter(id==2) %>% select(-id)
  
  aux = IntegrateTwoTrees(ExprsDataRef =in_data_A ,ExprsData = in_data_B,
                          ClusterColRef ="clusterr" ,
                          ClusterCol = "clusterr")
  
  return(list(aux,Imputed = FiTTed))
  
}



GetTreeVariableGenes <- function(mst,
                                 ExprsData,
                                 ClusterCol, 
                                 useWeight=TRUE, 
                                 Robust = TRUE,
                                 Model="GA",
                                 nCores=1){
  
  
  Idclust = sort(unique(ExprsData[,ClusterCol]))
  nclust = length(Idclust)
  nNodes = vcount(mst)
  
  if(nclust!=nNodes) stop("Nodes on data and graph are different")
  if(useWeight){  WeightGraphAdj = as_adjacency_matrix(mst,attr = "weight")
  }else{  WeightGraphAdj = as_adjacency_matrix(mst)}
  
  
  PrecWeightGraphAdj = Diagonal(nrow(WeightGraphAdj),x=rowSums(WeightGraphAdj!=0)) - WeightGraphAdj
  
  # Since the mst_graph_fromPriorKwldge was sorted according to in_data_grp_crossTab$cluster_A
  # We can link it to the orriginal data to get graph info into the data
  id = data.frame(cluster = Idclust, 
                  Graphid = seq_len(length(Idclust)))
  id[,ClusterCol] = id[,"cluster"]
  id = id %>% select(-cluster)
  
  # Update input data
  
  in_data_AA = ExprsData  %>% left_join(id,by = ClusterCol)
  
  # Calculate gene signal to noise ratio
  
  names_gene = in_data_AA %>% dplyr::select(-all_of(ClusterCol)) %>% select(-Graphid) %>% colnames()
  
  # Configure to hold predicted gene values
  in_data_AA.fitted = in_data_AA[,names_gene]
  if(is.na(sum(in_data_AA))) warning("Expression contains missing data")
  
  # Configured to hold Signal-to-noise ratio
  SNR = matrix(NA,ncol = length(names_gene))
  colnames(SNR) = names_gene
  
  #Log-like 
  BIC = matrix(NA,ncol = length(names_gene))
  colnames(BIC) = names_gene
  # Configure to hold tree-effect
  treeEffect = matrix(0,nrow =nrow(WeightGraphAdj) ,ncol = length(names_gene))
  colnames(treeEffect) = names_gene
  # Compute density to derive nodes weights
  
  # weight = in_data_AA %>% group_by_at(ClusterCol) %>% dplyr::summarise(n=n()) %>%
  #          group_by_at(ClusterCol) %>%
  #         dplyr::reframe(dense = n/nrow(in_data_AA),weight = 1/dense) %>% select(-dense)
  # weight$weight = weight$weight/sum(weight$weight)
  #weight = 1/dense
  #weight = weight/sum(weight)
  #  in_data_AA = in_data_AA%>% left_join(weight,by=ClusterCol)
  # Parallel computing
  
  
  cl <- makeCluster(nCores, outfile="")
  registerDoParallel(cl)
  
  
  ko_param <- foreach(i=seq_len(length(names_gene)),.errorhandling = "pass",
                      .packages = c("doParallel",
                                    "Rfast2",
                                    "tidyverse",
                                    "gamlss.spatial",
                                    "gamlss"
                      )
  ) %dopar% {
    
    cat("Modeling gene", (i/length(names_gene))*100,"%","\n" )
    
    
    model_dataFrame = in_data_AA %>% 
      dplyr::select(names_gene[i], Graphid) %>% 
      dplyr::rename(gene =names_gene[i]) %>% 
      dplyr::mutate(gene= gene+1e-02)
    
    if (sum(is.na(model_dataFrame$gene)) > 0 | sum(is.nan(model_dataFrame$gene))>0 ){
      warning("NA in Gene",names_gene[i], ". Hence will be neglected") 
      SNR[,names_gene[i]] = NA
      in_data_AA.fitted[, names_gene[i]] = NA
      treeEffect[,names_gene[i]]  = NA
      next
    }
    
    model_dataFrame$Graphid = as.factor(model_dataFrame$Graphid)
    
    
    if (i==length(names_gene)) cat(' Finalizing...',"\n")
    if(!Robust) {
      Model=Model
    }else{
      m = fitDist(model_dataFrame$gene,type = "realplus")
      Model = m$family[1]
    }
    
    
    m1 <- gamlss(gene ~ gmrf(Graphid,
                             precision=as.matrix(PrecWeightGraphAdj)),
                 data=model_dataFrame,
                 family = Model)
    
    
    # Get Singnal to Noise ratio
    sError = coef(m1$mu.coefSmo[[1]])[[1]] %>% exp
    sSpat = coef(m1$mu.coefSmo[[1]])[[2]] %>% exp
    
    SNR[,names_gene[i]] = sSpat/sError
    BIC[,names_gene[i]] =  m1$sbc #-gen.likelihood(m1)()
    # if(i%%500 ==0) write.csv(SNR[,1:i] ,file=paste0(getwd(),"/SNR.txt"))
    
    # Get smoothed/predicted gene expression
    in_data_AA.fitted[, names_gene[i]] = fitted(m1)
    # Get tree effects
    treeEffect[,names_gene[i]]  = getSmo(m1)$beta
    
    list(Gene= names_gene[i],
         SNR=SNR[,names_gene[i]],
         BIC = BIC[,names_gene[i]],
         fitted= in_data_AA.fitted[, names_gene[i]],
         treeEffect = treeEffect[,names_gene[i]]
    )
  }
  
  
  for (i in seq_len(length(names_gene))) {
    if(is.null(ko_param[[i]]$SNR)) next
    # Get smoothed/predicted gene expression
    SNR[,names_gene[i]] = ko_param[[i]]$SNR
    BIC[,names_gene[i]] = ko_param[[i]]$BIC
    # Get smoothed/predicted gene expression
    in_data_AA.fitted[, names_gene[i]] = ko_param[[i]]$fitted
    # Get tree effects
    treeEffect[,names_gene[i]]  = ko_param[[i]]$treeEffect
  }
  
  parallel::stopCluster(cl) 
  return(list(SNR = SNR,BIC=BIC, Fitted =in_data_AA.fitted, TreeEffect = treeEffect ))
}


GetTreeVariableGenesLocScale <- function(mst,
                                         ExprsData,
                                         ClusterCol, 
                                         useWeight=TRUE, 
                                         Robust = TRUE,
                                         Model="GA",
                                         nCores=1){
  
  
  Idclust = sort(unique(ExprsData[,ClusterCol]))
  nclust = length(Idclust)
  nNodes = vcount(mst)
  
  if(nclust!=nNodes) stop("Nodes on data and graph are different")
  if(useWeight){  WeightGraphAdj = as_adjacency_matrix(mst,attr = "weight")
  }else{  WeightGraphAdj = as_adjacency_matrix(mst)}
  
  
  PrecWeightGraphAdj = Diagonal(nrow(WeightGraphAdj),x=rowSums(WeightGraphAdj!=0)) - WeightGraphAdj
  
  # Since the mst_graph_fromPriorKwldge was sorted according to in_data_grp_crossTab$cluster_A
  # We can link it to the orriginal data to get graph info into the data
  id = data.frame(cluster = Idclust, 
                  Graphid = seq_len(length(Idclust)))
  id[,ClusterCol] = id[,"cluster"]
  id = id %>% select(-cluster)
  
  # Update input data
  
  in_data_AA = ExprsData  %>% left_join(id,by = ClusterCol)
  
  # Calculate gene signal to noise ratio
  
  names_gene = in_data_AA %>% dplyr::select(-all_of(ClusterCol)) %>% select(-Graphid) %>% colnames()
  
  # Configure to hold predicted gene values
  in_data_AA.fitted = in_data_AA[,names_gene]
  if(is.na(sum(in_data_AA))) warning("Expression contains missing data")
  
  # Configured to hold Signal-to-noise ratio
  SNR = matrix(NA,ncol = length(names_gene))
  SNRScale = matrix(NA,ncol = length(names_gene))
  colnames(SNR) = names_gene
  colnames(SNRScale) = names_gene
  # Configure to hold tree-effect
  treeEffect = matrix(0,nrow =nrow(WeightGraphAdj) ,ncol = length(names_gene))
  colnames(treeEffect) = names_gene
  # Compute density to derive nodes weights
  
  # weight = in_data_AA %>% group_by_at(ClusterCol) %>% dplyr::summarise(n=n()) %>%
  #          group_by_at(ClusterCol) %>%
  #         dplyr::reframe(dense = n/nrow(in_data_AA),weight = 1/dense) %>% select(-dense)
  # weight$weight = weight$weight/sum(weight$weight)
  #weight = 1/dense
  #weight = weight/sum(weight)
  #  in_data_AA = in_data_AA%>% left_join(weight,by=ClusterCol)
  # Parallel computing
  
  
  cl <- makeCluster(nCores, outfile="")
  registerDoParallel(cl)
  
  
  ko_param <- foreach(i=seq_len(length(names_gene)),.errorhandling = "pass",
                      .packages = c("doParallel",
                                    "Rfast2",
                                    "tidyverse",
                                    "gamlss.spatial",
                                    "gamlss"
                      )
  ) %dopar% {
    
    cat("Modeling gene", (i/length(names_gene))*100,"%","\n" )
    
    
    model_dataFrame = in_data_AA %>% 
      dplyr::select(names_gene[i], Graphid) %>% 
      dplyr::rename(gene =names_gene[i]) %>% 
      dplyr::mutate(gene= gene+1e-02)
    
    if (sum(is.na(model_dataFrame$gene)) > 0 | sum(is.nan(model_dataFrame$gene))>0 ){
      warning("NA in Gene",names_gene[i], ". Hence will be neglected") 
      SNR[,names_gene[i]] = NA
      in_data_AA.fitted[, names_gene[i]] = NA
      treeEffect[,names_gene[i]]  = NA
      next
    }
    
    model_dataFrame$Graphid = as.factor(model_dataFrame$Graphid)
    
    
    if (i==length(names_gene)) cat(' Finalizing...',"\n")
    if(!Robust) {
      Model=Model
    }else{
      m = fitDist(model_dataFrame$gene,type = "realplus")
      Model = m$family[1]
    }
    
    if("sigma" %in% m$parameters){
      m1 <- gamlss(gene ~ gmrf(Graphid,
                               precision=as.matrix(PrecWeightGraphAdj)),
                   ~   gmrf(Graphid,
                            precision=as.matrix(PrecWeightGraphAdj)),
                   data=model_dataFrame,
                   family = Model)
      
      # Get location Signal to Noise ratio
      sError = coef(m1$mu.coefSmo[[1]])[[1]] %>% exp
      sSpat = coef(m1$mu.coefSmo[[1]])[[2]] %>% exp
      
      SNR[,names_gene[i]] = sSpat/sError
      
      #  Get Scale Signal to Noise ratio
      sError = coef(m1$sigma.coefSmo[[1]])[[1]] %>% exp
      sSpat = coef(m1$sigma.coefSmo[[1]])[[2]] %>% exp
      
      SNRScale[,names_gene[i]] = sSpat/sError
      
    }else{
      
      m1 <- gamlss(gene ~ gmrf(Graphid,
                               precision=as.matrix(PrecWeightGraphAdj)),
                   data=model_dataFrame,
                   family = Model)
      
      # Get location Signal to Noise ratio
      sError = coef(m1$mu.coefSmo[[1]])[[1]] %>% exp
      sSpat = coef(m1$mu.coefSmo[[1]])[[2]] %>% exp
      
      SNR[,names_gene[i]] = sSpat/sError
      
    }
    
    
    # Get smoothed/predicted gene expression
    in_data_AA.fitted[, names_gene[i]] = fitted(m1)
    # Get tree effects
    treeEffect[,names_gene[i]]  = getSmo(m1)$beta
    
    list(Gene= names_gene[i],
         SNR=SNR[,names_gene[i]],
         SNRScale=SNRScale[,names_gene[i]],
         fitted= in_data_AA.fitted[, names_gene[i]],
         treeEffect = treeEffect[,names_gene[i]]
    )
  }
  
  
  for (i in seq_len(length(names_gene))) {
    if(is.null(ko_param[[i]]$SNR)) next
    # Get smoothed/predicted gene expression
    SNR[,names_gene[i]] = ko_param[[i]]$SNR
    SNRScale[,names_gene[i]] = ko_param[[i]]$SNRScale
    # Get smoothed/predicted gene expression
    in_data_AA.fitted[, names_gene[i]] = ko_param[[i]]$fitted
    # Get tree effects
    treeEffect[,names_gene[i]]  = ko_param[[i]]$treeEffect
  }
  
  parallel::stopCluster(cl) 
  return(list(SNR = SNR,SNRScale = SNRScale, Fitted =in_data_AA.fitted, TreeEffect = treeEffect ))
}


GetTreeVariableGenesDynamics <- function(mst,
                                         ExprsData,
                                         ClusterCol,
                                         TemporalCol,
                                         useWeight=TRUE, 
                                         Robust = TRUE,
                                         Model="GA",
                                         rho_tree = 0.9,
                                         rho_temp = 0.7,
                                         nCores=1){
  
  
  Idclust = sort(unique(ExprsData[,ClusterCol]))
  nclust = length(Idclust)
  nNodes = vcount(mst)
  
  if(nclust!=nNodes) stop("Nodes on data and graph are different")
  if(useWeight){  WeightGraphAdj = as_adjacency_matrix(mst,attr = "weight")
  }else{  WeightGraphAdj = as_adjacency_matrix(mst)}
  
  PrecWeightGraphAdj = Diagonal(nrow(WeightGraphAdj),x=rowSums(WeightGraphAdj!=0)) - rho_tree * WeightGraphAdj
  
  ## Temporal
  
  idtime = sort(unique(ExprsData[,TemporalCol]))
  CovTemporal = Q.AR1(length(idtime),1,rho_temp,vcov = T)
  
  TreeTemporaCov = kronecker(solve(PrecWeightGraphAdj),CovTemporal)
  TreeTemporaPrecision = solve(TreeTemporaCov) 
  
  # Tree & Temporal
  ExprsData[,"TreeTemp"] = (paste0(ExprsData[,ClusterCol],ExprsData[,TemporalCol]))
  
  idTreeTime = expand.grid(Idclust,idtime)
  idTreeTime = paste0(idTreeTime$Var1,idTreeTime$Var2)
  
  # Since the mst_graph_fromPriorKwldge was sorted according to in_data_grp_crossTab$cluster_A
  # We can link it to the original data to get graph info into the data
  id = data.frame(cluster = idTreeTime, 
                  Graphid = seq_len(length(idTreeTime)))
  
  missId  = id$Graphid[! id$cluster %in% unique(ExprsData[,"TreeTemp"])]
  id      = id[-missId,]
  id[,"TreeTemp"] = id[,"cluster"]
  id = id %>% select(-cluster)
  
  TreeTemporaPrecision = TreeTemporaPrecision[-missId,-missId]
  TreeTemporaPrecision  = forceSymmetric(TreeTemporaPrecision)
  # Update input data
  
  in_data_AA = ExprsData  %>% left_join(id, by = "TreeTemp")
  
  # Calculate gene signal to noise ratio
  
  names_gene = in_data_AA %>% dplyr::select(-all_of(ClusterCol)) %>%
    dplyr::select(-all_of(TemporalCol))  %>% 
    dplyr::select(-Graphid) %>% 
    dplyr::select(-TreeTemp) %>% 
    colnames()
  
  # Configure to hold predicted gene values
  in_data_AA.fitted =in_data_AA.fitted2= in_data_AA[,names_gene,drop=F]
  if(is.na(sum(in_data_AA.fitted))) warning("Expression contains missing data")
  
  # Configured to hold Signal-to-noise ratio
  SNR = SNR2 = matrix(NA,ncol = length(names_gene),nrow=2)
  colnames(SNR) = names_gene
  colnames(SNR2) = names_gene
  
  #Log-like 
  BIC = BIC2 = matrix(NA,ncol = length(names_gene))
  colnames(BIC) = names_gene
  colnames(BIC2) = names_gene
  # Configure to hold tree-effect
  treeEffect =treeEffect2 = matrix(NaN,nrow =nrow(TreeTemporaPrecision) ,ncol = length(names_gene))
  colnames(treeEffect) = names_gene
  colnames(treeEffect2) = names_gene
  rownames(treeEffect) = id$Graphid
  rownames(treeEffect2) = id$Graphid
  # Compute density to derive nodes weights
  
  # weight = in_data_AA %>% group_by_at(ClusterCol) %>% dplyr::summarise(n=n()) %>%
  #          group_by_at(ClusterCol) %>%
  #         dplyr::reframe(dense = n/nrow(in_data_AA),weight = 1/dense) %>% select(-dense)
  # weight$weight = weight$weight/sum(weight$weight)
  #weight = 1/dense
  #weight = weight/sum(weight)
  #  in_data_AA = in_data_AA%>% left_join(weight,by=ClusterCol)
  # Parallel computing
  
  
  cl <- makeCluster(nCores, outfile="")
  registerDoParallel(cl)
  
  
  ko_param <- foreach(i=seq_len(length(names_gene)),.errorhandling = "pass",
                      .packages = c("doParallel",
                                    "Rfast2",
                                    "tidyverse",
                                    "gamlss.spatial",
                                    "gamlss",
                                    "Matrix"
                      )
  ) %dopar% {
    
    cat("Modeling gene", (i/length(names_gene))*100,"%","\n" )
    
    ## Implement Zero inflation here
    
    model_dataFrame = in_data_AA %>% 
      dplyr::select(names_gene[i], Graphid) %>% 
      dplyr::rename(gene =names_gene[i]) %>% dplyr::mutate(geneBin= as.numeric(gene>0) )
    
    model_dataFrameNonBinary = model_dataFrame %>% filter(geneBin==1)
    missId = which(!sort(unique(model_dataFrame$Graphid)) %in% unique(model_dataFrameNonBinary$Graphid))
    # id required to capture thevtree effect below
    auxid=as.numeric(as.character(sort(unique(model_dataFrameNonBinary$Graphid))))
    
    if(length(missId)>0) {
      TreeTemporaPrecision_count = TreeTemporaPrecision
      TreeTemporaPrecision_count = TreeTemporaPrecision_count[-missId,-missId]
      TreeTemporaPrecision_count = forceSymmetric(TreeTemporaPrecision_count)
    }
    
    if (sum(is.na(model_dataFrame$gene)) > 0 | sum(is.nan(model_dataFrame$gene))>0 ){
      warning("NA in Gene",names_gene[i], ". Hence will be neglected") 
      SNR[,names_gene[i]] = NA
      in_data_AA.fitted[, names_gene[i]] = NA
      treeEffect[,names_gene[i]]  = NA
      next
    }
    
    model_dataFrame$Graphid = as.factor(model_dataFrame$Graphid)
    model_dataFrameNonBinary$Graphid=as.factor(model_dataFrameNonBinary$Graphid)
    
    #Check if zero inlated
    TruZero = 1-mean(model_dataFrame$geneBin)
    
    if(TruZero>0){
      if (i==length(names_gene)) cat(' Finalizing...',"\n")
      if(!Robust) {
        Model=Model
      }else{
        # m = fitDist(model_dataFrameNonBinary$gene,type = "realplus")
        m = fitDist(model_dataFrame$gene[model_dataFrame$gene>0],type = "realplus")
        Model = m$family[1]
      }
      
      
      #m1 <- gamlss(gene ~ gmrf(Graphid,
      #                         precision=as.matrix(TreeTemporaPrecision_count)),
      #             data=model_dataFrameNonBinary,
      #             family = Model)
      
      model_dataFrame$weight = 0
      model_dataFrame$weight[model_dataFrame$geneBin==0] = (1/sum(model_dataFrame$geneBin==0))*1e+3
      model_dataFrame$weight[model_dataFrame$geneBin==1] = 1/sum(model_dataFrame$geneBin==1)*1e+3
      model_dataFrame$gene = model_dataFrame$gene+1e-2
      model_dataFrameNonBinary$weight =  1/sum(model_dataFrame$geneBin==1)*1e+3
      
      #  m2 <- gamlss(geneBin ~ gmrf(Graphid,
      #                           precision=as.matrix(TreeTemporaPrecision)),
      #               data=model_dataFrame,weights = weight,
      #               family = "BI")
      
      
      
      m1 <-     gamlss(gene ~ gmrf(Graphid,
                                   precision=as.matrix(TreeTemporaPrecision_count)),
                       data=model_dataFrameNonBinary,#weights = weight,
                       family = Model)
      
      # Get Signal to Noise ratio
      sError = coef(m1$mu.coefSmo[[1]])[[1]] %>% exp
      sSpat = coef(m1$mu.coefSmo[[1]])[[2]] %>% exp
      
      SNR[1,names_gene[i]] = sSpat
      SNR[2,names_gene[i]] = sError
      
      BIC[,names_gene[i]] =  m1$sbc 
      
      # sError = coef(m2$mu.coefSmo[[1]])[[1]] %>% exp
      # sSpat = coef(m2$mu.coefSmo[[1]])[[2]] %>% exp
      
      # SNR2[1,names_gene[i]] = sSpat/sError
      #  BIC2[,names_gene[i]] =  m1$sbc 
      #-gen.likelihood(m1)()
      # if(i%%500 ==0) write.csv(SNR[,1:i] ,file=paste0(getwd(),"/SNR.txt"))
      
      # Get smoothed/predicted gene expression
      #in_data_AA.fitted[model_dataFrame$geneBin==1, names_gene[i]] = fitted(m1)
      #  in_data_AA.fitted[, names_gene[i]] = fitted(m1)
      #  in_data_AA.fitted2[, names_gene[i]] = fitted(m2)
      # Get tree effects
      
      #treeEffect[as.character(auxid),names_gene[i]]  = getSmo(m1)$beta
      #  treeEffect[ ,names_gene[i]]  = getSmo(m1)$beta
      #  treeEffect2[,names_gene[i]]  = getSmo(m2)$beta
      
    }else{
      model_dataFrame$gene = model_dataFrame$gene+1e-3
      
      if(!Robust) {
        Model=Model
      }else{
        m = fitDist(model_dataFrame$gene,type = "realplus")
        Model = m$family[1]
      }
      
      
      m1 <- gamlss(gene ~ gmrf(Graphid,
                               precision=as.matrix(TreeTemporaPrecision)),
                   data=model_dataFrame,
                   family = Model)
      
      # Get Singnal to Noise ratio
      sError = coef(m1$mu.coefSmo[[1]])[[1]] %>% exp
      sSpat = coef(m1$mu.coefSmo[[1]])[[2]] %>% exp
      
      SNR[,names_gene[i]] = sSpat/sError
      BIC[,names_gene[i]] =  m1$sbc 
      
      
      
      # Get smoothed/predicted gene expression
      # in_data_AA.fitted[, names_gene[i]] = fitted(m1)
      # Get tree effects
      
      # treeEffect[,names_gene[i]]  = getSmo(m1)$beta
    }
    
    
    list(Gene= names_gene[i],
         SNR=SNR[,names_gene[i]],SNR2=SNR2[,names_gene[i]],
         BIC = BIC[,names_gene[i]],BIC2 = BIC2[,names_gene[i]]
         #  fitted= in_data_AA.fitted[, names_gene[i]],fitted2= in_data_AA.fitted2[, names_gene[i]],
         #  treeEffect = treeEffect[,names_gene[i]],treeEffect2 = treeEffect2[,names_gene[i]]
    )
  }
  
  
  for (i in seq_len(length(names_gene))) {
    if(is.null(ko_param[[i]]$SNR)) next
    # Get smoothed/predicted gene expression
    SNR[,names_gene[i]] = ko_param[[i]]$SNR
    BIC[,names_gene[i]] = ko_param[[i]]$BIC
    SNR2[,names_gene[i]] = ko_param[[i]]$SNR2
    BIC2[,names_gene[i]] = ko_param[[i]]$BIC2
    # Get smoothed/predicted gene expression
    #in_data_AA.fitted[, names_gene[i]] = ko_param[[i]]$fitted
    # in_data_AA.fitted2[, names_gene[i]] = ko_param[[i]]$fitted2
    # Get tree effects
    # treeEffect[,names_gene[i]]  = ko_param[[i]]$treeEffect
    # treeEffect2[,names_gene[i]]  = ko_param[[i]]$treeEffect2
  }
  
  parallel::stopCluster(cl) 
  return(list(SNR = SNR,SNR2 = SNR2,
              BIC=BIC, BIC2=BIC2))
  # Fitted =in_data_AA.fitted, Fitted2 =in_data_AA.fitted2,
  # TreeEffect = treeEffect,TreeEffect2 = treeEffect2 ))
}



GetTreeVariableGenesDynamicsMisHandler <- function(mst,
                                                   ExprsData,
                                                   ClusterCol,
                                                   TemporalCol,
                                                   useWeight=TRUE, 
                                                   Robust = TRUE,
                                                   Model="GA",
                                                   rho_tree = 0.9,
                                                   rho_temp = 0.7,
                                                   nCores=1){
  
  
  Idclust = sort(unique(ExprsData[,ClusterCol]))
  nclust = length(Idclust)
  nNodes = vcount(mst)
  
  if(nclust!=nNodes) stop("Nodes on data and graph are different")
  if(useWeight){  WeightGraphAdj = as_adjacency_matrix(mst,attr = "weight")
  }else{  WeightGraphAdj = as_adjacency_matrix(mst)}
  
  rho_tree = 0.9
  PrecWeightGraphAdj = Diagonal(nrow(WeightGraphAdj),x=rowSums(WeightGraphAdj!=0)) - rho_tree * WeightGraphAdj
  
  ## Temporal
  rho_temp = 0.7
  idtime = sort(unique(ExprsData[,TemporalCol]))
  CovTemporal = Q.AR1(length(idtime),1,rho_temp,vcov = T)
  
  TreeTemporaCov = kronecker(solve(PrecWeightGraphAdj),CovTemporal)
  TreeTemporaPrecision = solve(TreeTemporaCov) 
  
  # Tree & Temporal
  ExprsData[,"TreeTemp"] = (paste0(ExprsData[,ClusterCol],ExprsData[,TemporalCol]))
  
  idTreeTime = expand.grid(Idclust,idtime)
  idTreeTime = paste0(idTreeTime$Var1,idTreeTime$Var2)
  
  # Since the mst_graph_fromPriorKwldge was sorted according to in_data_grp_crossTab$cluster_A
  # We can link it to the original data to get graph info into the data
  id = data.frame(cluster = idTreeTime, 
                  Graphid = seq_len(length(idTreeTime)))
  
  missId  = id$Graphid[! id$cluster %in% unique(ExprsData[,"TreeTemp"])]
  id      = id[-missId,]
  id[,"TreeTemp"] = id[,"cluster"]
  id = id %>% select(-cluster)
  
  TreeTemporaPrecision = TreeTemporaPrecision[-missId,-missId]
  TreeTemporaPrecision  = forceSymmetric(TreeTemporaPrecision)
  # Update input data
  
  in_data_AA = ExprsData  %>% left_join(id, by = "TreeTemp")
  
  # Calculate gene signal to noise ratio
  
  names_gene = in_data_AA %>% dplyr::select(-all_of(ClusterCol)) %>%
    dplyr::select(-all_of(TemporalCol))  %>% 
    dplyr::select(-Graphid) %>% 
    dplyr::select(-TreeTemp) %>% 
    colnames()
  
  # Configure to hold predicted gene values
  in_data_AA.fitted = in_data_AA[,names_gene]
  if(is.na(sum(in_data_AA.fitted))) warning("Expression contains missing data")
  
  # Configured to hold Signal-to-noise ratio
  SNR = matrix(NA,ncol = length(names_gene))
  colnames(SNR) = names_gene
  
  #Log-like 
  BIC = matrix(NA,ncol = length(names_gene))
  colnames(BIC) = names_gene
  # Configure to hold tree-effect
  treeEffect = matrix(0,nrow =nrow(TreeTemporaPrecision) ,ncol = length(names_gene))
  colnames(treeEffect) = names_gene
  rownames(treeEffect) = id$Graphid
  # Compute density to derive nodes weights
  
  # weight = in_data_AA %>% group_by_at(ClusterCol) %>% dplyr::summarise(n=n()) %>%
  #          group_by_at(ClusterCol) %>%
  #         dplyr::reframe(dense = n/nrow(in_data_AA),weight = 1/dense) %>% select(-dense)
  # weight$weight = weight$weight/sum(weight$weight)
  #weight = 1/dense
  #weight = weight/sum(weight)
  #  in_data_AA = in_data_AA%>% left_join(weight,by=ClusterCol)
  # Parallel computing
  
  
  cl <- makeCluster(nCores, outfile="")
  registerDoParallel(cl)
  
  
  ko_param <- foreach(i=seq_len(length(names_gene)),.errorhandling = "pass",
                      .packages = c("doParallel",
                                    "Rfast2",
                                    "tidyverse",
                                    "gamlss.spatial",
                                    "gamlss"
                      )
  ) %dopar% {
    
    cat("Modeling gene", (i/length(names_gene))*100,"%","\n" )
    
    
    model_dataFrame = in_data_AA %>% 
      dplyr::select(names_gene[i], Graphid) %>% 
      dplyr::rename(gene =names_gene[i]) %>% 
      dplyr::mutate(gene= gene+1e-02)
    
    if (sum(is.na(model_dataFrame$gene)) > 0 | sum(is.nan(model_dataFrame$gene))>0 ){
      warning("NA in Gene",names_gene[i], ". Hence will be neglected") 
      SNR[,names_gene[i]] = NA
      in_data_AA.fitted[, names_gene[i]] = NA
      treeEffect[,names_gene[i]]  = NA
      next
    }
    
    model_dataFrame$Graphid = as.factor(model_dataFrame$Graphid)
    
    
    if (i==length(names_gene)) cat(' Finalizing...',"\n")
    if(!Robust) {
      Model=Model
    }else{
      m = fitDist(model_dataFrame$gene,type = "realplus")
      Model = m$family[1]
    }
    
    
    for(k in 1: length(names(m$fits))){
      
      Model =  names(m$fits)[k]
      m1 =  tryCatch(
        expr = {
          m1 <- gamlss(gene ~ gmrf(Graphid,
                                   precision=as.matrix(TreeTemporaPrecision)),
                       data=model_dataFrame,
                       family = Model)
        },
        error = function(e){ 
          m1 = NULL
        },
        warning = function(w){
          
          m1 = m1
          
        },
        finally = {
          #m1
        }
      )
      
      if(!is.null(m1)){
        break
      }
    }
    # Get Singnal to Noise ratio
    sError = coef(m1$mu.coefSmo[[1]])[[1]] %>% exp
    sSpat = coef(m1$mu.coefSmo[[1]])[[2]] %>% exp
    
    SNR[,names_gene[i]] = sSpat/sError
    BIC[,names_gene[i]] =  m1$sbc #-gen.likelihood(m1)()
    # if(i%%500 ==0) write.csv(SNR[,1:i] ,file=paste0(getwd(),"/SNR.txt"))
    
    # Get smoothed/predicted gene expression
    in_data_AA.fitted[, names_gene[i]] = fitted(m1)
    # Get tree effects
    treeEffect[,names_gene[i]]  = getSmo(m1)$beta
    
    list(Gene= names_gene[i],
         SNR=SNR[,names_gene[i]],
         BIC = BIC[,names_gene[i]],
         fitted= in_data_AA.fitted[, names_gene[i]],
         treeEffect = treeEffect[,names_gene[i]]
    )
  }
  
  
  for (i in seq_len(length(names_gene))) {
    if(is.null(ko_param[[i]]$SNR)) next
    # Get smoothed/predicted gene expression
    SNR[,names_gene[i]] = ko_param[[i]]$SNR
    BIC[,names_gene[i]] = ko_param[[i]]$BIC
    # Get smoothed/predicted gene expression
    in_data_AA.fitted[, names_gene[i]] = ko_param[[i]]$fitted
    # Get tree effects
    treeEffect[,names_gene[i]]  = ko_param[[i]]$treeEffect
  }
  
  parallel::stopCluster(cl) 
  return(list(SNR = SNR,BIC=BIC, Fitted =in_data_AA.fitted, TreeEffect = treeEffect ))
}

# Dual Mesh

book.mesh.dual <- function(mesh) {
  if (mesh$manifold=='R2') {
    ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
      colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
    library(parallel)
    pls <- mclapply(1:mesh$n, function(i) {
      p <- unique(Reduce('rbind', lapply(1:3, function(k) {
        j <- which(mesh$graph$tv[,k]==i)
        if (length(j)>0) 
          return(rbind(ce[j, , drop=FALSE],
                       cbind(mesh$loc[mesh$graph$tv[j, k], 1] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 1], 
                             mesh$loc[mesh$graph$tv[j, k], 2] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 2])/2))
        else return(ce[j, , drop=FALSE])
      })))
      j1 <- which(mesh$segm$bnd$idx[,1]==i)
      j2 <- which(mesh$segm$bnd$idx[,2]==i)
      if ((length(j1)>0) | (length(j2)>0)) {
        p <- unique(rbind(mesh$loc[i, 1:2], p,
                          mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2]/2, 
                          mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2]/2))
        yy <- p[,2]-mean(p[,2])/2-mesh$loc[i, 2]/2
        xx <- p[,1]-mean(p[,1])/2-mesh$loc[i, 1]/2
      }
      else {
        yy <- p[,2]-mesh$loc[i, 2]
        xx <- p[,1]-mesh$loc[i, 1]
      }
      Polygon(p[order(atan2(yy,xx)), ])
    })
    return(SpatialPolygons(lapply(1:mesh$n, function(i)
      Polygons(list(pls[[i]]), i))))
  }
  else stop("It only works for R2!")
}

##################################
# Communication specific functions
#################################

# Function to project P_hat onto the space of ordered probabilities
project_stochastic_order <- function(P_hat, Sigma) {
  p <- length(P_hat)  # Number of dimensions
  
  # Compute the inverse of the covariance matrix
  Sigma_inv <- solve(Sigma)
  
  # Define the quadratic programming inputs
  Dmat <- 2 * Sigma_inv  # Quadratic term (must be positive definite)
  dvec <- 2 * Sigma_inv %*% P_hat  # Linear term
  
  # Constraint matrix for order restrictions P_1 <= P_2 <= ... <= P_p
  Amat <- matrix(0, nrow = p , ncol = p)
  for (i in 1:p) {
    # Amat[i, i] <- -1
    Amat[i, i] <- 1
  }
  
  # Constraint vector (enforcing increasing order)
  bvec <- rep(0.5, p)
  
  # Solve the quadratic programming problem
  result <- solve.QP(Dmat, dvec, t(Amat), bvec )
  
  return(result$solution)  # Projected vector P_bar
}

GetTreeVariableGenesDynamics.INLA <- function(
    mst,
    ExprsData,
    ClusterCol,
    TemporalCol,
    ConfoundFrame=NULL,
    useWeight = FALSE, 
    Model="gaussian",
    Transform = FALSE,
    Fun = function(x)x,
    rho_tree = 0.7,
    rho_temp = 0.9,
    IncZero= TRUE,
    DownSample = FALSE,
    nCores =1
){
  
  Idclust = as.vector(V(mst))#sort(unique(ExprsData[,ClusterCol]))
  nclust = length(Idclust)
  nNodes = vcount(mst)
  
  if(nclust!=nNodes) stop("Nodes on data and graph are different")
  if(useWeight){  WeightGraphAdj = as_adjacency_matrix(mst,attr = "weight")
  }else{  WeightGraphAdj = as_adjacency_matrix(mst)}
  
  PrecWeightGraphAdj = Diagonal(nrow(WeightGraphAdj),x=rowSums(WeightGraphAdj!=0)) - rho_tree * WeightGraphAdj
  
  ## Temporal
  
  idtime = sort(unique(ExprsData[,TemporalCol]))
  if(length(idtime)==1) {CovTemporal =matrix(1)
  }else{  
    CovTemporal = Q.AR1(length(idtime),1,rho_temp,vcov = T)
  }
  
  TreeTemporaCov = kronecker(solve(PrecWeightGraphAdj),CovTemporal)
  TreeTemporaPrecision = solve(TreeTemporaCov) 
  
  # Tree & Temporal
  
  ExprsData[,"TreeTemp"] = (paste0(ExprsData[,ClusterCol],ExprsData[,TemporalCol]))
  
  idTreeTime = expand.grid(Idclust,idtime)
  idTreeTime = paste0(idTreeTime$Var1,idTreeTime$Var2)
  
  # Since the mst_graph_fromPriorKwldge was sorted according to in_data_grp_crossTab$cluster_A
  # We can link it to the original data to get graph info into the data
  id = data.frame(cluster = idTreeTime, 
                  Graphid = seq_len(length(idTreeTime)))
  
  missId  = id$Graphid[! id$cluster %in% unique(ExprsData[,"TreeTemp"])]
  id$weight  = 1
  id$weight[id$Graphid %in% missId]  = 0
  id[,"TreeTemp"] = id[,"cluster"]
  id = id %>% dplyr::select(-cluster)
  
  #TreeTemporaPrecision = TreeTemporaPrecision[-missId,-missId]
  TreeTemporaPrecision  = forceSymmetric(TreeTemporaPrecision)
  # Update input data
  
  in_data_AA = ExprsData  %>% bind_cols(ConfoundFrame=ConfoundFrame)%>%  left_join(id, by = "TreeTemp")
  in_data_AA_aux = matrix(0,nrow = length(missId),ncol=ncol(in_data_AA))%>% as.data.frame()
  colnames(in_data_AA_aux) = colnames(in_data_AA)
  in_data_AA_aux$Graphid=id$Graphid[id$weight==0]
  in_data_AA_aux$TreeTemp=id$TreeTemp[id$weight==0]
  in_data_AA = bind_rows(in_data_AA,in_data_AA_aux)
  
  # Calculate gene signal to noise ratio
  ConfName = colnames(ConfoundFrame)
  
  names_gene = in_data_AA %>% dplyr::select(-all_of(ClusterCol)) %>% 
    dplyr::select(-all_of(ConfName)) %>%
    dplyr::select(-all_of(TemporalCol))  %>% 
    dplyr::select(-Graphid,-TreeTemp,-weight) %>% 
    colnames()
  
  SNR =matrix(NA,ncol = length(names_gene),nrow=3)
  colnames(SNR) = names_gene
  
  
  #Log-like 
  BIC = matrix(NA,ncol = length(names_gene),nrow=3)
  colnames(BIC) = names_gene
  
  
  #Log-like 
  Fam = matrix(NA,ncol = length(names_gene),nrow=2) %>%as.data.frame()
  colnames(Fam) = names_gene
  colnames(Fam) = names_gene
  
  # Tree Effect
  treeEffect = matrix(NaN,nrow =nrow(TreeTemporaPrecision) ,ncol = length(names_gene))
  colnames(treeEffect) = names_gene
  
  
  cl <- makeCluster(nCores, outfile="")
  registerDoParallel(cl)
  
  ko_param <- foreach(i=seq_len(length(names_gene)),.errorhandling = "pass",
                      .packages = c("doParallel",
                                    "Rfast2",
                                    "tidyverse",
                                    "gamlss.spatial",
                                    "gamlss",
                                    "Matrix",
                                    "INLA"
                      )
  ) %dopar% {
    
    cat("Modeling gene", (i/length(names_gene))*100,"%","\n" )
    
    ## Implement Zero inflation here
    
    model_dataFrame = in_data_AA %>% 
      dplyr::select(names_gene[i], Graphid,weight,all_of(ConfName)) %>%
      dplyr::rename(gene =names_gene[i]) %>% dplyr::mutate(geneBin= as.numeric(gene>0) )
    # model_dataFrame$Graphid = as.factor(model_dataFrame$Graphid)
    model_dataFrame$weight[model_dataFrame$gene==0] = model_dataFrame$weight[model_dataFrame$gene==0] * as.numeric(IncZero)
    
    # Downsample 
    Hcluster <- function(DataTocluster,thresholdGini=0.2,k,ClusterName){
      # Goal: Function to cluster a data using hierarchical clustering
      # INPUT: 
      # DataTocluster - > data matrix to cluster
      # thresholdGini - > Gini hyperparameter
      # OUTPUT: Cluste labled class
      
      cluster  <- genie::hclust2(objects=as.matrix(DataTocluster), thresholdGini=thresholdGini)
      clust    <- cutree(cluster, k = k)
      
      clust = as.matrix(clust) 
      colnames(clust) = ClusterName
      ClusteredData = cbind(clust,DataTocluster)
      
      return(ClusteredData)
    }
    if(DownSample){
      set.seed(23456)
      Aux_Downsamp = model_dataFrame %>% group_by(Graphid) %>% summarise(nn=n())
      Downsampled_model_dataFrame = NULL
      cat("Downsampling by node","\n")
      for (k in Aux_Downsamp$Graphid) {
        aux = model_dataFrame %>%filter(Graphid==k)
        aux$Graphid = as.numeric(aux$Graphid)
        if(nrow(aux)>100){
          # in_data = Hcluster(DataTocluster = (aux%>% select(gene)),
          #                    thresholdGini = 0.2,
          #                    k=50,ClusterName = "clusterr")
          #aux_down =   in_data %>% bind_cols(aux%>% select(-gene,-Graphid))%>% group_by(clusterr) %>% summarise_all(mean,na.rm=T) %>% select(gene,ConfName)  %>% mutate(Graphid = k)
          #aux_down$Graphid = as.numeric(aux_down$Graphid )
          aux_down = aux[sample(1:nrow(aux),100),]
          Downsampled_model_dataFrame = Downsampled_model_dataFrame%>% bind_rows(aux_down)
        }else{
          
          # Downsampled_model_dataFrame = Downsampled_model_dataFrame%>% bind_rows(aux%>% select(-weight, -geneBin))
          Downsampled_model_dataFrame = Downsampled_model_dataFrame%>% bind_rows(aux)
        }
        
      }
      model_dataFrame = Downsampled_model_dataFrame
      # model_dataFrame$Graphid = as.factor(model_dataFrame$Graphid)
    }
    if(Transform){
      
      model_dataFrame$gene = Fun(model_dataFrame$gene)
    }
    
    
    # This is where cancer model is different
    for (kk in seq_len(length(ConfName))) {
      aux = fastDummies::dummy_columns( model_dataFrame[,ConfName[kk]], remove_most_frequent_dummy = T)
      colnames(aux) = paste0(ConfName[kk],colnames(aux))
      model_dataFrame = model_dataFrame %>% bind_cols(aux)
    }
    
    ConfName_a = grep(".data_",colnames(model_dataFrame))
    ConfName_aux = colnames(model_dataFrame[,ConfName_a])
    
    
    formula = Y ~  1  + 
      f(s, model = "generic0", Cmatrix = TreeTemporaPrecision, constr = TRUE,initial = c(0))+
      f(s2, model = "generic0", Cmatrix = diag(nrow(TreeTemporaPrecision)), constr = TRUE,initial = c(0))
    
    
    
    if( length(ConfName_a)>0){
      for(kk in 1:length(ConfName_a)){
        f = as.formula(paste0("~.+", ConfName_aux[kk],sep=""))
        formula = update(formula,   f)
      }
    }
    model_dataFrame$Y = model_dataFrame$gene
    model_dataFrame$s = model_dataFrame$Graphid
    model_dataFrame$s2 =model_dataFrame$Graphid
    
    m1 <- inla(formula,
               data =model_dataFrame,
               family = Model,weights = weight,
               control.fixed = list(prec = 0.1), verbose=F, 
               control.compute = list(config = TRUE,dic=TRUE, cpo=TRUE, waic=TRUE))
    
    # Get Signal to Noise ratio
    sSpat = 1/m1$summary.hyperpar["Precision for s",1]
    sError = 1/m1$summary.hyperpar["Precision for s2",1]
    
    SNR[1,names_gene[i]] = sSpat
    SNR[2,names_gene[i]] = sError
    SNR[3,names_gene[i]] = sSpat/sError
    
    BIC[1,names_gene[i]] =  m1$dic$dic
    BIC[2,names_gene[i]] =  m1$waic$waic
    BIC[3,names_gene[i]] =  m1$mlik[1,1][[1]]
    
    Fam[1,names_gene[i]] = Model
    Fam[2,names_gene[i]] = Model
    
    treeEffect[,names_gene[i]]  = m1$summary.random$s$mean
    
    list(Gene= names_gene[i],
         SNR=SNR[,names_gene[i]],
         BIC = BIC[,names_gene[i]],
         Fam = Fam[,names_gene[i]],
         treeEffect = treeEffect[,names_gene[i]]
    )
  }
  
  
  for (i in seq_len(length(names_gene))) {
    if(is.null(ko_param[[i]]$SNR)) {Error = ko_param[[i]];next}
    # Get smoothed/predicted gene expression
    Error = 0
    SNR[,names_gene[i]] = ko_param[[i]]$SNR
    BIC[,names_gene[i]] = ko_param[[i]]$BIC
    Fam[,names_gene[i]] = ko_param[[i]]$Fam
    treeEffect[,names_gene[i]] = ko_param[[i]]$treeEffect
  }
  
  parallel::stopCluster(cl) 
  return(list(SNR = SNR,
              BIC=BIC,
              Fam = Fam,
              treeEffect = treeEffect,
              Error=Error)
  )
}



GetId <- function(
    mst,
    ExprsData,
    ClusterCol,
    TemporalCol,
    useWeight = FALSE
){
  
  Idclust = as.vector(V(mst))#sort(unique(ExprsData[,ClusterCol]))
  nclust = length(Idclust)
  nNodes = vcount(mst)
  
  if(nclust!=nNodes) stop("Nodes on data and graph are different")
  
  ## Temporal
  
  idtime = sort(unique(ExprsData[,TemporalCol]))
  
  # Tree & Temporal
  
  ExprsData[,"TreeTemp"] = (paste0(ExprsData[,ClusterCol],ExprsData[,TemporalCol]))
  
  idTreeTime = expand.grid(Idclust,idtime)
  idTreeTime = paste0(idTreeTime$Var1,idTreeTime$Var2)
  
  # Since the mst_graph_fromPriorKwldge was sorted according to in_data_grp_crossTab$cluster_A
  # We can link it to the original data to get graph info into the data
  id = data.frame(cluster = idTreeTime, 
                  Graphid = seq_len(length(idTreeTime)))
  
  missId  = id$Graphid[! id$cluster %in% unique(ExprsData[,"TreeTemp"])]
  id$weight  = 1
  id$weight[id$Graphid %in% missId]  = 0
  id[,"TreeTemp"] = id[,"cluster"]
  id = id %>% dplyr::select(-cluster)
  
  return(list(netwknode = id, TreeTemp_data = data.frame(TreeTemp=ExprsData[,"TreeTemp"])))
}

NetworkGet = function(metaData,datExpr){

un = unique(metaData$Slice_ID)
META=NULL
cout =0

for (i in un ) {
  
  ux = metaData[metaData$Slice_ID==i,]
  #AuxData = ExprData[metaData$Slice_ID==i,]
  aux = ux %>% select(x,y)
  
  
  coords <- data.frame(
    lon =aux$x,# 
    lat =aux$y#  
  )
  coords[,"lon"] = scale(coords[,"lon"])
  coords[,"lat"] = scale(coords[,"lat"])
  #########
  library(INLA)
  library(sf)
  mesh = inla.mesh.2d(loc.domain = coords, offset = c(1, 1), 
                      max.edge = c(2.0, 2.0), cutoff =2)
  dmesh    = book.mesh.dual(mesh=mesh)
  dmesh2 = st_as_sf(dmesh)
  coordsPoints = st_as_sf(coords, coords =c(1,2))
  aux = st_intersects(coordsPoints, dmesh2)
  aux = aux %>% unlist()
  
  coords_with_cells <- data.frame(ux, mesh_id = paste0(aux,Slice_ID=i,"^",metaData$Tier3[metaData$Slice_ID==i],"#",metaData$Leiden_neigh[metaData$Slice_ID==i]) )
  #coords_with_cells <- data.frame(ux, mesh_id = paste0(Slice_ID=i,"^",metaData$Tier3[metaData$Slice_ID==i],"#",metaData$Leiden_neigh[metaData$Slice_ID==i]) )
  #coords_with_cells <- data.frame(ux, mesh_id = paste0(metaData$Tier3[metaData$Slice_ID==i],"#",metaData$Leiden_neigh[metaData$Slice_ID==i]) )
  #rownames(coords_with_cells) = rownames(ux)
  META <- rbind(META,coords_with_cells)
  
  #AuxDATA <- rbind(AuxDATA,AuxData)
  
  cout = cout+1
  cat(cout,"\n")
}

ax_reduce_META = META
datExpr = datExpr[rownames(ax_reduce_META),]
#Metadat = META

#META$id = paste0(META$Tier1,META$mesh_id) # Downsampling by these parameters
#datExpr = AuxDATA
#un = unique(META$id)
#length(un)
#[1] 158083
############# Downsample
# 
# set.seed(1234)
# cout = 0
# reduce_META = reduce_datExpr = NULL
# for (i in un) {
#   
#   rex =  META[META$id==i,]
#   # rex_data = datExpr[META$id==i,]
#   if(nrow(rex)>20){
#     
#     id_rex = sample(1:nrow(rex),20)
#     reduce_META = rbind(reduce_META,rex[id_rex,])
#     # reduce_datExpr = rbind(reduce_datExpr,rex_data[id_rex,])
#     
#   }else{
#     reduce_META = rbind(reduce_META,rex)
#   }
#   if(cout%%500==0)cat(cout,"\n")
#   cout = cout+1
# }

#datExpr = ExprData[reduce_META$cell_IDs,]
#cl = colSums(datExpr==0)
### Plot
#Aux = data.frame(z=datExpr[,100],reduce_META) %>% 
#  filter(Slice_ID=="062221_D9_m3_2_slice_3")
#plotScatter(Aux$x,Aux$y,Aux$z,main="",size=.3,
#            legend.size = 1.5,ManualColor = F,cols = c25)
### Plot ends

#subset= c(6,7,15,16,19,20,21,22,23,26,31,32,33,37,41,42,47,49,50,53,55,56,57)



##############
#########################

#library(Seurat)
#seurat_obj <- CreateSeuratObject(counts = t(datExpr), meta.data = ax_reduce_META)

# Split object by sample
#seurat_list <- SplitObject(seurat_obj, split.by = "Slice_ID")

# Normalize each sample individually using SCTransform

#seurat_list <- lapply(seurat_list, NormalizeData, method = "LogNormalize", verbose = TRUE)

#features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)

#anchors <- FindIntegrationAnchors(
#  object.list = seurat_list,
#  anchor.features = SelectIntegrationFeatures(seurat_list, nfeatures = 3000),
#  verbose = TRUE
#)

#seurat_list <- lapply(seurat_list, function(x) {
#  x <- ScaleData(
#    object = x,
# vars.to.regress = c("Slice_ID"),
#    verbose = TRUE
#  )
#  return(x)
#})

#datExpr =ax_reduce_META= NULL

#for(i in 1:length(seurat_list)){
#  datExpr = cbind(datExpr,seurat_list[[i]]@assays$RNA$scale.data)
#  ax_reduce_META = rbind(ax_reduce_META, seurat_list[[i]]@meta.data)
#}
#datExpr = t(datExpr)
#############
datExpr_withNeigborhood = datExpr%>%as.matrix%>%as.data.frame()

#datExpr_withNeigborhood$NeighXCell = ax_reduce_META$mesh_id 

###
#aux = datExpr_withNeigborhood$NeighXCell%>%table
#nam_aux = names(aux[which(aux>500)])
#Mini_datExpr_withNeigborhood = datExpr_withNeigborhood %>%filter(NeighXCell%in%nam_aux)

#kmeans_result <- kmeans(Mini_datExpr_withNeigborhood%>%select(-NeighXCell), centers =20, nstart = 25,iter.max = 20)

#datExpr_withNeigborhood[rownames(Mini_datExpr_withNeigborhood),"NeighXCell"]= 
# paste0(datExpr_withNeigborhood[rownames(Mini_datExpr_withNeigborhood),"NeighXCell"],"##",kmeans_result$cluster)
####


####################
# Begin Modeling
###################
#library(umap)
#library(patchwork)


#################
#Modeling Day 4 (Control) to day 5 (Treatment)
#################

set.seed(12345)

#datExpr_withNeigborhood_sub = datExpr_withNeigborhood 

#umap_in_data = datExpr_withNeigborhood_sub%>%dplyr::select(-day.harvested,-NeighXCell,-NeighXCell_id)

#data_umap3d <- umap(umap_in_data, 
#                    n_neighbors = 20, min_dist = 1, 
#                    spread = 5,n_components = 5) #3D version

ax_reduce_META$NeighXCellXclust = paste0(ax_reduce_META$Leiden_neigh,"^",ax_reduce_META$Tier3)

###
aux = ax_reduce_META$NeighXCellXclust%>%table
nam_aux = names(aux[which(aux>500)])
Mini_ax_reduce_META = ax_reduce_META %>% filter(NeighXCellXclust %in% nam_aux)

A = datExpr_withNeigborhood[rownames(Mini_ax_reduce_META),]

kmeans_result <- kmeans(A, centers =15, nstart = 25,iter.max = 20)
remove(A)

ax_reduce_META[rownames(Mini_ax_reduce_META),"NeighXCellXclust"] = paste0(ax_reduce_META[rownames(Mini_ax_reduce_META),"NeighXCellXclust"],"^^",
                                                                          kmeans_result$cluster )

NeighXCell_id = as.numeric(as.factor(ax_reduce_META$NeighXCellXclust) )

#data_umap= data.frame(seurat_integrated@reductions$umap@cell.embeddings,NeighXCell_id=NeighXCell_id)

datExpr_withNeigborhood$NeighXCell_id = NeighXCell_id

Centers = datExpr_withNeigborhood%>%as.data.frame()%>% as.data.frame %>% group_by(NeighXCell_id) %>% summarise_all(mean,na.rm=T)%>%
  ungroup() %>% dplyr::select(-NeighXCell_id)

Centers2 = apply(Centers,2,scale)

return(list(Centers = Centers2,
       datExpr_withNeigborhood = datExpr_withNeigborhood,
       ax_reduce_META = ax_reduce_META)
       )
}

FindDEGs <- function(datExpr_withNeigborhood_sub,
                     ax_reduce_META,
                     mst,
                     day=1
                 ) {
  

Fun <- function(x)x+0.000
nam = colnames(datExpr_withNeigborhood_sub)


############ GetLigand-Receptor pairs
##################################
########### Get cell - cell talk#
##################################
Refined_data = datExpr_withNeigborhood_sub
DriverGenes_C = nam
library(SingleCellSignalR)
GeneEquivaalentInMouseHuman = mm2Hs
LiganReceiptorDatabaseHuman = LRdb
mouse_lr_pair = readRDS("~/Downloads/DSFMix-main/EMTVenosa/mouse_lr_pair.rds")
human_lr_pair = readRDS("~/Downloads/DSFMix-main/EMTVenosa/human_lr_pair.rds")

# From iTALK
ital_db <- iTALK:::database %>%select(Ligand.ApprovedSymbol,Receptor.ApprovedSymbol)


colnames(ital_db) = c("ligand", "receptor")
mm2Hs = GeneEquivaalentInMouseHuman
LRdb = rbind(LiganReceiptorDatabaseHuman[,c("ligand","receptor")],
             data.frame(ligand=human_lr_pair$ligand_gene_symbol,
                        receptor = human_lr_pair$receptor_gene_symbol))
LRdb = rbind(LRdb,ital_db)
LRdb = unique(LRdb)


findEqvGene = function(gene,mouse=T){
  
  load("/Users/egbonoa/Downloads/DSFMix-main/EMTVenosa/GeneEquivaalentInMouseHuman.Rdata")
  load("/Users/egbonoa/Downloads/DSFMix-main/EMTVenosa/LiganReceiptorDatabaseHuman.Rdata")
  
  Rex =NULL
  
  for (g in gene) {
    if(mouse){rex= mm2Hs[,2][which(mm2Hs[,1]%in%g)]
    
    }else{rex = mm2Hs[,1][which(mm2Hs[,2]%in%g)] }
    if(length(rex)==0)rex=NA
    Rex=c(Rex,rex[1])
  }
  
  return(Rex)
}

GetLiganRecptor <- function(data, genes,mouse=TRUE){
  
  if(mouse){
    
    TopGenes = data.frame(data_mouse = genes,data_human =  findEqvGene(gene=genes,mouse = T))
    TopGenes = na.omit(TopGenes)
    LiganInData = which(LRdb$ligand%in%TopGenes$data_human)
    LiganInData = LRdb[LiganInData,c("ligand","receptor")]
    LiganInData$IDFromData = "Ligan"
    ReceptorInData = which(LRdb$receptor %in%TopGenes$data_human)
    ReceptorInData =LRdb[ReceptorInData,c("ligand","receptor")]
    ReceptorInData$IDFromData = "Receptor"
    
    LigRecep = rbind(LiganInData,ReceptorInData) %>%as.data.frame()
    LigRecep$LiganInMouse =  findEqvGene(gene=LigRecep$ligand,mouse = F)
    LigRecep$RecepInMouse =  findEqvGene(gene=LigRecep$receptor,mouse = F)
    LigRecep = na.omit(LigRecep)
    dataNames = colnames(data)
    id_temp = which(LigRecep$LiganInMouse%in%dataNames & LigRecep$RecepInMouse%in%dataNames)
    LigRecep = LigRecep[id_temp,]
  }else{
    TopGenes = data.frame(data_human = genes)
    LiganInData = which(LRdb$ligand%in%TopGenes$data_human)
    LiganInData = LRdb[LiganInData,c("ligand","receptor")]
    LiganInData$IDFromData = "Ligan"
    ReceptorInData = which(LRdb$receptor %in%TopGenes$data_human)
    ReceptorInData =LRdb[ReceptorInData,c("ligand","receptor")]
    ReceptorInData$IDFromData ="Receptor"
    
    LigRecep = rbind(LiganInData,ReceptorInData) %>%as.data.frame()
    LigRecep$LiganInMouse =  findEqvGene(gene=LigRecep$ligand,mouse = F)
    LigRecep$RecepInMouse =  findEqvGene(gene=LigRecep$receptor,mouse = F)
    LigRecep = na.omit(LigRecep)
    dataNames = colnames(data)
    id_temp = which(LigRecep$LiganInMouse%in%dataNames & LigRecep$RecepInMouse%in%dataNames)
    LigRecep = LigRecep[id_temp,]
  }
  
  return(LigRecep)
}

Rex_axu = intersect(union(mouse_lr_pair$ligand_gene_symbol,
                          mouse_lr_pair$receptor_gene_symbol),
                    DriverGenes_C)

Rex_axu = c(which(mouse_lr_pair$ligand_gene_symbol%in%Rex_axu),which(mouse_lr_pair$receptor_gene_symbol%in%Rex_axu))
Rex_axu = unique(mouse_lr_pair[Rex_axu,])
Rex_axu = data.frame(LiganInMouse=Rex_axu$ligand_gene_symbol,RecepInMouse =Rex_axu$receptor_gene_symbol)
axs = Rex_axu$LiganInMouse%in%colnames(Refined_data)&Rex_axu$RecepInMouse%in%colnames(Refined_data)
LRInData_celltalkDB = Rex_axu[axs,]

LRInData = GetLiganRecptor(data = Refined_data, genes = c(DriverGenes_C),mouse = T)
LRInData = unique(LRInData[,c("LiganInMouse","RecepInMouse")])


LRInData = rbind(LRInData,LRInData_celltalkDB)
LRInData = unique(LRInData)

nam_sub  = union(LRInData$LiganInMouse,LRInData$RecepInMouse)
############################
############################


#day = 1

subdata = datExpr_withNeigborhood_sub %>%filter(day.harvested==day)
ctype_day1 = ax_reduce_META$Tier3[datExpr_withNeigborhood_sub$day.harvested==day]
ConfoundFrame = data.frame(ax_reduce_META$Slice_ID[datExpr_withNeigborhood_sub$day.harvested==day] %>%
                             as.numeric())
META_subb = ax_reduce_META[datExpr_withNeigborhood_sub$day.harvested==day,]


Result_Cancer_day2 = GetTreeVariableGenesDynamics.INLA(
  mst = mst,
  ExprsData = subdata[,c(nam_sub,"NeighXCell_id","day.harvested")] %>% as.data.frame(),
  ClusterCol ="NeighXCell_id",
  TemporalCol ="day.harvested",
  ConfoundFrame=NULL,
  useWeight = FALSE, 
  Model="normal",
  Transform = TRUE,
  Fun = Fun,
  rho_tree = 0.7,
  rho_temp = 0.9,
  IncZero= TRUE,
  DownSample = F,
  nCores =10
)

return(list(Result_Cancer_day= Result_Cancer_day2,
            LRInData = LRInData,
            subdata = subdata[,c(nam_sub,"NeighXCell_id","day.harvested")],
            ctype_day1=ctype_day1,
            META_subb = META_subb)
)

}


Teffect = function(Result_Cancer_day,
                                mst,
                                subdata,
                                 ctype_day1 
                                ){

###
##################
#treeEffect = Result_Cancer_day$treeEffect
treeEffect = Result_Cancer_day$treeEffect
#subdata$NeighXCell_id %>% unique()


DD = GetId(mst,
           ExprsData = subdata%>%as.data.frame(),
           ClusterCol ="NeighXCell_id" ,
           TemporalCol ="day.harvested",
           useWeight = FALSE
)

TreeTemp = DD$netwknode
TreeTemp_data = DD$TreeTemp_data

TreeTemp_effect = bind_cols(TreeTemp,treeEffect)

#ctype_day1
 
TreeTemp_data_effect = left_join(TreeTemp_data,TreeTemp_effect,by="TreeTemp")
TreeTemp_data_effect = bind_cols(subdata[,c("NeighXCell_id","day.harvested")],TreeTemp_data_effect)
TreeTemp_data_effect = bind_cols(Ctype = ctype_day1,TreeTemp_data_effect)

return(list(TreeTemp_data_effect= TreeTemp_data_effect,
            TreeTemp_data = TreeTemp_data)
      )
}

getNull  = function(Tree_effect,LRInData){

  treeEffect = Tree_effect
  uniq_ctype = treeEffect$Ctype %>% unique()

  snr = colnames(treeEffect)[-c(1:6)] # Remove non-genes
 
  LRInData = LRInData[LRInData$LiganInMouse %in% snr | LRInData$RecepInMouse %in% snr,]
  
  set.seed(1234)
  
  LR0_effect = NULL  
  
  # Get Null Random vector of ligand receptors
  
  for (k in 1:1000) {
    
    Rx = LRInData[sample(1:nrow(LRInData),1),c("LiganInMouse","RecepInMouse")]
    L_0 = Rx$LiganInMouse
    Rx = LRInData[sample(1:nrow(LRInData),1),c("LiganInMouse","RecepInMouse")]
    R_0 = Rx$RecepInMouse
    
    L0_effect = treeEffect %>%filter(Ctype == sample(uniq_ctype,1) ) %>%dplyr::select_at(L_0)
    R0_effect = treeEffect %>%filter(Ctype == sample(uniq_ctype,1)) %>%dplyr::select_at(R_0)
    
    nL = nrow(L0_effect)
    nR = nrow(R0_effect)
    nn = ceiling((nL+nR)/4)
    
    aux_exp = cbind(L0 = abs(L0_effect[sample(1:nL,nn,replace = T),1]),R0 = abs(R0_effect[sample(1:nR,nn,replace = T),1])) %>%as.matrix()
    LR0_effect = rbind(LR0_effect,aux_exp)
    if(k%%100==0)cat(k,"\n")
  }
  
  
  
  
DN = NULL



for (k in 1:1000) {
  
  
  # Compute P-hat
  X = LR0_effect%>%as.matrix()
  
  
  n11 = 1000
  n22 = 1000
  
  Xij_aux = X[sample(1:nrow(X),n11,replace = T),,drop=F]
  Xji_aux = X[sample(1:nrow(X),n22,replace = T),,drop=F]
  
  p1_hat_rex_1 = (expand_grid(Xij_aux[ ,1],Xji_aux[ ,1]))
  p1_hat_rex_2 = (expand_grid(Xij_aux[ ,2],Xji_aux[ ,2]))
  
  p1_hat_mean = c(mean(p1_hat_rex_1[,1]<=p1_hat_rex_1[,2]),
                  mean(p1_hat_rex_2[,1]<=p1_hat_rex_2[,2]))
  
  n1= nrow(Xij_aux)
  n2 =nrow(Xji_aux)
  n =n1+n2
  p1_Hat = NULL
  for (l in 1:30) {
    n11 = 30 # 
    n22 = 30 # 
    id_1 = sample(1:n1,n11,replace = T)
    id_2 = sample(1:n2,n22,replace = T)
    p1_hat_rex_1 = (expand_grid(Xij_aux[id_1,1],Xji_aux[id_2,1]))
    p1_hat_rex_2 = (expand_grid(Xij_aux[id_1,2],Xji_aux[id_2,2]))
    
    p1_hat = c(mean(p1_hat_rex_1[,1]<=p1_hat_rex_1[,2]),
               mean(p1_hat_rex_2[,1]<=p1_hat_rex_2[,2]))
    
    p1_Hat = rbind(p1_Hat,p1_hat)
    
  }
  
  Sigma0 = cov(p1_Hat)
  
  # Estimated probabilities
  
  P_hat <- p1_hat_mean 
  Sigma <- Sigma0
  
  aux_result <- try({
    
    P_bar <- project_stochastic_order(P_hat, Sigma)
    
  }, silent = TRUE) 
  
  if (inherits(aux_result, "try-error")) {
    
    next
  }
  P_bar <- aux_result
  #### Calculate Dn statistic
  P0 = matrix(c(0.5,0.5))
  Dn = drop(n*((t(P_hat-P0)%*%Sigma%*%(P_hat-P0)) - (t(P_hat-P_bar)%*%Sigma%*%(P_hat-P_bar))))
  DN = c(DN,Dn)
  if(k%%100==0)cat(k,"\n")
}
return(list(DN=DN, LR0_effect = LR0_effect))
}

#####################
# Estimate communication significance
#####################



ComputeCCS = function(TreeTemp_data_effect,
                      LRInData,
                      CTn,
                      DN,
                      LR0_effect){

treeEffect = TreeTemp_data_effect
uniq_ctype = as.character(unique(TreeTemp_data_effect$Ctype))
uniq_ctype = uniq_ctype[uniq_ctype %in% CTn]


CT = NULL

for(kprim1 in 1:length(uniq_ctype)){
  
  for (kprim2 in 1:length(uniq_ctype)) {
    CT = c(CT,paste0(uniq_ctype[kprim1],"^",uniq_ctype[kprim2]))
  }
}

ResultMatrix <- ResultMatrixPval <- array(NaN, dim = c(nrow(LRInData),length(uniq_ctype)*length(uniq_ctype)) )

rownames(ResultMatrix) = paste0(LRInData$LiganInMouse,"^",LRInData$RecepInMouse)
rownames(ResultMatrixPval) = paste0(LRInData$LiganInMouse,"^",LRInData$RecepInMouse)

colnames(ResultMatrix) = CT
colnames(ResultMatrixPval) = CT

nCores = 8

cl <- makeCluster(nCores, outfile="")
registerDoParallel(cl)

ko_param <- foreach(I=seq_len( nrow(LRInData)),.errorhandling = "pass",
                    .packages = c("doParallel",
                                  "Rfast2",
                                  "tidyverse",
                                  "Matrix",
                                  "quadprog"
                    )
) %dopar% {
  
  set.seed(1241453)
  ##############
  project_stochastic_order <- function(P_hat, Sigma) {
    p <- length(P_hat)  # Number of dimensions
    
    # Compute the inverse of the covariance matrix
    Sigma_inv <- solve(Sigma)
    
    # Define the quadratic programming inputs
    Dmat <- 2 * Sigma_inv  # Quadratic term (must be positive definite)
    dvec <- 2 * Sigma_inv %*% P_hat  # Linear term
    
    # Constraint matrix for order restrictions P_1 <= P_2 <= ... <= P_p
    Amat <- matrix(0, nrow = p , ncol = p)
    for (i in 1:p) {
      # Amat[i, i] <- -1
      Amat[i, i] <- 1
    }
    
    # Constraint vector (enforcing increasing order)
    bvec <- rep(0.5, p)
    
    # Solve the quadratic programming problem
    result <- solve.QP(Dmat, dvec, t(Amat), bvec )
    
    return(result$solution)  # Projected vector P_bar
  }
  ###################
  
  Rx = LRInData[I,c("LiganInMouse","RecepInMouse")]
  
  # Communication pair
  
  L = Rx$LiganInMouse
  R = Rx$RecepInMouse
  
  
  
  for(kprim1 in 1:length(uniq_ctype)){
    
    for (kprim2 in 1:length(uniq_ctype)) {
      
      L_effect = treeEffect %>%filter(Ctype == uniq_ctype[kprim1]) %>%dplyr::select_at(L)
      R_effect = treeEffect %>%filter(Ctype == uniq_ctype[kprim2]) %>%dplyr::select_at(R)
      LR_effect = expand_grid(L = abs(L_effect[,1]),R = abs(R_effect[,1]))
      
      RL = LR_effect %>%as.matrix()
      
      #if(nrow(L_effect)<5) next
      
      
      l = 500
      
      Xij = LR0_effect[sample(1:nrow(LR0_effect),l,replace = T),]
      Xji = RL[sample(1:nrow(RL),l,replace = T),] 
      
      
      n1 = nrow(Xij)
      n2 = nrow(Xji)
      n= n1+n2
      
      # Compute P-hat
      
      p1_hat_rex_1 = (expand_grid(Xij[ ,1],Xji[ ,1]))
      p1_hat_rex_2 = (expand_grid(Xij[ ,2],Xji[ ,2]))
      
      p1_hat_mean = c(mean(p1_hat_rex_1[,1]<=p1_hat_rex_1[,2]),
                      mean(p1_hat_rex_2[,1]<=p1_hat_rex_2[,2]))
      
      
      # Bootstrap samples to calculate the covariance of P_j
      p1_Hat = NULL
      for (l in 1:20) {
        
        
        n11 = 200
        n22 = 200
        
        
        id_1 = sample(1:n1,n11,replace = T)
        id_2 = sample(1:n2,n22,replace = T)
        
        p1_hat_rex_1 = (expand_grid(Xij[id_1,1],Xji[id_1,1]))
        p1_hat_rex_2 = (expand_grid(Xij[id_2,2],Xji[id_2,2]))
        
        p1_hat = c(mean(p1_hat_rex_1[,1]<=p1_hat_rex_1[,2]),
                   mean(p1_hat_rex_2[,1]<=p1_hat_rex_2[,2]))
        
        p1_Hat = rbind(p1_Hat,p1_hat)
      }
      
      Sigma0 = cov(p1_Hat)
      
      # estimated probabilities & variance
      
      P_hat <- p1_hat_mean  
      Sigma <- Sigma0
      
      # Estimate P_bar
      
      aux_result <- try({
        
        P_bar <- project_stochastic_order(P_hat, Sigma)
        
      }, silent = TRUE) 
      
      if (class(aux_result)=="try-error") {
        message(paste("Error at iteration",uniq_ctype[kprim1],"^",uniq_ctype[kprim2]))
        Dn_Derived = NA
        p_value = NA
        
        ResultMatrix[I,paste0(uniq_ctype[kprim1],"^",uniq_ctype[kprim2])] = Dn_Derived
        ResultMatrixPval[I,paste0(uniq_ctype[kprim1],"^",uniq_ctype[kprim2])] = p_value
        
        next
      }
      P_bar <- aux_result
      #### Calculate Dn statistic
      P0 = matrix(c(0.5,0.5))
      Dn_Derived = drop(n*((t(P_hat-P0)%*%Sigma%*%(P_hat-P0)) - (t(P_hat-P_bar)%*%Sigma%*%(P_hat-P_bar))))
      
      
      p_value = mean(DN>Dn_Derived)
      
      
      ResultMatrix[I,paste0(uniq_ctype[kprim1],"^",uniq_ctype[kprim2])] = Dn_Derived
      ResultMatrixPval[I,paste0(uniq_ctype[kprim1],"^",uniq_ctype[kprim2])] = p_value
      
    }
  }
  
  
  cat("",(I/nrow(LRInData))*100,"\n")
  
  ####################
  
  
  #  
  
  list(
    ResultMatrix = ResultMatrix[I,],
    ResultMatrixPval = ResultMatrixPval[I,]
  )
}

parallel::stopCluster(cl)

for (I in seq_len(nrow(LRInData))) {
  if(is.null(ko_param[[I]]$ResultMatrix)) {Error = ko_param[[I]];next}
  Error = 0
  ResultMatrix[I,] =  ko_param[[I]]$ResultMatrix
  ResultMatrixPval[I,] =  ko_param[[I]]$ResultMatrixPval
}

return(list(ResultMatrix= ResultMatrix,
            ResultMatrixPval = ResultMatrixPval,
            ko_param=ko_param)
       
       )

}



GetNetWk = function(ResultMatrix,ResultMatrixPval,
                    subdata,
                    mst,
                    cutoff=20){


ResultMatrixPval_adj = matrix(p.adjust(as.vector(ResultMatrixPval), method =  "BY"),nrow = nrow(ResultMatrixPval),byrow = F)

ResultMatrixPval_adj2 = matrix(as.numeric(as.vector(ResultMatrixPval_adj)<0.05),nrow = nrow(ResultMatrixPval_adj),byrow = F)

 
rownames(ResultMatrixPval_adj) = rownames(ResultMatrixPval_adj2) = rownames(ResultMatrixPval)
colnames(ResultMatrixPval_adj) = colnames(ResultMatrixPval_adj2) = colnames(ResultMatrixPval)

M  = ResultMatrix
SigComm = melt(ResultMatrixPval_adj2)
SigComm = SigComm %>% filter(value==1)

SigComm$value2 = unlist(lapply(1:nrow(SigComm), function(k) ResultMatrix[SigComm$Var1[k],SigComm$Var2[k] ] ))
 
a = SigComm$value2
names(a) = SigComm$Var2



df = data.frame(from = str_extract_part(names(a),"^",before = T),to = str_extract_part(names(a),"^",before = F),score = a)
 
df = df %>%arrange(desc(score)) 
df$id = 1:nrow(df)
df = df %>% filter(id<cutoff)

df$transparency <- (max(df$score) - df$score) / (max(df$score) - min(df$score))

# Unique senders and receivers
senders <- unique(df$from)
receivers <- unique(df$to)

library(RColorBrewer)
# Dynamically assign colors
sender_colors <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(senders)), senders)
receiver_colors <- setNames(colorRampPalette(brewer.pal(8, "Set3"))(length(receivers)), receivers)
grid.col <- c(sender_colors, receiver_colors)

# Draw the Circos plot with arrows
circos.clear()
circos.par(gap.degree = 5)

chordDiagram(
  x = df,
  grid.col = grid.col,
  col = colorRamp2(c(min(df$score), max(df$score)), c("lightblue", "darkred"))(df$score),
  transparency = df$transparency,
  directional = 1,
  # direction.type = "arrows",   # <- This adds arrowheads
  # link.arr.type = "big.arrow", # <- Big arrow style
  annotationTrack = "grid",
  preAllocateTracks = 1
)

# Add sector labels
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.name <- get.cell.meta.data("sector.index")
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], sector.name,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.6)
}, bg.border = NA)


# Show ligand receptor pairs on a table
sa = union(df$from,df$to)
nam_a = names(sort(a,decreasing = T))
nam_a = nam_a[(str_extract_part(nam_a,"^",before = T)%in% sa & str_extract_part(nam_a,"^",before = F)%in%sa)]
nam_a = nam_a[1:pmin(cutoff,length(nam_a))]
b_aux = NULL

for(kk in 1:length(nam_a)){
  
  aux = sort(M[,nam_a[kk]],decreasing = T)%>%head(n=1)
  
  b_aux = rbind(b_aux, data.frame(Celltype= nam_a[kk],Ligand = str_extract_part(names(aux),"^",before = T),Receptor = str_extract_part(names(aux),"^",before = F)))
  
}
b_aux$Celltype = paste(str_extract_part(b_aux$Celltype,"^",before = T),
                       str_extract_part(b_aux$Celltype,"^",before = F),sep=" --> ")

colnames(b_aux) = c("Sender->Receiver Cell", "Ligand", "Receptor")
b_aux = unique(b_aux)
tg <- tableGrob(b_aux, rows = NULL)

for (i in seq_along(tg$grobs)) {
  if (inherits(tg$grobs[[i]], "grob")) {
    tg$grobs[[i]]$gp$fill <- NA
  }
}

return(list(tg = tg, df=df))


}

 



Graphh = function(Result,
                  subdata,
                  META_subb,
                  mst,
                  tissue,
                  sa, # Cell types to view
                  antibody = "Itgb2"
                  ){
  Healthy = tissue

  
  
###################
  
  AA = GetId(
    mst,
    ExprsData = subdata%>%as.data.frame(),
    ClusterCol ="NeighXCell_id" ,
    TemporalCol ="day.harvested",
    useWeight = FALSE
  )
  
  #####
  netwknode = AA$netwknode
  netwknode = bind_cols(netwknode,Result$treeEffect)
  TreeTemp_data = AA$TreeTemp_data
  
  TreeTemp_data = left_join(TreeTemp_data,netwknode,by="TreeTemp")
  
  TreeTemp_data = bind_cols(META_subb,
                            subdata[,c("NeighXCell_id","day.harvested")],
                            TreeTemp_data)
  
###################


id = TreeTemp_data$Slice_ID%in%Healthy[1]
id2 = TreeTemp_data$Tier3 %in% sa
id = id&id2
 


if(is.null(antibody)){
  plotScatter(TreeTemp_data$x[id],TreeTemp_data$y[id],TreeTemp_data$Tier3[id],main="",size=.8,
              legend.size = 1.5,ManualColor = T,cols = c25)
}else{
  
  plotScatter(TreeTemp_data$x[id],TreeTemp_data$y[id],TreeTemp_data[id,antibody],main=antibody,size=.3,
                        legend.size = 1.5,ManualColor = F,cols = c25)
  
}


}


 

HeatMap = function(ResultMatrix,
                   ResultMatrixPval,
                   sa,
                   show.onlySig=TRUE,
                   fontsize_row = 4,
                   fontsize_col = 7,
                   zoom=0.5,
                   toCell = NULL){
  
  ResultMatrixPval_adj = matrix(p.adjust(as.vector(ResultMatrixPval), method =  "BY"),nrow = nrow(ResultMatrixPval),byrow = F)
  
  ResultMatrixPval_adj2 = matrix(as.numeric(as.vector(ResultMatrixPval_adj)<0.05),nrow = nrow(ResultMatrixPval_adj),byrow = F)
  
  
  rownames(ResultMatrixPval_adj) = rownames(ResultMatrixPval_adj2) = rownames(ResultMatrixPval)
  colnames(ResultMatrixPval_adj) = colnames(ResultMatrixPval_adj2) = colnames(ResultMatrixPval)
  M  = ResultMatrix
  
  aux_nam = str_extract_part(colnames(M),"^",before = T)%in%sa &  str_extract_part(colnames(M),"^",before = F)%in%sa
M2  = M[,aux_nam,drop=F]


if (show.onlySig){
M2[ResultMatrixPval_adj2[,aux_nam]==0] = NA
}

if(is.null(toCell)){
  
  M_F = M2
  
}else{
  
  M_F = M2[,str_extract_part(colnames(M2),"^",before = F)==toCell]
  
}



global_min <- 0
global_max <- max(M2,na.rm = T)

breaks <- seq(global_min, global_max, length.out = 50)
M_F2 = M_F
M_F2[is.na(M_F2)] = 0
 
M_F = M_F[rowSums(M_F2)!=0,]

M_F = M_F[rowSums(M_F,na.rm = T)>zoom,]

if (show.onlySig){
  pheatmap::pheatmap(M_F , 
                     treeheight_row = 0,
                     treeheight_col = 0,
                     fontsize = 7,
                     fontsize_row = fontsize_row,
                     fontsize_col = fontsize_col,
                     angle_col = 90,
                     cluster_rows = F,
                     cluster_cols = F,
                     na_col = "grey",
                     breaks = breaks,
                     color = inferno(49))
}else{
  pheatmap::pheatmap(M_F , 
                     treeheight_row = 0,
                     treeheight_col = 0,
                     fontsize = 7,
                     fontsize_row = fontsize_row,
                     fontsize_col = fontsize_col,
                     angle_col = 90,
                     cluster_rows = T,
                     cluster_cols = T,
                     na_col = "grey",
                     breaks = breaks,
                     color = inferno(49))
  
}

}



StochasticCC = function(metaData,
                        datExpr,
                        day_var ="Sample_type",
                        CTn = CTn
){
  
  HH = NetworkGet(metaData,datExpr)
  
  # Extract Processed data
  
  datExpr = HH$datExpr_withNeigborhood
  Centers2 = HH$Centers
  ax_reduce_META = HH$ax_reduce_META
  datExpr$day.harvested = as.numeric(ax_reduce_META[,day_var])
  
  # Estimate network
  set.seed(12345)
  mst = ClusterToTree(Centers = Centers2)
  
  
  # Estimate phi for day, say 3
  ResultMatrix =TreeTemp_data_effect2= ResultMatrixPval=BB=LR0_effect=DN=TreeTemp_data=TreeTemp_data_effect=META_subb=LRInData=ctype_day1=subdata=Result=DEGs=list()
  
  for (lk in unique(datExpr$day.harvested )) {
    
    DEGs[[lk]] = FindDEGs(datExpr,
                          ax_reduce_META,
                          mst,
                          day=lk 
    )
    
    
    # Collect outputs
    Result[[lk]]  = DEGs[[lk]]$Result_Cancer_day
    subdata[[lk]] = DEGs[[lk]]$subdata 
    ctype_day1[[lk]]  =  DEGs[[lk]]$ctype_day1
    LRInData[[lk]] = DEGs[[lk]]$LRInData
    META_subb[[lk]] = DEGs[[lk]]$META_subb
    
    # Extract effect
    
    TreeTemp_data_effect = Teffect(Result[[lk]],
                                   mst,
                                   subdata[[lk]],
                                   ctype_day1[[lk]])
    
    TreeTemp_data_effect2[[lk]] = TreeTemp_data_effect$TreeTemp_data_effect
    TreeTemp_data[[lk]] = TreeTemp_data_effect2[[lk]]$TreeTemp_data
    
    
    # Get Null distribution of communication scores
    
    AA = getNull(TreeTemp_data_effect2[[lk]],LRInData[[lk]])
    
    DN[[lk]] = AA$DN
    LR0_effect[[lk]] = AA$LR0_effect
    
    # Estimate cell-cell communication
    BB[[lk]] = ComputeCCS(TreeTemp_data_effect2[[lk]],
                          LRInData[[lk]],
                          CTn,
                          DN[[lk]],
                          LR0_effect[[lk]])
    
    
    ResultMatrix[[lk]]     =     BB[[lk]]$ResultMatrix
    ResultMatrixPval[[lk]] =     BB[[lk]]$ResultMatrixPval
    
  }
  
  
  return(list(datExpr=datExpr,
              ax_reduce_META= ax_reduce_META,
              Result = Result,
              ResultMatrix=ResultMatrix,
              ResultMatrixPval = ResultMatrixPval,
              mst=mst,
              subdata=subdata,
              META_subb= META_subb,
              LRInData = LRInData
              
  )
  )
}
###################################################
############### DYNAMIC FUNCTIONS ##################
####################################################


StochasticCCDynm = function(Re= R,
                            LRInData = NULL,
                            ctype_from = "Colonocytes",
                            ctype_to = "Stem cells",
                            uniq_day = c(1,2,3),
                            nCores = 9,
                            usePvalue=TRUE
){
  
  datExpr_withNeigborhood_subdata = Re$datExpr
  ax_reduce_META_subb = Re$ax_reduce_META
  
  ConfoundFrame = data.frame(ax_reduce_META_subb$Slice_ID %>%
                               as.numeric())
  
  if(is.null(LRInData)){
    LRInData = Re$LRInData[[1]]
  }
  
  
  TreeEffect = NULL
  
  for(day in uniq_day ){
    treeEffect = Re$Result[[day]]$treeEffect %>% as.data.frame()
    treeEffect$NeighXCellXclust_id =  unique(as.numeric(as.factor(ax_reduce_META_subb$NeighXCellXclust) )) %>% sort
    ax_reduce_META_subb$NeighXCellXclust_id = as.numeric(as.factor(ax_reduce_META_subb$NeighXCellXclust))
    
    treeEffect = left_join(ax_reduce_META_subb, treeEffect,by="NeighXCellXclust_id")
    treeEffect$day.harvested = day
    
    TreeEffect = rbind(TreeEffect,treeEffect)
  }
  
  datExpr_withNeigborhood_subdata = TreeEffect
  
  bivariate_monotonic_test <- function(groups) {
    k <- length(groups)
    ns <- sapply(groups, nrow)
    pooled <- do.call(rbind, groups)
    
    # Grid for ECDF
    x_grid <-  seq(min(sort(unique(pooled[,1]))),max(sort(unique(pooled[,1]))),length=100) # sort(unique(pooled[,1]))
    y_grid <-  seq(min(sort(unique(pooled[,2]))),max(sort(unique(pooled[,2]))),length=100) #sort(unique(pooled[,2]))
    
    m1 <- length(x_grid)
    m2 <- length(y_grid)
    
    # Empirical joint CDFs
    Fhat_array <- Fhat_array_scale <- array(0, dim = c(k, m1, m2))
    for(a in 1:m1){
      for(b in 1:m2){
        for(i in 1:k){
          Xi <- groups[[i]]
          Fhat_array[i,a,b] <- mean(Xi[,1] <= x_grid[a] & Xi[,2] <= y_grid[b])
        }
        Fhat_array_scale[,a,b] = (1-Fhat_array[,a,b])/sum(1-Fhat_array[,a,b])
      }
    }
    Simplex =  cbind(expand_grid(b=x_grid, a=y_grid), T1 = as.numeric(Fhat_array_scale[1,,]),
                     T2 = as.numeric(Fhat_array_scale[2,,]),
                     T3 = as.numeric(Fhat_array_scale[3,,]))         
    
    # Isotonic fit via PAVA
    Fiso_array <- array(0, dim = c(k, m1, m2))
    for(a in 1:m1){
      for(b in 1:m2){
        Fiso_array[,a,b] <- pava(1-Fhat_array[,a,b], w = ns) # Fhat_array[,a,b] is decreasing across time, however, 1-Fhat_array[,a,b] is increasing
      }
    }
    
    # Weights for integration
    w_matrix <- matrix(1/(m1*m2), nrow = m1, ncol = m2)
    
    
    
    # Test statistic
    T_obs <- sum( w_matrix * apply((Fhat_array - Fiso_array)^2, c(2,3), function(z) sum(ns * z)) )
    
    # Permutation test
    T_obs
    B <- 100
    Tperm <- numeric(B)
    group_id <- rep(1:k, times = ns)
    
    for(b in 1:B){
      perm_labels <- sample(group_id)
      perm_groups <- lapply(1:k, function(i) pooled[perm_labels == i, , drop = FALSE])
      
      Fh_perm <- array(0, dim = c(k, m1, m2))
      for(i in 1:k){
        Xi <- perm_groups[[i]]
        for(a in 1:m1){
          for(bb in 1:m2){
            Fh_perm[i,a,bb] <- mean(Xi[,1] <= x_grid[a] & Xi[,2] <= y_grid[bb])
          }
        }
      }
      
      Fiso_perm <- array(0, dim = c(k, m1, m2))
      for(a in 1:m1){
        for(bb in 1:m2){
          Fiso_perm[,a,bb] <- pava(Fh_perm[,a,bb], w = ns)
        }
      }
      
      Tperm[b] <- sum(w_matrix * apply((Fh_perm - Fiso_perm)^2, c(2,3), function(z) sum(ns * z)))
    }
    
    pval <- (1 + sum(Tperm <= T_obs)) / (1 + B)
    list(pval=pval,T_obs =T_obs, Fhat_array=Fhat_array, Simplex=Simplex,Fiso_perm=Fiso_perm)
  }
  
  detected =detectedObs =NULL
  
  
  
  library(doParallel)
  nCores = nCores
  cl <- makeCluster(nCores, outfile="")
  registerDoParallel(cl)
  
  replic = length(LRInData$RecepInMouse)
  detected = detectedObs = L = R = vector("numeric", replic)
  P = L=R= matrix(NA,nrow = length(uniq_day), ncol = replic)
  
  
  ko_param <- foreach(I=seq_len(replic),.errorhandling = "pass",
                      .packages = c("doParallel",
                                    "Rfast2",
                                    "tidyverse",
                                    "Iso"
                                    
                      )
  ) %dopar% {
    
    
    set.seed(12353)
    
    if(I%%100==0)cat(I,"\n")
    mono_data = list()
    
    for (day in uniq_day) {
      
      id1l = which(datExpr_withNeigborhood_subdata$day.harvested==day & datExpr_withNeigborhood_subdata$Tier3==ctype_from) %>% sample(1500, replace = T)
      id1r = which(datExpr_withNeigborhood_subdata$day.harvested==day &datExpr_withNeigborhood_subdata$Tier3==ctype_to)%>% sample(1500, replace = T)
      
      mono_data[[day]] = cbind(Ligand = abs(datExpr_withNeigborhood_subdata[id1l,LRInData$LiganInMouse[I]]),
                               Receptor = abs(datExpr_withNeigborhood_subdata[id1r,LRInData$RecepInMouse[I]] ))
      
    }
    
    pval_mono <- bivariate_monotonic_test(mono_data)
    
    detected[I] <- pval_mono$pval
    detectedObs[I] <- pval_mono$T_obs
    Fhat_array <- pval_mono$Fhat_array
    Fiso_perm <- pval_mono$Fiso_perm
    Simplex = pval_mono$Simplex
    
    aL= list()
    aR= list()
    
    for (day in uniq_day) {
      
      aL[[day]] = mono_data[[day]][,1,drop=F] %>% colMeans()
      aR[[day]] = mono_data[[day]][,2,drop=F] %>% colMeans()
    }
    
    L[,I] =  unlist(aL)        #c(a1[1],a2[1],a3[1])
    R[,I]  = unlist(unlist(aR))# c(a1[2],a2[2],a3[2])
    
    aux1 = Fhat_array[,sample(1:dim(Fhat_array)[2],10),sample(1:dim(Fhat_array)[3],1)]
    aux2 = Fhat_array[,sample(1:dim(Fhat_array)[2],1),sample(1:dim(Fhat_array)[3],10)]
    
    ranked_aux1 <- apply(1-aux1, 2, rank)
    ranked_aux2 <- apply(1-aux2, 2, rank)
    
    rownames(ranked_aux1) = c("T1","T2","T3")
    rownames(ranked_aux2) = c("T1","T2","T3")
    
    
    P[,I] = colMeans(rbind(rowMeans(ranked_aux1),rowMeans(ranked_aux2)))
    
    list(
      detected  =  detected[I],
      detectedObs   =  detectedObs[I],
      L =  L[,I],
      R =  R[,I],
      P = P[,I],
      Simplex = Simplex
    )
  }
  
  Simplex =list()
  for (i in 1:replic) {
    if(is.null(ko_param[[i]]$detected)) {Error = ko_param[[i]];next}
    # Get smoothed/predicted gene expression
    Error = 0
    detected[i] = ko_param[[i]]$detected
    detectedObs[i] = ko_param[[i]]$detectedObs
    L[,i] = ko_param[[i]]$L
    R[,i] = ko_param[[i]]$R
    P[,i] = ko_param[[i]]$P
    Simplex[[i]] = ko_param[[i]]$Simplex
  }
  
  parallel::stopCluster(cl) 
  
  
  colnames(L) = paste0(LRInData$LiganInMouse,"^",LRInData$RecepInMouse)
  colnames(R) = paste0(LRInData$LiganInMouse,"^",LRInData$RecepInMouse)
  colnames(P) = paste0(LRInData$LiganInMouse,"^",LRInData$RecepInMouse)
  names(Simplex) = paste0(LRInData$LiganInMouse,"^",LRInData$RecepInMouse)
  
  rownames(L) = c("T1","T2","T3")
  rownames(R) = c("T1","T2","T3")
  rownames(P) = c("T1","T2","T3")
  
  
  names(detected) = paste0(LRInData$LiganInMouse,"^",LRInData$RecepInMouse)
  names(detectedObs) = paste0(LRInData$LiganInMouse,"^",LRInData$RecepInMouse)
  
  return(
    list(detected = detected,
         detectedObs = detectedObs,
         Simplex= Simplex,
         L = L,
         R= R,
         P = P,
         Dirtn = uniq_day)
  )
}



GetDynamicNtwk = function(Re= Re,
                          LRInData = NULL,
                          ctype_from = "Colonocytes",
                          ctype_to = "Stem cells",
                          uniq_day = Test,
                          nCores = 9,
                          usePvalue =TRUE,
                          NoG = 5
){
  
  ORR =  NULL
  Simplex = DetectedObs = Detected =list()
  
  for(kl in 1:length(Test) ){
    
    R2 = StochasticCCDynm(Re= Re,
                          LRInData = LRInData,
                          ctype_from = ctype_from,
                          ctype_to = ctype_to,
                          uniq_day = uniq_day[[kl]],
                          nCores = nCores,
                          usePvalue=usePvalue
    )
    
    
    detectedObs = R2$detectedObs
    detected =  R2$detected
    Simplex[[paste0("T",uniq_day[[kl]],collapse = "")]] = R2$Simplex
    
    NoG = pmin(NoG, length(detected))
    
    nam = sort(detectedObs) %>% names() #names(a)
    
    if(usePvalue){
      
      auX = detected[nam][detected[nam] <0.05]
      
      if(length(auX)==0){
        nam = nam[1:NoG]
      }else{
        nam = names(auX)
      }
      
    }else{
      nam = nam[1:NoG]
    }
    
    
    
    OR=NULL
    for(k in 1:length(nam)){
      OR = rbind(OR ,order = paste0("T",uniq_day[[kl]]) )
    }
    
    
    rownames(OR) = nam
    OR  = cbind(OR,detectedObs[nam])
    
    ORR = rbind(ORR,OR)
    
  }
  
  colnames(ORR) = c("V1","V2","V3","V4")
  
  return( list(OR =ORR,
               Simplex = Simplex)
  )
  
}
