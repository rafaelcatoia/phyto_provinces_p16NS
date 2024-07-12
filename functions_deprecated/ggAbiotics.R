ggAbiotics <- function(
    dfLonger,
    sampleIdVar='SampleID',
    clusterVar='',
    abioticVars = c(
      'Temperature','Salinity','Oxygen',
      'Silicate','NO2','NO3','NOx','PO4')
    ){
  
  library(GGally)
  out <- list()
  
  dfLonger$Cluster = dfLonger %>% dplyr::select(one_of(clusterVar)) %>% pull()
  dfLonger$ID = dfLonger %>% dplyr::select(one_of(sampleIdVar)) %>% pull()
  
  
  df_Abio <- dfLonger %>% select(ID,Cluster,Latitude,Depth,any_of(abioticVars)) %>% 
    distinct()

  ############# ggpairs

  abioticPairs <- GGally::ggpairs(df_Abio %>% select(-ID,-Latitude,-Depth),
                                  columns =2:9,upper = NULL)+
    theme_minimal()
  
  ############# ggpairs colored
  
  ## Function to plot only the Scatterplots
  gpairs_lower <- function(g){
    g$plots <- g$plots[-(1:g$nrow)]
    g$yAxisLabels <- g$yAxisLabels[-1]
    g$nrow <- g$nrow -1
    
    g$plots <- g$plots[-(seq(g$ncol, length(g$plots), by = g$ncol))]
    g$xAxisLabels <- g$xAxisLabels[-g$ncol]
    g$ncol <- g$ncol - 1
    
    g
  }
  
  abioticPairsColored = GGally::ggpairs(
    aes(label = Cluster,color= Cluster),
    data = df_Abio %>% select(-ID,-Latitude,-Depth),
    columns =2:9,upper = NULL,diag = list(continuous = "blankDiag"),
    lower=list(continuous = wrap("points", size = 0,alpha=0.00001)))+
    #geom_point(color='white')+
    geom_text(aes(label =Cluster), alpha = 0.5)+
    theme_minimal()
  
  
  df_Abio_longer = df_Abio %>%
    tidyr::pivot_longer(cols = any_of(abioticVars),names_to = 'AbioticVar' )
  
  out$abioticMatrix = df_Abio_longer %>% 
    group_split(AbioticVar) %>%
    purrr::map(
      ~ggplot(., aes(Latitude, Depth, color = value,label=Cluster)) + 
        geom_label(alpha=0.3) +
        scale_colour_gradient2(
          low = "#C1FFFF", 
          mid = "#007CE3", 
          high = "#06004D", 
          midpoint = median(.$value)
        ) +    
        scale_y_reverse() + theme_minimal()+
        facet_wrap(~ AbioticVar, labeller = function(x) label_value(x, multi_line = FALSE))
    ) %>% 
    cowplot::plot_grid(plotlist = ., align = 'hv', ncol = 2)
  
  ############### MDS plot with only abiotic factors ----------------
  
  mds_obj <- cmdscale(d = dist(df_Abio[,-c(1,2)] %>% scale(),
                               #method = 'manhattan'),k = 3,eig = T,add=T)
                               method = 'euclid'),k = 3,eig = T,add=T)
  pct_explained = round(100 * mds_obj$eig/sum(mds_obj$eig),1)
  
  df_Abio <- df_Abio %>% mutate(
    MDS1 = mds_obj$points[,1],
    MDS2 = mds_obj$points[,2],
    MDS3 = mds_obj$points[,3]
  )
  
  outputMDS2dAbiotic = ggplot(
    data = df_Abio,
    mapping = aes(x=MDS1,y=MDS2,
                  color=Cluster,
                  label=Cluster))+
    geom_label(alpha=0.3) +
    theme_minimal() +
    xlab(paste('MDS1 -',pct_explained[1],'%'))+
    ylab(paste('MDS2 -',pct_explained[2],'%'))+
    theme(legend.position = '')
  
  #### Plot 3D
  axx <- list(
    title = "MDS1"
  )
  axy <- list(
    title = "MDS2"
  )
  axz <- list(
    title = "MDS3"
  )
  
  outputMDS3dAbiotic <- plotly::plot_ly(
    x=df_Abio$MDS1,
    y=df_Abio$MDS2,
    z=df_Abio$MDS3,
    type="scatter3d",
    color = as.numeric(df_Abio$Cluster),
    mode="markers") %>%
       plotly::layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz)) %>% 
    plotly::hide_colorbar() 
  
  ############ ----- princomp
  
  pca_obj = princomp(scale(df_Abio %>% select(one_of(abioticVars))))
  
  summary(pca_obj)
  
  df_Abio = df_Abio %>% mutate(
    PC1=pca_obj$scores[,1],
    PC2=pca_obj$scores[,2],
    PC3=pca_obj$scores[,3]
  )
  
  pct_explained_PCA  = round(100*pca_obj$sdev^2/sum(pca_obj$sdev^2),2)
  
  outputPCA2dAbiotic = ggplot(
    data = df_Abio,
    mapping = aes(x=PC1,y=PC2,
                  color=Cluster,
                  label=Cluster))+
    geom_label(alpha=0.3) +
    theme_minimal() +
    xlab(paste('PC1 -',pct_explained_PCA[1],'%'))+
    ylab(paste('PC2 -',pct_explained_PCA[2],'%'))+
    theme(legend.position = '')
  
  #### Plot 3D
  axx <- list(
    title = "PC1"
  )
  axy <- list(
    title = "PC2"
  )
  axz <- list(
    title = "PC3"
  )
  
  outputPCA3dAbiotic <- plotly::plot_ly(
    x=df_Abio$PC1,
    y=df_Abio$PC2,
    z=df_Abio$PC3,
    type="scatter3d",
    color = as.numeric(df_Abio$Cluster),
    mode="markers") %>%
    plotly::layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz)) %>% 
    plotly::hide_colorbar() 
  
  
  
return(list(abioticPairs=abioticPairs,
            abioticPairsColored=gpairs_lower(abioticPairsColored),
            MDS2d = outputMDS2dAbiotic,
            MDS3d = outputMDS3dAbiotic,
            pca_obj = pca_obj,
            PCA2d = outputPCA2dAbiotic,
            PCA3d = outputPCA3dAbiotic
            ))  
}

#example::
# dfZoop<-data.table::fread('/Users/rafaelcatoia/Desktop/repos/Capstone/data/grump_Zoop_16NS_longer.csv')
# dfTest <- grump2comp(dfZoop,sampleIdVar = 'SampleID',bioGranVar = 'ASV_name',
#                      countingVar = 'Raw.Sequence.Counts',fillPercentage =T)
# adist <- robCompositions::aDist(x = dfTest[,-1])
# plot(hclust(adist,method = 'ward.D2'))
# dfClusters = dfTest %>% select(SampleID) %>%
#   bind_cols( Clust = factor(cutree(hclust(adist,method = 'ward.D2'),k = 6),levels=1:6))
# 
# ggAbiotics_obj <- ggAbiotics(dfLonger = dfZoop %>% left_join(dfClusters),
#                 sampleIdVar = 'SampleID',clusterVar = 'Clust')

