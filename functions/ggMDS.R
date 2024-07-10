############## MDS plots ----- 

ggMDS <- function(
    dfLonger, ## dflonger = stacked dataset
    sampleIdVar='SampleID',
    bioGranVar ='ASV_name',
    countingVar='Raw.Sequence.Counts',
    clusterVar=''
  ){
  
  #### 
  dfLonger$IdVar = dfLonger %>% select(one_of(sampleIdVar)) %>% pull()
  
  #### Firt step is to create the composition
  dfComposition = grump2comp(dfLonger,
                             sampleIdVar = sampleIdVar,
                             bioGranVar =  bioGranVar,
                             countingVar = countingVar,
                             fillPercentage =T)
  
  # creating the atichison distance
  aitDist <- vegan::vegdist(dfComposition[,-1],method='aitchison')
  
  # finding the MDS 
  mds_obj<-cmdscale(
    d = aitDist,
    k = 3,eig = T,add=T
  )
  pct_explained = round(100 * mds_obj$eig/sum(mds_obj$eig),1)

  dfComposition = dfComposition %>% select(one_of(sampleIdVar)) %>% mutate(
    MDS1=mds_obj$points[,1],
    MDS2=mds_obj$points[,2],
    MDS3=mds_obj$points[,3]
  )
  
  if(clusterVar!=''){
    dfLonger$ClusterAux = dfLonger %>% dplyr::select(one_of(clusterVar)) %>% pull()
    dfCLusterAux = dfLonger %>% select(IdVar,ClusterAux) %>% distinct()
    colnames(dfCLusterAux) = c(sampleIdVar,'Cluster')
    dfComposition = dfComposition %>% left_join(dfCLusterAux)
    
    #### Plot 2D
    outputAit2d <- ggplot(
      data = dfComposition,
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
    
      outputAit3d <- plotly::plot_ly(
        x=dfComposition$MDS1,
        y=dfComposition$MDS2,
        z=dfComposition$MDS3,
        type="scatter3d",
        color = as.numeric(dfComposition$Cluster),
        mode="markers") %>%
        #plotly::layout(
        #  scene = list(aspectratio = list(x = 4, y = 4, z = 3))) %>% 
        #list(xaxis=list(title='MDS1'),yaxis=list(title='MDS2'),zaxis=list(title='MDS3'))) %>% 
        plotly::layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz)) %>% 
        plotly::hide_colorbar() 
    
  }else{
    outputAit2d=ggplot(
      data = dfComposition,
      mapping = aes(x=MDS1,y=MDS2))+
      geom_label(alpha=0.3) +
      theme_minimal() +
      xlab(paste('MDS1 -',pct_explained[1],'%'))+
      ylab(paste('MDS2 -',pct_explained[2],'%'))+
      theme(legend.position = '') 
    
    axx <- list(
      title = "MDS1"
    )
    
    axy <- list(
      title = "MDS2"
    )
    
    axz <- list(
      title = "MDS3"
    )
    
    outputAit3d <- plotly::plot_ly(
      x=dfComposition$MDS1,
      y=dfComposition$MDS2,
      z=dfComposition$MDS3,
      type="scatter3d",
      mode="markers") %>%
      #plotly::layout(
      #  scene = list(aspectratio = list(x = 4, y = 4, z = 3))) %>% 
      #list(xaxis=list(title='MDS1'),yaxis=list(title='MDS2'),zaxis=list(title='MDS3'))) %>% 
      plotly::layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz)) %>% 
      plotly::hide_colorbar() 
  }
  
  
  ##### Now using Bray Curtiss
  #### 
  dfLonger$IdVar = dfLonger %>% select(one_of(sampleIdVar)) %>% pull()
  
  #### Firt step is to create the composition
  dfComposition = grump2comp(dfLonger,
                             sampleIdVar = sampleIdVar,
                             bioGranVar =  bioGranVar,
                             countingVar = countingVar,
                             fillPercentage =F)
  
  # creating the atichison distance
  bcDist <- vegan::vegdist(dfComposition[,-1],method='bray')
  
  # finding the MDS 
  mds_obj<-cmdscale(
    d = bcDist,
    k = 3,eig = T,add=T
  )
  pct_explained = round(100 * mds_obj$eig/sum(mds_obj$eig),1)
  
  dfComposition = dfComposition %>% select(one_of(sampleIdVar)) %>% mutate(
    MDS1=mds_obj$points[,1],
    MDS2=mds_obj$points[,2],
    MDS3=mds_obj$points[,3]
  )
  
  if(clusterVar!=''){
    dfLonger$ClusterAux = dfLonger %>% dplyr::select(one_of(clusterVar)) %>% pull()
    dfCLusterAux = dfLonger %>% select(IdVar,ClusterAux) %>% distinct()
    colnames(dfCLusterAux) = c(sampleIdVar,'Cluster')
    dfComposition = dfComposition %>% left_join(dfCLusterAux)
    
    #### Plot 2D
    outputBC2d <- ggplot(
      data = dfComposition,
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
    
    outputBC3d <- plotly::plot_ly(
      x=dfComposition$MDS1,
      y=dfComposition$MDS2,
      z=dfComposition$MDS3,
      type="scatter3d",
      color = as.numeric(dfComposition$Cluster),
      mode="markers") %>%
      #plotly::layout(
      #  scene = list(aspectratio = list(x = 4, y = 4, z = 3))) %>% 
      #list(xaxis=list(title='MDS1'),yaxis=list(title='MDS2'),zaxis=list(title='MDS3'))) %>% 
      plotly::layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz)) %>% 
      plotly::hide_colorbar() 
    
  }else{
    outputBC2d=ggplot(
      data = dfComposition,
      mapping = aes(x=MDS1,y=MDS2))+
      geom_label(alpha=0.3) +
      theme_minimal() +
      xlab(paste('MDS1 -',pct_explained[1],'%'))+
      ylab(paste('MDS2 -',pct_explained[2],'%'))+
      theme(legend.position = '') 
    
    axx <- list(
      title = "MDS1"
    )
    
    axy <- list(
      title = "MDS2"
    )
    
    axz <- list(
      title = "MDS3"
    )
    
    outputBC3d <- plotly::plot_ly(
      x=dfComposition$MDS1,
      y=dfComposition$MDS2,
      z=dfComposition$MDS3,
      type="scatter3d",
      mode="markers") %>%
      #plotly::layout(
      #  scene = list(aspectratio = list(x = 4, y = 4, z = 3))) %>% 
      #list(xaxis=list(title='MDS1'),yaxis=list(title='MDS2'),zaxis=list(title='MDS3'))) %>% 
      plotly::layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz)) %>% 
      plotly::hide_colorbar() 
  }
  
  return(list(
    MDS_Bray2d = outputBC2d,
    MDS_Bray3d = outputBC3d,
    MDS_Ait2d = outputAit2d,
    MDS_Ait3d = outputAit3d
  ))
}

########### Example 
#example::

# dfZoop<-data.table::fread('/Users/rafaelcatoia/Desktop/repos/Capstone/data/grump_Zoop_16NS_longer.csv')
# dfTest <- grump2comp(dfZoop,sampleIdVar = 'SampleID',bioGranVar = 'ASV_name',
#                      countingVar = 'Raw.Sequence.Counts',fillPercentage =T)
# adist <- robCompositions::aDist(x = dfTest[,-1])
# plot(hclust(adist,method = 'ward.D2'))
# dfClusters = dfTest %>% select(SampleID) %>%
#   bind_cols( Clust = factor(cutree(hclust(adist,method = 'ward.D2'),k = 6),levels=1:6))
# 
# ggMDS_obj <- ggMDS(dfZoop %>% left_join(dfClusters),
#                 sampleIdVar = 'SampleID',clusterVar = 'Clust')

