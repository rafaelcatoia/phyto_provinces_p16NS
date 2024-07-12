##  ggLatitudeDepth
##  Function to create a scatterplot containing the Latitude and Depth
##  as input, recieves a stacked dataset containing the sample ID,
####sample cluster membership, latitude, depht 

ggLatitudeDepth <- function(
  dfLonger, ## dflonger = stacked dataset
  sampleIdVar, ## SampleIdVar = name of the variable containing the sample ID
  clusterVar, ## Cluster Identifier for each sample
  title ='',
  CHull=F,
  OceanLayers=F,
  FacetClusters=F,
  LatitudeRegions=T,
  labelSize=12,
  baseSize=16
){
  
  dfLonger$Cluster = dfLonger %>% dplyr::select(one_of(clusterVar)) %>% pull()
  
  df_geo <- dfLonger %>% 
    select(SampleID=one_of(sampleIdVar),Latitude,Depth,Cluster) %>%
    distinct()
  
  output <- 
    ggplot(df_geo,
           aes(x=Latitude,y=Depth,color=Cluster,label=Cluster))+
    geom_vline(xintercept = c(-60,-45,-5,5,35),col='gray75')+
    #geom_hline(yintercept =0,col='gray75')+
    geom_label(alpha=0.3,size = labelSize/.pt)+
    scale_y_reverse() +
    geom_hline(yintercept = c(0),col='gray55')+
    theme_minimal(base_size = baseSize)+ 
    ggtitle(title)+
    theme(legend.position = 'none',
          panel.grid.major = element_blank())
  

  if(LatitudeRegions){
    output = output+
      annotate("text",size = labelSize*0.3, x = -67, y = -50, label = "Southern Ocean",hjust=0.5)+
      annotate("text",size = labelSize*0.3, x = -52, y = -50, label = "Subantarctic",hjust=0.5)+
      annotate("text",size = labelSize*0.3, x = -25, y = -50, label = "South Pacific Gyre",hjust=0.5)+
      annotate("text",size = labelSize*0.3, x = 0,   y = -50, label = "Equatorial",hjust=0.5)+
      annotate("text",size = labelSize*0.3, x = 20,  y = -50, label = "North Pacific Gyre",hjust=0.5)+
      annotate("text",size = labelSize*0.3, x = 45,  y = -50, label = "Subarctic",hjust=0.5)
  }
  
  if(CHull){
    chull_ <- df_geo %>% 
      group_by(Cluster) %>%
      slice(chull(Latitude, Depth))
    output =  output +
      geom_polygon(data = chull_, alpha = 0.2,aes(fill=Cluster))
  }
  
  return(output)
}

#example::
# dfZoop<-data.table::fread('/Users/rafaelcatoia/Desktop/repos/Capstone/data/grump_Zoop_16NS_longer.csv')
# dfTest <- grump2comp(dfZoop,sampleIdVar = 'SampleID',bioGranVar = 'ASV_name',
#                      countingVar = 'Raw.Sequence.Counts',fillPercentage =T)
# adist <- robCompositions::aDist(x = dfTest[,-1])
# plot(hclust(adist,method = 'ward.D2'))
# dfClusters = dfTest %>% select(SampleID) %>%
#   bind_cols( Clust = factor(cutree(hclust(adist,method = 'ward.D2'),k = 15),levels=1:15))
# 
# ggLatitudeDepth(dfZoop %>% left_join(dfClusters),
#                 sampleIdVar = 'SampleID',clusterVar = 'Clust')

                