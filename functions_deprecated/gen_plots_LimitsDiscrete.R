######## Genete plots ---------------------------------
# list_cluster_membership_and_bounderies = readRDS(file = paste0(savingdir,'/','list_cluster_membership_and_bounderies'))
# list_cluster_membership_and_bounderies_mirroredLat = readRDS(file = paste0(savingdir,'/','list_cluster_membership_and_bounderies_mirroredLat'))
# grid_base = readRDS(file = paste0(savingdir,'/','grid_base'))
# clust_member_limits = unlist_coloring_obj(list_cluster_membership_and_bounderies_mirroredLat)
# df_GeoAbio = readRDS(file = paste0(savingdir,'/','df_geo_abiotics'))

#gridBase = grid_base %>% rename('lat_grid'=1,'depth_grid'=2)
#clustMemLim = clust_member_limits
#df_GeoAbio =df_geo_abiotics

gen_plots <- function(gridBase,clustMemLim,df_GeoAbio){
  
  ## first lets bring LAT and Depht to it's original scale
  mean_Lat = df_GeoAbio$Latitude %>% mean()
  sd_Lat = df_GeoAbio$Latitude %>% sd()
  
  mean_Depth = df_GeoAbio$Depth %>% mean()
  sd_Depth = df_GeoAbio$Depth %>% sd()
  
  
  gridBase = gridBase %>% 
    mutate(Latitude= lat_grid*sd_Lat + mean_Lat,
           Depth = depth_grid*sd_Depth + mean_Depth)

  
  #################### First the summary plots ----------------------------------------------
  meanLimits = lapply(clustMemLim$regionLimits,rowMeans)
  meanLimitsDiscrete = lapply(meanLimits,function(x){
    cut(x,breaks = seq(0,1,0.1),#labels = c(<0.2,0.2,0.4,0.6,0.8,1) 
        include.lowest = T)})
  
clustRegion = parallel::mclapply(clustMemLim$regionCluster,most_frequent_random,mc.cores = 12)
  number_replicates = ncol(clustMemLim$regionCluster[[1]])
  
  
  #limitsDiscrete = 
  ##### Plotting limits
    
  gridBase %>% data.frame %>% dplyr::select(Latitude , Depth) %>% 
    mutate(limits = meanLimitsDiscrete[[1]]) %>% 
    ggplot(aes(x=Latitude,y=Depth,fill=limits))+
    geom_tile()+
    scale_y_reverse()+
    theme_minimal()+
    scale_fill_brewer(palette = 'Greys',
    )
  
  
  list_summary_Limits = lapply(seq_along(meanLimits), function(i){
    gridBase %>% select(Latitude,Depth) %>% 
      mutate(limits = meanLimits[[i]]) %>% #filter(limits>0) %>% 
      ggplot(aes(x=Latitude,y=Depth,alpha=limits))+
      geom_tile()+
      scale_y_reverse()+
      theme_minimal()+
      ggtitle(names(meanLimits)[i])
    
  })
  
  limits_faceted <- ggpubr::ggarrange(plotlist = list_summary_Limits,ncol = 4,nrow=2,
                                      common.legend = F, legend="bottom")
  
  
  list_summary_clustRegion = lapply(seq_along(clustRegion), function(i){
    
    gridBase %>% select(Latitude,Depth) %>% 
      mutate(clust_Region = factor(clustRegion[[i]]$value),
             pct = clustRegion[[i]]$frequency/number_replicates) %>% 
      ggplot(aes(x=Latitude,y=Depth,fill=clust_Region,alpha=pct))+
      geom_tile()+
      scale_y_reverse()+
      theme_minimal()+
      ggtitle(names(clustRegion)[i])
    
  })
  
  clustRegion_facetd <- ggpubr::ggarrange(
    plotlist = list_summary_clustRegion,
    ncol = 4,nrow=2,
    common.legend = F, legend="bottom")
  
  
  return(list(
    limits_faceted=limits_faceted,
    clustRegion_facetd=clustRegion_facetd
  )
  )  
  
}

most_frequent_random <- function(matrix) {
  get_most_frequent <- function(row) {
    tab <- table(row)
    max_count <- max(tab)
    most_frequent_values <- names(tab[tab == max_count])
    selected_value <- sample(most_frequent_values, 1)
    return(c(as.numeric(selected_value), max_count))
  }
  
  result <- t(apply(matrix, 1, get_most_frequent))
  df <- data.frame(value = result[, 1], frequency = result[, 2])
  return(df)
}

## Example usage
#matrix <- matrix(c(1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5), nrow = 4, byrow = TRUE)
#print(matrix)
#print(most_frequent(matrix))