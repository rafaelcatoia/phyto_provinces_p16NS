############## ------ function to get the regions

## here we have a grid that contains the closest 5 samples to each lat,depth
# grid_base <- readRDS(paste0(savingdir,'/','grid_base'))
# D = list_AitDist[[1]]
# 

## here are the clusters that we will use as the TRUE labels
coloring_map_matching <- function(
    D,
    trueClusterMembership,
    gbase=grid_base,
    list_normalized_dist=list_normalized_geo_abiotics_dists, # those dists can be used to create a convex mixture of distance matrices
    list_dist_toMatch=list_geo_abiotics_dists,
    vet_alpha=c(0.05,0.10,0.20),
    latMirrored=T,
    nclusters=10){
  
  D=as.matrix(D)
  
  nrows = dim(D)[1]
  ####################################                ####################################
  #################################### Creating Phase ####################################
  ####################################                ####################################
  list_normalized_dist$abioticDist
  
  if(!latMirrored){
    ### Settind dist matrices 
    D1 = as.dist(D)
    D2= as.dist( (1-as.numeric(vet_alpha[1]))*D + as.numeric(vet_alpha[1])*list_normalized_dist$geo_Dist)
    D3= as.dist( (1-as.numeric(vet_alpha[2]))*D + as.numeric(vet_alpha[2])*list_normalized_dist$geo_Dist)
    D4= as.dist( (1-as.numeric(vet_alpha[3]))*D + as.numeric(vet_alpha[3])*list_normalized_dist$geo_Dist)
  }else{
    ### Settind dist matrices 
    D1 = as.dist(D)
    D2 = as.dist( (1-as.numeric(vet_alpha[1]))*D + as.numeric(vet_alpha[1])*list_normalized_dist$geo_Dist_mirrored)
    D3 = as.dist( (1-as.numeric(vet_alpha[2]))*D + as.numeric(vet_alpha[2])*list_normalized_dist$geo_Dist_mirrored)
    D4 = as.dist( (1-as.numeric(vet_alpha[3]))*D + as.numeric(vet_alpha[3])*list_normalized_dist$geo_Dist_mirrored)
  }
  
  ####################################   Ward    ####################################
  #################################### Creating hclust objs
  hclust_0 = hclust(d=D1,method='ward.D')
  hclust_0.05 = hclust(d=D2,method='ward.D')
  hclust_0.10 = hclust(d=D3,method='ward.D')
  hclust_0.20 = hclust(d=D4,method='ward.D')
  
  #################################### With geodist 
  mat_cluster_membership <- cbind.data.frame(
    v1_hclust=cutree(hclust_0,k = nclusters),
    v2_hclust=cutree(hclust_0.05,k = nclusters),
    v3_hclust=cutree(hclust_0.10,k = nclusters),
    v4_hclust=cutree(hclust_0.20,k = nclusters),
    ####################################   pam    ####################################
    v1_pam=cluster::pam(x=D1,k = nclusters,cluster.only = T),
    v2_pam=cluster::pam(x=D2,k = nclusters,cluster.only = T),
    v3_pam=cluster::pam(x=D3,k = nclusters,cluster.only = T),
    v4_pam=cluster::pam(x=D4,k = nclusters,cluster.only = T)
  )
  
  ####################################################################################
  ## here we will match the clusters -------------------------------------------------
  ####################################################################################
  for(i in 1:ncol(mat_cluster_membership)){
    mat_cluster_membership[,i]<-
      matching_clusters(
      DistMatrix = list_dist_toMatch$geo_Dist_mirrored,
      trueLabel = trueClusterMembership[,i],
      proposedLabel = mat_cluster_membership[,i])
  }
  ####################################################################################
  out_list <- list()
  for(i in 1:ncol(mat_cluster_membership)){
    nneigh_colnames <- gbase %>% select(starts_with('n_neig')) %>% colnames()
    aux = 1:length(nneigh_colnames)
    custDominates =  ifelse(aux<10,
                            paste0('ClustDom0',aux),
                            paste0('ClustDom',aux))
    matCluster = matrix(NA,nrow=nrow(gbase),ncol=length(nneigh_colnames))
    
    
    for(j in 1:length(nneigh_colnames)){
      clustDom = grid_base %>% select(any_of(nneigh_colnames[j])) %>% pull()
      matCluster[,j] <- getClustMembership(clustDom,vector_cluster_mebership = mat_cluster_membership[,i])
    }
    
    clustRegion = most_frequent_value(matCluster)
    gbase = gbase %>% mutate(ClustRegion = clustRegion)
    grid_base_matrix = gbase %>% 
      tidyr::pivot_wider(id_cols = lat_grid, names_from = depth_grid,
                         values_from = ClustRegion) %>% 
      as.matrix()
    grid_base_matrix_limits = check_adjacent_values(grid_base_matrix[,-1]) %>% data.frame() %>% 
      tidyr::pivot_longer(everything())
    out_list[[i]]<-data.frame(clustRegion=clustRegion,limits = grid_base_matrix_limits$value)
  }
  
  names(out_list) <- c(paste0('ward_',c(0,0.05,0.1,0.2)),paste0('pam_',c(0,0.05,0.1,0.2)))
  
  return(out_list)
}


