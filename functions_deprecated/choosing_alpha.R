############# -------------------------------------------
## script to create a function to evaluate the clustering
############# -------------------------------------------
 
# df_geo_abiotics = readRDS(file = paste0(savingdir,'/','df_geo_abiotics'))
# df_geo_abiotics = df_geo_abiotics %>% arrange(SampleID)
# list_normalized_geo_abiotics_dists = readRDS(file = paste0(savingdir,'/','list_normalized_geo_abiotics_dists'))
# list_geo_abiotics_dists = readRDS(file = paste0(savingdir,'/','list_geo_abiotics_dists'))
# list_AitDist = readRDS(file = paste0(savingdir,'/','list_AitDist'))
# D = list_AitDist[[1]]

choosing_alpha <- function(
  D,
  list_normalized_dist=list_normalized_geo_abiotics_dists, # those dists can be used to create a convex mixture of distance matrices
  list_dist = list_geo_abiotics_dists, # those will be use to perform cluser evaluation
  latMirrored=T,vet_alpha=seq(0,1,0.1),K=10){
  
  D=as.matrix(D)
  
  nrows = dim(D)[1]
  ####################################                ####################################
  #################################### Creating Phase ####################################
  ####################################                ####################################
  
  ### Settind dist matrices 
  
  if(!latMirrored){
    D_geo = list_normalized_dist$geo_Dist
  }else{
    D_geo = list_normalized_dist$geo_Dist_mirrored
  }
  
  mat_hclust <- matrix(NA,nrow=nrows,ncol=length(vet_alpha))
  mat_pam <- matrix(NA,nrow=nrows,ncol=length(vet_alpha))
  list_df_eval = list()
  for( i in 1:length(vet_alpha)){
    D_mix = (1-vet_alpha[i]) * D + vet_alpha[i]*D_geo
    mat_hclust[,i] <- cutree(hclust(d=as.dist(D_mix),method='ward.D'),k = K)
    mat_pam[,i] <- cluster::pam(x=as.dist(D_mix),k = K,cluster.only = T)
    list_df_eval[[i]]<- bind_rows(
      clustEvaluation(distMatrixEval = list_dist$geo_Dist,
                      vetClusterMembership = mat_hclust[,i]) %>%
        mutate(method='hclust',alpha=vet_alpha[i],DistMetric='Geo'), 
      clustEvaluation(distMatrixEval = list_dist$geo_Dist_mirrored,vetClusterMembership = mat_hclust[,i]) %>%
        mutate(method='hclust',alpha=vet_alpha[i], DistMetric='Geo_Mirrored') , 
      clustEvaluation(distMatrixEval = list_dist$abioticDist,vetClusterMembership = mat_hclust[,i]) %>%
        mutate(method='hclust',alpha=vet_alpha[i], DistMetric='Abiotics'), 
      clustEvaluation(distMatrixEval = list_dist$aitDist,vetClusterMembership = mat_hclust[,i]) %>%
        mutate(method='hclust',alpha=vet_alpha[i], DistMetric='Aitchison') ,
      clustEvaluation(distMatrixEval = list_dist$geo_Dist,vetClusterMembership = mat_pam[,i]) %>%
        mutate(method='pam',alpha=vet_alpha[i], DistMetric='Geo'), 
      clustEvaluation(distMatrixEval = list_dist$geo_Dist_mirrored,vetClusterMembership = mat_pam[,i]) %>%
        mutate(method='pam',alpha=vet_alpha[i], DistMetric='Geo_Mirrored'), 
      clustEvaluation(distMatrixEval = list_dist$abioticDist,vetClusterMembership = mat_pam[,i]) %>%
        mutate(method='pam',alpha=vet_alpha[i], DistMetric='Abiotics'), 
      clustEvaluation(distMatrixEval = list_dist$aitDist,vetClusterMembership = mat_pam[,i]) %>%
        mutate(method='pam',alpha=vet_alpha[i], DistMetric='Aitchison')
    )
  }
  df_out <- data.table::rbindlist(list_df_eval)
  return(df_out) 
}
