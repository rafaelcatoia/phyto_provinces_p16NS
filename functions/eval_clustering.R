############# -------------------------------------------
## script to create a function to evaluate the clustering
############# -------------------------------------------
 
# df_geo_abiotics = readRDS(file = paste0(savingdir,'/','df_geo_abiotics'))
# df_geo_abiotics = df_geo_abiotics %>% arrange(SampleID)
# list_normalized_geo_abiotics_dists = readRDS(file = paste0(savingdir,'/','list_normalized_geo_abiotics_dists'))
# list_geo_abiotics_dists = readRDS(file = paste0(savingdir,'/','list_geo_abiotics_dists'))
# list_AitDist = readRDS(file = paste0(savingdir,'/','list_AitDist'))
# D = list_AitDist[[1]]

eval_clustering <- function(
  D,
  list_normalized_dist=list_normalized_geo_abiotics_dists, # those dists can be used to create a convex mixture of distance matrices
  list_dist = list_geo_abiotics_dists, # those will be use to perform cluser evaluation
  latMirrored=T,
  vet_alpha=c(0.05,0.10,0.20),
  
  Kmax=25){
  
  D=as.matrix(D)
  
  nrows = dim(D)[1]
  ####################################                ####################################
  #################################### Creating Phase ####################################
  ####################################                ####################################
  
  ### Settind dist matrices 
  
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
  mat_hclust0 <- matrix(NA,nrow=nrows,ncol=Kmax)
  mat_hclust0[,1]<-1
  
  mat_hclust_0.05<- matrix(NA,nrow=nrows,ncol=Kmax)
  mat_hclust_0.05[,1]<-1
  
  mat_hclust_0.1<- matrix(NA,nrow=nrows,ncol=Kmax)
  mat_hclust_0.1[,1]<-1
  
  mat_hclust_0.2<- matrix(NA,nrow=nrows,ncol=Kmax)
  mat_hclust_0.2[,1]<-1
  
  for(ii in 2:Kmax){
    mat_hclust0[,ii]=cutree(hclust_0,k = ii)
    mat_hclust_0.05[,ii]=cutree(hclust_0.05,k = ii)
    mat_hclust_0.1[,ii]=cutree(hclust_0.10,k = ii)
    mat_hclust_0.2[,ii]=cutree(hclust_0.20,k = ii)
  }
  
  ####################################   pam    ####################################
  mat_pam0 <- matrix(NA,nrow=nrows,ncol=Kmax)
  mat_pam0[,1]<-1
  
  mat_pam_0.05<- matrix(NA,nrow=nrows,ncol=Kmax)
  mat_pam_0.05[,1]<-1
  
  mat_pam_0.1<- matrix(NA,nrow=nrows,ncol=Kmax)
  mat_pam_0.1[,1]<-1
  
  mat_pam_0.2<- matrix(NA,nrow=nrows,ncol=Kmax)
  mat_pam_0.2[,1]<-1

  
  for(ii in 2:Kmax){
    mat_pam0[,ii]     = cluster::pam(x=D1,k = ii,cluster.only = T)
    mat_pam_0.05[,ii] = cluster::pam(x=D2,k = ii,cluster.only = T)
    mat_pam_0.1[,ii]  = cluster::pam(x=D3,k = ii,cluster.only = T)
    mat_pam_0.2[,ii]  = cluster::pam(x=D4,k = ii,cluster.only = T)
  }
  
  
  ####################################                  ####################################
  #################################### Evaluating Phase ####################################
  ####################################                  ####################################
  Kmax = max(mat_hclust0)
  
  methodId = 'ward'
  df_hclust = bind_rows(
    bind_rows(
      evaluate_cluster(distMatrixEval = list_dist$geo_Dist %>% as.matrix,df_clusters = mat_hclust0,Kmax = Kmax) %>% mutate(method=methodId,alpha=0),
      evaluate_cluster(distMatrixEval = list_dist$geo_Dist %>% as.matrix,df_clusters = mat_hclust_0.05,Kmax = Kmax) %>% mutate(method=methodId,alpha=0.05),
      evaluate_cluster(distMatrixEval = list_dist$geo_Dist %>% as.matrix,df_clusters = mat_hclust_0.1,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.1),
      evaluate_cluster(distMatrixEval = list_dist$geo_Dist %>% as.matrix,df_clusters = mat_hclust_0.2,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.2)
    ) %>% mutate(DistMetric='Geo'),
    bind_rows(
      evaluate_cluster(distMatrixEval = list_dist$geo_Dist_mirrored %>% as.matrix,df_clusters = mat_hclust0,Kmax = Kmax) %>% mutate(method=methodId,alpha=0),
      evaluate_cluster(distMatrixEval = list_dist$geo_Dist_mirrored %>% as.matrix,df_clusters = mat_hclust_0.05,Kmax = Kmax) %>% mutate(method=methodId,alpha=0.05),
      evaluate_cluster(distMatrixEval = list_dist$geo_Dist_mirrored %>% as.matrix,df_clusters = mat_hclust_0.1,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.1),
      evaluate_cluster(distMatrixEval = list_dist$geo_Dist_mirrored %>% as.matrix,df_clusters = mat_hclust_0.2,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.2)
    ) %>% mutate(DistMetric='Geo_Mirrored') ,
    
    bind_rows(
      evaluate_cluster(distMatrixEval = list_dist$abioticDist %>% as.matrix,df_clusters = mat_hclust0,Kmax = Kmax) %>% mutate(method=methodId,alpha=0),
      evaluate_cluster(distMatrixEval = list_dist$abioticDist %>% as.matrix,df_clusters = mat_hclust_0.05,Kmax = Kmax) %>% mutate(method=methodId,alpha=0.05),
      evaluate_cluster(distMatrixEval = list_dist$abioticDist %>% as.matrix,df_clusters = mat_hclust_0.1,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.1),
      evaluate_cluster(distMatrixEval = list_dist$abioticDist %>% as.matrix,df_clusters = mat_hclust_0.2,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.2)
    ) %>% mutate(DistMetric='Abiotics'),
    bind_rows(
      evaluate_cluster(distMatrixEval = list_dist$aitDist %>% as.matrix,df_clusters = mat_hclust0,Kmax = Kmax) %>% mutate(method=methodId,alpha=0),
      evaluate_cluster(distMatrixEval = list_dist$aitDist %>% as.matrix,df_clusters = mat_hclust_0.05,Kmax = Kmax) %>% mutate(method=methodId,alpha=0.05),
      evaluate_cluster(distMatrixEval = list_dist$aitDist %>% as.matrix,df_clusters = mat_hclust_0.1,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.1),
      evaluate_cluster(distMatrixEval = list_dist$aitDist %>% as.matrix,df_clusters = mat_hclust_0.2,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.2)
    ) %>% mutate(DistMetric='Aitchison')
  )
  
  methodId ='pam'
  df_pam = bind_rows(
    bind_rows(
      evaluate_cluster(distMatrixEval = list_dist$geo_Dist %>% as.matrix,df_clusters = mat_pam0,Kmax = Kmax) %>% mutate(method=methodId,alpha=0),
      evaluate_cluster(distMatrixEval = list_dist$geo_Dist %>% as.matrix,df_clusters = mat_pam_0.05,Kmax = Kmax) %>% mutate(method=methodId,alpha=0.05),
      evaluate_cluster(distMatrixEval = list_dist$geo_Dist %>% as.matrix,df_clusters = mat_pam_0.1,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.1),
      evaluate_cluster(distMatrixEval = list_dist$geo_Dist %>% as.matrix,df_clusters = mat_pam_0.2,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.2)
    ) %>% mutate(DistMetric='Geo'),
    bind_rows(
      evaluate_cluster(distMatrixEval = list_dist$geo_Dist_mirrored %>% as.matrix,df_clusters = mat_pam0,Kmax = Kmax) %>% mutate(method=methodId,alpha=0),
      evaluate_cluster(distMatrixEval = list_dist$geo_Dist_mirrored %>% as.matrix,df_clusters = mat_pam_0.05,Kmax = Kmax) %>% mutate(method=methodId,alpha=0.05),
      evaluate_cluster(distMatrixEval = list_dist$geo_Dist_mirrored %>% as.matrix,df_clusters = mat_pam_0.1,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.1),
      evaluate_cluster(distMatrixEval = list_dist$geo_Dist_mirrored %>% as.matrix,df_clusters = mat_pam_0.2,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.2)
    ) %>% mutate(DistMetric='Geo_Mirrored'),
    bind_rows(
      evaluate_cluster(distMatrixEval = list_dist$abioticDist %>% as.matrix,df_clusters = mat_pam0,Kmax = Kmax) %>% mutate(method=methodId,alpha=0),
      evaluate_cluster(distMatrixEval = list_dist$abioticDist %>% as.matrix,df_clusters = mat_pam_0.05,Kmax = Kmax) %>% mutate(method=methodId,alpha=0.05),
      evaluate_cluster(distMatrixEval = list_dist$abioticDist %>% as.matrix,df_clusters = mat_pam_0.1,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.1),
      evaluate_cluster(distMatrixEval = list_dist$abioticDist %>% as.matrix,df_clusters = mat_pam_0.2,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.2)
    ) %>% mutate(DistMetric='Abiotics'),
    bind_rows(
      evaluate_cluster(distMatrixEval = list_dist$aitDist %>% as.matrix,df_clusters = mat_pam0,Kmax = Kmax) %>% mutate(method=methodId,alpha=0),
      evaluate_cluster(distMatrixEval = list_dist$aitDist %>% as.matrix,df_clusters = mat_pam_0.05,Kmax = Kmax) %>% mutate(method=methodId,alpha=0.05),
      evaluate_cluster(distMatrixEval = list_dist$aitDist %>% as.matrix,df_clusters = mat_pam_0.1,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.1),
      evaluate_cluster(distMatrixEval = list_dist$aitDist %>% as.matrix,df_clusters = mat_pam_0.2,Kmax = Kmax)  %>% mutate(method=methodId,alpha=0.2)
    ) %>% mutate(DistMetric='Aitchison')
  )
  
  df_out = bind_rows(df_hclust,df_pam)

  return(df_out) 
}

