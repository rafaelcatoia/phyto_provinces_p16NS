################################## 
library(dplyr) ; library(tidyr)

root <- rprojroot::has_file(".git/index")
datadir = root$find_file("data")
funsdir = root$find_file("functions")
savingdir = root$find_file("saved_files")

files_vec <- list.files(funsdir)

for( i in 1:length(files_vec)){
  source(root$find_file(paste0(funsdir,'/',files_vec[i])))
}

###############################
## Loading
###############################
df_geo_abiotics = readRDS(file = paste0(savingdir,'/','df_geo_abiotics'))
dfGrump_longer_filtered = readRDS(paste0(savingdir,'/','dfGrump_longer_filtered'))
list_AitDist = readRDS(paste0(savingdir,'/','list_AitDist'))
list_normalized_geo_abiotics_dists = readRDS(file = paste0(savingdir,'/','list_normalized_geo_abiotics_dists'))
list_geo_abiotics_dists = readRDS(file = paste0(savingdir,'/','list_geo_abiotics_dists'))
grid_base = readRDS(file = paste0(savingdir,'/','grid_base'))

#############################
######## ---------------------------------------------------------------------
## Creating cluster label ####### 
nclusters = 10
D = list_geo_abiotics_dists$aitDist %>% normalizeMatrix()
vet_alpha = c(0,0.05,0.1,0.2)
D1 = as.dist(D)
D2 = as.dist( (1-as.numeric(vet_alpha[1]))*D + as.numeric(vet_alpha[1])*list_normalized_geo_abiotics_dists$geo_Dist_mirrored)
D3 = as.dist( (1-as.numeric(vet_alpha[2]))*D + as.numeric(vet_alpha[2])*list_normalized_geo_abiotics_dists$geo_Dist_mirrored)
D4 = as.dist( (1-as.numeric(vet_alpha[3]))*D + as.numeric(vet_alpha[3])*list_normalized_geo_abiotics_dists$geo_Dist_mirrored)

####################################   Ward    ####################################
#################################### Creating hclust objs
hclust_0 = hclust(d=D1,method='ward.D')
hclust_0.05 = hclust(d=D2,method='ward.D')
hclust_0.10 = hclust(d=D3,method='ward.D')
hclust_0.20 = hclust(d=D4,method='ward.D')

#################################### With geodist 
mat_cluster_membership_label <- cbind.data.frame(
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
######## ---------------------------------------------------------------------

#list_cluster_membership_and_bounderies = lapply(list_AitDist_clean[1:100],coloring_map)
#saveRDS(list_cluster_membership_and_bounderies,file = paste0(savingdir,'/','list_cluster_membership_and_bounderies'))

list_AitDist <- list_AitDist[1:100]

##
## ----------------------------------------------------------------------------- 10 neigh
##


list_cluster_membership_and_bounderies_mirroredLat_matched = parallel::mclapply(
  list_AitDist,
  function(x){coloring_map_matching(
    D = x,latMirrored = T,nclusters = 10,
    trueClusterMembership = mat_cluster_membership_label,
    list_dist_toMatch=list_geo_abiotics_dists,
    gbase = grid_base,
    list_normalized_dist = list_normalized_geo_abiotics_dists)},mc.cores = 10)

saveRDS(list_cluster_membership_and_bounderies_mirroredLat_matched,
        file = paste0(savingdir,'/','list_cluster_membership_and_bounderies_mirroredLat_matched'))

##
## ----------------------------------------------------------------------------- 5 neigh
##


list_cluster_membership_and_bounderies_mirroredLat_matched_n5 = parallel::mclapply(
  list_AitDist,
  function(x){coloring_map_matching(
    D = x,latMirrored = T,nclusters = 10,
    trueClusterMembership = mat_cluster_membership_label,
    list_dist_toMatch=list_geo_abiotics_dists,
    gbase = grid_base %>% select(lat_grid,depth_grid,n_neighs01,n_neighs02,n_neighs03,n_neighs04,n_neighs05),
    list_normalized_dist = list_normalized_geo_abiotics_dists)},mc.cores = 10)

saveRDS(list_cluster_membership_and_bounderies_mirroredLat_matched_n5,
        file = paste0(savingdir,'/','list_cluster_membership_and_bounderies_mirroredLat_matched_n5'))


##
## ----------------------------------------------------------------------------- 3 neigh
##

## three neigh
list_cluster_membership_and_bounderies_mirroredLat_matched_n3 = parallel::mclapply(
  list_AitDist,
  function(x){coloring_map_matching(
    D = x,latMirrored = T,nclusters = 10,
    trueClusterMembership = mat_cluster_membership_label,
    list_dist_toMatch=list_geo_abiotics_dists,
    gbase = grid_base %>% select(lat_grid,depth_grid,n_neighs01,n_neighs02,n_neighs03),
    list_normalized_dist = list_normalized_geo_abiotics_dists)},mc.cores = 10)

saveRDS(list_cluster_membership_and_bounderies_mirroredLat_matched_n3,
        file = paste0(savingdir,'/','list_cluster_membership_and_bounderies_mirroredLat_matched_n3'))
