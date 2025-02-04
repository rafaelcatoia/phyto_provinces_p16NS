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


##
## ----------------------------------------------------------------------------- 10 neigh
##

list_cluster_membership_and_bounderies_mirroredLat = parallel::mclapply(
  list_AitDist[1:100],
  function(x){coloring_map(D = x,latMirrored = T,nclusters = 10,
                           gbase = grid_base,
                           list_normalized_dist = list_normalized_geo_abiotics_dists)},mc.cores = 10)

saveRDS(list_cluster_membership_and_bounderies_mirroredLat,
        file = paste0(savingdir,'/','list_cluster_membership_and_bounderies_mirroredLat'))


##
## ----------------------------------------------------------------------------- 5 neigh
##

list_cluster_membership_and_bounderies_mirroredLat_n5 = parallel::mclapply(
  list_AitDist[1:100],
  function(x){coloring_map(D = x,latMirrored = T,nclusters = 10,
                           gbase = grid_base %>% select(lat_grid,depth_grid,n_neighs01,n_neighs02,n_neighs03,n_neighs04,n_neighs05),
                           list_normalized_dist = list_normalized_geo_abiotics_dists)},mc.cores = 10)

saveRDS(list_cluster_membership_and_bounderies_mirroredLat_n5,
        file = paste0(savingdir,'/','list_cluster_membership_and_bounderies_mirroredLat_n5'))

##
## ----------------------------------------------------------------------------- 3 neigh
##

list_cluster_membership_and_bounderies_mirroredLat_n3 = parallel::mclapply(
  list_AitDist[1:100],
  function(x){coloring_map(D = x,latMirrored = T,nclusters = 10,
                           gbase = grid_base %>% select(lat_grid,depth_grid,n_neighs01,n_neighs02,n_neighs03),
                           list_normalized_dist = list_normalized_geo_abiotics_dists)},mc.cores = 10)

saveRDS(list_cluster_membership_and_bounderies_mirroredLat_n3,
        file = paste0(savingdir,'/','list_cluster_membership_and_bounderies_mirroredLat_n3'))
