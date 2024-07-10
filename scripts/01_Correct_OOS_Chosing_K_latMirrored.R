library(dplyr) ; library(tidyr)

root <- rprojroot::has_file(".git/index")
datadir = root$find_file("data")
funsdir = root$find_file("functions")
savingdir = root$find_file("saved_files")

library(parallel) ; library(dplyr)

files_vec <- list.files(funsdir)

for( i in 1:length(files_vec)){
  source(root$find_file(paste0(funsdir,'/',files_vec[i])))
}

###############################
## Loading
###############################
df_geo_abiotics = readRDS(file = paste0(savingdir,'/','df_geo_abiotics'))
dfGrump_longer_filtered = readRDS(paste0(savingdir,'/','dfGrump_longer_filtered'))
list_AitDist = readRDS(file = paste0(savingdir,'/','list_AitDist_inSample_OOS'))
list_normalized_geo_abiotics_dists = readRDS(file = paste0(savingdir,'/','list_normalized_geo_abiotics_dists'))
list_geo_abiotics_dists = readRDS(file = paste0(savingdir,'/','list_geo_abiotics_dists'))

## there are two positions for which the dimension is not 165
df_evaluation <- data.table::rbindlist(parallel::mclapply(
  list_AitDist,function(x){eval_clustering_OOS(
    D = x,
    list_normalized_dist = list_normalized_geo_abiotics_dists,
    list_dist = list_geo_abiotics_dists,
    latMirrored = T)},mc.cores = 10))

###### THE WARNINGS ARE SUPOSED TO HAPPEN!


saveRDS(df_evaluation,paste0(savingdir,'/','df_evaluation_lat_mirrored_OOS'))
#df_evaluation = readRDS(paste0(savingdir,'/','df_evaluation'))

