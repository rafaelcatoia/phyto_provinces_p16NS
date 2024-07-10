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
list_AitDist = readRDS(paste0(savingdir,'/','list_AitDist'))
list_normalized_geo_abiotics_dists = readRDS(file = paste0(savingdir,'/','list_normalized_geo_abiotics_dists'))
list_geo_abiotics_dists = readRDS(file = paste0(savingdir,'/','list_geo_abiotics_dists'))

## there are two positions for which the dimension is not 165
df_evaluation <- data.table::rbindlist(parallel::mclapply(
  list_AitDist,function(x){
    choosing_alpha(
      D = x,
      list_normalized_dist = list_normalized_geo_abiotics_dists,
      list_dist = list_geo_abiotics_dists,
      latMirrored = T,
      vet_alpha = seq(0,1,0.01),
      K = 10)},
  mc.cores = 10))

saveRDS(object = df_evaluation,file = paste0(savingdir,'/','df_evaluation_alpha'))
#df_evaluation_alpha = readRDS(paste0(savingdir,'/','df_evaluation_alpha'))

######### 
# df_summary = bind_rows(
#   df_evaluation_lat_mirrored %>% 
#     group_by(ncluster,method,alpha,DistMetric) %>% 
#     summarise_all(.funs = mean) %>% mutate(summary_metric='Mean') %>% ungroup() %>% 
#     pivot_longer(cols = -c(ncluster,method,alpha,DistMetric,summary_metric),names_to = 'Metric'),
#   
#   df_evaluation_lat_mirrored %>% 
#     group_by(ncluster,method,alpha,DistMetric) %>% 
#     summarise_all(.funs = sd) %>% mutate(summary_metric='SD') %>% ungroup() %>% 
#     pivot_longer(cols = -c(ncluster,method,alpha,DistMetric,summary_metric),names_to = 'Metric')
# )