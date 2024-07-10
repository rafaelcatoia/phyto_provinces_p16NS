library(dplyr)
#install.packages('rlist')

root <- rprojroot::has_file(".git/index")
datadir = root$find_file("data")
funsdir = root$find_file("functions")
savingdir = root$find_file("saved_files")

files_vec <- list.files(funsdir)

for( i in 1:length(files_vec)){
  source(root$find_file(paste0(funsdir,'/',files_vec[i])))
}

`%not_in%` <- purrr::negate(`%in%`)
######################## -------------------------------------------
## Creating the split for evaluating Cluster Size ------------------
######################## -------------------------------------------
dframe = readRDS(paste0(savingdir,'/','dfGrump_longer_filtered'))

### --------------------------------------------------------------
### Creating the B replicates of aitchison distances 
### --------------------------------------------------------------

B=500

list_AitDist = list()

idASVs = dframe %>% select(ID_ASV) %>% distinct() %>% pull
pct_colSubSample = 0.5 

min_raw_count = dframe %>% select(Raw.Sequence.Counts) %>% min()
min_raw_count = min_raw_count/1000

accepted_sample <- 1 
seedI = 5757

while(accepted_sample<=B){
  asvSubset=sample(idASVs,size = pct_colSubSample*length(idASVs))
  mat_repB <- dframe %>%
    filter(ID_ASV %in% asvSubset) %>% 
    mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
    select(SampleID,ID_ASV,Raw.Sequence.Counts) %>% distinct() %>% 
    group_by(SampleID,ID_ASV) %>% 
    summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
    tidyr::pivot_wider(id_cols = SampleID,names_from = ID_ASV ,
                       values_from = Sum_RawCounts,
                       values_fill = min_raw_count) %>%
    data.frame() %>% 
    mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
    relocate(SampleID,sort(names(.))) %>% 
    arrange(SampleID) %>% 
    select(-SampleID) %>% 
    vegan::vegdist(method = 'aitchison') %>% as.matrix() %>% 
    normalizeMatrix()
  
  mat_evalB <- dframe %>%
    filter(ID_ASV %not_in% asvSubset) %>% 
    mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
    select(SampleID,ID_ASV,Raw.Sequence.Counts) %>% distinct() %>% 
    group_by(SampleID,ID_ASV) %>% 
    summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
    tidyr::pivot_wider(id_cols = SampleID,names_from = ID_ASV ,
                       values_from = Sum_RawCounts,
                       values_fill = min_raw_count) %>%
    data.frame() %>% 
    mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
    relocate(SampleID,sort(names(.))) %>% 
    arrange(SampleID) %>% 
    select(-SampleID) %>% 
    vegan::vegdist(method = 'aitchison') %>% as.matrix()
  
  if( sum(c(dim(mat_repB),dim(mat_evalB))!=174)  == 0){
    list_AitDist = rlist::list.append(list_AitDist,list(INS=mat_repB,OOS=mat_evalB))
    accepted_sample = accepted_sample+1
  }
  cat(paste0('Iteration number ------ ::',accepted_sample,'\n'))
}

####### here we have a list of "bootstraped" normalized aitDist
####### now lets save this and move to the next step

saveRDS(list_AitDist,file = paste0(savingdir,'/','list_AitDist_inSample_OOS'))

############################### --------------------------------------------------------------------------------

